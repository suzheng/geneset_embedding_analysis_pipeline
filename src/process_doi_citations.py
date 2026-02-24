#!/usr/bin/env python3
"""
Process inline DOI citations in text/markdown and export RIS.

Input citation format (within the text):
  [https://doi.org/DOI_1, https://doi.org/DOI_2]

This script will:
  1) Replace each bracketed DOI list with:
     {Author1, Year1, first three words of Title1; Author2, Year2, first three words of Title2}
  2) Resolve DOIs via Crossref API (https://api.crossref.org/works/{doi})
  3) Write a RIS file containing all successfully resolved DOIs (deduplicated)

Notes:
  - Some DOIs may be invalid/unresolvable; those will be kept as placeholders in the text output
    and skipped in RIS output.
  - Uses tqdm for progress display.
"""

from __future__ import annotations

import argparse
import html
import json
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote
from urllib.request import Request, urlopen

try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm is not available
    def tqdm(iterable, *args, **kwargs):
        return iterable


DOI_CITATION_BLOCK_RE = re.compile(
    r"""
    \[
      (?P<inner>
        \s*https?://doi\.org/[^\s,\]]+\s*
        (?:,\s*https?://doi\.org/[^\s,\]]+\s*)*
      )
    \]
    """,
    re.VERBOSE,
)

DOI_URL_RE = re.compile(r"https?://doi\.org/(?P<doi>[^\s,\]]+)", re.IGNORECASE)


@dataclass(frozen=True)
class ResolveResult:
    doi: str
    ok: bool
    message: str
    meta: Optional[Dict[str, Any]] = None


def _load_cache(cache_path: Optional[Path]) -> Dict[str, Dict[str, Any]]:
    if not cache_path:
        return {}
    if not cache_path.exists():
        return {}
    try:
        with cache_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict):
            # Expect: {doi: crossref_message_dict}
            return {str(k): v for k, v in data.items() if isinstance(v, dict)}
    except Exception:
        # Cache is optional; ignore corruption and proceed without it.
        return {}
    return {}


def _save_cache(cache_path: Optional[Path], cache: Dict[str, Dict[str, Any]]) -> None:
    if not cache_path:
        return
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(cache_path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)
    tmp.replace(cache_path)


def normalize_doi(raw: str) -> str:
    """
    Normalize a DOI extracted from a https://doi.org/... URL.
    """
    doi = unquote(raw).strip()
    # Strip trailing punctuation that commonly appears in prose.
    doi = re.sub(r"[.,;:]+$", "", doi)
    return doi


def extract_dois_from_block(inner: str) -> List[str]:
    dois: List[str] = []
    for m in DOI_URL_RE.finditer(inner):
        doi = normalize_doi(m.group("doi"))
        if doi:
            dois.append(doi)
    return dois


def crossref_resolve(
    doi: str,
    *,
    mailto: Optional[str],
    timeout_s: float,
) -> ResolveResult:
    """
    Resolve DOI metadata via Crossref.
    Returns Crossref "message" dict on success.
    """
    # Crossref expects the DOI in the path; must be URL-encoded.
    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    ua = "sports-meta-analysis/0.1 (citation processing script)"
    if mailto:
        ua = f"{ua} (mailto:{mailto})"
    req = Request(url, headers={"User-Agent": ua, "Accept": "application/json"})
    try:
        with urlopen(req, timeout=timeout_s) as resp:
            payload = resp.read().decode("utf-8", errors="replace")
        data = json.loads(payload)
        if not isinstance(data, dict) or "message" not in data:
            return ResolveResult(doi=doi, ok=False, message="Crossref response missing 'message'")
        msg = data.get("message")
        if not isinstance(msg, dict):
            return ResolveResult(doi=doi, ok=False, message="Crossref 'message' is not an object")
        return ResolveResult(doi=doi, ok=True, message="ok", meta=msg)
    except HTTPError as e:
        return ResolveResult(doi=doi, ok=False, message=f"HTTP {e.code}: {e.reason}")
    except URLError as e:
        return ResolveResult(doi=doi, ok=False, message=f"URL error: {e.reason}")
    except json.JSONDecodeError:
        return ResolveResult(doi=doi, ok=False, message="Failed to decode JSON from Crossref")
    except Exception as e:
        return ResolveResult(doi=doi, ok=False, message=f"Unexpected error: {e}")


def _pick_year(msg: Dict[str, Any]) -> Optional[int]:
    for key in ("issued", "published-print", "published-online", "created"):
        obj = msg.get(key)
        if isinstance(obj, dict):
            dp = obj.get("date-parts")
            if isinstance(dp, list) and dp and isinstance(dp[0], list) and dp[0]:
                y = dp[0][0]
                if isinstance(y, int):
                    return y
                if isinstance(y, str) and y.isdigit():
                    return int(y)
    return None


def _first_author_family(msg: Dict[str, Any]) -> str:
    authors = msg.get("author")
    if isinstance(authors, list) and authors:
        a0 = authors[0]
        if isinstance(a0, dict):
            fam = a0.get("family")
            if isinstance(fam, str) and fam.strip():
                return fam.strip()
            given = a0.get("given")
            if isinstance(given, str) and given.strip():
                return given.strip()
    return "Unknown"


def _first_three_words(title: str) -> str:
    t = re.sub(r"<[^>]+>", "", title or "").strip()
    parts = t.split()
    if not parts:
        return "Untitled"
    return " ".join(parts[:3])


def format_inline_citation(msg: Dict[str, Any]) -> str:
    author = _first_author_family(msg)
    year = _pick_year(msg)
    title_list = msg.get("title")
    title = title_list[0] if isinstance(title_list, list) and title_list and isinstance(title_list[0], str) else ""
    three = _first_three_words(title)
    year_str = str(year) if year is not None else "n.d."
    return f"{author}, {year_str}, {three}"


def _ris_type(msg: Dict[str, Any]) -> str:
    t = msg.get("type")
    if t == "journal-article":
        return "JOUR"
    if t in ("book-chapter", "book-section"):
        return "CHAP"
    if t == "book":
        return "BOOK"
    if t == "proceedings-article":
        return "CPAPER"
    return "GEN"


def _ris_escape(v: str) -> str:
    # Decode HTML entities (e.g., &amp; -> &) for RIS compatibility with importers like EndNote
    v = html.unescape(v or "")
    return v.replace("\r", " ").replace("\n", " ").strip()


def to_ris(msg: Dict[str, Any], doi: str) -> str:
    lines: List[str] = []
    lines.append(f"TY  - {_ris_type(msg)}")

    # Title
    title_list = msg.get("title")
    if isinstance(title_list, list) and title_list and isinstance(title_list[0], str):
        lines.append(f"TI  - {_ris_escape(title_list[0])}")

    # Authors
    authors = msg.get("author")
    if isinstance(authors, list):
        for a in authors:
            if not isinstance(a, dict):
                continue
            family = a.get("family")
            given = a.get("given")
            if isinstance(family, str) and isinstance(given, str) and family.strip() and given.strip():
                lines.append(f"AU  - {_ris_escape(family)}, {_ris_escape(given)}")
            elif isinstance(family, str) and family.strip():
                lines.append(f"AU  - {_ris_escape(family)}")

    # Year
    year = _pick_year(msg)
    if year is not None:
        lines.append(f"PY  - {year}")

    # Journal / container
    ct = msg.get("container-title")
    if isinstance(ct, list) and ct and isinstance(ct[0], str) and ct[0].strip():
        lines.append(f"JO  - {_ris_escape(ct[0])}")

    # Volume / issue
    vol = msg.get("volume")
    iss = msg.get("issue")
    if isinstance(vol, str) and vol.strip():
        lines.append(f"VL  - {_ris_escape(vol)}")
    if isinstance(iss, str) and iss.strip():
        lines.append(f"IS  - {_ris_escape(iss)}")

    # Pages
    page = msg.get("page")
    if isinstance(page, str) and page.strip():
        # Best-effort split "123-130"
        if "-" in page:
            sp, ep = page.split("-", 1)
            if sp.strip():
                lines.append(f"SP  - {_ris_escape(sp)}")
            if ep.strip():
                lines.append(f"EP  - {_ris_escape(ep)}")
        else:
            lines.append(f"SP  - {_ris_escape(page)}")

    # DOI & URL
    lines.append(f"DO  - {_ris_escape(doi)}")
    lines.append(f"UR  - https://doi.org/{_ris_escape(doi)}")
    lines.append("ER  -")
    return "\n".join(lines) + "\n"


def process_text(
    text: str,
    *,
    cache: Dict[str, Dict[str, Any]],
    mailto: Optional[str],
    timeout_s: float,
    sleep_s: float,
    show_unresolved: bool = False,
) -> Tuple[str, Dict[str, ResolveResult]]:
    """
    Returns:
      - processed text
      - per-doi resolve results (deduped by DOI)
    
    Args:
      show_unresolved: If True, show "Unresolved, n.d., {doi}" for failed DOIs.
                       If False (default), remove unresolved citations entirely.
    """
    # First, extract all unique DOIs from the text
    all_dois: List[str] = []
    for m in DOI_CITATION_BLOCK_RE.finditer(text):
        inner = m.group("inner")
        dois = extract_dois_from_block(inner)
        all_dois.extend(dois)
    
    # Deduplicate while preserving order
    unique_dois: List[str] = []
    seen: Set[str] = set()
    for doi in all_dois:
        if doi not in seen:
            unique_dois.append(doi)
            seen.add(doi)
    
    # Resolve all DOIs with progress bar
    resolved: Dict[str, ResolveResult] = {}
    
    print(f"Found {len(unique_dois)} unique DOIs to resolve...")
    for doi in tqdm(unique_dois, desc="Resolving DOIs", unit="DOI"):
        if doi in cache:
            rr = ResolveResult(doi=doi, ok=True, message="ok (cache)", meta=cache[doi])
            resolved[doi] = rr
        else:
            rr = crossref_resolve(doi, mailto=mailto, timeout_s=timeout_s)
            resolved[doi] = rr
            if rr.ok and rr.meta is not None:
                cache[doi] = rr.meta
            if sleep_s > 0:
                time.sleep(sleep_s)
    
    # Now replace citations in the text
    def repl(m: re.Match) -> str:
        inner = m.group("inner")
        dois = extract_dois_from_block(inner)
        parts: List[str] = []
        for doi in dois:
            rr = resolved.get(doi)
            if rr and rr.ok and rr.meta is not None:
                parts.append(format_inline_citation(rr.meta))
            elif show_unresolved:
                # Keep a readable placeholder in the rewritten manuscript.
                parts.append(f"Unresolved, n.d., {doi}")
        # If no parts (all unresolved and show_unresolved=False), return empty string
        if not parts:
            return ""
        return "{" + "; ".join(parts) + "}"

    processed = DOI_CITATION_BLOCK_RE.sub(repl, text)
    return processed, resolved


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Replace inline DOI citation blocks with {Author, Year, Title...} and export RIS via Crossref."
    )
    p.add_argument("--input", required=True, help="Input .md/.txt file path")
    p.add_argument("--output", required=True, help="Output processed file path")
    p.add_argument("--ris", default="", help="Output RIS file path (resolved DOIs only) - DEPRECATED, use --ris-folder instead")
    p.add_argument("--ris-folder", default="", help="Output folder for individual RIS files (one per citation)")
    p.add_argument("--cache", default="", help="Optional JSON cache path for Crossref 'message' objects")
    p.add_argument("--mailto", default="", help="Optional email for Crossref polite User-Agent (recommended)")
    p.add_argument("--timeout", type=float, default=20.0, help="HTTP timeout seconds (default: 20)")
    p.add_argument("--sleep", type=float, default=0.1, help="Sleep seconds between Crossref calls (default: 0.1)")
    p.add_argument(
        "--show-unresolved",
        action="store_true",
        help="Show 'Unresolved, n.d., {doi}' for failed DOIs. Default: remove unresolved citations entirely.",
    )
    args = p.parse_args(argv)

    in_path = Path(args.input)
    out_path = Path(args.output)
    ris_path = Path(args.ris) if args.ris else None
    ris_folder = Path(args.ris_folder) if args.ris_folder else None
    cache_path = Path(args.cache) if args.cache else None
    mailto = args.mailto.strip() or None

    if not in_path.exists():
        print(f"ERROR: input file not found: {in_path}", file=sys.stderr)
        return 2

    if not ris_path and not ris_folder:
        print(f"ERROR: must provide either --ris or --ris-folder", file=sys.stderr)
        return 2

    cache = _load_cache(cache_path)

    print(f"Reading input file: {in_path}")
    text = in_path.read_text(encoding="utf-8")
    
    processed, resolved = process_text(
        text,
        cache=cache,
        mailto=mailto,
        timeout_s=float(args.timeout),
        sleep_s=float(args.sleep),
        show_unresolved=args.show_unresolved,
    )

    print(f"\nWriting output file: {out_path}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(processed, encoding="utf-8")

    # Write RIS files
    ok_count = 0
    bad_count = 0

    if ris_folder:
        print(f"Writing individual RIS files to: {ris_folder}")
        ris_folder.mkdir(parents=True, exist_ok=True)
        for doi in tqdm(sorted(resolved.keys()), desc="Writing RIS files", unit="file"):
            rr = resolved[doi]
            if rr.ok and rr.meta is not None:
                # Create safe filename from DOI
                safe_filename = doi.replace("/", "_").replace(".", "_") + ".ris"
                ris_file = ris_folder / safe_filename
                ris_file.write_text(to_ris(rr.meta, doi), encoding="utf-8")
                ok_count += 1
            else:
                bad_count += 1
    elif ris_path:
        print(f"Writing combined RIS file: {ris_path}")
        ris_path.parent.mkdir(parents=True, exist_ok=True)
        with ris_path.open("w", encoding="utf-8") as f:
            for doi in tqdm(sorted(resolved.keys()), desc="Writing RIS", unit="entry"):
                rr = resolved[doi]
                if rr.ok and rr.meta is not None:
                    f.write(to_ris(rr.meta, doi))
                    ok_count += 1
                else:
                    bad_count += 1

    _save_cache(cache_path, cache)

    print(f"\n✓ Processed: {in_path} -> {out_path}")
    if ris_folder:
        print(f"✓ RIS files written to: {ris_folder}")
    elif ris_path:
        print(f"✓ RIS written: {ris_path}")
    print(f"✓ DOIs resolved: {ok_count}, unresolved: {bad_count}, total seen: {len(resolved)}")

    # Print unresolved details to stderr for easier debugging.
    if bad_count:
        print("\nUnresolved DOIs:", file=sys.stderr)
        for doi in sorted(resolved.keys()):
            rr = resolved[doi]
            if not rr.ok:
                print(f"- {doi}: {rr.message}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

