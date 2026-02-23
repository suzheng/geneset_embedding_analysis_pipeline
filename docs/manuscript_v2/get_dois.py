import sys, json, urllib.request, urllib.parse
query = urllib.parse.quote(" ".join(sys.argv[1:]))
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term={query}"
try:
    req = urllib.request.urlopen(url)
    data = json.loads(req.read())
    idlist = data.get("esearchresult", {}).get("idlist", [])[:5]
    if not idlist:
        print("No results")
        sys.exit(0)
    ids = ",".join(idlist)
    sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id={ids}"
    sum_req = urllib.request.urlopen(sum_url)
    sum_data = json.loads(sum_req.read())
    for i in idlist:
        article = sum_data["result"][i]
        title = article.get("title", "No Title")
        dois = [a["value"] for a in article.get("articleids", []) if a["idtype"] == "doi"]
        doi = dois[0] if dois else "No DOI"
        print(f"Title: {title}\nDOI: {doi}\n")
except Exception as e:
    print(e)
