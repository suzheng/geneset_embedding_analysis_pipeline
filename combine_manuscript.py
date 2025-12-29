import re
import os

def md_to_html(text):
    # Handle <it> tag
    text = text.replace('<it>', '<em>').replace('</it>', '</em>')
    
    # Handle Headers
    lines = text.split('\n')
    html_lines = []
    
    in_paragraph = False
    
    for line in lines:
        line = line.strip()
        if not line:
            if in_paragraph:
                html_lines.append('</p>')
                in_paragraph = False
            continue
            
        # Headers
        if line.startswith('# '):
            if in_paragraph: html_lines.append('</p>'); in_paragraph = False
            html_lines.append(f'<h1>{line[2:]}</h1>')
        elif line.startswith('## '):
            if in_paragraph: html_lines.append('</p>'); in_paragraph = False
            html_lines.append(f'<h2>{line[3:]}</h2>')
        elif line.startswith('### '):
            if in_paragraph: html_lines.append('</p>'); in_paragraph = False
            html_lines.append(f'<h3>{line[4:]}</h3>')
        else:
            # Inline formatting
            # Bold
            line = re.sub(r'\*\*(.*?)\*\*', r'<strong>\1</strong>', line)
            # Italic (simple check)
            line = re.sub(r'\*([^\*]+)\*', r'<em>\1</em>', line)
            
            # Links [[url]] -> <a href="url">[Citation]</a> or similar? 
            # The user text has [[url]]. Let's make it clickable.
            # Pattern: [[url]]
            def replace_link(match):
                url = match.group(1)
                return f'<a href="{url}">[{url}]</a>'
            
            line = re.sub(r'\[\[(https?://[^\]]+)\]\]', replace_link, line)
            
            # Standard markdown links [text](url)
            line = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', line)

            if not in_paragraph:
                html_lines.append('<p>')
                in_paragraph = True
            
            html_lines.append(line)
            
    if in_paragraph:
        html_lines.append('</p>')
        
    return '\n'.join(html_lines)

def read_file(path):
    with open(path, 'r') as f:
        return f.read()

def main():
    base_dir = "docs/manuscript"
    
    # Read files
    abstract_intro = read_file(os.path.join(base_dir, "abstract_and_introduction.md"))
    methods = read_file(os.path.join(base_dir, "methods.md"))
    results = read_file(os.path.join(base_dir, "results.html"))
    discussion = read_file(os.path.join(base_dir, "discussion.md"))
    
    # Convert MD to HTML
    abstract_intro_html = md_to_html(abstract_intro)
    methods_html = md_to_html(methods)
    discussion_html = md_to_html(discussion)
    
    # Process Results HTML
    # Extract CSS
    style_match = re.search(r'<style>(.*?)</style>', results, re.DOTALL)
    css = style_match.group(1) if style_match else ""
    
    # Extract Body content (between <body> and </body>)
    body_match = re.search(r'<body>(.*?)</body>', results, re.DOTALL)
    results_body = body_match.group(1) if body_match else results
    
    # We might want to remove the <h1>Results</h1> from results_body if we are restructuring,
    # but likely keeping it is fine as other sections have headers too.
    # Abstract/Intro MD has # Abstract and # Introduction -> <h1>.
    # So Results having <h1>Results</h1> is consistent.
    
    # Combine
    final_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Identifying Immune Pathways and Drug Targets in Sepsis Through Gene Expression Profiling</title>
    <style>
        {css}
        /* Add some extra spacing for distinct sections if needed */
        section {{
            margin-bottom: 0.5in;
        }}
        .main-title {{
            font-size: 24pt;
            font-weight: bold;
            text-align: center;
            margin-bottom: 0.5in;
            margin-top: 0.5in;
        }}
    </style>
</head>
<body>

<h1 class="main-title">Identifying Immune Pathways and Drug Targets in Sepsis Through Gene Expression Profiling</h1>

<section id="abstract-intro">
{abstract_intro_html}
</section>

<section id="methods">
{methods_html}
</section>

<section id="results">
{results_body}
</section>

<section id="discussion">
{discussion_html}
</section>

</body>
</html>"""

    output_path = os.path.join(base_dir, "manuscript_combined.html")
    with open(output_path, 'w') as f:
        f.write(final_html)
    
    print(f"StartLine:1")
    print(f"EndLine:{len(final_html.splitlines())}") # Just for info
    print(f"Created {output_path}")

if __name__ == "__main__":
    main()
