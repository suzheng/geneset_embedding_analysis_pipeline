search_epmc() {
  query=$1
  curl -s -G --data-urlencode "query=$query" \
    --data-urlencode "format=json" \
    --data-urlencode "resultType=core" \
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search" | \
    python3 -c "
import sys, json
data = json.load(sys.stdin)
for res in data.get('resultList', {}).get('result', [])[:5]:
    print('Title:', res.get('title'))
    print('DOI:', res.get('doi'))
    print('Abstract:', res.get('abstractText', 'No abstract').replace('\n', ' '))
    print('-'*80)
"
}
search_epmc "$1"
