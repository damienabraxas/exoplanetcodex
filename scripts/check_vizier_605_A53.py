"""
scripts/check_vizier_605_A53.py
================================
Quick check: does VizieR J/A+A/605/A53 have machine-readable Ca NLTE tables?
(Mashonkina, Sitnova & Belyaev 2017, A&A 605, A53)
"""
import requests

url = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/605/A53/"
r = requests.get(url, timeout=30)
print(f"Status: {r.status_code}")
if r.status_code == 200:
    print(r.text[:2000])
else:
    print("Not found at FTP path — trying viz-bin/cat:")
    alt = "https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/605/A53"
    r2 = requests.get(alt, timeout=30)
    print(f"Alternate status: {r2.status_code}  len={len(r2.text)}")
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(r2.text, 'html.parser')
    print(soup.title.text if soup.title else r2.text[:500])
