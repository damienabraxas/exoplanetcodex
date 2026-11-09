import re, sys, csv
S = sys.argv[1]
NUM = r"[−\-+]?\d+\s?\.\d+"
def f(x): return float(x.replace("−", "-").replace(" ", ""))
rows = []
# Paper I (Na-Ca): element headers lose the ion numeral; 1st occurrence = I, 2nd = II
seen = {}; el = None; ion = None
for line in open(f"{S}/s15a.txt"):
    line = re.sub(r"[^\x20-\x7e\u2212\u27e8\u27e9]", " ", line).strip()
    if line.startswith("Table 3"): break
    m = re.fullmatch(r"([A-Z][a-z]?)", line)
    if m:
        el = m.group(1); seen[el] = seen.get(el, 0) + 1; ion = "I" if seen[el] == 1 else "II"; continue
    m = re.match(rf"^(\d{{3,4}}\.\d{{4}})[a-z]?\s+({NUM})\s+({NUM})\s+\d+\s+({NUM})\s+\d+\*?\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})", line)
    if m and el:
        g = m.groups()
        rows.append(dict(paper="Scott+2015a", element=el, ion=ion, wavelength_nm=g[0], eexc=f(g[1]), loggf=f(g[2]),
                         A_3D=f(g[4]), A_mean3D=f(g[5]), A_HM=f(g[6]), A_MARCS=f(g[7]), A_MISS=f(g[8])))
# Paper II (Sc-Ni): explicit "Ti i" headers; lambda on its own line, numbers on a later line
el = ion = None; lam = None
for line in open(f"{S}/s15b.txt"):
    line = re.sub(r"[^\x20-\x7e\u2212\u27e8\u27e9]", " ", line).strip()
    m = re.fullmatch(r"([A-Z][a-z]?) (i|ii)", line)
    if m:
        el, ion = m.group(1), m.group(2).upper(); lam = None; continue
    m = re.match(r"^(\d{3,4}\.\d{4})[a-z]?(\s|$)", line)
    if m: lam = m.group(1)
    m = re.search(rf"({NUM})\s+({NUM})\*?\s+\d+\s+({NUM})\s+\d+\*?\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+({NUM})(\s|$)", line)
    if m and el and lam:
        g = m.groups()
        rows.append(dict(paper="Scott+2015b", element=el, ion=ion, wavelength_nm=lam, eexc=f(g[0]), loggf=f(g[1]),
                         A_3D=f(g[3]), A_mean3D=f(g[4]), A_HM=f(g[5]), A_MARCS=f(g[6]), A_MISS=f(g[7])))
        lam = None
w = csv.DictWriter(open(f"{S}/scott2015.csv", "w"), fieldnames=list(rows[0]))
w.writeheader(); w.writerows(rows)
from collections import Counter
print(len(rows), sorted(Counter((r["element"], r["ion"]) for r in rows).items()))
