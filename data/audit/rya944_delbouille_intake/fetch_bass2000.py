"""RYA-944: acquire the BASS2000 visible arm (Jungfraujoch/Delbouille) 3000-10000 A.

Chunked pull with retry. Each chunk is kept as its own raw file so the download is
resumable and every byte is traceable to the request that produced it.
"""
import sys, time, urllib.request, urllib.error
from pathlib import Path

OUT = Path("/private/tmp/claude-501/-Users-ryanschmitt/07bd9c77-6066-4985-9818-751ff0c31e9c/scratchpad/delb")
ENDPOINT = "https://bass2000.obspm.fr/php/getSolarSpectrumDB.php?WL={start}&DW={dw}&resol=0.002&fmt=txt"
START, STOP, DW = 3000, 10000, 250

def get(start, dw, tries=4):
    url = ENDPOINT.format(start=start, dw=dw)
    for a in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.read()
        except Exception as e:
            print(f"  retry {a+1}/{tries} @{start}: {e}", flush=True)
            time.sleep(5 * (a + 1))
    raise SystemExit(f"FAILED permanently at {start}")

for s in range(START, STOP, DW):
    f = OUT / f"chunk_{s:05d}.csv"
    if f.exists() and f.stat().st_size > 1000:
        print(f"skip {s} ({f.stat().st_size} B)", flush=True); continue
    b = get(s, DW)
    f.write_bytes(b)
    nl = b.count(b"\n")
    print(f"got {s}-{s+DW}: {len(b)} B, {nl} lines", flush=True)
print("DOWNLOAD COMPLETE", flush=True)
