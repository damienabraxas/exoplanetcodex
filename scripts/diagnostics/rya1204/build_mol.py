"""RYA-1204 Lever B — write the 6,223 excluded near-UV VALD molecular lines as
Turbospectrum .bsyn lists.

🔴 THE ISOTOPOLOGUE, WHICH IS WHY THEY WERE EXCLUDED. VALD's extract says 'OH 1' and
names no isotopologue; TS's molecular row needs one. The near-UV list builder refused to
guess, correctly citing RYA-684. The assignment here is the DOMINANT isotopologue, and
the reason that is not the RYA-684 guess is that the error it can introduce is BOUNDED
AND TINY: 16O is 99.76% of solar O, 14N 99.64%, 12C 98.93%, 1H 99.99%. If VALD's gf is
the all-isotopologue total, assigning it to the dominant one over-counts that species by
<=0.25%/1.1%; if VALD's gf is already per-dominant-isotopologue, the assignment is exact.
Either way the error is at most ~1%, against a status quo of 100% -- the band currently
carries NO molecular opacity at all.

The isotopologue also cannot move a line: it sets the reduced mass, and these wavelengths
come from VALD, not from a computed term value. It enters only through the abundance.
"""
import re, numpy as np, collections
RAW='/Users/ryanschmitt/codex/rya1195/data/linelists/vald_solar_nearuv_2000_3780_hfson_raw.txt'
CODE={'OH':'0108.000016','NH':'0107.000014','CN':'0607.012014','CH':'0106.000012'}
LO,HI=3000.0,3780.0
rows=collections.defaultdict(list)
pat=re.compile(r"'(OH|NH|CN|CH) (\d+)',\s*([\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),\s*(-?[\d.]+),")
for l in open(RAW,errors='replace'):
    m=pat.match(l)
    if not m: continue
    sp,ion,wl,gf,elo,jlo,eup=m.group(1),int(m.group(2)),float(m.group(3)),float(m.group(4)),float(m.group(5)),float(m.group(6)),float(m.group(7))
    if not (LO<=wl<=HI): continue
    ju=float(l.split(',')[6])
    rows[sp].append((wl,elo,gf,2*ju+1))
tot=0
for sp,v in sorted(rows.items()):
    v.sort()
    fn=f'molecules/{sp}_300-378.bsyn'
    with open(fn,'w') as f:
        f.write(f"'         {CODE[sp]}'    1   {len(v)}\n")
        f.write(f"'VALD near-UV {sp} A-X (RYA-1204); dominant isotopologue'\n")
        for wl,elo,gf,gu in v:
            # wavelength, Elow(eV), loggf, fdamp, gu, gamma_rad -- VALD gives 99.000 for
            # unknown damping, so 0.0 is written and Turbospectrum's own default applies,
            # exactly as the atomic builder does for a missing gamma_rad.
            f.write(f"{wl:12.3f} {elo:9.5f} {gf:7.3f}    0.000 {gu:6.1f}  0.00E+00 'X' 'X'   0.0    0.0\n")
    print(f'  {fn}: {len(v):5d} lines  {v[0][0]:.1f}-{v[-1][0]:.1f} A  code {CODE[sp]}')
    tot+=len(v)
print(f'total {tot} molecular lines written (expected 6223)')
