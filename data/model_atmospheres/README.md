# Model Atmospheres

This directory holds stellar atmosphere models used for abundance analysis.

## ATLAS9 Grid (Kurucz)

Model type: `ATLAS9` — 1D, plane-parallel, LTE, opacity-sampled

These files are large binary/text files and are **not included in the repo**.
Download instructions below.

### Acquiring the ATLAS9 Grid

**Option A: Kurucz website**
1. Visit [Kurucz ATLAS9 models](http://kurucz.harvard.edu/grids.html)
2. Download the grid matching 55 Cancri's parameters:
   - Teff range: 4750–5750 K (step 250 K)
   - log g range: 4.0–5.0 (step 0.5)
   - [Fe/H] range: 0.0–+0.5 (step 0.1)

**Option B: MOOG companion grid**
If using MOOG for radiative transfer, download the standard `odfnew` grid:
```
wget http://kurucz.harvard.edu/grids/gridm01odfnew/
```

**Option C: iSpec bundled grid**
iSpec ships with a pre-interpolated ATLAS9 grid.
See [iSpec documentation](https://www.blancocuaresma.com/s/iSpec).

## Model Parameters for 55 Cancri A

Interpolation target (from `config/constants.py` → `STAR_55CNC`):

| Parameter | Value | Source |
|-----------|-------|--------|
| Teff | 5196 K | von Braun et al. 2011 (CHARA) |
| log g | 4.41 | von Braun et al. 2011 |
| [Fe/H] | +0.32 | Initial estimate |
| vturb | 0.9 km/s | Initial estimate |

These are starting values; the final parameters are determined iteratively
through excitation/ionization equilibrium in `pipeline/05_abundances.py`.
