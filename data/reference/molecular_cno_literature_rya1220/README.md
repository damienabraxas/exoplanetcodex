# Molecular CNO literature reference (RYA-1220 Work Package A)

How molecular C/N/O has actually been measured, what those authors got, and which
transitions they used. Assembled from holdings already in this repository: the ingested
Amarsi et al. 2021 Table 2 (CDS `J/A+A/656/A113`, 408 used lines) and the RYA-1220
evidence matrix. Nothing here is a new literature claim; it is the published record
put next to our own diagnostics.

## 1. The modern 3D reference — Amarsi et al. 2021

Homogeneous 3D LTE solar analysis, 408 lines. Per-line abundances are published under
five model atmospheres, which is what makes the model-form term measurable rather than
assumed.

| molecule | system | n | λ_vac (nm) | A(3D) | A(MARCS 1D) | 3D−MARCS | 3D−⟨3D⟩ | line-to-line sd |
|---|---|---:|---|---:|---:|---:|---:|---:|
| C₂ | Swan | 39 | 472.8–562.4 | 8.454 | 8.431 | **+0.023** | −0.018 | 0.028 |
| CH | X-X | 54 | 425.4–3795.5 | 8.471 | 8.436 | **+0.035** | −0.028 | 0.023 |
| ¹²C¹⁶O | X-X | 80 | 2295.4–6329.2 | 8.470 | 8.610 | **−0.140** | −0.211 | 0.019 |
| CN | A-X | 59 | 1087.5–1320.8 | 7.864 | 7.925 | **−0.061** | −0.098 | 0.025 |
| NH | X-X | 31 | 2891.0–15023.0 | 7.908 | 7.941 | **−0.033** | −0.093 | 0.032 |
| OH | X-X | 145 | 1528.3–12280.3 | 8.700 | 8.765 | **−0.065** | −0.093 | 0.030 |

Per element: C 8.468 (3D−MARCS −0.009, n=173), N 7.876 (−0.052, n=90), O 8.700 (−0.065, n=145).

Two things follow directly:

- **The 1D→3D correction is molecule-specific, not a single CNO offset.** It spans +0.035
  (CH) to −0.140 (CO). Any "molecular 3D term" applied uniformly is wrong by up to 0.18 dex.
- **Line-to-line scatter is small (0.019–0.032 dex) in every molecule.** Our molecular
  routes run red_chi2 33–213. That gap is ours, not the method's.

🔴 **NH UV electronic lines are NOT retained** by this analysis (its own stated limitation).
Their NH is vibration-rotation in the infrared, 2.9–15 µm. The near-UV NH A-X band is
absent from the modern reference set by choice.

## 2. What the rest of the field does, by regime

| study | regime | transitions used | what they got / the limit they hit |
|---|---|---|---|
| Amarsi 2020 | Sun | N I 7442.29, 8216.33, 8629.23, 8683.40, 10108.90 Å | 0.05 dex systematic vs 0.002 dex mean SE; 4 of the lines are CN-blended and corrected empirically |
| Amarsi 2021 | Sun | 59 CN A-X (0-0); 13 NH rotational + 18 NH rovibrational | see table above; temperature sensitivity is the stated failure mode |
| AGSS21 | Sun | the 5 N I + the molecular set | **atomic 7.77 vs molecular 7.90; recommended 7.83 ± 0.07** |
| Mashonkina 2024 | Sun + A-F | same five N I, plus stellar N I | solar 7.92 ± 0.03 scatter; 0.12/0.05 dex continuum shifts |
| Botelho 2020 | solar twins | CN B-X 4180.02, 4192.94, 4193.40, 4195.95, 4212.25 Å | differential precision only; does not establish absolute gf |
| Smith 2013 | K/M giants | H-band CN A-X with CO/OH, ¹³CO | coupled CNO iteration, CO then OH; dredge-up moves photospheric C/N |
| Suárez-Andrés 2016 | planet-host FGK | NH 3345–3375 Å | continuum allowance 0.1 dex; **solar-calibrated (tuned) gf** |
| Spite 2005 | EMP giants | NH 336 nm; CN 388.8 nm | **NH/CN offset 0.4 dex**, patched by an empirical correction |
| Spite 2022 | metal-poor dwarfs | NH A-X near 336 nm | **line-list choice alone moves NH by 0.44–0.52 dex over 3357–3365 Å** |

The atomic-vs-molecular N tension is **0.13 dex in the published literature itself**
(7.77 vs 7.90). It is not an artifact we introduced, and no reconciliation should be
expected to close below it.

## 3. Our diagnostics against that record

| our key | molecule | our window (air Å) | reference used-range (Å) | in the reference set? | our red_chi2 |
|---|---|---|---|---|---:|
| `CN_AX_J` | CN | 11640–12076 | 10875–13208 | **YES** | 56 |
| `CN_AX_IR` | CN | 10871–11801 | 10875–13208 | **YES** | 195–213 (KP) |
| `CN_red` | CN | 6125–6200 | 10875–13208 | **NO** | 33 |
| `NH_AX` | NH | 3358–3373 | 28910–150230 | **NO** | 37–44 |

🔴 **Two of our four molecular N diagnostics sit on bands the modern reference analysis
does not use**, and they are the two carrying our worst pathologies:

- `CN_red` (6125–6200 Å) is the diagnostic with the **2.404 dex** target−Sun difference.
  The CN red system is a different vibrational band from the A-X (0-0) set at 1.09–1.32 µm.
- `NH_AX` (3358–3373 Å) is exactly the band Amarsi 2021 declined to retain and where
  Spite 2022 measured a **0.44–0.52 dex** spread from line-list choice alone. Our own
  near-UV opacity deficit (RYA-1189/1190/1204/1207) sits on top of that.

This separates two different failures that the uniform HOLD had merged:

- **wrong band** — `CN_red`, `NH_AX`. No amount of covariance work rescues these.
- **right band, bad fit** — `CN_AX_IR` on Kitt Peak (red_chi2 195–213) while the same
  band on IAG and CRIRES+ J is the one that already gives A(N) ≈ 8.00.

## Provenance

- `amarsi2021_molecular_reference_by_species.csv` — derived from the ingested Table 2.
- `our_diagnostics_vs_amarsi2021.csv` — the cross-match above.
- Source rows: `data/reference/amarsi2021_cno/`, `data/reference/nitrogen_method_rya1220/evidence_matrix.csv`.

Per-line rotational identities are **not** in Table 2, so those rows stay
`CROSSMATCH_REVIEW` until joined to the upstream molecular releases on physical identity.
A wavelength match alone is not a transition match.
