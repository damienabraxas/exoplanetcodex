# RYA-1395: O I s0 (spectroscopic, Dirac-Fock) stage complete (v72 even + v75 odd)

Status: **HOLD, report only.** This stage has no correlation, and nothing here is adopted or graded.

## What this stage is
- Our own GRASP MCDHF-EOL model of O I, following the Li et al. 2023 scheme.
- Orbitals are added one per step: half grid, strict SCF criterion (ACCY 6e-8).
- Every step must pass three gates:
  - identity (the leading LSJ families),
  - radial node count = n − l − 1,
  - orbital energy within 30% of the NIST binding.
- The final fixed-orbital CI list of each parity equals the v56 s0 list.

## Run history
| Version | Result | What we learned |
|---|---|---|
| v71 | FAILED_RETAINED | GRASP "screened hydrogenic" = bare Z = 8 for core-free lists; odd 3d went to the continuum |
| v72 | **even COMPLETE**; odd FAILED at 5s | Own hydrogenic starts (Z_eff = n/n*); 5s took a 6-node root under method 3 |
| v73 | odd FAILED at 5d | Method 1 (node counting) fixes 5s and 6s; 5d stagnates at self-consistency 5e-7 |
| v74 | odd FAILED at 5d | Method 1 on 5d gives the **identical** fixed point, so the residual is a numerical floor |
| v75 | **odd COMPLETE** | Declared measured-floor admission for 5d only (≤ 1e-6, 500 identical iterations) |

## Orbitals vs NIST binding
- p series and d series: within 1%.
- 3s −7%, 4s −4%, 5s −3%, 6s −2%.
- 4f and 5f: hydrogenic Z_eff = 1, −0.2%.

## E1 s0 log gf vs Li 2023 Table A2 (Babushkin gauge, NIST energies)
989 of 994 lines are matched. `d` = ours − Li:

| Li class | n | median d | 68% of abs(d) |
|---|---|---|---|
| A | 79 | −0.008 | 0.117 |
| B | 126 | +0.006 | 0.136 |
| C | 316 | −0.162 | 0.639 |
| D | 270 | −0.611 | 2.206 |
| E | 198 | −0.452 | 1.701 |

Key solar lines (d, dex):

| Line | d |
|---|---|
| 777 triplet | +0.029 |
| 844.6 | +0.005 |
| 926 | −0.097 |
| 1129 | +0.005 |
| 1316 | +0.007 |
| 615.8 | +0.095 |
| 645.6 | −0.145 |
| 130.2–130.6 | −0.050 |
| 135.6 (intercombination) | −0.14 / −0.10 |

## Known limits (they motivate the correlation stage)
- **Gauge disagreement is large.** Median abs(B − C) is 0.37 dex for class A lines and 1.39 dex for class B. The Coulomb gauge is driven by the model transition energy.
- **The cross-parity energy scale is off.** Even and odd are optimised separately, and the model puts 3d below 3p, which flips 309 line orientations in the report. Odd 3s ⁵S sits about 12,000 cm⁻¹ low against the ground state, and even 3p ⁵P about 6,600 cm⁻¹ high.
- **Class C/D/E lines** (weak, cancellation, intercombination) are not reliable at the s0 level.
- **Agreement with Li at class A/B is partly fortuitous.** Babushkin S is robust at the DF level, but a reproduction without correlation is not a validation.

## Files
- **Drivers:**
  - `scripts/rya1395_mcdhf_s0_v7[2-5].py`
  - `scripts/rya1395_hydrogenic_rwfn.py` (control: `scripts/rya1395_hydrogenic_control.py`)
  - `scripts/rya1395_s0_e1_report.py`
- **Data** (`data/audit/rya1395/mcdhf_s0_v75/`): `s0_loggf_vs_li2023.csv`, `s0_report.txt`, `lower.upper.t`, and the receipts.
- **Native runs** (Sirius, `/mnt/codex-data/codex/rya1354/runs/`):
  - `mcdhf_s0_v72_even`
  - `mcdhf_s0_v75_odd`
  - `mcdhf_s0_v75_E1`
  - kept failures: v71, v72_odd, v73_odd and v74_odd
