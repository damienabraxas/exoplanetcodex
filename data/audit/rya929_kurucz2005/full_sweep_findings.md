# RYA-929 full-sweep diagnostic

Comparison of the 1984 Kitt Peak Flux Atlas, the Kurucz (2005, revised 2010)
irradiance file, and the IAG/Baker 2020 telluric-free solar atlas. The sweep is
diagnostic only: it does not replace the EW or synthesis products and does not
produce an abundance.

## Window-level result

`fraction_below_0p5` is the fraction of sampled pixels below half of the local
continuum. The Kurucz and IAG values agree closely in the strong telluric
regions, while the original Kitt Peak atlas retains substantial absorption.

| Region (Å) | Kitt Peak 1984 | Kurucz 2005 | IAG/Baker 2020 |
|---|---:|---:|---:|
| O2 B 6867–6884 | 0.1493 | 0.0000 | 0.0000 |
| H2O 7160–7340 | 0.0287 | 0.0024 | 0.0018 |
| O2 A 7594–7685 | 0.3354 | 0.0029 | 0.0018 |
| H2O 8100–8400 | 0.0353 | 0.0032 | 0.0029 |
| H2O 9280–9600 | 0.2927 | 0.0003 | 0.0000 |
| Clean control 6600–6650 | 0.0014 | 0.0014 | 0.0015 |
| Clean control 6690–6702 | 0.0000 | 0.0000 | 0.0000 |
| Clean control 7690–7705 | 0.0163 | 0.0101 | 0.0101 |

This proves the broad correction behavior: Kitt Peak 1984 is not telluric-free;
Kurucz 2005 behaves like the corrected IAG reference at window scale. The
continuum normalization is a local diagnostic normalization, not a production
continuum.

## Line-level diagnostic

Depth is one minus the locally normalized flux minimum. These values are useful
for identifying failure modes, not for selecting a final abundance line.

| Line | Kitt Peak 1984 | Kurucz 2005 | IAG/Baker 2020 | Interpretation |
|---|---:|---:|---:|---|
| Al I 6631.218 | 0.0269 | 0.0142 | 0.0311 | broadly consistent clean control |
| Al I 6696.185 | 0.2544 | 0.0081 | 0.2572 | Kurucz line-level mismatch; investigate before use |
| Fe I 6872.162 | 1.0010 | 0.0286 | 0.0980 | Kitt Peak saturated by O2-B contamination |
| Fe I 6875.445 | 1.0009 | 0.0127 | 0.0502 | Kitt Peak saturated by O2-B contamination |
| Si I 6876.359 | 1.0011 | 0.0058 | 0.0756 | Kitt Peak saturated; line assignment/continuum needs review |
| Fe I 6881.442 | 0.2137 | 0.0542 | 0.2016 | corrected references agree more closely than Kitt Peak |
| K I 7665.000 | 1.0015 | 0.0263 | 0.6425 | Kurucz 7665 mismatch; do not accept as a corrected line result |
| K I 7698.964 | 0.8078 | 0.0392 | 0.8074 | Kurucz mismatch despite agreement of Kitt Peak/IAG |
| N I 8216.000 | 0.0685 | 0.1142 | 0.0651 | Kurucz line-level discrepancy in H2O region |

The broad-window result is therefore a GO for framework calibration, not a GO
for blindly substituting Kurucz line measurements. The K I 7665/7699 behavior,
Al I 6696, and N I 8216 require wavelength-medium, sampling, and line-window
forensics before any EW product can use Kurucz 2005. The 1984 Kitt Peak values
remain useful as the raw/uncorrected comparison arm.

Artifacts:

- `full_sweep_baker.csv`
- `line_sweep_baker.csv`
- `../diagnostic_line_plan.csv`
- `../../../../scripts/rya929_full_sweep.py`
