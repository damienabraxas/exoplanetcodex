# RYA-546 — Model-atom H-collision vintage audit + production-TS absolute-offset check

**Date:** 2026-07-12 · Branch `ryandamienschmitt/rya-546-...` (off origin/main 73088bc). **Audit/recommendation ticket — NO MERGE, no production change.**
Diagnostic script: `scripts/rya546_ts_offset.py` (Sirius, iSpec-TS). Sources: registry `config/constants.py` NLTE grids + the per-element NLTE papers + Gerber 2023 (A&A 669 A43).

## The pattern Ti found (RYA-542/544)

The older **MAFAGS-OS Bergemann-group** model atoms compute inelastic H collisions with the **Drawin
(1969)** formula (scaled by an SH factor, or set to SH=0 = no H collisions). Drawin under-thermalizes
trace neutrals → over-estimates the UV over-ionization → **inflates the positive NLTE abundance
correction**. Modern **ab-initio** rates (Barklem; Grumer & Barklem 2020; Amarsi/GALAH; Mallinson-2024)
thermalize the neutral → **smaller, ionization-balance-validated** corrections. Ti I: Bergemann-2011
Drawin gave +0.108 (Engine-A) / +0.20 (Engine-B); ab-initio Mallinson-2024 gives ~+0.05 (RYA-544).

## Part A — vintage audit of the registered Engine-A NLTE corrections

| element | Engine-A production grid (`nlte_grids/`) | H-collision recipe | solar δ applied | vintage verdict | modern ab-initio grid | recommendation |
|---|---|---|---|---|---|---|
| **Ti I** | `Ti_Bergemann2011_MPIA.csv` | Drawin-scaled (Bergemann 2011) | +0.108 | **CONFIRMED too large** (ab-initio ~+0.05, RYA-544) | **Mallinson-2024, Zenodo 10753497** ✓ acquired | in progress RYA-544/545 |
| **Mn I** | `Mn_Bergemann_MPIA.csv` | Drawin-scaled (Bergemann & Gehren 2008) | +0.11–0.17 | **SUSPECT — too large** | GALAH/Amarsi-2020 Mn + Gerber `atom.mn281kbc` + Grumer-Barklem-2020 Mn+H | **HIGH — re-derive ab-initio** |
| **Cr I** | `Cr_Bergemann2010_MPIA.csv` | Drawin, **SH=0** (no H collisions → maximal NLTE; Bergemann & Cescutti 2010) | +0.05–0.10 | **SUSPECT — likely spurious** (modern view: Cr I ≈ LTE) | none (not in Gerber-2023 / GALAH-13) | **re-examine: drop→LTE or acquire ab-initio Cr+H** |
| Ca I | `Ca_Mashonkina2017.csv` | Mashonkina 2017 (recipe unconfirmed) | small (~+0.02) | low impact | Gerber `atom.ca105b` | note; low priority |
| Na/Mg/Si/Al/S/K/N | `*_Amarsi2020/2025_PySME.csv` | **ab-initio (Barklem)** | (various) | **OK — modern** | (in use) | no action |
| C/O | `nlte_cno` (Amarsi 2019) | **ab-initio** | (various) | **OK — modern** | (in use) | no action |
| Ba II | `Ba_Korotin2015.csv` | Korotin 2015 | (ion) | ion — H-collision over-ionization mechanism N/A | Gerber `atom.ba111` | low priority |
| Sr II | `Sr_Bergemann2012_INSPECT.csv` | Bergemann 2012 | (ion) | ion — low H-collision sensitivity | Gerber `atom.sr191` | low priority |
| Co / Ni | **LTE** (no Engine-A NLTE) | — | 0 | no scaled-Drawin issue (just missing small NLTE) | Gerber `atom.co247qm`/`ni538qm` (**ab-initio "qm"**) | optionally ADD ab-initio NLTE (improvement) |

**The class Ti found:** the three MAFAGS-OS Bergemann-group **neutrals — Ti I, Mn I, Cr I — all use
the Drawin recipe and all inflate the over-ionization correction.** The Amarsi-PySME elements
(Na/Mg/Si/Al/S/K/N + C/O) are already ab-initio and fine. The ions (Ba II/Sr II) are insensitive.
Co/Ni are LTE in production (no vintage problem; the Gerber "qm" atoms would add correct small NLTE).

**Corroborating evidence for Mn (the same Ti pattern):** the RYA-534 Engine-B TS-Gerber Mn gate gave
δ **+0.043** — about HALF the Engine-A Drawin +0.107 — exactly as ab-initio should relative to
scaled-Drawin. This is the independent-engine cross-check (RYA-542 tool) already pointing at Mn. It
also explains the long-standing RYA-411 "Mn ~½" question: scaled-Drawin (Engine-A) vs ab-initio.

## Part B — production TS-synth-EW absolute-scale offset (RYA-545 addendum)

RYA-545: the production Turbospectrum synth-EW path (RYA-285 `_bisect_synth_abundance`) ran Ti ~+0.3
dex high in ABSOLUTE terms vs Route 1 (PySME, on solar). Two checks:

**B.1 — flat or line-dependent?** Per-line TS−PySME offset on the shared Ti lab lines:
+0.147/+0.325/+0.199/+0.228/+0.171 → **mean +0.21, stdev 0.06** → dominantly FLAT, ~0.06 line-dependent.

**B.2 — does it ride along on the audit elements?** Same production TS path, non-K10 lab-gf weak
solar lines, absolute A(X) vs Asplund-2021 (`ts_offset.log`):

| element | production TS A(X) | Asplund | offset |
|---|---|---|---|
| Ti I | 5.160 | 4.97 | **+0.19** |
| Cr I | 5.960 | 5.62 | **+0.34** |
| Ni I | 6.461 | 6.20 | **+0.26** |
| Mn I | — (HFS-limited: needs HFS synthesis, RYA-473) | 5.42 | — |

(With the K10 lines INCLUDED the raw offsets blow up — Ti +0.93, Cr +0.56 — because the K10 gf-scale
adds on top; excluding K10 collapses them to the ~+0.2–0.3 band. The K10 gf-scale is the RYA-521/545
problem, orthogonal to this.)

**Verdict:** the production TS-synth-EW absolute offset is **~+0.26 ± 0.08 dex, roughly CONSTANT across
elements** (Ti/Cr/Ni) and dominantly flat per-line. It is a **flat zero-point of the ew_floor
inversion path** (likely the blend-floor attribution + GES-synthesis normalization vs PySME), and it
**CANCELS in [X/H] = star − sun measured through the same code** (differential). The residual line-
(~0.06) and element- (~0.08) variation adds minor scatter, not a bias. **The RYA-545 balance
corroboration is untouched** (it is internal/differential). CAVEAT: **absolute** solar A(X) taken
from the production TS-synth-EW path runs ~+0.25 high — relevant only if an absolute (not
differential) A(X) is ever read off this path.

## Recommendations

1. **Mn I (HIGH)** — re-derive the NLTE correction via the ab-initio path (RYA-544 template): acquire
   the GALAH/Amarsi-2020 Mn departure grid (PySME `.grd`) OR adopt the Gerber `atom.mn281kbc` Engine-B
   value (+0.043). Mn is a registered/gold PASS element, so the ~2× inflation matters. This closes
   RYA-411's "Mn ~½".
2. **Cr I** — re-examine: the SH=0 Bergemann-2010 correction (+0.05–0.10) is likely too large; the
   modern view is Cr I ≈ LTE. Decide **drop Cr toward LTE** vs acquire an ab-initio Cr+H grid (none in
   Gerber-2023/GALAH — would need Barklem Cr+H if it exists).
3. **Ti I** — in progress (RYA-544/545; Mallinson ab-initio ~+0.05, grid banked).
4. **Co/Ni** — optionally adopt the Gerber ab-initio ("qm") NLTE (adds the correct small NLTE; today LTE).
5. **Flag** — the production TS-synth-EW ABSOLUTE scale is ~+0.25 high (flat; cancels in [X/H]). No
   action for differential abundances; note it if absolute A(X) is ever read from this path.

## Sources
- Ti: Bergemann 2011 (MNRAS 413 2184, Drawin); Mallinson 2022/2024 (ab-initio Grumer-Barklem 2020).
- Mn: Bergemann & Gehren 2008 (A&A 492 823, Drawin); Grumer & Barklem 2020 (A&A 637 A68, ab-initio Mn+H); GALAH/Amarsi 2020 (A&A 642 A62).
- Cr: Bergemann & Cescutti 2010 (A&A 522 A9, Drawin SH=0).
- Gerber et al. 2023 (A&A 669 A43) — Engine-B TS-native atom set.
- Amarsi et al. 2020 (A&A 642 A62) / 2025 (A&A 703 A35) — ab-initio PySME grids (Na/Mg/Si/Al/K/S/N).
