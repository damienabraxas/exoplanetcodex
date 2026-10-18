# Solar Mg: Asplund, Amarsi & Grevesse (2021)

RYA-1250. Source: A&A 653, A141, DOI 10.1051/0004-6361/202140445,
local `aa40445-21.pdf`, Table A.1 on printed page 27, visually checked.
The magnesium discussion is **Sect. 2.4, page 9**, not Sect. 3.

`mg_agss21_lines.csv` transcribes all eight Mg I and six Mg II rows, including
the authors' log gf, lower/upper energies, equivalent widths and 3D-NLTE
abundances. Wavelengths are the printed air wavelengths in nm multiplied by
10; EW remains in pm. No missing significant figures have been invented.

Existing canonical IDs were checked on species, wavelength (within 0.1 A of
the low-precision published wavelength), and lower energy (within 0.001 eV).
These are explicit joins, not new automatic nearest-wavelength matches.
The generated set keeps published and canonical wavelengths separately.
Mg I 9983.20 and Mg II 10092.20 have no canonical row in that window.
The 8712/8717 and 15879 features have fine-structure components; a canonical
join identifies the existing feature but does not validate its component gf
allocation. That must be resolved before adoption.

## Application update (2026-10-18)

Table A.1's full eight model columns are transcribed. The shared solar-table
registry declares cited MARCS-NLTE minus MARCS-LTE and 3D-NLTE minus MARCS-LTE
corrections. These are additive corrections, not our own 3D transfer calculation.

The saved official NIST response `nist_mgii_asd_20261018.csv` verifies A+ for
7877/7896/9218/9244 and A for 10914 using species, endpoint energies, wavelength
and oscillator-strength consistency. Keep AGSS21 gf exactly; the uncertainty
is the existing NIST-grade conversion and is not an AGSS21 quoted sigma.
Five Mg II rows are now adopted and graded. The missing 10092
feature and all eight Mg I published-set uncertainties remain unresolved.
The conservative source-error model treats these five NIST errors as correlated.

The six saturated Mg I candidates are retained for synthesis, not culled.
The Mg I production correction table now contains 66 rows after extension.
The following sections retain the initial investigation and its earlier holds.

## Oscillator strengths and uncertainty holds

Keep the exact Table A.1 gf. It supplies no per-line gf uncertainties or
per-row source labels. Sect. 2.4 cites Pehlivan Rhodin et al. (2017) for Mg I;
it describes all six Mg II gf as NIST A/A+ with errors <= 0.01 dex, but does
not identify which grade belongs to which line. The CSV retains that bound
as a note, not an invented per-line standard deviation.

Pehlivan Rhodin et al., A&A 598 A102, DOI 10.1051/0004-6361/201629849,
author manuscript https://arxiv.org/abs/1611.07000, Tables 2 and 3 inspected:
theoretical Table 3 gives -2.020, -2.242, -1.679, -2.177 and -1.718 for the
6318, 6319, 8923, 9983 and 10312 features, matching AGSS21. Table 3 has no
per-line uncertainty column. Summing its 8712 and 8717 component gf gives
the AGSS21 feature gf to printed precision. Do not borrow the experimental
Table 2 errors for these theoretical values. At 15879, Table 2 gives
-1.326 +/- 0.06 for the 15879.567 component while Table 3 gives -1.194;
neither is the AGSS21 feature value -1.226. Its adopted feature gf and error
need clarification; retain -1.226 exactly and leave sigma unknown.

Consequently the first build has **8/6 full rows and 0/0 graded rows**.
The prepare code only proposes gf adoption when both gf and sigma exist.
Registering this set alone therefore does not put the authors' gf into the
measurement. This is an explicit hold, not a completed reproduction.

## Combination recipe

Sect. 2.4: first take an unweighted mean within each ion. Mg I is
7.56 +/- 0.01 statistical +/- 0.02 systematic; Mg II is
7.52 +/- 0.03 statistical +/- 0.03 systematic. Combine ion means using their
total statistical and systematic errors. Interpreting those weights as
inverse total variance gives weights proportional to 1/0.0005 and 1/0.0018,
and 7.5513 from the rounded ion means. This reproduces the published 7.55;
the published total uncertainty is 0.03. A website implementation should
use measured ion budgets, preserve both ion results, and must not label an
LTE-only Mg II result as 3D-NLTE or silently substitute Mg I alone.

## Culls awaiting Ryan

None of the six initial candidates is in Table A.1: 4571.096, 5167.322,
5172.684, 5183.604, 5528.405 and 5711.088 A. However, Lodders et al. (2025),
local `AbundancesLodders.pdf`, pages 10-12, explicitly discusses 5711 A as
the clean solar Mg feature. Absence from the registered AGSS21 set is not
proof that no solar analysis uses it. No culls have been applied.

## Sirius solar-node NLTE check

The dedicated PySME interpreter completed `--solar-only`. Existing corrections
reproduce within 0.0002 dex: 4730.029 is -0.0088 (banked -0.0089), and
5711.088 is -0.0345 (banked -0.0343). New solar-node corrections are
4571.096 -0.0166, 5167.322 +0.0024, 5172.684 +0.0025, 5183.604 +0.0032,
and 5528.405 -0.0326 dex. The full-node `--write` is deferred at the
checkpoint; these five lines are also proposed culls. Mg II remains absent
from the registered 22-row table. Log: `data/output/rya1250/nlte_solar_check_pysme.log`.

## Current scope and Mg II readiness

Ryan confirmed both Mg I and Mg II, with Asplund as the guiding reference.
The earlier Mg-I-only interpretation is superseded. The staged Amarsi/PySME
1D departure grid contains levels for all six Mg II features, and shared ion-aware adapter/extension support is now implemented in the
unmerged branch. Five lines pass departure-engagement smoke checks; only two
yield finite diagnostic corrections. The Gerber atom has no Mg II
bound-bound transitions. Full 3D NLTE equivalence is not established.
See `docs/audit/rya1250_mgii_nlte_readiness.md`. No gf adoptions, culls, production correction-table writes or publication were applied.

## Neutral NIST follow-up

The 2026-10-18 Mg I ASD query is saved under `data/output/rya1250/`. It
does not resolve the adopted PR2017 uncertainties. In particular, 8923 has a
NIST B-rated oscillator strength numerically close to AGSS21, but its T5539
reference is Chang & Tang (1990), DOI 10.1016/0022-4073(90)90053-9, not
Pehlivan Rhodin et al. (2017). That grade is not transferred to the adopted
PR2017 value merely because the two calculations agree. The other features
include different values and unresolved component combinations.
