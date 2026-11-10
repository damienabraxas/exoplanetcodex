# CN source interface audit — RYA-1281 / RYA-1278

2026-11-08. Scientific HOLD. Worktree codex/rya1279, branch
ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, stacked base4635c742,
starting5899d41c. Predeclared plan/implementation906e2b51. Mac only, no RT.

## Source semantics established

Native Fortran pin21269948750e2177fe32c5df38e1d9bb81363374 and installed inputs
retain separate authorities. Eleven source/data/input receipts pass before and
after reading; model and exact original CN input also match. Source excerpts with
line numbers and fingerprints are in the report, not inferred from molecule names.

- bsyn.f:1281–1284 forms total CN molecules per gram from full-network partial
  pressure divided by kT and mass density, then obtains molecular partition Q.
- bsyn.f:1350–1351 multiplies one isotope fraction per constituent. The special
  symmetry factor at1250ff is for C2, not an automatic CN multiplicity.
- partffordepth.f:76 returns default-REAL Q from elemental molecule matching;
  no separate isotope-specific CN Q is selected by that interface.
- molecpartf.f:161–167 evaluates exp(Horner(coefficients,ln(TEMP))). TEMP/log
  expression is default REAL before assignment to double; coefficients and
  accumulator are double. File IRWIN_molecules_v15.1.dat:345 has CN code607.00,
  file energy -7.724 and six coefficients. No independent physical energy or
  temperature-validity qualification is inferred from the row.
- die_pe_lu.f:167 uses reciprocal molecular Q in the equilibrium calculation;
  changing Q only in the line opacity would not establish coherent chemistry.
- bsyn.f:949 sets g=1 and1235 computes f=10**loggf; depth.f:163 divides the
  Boltzmann-weighted population by Q, then165 applies stimulated-emission factor.
  This identifies conventions to match, not an independently validated coefficient.

The source polynomial projected onto56unchanged model temperatures4045.9–9909.1K
produces Q6430.37793–51498.5625. CSV retains each depth. This is a Python
source-semantics projection, **not a call to compiled native partition code**, not
an independent partition benchmark and not proof of fit validity across this range.
No full-network pressures, charge state, CN number densities or line extinction
coefficients were calculated/exported. No nuclear-spin factor added or removed.

## Isotope serialization result

Source fractions: C12=.9893,C13=.0107,N14=.99632,N15=.00368. The original native
input serializes .989/.011/.996/.004, exactly confirmed against the accepted input.
Products are evaluated once per constituent, without renormalizing or changing inputs.

| Species | Source fraction product | Serialized product | Relative change |
|---|---:|---:|---:|
|12C14N|0.985659376|0.985044|−0.062433%|
|13C14N|0.010660624|0.010956|+2.770720%|
|12C15N|0.003640624|0.003956|+8.662691%|
|13C15N|0.000039376|0.000044|+11.743194%|

The first three have delivered750–800nm routed files; no13C15Nfile is included
in the predeclared three-file CN route. These fractions do not determine the
spectral effect or line-list completeness error. Isotope masses and line positions
also matter. The existing spectrum remains unchanged; this audit neither repairs
its wrapper-compatible serialization nor claims a measured spectral rounding error.
Source labels such as Sneden web do not supply the missing typed primary release,
license, frame, partition or isotope/statistical provenance. Existing canonical
BLOCKED_PROVENANCE status remains authoritative. Legacy native LTE sensitivity
must stay distinct from a qualified molecular model in Bride NLTE transfer.

## Attempts, verification and reproducibility

First audit command failed before creating any result directory: it incorrectly
asked the Fortran Git tree to own installed SPECTRUM.lst. Error128from git show,
Python exit1, log /private/tmp/rya1281-cn-interface.log preserved. Fixfb590dd2
checks separately held isotope/list receipts against accepted original input and
inventory records; Fortran/DATA still byte-match pinned native Git. This is an
auditor provenance-contract defect, not a defect in the isotope source or native
physics. Run-bug-ledger row appended. Fresh outputv2, no overwritten results.

Command (from worktree, after source config/hosts/mac.env):

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_cn_interface.py --output <fresh-output>
```

Accepted registered output experiments.bride_cn_interface. Two results are copied
byte for byte to data/results/rya1281/cn_interface_v1 and generator-registered.
stdout /private/tmp/rya1281-cn-interface-v2.log, exit0. Source before/after checks
are additional runtime evidence, not extra independent tests.

Preflight2tests passed0.12s; final16tests passed0.19s, exit0:

```
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q tests/test_bride_cn_interface_rya1281.py tests/test_bride_native_lte_build_rya1281.py tests/test_bride_native_molecules_input_rya1281.py tests/test_molecular_canonical_rya1130.py
```

Tests cover ambiguous/malformed neutral-CN rows, natural-log Horner projection,
invalid temperatures, duplicate isotope rows and exactly-once serialized products,
plus existing native build/input and canonical intake contracts. They do not prove
native runtime partition parity or primary physical correctness.

Fetched current main3baabba9 (PR610Ca inputs/products and PR611completeness).
Current-main canonical_gf changed; this audit does not consume it or rerun atomic
spectra. Previously pulled RYA1315 ORCH continuum/gf-resolver/molecular owners are
unchanged in the inspected delta; retained earlier experiment stays pinnedde96ae51.
No claim that the old full canonical table equals current main. No merge/rebase.

Library audit196documents/127accounted/69unindexed/12off-Mac, exit1. No new
literature acquired, primary claims fabricated or bibliography entries inferred.
Existing CNO/Fe/refinement debts are retained; no phase-close or gate admission.

## Next bounded task

Prepare an isolated source-receipted **reporting** build that exports CN pressure,
Q, isotope fraction, mass and actual line coefficient at all56depths. Prove its
normal spectra unchanged against the accepted native binary before comparing the
exported coefficient with a separately implemented equation under identical state.
Carry EOS pressure/density, abundance and complete equilibrium network explicitly;
never approximate CN as undepleted C times N or add a fitted spin factor. Retain
source/table validity and primary line-list qualification as separate open gates.
Only then evaluate a coherent same-state Bride adapter and fresh O NLTE response.
No multiplication/addition of normalized spectra, abundance tuning, production,
Sirius, publication, adopted product change or merge is authorized by this audit.

Final checks: MIRROR/version/sequence/whitespace pass. New main integration was
recorded explicitly for merged RYA1311/1349 before sequence reconciliation.
Generator diagnostic set exactly2436historical, zero added/removed. Historical
RYA587blank science_impact and tracker/gold debts remain; no global-green claim.
