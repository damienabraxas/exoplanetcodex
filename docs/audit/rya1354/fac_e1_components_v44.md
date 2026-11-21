# RYA-1354 — signed native E1 component audit, v44

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. This milestone decomposes the EXISTING sealed v41 CI/rate models. It runs no new SCF, CI, Hamiltonian, expansion, rates, lifetime or grade.

Continue from `2b44a4cf05e286fd0406a30450c409f041707575` on the RYA-1354 research branch and draft PR616; worktree `/Users/ryanschmitt/codex/rya1354`. Resumed by Claude at Ryan's request after the 2026-11-17 handoff. Prospective declaration: Linear `d2b59f6f-c379-43df-899b-9ce5a8288f86`, before any build or native run. Plan `config/atomic/rya1354/lifetime_plan.v44.json`.

## The question

v43 showed that similar Babushkin (B, length) rates across the ground and equal-reference potentials hide compensating energy and amplitude changes. It also showed that the full many-electron Coulomb/Babushkin (C/B) amplitude ratios (0.892 ground, 0.702 equal) are not explained by the isolated reference 3s→2p orbital pair. v44 asks where inside the native E1 amplitude the gauge and potential differences actually live.

## What FAC actually sums (pinned source)

`TRMultipole` adds `AngularZMix` coefficients × `MultipoleRadialFR(aw)`, with `aw = α(E_up − E_low)`. The radial integrals are cached on, and interpolated across, the 2–3-point photon grid that the preceding `TransitionTable` sets from all requested E1 channels.

`AngularZMix` builds the CSF loop (`structure.c`) as follows:
- Stored level mixing is multiplied by the `AngularZMixStates` CSF-pair angular coefficients.
- Terms are skipped when |mix|, |mix₁·mix₂| or |term| is below `ANGZCUT = 1e-8`.
- Kept terms are merged by orbital pair after a `qsort`, and the bra/ket is swapped when the lower level's Hamiltonian index is higher.
- Before any of that, stored mixing is already truncated at `|mix| < MIXCUT·|principal| = 1e-5·|principal|` and renormalised (`AddToLevels`).

The v41 input sets no AngZ, mix-cut, PrepAngular or max_nbs control, so the source defaults apply. The driver rejects any input that would set them.

## Probe and controls

`scripts/rya1354_fac_e1_components.c` adds the observability-only SFAC command `RawE1Components` to the exact v27 rate frontend that v41 used (`toolchains/fac-config-probe-v27/sfac.c`; RawTransition/ConfigDiagnostics unchanged), linked to the unchanged `libfac.a`. It is called right after the last E1 `RawTransition` and before the next `TransitionTable`, so it sees the native grid and radial cache. It writes four record types:
- **A** — the native `AngularZMix` list in native order, with cached radial, product and running sum.
- **B** — every CSF-pair angular term, including the ones the native cut drops. Those are flagged and never summed.
- **L/N** — both levels' stored basis, mixing, configuration labels, and stored vs Hamiltonian basis counts and norm.
- **G** — the photon grid.

It does **not** call `TRMultipole`, whose optional pair selector reads an uninitialised variable (v43 finding). The native amplitude for comparison is the preceding unchanged `RawTransition` row.

Before declaration, a disposable four-configuration O I smoke calculation (throwaway build, `~/scratch` on Sirius; fixture `tests/fixtures/rya1354_e1_components_smoke`) reproduced the native amplitude bit-for-bit in both gauges, through both the swapped and unswapped bra/ket branches.

Measured null for byte identity: the v36 CI directory vs its v41 regeneration differ only in `input.sf`, `native.time`, `o.lev` and `o.lev.b` (session stamp and timing). `rates.tr(.b)` carry the same stamp and are checked by the unchanged independent rate audit. Every other file, including `sfac.stdout`/`sfac.stderr` and the 4M-line `o.basis`, must be byte-identical.

## Run

```sh
# Sirius, /mnt/codex-data/codex/rya1354
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 nice -n 19 \
  /usr/bin/time -v -o evidence/fac_e1_components_driver_v44/batch.time \
  python3 evidence/fac_e1_components_driver_v44/rya1354_fac_e1_components.py \
  --root . --source toolchains/fac-9ab6def \
  --index evidence/fac_e1_components_driver_v44/index.json \
  --index-sha256 a6bf301ccf6026e9d8fad77c070afd4c543089c2862c8de42fe7787469394644 \
  --probe toolchains/fac-e1-components-probe-v44 --out runs/fac_e1_components_v44
```

Index preparation took 18.82 s / 373180 KiB. The batch exited 0 in 4:02.10 at 1065352 KiB, with no swaps.
- The canary (new frontend on the unchanged v41 input) reproduced all 56 stable v41 files byte-for-byte.
- In all 8 cases (4 models × 2 gauges), the CI replay matches, all 56 stable files including `raw.tsv` are byte-identical to v41, and the independent rate audit equals v41.
- All 40 E1 contexts (5 channels per case) pass every exact identity:
  - the A-order running sum equals the native amplitude exactly;
  - products, cut flags and kept-term radial cache values are exact;
  - per-orbital-pair CSF sums fall within the recursive-summation bound γₙΣ|x| with u = 2⁻⁵³ (Higham eq. 4.4), since only the merge order differs. This is a floating-point bound, not a physical tolerance.

## Evidence (new lossless reference format)

`data/audit/rya1354/fac_e1_components_complete_v44/` holds 882 native paths:
- **138** new files are stored in 17 parts (12702050 compressed bytes, every part ≤5 MiB).
- **744** files are byte-identical to files already committed in `fac_radial_rates_complete_v41` (3.46 GB raw). They are recorded as **references**, not stored a second time.

Schema `rya1354.native-evidence-parts.v3` is v2 plus `referenced_archives` (pinned by manifest SHA-256) and `referenced_files`. The same decoder (`rya1354_native_parts.load_parts`, with `repo_root=`) decodes the earlier archive, restores each referenced path, and checks its hash before verifying every receipt. v1/v2 behaviour is unchanged.

Packing plus independent decoding took 17.11 s at 3.5 GB RSS on the Mac; log in `source_curation_v2/fac_e1_components_packing_v44.time.txt`. The machine report `source_curation_v2/fac_e1_components_v44.json` is regenerated from the decoded archives by `scripts/rya1354_fac_e1_components_report.py`, and a test requires exact equality.

## Findings

**Numerical side-channels are negligible for this line.** For the reference 3S1→3P2 line, amplitudes on the native photon grid equal the v43 exact own-energy radials to ≤5.0e-16 relative. 3P2 sits on, or one ulp from, a grid point. The 3P1/3P0 channels lie about 2–3e-3 (relative) inside the grid and are interpolated. 7154–7754 angular-cut terms drop at most 3.8e-6 of the 3P2 amplitude. Coherent10 and coherent100 give identical amplitudes, consistent with v41.

**What is not quantified post hoc.** The mixing cut leaves 577/370 (ground) and 612/395 (equal) Hamiltonian components out of the stored lower/upper eigenvectors, and FAC renormalises after dropping them. Their weight is not recoverable from the stored levels, so it remains an unquantified numerical axis.

**Dominance is not configuration probability.** Reference 3P2 amplitude by configuration pair (lower ← upper), coherent100:

| Model / gauge | Total | 2p4 ← 2p3 3s | 2p3 3p ← 2p3 3s | Lower 2p3 3p weight |
| --- | ---: | ---: | ---: | ---: |
| ground B | −1.40464 | −1.43748 | +0.07674 | 0.44% |
| ground C | −1.25313 | −1.21814 | −0.01029 | 0.44% |
| equal B | −2.02034 | −1.09804 | **−0.91508** | 1.01% |
| equal C | −1.41811 | −1.61544 | +0.16146 | 1.01% |

In the equal-reference model, a **1.0%** 2p3 3p admixture in the lower 2p4 3P2 level carries **45%** of the B amplitude. It reaches the upper 2p3 3s level through the one-electron 3s↔3p radial integral.

That integral is **−8.80 au** in length form but **+1.55** in velocity form at the line energy. Their ratio, 0.1765, equals the 3s–3p orbital gap divided by the photon energy (0.0519/0.2938) — the orbital-gap scaling v43 isolated. Ground gives −9.05 vs +1.21, ratio 0.134 = 0.0503/0.3751.

The admixture's angular-times-mixing coefficient is +0.0534 in the equal model and −0.0032 in the ground model. That is 17× larger and of opposite sign, far more than the √(weight) ratio of about 1.5, so the ground model's CSF terms in this channel nearly cancel.

The configuration-pair shifts account for each total change exactly:
- **ground→equal:** B −0.6157, made of −0.9918 from the 3p admixture and +0.3394 from the 2p4←3s principal channel. C −0.1650 = −0.3973 + 0.1718.
- **equal, B→C:** +0.6022, made of +1.0765 from the admixture and −0.5174 from the principal channel.
- **ground, B→C:** +0.1515, mostly the principal channel (+0.2193).

The same pattern holds for the 3P1 and 3P0 branches. The 1D2 and 1S0 intercombination channels have cancellation indices of 190–1200, so their tiny amplitudes are extremely fragile. They carry negligible rate.

**Interpretation (diagnostic, not a calibration).** The residual gauge split is set by a weakly weighted lower-level 2p→3p single-excitation admixture whose length and velocity radials differ by ω/Δε ≈ 5.7–7.5. Gauge agreement needs the correlation space to represent this channel consistently, not just more weight. In jj coupling the "principal CSF" share (18–38%) is not meaningful because the LS 3P2 is spread over several jj CSFs. The configuration pair is the meaningful unit. Nothing here selects a gauge, rescales energies or validates either model.

## Validation

```sh
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q \
  tests/test_rya1354*.py tests/test_atomic_benchmark_rya1354.py \
  tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py \
  tests/test_line_match_rya1033.py
```

New tests:
- 12 in `tests/test_rya1354_fac_e1_components.py`: insertion context, measured-null file set, smoke-fixture reconstruction in both branches, tamper rejection, exact-amplitude equality, parser structure, the summation bound, and the actual-archive replay with machine-report equality and the mechanism.
- 2 new v3 tests in `tests/test_rya1354_pack_native_stream.py`: reference round-trip, and rejection of a wrong reference or a changed earlier archive.

Final result: **260 passed in 234.92 s** (the previous 246 plus 14 new), 4984946688-byte max RSS, Mac Python 3.12.13 (source/metadata/report checks only).

Original FAC source diff is empty; the original library and executable hashes are unchanged; no owned native process remains. The throwaway compile-check directory `~/scratch/rya1354_v44_compilecheck` on Sirius is outside the repository and the run tree. Canonical scientific surfaces (grades, products, linelists, reference, catalog, stars, source spectra) are unchanged. No publication, merge or Done.

## Next bounded investigation (proposal — needs its own declaration)

1. **Numerical control first:** repeat one v41 case with `SetMixCut(0,…)` (separately declared) to bound the stored-mixing truncation that v44 cannot recover post hoc.
2. **Physical axis indicated by v44:** expand the lower-level 2p→n′p single-excitation and pair-correlation space, with p correlation orbitals beyond the spectroscopic 3p/4p. Then test whether the 2p3 n′p ← 2p3 3s channel converges and gauge-balances, with both potentials kept.

The v44 decomposition is the comparison instrument for each step. Core or radial expansion is not indicated as the leading axis by this evidence. Primary experimental provenance/covariance review and independent physical-convergence review remain owed before any theoretical grade or adoption.
