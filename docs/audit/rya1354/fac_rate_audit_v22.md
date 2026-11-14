# Minimal FAC model: complete one-photon rate accounting

This audit resolves the two native E3 refusals retained at v21. Scientific HOLD and EXPOSED_NONBLIND remain. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. These are exploratory totals in the declared minimal model, not qualified atomic data or grades.

The independently requalified FAC wavefunctions contain exactly15 CSFs from 1s2 2s2 2p4 and 1s2 2s2 2p3 3s1. The auditor verifies every CSF occupation and parity, not only the orbital list. Between these opposite-parity configurations the one-electron transition is 3s -> 2p. Its largest one-electron j sum is 1/2+3/2=2, so tensor rank3 fails the angular triangle. Thus both total-J-allowed E3 channels are exact angular zeros within this restricted basis. Native AngularZMix returns nz0/nmk0, and TRMultipole's corresponding refusal branch returns-1. Raw NaNs and native statuses remain unchanged; the derived report records a separately justified model zero. This argument does not carry over to an expanded configuration space.

The remaining11 channels have finite native rates. Replaying pinned dbase.c OscillatorStrength, using archived consts.h values, reproduces raw gf/A within double roundoff. Each of the10 standard printed rows then matches exactly at native %.6E precision after independently rounding the amplitude to IEEE float32, as required by TR_RECORD.strength. Header electron count/gauge/mode, duplicate channels and block counts are checked. E3 produces no table block. The only nonzero omission is M3 to level5: A=9.463351238284966e-20/s at coherent100. Dividing by RATE_AU gives a rate below the explicit SaveTransition0 trd<1e-30 gate. It is retained in the total from the raw diagnostic; absent printed rows are never generally assumed zero.

Upper level6 has J1. All lower model levels are the five even states0..4 and odd quintet5, with J<=2. Total-J triangles therefore cap rank at3. Parity selects E1/M2/E3 to the even states and M1/E2/M3 to the quintet:13 allowed requests,11 finite rates and2 explicit model angular zeros. This is complete one-photon accounting for this15-level model, not a claim about omitted configurations, continuum effects, multiphoton processes or physical accuracy.

| Numerical setting | Babushkin lifetime (ns) | Coulomb lifetime (ns) |
|---|---:|---:|
| coherent10 | 2.4121056376677217 | 2.3169782257545464 |
| coherent100 | 2.4121030293305274 | 2.316978563137304 |

Numerical stability across these two controls is not a model uncertainty. Both gauges are retained; neither is selected against the exposed experimental lifetime. FAC central-potential, nuclear defaults, minimal correlation space and ground-only optimization differ from the existing correlated n5 GRASP model. Their lifetimes are not a matched-method accuracy test. Shared GRASP angular conversion remains declared.

Reproduce with `python scripts/rya1354_fac_rate_audit.py`. Inputs are the immutable fac_rates_complete_v21, fac_coherent_complete_v18 and fac_recoupling_complete_v19 archives. No new native calculation or source modification occurred. Machine report: source_curation_v2/fac_rate_audit_v22.json. Regression tests reject a changed printed rate, unexplained native refusal, changed basis and missing channel. Full native raw/formatted evidence, including failed v20, remains preserved.

Next: declare a bounded FAC correlation-space expansion independently of exposed reference values, preserve both numerical settings and gauges, re-establish state/LS identities and complete lower-state inventory at every new space, and assess model-space stability. The existing minimal model does not yet justify known-grade scoring, calibrated uncertainty, adoption, publication, merge or Done.

Validation: all111 tests selected by tests/test_rya1354*.py passed in93.26s. The final six rate-audit cases pass after strengthening the exact CSF occupation guard and adding the occupation-mutation case. Whitespace and unchanged protected-surface comparison against the task base pass. The bibliography audit still stops on existing civis2013/mallinson2024 verification-state values before library traversal; no library audit success is claimed. No new upstream defect was discovered: the E3 selection rule and M3 writer cutoff are source behavior, now explicitly accounted for.
