# RYA-1354 v53–v60: Li-scheme MCDHF-EOL, spectroscopic (s0) stage

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. No rate from these runs is a result.

**Goal.** Build Li et al. 2023's GRASP MCDHF-EOL orbitals for O I, layer by layer, starting with the Dirac-Fock spectroscopic stage (s0).
- **Even targets** (45 states): 2p4; 2p3(4S)np, n=3–6; 2p3(4S)nf, n=4,5.
- **Odd targets** (36 states): 2p3(4S)ns, n=3–6; 2p3(4S)nd, n=3–5; 2p3(2D)3s.
- **Orbital sets:** separate for even and odd parity, as in Li Sect. 2.3.

Every orbital must meet the strict native criterion: the self-consistency of each varied orbital must reach ACCY = 5.96e-8 under the `grasp-strict-scf-v1` build, which has no energy-only stop.

Each version had a prospective Linear declaration, a frozen driver and a fresh directory. Each was stopped and preserved at its first failure, with a ledger row and archive `data/audit/rya1354/mcdhf_*_v5x`.

## 1. What each run established

| Run | Design | Outcome |
|---|---|---|
| v53 | all MR orbitals varied together, default method | FAILED_RETAINED (32 s). Method 1 cannot solve 4f/5f; 5f (j=7/2) sits at 0.166 and 2s at 1.7e-5. The EOL energy is stationary to 1e-10, so an energy-only stop would have reported "converged". |
| v54 | s0 integration-method control matrix | FAILED_RETAINED: my stdin omitted the SETDBG "Generate debug output?" answer. |
| v55 | same matrix, corrected stdin | FAILED_RETAINED: native IMPROV stop (6p) counted as a refusal by my own guard. |
| v56 | same matrix; a native stop counts as a measured outcome | **Complete.** No variant converges in either parity. Every stalled orbital lies in a set whose rotations leave the s0 CSF space invariant ({1s,2s}; complete Rydberg manifolds). |
| v57 | test that hypothesis with the standard build from two starts | FAILED_RETAINED: the all-hydrogenic start diverged. The other start permuted the 4S 4p/5p leading labels, showing labels are not invariant under those rotations. |
| v58 | sequential: one spectroscopic orbital per step, earlier ones frozen | FAILED_RETAINED: core step 2s/2p converge, but **1s is stuck at 2.16e-7 on the default grid**. |
| v59 | v58 plus the lifetime-ladder **half grid** (H=0.03125, N=439, same RNT and ACCY) | Core and **3p, 4p, 5p, 6p all converge strictly with correct Li families**. 4f (j=7/2) collapses (e = 0.62 Ha) because the 2p4 ground states are in its EOL. |
| v60 | each Rydberg orbital on its own l-series only, core fixed | Core and np reproduced. In the 4f-only step the selected targets carry **no f7/2** (occupation 4e-8), so the f7/2 orbital is unweighted and wanders to the continuum. |

## 2. Findings

1. **Do not vary all MR orbitals together.** In Li's MCDHF MR without excitations, varying every MR orbital at once leaves rotations within {1s,2s} and within each complete Rydberg manifold (3p–6p, 4f–5f, 3s–6s, 3d–5d) free. The strict criterion cannot converge there, and GRASP's energy stop hides it (v53, v56).
2. **1s needs the half grid.** The default grid (H=1/16, N=220) leaves 1s at about 2e-7, above ACCY, in every run. The half grid fixes it. This was already documented in `lifetime_continuation_v4.md`, and I should have reused it from v53 (ledger lesson).
3. **The sequential approach works for core plus p series.** With the half grid and method 3, the even core and all four np orbitals converge strictly in 8–10 cycles each, twice (v59, v60), with the correct leading LS families.
4. **The open problem is nf.** With the lowest-K-per-J target rule, the 4f step does not select the 4S-parent f7/2 states: the J=5 target is a single CSF with no f7/2. With ground states included instead, 4f collapses. A different target rule is needed for the f series; see section 3. The odd series (ns, nd) has not been run yet, because every run stops at its first failure and even parity runs first.

## 3. Choices before v61

The 2p3 nf states are not in any line we currently grade (130 nm, 135.6 nm, 777 nm, 844.6 nm, 926 nm, 1.13/1.32 µm). In Li's scheme their orbitals feed the MR and the correlation expansion.
- **(a) Explicit f targets.** Select the 4S 4f/5f targets by their leading jj2lsj family from a fixed-orbital CI, not by lowest K per J, and optimise f on exactly those levels. Closest to Li.
- **(b) Fixed hydrogenic nf.** Keep nf in the MR but not optimised: screened hydrogenic, fixed. A documented deviation that is cheap and stable, because Rydberg nf is nearly hydrogenic.
- **(c) Drop nf.** Remove nf from the MR and targets. The largest deviation from Li.

Recommendation: (a), with (b) as the declared fallback if (a) cannot pin f7/2.

## 4. Evidence

- **Archives:**
  - `mcdhf_eol_failed_v53`
  - `mcdhf_s0_controls_failed_v54`, `_failed_v55`, `mcdhf_s0_controls_v56`
  - `mcdhf_s0_invariance_failed_v57`
  - `mcdhf_s0_sequential_failed_v58`, `_v59`, `_v60`
- **Plans:** `config/atomic/rya1354/lifetime_plan.v53.json` … `v60.json`, each with its outcome.
- **Drivers:** `scripts/rya1354_mcdhf_eol.py`, `rya1354_mcdhf_s0_controls.py`, `rya1354_mcdhf_s0_invariance.py`, `rya1354_mcdhf_s0_sequential.py`.
- **Tests:** each archive is replayed by `tests/test_rya1354_mcdhf_*.py`. The tests pin the rmcdhf and rwfnestimate prompt paths to the lifetime ladder's natively proven inputs.
