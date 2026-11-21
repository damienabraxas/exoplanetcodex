# RYA-1354 v48–v52: RCI reduction and cost feasibility (Li 2023 scheme, GRASP)

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. These runs used Thomas-Fermi or screened-hydrogenic orbitals. They measure CSF counts, wall time, memory and angular disk only. No energy, rate, lifetime or grade from them is a result.

**Why.** v47 sized Li et al. 2023 Table 1 at 0.3–2.1M CSFs per full RCI block, which is too large for a serial GRASP build with 15 GB of memory. Before choosing the scaled scheme, we measured two things on Sirius:
- **A.** How much the literature's own reduction removes. This is `rcsfinteract` (Gaigalas & Froese Fischer), in its DC and DCB modes.
- **B.** What one reduced block costs in `rangular` and `rci`.

Each attempt had its own plan, `config/atomic/rya1354/lifetime_plan.v48.json` … `v52.json`, and its own prospective Linear declaration (a595852c, 5bcec254, 305d1629, 02cf4798, d2dd1f15). The driver is `scripts/rya1354_rci_feasibility.py`, frozen per version in `evidence/rci_feasibility_driver_vNN` on Sirius.

## 1. Retained failures (FAILED_RETAINED, each with a run-bug ledger row)

| Run | Stopped at | Cause | Fix in next version |
|---|---|---|---|
| v48 | first `rcsfgenerate` | GNU `time -o` relative path resolved against the child's cwd | absolute paths |
| v49 | guard before `rcsfinteract` | a 0-excitation MR lists only occupied orbitals, so its peel header is not a prefix of the full list's header | `rcsfmr.inp` = full header + native MR body, after an ordered-subset check |
| v50 | section B n4 `rangular` link | `rwfnestimate` writes `rwfn.inp`, not `rwfn.out` | `ESTIMATE_OUTPUT='rwfn.inp'` |
| v51 | section B n5 `rwfnestimate` | Thomas-Fermi cannot build 5g, and the next prompt hit EOF (rc=2) | screened-hydrogenic fallback for the remaining subshells |

The archives are `data/audit/rya1354/rci_feasibility_failed_v48` … `_v51`. v50 and v51 both completed section A and gave identical results.

## 2. Section A: what rcsfinteract keeps (v50 = v51)

The table gives the CSFs kept by DC / DCB as a fraction of the v47 list. Every regenerated full list matched its v47 sha.

| n_max | even 0+ | even 1+ | even 2+ | odd 1− | odd 2− |
|---|---|---|---|---|---|
| 6 | 0.70 / 0.73 | 0.79 / 0.82 | 0.79 / 0.82 | 0.58 / 0.62 | 0.58 / 0.62 |
| 7 | 0.64 / 0.66 | 0.74 / 0.77 | 0.75 / 0.79 | 0.57 / 0.60 | 0.57 / 0.61 |

The largest DCB-reduced blocks are even 2+: 253 482 at n6 and 559 970 at n7. Each reduction ran in 0.5–2 min at about 4 MB RSS.

**Finding: `rcsfinteract` removes only 20–40%.** It discards CSFs that cannot interact with the MR. Because the Li MR is large, almost everything interacts with it.

## 3. Section B: per-block cost (v52, complete)

The block is odd J=1, which holds the 3S1 upper level, reduced with DCB. Every native step ran under `prlimit --as=12 GiB`. No budget or memory stop occurred.

| n_max | CSFs | rangular wall / peak RSS | rci wall / peak RSS | angular files |
|---|---|---|---|---|
| 4 | 7 596 | 16.4 s / 0.56 GB | 40.9 s / 58 MB | 0.11 GB |
| 5 | 33 575 | 282 s / 6.31 GB | 495 s / 0.36 GB | 1.12 GB |
| 6 | 87 723 | 2 293 s / 6.68 GB | 4 165 s / 1.47 GB | 5.02 GB |

The v51 run of the same n4 block took 32.7 s in `rangular`, so wall times carry roughly 2× noise from machine load.

**Scaling from n5 to n6.** `rangular` and `rci` times both grow as about N^2.2. Angular bytes grow as N^1.56. The `rangular` peak RSS levels off near 6.7 GB, so memory is not the binding limit; time and angular disk are.

**Extrapolations** (not measured). Each is for one block, with one rci eigenvalue and a serial build:

| Block | CSFs | rangular | rci | angular disk |
|---|---|---|---|---|
| even 2+ n6, DCB | 253k | ~7 h | ~12 h | ~26 GB |
| even 2+ n7, DCB | 560k | ~38 h | ~71 h | ~90 GB |
| full Li blocks | 2.1M | months | months | ~0.7 TB |

Sirius had 108 GB free.

## 4. Conclusion for the O I model

1. The Li 2023 RCI MR with SD excitations is feasible here up to about n6, at roughly half a day per large block. It is not feasible at n7 or above as generated, and `rcsfinteract` does not change that.
2. The literature's next lever is a **weight-selected MR**. GRASP `rcsfmr` keeps the configurations above a weight threshold from a converged calculation; Li N I uses weights above 0.05. That needs physical orbitals first.
3. Next declared step (v53): GRASP **MCDHF-EOL** with the Li MCDHF MR, layers n=3→6, both gauges, plus a component audit at each layer. Its converged n6 mixing gives the weight-selected RCI MR. The section-A/B machinery here is reused for sizing that list before any RCI.

## 5. Evidence

- **v52 archive.** `data/audit/rya1354/rci_feasibility_complete_v52` holds 9 parts, 158 stored files and 8.4 MB.
  - 37 files of 4 MiB or more (8.1 GB) stay on Sirius under `runs/rci_feasibility_v52`. These are lists, angular files and `rci.res`.
  - The manifest's `off_archive_regenerable` gives each one's sha256/bytes and the archived stdin that regenerates it. The receipt's `files` list vouches for the same hashes.
- **Packer.** `scripts/rya1354_pack_rci_feasibility.py`.
- **Tests.** `tests/test_rya1354_rci_feasibility.py` replays every archive. The v52 test cross-checks the n6 J=1 count against the v51 section-A receipt, which is pinned by its sha.
