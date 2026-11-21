# RYA-1354 v47 — literature check, candidate theoretical grade, and Li-scheme sizing

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false. Accepted lifetime, uncertainty, covariance and calibrated sigma remain null. v47 generates configuration lists only: no orbitals, energies, rates or grades of any product.

Ryan (2026-11-17): continue the correlation expansion "along with any other parameters for making sure we have a solid Oxygen model". The end use is theoretical grading of lines, especially UV and IR, with Fe and other elements next. Plan `config/atomic/rya1354/lifetime_plan.v47.json`; prospective declaration Linear `3450dd67-8cf7-440d-9203-a5534c9b3256`. Continues from local `ebd4dc2a` (v46).

## 1. Literature check (before choosing the expansion)

Local reference PDFs were checked first. Li et al. 2023 (O I) is held on Sirius as `sources_curation_v2/rya1354-li2023-arxiv.pdf`.

**How the reference O I calculation is built.** Li, Jönsson, Amarsi, Li & Grumer 2023 (arXiv 2304.02310) use GRASP MCDHF in the extended-optimal-level (EOL) scheme, followed by RCI. Their Table 1:
- **MCDHF multireference (MR).**
  - Even: 2p3 np (n=2–6) and 2p3 nf (n=4,5); n=2 is the 2p4 ground configuration.
  - Odd: 2p3 ns (n=3–6) and 2p3 nd (n=3–5).
  - SD substitutions from 2p (and the outer electron); 1s and 2s closed.
- **RCI MR.** Extends the MCDHF MR: np up to n=8, ns up to n=7, and 2s2p4 ns (n=3–7) and 2s2p4 nd (n=3–6). SD substitutions from every n=2 subshell.
- **Core–valence correlation.** Opening 1s was tested and found negligible.
- **Size.** Final active set {12s,12p,11d,11f,9g,7h}: 12 641 532 even and 7 683 274 odd CSFs.
- **Result.** The 3S1 lifetime is 1.70/1.69 ns (B/C). B/C lifetimes agree to 6% for all selected levels, and ⟨dT⟩ is 2.29% for A ≥ 10⁶ s⁻¹.

The same scheme is used for C I–IV (Li 2021, `2101.09114v1.pdf`) and N I (Li 2023, `li2023_n_nitrogen_rya1220.pdf`). Papoulia et al. 2019 (Al, `1808.09478v1.pdf`) and Pehlivan Rhodin et al. 2024 (Si, `aa45686-22.pdf`; theory lifetimes 5–10%) use the same dT practice.

**How theoretical accuracy is graded** (Li 2023 O I, Sect. 2.2 and Tables 4–5):
- **dT** = |A_B − A_C| / max(A_B, A_C). The authors stress it is an estimate "in a statistical sense for a group of transitions".
- **CF** = (|Σ c_j⟨j‖T‖k⟩c′_k| / Σ|c_j⟨j‖T‖k⟩c′_k|)², taken in the Babushkin gauge (Cowan 1981).
- **d̃T** = max(dT, ⟨dT⟩ of the line's A-magnitude group), using six groups.
- **Class.** NIST-style A ≤ 3%, B ≤ 10%, C ≤ 25%, D ≤ 50%, E > 50%. CF < 0.1 moves a line down one class.
- **Their outcome.** 205 of 989 O I transitions are in classes AA–B.

**Our method compared with theirs.** v44 showed the gauge split is carried by a ~1% 2p3 3p admixture in the lower level. Li puts exactly those configurations (2p3 np) into the MR, with SD correlation built on top of them, and adds the 2s2p4 families. Neither our FAC models nor our GRASP n5 VV model have them in the MR. That GRASP model stopped at 2.284/1.713 ns, or dT ≈ 25%.

FAC builds every orbital from one local central potential, so it cannot provide variationally contracted correlation orbitals. The expansion therefore goes to GRASP, which is already the plan's declared primary method (`benchmark_plan.v1`). FAC remains the independent cross-check, and the v44 component instrument remains the diagnostic.

Two earlier positions from this ticket still stand. Li's energy-adjusted rates are not adopted (`lifetime_plan.v4`), and no gauge spread becomes a calibrated sigma.

## 2. Candidate theoretical grade (`scripts/rya1354_theory_grade.py`)

This implements dT, the six-group d̃T, Babushkin-gauge CF, and the Table 5 class rule. CF is computed from v44 component records, with terms summed per CSF pair before taking magnitudes, as in Li's Eq. 4. It is always ≥ 1/cancellation_index², because v44's index sums finer terms. Status: `CANDIDATE_HOLD_NOT_ADOPTED`. The RYA-1352 THEORETICAL_CODEX thresholds and adoption need Ryan's ruling. Nothing writes a grade into a product.

**Baseline on the v46 FAC models (3S1 → lower, Babushkin CF), diagnostic only:**

| Lower | ground dT | ground CF | class | equal dT | equal CF | class |
| --- | ---: | ---: | :---: | ---: | ---: | :---: |
| 3P2 | 20.4% | 0.561 | C | 50.7% | 0.863 | E |
| 3P1 | 20.3% | 0.507 | C | 50.7% | 0.860 | E |
| 3P0 | 20.3% | 0.498 | C | 50.6% | 0.852 | E |
| 1D2 | 0.6% | 4.5e-6 | **B** | 41.3% | 7.4e-6 | E |
| 1S0 | 45.5% | 2.6e-5 | E | 79.3% | 2.8e-5 | E |

**Finding for the grading ruling.** Under Table 5 as published, a cancellation-dominated line with a coincidentally small dT can be class B: CF = 4.5×10⁻⁶ costs only one class. Li's group-averaged d̃T protects against this only when the A-magnitude group holds many lines. A CF floor (for example CF < 10⁻³ → E), or a minimum group size, would be a deliberate departure from Li. It is not adopted here. The table also gives the target for the expansion: the dominant 3P lines must move from C/E to A/B.

## 3. v47 sizing (Sirius, native `rcsfgenerate` only)

Driver `rya1354_li_scheme_sizing.py`, sha256 `679d8a3a…320c630b`; GRASP pinned `9006157730a8`. The build is serial: there are no MPI binaries and no `mpirun` on Sirius (4 cores, 15 GB). The batch exited 0 in 15:08 at 739952 KiB, status `SIZING_COMPLETED_NO_ATOMIC_PREDICTION`.

All 16 cases passed the streamed native inventory: 8 electrons, closed 1s (closed 2s for the MCDHF stage), parity, no duplicates. MR orbitals above n_max are omitted and recorded per case: 3 even and 1 odd at n6, 1 even at n7. CSFs per J block:

| n_max | MCDHF even 0+/1+/2+ | MCDHF odd 1−/2− | RCI even 0+/1+/2+ | RCI odd 1−/2− |
| ---: | --- | --- | --- | --- |
| 6 | 7 688 / 21 193 / 30 468 | 17 977 / 25 617 | 75 419 / 211 298 / 308 729 | 142 302 / 205 391 |
| 7 | 14 521 / 40 403 / 58 825 | 34 084 / 49 327 | 170 953 / 482 042 / 711 992 | 288 437 / 421 220 |
| 8 | 23 763 / 66 450 / 97 391 | 55 625 / 81 113 | 315 856 / 893 259 / 1 325 838 | 477 181 / 700 902 |
| 9 | 35 414 / 99 334 / 146 166 | 82 600 / 120 975 | 499 982 / 1 416 230 / 2 107 505 | 713 955 / 1 052 306 |

**Reading.**
- **MCDHF stage.** Every layer through n9 is at most 146k CSFs per block, comparable to blocks GRASP handles serially. Angular-data size and SCF time are still unmeasured.
- **Full RCI stage.** Already 0.3M CSFs per block at n6 and 2.1M at n9. A serial build with 15 GB cannot treat this as it stands.
- **Feasibility levers to measure next, in the literature's own toolkit:**
  - `rcsfinteract`: keep only CSFs that interact with the MR.
  - A weight-selected RCI MR: Li N I keeps configurations with LSJ weight > 0.05.
  - A per-block cost benchmark of `rangular`/`rci` before any large run.

## 4. Evidence

`data/audit/rya1354/li_scheme_sizing_v47/`: 110 parts, 13 870 234 compressed bytes, all ≤ 5 MiB.
- **Archived:** 114 native files, including every stdin/stdout/stderr, the receipts, the frozen driver and all eight MCDHF-stage CSF lists.
- **Off-archive:** the eight RCI-stage lists (2 139 542 740 bytes). They stay on Sirius in `runs/li_scheme_sizing_v47/`. Each is declared in the manifest with sha256, bytes and its archived `rcsfgenerate.stdin`, so it is regenerable with the pinned executable.

**Defect, recorded and fixed** (run-bug ledger `li_scheme_sizing_v47_receipt_files`). The new driver did not write the receipt `files` inventory that earlier drivers write, so the decoder refused the pack.
- The native receipt was not edited.
- A read-only post-run sidecar (`evidence/li_scheme_sizing_inventory_v47/files.json`, 115 files) is archived. It independently matches all 16 rcsf hashes the driver recorded during the run.
- The decoder gained an opt-in, manifest-declared `receipt_inventories` / `off_archive_regenerable` path. It verifies every archived file against the sidecar, every off-archive file against both the manifest and the sidecar, and that each off-archive file's regeneration input is archived.
- A receipt with no `files` list and no declared sidecar still fails. Earlier archives are unaffected.
- Future drivers must write the receipt files list.

## 5. Validation

`pytest -q tests/test_rya1354*.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py` → **291 passed in 269.48 s**. New tests:
- the streaming inventory equals the tested `inventory()` on all 8 native v5 lists, and keeps its refusals;
- Li Table 1 families and layer truncation;
- the Table 5 class boundaries, d̃T groups, and CF summed per CSF pair;
- CF ≥ 1/index² on the actual v46 components;
- sidecar acceptance plus six refusals, and that a receipt with no files list still fails;
- the actual v47 archive replay.

The ledger row passes the scoped schema checks; the whole-ledger rows that already failed are not rewritten. No canonical surface, publication, merge or Done.

## 6. Next bounded step (proposal)

**v48 — RCI feasibility levers**, CSF lists and cost only:
- `rcsfinteract` reduction of the n6/n7 RCI lists against the Li RCI MR;
- one-block `rangular` + `rci` time and memory benchmark at increasing size, so the scaled scheme is chosen from measurements.

**v49 — the first physical layer.** GRASP MCDHF-EOL with the Li MCDHF MR, layers n=3→6, strict SCF, both gauges. The v44-style component audit tracks the 2p3 np ← 2p3 3s pair and the B/C ratio at each layer.
