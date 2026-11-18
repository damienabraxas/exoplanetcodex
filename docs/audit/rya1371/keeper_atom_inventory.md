# Keeper model-atom inventory for OpenMULTI3D products (RYA-1371)

2026-11-18. Ryan's direction: OpenMULTI3D is a publication-pipeline engine. It becomes a
product for every element the ORCH has already run, and the Bride fills only gaps where no
3D model exists. This inventory is the Phase 1 input check: which model atoms exist, which
version we hold, and what is missing.

## Source and snapshot

- Source: MPG Keeper share `6eaecbf95b88448f98a4`, `/dep-grids/<El>/atom.*`. The share has
  no DOI and is mutable. There is no `/model_atoms` folder; the atoms sit inside each
  element's `dep-grids` folder.
- Snapshot: Sirius `grids.keeper_atoms_20261118` = `/mnt/codex-data/grids/nlte/keeper_atoms_2026-11-18/`,
  fetched 2026-11-18T07:31:28Z. It holds 19 files, made read-only, with full SHA-256 sums in
  `SHA256SUMS`. The existing `grids.gerber_ts` deck (RYA-534) was not touched.
- Compared against: the Gerber deck atoms (RYA-534) and the four atoms shipped with OpenMULTI3D
  at pin `32085d26` (`input/atoms/`).

## Atoms per element

Counts are levels / bound-bound lines / continua, from each file's header.

| El | Keeper atom (date on Keeper) | Levels/BB/BF | vs Gerber deck (Sirius) | ORCH run? |
|---|---|---|---|---|
| H | h20 (2021-08-02) | 20/171/5 | identical | — |
| O | **o41f (2026-11-17, CHANGED)**, o41f_ref, **o41j (NEW)** | 41/256/22 | o41f_ref = held o41f, byte-identical | CNO |
| Na | na_qmh (2023-07-18) | 290/584/289 | identical | yes |
| Mg | **mg86g_nq_qmax (2026-10-26, NEW)** | 86/245/85 | deck holds older mg86d | yes |
| Al | al_qmh (2023-07-27) | 354/2543/353 | identical | yes |
| Si | si340 (2022-03-10) | 340/2956/289 | identical | yes |
| Ca | ca105b (2021-08-02) | 105/624/104 | identical | yes |
| Ti | ti503b (2026-03-16) | 503/10763/502 | identical | yes |
| Cr | cr374 (2024-12-12) | 374/5879/373 | identical | yes |
| Mn | mn281kbc (2021-08-02) | 281/1689/198 | identical | — |
| Fe | fe607a (2021-08-02) | 607/12635/606 | identical | yes |
| Co | co247qm (2024-11-14) | 247/6027/246 | identical | — |
| Ni | ni538qm (2022-06-02) | 538/3357/196 | identical | — |
| Sr | sr191 (2022-03-03) | 191/549/190 | identical | — |
| Y | y423qm_pf (2024-01-15) | 424/11819/187 | identical | — |
| Ba | ba111 (2021-08-02) | 111/284/110 | identical | — |
| Eu | eu662qmkhfs (2024-01-15) | 662/162/661 | identical | — |

The `CH` folder on Keeper holds no atom. There is **no Keeper atom for C, N, K, S, Sc, V, Cu
or Zn.**

## What changed, and why it matters

1. **O, updated on Keeper 2026-11-17.** The old file was renamed `o41f_ref`; it is
   byte-identical to our held atom (`d058b75b…`).
   - The new `o41f` (`5cc8a84d…`) adds a header: "Barklem + Kaulakys data (co-added),
     closed by Drawin". It also zeroes one broadening value on line 4→1 ([O I] 6300,
     2p4 ¹D₂→³P₂).
   - The new `o41j` (`9dc06ba0…`) is "Belyaev data, Kaulakys data, closed by Drawin".
   - H collisions are what separates the 3D NLTE literature values (A21 8.69 vs
     Steffen/Bergemann/Magg 8.75–8.77; Amarsi & Grevesse 2026), so both files go into the
     oxygen run as named variants.
2. **Mg.** The deck's `mg86d` and OpenMULTI3D's `mg86g` are both superseded by
   `mg86g_nq_qmax`: "MB 01.06.2026: updated RBF to cover the full range of x-sections".
3. **Copies shipped with OpenMULTI3D are not the Keeper files.** Its `fe607a` (`9715c0de…`)
   differs from Keeper's (`53446fe1…`) only in the element-name case and the frequency
   sampling (NQ/QMAX) of 51 lines; the physical data are unchanged. Its `o41f` differs the
   same way, on one line (NQ 9→101). Products must name which file they used.

## Coverage of ORCH elements already run

- **Atoms available:** Na, Mg, Al, Si, Ca, Ti, Cr, Fe, O. Every one has a Keeper atom, so
  OpenMULTI3D products are possible for all of them.
- **Missing: C and N.** No Keeper atom exists. These are the first Bride gap targets, or
  need another published atom (to be checked before building one).
- **Resource load (number of levels):**

  | Atoms | Levels |
  |---|---|
  | Fe | 607 |
  | Ti | 503 |
  | Cr | 374 |
  | Al | 354 |
  | Si | 340 |
  | Na | 290 |

  These need native OpenMULTI3D runs, not the Bride (20,000-wavelength cap).

## Licensing

The Keeper share carries no explicit data licence. We cite the source papers in each atom's
header and do not redistribute raw atoms in Git; atoms stay Sirius-only.
