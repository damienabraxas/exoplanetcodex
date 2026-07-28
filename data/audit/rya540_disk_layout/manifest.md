# RYA-540 STEP-0 — Sirius disk-layout manifest (read-only audit)

**Date:** 2026-07-09 · **Host:** Sirius · **Scope:** audit + manifest ONLY. No cache module,
no file moves, no policy change yet (per "report the manifest before touching the cache policy").

## Drives (two SSDs)

| Drive | Model / form factor | Mount | Total | Used | **Free** | Intended role |
|---|---|---|---|---|---|---|
| **sdb** (sdb2) | NT-512 **2280 = M.2** | `/` (OS) | 468 G | 53 G | **391 G** | **M.2 (OS + software):** grids, model atmospheres, atlases, engines |
| **sda** (sda1) | WD Blue SA510 **2.5″ SATA SSD** | `/mnt/codex-data` | 458 G | 207 G | **228 G** | **SSD (data):** spectra (FITS) |

The M.2 (`/`, 391 G free) holds **zero pipeline artifacts** — only OS/home (ollama LLM blobs
~10 G, chrome cache, `/swap.img`). **The entire pipeline working set is on the SSD (sda).**

## Artifact trees (all currently on sda `/mnt/codex-data`)

| Tree | Size | Intended drive | Correct? | Action |
|---|---|---|---|---|
| `grids/nlte/amarsi_galah` (13 PySME `.grd`, Engine-A) | **69 G** | M.2 | ✗ on SSD | → move to M.2 |
| `grids/nlte/gerber_ts` (atoms + aux; `.bin` FREED by RYA-534) | 0.7 G | M.2 | ✗ on SSD | → move to M.2 |
| `grids/model_atmospheres/marcs_standard_comp` | 0.5 G | M.2 | ✗ on SSD | → move to M.2 |
| `engines/nirps_drs` (NIRPS DRS + bundled data, RYA-498) | **131 G** | M.2? (software) — **FLAG** | ✗ on SSD | DECISION (see capacity) |
| `engines/{ispec_src, Turbospectrum_NLTE, TSFitPy, pysme_src}` | 0.46 G | M.2 | ✗ on SSD | → move to M.2 |
| `solar_reference` (IAG Baker-2020/Reiners-2016 + Wallace KPNO atlases) | 0.7 G | M.2 (atlases) | ✗ on SSD | → move to M.2 |
| `spectra` (tau_boo 1.9 G, procyon 431 M, _quarantine 50 M) | 2.3 G | SSD | ✓ correct | keep on SSD |
| `codex/` worktrees (repo 729 M, rya519 274 M, …) + `venv312`/`venv_pysme` | ~4.4 G | software | — | keep (see flags) |

Per-grid `.grd` (amarsi_galah, Engine-A PySME set): Ca 9.6 G, Mn 8.5 G, Si 8.3 G, C 7.4 G,
Li 7.1 G, N 5.8 G, Mg 5.8 G, O 4.4 G, K 3.5 G, Ba 3.3 G, Na 2.6 G, Al 2.2 G, H 0.3 G = **69 G**.

## Flags

1. **WRONG DRIVE (the headline):** the NLTE working set — grids (70 G) + model atmospheres (0.5 G)
   + atlases (0.7 G) + synth engines (0.46 G) — is entirely on the **SSD (sda)**, but the intended
   drive is the **M.2 (sdb, `/`)**, which is nearly empty (391 G free). Spectra (2.3 G) ARE correctly
   on the SSD. So the fix is a **relocation to the M.2**, then build the cache against the M.2 paths.
2. **`engines/nirps_drs` = 131 G** dominates the footprint. It is the NIRPS **reduction DRS + bundled
   calibration data** (RYA-498), not an NLTE-synthesis working-set item. Whether it belongs on the M.2
   (as "software") or stays on the SSD / is pruned is a **decision** — and it is the capacity swing
   factor (below).
3. **No duplicates:** the worktree `.grd` under `codex/rya519/...` are **symlinks** into canonical
   `grids/nlte/amarsi_galah/` (verified), not copies. (Minor: `codex/data/spectra` is a small 51 M
   worktree data dir, not a copy of the 2.3 G spectra.)
4. **No orphaned partial/truncated `.bin`:** the RYA-534-freed Gerber `.bin` grids are cleanly gone
   (only atoms/aux remain in `gerber_ts`); no leftover concurrent-unzip-truncation artifacts.
5. **The Gerber TS `.bin` grids are ABSENT** (freed) — exactly what RYA-540's cache must retain.
   Re-acquiring the full 11-element set is large: known Ti 55 G, Sr 38 G, Si 24 G, Na ~16 G, O 6.3 G;
   remainder (Ba/Ni/Co/Mn/Ca/Mg) ~5–15 G each → **full Gerber set ≈ 150–200 G** to hold locally.

## Capacity check (the governing "add storage, don't re-download" rule)

M.2 free = **391 G**. Prospective working sets to hold there:

| Scenario | M.2 footprint | Fits in 391 G free? |
|---|---|---|
| NLTE set only (both grid families + atmospheres + atlases + synth engines), **excluding** nirps_drs — current 71 G + full Gerber ~190 G | **~261 G** | ✅ yes, ~130 G headroom |
| + nirps_drs (131 G) also on M.2 | **~392 G** | ❌ **over** the 391 G free (0 headroom) |

**Verdict:** the NLTE grid working set (both families kept locally per RYA-540) **fits comfortably on
the M.2 with ~130 G headroom IF nirps_drs (131 G) is NOT co-located there** (it can stay on the SSD or
be pruned — it's a reduction pipeline, not NLTE synthesis). **If nirps_drs must also live on the M.2,
the combined set exceeds M.2 free → this is the "add capacity" trigger** (external SSD now / Orion
later, per RYA-477), NOT a reason to keep freeing/re-downloading. **Number for Ryan: ~392 G needed vs
391 G free with nirps_drs included; ~261 G vs 391 G without it.**

## Next (after Ryan's steer — NOT done here)

1. Decide nirps_drs placement (drives the capacity answer).
2. Relocate the misplaced NLTE working set (grids + atmospheres + atlases + engines) from the SSD to
   the M.2; keep spectra on the SSD; update the deck's grid/engine paths (currently hardcoded to
   `/mnt/codex-data/...`) to the M.2 location.
3. Build the md5-pinned persistent cache (rest of RYA-540) against the corrected M.2 paths.
