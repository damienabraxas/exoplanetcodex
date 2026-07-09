# RYA-540 — relocation + persistent md5-pinned grid cache (execution record)

**Date:** 2026-07-08 · **Host:** Sirius · Companion to `manifest.md` (STEP-0 audit).
Governing rule (Ryan, standing): *Sirius holds the entire pipeline working set locally.
Re-downloading an already-acquired artifact is a DEFECT.* Capacity is solved with STORAGE,
not re-downloads. Eviction is an explicit admin action, never a workflow step.

Drive policy (Ryan): **stars / spectra → SSD (sda, `/mnt/codex-data`); grids, software,
models, atlases, engines → M.2 (sdb, `/`).**

## 1. Relocation (SSD → M.2), symlink-compat, verify-before-delete

The whole NLTE working set was moved from the SSD to the M.2 under `/srv/codex/`, then a
**compat symlink** was left at each old `/mnt/codex-data/...` path so the deck's hardcoded
paths and the worktree `.grd` symlinks keep working with **zero code change**. Each move was
`rsync` → byte-count + md5 verify → rename original to `*_OLD_rya540` (kept) → symlink.

| Moved tree | Size | New home (M.2) | Compat symlink |
|---|---|---|---|
| `grids/` (amarsi_galah 69 G + gerber_ts + model_atmospheres) | 70 G | `/srv/codex/grids` | `/mnt/codex-data/grids` → |
| `solar_reference/` (IAG + Wallace/KPNO atlases) | 729 M | `/srv/codex/solar_reference` | `/mnt/codex-data/solar_reference` → |
| `engines/Turbospectrum_NLTE` | 99 M | `/srv/codex/engines/Turbospectrum_NLTE` | symlink → |
| `engines/ispec_src` | 252 M | `/srv/codex/engines/ispec_src` | symlink → |
| `engines/TSFitPy` | 83 M | `/srv/codex/engines/TSFitPy` | symlink → |
| `engines/pysme_src` | 21 M | `/srv/codex/engines/pysme_src` | symlink → |

Verification: `grids` rsync SRC==DST==74711442532 bytes, 15765 files both sides; Ca grid
md5-match; TS_NLTE `bsyn_lu` + `solar_reference` read-throughs resolve via symlink. The deck
path now backs onto **`/dev/sdb2` (M.2)** — confirmed by `grid_cache.py --status`
(`backing device: /dev/sdb2`, 343.5 GB free).

**`spectra/` stays on the SSD** (correct per policy). **`engines/nirps_drs` (131 G) stays on
the SSD** — it is the NIRPS reduction DRS + bundled star frames/cal data (RYA-498), i.e.
star-data + reduction scratch, not an NLTE-synthesis working-set item; co-locating it on the
M.2 would overflow (manifest capacity check: ~392 G needed vs 391 G free with it included).

### `*_OLD_rya540` backups — reclaim PENDING explicit confirm
Originals are retained on the SSD (verify-before-delete): `grids_OLD_rya540` 70 G,
`solar_reference_OLD_rya540` 729 M, `engines/{ispec_src,TSFitPy,pysme_src,Turbospectrum_NLTE}_OLD_rya540`
~455 M — **~71 G total**. Reclaiming (`rm -rf`) frees the SSD but is **not done here**; it is
an explicit, owner-confirmed step (nothing irreversible without a go).

Post-move `df`: SSD (sda) 207 G used / 228 G free; M.2 (sdb) 124 G used / 320 G free.

## 2. `scripts/grid_cache.py` — the persistent cache (RYA-540 core)

Policy (replaces RYA-534's free-after-gate):
- **Acquire once:** on a cache MISS, download the `.bin.zip` → verify against the ZIP md5
  pinned in the element's provenance JSON → unzip → record the unzipped `.bin` md5 in the
  cache index (`gerber_ts/_cache_index.json`) → **keep the `.bin`, delete the `.zip`.**
- **Reuse forever:** on load, md5 the cached `.bin` vs the recorded index md5. **HIT →
  skip the download entirely.** MISMATCH (corrupt/truncated — the RYA-534 concurrent-unzip
  class) → RAISE + re-fetch **that grid only**.
- **No free-after-gate:** the `.bin` is never deleted by the pipeline. `evict(el)` is an
  explicit admin-only reclaim.
- **Capacity = loud-fail, never thrash:** before a download, check free space; a shortfall
  RAISES with the number and points at RYA-477 (external SSD / Orion). A full disk during
  download fails loudly (curl no-space); the post-download unzip check RAISES before unzip.

Provenance is the **single source of truth** for the `.bin` filename + md5. Schema note: the
key is inconsistent across tickets — RYA-533 (Na) used `files.grid_1d_bin`, RYA-534 (the
other 10) used `files.grid_1d_bin_zip`; the loader accepts either. **Known gap:** the 10
RYA-534 prov entries omit a `bytes` field, so the *pre*-download capacity guard falls back to
a floor for them (post-download md5 + unzip-capacity guards still fire). Backfill exact
`bytes` into those prov JSONs after the next real pull.

Self-contained smoke on Sirius (`/tmp/rya540_smoke.py`, synthetic tiny grid via `file://`,
no Keeper pull): **10/10 pass** — miss→download→md5-verify→unzip→keep; HIT→no re-download;
corrupt→mismatch→re-fetch; capacity→loud-fail (RYA-477); no-prov→refuse.

## 3. Deck wiring (`ts_gerber_gate.py::make_departure`)

The deck previously did `binf = f"{GT}/{cfg['grid']}"` then
`if not os.path.exists(binf): raise SystemExit("grid .bin missing (provision first)")` — it
never downloaded, and the RYA-534 runner freed grids after each gate, forcing re-provision.
RYA-540 replaces that with:

```python
from grid_cache import ensure_grid   # RYA-540
...
    binf = str(ensure_grid(el))       # cache-or-download, md5-verified, retained
    aux = f"{GT}/{cfg['aux']}"
```

This makes the prov JSON (not a hardcoded `cfg['grid']` string) the SSOT for the `.bin` name.
Cross-check: all **10/10** RYA-534 prov `bin_name`s match the deck `cfg['grid']` exactly, so
the swap is behaviour-preserving on a HIT and download-once on a MISS.

Applied **live on Sirius** (`/mnt/codex-data/codex/rya534/ts_gerber_gate.py`, the standalone
working deck) so the running pipeline stops re-downloading immediately; verified: deck parses,
`ensure_grid` resolves, and with the departure file absent `make_departure` routes grid
loading through `ensure_grid(el)`. The same 2-line diff is tracked here to apply to the deck
on the RYA-534 branch when 534/540 integrate (the deck is not on the rya540 branch).

## Follow-ups
- Reclaim `*_OLD_rya540` (~71 G) after owner confirm.
- Backfill exact zip `bytes` into the 10 RYA-534 prov JSONs (strengthens the pre-download guard).
- The one-time real re-pulls (Ti 55 G etc., RYA-535) now land in the persistent cache and are
  never re-downloaded again — the defect RYA-540 was filed to end.
