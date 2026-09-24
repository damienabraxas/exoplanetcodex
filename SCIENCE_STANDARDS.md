# Codex data-stewardship standards

## 1. What GitHub carries — and what it never carries

**GitHub carries the Codex's own work.** Measurements and derived products (built line
lists, product feeds, results artifacts), code, config, and small provenance manifests.

🔴 **Raw external data that is NOT our measurement is NEVER committed or LFS-tracked.**
It lives under `~/Documents/Exoplanet Codex` (RYA-461) and/or Zenodo (RYA-1124), and is
referenced from a provenance manifest in the repository. This covers VALD extractions,
HARPS/ESO/any-archive spectra (FITS), model-atmosphere grids, and any archive-sourced or
re-fetchable input.

*Standing policy, Ryan, 2026-09-18 (RYA-1228).*

### Why this rule exists, stated as what actually happened

The GitHub LFS quota was exhausted by a single `.gitattributes` line. `data/linelists/
vald_*_raw.txt` was the only LFS-tracked pattern in the repository and it accounted for
the entire footprint: **53 objects, 853.9 MiB**, of which **30 were reachable only from
history** and invisible to `git lfs ls-files`.

⚠️ **It was never a `.gitignore` failure.** Spectra and model grids were already ignored
correctly. VALD came in through `.gitattributes` LFS tracking — a *separate route* that
`.gitignore` does not see. Closing one route does not close the other, which is why the
enforcement below checks the **committed tree**, where both routes converge.

⚠️ **And none of it was re-fetchable.** VALD3 is a manual web extraction with no API, so
every object had to be preserved with a verified checksum *before* it could be purged.
"We can re-download it later" was never available.

## 2. Enforcement — structural, not documentary

| layer | mechanism |
|---|---|
| `.gitignore` | `data/linelists/vald_*_raw.txt`, `vald_*.vald`, `vald_procyon_raw/`; FITS and model grids already covered |
| `.gitattributes` | **no `filter=lfs` pattern may exist.** The repository tracks no data in LFS |
| pre-commit | `scripts/check_stewardship.py` → `check_no_raw_external_data()` refuses the commit |
| provenance | `data/linelists/VALD_PROVENANCE.md` — one row per raw delivery, with its sha256 and archive location |

The guard applies **two independent tests**, because each alone has a known blind spot:

* **by pattern** — catches the file types we already know about, and by construction
  misses the next one;
* **by size ceiling** (5 MiB, non-CSV/non-code) — catches the next one by shape, and would
  misfire on a legitimate product, so tabular and code suffixes are exempt.
  `linelist_solar.csv` is 25 MB and entirely legitimate.

## 3. What stays

`linelist_*.csv`, `canonical_gf.csv` and every product feed are **Codex measurements**.
They are plain committed files, they were never LFS, and nothing in this standard touches
them. The distinction is authorship, not size: we commit what we measured.
