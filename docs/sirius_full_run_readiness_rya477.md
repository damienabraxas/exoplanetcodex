# Sirius full-run readiness manifest — RYA-477

**Goal:** one place that says everything a full benchmark run needs, where it lives
(Mac vs Sirius), and the gaps — the reproducibility checklist feeding RYA-337 and the
RYA-259 cross-platform-verification prerequisite.

**Sirius:** HP ProBook 450, Ubuntu (kernel 7.0), Python 3.14.4, reached over `ssh sirius`.
**Data drive (RYA-419):** `/mnt/codex-data` — 458 GB ext4, UUID-mounted via fstab with the
`.codex_mounted` sentinel; **419 GB free**. Everything big lands here (grids, linelists,
engine source, spectra) — never the system disk.

> Status legend: ✅ on Sirius / verified · ⏳ downloading now (Part B) · 📦 in-repo,
> arrives via Part D repo transfer · ⚠️ gap/owed · ⛔ externally blocked.

---

## Part A — Full-run dependency manifest

### A.1 Software / engines
| Dependency | Version / source | Size | On Sirius |
|---|---|---|---|
| iSpec (bundles Turbospectrum / MOOG / SPECTRUM backends) | github.com/marblestation/iSpec | source 109 MB (bundle w/ data ~9.3 GB on Mac) | ✅ source cloned `/mnt/codex-data/engines/ispec_src`; compile = RYA-172 |
| PySME (SME) | github.com/AWehrhahn/SME | 21 MB | ✅ source cloned `/mnt/codex-data/engines/pysme_src` |
| MOOG / Turbospectrum / SPECTRUM | vendored inside iSpec | — | ✅ via iSpec source; **Fortran compile owed (RYA-172)** |
| molecfit (esorex 4.4.2) | ESO / Apptainer | — | ⚠️ RYA-375 route (not fetched here — ESO recipe, separate) |
| Python env + pinned reqs | `requirements.txt` (Mac py3.9.6) | — | ⚠️ **version skew: Sirius py3.14.4 vs Mac 3.9.6** — pin/venv owed for the RYA-172 compile |
| Repo working tree | github (private) | — | 📦 Part D transfer (Sirius lacks `gh`/verified github key — RYA-80) |

### A.2 NLTE / abundance grids (the big-disk items)
| Grid | Source | Size | On Sirius |
|---|---|---|---|
| **Amarsi GALAH PySME `.grd` set (13 elements: Al Ba C Ca H K Li Mg Mn N Na O Si)** | Zenodo 3982506 (Amarsi+2020, A&A 642 A62) | 5.34 GB archives → ~35 GB extracted | ⏳ **Part B, md5-verified per item** → `/mnt/codex-data/grids/nlte/amarsi_galah/` |
| ↳ incl. `nlte_Mn_scatt_pysme.grd` (8.4 GB, PASS-load-bearing, RYA-476) | Zenodo 3982506, md5 `ba0159…ebf72` | 1.40 GB → 8.4 GB | ⏳ in the set |
| Amarsi 2019 C I/O I CNO grid | CDS (RYA-359), vendored | 7.2 MB | 📦 in-repo `data/nlte_grids/amarsi2019_cno/` → Part D |
| MPIA/Bergemann Fe NLTE grid | RYA-319, committed CSV | 0.9 MB | 📦 in-repo `Fe_Bergemann_MPIA.csv` → Part D |
| MPIA Ca / Ti / Cr / Mn / Mg / Si grids | RYA-245/411, committed CSVs | <0.5 MB | 📦 in-repo → Part D |
| Other repo NLTE CSVs (Ba Korotin, Na Lind, Sr, Al/K/Mg/Na/Si/S Amarsi-PySME-CSV) | committed | small | 📦 in-repo → Part D |
| STAGGER / MARCS.GES model-atmosphere grid | RYA-457 | — | ⛔ externally blocked (RYA-457) — skeleton dir `grids/model_atmospheres` ready |
| Cu (Caliskan 2024) / S (Amarsi 2025) PySME `.grd` | Zenodo 15062813 / 17064337 (separate records) | — | ⚠️ not in the 3982506 set — add to a follow-up fetch if needed |

### A.3 Reference data
| Data | Source | On Sirius |
|---|---|---|
| Frozen solar gold reference `solar_abundances_v1.csv` (+ hash manifest) | RYA-469, committed | 📦 in-repo → Part D |
| canonical_gf.csv (single-source gf, incl. Den Hartog Mn) | committed (15 MB) | 📦 → Part D |
| Built VALD linelists (master + per-system: full, 55cnc, αCenA, αCenB) | committed (~80 MB) | 📦 → Part D |
| **VALD raw deliveries** (account-gated extractions) | VALD3 (Ryan-pulled) | ⚠️ account-gated — stage from Mac (Part D), not a clean URL fetch |
| Kitt Peak NSO solar flux atlas (251 segs) | NSO (RYA-459), EXTERNAL ~44 MB | ⚠️ external download (NSO) — not in-repo; fetch or stage |
| CALSPEC UV composite + IR atlases | RYA-459, committed (3 MB) | 📦 → Part D |
| Den Hartog / Meléndez gf tables | committed (audit dirs) | 📦 → Part D |

### A.4 Constants / config
| `config/constants.py`, `config/stars.yaml`, `config/physics_regime_rya400.yaml`, regions files — committed → 📦 Part D. |

### A.5 Spectra (per audited target — see Part C)
| Target | On Sirius `/mnt/codex-data/spectra/` |
|---|---|
| Solar (Vesta/HARPS) | ⚠️ `sol/` empty — stage (Part D) |
| Procyon | ✅ staged (185 files, 431 MB, `.staged.sha256`) — RYA-420 |
| τ Boo | ✅ staged (125 files, 1.9 GB, `.staged.sha256`) — RYA-420 |
| α Cen A / B | ⚠️ `alpha_cen_a/_b` empty — audited keepers staged in Part D |
| 55 Cnc A | ⚠️ `55cnc_a` empty — audited keepers staged in Part D |

---

## Part B — Downloads (environment-independent, in progress NOW)

Driven from Mac over SSH onto `/mnt/codex-data`, **md5-verified, no silent partials**
(`scripts/sirius_fetch_grids_rya477.py` — size-gate then md5-gate, extract only on PASS,
per-item provenance JSON; the RYA-359/461 stewardship pattern). Manifest:
`data/sirius_manifest/grids_zenodo_3982506.json`.

- **Zenodo 3982506 PySME grid set (13 items, 5.34 GB archives):** ✅ **COMPLETE — 13/13
  verified** (finished 2026-06-29T20:49Z) → `grids/nlte/amarsi_galah/`. Each archive
  md5-checked against the Zenodo checksum before its `.grd` was extracted; archives freed
  post-extract (re-fetchable by md5). **~78 GB extracted** (Al 2.17 · Ba 3.26 · C 7.37 ·
  Ca 9.59 · H 0.28 · K 3.47 · Li 7.06 · Mg 5.78 · **Mn 8.40** · N 5.76 · Na 2.57 · O 4.38 ·
  Si 8.29 GB); Mn `.grd` size bit-identical to the RYA-476 Mac extraction. Status:
  `grids/nlte/_fetch_status_sirius_grids_zenodo_3982506.json`.
- **Engine source:** ✅ iSpec + PySME cloned to `/mnt/codex-data/engines/` (for the RYA-172/375 compile).
- **Committed substrate (Amarsi CNO, MPIA Fe/Ca/Ti/Cr, built VALD linelists, atlases, gold
  ref, constants):** routes via the **Part D repo/stager transfer** — Sirius cannot clone the
  private repo (no `gh`/verified github key; RYA-80). Documented, not silently skipped.
- **VALD raw + Kitt Peak NSO atlas + STAGGER atmospheres:** account-gated / external / blocked
  — see A.2/A.3 notes; not clean URL fetches.

---

## Part C — Benchmark data-audit sweep (one-line per-star verdict)

> The `codex-data-audit` skill named in the ticket is **not available in this
> environment**; the audit was done inline from the committed substrate (RYA-384 / RYA-303
> / RYA-382). Verdicts are real gates, not rubber stamps.

| Star | Audited? | Keeper manifest | One-line verdict |
|---|:--:|---|---|
| **Solar** | ✅ | `data/reference/solar/solar_abundances_v1.csv` (frozen, RYA-469) | Vetted — gold reference frozen. |
| **Procyon** | ✅ | RYA-351; staged `/mnt/codex-data/spectra/procyon/.staged.sha256` | Vetted multi-instrument; on Sirius. |
| **τ Boo** | ✅ | RYA-415; staged `…/tau_boo/.staged.sha256` | Vetted per-arm; on Sirius. |
| **α Cen A** | ✅ | `data/audit/acen_holdings_rya384/spectral_inventory.csv` + `data/audit/alpha_cen_stis_uv/alpha_cen_A_*_files.txt` | **AUDITED** — optical keeper HARPS (88 files, all HD128620, clean); UV keeper = RYA-303 header-truth STIS-A lists (E140M 25 / E140H 3 / E230H 20, vac→air flagged); K-band absent (owed, non-blocking). No blocking gap for an optical+UV run. |
| **α Cen B** | ⚠️ | same dirs, `…_B_*_files.txt` | **AUDITED w/ caveat** — UV keeper = RYA-303 header-truth STIS-B lists (clean); optical primary = CHIRON (76) but **format-aware read confirmation owed (RYA-384)** before optical staging. |
| **55 Cnc A** | ⚠️ | `RYA382_55cnc_audit.md` — **in unmerged worktree rya-382, NOT on main** | **AUDITED (RYA-382)** — optical GO (HARPS 88 / HARPS-N 126 / ESPRESSO 4), UV GO (HST COS G130M+E230M), IR SPIRou GO-caveats; NO-GO = KOA raw (NIRSPEC/KPF/HIRES, some are nph-getKOA stubs); **K-band IR arm blocked on the HFS-on K-band line list (RYA-378/389)**. Optical+UV ready; **the audit doc needs to land on main**. |

**Audit findings worth surfacing:**
1. α Cen STIS-UV is **resolved** by RYA-303 (header-truth dedup → per-star/per-mode keeper
   lists); the RYA-426 INDETERMINATE was folder-name ambiguity, addressable via these lists.
2. The **55 Cnc audit is not committed to main** (lives in the rya-382 worktree) — a
   reproducibility gap to close before Part D.
3. α Cen B optical (CHIRON) needs a format-read confirmation (the one real open audit item).

---

## Part D — Transfer (DEFERRED — gated on RYA-348 + audits green)

Not started (per ticket). Once Procyon/348 is done and the audits above are green:
migrate the committed substrate + audited spectra keepers Mac→Sirius via the RYA-420 stager
(extended beyond Procyon/τ Boo), md5-verify both ends, then free the Mac. **Rescue-before-
clean (RYA-461): nothing deleted from the Mac until confirmed on Sirius + backup (Sirius
drive + Google Drive "Exoplanet Codex Data").**

Pre-transfer checklist: (a) grid fetch 13/13 verified; (b) 55 Cnc audit doc on main;
(c) α Cen B CHIRON format-read confirmed; (d) engine compile (RYA-172) so Sirius can
actually run, not just store.
