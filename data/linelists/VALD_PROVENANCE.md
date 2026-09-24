# VALD raw-extraction provenance (RYA-1228)

**Why this file exists.** `data/linelists/vald_*_raw.txt` are the only LFS-tracked data in
this repository: **53 LFS objects, 853.9 MiB** across all refs. One further raw extraction
(**`vald_55cnc_raw.txt` @ `a7d016acd5…`, 7.28 MiB**) was committed as an ORDINARY git blob
before this repository adopted LFS, bringing the census to **54 objects, 861.2 MiB**.

🔴 **A plain blob is invisible to `git lfs ls-files --all`.** The Phase 1 inventory was
built from that command and therefore missed it — it would have been purged with no
preserved copy, the one outcome this ticket calls CRITICAL. The census is now taken by
PATH (`git rev-list --all --objects -- data/linelists/vald_*_raw.txt`) and each object is
then classified by storage, which is the only enumeration that cannot miss this class.

⚠️ **The GitHub LFS quota is 43 objects / 692.7 MiB, not 53 / 853.9 MiB.** Ten objects
(161.2 MiB, including all five tau Boo extractions) exist only on local branches that were
never pushed, so they never consumed quota. The plain blob is ordinary pack, not quota. They are raw VALD3 *Extract Stellar* deliveries, not Codex measurements.

🔴 **NONE OF THESE IS RE-FETCHABLE.** VALD3 is a manual web extraction with no API. Losing
one means a manual re-extraction per the `codex-vald-extraction` skill, against a service
whose output can change. Every object below is preserved OUTSIDE git with its sha256 before
anything is removed from history (RYA-461 archive policy).

**What this file guarantees.** After the raw deliveries leave the repository, the built line
lists (`linelist_*.csv`, `canonical_gf.csv`) remain reproducible: each row below names the
exact delivery, its content hash, its wavelength span and the archive path holding it.
The built lists are plain committed files and are NOT affected by this ticket.

## Archive

```
~/Documents/Exoplanet Codex/vald_raw_archive/
  current/      23 objects, 440.2 MiB   -- the version each path resolves to at HEAD
  superseded/   30 objects, 413.7 MiB   -- history-only versions, named <stem>.<oid10>.txt
  MANIFEST.json         full record, one entry per object
  MANIFEST.sha256       `shasum -a 256 -c` format
```

⚠️ **The LFS oid IS the sha256 of the file content** — verified directly on two objects
including a 51 MiB one — so `archive filename hash == oid` is an exact identity check, not
a re-derivation. The plain blob is hashed from its git content directly. All copies were
verified this way: **54/54 pass, 0 failures** (53 LFS + 1 plain blob).

## CURRENT — the version a build resolves today (23 objects, 440.2 MiB)

| file | MiB | requested Å | actual Å | HFS | elts | added | sha256 (LFS oid) |
|---|---:|---|---|---|---:|---|---|
| `vald_55cnc_nir_raw.txt` | 2.1 | 6910.0–17000.0 | 6911.5–16997.2 | unknown | 25.0 | 2026-06-07 | `f7c76015b3331cdf…` |
| `vald_55cnc_raw.txt` | 51.0 | 3780.0–6910.0 | 3780.0–6909.8 | unknown | 71.0 | 2026-05-30 | `3db26f88ea3a2b94…` |
| `vald_55cnc_uv_a_raw.txt` | 9.9 | 1150.0–2000.0 | 1150.3–2000.0 | unknown | 46.0 | 2026-06-13 | `125bcf41756411ab…` |
| `vald_55cnc_uv_b_raw.txt` | 36.6 | 2000.0–3780.0 | 2000.1–3780.0 | unknown | 72.0 | 2026-06-13 | `8ecb6a3c852ba1af…` |
| `vald_55cnc_uv_raw.txt` | 31.4 | 1150.0–3780.0 | 1150.3–3780.0 | unknown | 77.0 | 2026-06-12 | `1b5b17204ce3df1d…` |
| `vald_alpha_cen_a_nir_raw.txt` | 1.5 | 6910.0–17000.0 | 6911.5–16997.2 | unknown | 24.0 | 2026-06-14 | `d865563ea3195a7a…` |
| `vald_alpha_cen_a_optical_raw.txt` | 11.3 | 3780.0–6910.0 | 3780.3–6907.7 | unknown | 53.0 | 2026-06-14 | `109978735a86578d…` |
| `vald_alpha_cen_a_uv1_raw.txt` | 11.3 | 1150.0–2000.0 | 1150.3–2000.0 | unknown | 48.0 | 2026-06-14 | `37568699d4fd00bb…` |
| `vald_alpha_cen_a_uv2_raw.txt` | 33.4 | 2000.0–3780.0 | 2000.1–3779.7 | unknown | 72.0 | 2026-06-14 | `f5282dd04047e421…` |
| `vald_alpha_cen_b_nir_raw.txt` | 2.1 | 6910.0–17000.0 | 6911.5–16997.2 | unknown | 25.0 | 2026-06-14 | `c8a6051daade6317…` |
| `vald_alpha_cen_b_optical_raw.txt` | 15.0 | 3780.0–6910.0 | 3780.3–6907.7 | unknown | 55.0 | 2026-06-14 | `eed3ba39471e80c3…` |
| `vald_alpha_cen_b_uv1_raw.txt` | 9.8 | 1150.0–2000.0 | 1150.3–2000.0 | unknown | 46.0 | 2026-06-14 | `ed03c7263e64ac3f…` |
| `vald_alpha_cen_b_uv2_raw.txt` | 36.3 | 2000.0–3780.0 | 2000.1–3780.0 | unknown | 72.0 | 2026-06-14 | `68c30fd8a917c686…` |
| `vald_procyon_nir_hfson_raw.txt` | 0.5 | 6910.0–11000.0 | 6914.5–10984.5 | on | 22.0 | 2026-06-14 | `b3008e27601b2ba2…` |
| `vald_procyon_optical_hfson_raw.txt` | 15.4 | 3000.0–6910.0 | 3000.1–6907.7 | on | 64.0 | 2026-06-14 | `1f5acca42cbe14f5…` |
| `vald_procyon_uv_hfson_raw.txt` | 28.2 | 1150.0–3000.0 | 1150.3–3000.0 | on | 67.0 | 2026-06-14 | `cbaf699bb29c7080…` |
| `vald_solar_fuv_1150_2000_hfson_raw.txt` | 17.3 | 1150.0–2000.0 | 1150.3–2000.0 | on | 51.0 | 2026-06-20 | `e0341f3247bb84d1…` |
| `vald_solar_ir_17000_25000_hfson_raw.txt` | 4.0 | 17000.0–25000.0 | 17000.6–24998.0 | on | 27.0 | 2026-06-20 | `5b39cfa9058cfc94…` |
| `vald_solar_ir_6910_14000_raw.txt` | 1.0 | 6910.0–9500.0 | 6911.5–13997.7 | unknown | 24.0 | 2026-08-11 | `6097a36e67844453…` |
| `vald_solar_ir_9500_17000_hfson_raw.txt` | 6.8 | 9500.0–17000.0 | 9500.3–16998.5 | on | 40.0 | 2026-06-20 | `57d31b77176c2628…` |
| `vald_solar_nearuv_2000_3780_hfson_raw.txt` | 65.3 | 2000.0–3780.0 | 2000.0–3780.0 | on | 76.0 | 2026-06-20 | `a39ffc7b844cbac2…` |
| `vald_solar_raw.txt` | 44.3 | 3780.0–6910.0 | 3780.0–6909.8 | unknown | 73.0 | 2026-05-30 | `5cba963c5cc89113…` |
| `vald_solar_redopt_6910_9500_hfson_raw.txt` | 5.8 | 6910.0–9500.0 | 6910.2–9499.9 | on | 45.0 | 2026-06-20 | `1a834a7616c9778a…` |

## SUPERSEDED — present only in history (30 objects, 413.7 MiB)

Three distinct reasons, each measured rather than assumed:

| reason | n | MiB |
|---|---:|---:|
| older version of a path that still exists at HEAD | 15 | 232.0 |
| star not in the roster (tau Boo — no built line list) | 5 | 136.1 |
| deleted narrow sub-range, superseded by a consolidated extraction | 10 | 45.5 |

⚠️ **A citation names a PATH, not a version.** Where a superseded object shares a filename
with a live one, the code references resolve to the HEAD object; the historical version is
referenced by nothing. Reference *counts* alone would wrongly credit it.

| file (archived as) | MiB | reason | sha256 (LFS oid) |
|---|---:|---|---|
| `vald_tauboo_uv_b_raw.3c1912ae6f.txt` | 62.8 | star 'tauboo' has no built line list in data/linelists/ (r | `3c1912ae6f112cf2…` |
| `vald_tauboo_optical_raw.6320c92da4.txt` | 41.6 | star 'tauboo' has no built line list in data/linelists/ (r | `6320c92da4e87abc…` |
| `vald_55cnc_uv_b_raw.1abf9e78fb.txt` | 38.1 | older version of a path that still exists at HEAD -- the c | `1abf9e78fb38c0c4…` |
| `vald_alpha_cen_b_uv2_raw.ed479edced.txt` | 37.0 | older version of a path that still exists at HEAD -- the c | `ed479edceda02d65…` |
| `vald_alpha_cen_a_uv2_raw.a208050b2f.txt` | 36.4 | older version of a path that still exists at HEAD -- the c | `a208050b2fea9083…` |
| `vald_solar_nearuv_2000_3780_hfson_raw.2575fdc2a7.txt` | 33.8 | older version of a path that still exists at HEAD -- the c | `2575fdc2a7368b50…` |
| `vald_55cnc_optical_raw.5e0c8ce344.txt` | 19.0 | deleted from the tree | `5e0c8ce344d258d0…` |
| `vald_tauboo_uv_a_raw.426f9aa149.txt` | 18.3 | star 'tauboo' has no built line list in data/linelists/ (r | `426f9aa1493ecb12…` |
| `vald_alpha_cen_b_optical_raw.a6818b4671.txt` | 17.2 | older version of a path that still exists at HEAD -- the c | `a6818b4671190156…` |
| `vald_alpha_cen_a_optical_raw.eb959e84c0.txt` | 13.9 | older version of a path that still exists at HEAD -- the c | `eb959e84c084e103…` |
| `vald_solar_fuv_1150_2000_hfson_raw.7a95064a02.txt` | 11.4 | older version of a path that still exists at HEAD -- the c | `7a95064a02116243…` |
| `vald_tauboo_nir_raw.90f64669dc.txt` | 10.5 | star 'tauboo' has no built line list in data/linelists/ (r | `90f64669dcc305c5…` |
| `vald_alpha_cen_a_uv1_raw.1dd15f0978.txt` | 10.4 | older version of a path that still exists at HEAD -- the c | `1dd15f09781d6adb…` |
| `vald_55cnc_nir_raw.03d171ee11.txt` | 8.7 | older version of a path that still exists at HEAD -- the c | `03d171ee11361a77…` |
| `vald_alpha_cen_b_uv1_raw.950be5adc0.txt` | 8.6 | older version of a path that still exists at HEAD -- the c | `950be5adc0025221…` |
| `vald_55cnc_uv_a_raw.d4b8448fec.txt` | 8.2 | older version of a path that still exists at HEAD -- the c | `d4b8448fec0aedf1…` |
| `vald_55cnc_ir_17000_25000_hfson_raw.4ab48723da.txt` | 7.5 | deleted from the tree | `4ab48723da14dffa…` |
| `vald_alpha_cen_b_ir_17000_24000_raw.afc4ef94c8.txt` | 6.0 | deleted from the tree | `afc4ef94c854499a…` |
| `vald_alpha_cen_a_ir_17000_24000_raw.24b6456a80.txt` | 5.0 | deleted from the tree | `24b6456a809afbec…` |
| `vald_procyon_nir_9500_17000_raw.a6f96433c3.txt` | 4.3 | deleted from the tree | `a6f96433c3fc4b9d…` |
| `vald_55cnc_nir_raw.0260721b64.txt` | 3.1 | older version of a path that still exists at HEAD -- the c | `0260721b6487ae43…` |
| `vald_tauboo_ir_17000_24000_raw.4c3ce9740d.txt` | 3.0 | star 'tauboo' has no built line list in data/linelists/ (r | `4c3ce9740da4c8f2…` |
| `vald_alpha_cen_b_nir_raw.8bd0902bbb.txt` | 2.7 | older version of a path that still exists at HEAD -- the c | `8bd0902bbb95b6a7…` |
| `vald_procyon_ir_17000_24000_raw.a3a9d4b06d.txt` | 2.2 | deleted from the tree | `a3a9d4b06d8ee50c…` |
| `vald_alpha_cen_a_nir_raw.258df92c63.txt` | 2.0 | older version of a path that still exists at HEAD -- the c | `258df92c6304a689…` |
| `vald_solar_ir_14000_18000_hfson_raw.7236a2fc71.txt` | 0.7 | deleted from the tree | `7236a2fc7176935f…` |
| `vald_solar_redopt_6910_9500_hfson_raw.54db11875e.txt` | 0.6 | older version of a path that still exists at HEAD -- the c | `54db11875ee4f367…` |
| `vald_solar_ir_18000_25000_hfson_raw.7dca0c6f71.txt` | 0.4 | deleted from the tree | `7dca0c6f71a97855…` |
| `vald_solar_ir_9500_14000_hfson_raw.e7fbd9f665.txt` | 0.4 | deleted from the tree | `e7fbd9f665949d62…` |
| `vald_55cnc_supp_019572_raw.c99e147bf5.txt` | 0.0 | deleted from the tree | `c99e147bf569cca2…` |

## Not in scope

* `data/linelists/linelist_*.csv`, `canonical_gf.csv` — Codex-built products, plain
  committed files, never LFS. Untouched.
* `data/linelists/vald_procyon_raw/*.vald` — raw VALD deliveries that are NOT LFS-tracked
  (the `.gitattributes` filter covers `vald_*_raw.txt` only), so they cost no LFS quota.
  They are still raw external data and fall under the Phase 2 standing policy.

## Status

Phase 1 (this file): inventory + preservation only. **Nothing has been removed from history.**
Phase 2 — the history rewrite, the `.gitattributes`/`.gitignore` changes and the pre-commit
guard — awaits Ryan's approval of the Phase 1 report, and Ryan performs the force-push.

## PRE-LFS — committed as an ordinary git blob (1 object, 7.3 MiB)

Not an LFS object, so absent from every `git lfs ls-files` listing; destroyed by the same
history rewrite all the same. Preserved at `superseded/vald_55cnc_raw.a7d016acd5.txt`.

| file | MiB | storage | sha256 (of content) |
|---|---:|---|---|
| `vald_55cnc_raw.txt` | 7.28 | plain git blob (pre-LFS) | `a7d016acd567dbca…` |

Two further plain blobs exist on these paths (`vald_55cnc_raw.txt` @ `3db26f88ea3a…` and
`vald_solar_raw.txt` @ `5cba963c5cc8…`); both are byte-identical to an LFS object already
listed above, so they are covered by an existing preserved copy and are not re-listed.
