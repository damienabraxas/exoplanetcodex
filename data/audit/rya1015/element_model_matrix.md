# Element x model availability (RYA-1015)

Generated 2026-08-23 by `pipeline/model_availability_matrix.py (RYA-1015)`.
Disk half: `data/audit/rya1015/sirius_scan_raw.txt` (Sirius `find -L`, control PASS).

| element | 1D_LTE | 1D_NLTE | MEAN3D_NLTE | FULL_3D_NLTE |
|---|---|---|---|---|
| Li | HAVE | **DISK_ONLY** | REQUEST_ONLY | NONE |
| C | HAVE | HAVE | REQUEST_ONLY | NONE |
| N | HAVE | HAVE | REQUEST_ONLY | NONE |
| O | HAVE | HAVE | REQUEST_ONLY | NONE |
| Na | HAVE | HAVE | REQUEST_ONLY | NONE |
| Mg | HAVE | HAVE | REQUEST_ONLY | NONE |
| Al | HAVE | HAVE | **DISK_ONLY** | NONE |
| Si | HAVE | HAVE | REQUEST_ONLY | NONE |
| P | HAVE | NONE | NONE | NONE |
| S | HAVE | HAVE | REQUEST_ONLY | NONE |
| K | HAVE | HAVE | REQUEST_ONLY | NONE |
| Ca | HAVE | HAVE | REQUEST_ONLY | NONE |
| Sc | HAVE | NONE | REQUEST_ONLY | NONE |
| Ti | HAVE | HAVE | REQUEST_ONLY | NONE |
| V | HAVE | NONE | NONE | NONE |
| Cr | HAVE | HAVE | **DISK_ONLY** | NONE |
| Mn | HAVE | HAVE | NONE | NONE |
| Fe | HAVE | HAVE | REQUEST_ONLY | NONE |
| Fe II | HAVE | HAVE | REQUEST_ONLY | NONE |
| Co | HAVE | NONE | REQUEST_ONLY | NONE |
| Ni | HAVE | **DISK_ONLY** | NONE | NONE |
| Cu | HAVE | HAVE | REQUEST_ONLY | NONE |
| Zn | HAVE | NONE | REQUEST_ONLY | NONE |
| Sr | HAVE | HAVE | REQUEST_ONLY | NONE |
| Y | HAVE | **DISK_ONLY** | **DISK_ONLY** | NONE |
| Zr | HAVE | NONE | NONE | NONE |
| Ba | HAVE | HAVE | REQUEST_ONLY | NONE |
| Eu | HAVE | **DISK_ONLY** | **DISK_ONLY** | NONE |

**PROBLEM cells: 8**

- `Li` / `1D_NLTE` -> **DISK_ONLY**: 1 departure grid(s) on Sirius but NO code entry in NLTE_CORRECTION_ELEMENTS -- unwired.
- `Al` / `MEAN3D_NLTE` -> **DISK_ONLY**: <3D> STAGGERmean3D deck present on Sirius (1) but no code path consumes a <3D> departure deck -- unwired capability.
- `Cr` / `MEAN3D_NLTE` -> **DISK_ONLY**: <3D> STAGGERmean3D deck present on Sirius (1) but no code path consumes a <3D> departure deck -- unwired capability.
- `Ni` / `1D_NLTE` -> **DISK_ONLY**: 1 departure grid(s) on Sirius but NO code entry in NLTE_CORRECTION_ELEMENTS -- unwired.
- `Y` / `1D_NLTE` -> **DISK_ONLY**: 1 departure grid(s) on Sirius but NO code entry in NLTE_CORRECTION_ELEMENTS -- unwired.
- `Y` / `MEAN3D_NLTE` -> **DISK_ONLY**: <3D> STAGGERmean3D deck present on Sirius (1) but no code path consumes a <3D> departure deck -- unwired capability.
- `Eu` / `1D_NLTE` -> **DISK_ONLY**: 1 departure grid(s) on Sirius but NO code entry in NLTE_CORRECTION_ELEMENTS -- unwired.
- `Eu` / `MEAN3D_NLTE` -> **DISK_ONLY**: <3D> STAGGERmean3D deck present on Sirius (1) but no code path consumes a <3D> departure deck -- unwired capability.
