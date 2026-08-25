# Fe ⟨3D⟩ Gerber deck — aux tables (RYA-1035)

The two aux tables shipped with `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin`, committed
here **as evidence and as a test fixture**, not as a pipeline input. The deck binary itself
(92.9 MB unpacked) stays Sirius-only like every other Gerber grid; these two files are 19 KB
and 31 KB and are the whole reason the route needs a decision rather than a registry line.

| file | vendor name convention | usable? |
| --- | --- | --- |
| `auxData_Fe_STAGGERmean3D_May-21-2021.txt` | STAGGER (`p5777g44m00`) | ✅ **register against this one** |
| `auxData_Fe_STAGGERmean3D_May-21-2021_marcs_names.txt` | MARCS-style alias | ❌ unrecoverable — see below |

## 🔴 The vendor's `[Fe/H]` column is zeroed on the seven Teff = 5777 rows

The deck holds 189 atmosphere nodes. Seven of them are the solar-Teff sequence
`p5777g44m00 / m05 / m10 / m20 / m30 / m40 / p05` — the full metallicity axis at Teff 5777 —
but the file's `[Fe/H]` column reads **+0.00 for all seven**. The other 182 rows are correct:
name and column agree exactly, 182/182. So the column is wrong and the name is right.

Left alone this is not a crash, it is a **wrong star**. All seven tie at the solar node,
the tie-break is on A(X) — identical at 7.50 across all seven — and the first wins, which is
`p5777g44m10`, **[Fe/H] = −1.0**. The true solar record `p5777g44m00` is sixth in file order
and unreachable.

`pipeline.gerber_nlte._parse_aux_text` refereeds the column with the model name and records
the override (`feh_aux`, `feh_from_name`), which restores all six and resolves the solar node
to exactly `p5777g44m00`.

## 🔴 Why the `_marcs_names` variant cannot be used — even though it is the canonical one

TSFitPy's own downloader (`utilities/nlte_grids_links.cfg`, `[Fe] 3d_aux_link`) points at the
`_marcs_names.txt` file, and the Al ⟨3D⟩ deck is registered against its `_marcs_names` aux.
For Fe that file is **unrecoverable**: the vendor's `convert_3d_grid_to_marcs_names.py` builds
the new name *from the `[Fe/H]` column*, so it propagated the zeroing into the name. All seven
rows come out as the byte-identical string
`p5777_g+4.4_m0.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00`.

Name and column now agree — and both are wrong for six of the seven. There is no signal left
to referee with, so the name check finds nothing to correct and the reader falls back to the
7-way tie. Measured end-to-end: reading the solar node through the `_marcs_names` aux is
**refused** by the RYA-821 record-vs-aux check; through the plain aux it returns
`p5777g44m00`, ndep 101, nlev 607, log τ −5.000…+5.000, all finite.

**⇒ `Fe@mean3D` must be registered against the PLAIN aux — the opposite of `Al@mean3D`.**

## 🔴 Fe is not alone — Mn has the identical defect (swept 2026-08-25)

`scripts/rya1035_mean3d_aux_defect_sweep.py` scored all **17** published ⟨3D⟩ aux tables
with `gerber_nlte`'s own parser. **Two are defective, and both are May-2021 vintage:**

| deck | rows | `[Fe/H]` overridden | where | `_marcs_names` |
| --- | --- | --- | --- | --- |
| **Fe** (May-21-2021) | 189 | **6** | (5777, 4.44) | **collapsed** 7 rows → 1 name |
| **Mn** (May-17-2021) | 4769 | **150** | (5777, 4.44) | **collapsed** 175 rows → 1 name |
| the other 15 | — | 0 | — | faithful |

**The same defect, differently sized.** Six metallicities are wrong on each; Fe resolves
ONE abundance per node so that is 6 rows, Mn resolves 25 so it is 150. Counting rows rather
than nodes would make Mn look like a worse problem and Fe like a rounding error.

**Al is the positive control**: its conversion keeps all 7 names at Teff=5777, so the
collapse belongs to the defective input rather than to the converter in general.

⇒ **Fe and Mn must be registered against the PLAIN aux. The other 15 may use either**, and
`Al@mean3D` correctly uses `_marcs_names`. The per-element verdict is a committed row in
`data/results/rya1035/mean3d_aux_defect_sweep.csv` — read it before writing a registry line.

## Provenance

- Source: MPG Keeper (Seafile) share `https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/`,
  path `/dep-grids/Fe/`. Fetched 2026-08-24 under RYA-1035.
- Citation: Gerber, Bergemann et al. 2023, A&A 669, A43 (arXiv 2206.00967).
- md5 — plain `302769ca102530e9a209e1887df8c3cd`, marcs_names `4ceb02f997dff488514f229f6422bbc6`.
- No DOI; the Seafile share is mutable, hence the md5 pin (RYA-540 convention).
- The full file set, md5s and the canonical per-file links are in the sibling
  `../Fe_gerber2023.prov.json` under `files.grid_3d_*`.
