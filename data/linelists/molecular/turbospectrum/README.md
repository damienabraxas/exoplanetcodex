# Vendored C/N/O molecular line lists (RYA-360)

The git-tracked **secure record of truth** for the Turbospectrum molecular line lists
the C/N/O synthesis keystone (RYA-237) depends on. Machine-readable provenance +
per-molecule baselines live in [`MOLECULAR_MANIFEST.json`](MOLECULAR_MANIFEST.json);
this file is the prose.

## Reconciliation (the RYA-237 false-absence)

RYA-237's recon reported "no CH/CN Turbospectrum lists" — it looked in the
**raw-download dir** `data/linelists/molecular/exomol/` (the `.bz2` archives). The
usable Turbospectrum `.bsyn` lists actually live in the **iSpec tool dir**
`ispec/input/linelists/turbospectrum/molecules/` (the RYA-236 audit's authoritative
location). They were present and verified all along; the recon looked in the wrong
place. RYA-360 re-verified at the correct path and reproduced the RYA-236 HARPS-window
counts **exactly**:

| diagnostic | window (Å) | isotopologue | lines | RYA-236 |
|---|---|---|---|---|
| CH A-X G-band | 4290–4315 | 12CH | 583 | ≈583 ✓ |
| CN red A-X | 6000–6200 | 12C14N | 3534 | ≈3534 ✓ |
| C2 Swan (0,0) | 5130–5170 | 12C12C | 1019 | ≈1019 ✓ |

## The problem this closes — security, not absence

The project repo tracked **zero** `.bsyn` / molecular artifacts; they lived only in the
iSpec **install** tree. An iSpec reinstall/rebuild/relocation would silently reset them
to the stock bundle and **wipe the RYA-236 CO addition** (`CO_IR_Li2015.dat`), with
nothing to catch it — the same failure class as the gf / blend_flag / STAR_PARAMS arc.
These vendored copies + the RYA-355 `[molecular]` stewardship invariant make it
mechanical: a missing/emptied list, or an iSpec dir that no longer matches this vendored
record, is now a **loud CI failure**, not a silent RYA-237-style false-absence.

## What is vendored

`<MOL>/` per molecule (isotopologues + wavelength-region chunks, iSpec filenames kept):

- **CH** — `12CH`, `13CH` (A-X; Masseron et al. 2014/2022)
- **CN** — `12C14N`, `13C14N`, `12C15N` (red A-X; Brooke/Sneden et al. 2014, isotopologues)
- **C2** — `12C12C`, `12C13C`, `13C13C` (Swan)
- **OH** — `16OH` (A-X electronic)
- **NH** — `14NH` (A-X electronic)
- **CO** — `CO_IR_Li2015.dat` (ExoMol Li et al. 2015; RYA-236 conversion to Turbospectrum)

Distribution/compilation of the `.bsyn` set: the Masseron VALD-based compilation shipped
with iSpec/Turbospectrum (Gerber et al. 2023, A&A 669, A43; uploaded 2022 to
`keeper.mpdl.mpg.de`, split per wavelength region — see the molecules-dir `README.md`).

## Band regime — electronic vs mid-IR (RYA-499 / RYA-360, measured not inferred)

The `.bsyn` lists (CH/CN/C2/OH/NH) span ≈4200–9200 Å — **400–950 nm ELECTRONIC bands
only**. Measured row counts in the mid-IR ro-vibrational **fundamental** windows are
**zero** for all of OH/NH/CH:

| species | fundamental window | Å window | rows | verdict |
|---|---|---|---|---|
| OH | 2600–3600 cm⁻¹ | 27778–38462 | 0 | mid-IR ABSENT (RYA-503 unblocked) |
| NH | 3000–3500 cm⁻¹ | 28571–33333 | 0 | mid-IR ABSENT (RYA-503 unblocked) |
| CH | 2650–3100 cm⁻¹ | 32258–37736 | 0 | mid-IR ABSENT (RYA-503 unblocked) |

**CO is the sole mid-IR exception**: `CO_IR_Li2015.dat` spans ≈4515–99958 Å — NIR
overtones through the ~4.6 µm fundamental (a genuine mid-IR ro-vibrational list).

## Pipeline resolution

Turbospectrum (iSpec, `use_molecules=True`) reads from the iSpec molecules dir; that
path is hard-coded in iSpec, so repointing it mid-stream is risky. The vendored copy is
therefore the **tracked backup-of-record**: the `[molecular]` guard gates on the vendored
copy (present + non-empty + provenance) and, when iSpec is present, additionally asserts
the iSpec dir still matches it (drift = a reinstall reset/wiped a list). Re-vendor with
`scripts/vendor_molecular_lists_rya360.py`.
