# Model and correction assets

The live synthesis code consumes iSpec-formatted packs. Downloading a similarly
named upstream grid is not sufficient unless it is converted to the layout iSpec
expects. `python scripts/validate_installation.py --full` validates paths and
executables; an interpolation test in `tests/test_abundances_derive.py` validates
usable atmosphere packs when iSpec is available.

## Asset matrix

| Asset | Classification | Live purpose | Expected path |
|---|---|---|---|
| `ATLAS9.Castelli` | Mandatory for FGK full runs | 1D LTE plane-parallel atmosphere | `$ISPEC_DIR/input/atmospheres/ATLAS9.Castelli/` |
| `MARCS.GES` | Target-specific/experimental | Cool-star/M-dwarf atmosphere | `$ISPEC_DIR/input/atmospheres/MARCS.GES/` |
| Turbospectrum executable/data | Mandatory for synthesis | Radiative transfer, opacity/EOS | `$ISPEC_DIR/synthesizer/turbospectrum/` |
| iSpec atomic/molecular inputs | Mandatory for synthesis | Lines, isotopes, molecules | `$ISPEC_DIR/input/` |
| MPIA/INSPECT/CDS NLTE CSVs | Vendored; mandatory when registry applies them | Fe, Ca, Ti, Cr, Na, Mg, Ba, Mn, Si corrections | `data/nlte_grids/` |
| Amarsi C/O tables | Vendored; mandatory for C/O correction path | 3D-NLTE/1D-NLTE C I and O I corrections | `data/nlte_grids/amarsi2019_cno/` |
| `vendor/1L-3NErrors` | Archived/optional | Fe 3D-NLTE cross-check, not live Fe source | `vendor/1L-3NErrors/` |
| Standalone solar atmosphere | Not used | No separate live asset | — |

## Atmosphere and radiative-transfer sources

ATLAS9.Castelli is the FGK choice hardcoded in
`pipeline/abundances_derive.py`. The scientific source is Castelli & Kurucz
(2004); upstream grids are available from
[Kurucz](http://kurucz.harvard.edu/grids.html) and the
[STScI Castelli/Kurucz page](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/castelli-and-kurucz-atlas).
The exact runtime artifact is the atmosphere pack distributed for the recorded
iSpec commit, not an arbitrary raw Kurucz directory.

MARCS.GES is selected for cool stars and evolved stars. MARCS describes 1D,
hydrostatic LTE atmospheres; cite Gustafsson et al. (2008). Source and usage
information are on the [MARCS site](https://marcs.astro.uu.se/) and its
[download page](https://marcs.astro.uu.se/data.html). The reference iSpec pack
is about 39 MB; the upstream archive sizes differ. M-dwarf reproduction remains
experimental even when interpolation passes.

Turbospectrum is the live synthesis radiative-transfer code
(`RADIATIVE_TRANSFER_CODE='turbospectrum'`). It is shipped/compiled through the
tested iSpec tree. Cite Alvarez & Plez (1998), Plez (2012), and the versioned
[Turbospectrum source](https://github.com/bertrandplez/Turbospectrum2019).
Its continuous opacity, equation-of-state, isotope, and molecular assets must
match the executable/iSpec bundle. Do not mix bundles without a new validation.

### Version and checksum policy

The repository does not redistribute the iSpec input archive and currently has
no upstream archive checksum. Record the iSpec Git commit, input-archive name,
download date/URL, byte size, and locally computed SHA-256:

```bash
git -C "$ISPEC_DIR" rev-parse HEAD
find "$ISPEC_DIR/input/atmospheres/ATLAS9.Castelli" -type f -print0 \
  | sort -z | xargs -0 shasum -a 256 > atlas9-files.sha256
find "$ISPEC_DIR/input/atmospheres/MARCS.GES" -type f -print0 \
  | sort -z | xargs -0 shasum -a 256 > marcs-files.sha256
python scripts/validate_installation.py --full
python -m pytest tests/test_abundances_derive.py -q
```

Licenses and redistribution rights are upstream-controlled. A downloadable
asset is not automatically redistributable; retain upstream notices and cite it.

## NLTE and 3D corrections

The live Fe correction source is `Fe_Bergemann_MPIA.csv`: a 1D-NLTE
MPIA/SpectrumTools grid interpolated by line wavelength and
`(Teff, logg, [Fe/H])`. A query outside the convex hull or without a wavelength
match returns unavailable; it must not silently become LTE. The archived Amarsi
MLP under `vendor/1L-3NErrors` is only a 3D-versus-1D cross-check.

`NLTE_CORRECTION_ELEMENTS` in `config/constants.py` is authoritative:

- Ca I: Mashonkina et al. (2017), MPIA MAFAGS-OS.
- Ti I: Bergemann (2011), MPIA MAFAGS-OS.
- Cr I: Bergemann & Cescutti (2010), MPIA MAFAGS-OS.
- Na I: Lind et al. (2011), INSPECT.
- Mg I, Mn I, Si I: MPIA/SpectrumTools sources recorded in `.prov.json`.
- Ba II: Korotin et al. (2015), CDS/VizieR.

Every grid is line-specific. Provenance sidecars provide source URLs where
available; builder/fetch scripts under `scripts/` document reconstruction.
Missing/out-of-grid corrections return NaN or a flagged result; they are not
permission to label an LTE result “NLTE”.

For C I/O I, `pipeline/nlte_cno.py` uses Amarsi, Nissen & Skúladóttir (2019):
the 3D-NLTE minus 1D-LTE tables at/below 6500 K and the 1D-NLTE leg above that
ceiling. Parameter interpolation is over Teff, log g, metallicity, microturbulence,
and abundance. The tables, CDS ReadMe, and
`provenance.json` are vendored. C/O corrections are automated only for supported
lines and in-grid parameters.

The constant `CORRECTIONS_3D` contains literature offsets for forbidden oxygen
lines, but not every older design note describes a currently executed path.
Treat hardcoded/manual corrections as method-specific and verify their call site
before reporting them as applied. No general 3D atmosphere synthesis grid is
currently operational.
