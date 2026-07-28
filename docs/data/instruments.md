# Multi-wavelength instruments and data sources

The Codex is a UV–IR program, not a HARPS-only pipeline. Its observational
strategy combines high-resolution space UV, stabilized optical echelles,
ground-based NIR/IR echelles, and solar/reference atlases. The canonical source
is [`data/catalog/instrument_catalog.csv`](../../data/catalog/instrument_catalog.csv);
this page explains how to interpret it. Do not copy this inventory into another
static table.

## Capability, support, and holdings are different

An instrument row answers what the source can provide and how the project treats
it. `codex_status=active` means a scientifically supported project source;
`experimental` means real repository evidence exists but the full science route
is incomplete; `candidate` means capability is registered but holdings/use are
unverified; `rejected` records a deliberate exclusion.

None of those states proves a spectrum exists for a particular star.
Per-system evidence belongs in a holdings manifest:

```text
system_catalog.star_params_key
        ↓
per-system holdings manifest
        ↓ instrument_id
instrument_catalog.instrument_id
```

Element/status questions join through
[`data/audit/element_status_tracker.csv`](../../data/audit/element_status_tracker.csv).

## Wavelength strategy

| Regime | Principal sources | Scientific role | Dominant caveats |
|---|---|---|---|
| FUV/NUV | HST/STIS, HST/COS; UVES blue edge | C/N/O/S/P access and ionization checks | vacuum wavelengths, grating-dependent resolution, scattered light, chromospheric masks |
| Optical/VIS | HARPS, ESPRESSO, UVES; HARPS-N/HIRES/KPF/GHOST candidates | Fe backbone, refractories, C I, [O I], O I 777 | product-specific BERV/frame, iodine/calibration rejection, continuum/blends |
| Red optical/NIR | UVES, NIRPS, SPIRou, Kitt Peak solar atlas | N I, P I, OH/CN and cross-arm checks | tellurics, vacuum-to-air conversion, mode gaps |
| H/K/IR | CRIRES+, SPIRou, iSHELL/NIRSPEC candidates, solar FTS atlases | CO/OH/CN, C/O, isotopes | strong tellurics, per-setting coverage, thermal/blaze structure, custom reduction |
| Reference | Kitt Peak solar FTS, CALSPEC composite, ACE/NSO IR atlases | solar differential anchor and absolute-flux context | CALSPEC solar product is cited composite—not a direct HST solar observation |

Exact coverage and resolving-power fields are in the CSV in nanometers.
Mode/grating-specific metadata in live loaders overrides a broad instrument
envelope when processing a file.

## Archive coverage

- ESO Science Archive: HARPS, ESPRESSO, UVES, FEROS, CRIRES+, NIRPS, FLAMES.
- MAST: HST/STIS, HST/COS, CALSPEC/reference products.
- TNG/IA2: HARPS-N.
- Keck Observatory Archive: HIRES, NIRSPEC, KPF.
- CADC and PolarBase: SPIRou, ESPaDOnS, NARVAL.
- Gemini Observatory Archive: GHOST and GNIRS.
- NOIRLab/SMARTS: Phoenix and CHIRON evidence.
- IRTF/IRSA: iSHELL and IRTF reference spectra.
- NSO/Kitt Peak: disk-integrated solar FTS atlases.

The catalog records archive capability. Check the System Register and a
target-specific manifest before claiming verified holdings.

## Preprocessing rules

- Tellurics: NIR/IR products require a declared correction or forward model.
  SPIRou ingestion accepts only APERO telluric-corrected `_t.fits`; CRIRES+
  remains per-night/molecfit gated.
- Wavelength frames: HST and most IR sources are native vacuum. The pipeline's
  optical convention is air at/above 200 nm; conversion happens once at the
  loader boundary. Never infer units from array magnitude when a DRS contract is
  available.
- Velocity frames: inspect `SPECSYS`, BERV keywords, and the instrument contract.
  NIRPS products are already corrected in the verified path; UVES/CRIRES products
  may require explicit correction. Double-BERV is a hard failure.
- Iodine/calibration: reject iodine-cell spectra for abundance work unless an
  explicit contamination model is part of the method. Never ingest calibration,
  stacked velocity-smeared, polarimetric, or raw products as if they were a
  normalized stellar S1D.
- Custom reduction: raw Phoenix, HIRES, NIRSPEC, MIKE, and some survey/archive
  holdings are not science-ready merely because they are downloadable.

## Minimal and full reproduction

The repository-native quick start uses a vendored correction grid and no
instrument data. The smallest real-spectra benchmark path uses one approved
HARPS solar or Procyon product plus its manifest and the full synthesis assets.
A full multi-wavelength reproduction adds, as applicable, STIS/COS UV, UVES or
ESPRESSO optical/red, NIRPS/SPIRou NIR, CRIRES+ K-band, and the solar atlases.
Each arm is independently conditioned and validated before cross-arm comparison;
values from disagreeing arms are flagged, not silently averaged.

## Validation and update discipline

```bash
python -m data.catalog.instruments
python -m pytest tests/test_instrument_catalog.py -q
```

Update the CSV first, cite the authoritative instrument/DRS documentation, add
or adjust tests, then update prose. Candidate/rejected rows should remain in the
register so future agents do not repeat the same archive evaluation.
