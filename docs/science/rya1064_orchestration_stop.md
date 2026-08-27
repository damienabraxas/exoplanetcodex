# RYA-1064 — Alpha Cen A HARPS/VIS Fe orchestration experiment

Date: 2026-08-26 (America/Denver)  
Repository base: `origin/main` at `c9a76ec`  
Disposition: **SCIENTIFIC STOP — no abundance execution permitted**

## Requested run

`target=alpha_cen_a | instrument=harps | region=VIS | element=Fe | ion=I | line_tier=graded | engines=all_valid`

The requested science result was not run. This is the hard-stop outcome required by
RYA-1064: the repository cannot prove a current, approved, continuum-conditioned Alpha
Cen A HARPS product, and the measurement harness cannot address the registered holding.
No legacy spectrum, alternate instrument, ungraded line, solar-calibrated gf, or engine
substitute was used.

## Read-first reconstruction

The repository state register, target registry, element/method standards, holdings
registry, `STAR_PARAMS`, run descriptor, measurement-harness holding table, and current
Fe graded-gf implementation were inspected. Linear lineage inspected: RYA-609, RYA-384,
RYA-323, RYA-302, RYA-308, RYA-824, RYA-836, RYA-850, RYA-945, RYA-435, and RYA-1063.

Relevant current state:

- `config/stars.yaml` is the authoritative parameter entry. Alpha Cen A currently has
  Teff 5792 +/- 16 K, logg 4.30 +/- 0.01, reference [Fe/H] +0.20 +/- 0.04,
  xi 1.1 km/s, vsini 1.9 km/s, and vmac 2.3 km/s. Teff/logg/xi are pinned and
  [Fe/H] is solved. RYA-435 remains open over the cited Jofre value mismatch; no value
  was changed here.
- RYA-323 remains open: the corrected-metallicity target VALD re-extraction is not yet
  complete. This would be a later line-list completeness warning, not authority to use
  a solar-default list.
- RYA-1063 remains open, so the requested post-measurement literature comparison is not
  ready. Because measurement stopped first, literature was not used to tune or reject
  anything.

## Hard data/provenance gate

| required field | resolved value |
|---|---|
| holding_id | `alpha_cen_a_harps` |
| current local path | none on the workstation; data moved to Sirius |
| Sirius path | `/mnt/codex-data/spectra/_mac_import_20260816/Alpha Centauri (vetted)/Alpha Cen A/HARPS` |
| files | 205 FITS products |
| target identity | manifest rows resolve `HD128620` to `alpha_cen_a`; the historical folder label is not authoritative |
| source archive/product ID | individual ESO `ADP.*.fits` products; no single run product is selected |
| wavelength coverage | not proven for a selected product through a supported reader |
| resolving power/product type | registry notes `PRO CATG=SCIENCE.SPECTRUM`; exact resolving power is not proven for a selected product |
| continuum/normalization provenance | **UNKNOWN** |
| rest-frame status | manifest records `BARYCENT` for inspected HARPS rows, but no supported reader verifies the selected run product |
| checksum | **UNRESOLVED** because no exact science product/coadd is selected |
| source ticket | RYA-479 / RYA-806; measurement reachability finding RYA-1030 |

Committed `data/catalog/holdings_manifest_registry.csv` states:

- `evidence_state=verified` for the inventory/identity evidence;
- `telluric_applied=unknown` because the headers contain no reduction-chain,
  transmission-extension, or telluric-column evidence;
- `normalization_state=unknown`;
- there is no `HoldingSpec`/reader for `alpha_cen_a_harps` in
  `scripts/measure_band_ew.py::_INSTRUMENT_HOLDINGS`.

The presence of 205 files proves a holding exists. It does **not** prove which one is the
approved current product, its continuum contract, or that the abundance harness can read
it. RYA-1064 explicitly requires a STOP if the product has not reached the current
continuum/conditioning standard; therefore no checksum was manufactured by arbitrarily
choosing the first FITS file and no 205-file coadd was invented.

## Line-pool and engine disposition

The data gate precedes the line-pool and engine gates. Consequently:

| deliverable | status | reason |
|---|---|---|
| graded Fe I pre-run manifest | NOT EMITTED | normal production path cannot resolve a readable Alpha Cen A HARPS holding |
| Fe II diagnostic | NOT ATTEMPTED | same upstream gate; Fe II cannot redefine the requested Fe I headline |
| all valid engines | NOT ATTEMPTED | no approved common observed spectrum/line measurement exists |
| matched-line comparisons | NOT COMPUTABLE | no engine was served any line |
| A(Fe), uncertainty, trends, residuals, influence checks | NOT COMPUTABLE | scientifically correct STOP before measurement |
| literature comparison | DEFERRED | must occur only after measurement; RYA-1063 is also open |

No engine is labelled “unavailable” here: engine availability was deliberately not
conflated with an upstream data-product refusal.

## Orchestration sequence and classification

| step | classification | outcome |
|---|---|---|
| identify project/target/species scope | required Linear/history knowledge | ticket and parent define it |
| resolve `STAR_PARAMS` | required repository knowledge | automatic after knowing `config/stars.yaml` is authoritative |
| resolve registered HARPS holding | automatic | `alpha_cen_a_harps` found |
| locate bytes | required repository knowledge | path-register notes reveal the Mac-to-Sirius migration |
| choose exact approved product | failed/ambiguous | 205 source products, no declared run product/coadd |
| prove continuum state | failed/ambiguous | registry says `unknown` |
| resolve harness reader | automatic refusal | no `HoldingSpec` for the holding |
| enumerate graded Fe I lines | not reached | must follow approved spectrum selection |
| enumerate/run all valid engines | not reached | must share the same approved spectrum and lines |

Commands used for the reproducible preflight were read-only except for creating this
ticket worktree and report:

```text
git fetch origin
git worktree add -b ryandamienschmitt/rya-1064-experiment-mr-codex-orchestration-test-alpha-cen-a-harps-vis /Users/ryanschmitt/codex/rya1064 origin/main
rg -n -i "alpha cen|alpha_cen|hd 128620|71683|graded.gf|line.tier|engines=all|all.valid" README.md SEQUENCE.md CODEX_STATE_REGISTER.md config docs pipeline scripts
rg -n -C 3 "alpha_cen_a_harps" data config pipeline docs scripts
python3 scripts/derive_band_products.py --help
ssh sirius 'find "/mnt/codex-data/spectra/_mac_import_20260816/Alpha Centauri (vetted)/Alpha Cen A/HARPS" -maxdepth 1 -type f -name "*.fits"'
```

The nominal descriptor attempted was:

```python
RunDescriptor(element="Fe", ion="I", instrument="harps",
              holding="alpha_cen_a_harps", lo_A=3782.6, hi_A=6910.0,
              engine_deck="ts-lte")
```

On the workstation, importing the measurement holding table while resolving that
descriptor exits first on a missing **solar Kitt Peak** directory. Thus even a pure
Alpha Cen preflight has an unrelated solar-data import side effect before it can report
the intended `alpha_cen_a_harps is not wired` refusal.

## What worked automatically

- The canonical system ID, HARPS holding ID, registry evidence state, target identity,
  and Sirius storage root are discoverable from committed registries.
- The registry distinguishes an existing-but-unmeasured holding from missing data.
- The telluric policy preserves `unknown` rather than silently treating it as corrected
  or not required.
- The harness holding table makes clear that only the two solar HARPS holdings have
  readers; it does not silently fall back to one of them.
- The graded Fe policy and uncertainty machinery are present in committed code, but were
  correctly left downstream of the failed data gate.

## Friction / defects

### Critical

1. **No selected, approved Alpha Cen A HARPS science product.** The registry points to a
   205-row/source-product inventory, not one continuum-conditioned abundance input.
   Choosing or coadding one by agent judgment could silently change the science.
2. **Continuum state is unknown.** RYA-1064 specifically disallows legacy pre-fix
   normalization as validation evidence.
3. **The registered holding is unreachable by the measurement harness.** No
   `HoldingSpec`, reader, span, continuum contract, or reference-continuum declaration
   exists for `alpha_cen_a_harps`.

### Warning

1. `telluric_applied=unknown`. HARPS presently uses line-selection policy rather than a
   universal correction, but the product state still must be measured and recorded.
2. The current run descriptor has no target/star field and its `engine_deck` names one
   deck, so the requested `target=alpha_cen_a` and `engines=all_valid` cannot be expressed
   by the highest-level run object.
3. Resolving the descriptor imports a harness module with eager solar-atlas discovery;
   on the workstation it stops on absent Kitt Peak data before evaluating the named
   Alpha Cen holding.
4. RYA-323 and RYA-435 remain open, so corrected-metallicity VALD completeness and the
   [Fe/H] citation/value mismatch must remain visible in any later run.

### Suggestion

1. `python3 scripts/derive_band_products.py --help` crashes with
   `ValueError: unsupported format character 'b'` because an argparse help string
   contains an unescaped percent sign. The supported interface cannot currently explain
   itself to an operator.
2. A standard repository-only workflow is insufficient: the intended run order and
   several gates required historical Linear tickets to interpret.

## Follow-up tickets recommended

1. **Condition/register Alpha Cen A HARPS production product.** Select and document the
   exact exposure/coadd policy; prove source ADP IDs, identity, coverage, resolution,
   rest frame, continuum method/version, telluric state, checksums, and supersession of
   legacy products.
2. **Wire the approved product into the measurement harness.** Add a star-aware
   `HoldingSpec` and reader with declared span, `pre_normalised`, reference-continuum,
   and refusal tests. Do this only after follow-up 1 establishes the product contract.
3. **Add a declarative all-engine executor.** Accept the requested run spec, query valid
   engines, run a common line manifest, retain exact served subsets, and emit explicit
   unavailable statuses. The scientific modules remain independently runnable.
4. **Make preflight side-effect free and repair CLI help.** Registry resolution must not
   require unrelated solar bytes, and every supported CLI must render `--help` in a
   data-free environment.

## Closeout

Science result: **STOP, no abundance quoted.**  
Code changes: **zero.**  
Artifact added: this report only.  
Branch: `ryandamienschmitt/rya-1064-experiment-mr-codex-orchestration-test-alpha-cen-a-harps-vis`.

