# Model-atmosphere pointer

Large model grids are not committed here. The live abundance code uses the
iSpec-formatted `ATLAS9.Castelli` and `MARCS.GES` packs below `ISPEC_DIR`, not a
partial grid copied into this directory.

See [the model-assets guide](../../docs/models/assets.md) for classifications,
authoritative sources, expected paths, checksum policy, licensing caveats, and
validation commands. Run:

```bash
export ISPEC_DIR=/absolute/path/to/ispec
python scripts/validate_installation.py --full
python -m pytest tests/test_abundances_derive.py -q
```

Do not populate this directory with an ad hoc model subset and assume the live
iSpec interpolation path will discover it.
