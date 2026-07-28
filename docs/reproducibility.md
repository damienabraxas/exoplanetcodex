# Reproducibility and provenance

Every result bundle should identify enough state for another researcher to
reconstruct the calculation:

```text
repository_commit
git_status (must explain local changes)
python_version, platform, compiler
pip_freeze or conda_explicit_spec
iSpec_commit and input-bundle identity
Turbospectrum/MOOG executable checksum
atmosphere-pack file checksums
atomic/molecular-data versions and checksums
spectrum archive IDs, reduction version, frame, and SHA-256
config/stars.yaml and data/method_policy.yaml SHA-256
command line and working directory
random seeds (or “none used”)
start/end UTC, CPU/thread count
output file SHA-256 and scientific validation verdicts
```

Suggested capture:

```bash
mkdir -p results/tables/provenance
git rev-parse HEAD > results/tables/provenance/repository-commit.txt
git status --short > results/tables/provenance/git-status.txt
python -m pip freeze > results/tables/provenance/pip-freeze.txt
git -C "$ISPEC_DIR" rev-parse HEAD > results/tables/provenance/ispec-commit.txt
shasum -a 256 config/stars.yaml data/method_policy.yaml \
  > results/tables/provenance/config.sha256
python scripts/validate_installation.py --full \
  --json results/tables/provenance/install-validation.json
```

For stochastic inference, explicitly pass and record a seed. For deterministic
EW/synthesis paths, record that no random sampling was used. Do not overwrite a
previous run directory; include target, date, commit, and method in its name.

## Verification

The quick-start manifest contains input values, output, grid hashes, commit, and
citation. A successful reference evaluation currently gives O I 777.194 nm
`delta ≈ -0.169 dex` and `A(O) ≈ 8.521`. This tolerance check validates grid
parsing/interpolation only.

Full science verification must additionally compare output row counts, line
identities, parameter values, abundance-scale labels, per-line scatter, Fe
ionization balance, reduced-EW slopes, synthesis residuals/chi-square, correction
coverage, and output hashes or declared floating-point tolerances. A byte hash
alone cannot establish scientific validity.

## Citation template

Replace brackets with the exact run assets:

> Exoplanet Codex science pipeline, commit `[SHA]`, MIT license,
> https://github.com/damienabraxas/exoplanetcodex. Methodology and target
> parameters: `[repository documentation and primary stellar references]`.
> Spectra: `[archive, program/observation IDs, instrument, reduction]`.
> Atmospheres: `[Castelli & Kurucz 2004 / Gustafsson et al. 2008, pack version]`.
> Radiative transfer: `[Turbospectrum / MOOG version and primary citations]`.
> Atomic/molecular data: `[VALD3, NIST ASD, ExoMol/other releases]`.
> NLTE/3D corrections: `[only the grids actually applied]`.

Also cite Asplund et al. (2021) for the repository's solar abundance reference.
Each data provider may impose additional acknowledgement or usage language;
consult the upstream license at retrieval time.

Use the [categorized project bibliography](references.md) to distinguish
methodology, software, model/grid, instrument, archive, atomic/molecular, and
reference-dataset citations. Cite only resources actually used by the run.
