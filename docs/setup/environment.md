# Environment and prerequisites

## Compatibility matrix

| Item | Supported/tested | Notes |
|---|---|---|
| macOS | Tested: macOS 15, Apple silicon | Xcode command-line tools needed to compile external engines |
| Linux | Expected: current x86_64 distributions | Not exercised for this documentation revision |
| Windows | Unsupported natively | Use a Linux VM/WSL experimentally; iSpec engines are Unix-oriented |
| Python | Supported 3.9–3.11; tested 3.9.6 | 3.12+ is not validated with the pinned scientific stack/iSpec |
| Environment | `venv` documented | Conda/mamba may work but no lockfile is committed |
| Compiler | C/C++ and Fortran for full iSpec engines | Core correction example uses binary wheels only |

The core checkout is about 0.7 GB because substantial line-list material is
versioned. Budget 2 GB for clone, environment, and test outputs. A tested full
iSpec checkout with its input bundle occupies about 9.3 GB (ATLAS9.Castelli
about 52 MB, MARCS.GES about 39 MB, and molecular lists about 4.3 GB on the
reference workstation). Keep at least 15 GB free; archive spectra and synthesis
caches may require considerably more.

The core smoke example needs one CPU and under 2 GB RAM. Full synthesis is
line/window dependent; start with 4 cores and 8 GB RAM, avoid unvalidated
multiprocessing, and benchmark before scheduling an HPC array.

## Python environment

`requirements.txt` declares supported minimums.
`requirements-lock-py39.txt` records the exact Python 3.9 environment tested
for RYA-642. Python 3.10/3.11 users must resolve the minimum declaration and
record `python -m pip freeze` with every science run. The tested setup is:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock-py39.txt
python scripts/validate_installation.py
```

Core packages are NumPy, SciPy, pandas, Matplotlib, Astropy, Specutils,
Astroquery, PyYAML, tqdm, and pytest. `lockfile` is required by the tested iSpec
tree; Cython, dill, h5py, hdf5plugin, and statsmodels support its imported
runtime. There is no supported editable/package install: run from the repository
root so `config` and `pipeline` resolve.

## Full iSpec and radiative-transfer environment

This repository imports iSpec from a sibling checkout by default or from
`ISPEC_DIR`. The reference installation is upstream
[`marblestation/iSpec`](https://github.com/marblestation/iSpec), commit
`21269948750e2177fe32c5df38e1d9bb81363374`. That commit is **known to work**,
not a promise that newer commits are incompatible.

```bash
git clone https://github.com/marblestation/iSpec.git /absolute/path/to/ispec
cd /absolute/path/to/ispec
# Follow iSpec's installation instructions and unpack its matching input bundle.
cd /absolute/path/to/exoplanetcodex
export ISPEC_DIR=/absolute/path/to/ispec
python scripts/validate_installation.py --full
```

iSpec supplies the integration layer, compiled Turbospectrum and MOOG binaries,
solar abundance tables, atmosphere packs, line regions, opacity data, and
molecular lists used by the live synthesis modules. Do not install the unrelated
PyPI package named `ispec` and assume it is equivalent. Follow the upstream
[iSpec documentation](https://www.blancocuaresma.com/s/iSpec) and
[Turbospectrum source](https://github.com/bertrandplez/Turbospectrum2019).

Full validation expects:

```text
$ISPEC_DIR/input/atmospheres/ATLAS9.Castelli/
$ISPEC_DIR/input/atmospheres/MARCS.GES/
$ISPEC_DIR/input/linelists/turbospectrum/molecules/
$ISPEC_DIR/synthesizer/turbospectrum/exec-gf/bsyn_lu
$ISPEC_DIR/synthesizer/moog/MOOGSILENT
```

These assets are mandatory only for their corresponding full workflows. They
are not required by the repository-native quick start.

## Platform notes

- On macOS, install Xcode command-line tools before compiling engines. The
  system Python 3.9 may emit a LibreSSL warning through urllib3; use a
  python.org/Homebrew Python linked to OpenSSL for archive downloads.
- On Linux, install the distribution equivalents of a C compiler, a Fortran
  compiler, BLAS/LAPACK development libraries, and FITS libraries if wheels are
  unavailable.
- Run synthesis serially first. Forked workers can inherit non-thread-safe
  Fortran state and exhaust temporary storage.
- HPC deployments should place read-only model grids on shared storage and
  give each job a private temporary/output directory. Record scheduler,
  compiler, CPU, and module versions.
