# RYA-517 — reference compute stack: compatibility PROBE (no build/migrate)

Probe only, per the scope: determine the newest Python + numpy each engine supports
(compiles + tests green), report the per-engine ceiling + the recommended stack. NO build,
NO migration. Branch `ryandamienschmitt/rya-517-compat-probe`.

## Headline — iSpec is NOT the low cap the ticket feared; it's the *floor-raiser*
The current reference stack (Mac **Python 3.9.6 + numpy 1.26.4 + scipy 1.13.1 + astropy 6.0.1**)
is **older than iSpec's own declared requirements** (`numpy>=2.2.5, scipy>=1.15.3, astropy>=7.1.0,
pandas>=2.2.3, Cython>=3.0.12`). iSpec (v2023.08.04-20-g2126994) is a modern-stack codebase; the
accidental EOL environment is behind it, not ahead of it.

## Per-engine ceiling (evidence-based)

| engine | newest supported | numpy | evidence |
|---|---|---|---|
| **iSpec (C synthesizer)** | **Python 3.12 PROVEN**; 3.13 likely; 3.14 unverified | **2.x** (req ≥2.2.5) | ships compiled `synthesizer.cpython-39-darwin.so` AND `synthesizer.cpython-312-darwin.so` → builds on 3.9 (vs numpy 1.x) AND 3.12 (vs numpy 2.x); `requirements.txt` numpy≥2.2.5; **zero numpy-2.0-removed APIs** in its Python code (grep clean: no np.float/np.NaN/np.infty/np.round_/np.in1d/…) |
| **Turbospectrum** | Python-AGNOSTIC | n/a | external Fortran binary under `synthesizer/turbospectrum/`, invoked by path via subprocess; no Python/numpy cap. Build needs gfortran. |
| **MOOG** | Python-AGNOSTIC | n/a | external Fortran binary `synthesizer/moog/MOOGSILENT` (`Abfind.f` …), invoked by path (`ispec/common.py:59`); no Python/numpy cap. Build needs gfortran. |
| **scipy** | ≥ Python 3.14 | 2.x | scipy 1.18.0 installs cleanly on Sirius py3.14 |
| **astropy** | ≥ Python 3.14 | 2.x | astropy 8.0.0 installs cleanly on Sirius py3.14 |
| **numpy** | ≥ Python 3.14 | — | numpy 2.5.0 installs cleanly on Sirius py3.14 |
| **pandas** (also required) | ≥ Python 3.14 | 2.x | pandas 3.0.3 installs cleanly on Sirius py3.14 |

Data points:
- Sirius (Python 3.14.4), fresh venv: **numpy 2.5.0 / scipy 1.18.0 / astropy 8.0.0 / pandas 3.0.3**
  all install (wheels exist for 3.14). So the pure-Python-lib ceiling is ≥ 3.14.
- The Mac's pip only *offers* up to numpy 2.0.2 / scipy 1.13.1 / astropy 6.0.1 — that is the py3.9
  EOL cap showing itself, not the libraries' real ceiling.

**The single binding constraint is iSpec's C-extension COMPILE**, and it is proven at **Python 3.12
+ numpy 2.x** (the highest with an existing compiled `.so`). 3.13 is a low-risk step (all libs
support it; iSpec's `Cython>=3.0.12` covers 3.13; no removed APIs). 3.14 is higher-risk: all the
Python libs install, but iSpec's pinned Cython (3.0.12) predates 3.14 support (needs Cython ≥3.1),
and no iSpec 3.14 build exists — it must be compiled to confirm (deferred; this is a probe).

## Recommended reference stack
**Python 3.12 + numpy 2.2.x (≥2.2.5) + scipy ≥1.15.3 + astropy ≥7.1.0 + pandas ≥2.2.3.**

Rationale: the **newest stack with a PROVEN iSpec C build** (the binding constraint) where every
other dependency is mature and fully wheel-supported. It clears the EOL problem (3.9→3.12, EOL Oct
2028 vs 3.9 EOL Oct 2025), moves onto numpy 2.x as iSpec itself requires, and carries zero unknowns
— i.e. genuinely "compiles + tests green today" without a speculative build.

Upside, pending one confirmatory iSpec recompile (low risk, deferred): **Python 3.13** — same
numpy 2.x, all libs support it, iSpec code is clean and its Cython pin covers 3.13. Worth
confirming when the migration runs.

Not recommended yet: **Python 3.14** (Sirius's system default) — the Python libraries all install,
but the iSpec C build on 3.14 is unproven AND its Cython pin needs bumping to ≥3.1; adopting it
means a real build + a Cython bump, so it is not a "tests-green-today" choice.

## Guard note (per the ticket)
iSpec does NOT cap the stack below current — it *raises* the floor (numpy 1.x → 2.x, py3.9 → ≥3.12).
No OPEN_QUESTIONS "modernize the synthesizer" item is warranted; the synthesizer is already modern.
The only thing NOT provable without a build is 3.13/3.14 iSpec compilation — flagged for the
migration step (RYA-517 step 2 / RYA-511 Phase 0), not this probe.

RYA-514's force-fork + single-thread BLAS carry to any of these stacks unchanged (platform-agnostic).
