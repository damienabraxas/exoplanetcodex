# RYA-1202 — near-UV continuous opacity: what Fe I bound-free actually carries

Diagnostic only. Nothing here is a product and nothing writes to `data/`.

The question was whether our synthesis uses proper Fe I photoionization opacity or a
hydrogenic approximation — the ~50% lever Bell, Balachandran & Bautista (2001) identify
for the missing near-UV opacity. **It uses the Opacity Project R-matrix data, not a
hydrogenic approximation**, so that premise is refuted; the measurements here quantify
what the component carries and what refining it could recover.

## The chain that was verified, not assumed

`pipeline/abundances_derive` calls `ispec.generate_spectrum(code='turbospectrum')`;
iSpec's `atmospheres.calculate_opacities` runs `babsma_lu` in a temp dir with its bundled
`synthesizer/turbospectrum/DATA` symlinked in, and passes **no** `CONTINOPAC` key, so
babsma's `input.f` default `DATA/jonabs_vac_v19.2.dat` is what our near-UV runs on. That
file is **byte-identical** (sha256 `867f5a3f…`) to the one in `engines/Turbospectrum_NLTE`.

## Method

`babsma_pair.py` runs babsma twice over 3000–3800 Å on a solar MARCS model
(5772 / 4.44 / 0.0, ξ=1.0 — the feed's own pin), changing **only** the Fe I bf
cross-sections, which `compare_opac.py` then reads out of the `MODELOPAC` file.

The control is the point: the `stock` arm uses a copy that is byte-identical to the
shipped table, so any `stock`-vs-production difference would be the harness moving
something, and a `stock`-vs-`fezero` difference is attributable to Fe I bf alone.

`synth_pair.py` repeats the swap one level up, monkeypatching iSpec's own
`calculate_opacities` so that model, abundances, line list and wave grid are the
production call unchanged. `bracket.py` and `feedback.py` interrogate the response.

## What it found

* Fe I bf carries **11% of the continuous opacity** at the continuum-forming layer across
  3000–3800 Å — 20.8% at 3000 Å falling to ~2–4% by 3700 Å, and up to 48% in the outer
  layers. It is not a small correction.
* **The window's abundance response turns over.** With Fe I bf on, mean depth stops
  responding to A(Fe) above ~+0.4 dex and then reverses. With Fe I bf off the response is
  monotonic. One-variable test, `feedback.py` — the turnover is caused by Fe I bf, because
  Fe I bound-free opacity **scales with the Fe abundance being measured**.
* So a dex equivalent obtained by extrapolating a linear fit is an artifact.
  `bracket.py` refuses to report one: the effect is not bracketed even at +1.5 dex.

## Reproduce

    python scripts/diagnostics/rya1202/babsma_pair.py   # needs solar_marcs.mod alongside
    python scripts/diagnostics/rya1202/compare_opac.py

Both need `PYTHONPATH` to include the iSpec checkout. They write only into their own
run directories.
