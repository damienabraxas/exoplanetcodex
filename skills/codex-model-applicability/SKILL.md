---
name: codex-model-applicability
description: Audit which EW, synthesis, LTE, NLTE, 3D, or mean-3D products are genuinely applicable and in-domain for a scientific task. Use before an all-model sweep, when a model product is missing/NaN, or when reporting unsupported and out-of-domain paths. Does not run abundance science or select the preferred answer.
---

# Model-domain and applicability audit

## Purpose and boundary

Produce the complete matrix of model/route products that **may** run for the task and
the evidence for each disposition. Do not execute fits, extrapolate, convert one model
class into another, or choose the result closest to literature. Start with the authority
map in `LEDGERS.md`; scientific rules come from `SCIENCE_STANDARDS.md`, current model and
asset state from the ledgers/config/registries, and executable support from committed
code/tests.

## Build the candidate matrix

Enumerate candidates as independent stored axes, never as parsed display labels:

- route (`ew` or `synth`);
- scale (`1D-LTE`, `1D-NLTE`, `3D-NLTE`, or an explicitly registered class such as
  `MEAN3D_NLTE`);
- model/atom/grid family and vintage;
- atmosphere family;
- oscillator-strength source/grade;
- instrument/holding, band, species/ion, and line.

Derive display names using the current naming implementation. Preserve legacy labels
only as provenance. Do not infer route from a treatment string.

## Prove applicability per candidate

For every row record:

1. The task and protocol require or permit the route for this wavelength/line class.
2. A verified holding covers the line with usable pixels and passed data audit.
3. The line/species/ion exists in the required line list and model atom.
4. Stellar parameters, wavelength/transition, abundance, atmosphere, and other model
   dimensions are inside the verified domain. Inspect the underlying model domain; a
   derived extract's endpoints do not prove the physical model boundary.
5. Required assets, dependency versions, conversion products, handlers, and call paths
   are present and reachable. A green skipped test or swallowed import error is not
   support.
6. Applicable curation, telluric, frame, gf-grade, and ratified constraints can be
   applied without bypass.

Classify each row `APPLICABLE`, `UNSUPPORTED`, `OUT-OF-DOMAIN`,
`NOT-APPLICABLE`, `DATA-GAP`, or `AMBIGUOUS-STOP`, with evidence and the action needed.
A missing or invalid dependency is an implementation failure, not a physics-domain NaN.
Expected scientific absence may yield a dispositioned NaN; unexpected NaN fails loudly.

## Dimensionality honesty

Never label a mean-3D (`<3D>`) correction or atmosphere as full 3D radiative transfer.
Never call a correction grid a synthesis engine. Do not extrapolate beyond a verified
hull or silently clamp; only a specifically ratified bounded policy may do so, with its
flag and uncertainty. Keep LTE, NLTE, 3D correction, and 3D synthesis products separate.

## Deliverable

Return one row per independent product/line combination with:

`holding | band | line | route | scale | model | atmos | gf | domain evidence | code/asset evidence | disposition | reason`

Then summarize attempted/applicable/unsupported/out-of-domain counts. Absence from the
matrix is a defect. Any unresolved scientifically material ambiguity is a STOP before
execution and must be reported in the active Linear issue.
