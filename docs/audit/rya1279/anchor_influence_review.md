# RYA-1281: influential continuum anchors and source availability

Predeclared2026-11-05 before evaluation. Mac only, no RT. Start0c9a93b9,
mainf8c723a8, stacked base4635c742/draftPR598. Frozen canonical/lower-g spectra,
original corrected Kitt holding and both reference legs remain unchanged.

Manifest config/experiments/rya1281_anchor_influence.json declares0.1AIR-A bins
relative to each O777 line centre. Group every actual coefficient-support pixel
for both masks, with no reselection or clip. Sum signed centre-prediction weights
times (observed/IAG minus observed/synthesis), independently verify with linear
fits and the previous fixed-mask report. Keep every bin, including cancellation.
Rank only for plot readability (12largest absolute bins); full tables untruncated.
No noise/physical-opacity/independent-error budget inferred.

Read original canonical/GES source rows with exact species/wavelength/excitation
identity, retaining zero/multiple matches and gf differences. Inventory candidate
centres within declared bins, distinguish modeled41FeI, excluded combined FeI,
active-OI candidates and other candidates. This is availability, not a blend
identification, predicted line strength, completeness, molecular treatment, wing
coverage or permission to add opacity. No transition chosen from proximity alone.
Physical attribution requires separate source-qualified mixed-transfer controls.

## Findings and interpretation

Plan/module/script committed before evaluation at09e077c3 (full SHA in Git).
The accounting reconstructs all12case/line/mask centre differences within
1.485e-15, independently checked against fresh linear fits and prior fixed-mask
levels. All old source/pixel identities hash-match; no native fit or spectrum is
changed. These are the prior no-more-clipping counterfactual supports, not a new
adopted continuum. Signed bin contributions account for the centre level only;
finite-window area is nonlinear in the continuum and is cited separately from
the prior report, not claimed to be a sum of linear area errors.

Canonical dimensionless centre-level changes (IAG-denominator minus synthesis):

| AIR centre A | Frozen mask | Net level change | Sum of absolute pixel contributions |
|---|---|---:|---:|
|7771.944|synthesis|0.009035861|0.098676428|
|7771.944|IAG|0.006446399|0.018115321|
|7774.166|synthesis|0.007542229|0.109373607|
|7774.166|IAG|0.006607031|0.006607031|
|7775.388|synthesis|0.011871433|0.011871433|
|7775.388|IAG|0.006803349|0.006803349|

Synthesis7771/7774 absolute-to-net ratios10.92/14.50 expose strong cancellation
from one-sided extrapolated linear predictions. These sums are not uncertainties,
opacity fractions, independent errors or fractions of a physical line deficit.
Other supports bracket centre or have much less amplification; IAG7771 still
shows some cancellation. Lower-g and canonical contexts share observations;
repetitions are not independent observational measurements.

For canonical synthesis7771, the largest signed0.1A bins are
[7767.744,7767.844):+0.0240609,
[7766.944,7767.044):-0.0189535,
[7767.644,7767.744):+0.0178414.
For synthesis7774 they are
[7769.966,7770.066):+0.0248543,
[7769.166,7769.266):-0.0225471,
[7769.866,7769.966):+0.0160861.
None of these bins contains a canonical candidate centre. This does not establish
absence of line wings, molecular/uncataloged opacity, or a normalization offset.
A centre-only join cannot measure opacity coverage. No missing blend is identified.
At7775, the leading synthesis bin [7777.088,7777.188) contains an FeII candidate,
but presence alone supplies no predicted strength or permission to add it.

## Source availability, not opacity admission

All388canonical rows in7763-7783A retained. Treatments:41modeled FeI entries,
6explicit combined-FeI exclusions,3active-OI candidates,338other candidates.
Exact original GES identity uses species, printed AIR wavelength and excitation;
gf is compared separately so an adjudicated canonical value cannot be silently
replaced.354rows have one exact identity match,18have two,1has three,15have none.
Six uniquely matched rows have differing canonical/source gf values; both retained,
no fitted replacement. The inventory alone does not settle each adjudication.
Multiple/absent matches remain explicit; no arbitrary source component selected.

Original GES rad/Stark/waals fields are availability evidence only; their treatment
and units are not newly qualified here. Atomic centre matches do not provide
ionization/partition populations, wing profiles, HFS/isotopes, molecules or a
consistent mixed atmosphere/EOS. The original41FeI packet remains unchanged.
No new line opacity is installed. Source-qualified differential strengths are
still required before attributing any residual to these candidates.

## Reproduction, verification and handoff

```sh
source config/hosts/mac.env
/Users/ryanschmitt/codex/venv/bin/python scripts/audit_bride_anchor_influence.py --output /Users/ryanschmitt/codex/bride_mac_checks/rya1281-anchor-influence-v1
```

Exit0. Immutable registered experiments.bride_anchor_influence output, with three
results/two visually inspected plots copied byte-for-byte to
 data/results/rya1281/anchor_influence_v1 and data/output/rya1281/anchor_influence_v1.
Largest file520433bytes. Derived ratios were already archived; no external raw
spectral flux arrays copied. Source fragments preserve original field values.
Full blocks.csv retains every bin and candidate ID; plots show at most12largest
absolute bins per context, never a new support selection.

[Canonical influence](../../../data/output/rya1281/anchor_influence_v1/canonical_refined_influence.png),
[lower-g influence](../../../data/output/rya1281/anchor_influence_v1/lower_gravity_refined_influence.png).

New focused tests3passed in0.11s; expanded science-env suite71passed in7.23s,
292existing backend printUpdate warnings. Exact expanded command is the prior
19-module suite in reference_resolution_review.md plus
 tests/test_bride_anchor_influence_rya1281.py. Tests independently reconstruct
linear-fit differences, prove constant/slope controls at centre, preserve signed
cancellation/bin closure and refuse invalid/unordered supports or bin widths.

Self-review: masks are actual coefficient supports rather than postclip survivors;
ratios share the exact observed grid; source/pixel/report/packet hashes validated.
Each canonical candidate is accounted for, source ambiguity retained, gf not
silently harmonized. Bins do not create a physical opacity/uncertainty inference.
No new bug, silent fallback or source alteration discovered in this phase; earlier
clipping, fractional-flag and numerical-parity ledger entries retained.

Next concrete work is a source-qualified mixed-opacity/continuum differential
covering the audited regions and line wings, with independently justified line
strengths/partition/EOS/damping/components and molecular coverage. Preserve the
fixed O/background contract and source-only negative control; do not add a
candidate merely because it lies nearby. Independent O atom/EOS/reference and
native FTS calibration, broader xi/T/g/composition interactions, additional O
lines, C/N then Fe then reference star and adaptive3D/all-band goals remain OPEN.
No RT, Sirius, fitting/adoption, default continuum, production, merge or publication
change. RYA-1281 In Progress / RYA-1278 In Review; scientific HOLD.

Final checks: MIRROR targets/v214 pointer/whitespace pass. Generator guard retains exactly2436existing diagnostics, zero added/removed versus prior resolution audit. Fresh library audit194documents/127accounted/67unindexed/12off-Mac registered paths. Existing historical RYA-587 ledger failure and gold/tracker debts remain documented; no green CI claim. Three results/two figures tracked, largest520433bytes. No new bug discovered; no empty ledger placeholder added.
