# RYA-1281: FeII feature comparison and native broadening source audit

Predeclared2026-11-06 before new measurements. Start2235bab9, main9a74d843
unchanged, stacked base4635c742/draftPR598. Mac only; no new RT in this phase.
Historical zero/FeII spectra and both source columns remain frozen, no fitted
gf/width/RV/LSF, source/continuum adoption or abundance inference.

Manifest config/experiments/rya1281_feii_profile.json fixes core+/-0.1A and
symmetric shoulder bands0.35-0.55A from each unchanged fixed-RV observer centre.
Retain raw shipped Kitt and independent IAG native/Gaussian-reference views.
Fixed shoulder contrast cancels a multiplicative level; it is shape bookkeeping,
not an isolated EW, new ORCH continuum, profile fit or physical noise budget.
Interpolation0.005A adds no resolution or independent samples. Preserve both
no-gravity/plus-gravity conventions; no post-hoc recentring or shoulder movement.
Companion lines/molecules can affect both core and shoulders; don't assign the
whole observed feature to FeII. Actual observation arrays remain read-only.

Inspect pinned source21269948750e2177fe32c5df38e1d9bb81363374: original
anstee/bsyn/depth/IRWIN/GES inputs and iSpec export. Execute actual FeII pressure
priority via isolated original anstee routine, then exact native width fragments.
Test stage2ionization and provided versus exported-zero/Cowley Stark separately
across3T(4000/5772/10000),3H(1e19/1e21/1e23m^-3),3ne(0/1e18/1e20m^-3).
Positive H avoids native PH denominators; H2 explicitly0. Compare all54total
widths under the existing128*float32eps arithmetic guard, not a physical-error
criterion. Then evaluate every native depth with saved actual H0/He0/ne values,
without changing opacity or spectra. Never silently use a neutral-Fe Cowley
coefficient for FeII or promote a supplied coefficient over a source lookup.

Existing bibliography barklem2005_fe2 found before searching/download. Held
library scan found no matching local PDF. A&A web page403; author arXiv abstract
https://arxiv.org/abs/astro-ph/0502098 accessible. Abstract describes24188FeII
H-collision calculations and a Kurucz-list gf>-5 applicability selection; target
RU-tag gf=-5.33 does not itself prove or disprove original table coverage. Actual
native lookup is required. Abstract only inspected, no full-paper validation or
new document acquisition claimed. DOI10.1051/0004-6361:20042469.
