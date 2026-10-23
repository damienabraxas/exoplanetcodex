# RYA-1250 investigation: scope, uncertainty and culls

Current scope: Ryan confirmed both Mg I and Mg II, following Asplund as closely as possible. The earlier Mg-I-only recommendation below is superseded. See [Mg II readiness](rya1250_mgii_nlte_readiness.md).

Historical investigation snapshot, superseded by Ryan's subsequent authorization to apply the models and continue. The current implementation uses shared species/window preparation and declared solar corrections; six saturated lines are retained for synthesis.

2026-10-18. Ryan initially authorized investigation, not one-off orchestrator fixes.
No orchestrator code, canonical gf, culls or products were changed in this investigation.

## Earlier scope interpretation (superseded)

Mg II came from the original ticket, which explicitly requested Mg I + Mg II and
AGSS21's combined-ion headline. Ryan's subsequent instruction supersedes that:
continue with neutral Mg I; ionized Mg is not a prerequisite. Retain the Mg II
transcription as literature evidence, not required measurement work.

`run_matrix.graded_ions("Mg")` currently returns `["I"]`. Explicit
`expand(..., ions=["I"])` also produces only I. The existing run CLI supports
`--element Mg --ion I`. No Mg-specific switch is needed. However `--prepare`
currently does not pass `args.ion`, so it considers all declared sets for the
symbol. Scope propagation into prepare is a general interface gap, not evidence
that the science requires a second ion. Publication selects units from the actual
run report, while coverage and withdrawal operate at element scope; preserve this
limitation for any future mixed-ion run rather than assuming every stage is scoped.

## The uncertainty recipe already exists

The actual call chain is `element_publish.publish` ->
`scripts/rya1230_cno_budget_legs.py` -> `scripts/rya1230_cno_budget.py` ->
`pipeline.uncertainty_contract`. The historical CNO filename is misleading;
the assembler takes species from each product and is used for other elements.
It handles feature scatter/SEM or single-line fit curvature, paired nuisance
legs, source-correlated gf errors, and NLTE/3D terms. The implemented ordinary
NLTE term is half the pool's correction with a 0.03 dex floor; no Mg-only recipe
needs to be invented. `gf_error_model.csv` already declares source correlation.

The earlier checkpoint's wording was too broad: this is an input/admission
problem, not absence of an uncertainty recipe. Three distinct gates exist:

- prepare adopts a published gf only when the source row also supplies sigma;
- published-set execution reads the graded subset, not the full transcription;
- budget assembly uses stored sigma or the cited NIST class, and explicitly holds
  a product when a line has neither.

Older grading code has an ungraded systematic allowance, but the current
published-set/publishing path does not use it as a general fallback. Its Kurucz
semi-empirical allowance cannot simply be relabelled as a measured uncertainty
for the AGSS21 theoretical Mg gf. Any change to that policy is cross-element
scientific policy, not a numerical fix to get Mg through.

## Cull finding

The six proposed culls are exactly all six currently graded Mg I lines:
4571.096, 5167.322, 5172.684, 5183.604, 5528.405, 5711.088 A. Applying all six
would leave zero graded Mg I lines; both AGSS21 subsets currently have zero too.
The 4730.029 NLTE table entry does not imply that line is currently gf-graded.

`element_prepare` measures reduced EW and depth, checks absence from REGISTERED
solar-literature sets, then proposes `exclude/active`. That action removes a
line from synthesis aggregates as well as EW. The broad pre-existing Mg registry
row instead says the b triplet and 5528 require synthesis. Its shorthand scope
(`b-triplet 5167/72/83 + 5528 (strong)`) does not match any of the exact lines in
`disposition_for_line`; checked all six. Prepare loads the disposition table but
does not consult it in the candidate loop. Thus the older science routing is
not protecting these lines. Saturation alone establishes difficulty for EW;
it is not proof that profile synthesis must discard a line.

## A data-only route is available

Visually inspected Magg et al. 2022, local author manuscript
`Magg2022_solar_standard_composition_2203.02255.pdf`, Table 1 page 4, Sect. 3.1.4
page 11, and per-line Table A.1 page 18. Both 5528 and 5711 are explicitly used
for solar Mg I synthesis. Their Table 1 adopted inputs are:

| Air wavelength (A) | Lower energy (eV) | log gf | Stated error (dex) |
|---|---|---|---|
| 5528.405 | 4.346 | -0.547 | 0.02 |
| 5711.088 | 4.346 | -1.742 | 0.05 |

Their per-line solar results are 7.51 and 7.59, not a reproduction claimed here.
The paper cites Pehlivan Rhodin 2017, but the inspected 2016 arXiv manuscript
of that atomic paper gives experimental 5711 gf -1.842 and theoretical -1.742.
Therefore cite the exact adopted values/errors to Magg Table 1; do not claim
that the different PR2017 experimental value was independently verified to
match. This provenance discrepancy remains visible.

A temporary in-memory declaration probe uses the existing prepare code with
these two source rows and writes only under `data/output/rya1250/`. It does
not register a new production set or apply gf changes. The probe returned four
cull candidates (4571 and the b triplet), two gf adoptions with the exact errors
above, and one unclassified source-correlation entry. This confirms the data-only
route mechanically; it is not an approved scientific adoption. Registering this
literature selection and its source-correlation declaration is ordinary input
data, not a special Mg branch in the orchestrator. It would preserve both
solar lines and give the existing budget published error inputs.

Recommended direction: Mg I only, reuse the existing budget, resolve published
line selections through data declarations, and keep the remaining culls pending
method-specific evidence. Do not silently replace the requested AGSS21 set with
Magg or manufacture sigmas for AGSS21. The choice of primary published selection
must be explicit. Any infrastructure work should address cross-element scope,
readiness for empty required pools, and route-specific cull decisions together.

Evidence: `data/output/rya1250/investigation.json`,
`data/output/rya1250/magg_declaration_probe/solar_Mg.json`, and source images in
`data/output/rya1250/`. No measurement, full feed test suite or publication was
run for this read-only investigation.
