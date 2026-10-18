# RYA-1250: Mg II NLTE readiness

2026-10-18. Ryan authorized model application and continued orchestration for
both ions. The sections below retain the initial investigation as history.

## Current application status

**Correction to the initial inventory:** AGSS21 Table A.1 already publishes
all eight model abundance columns, including MARCS LTE/NLTE and 3D LTE/NLTE.
Those columns are now transcribed for all 14 Mg features. The declared solar
correction route uses MARCS NLTE minus MARCS LTE for the 1D leg, and 3D NLTE
minus MARCS LTE for the 3D leg. These are cited additive corrections on our
LTE fits, not a new in-transfer 3D calculation; our atmosphere can differ.
Full reproduction of the authors' transfer calculation remains unestablished.

Five Mg II lines have current NIST A/A+ grades verified against their level
energies, wavelength and rounded oscillator strength. The exact AGSS21 gf is
retained. Errors use the existing project conversion from NIST accuracy grades,
not invented AGSS21 per-line standard deviations. The 10092 feature remains
unresolved; all eight AGSS21 Mg I uncertainties remain unresolved.

The Mg I PySME grid extension completed on Sirius: 22 retained + 44 new = 66
rows. The validated RYA-534 Gerber Mg deck is registered in the shared adapter;
its aux table resolves 21 abundances from 6.55 to 8.55. Its atom supports Mg I
only, so Mg II Gerber-NLTE cells are explicitly held. The six saturated Mg I
lines are retained with synthesis-required dispositions; no exclusions applied.
Preparation now honors these dispositions and scopes unresolved line inputs to
the relevant ion and wavelength window. The existing shared uncertainty recipe
is reused by the orchestrator publication stage.

The published-set source declares `aggregation: mean`; the shared product builder
honors that recipe for LTE and both correction legs while other pools retain
their existing median. Coverage accounts for unavailable cited 3D corrections
as explicit holds after the underlying LTE pool has been handled. Canonical gf
adoption preserves unrelated high-precision fields; only the five intended rows
change in this branch. The full run was restarted on this implementation.

The same declared statistic now reaches budget restatistics, paired perturbations,
and NLTE/3D effects. Citation provenance is written into the actual product CSV,
and names a cited table only when that table supplied a correction to the pool.
Declared solar sets with no graded lines appear as held coverage cells.

The resumed production run uses `orchestrator_run_v3.log` on Sirius. The latest
capability/grid/recipe regression pass is 65/65. The full-suite comparison against
the same main commit found 35 shared failures; grid regressions were corrected and
rechecked, the old "Mg has no deck" fixture now uses an unregistered element, and
the seven historical clean-working-tree checks must be rerun after the reviewed
canonical adoption is committed. Production results remain pending.

## Initial investigation (superseded where noted above)

The following records describe the pre-application snapshot. References to
22 rows, zero graded Mg II rows, unapplied gf changes, or missing published 3D
correction data are historical, not the current status.

## Assets inspected on Sirius

| Asset | Verified content | Consequence |
| --- | --- | --- |
| `amarsi_galah/nlte_Mg_scatt_pysme.grd` | 229 Mg I, 86 Mg II, 1 Mg III exported fine-structure levels | Candidate Mg II 1D NLTE substrate exists |
| Same grid, nearest solar node `t5750_g+4.50_m+0.00_a-0.00` | Finite positive, nonunity departure coefficients for all six features' endpoint candidates | Actual departure data present; not yet an interpolated synthesis validation |
| `gerber_ts/atom.mg86d` and associated binary/aux files | 85 Mg I levels, one Mg II ground level; all 247 bound-bound transitions are Mg I | This staged TS atom cannot synthesize Mg II in NLTE |
| `data/nlte_grids/Mg_Amarsi2020_PySME.csv` | 22 rows, Mg I only (4730/5711) | No production Mg II corrections yet |

Grid root is `/mnt/codex-ext/codex-grids/nlte` (registered alias
`/mnt/codex-data/grids/nlte`). Audits copied to
`data/output/rya1250/mgii_amarsi_grid_audit.json` and
`data/output/rya1250/mg_atom_ion_reach.json`. Exported fine-structure counts
must not be equated to model-atom superlevel counts.

Mg II endpoint matches were found for 7877.0, 7896.4, 9218.2, 9244.3,
10092.2 and 10914.2 A within 0.00002 eV of the transcribed AGSS21 energies.
The 10092 feature has multiple fine-structure candidates; energy matching
alone does not establish the allowed components or their gf allocation.

## Initial readiness gaps (implementation update below)

1. **General ion-aware synthesis and derivation.** `pipeline/pysme_nlte.py`
   hardcodes neutral species and ionization in `_linelist_rows`; level-label
   matching also needs ion filtering because energy origins overlap.
   `scripts/extend_nlte_grid.py` selects neutral lines and writes ion 1.
   Prepare's extension discovery is neutral-only too. Fix these as shared
   ion-aware capabilities, not Mg-specific exceptions. Test ion identity,
   endpoint selection, actual departure engagement and explicit unsupported
   route failures before deriving or publishing corrections.
2. **Atomic inputs.** Resolve the absent canonical 10092.2-A feature,
   component allocations, radiative/electron/H broadening and line-specific
   gf provenance/errors. AGSS21 describes NIST A/A+ accuracy collectively;
   this is not a per-line standard deviation. Both published graded subsets
   currently have zero rows. No arbitrary errors or broadening defaults.
3. **Observational validation.** Verify spectral windows, blends and tellurics
   for all six red/NIR features; compare solar synthesized profiles/EWs and
   derived abundances to the source. The nearest-node binary inventory is
   not a solar abundance or correction validation.
4. **Full 3D NLTE inputs.** AGSS21 used a 48 Mg I + 27 Mg II + Mg III ground
   model atom with 482 bound-bound and 71 bound-free transitions. The staged
   1D GALAH grid is not evidence that this same calculation is reproduced.
   No matching Mg 3D input was found in the inspected grid inventory.
   Obtain/identify the authors' solar synthesis/correction data for these
   six lines or a reproducible equivalent, with collision and broadening
   prescriptions checked. Do not equate 1D NLTE with full 3D NLTE.
5. **Uncertainty and combination.** Reuse the existing orchestrator budget
   infrastructure; resolve its source inputs and compare its recipe with
   AGSS21. Preserve separate ion results and combine using their total
   statistical/systematic errors. Culls remain unapplied and need
   method-specific review; selecting AGSS21 lines does not justify deleting
   other scientifically usable lines from every route.

## Primary references

- [Asplund et al. 2021](https://doi.org/10.1051/0004-6361/202140445),
  Sect. 2.4 and Table A.1: eight Mg I and six Mg II lines, full 3D NLTE;
  Mg II mean 7.52 versus 7.58 in 3D LTE. Local source was visually checked.
- [Amarsi et al. 2020](https://arxiv.org/abs/2008.09582): GALAH departure
  coefficients on 1D MARCS atmospheres; the staged candidate grid family.
- [Matsuno et al. 2024](https://arxiv.org/html/2405.13486v1), Sect. 2.1:
  Balder calculations with the AGSS21 model atom. Its Table 1 lists the
  Mg II 448.1-nm multiplet, not the six AGSS21 red/NIR features; this paper's
  existence alone does not establish correction coverage for our selection.

## Continuation: implementation and real synthesis checks

The shared functions now accept an ion stage (neutral remains the default):
`auto_labels`, `labels_by_selection_rules`, `_linelist_rows`, `_synth_ew`, and
`nlte_delta`. Grid energy matches are restricted to the requested species.
`extend_nlte_grid.py --ion II` carries that stage through discovery, canonical
lookup, synthesis-list caching, workers, CSV rows and per-ion failure provenance.
Prepare coverage uses ion plus wavelength. Historical wavelength-only failure
provenance is explicitly neutral, matching its old writer.

Ionized rows require explicit radiative/Stark damping and matching component J
values. A mixed/unidentified component set raises instead of choosing its first
row. Native PySME flags must confirm departure application for every component.
The environment must reset departures before the LTE comparison. The returned
result now includes its LTE curve and target NLTE EW for diagnosis. This is shared
infrastructure; no Mg-specific branch or new uncertainty recipe was introduced.

Synthetic unit/regression checks: **97 passed, 4 environment-dependent skips**
across ion matching/coverage/writing, HFS, Al level reach and prepare, departure
reset, b-factor synthesis, hull bounds, orchestrator sweep, Si and coverage tests.
The dedicated Sirius PySME interpreter has no pytest installed; real synthesis
below ran there and passed the live reset and departure-engagement guards.

Solar technical smoke, using existing canonical gf and catalog broadening:

| Published feature (A) | Outcome |
| --- | --- |
| 7877.0 | NLTE engaged; diagnostic delta -0.02176 dex |
| 7896.4 | NLTE engaged; diagnostic delta -0.03036 dex |
| 9218.2 | NLTE engaged; no valid correction at either +/-0.2 or +/-0.8 dex brackets |
| 9244.3 | NLTE engaged; no valid correction at either bracket |
| 10092.2 | Held before synthesis: no canonical feature and unresolved components |
| 10914.2 | NLTE engaged; no valid correction at either bracket |

For 9218, the wider LTE curve peaks at 60.897 mA then declines, while the NLTE
EW is 65.033 mA. Thus the failed inversion is not solved by a wider bracket;
extrapolating or clipping would invent a result. Its physical/numerical cause
remains unresolved. Do not infer this exact cause for the other two failed lines
without their own curves. [Diagnostic plot](../../data/output/rya1250/mgii_9218_cog.png).
The neutral reference lines reproduce the earlier check: 4730.029 = -0.0088448
and 5711.088 = -0.0344944 dex. No production table was updated.

The candidate Mg II inputs have catalog zero damping fields (all five radiative;
10914 also Stark), and gf values are not yet the adopted AGSS21 feature values.
Finite smoke deltas therefore do not validate Asplund replication or justify
product admission. Broadening/source errors and collision-model differences
still need resolution.

### 10092 component and source findings

The existing published-line helper's +/-0.10-A window finds 10092.217 A
(J=3.5 -> 4.5, log gf=1.02). A wider read-only inspection of the held VALD
archive also finds 10092.095 A (J=2.5 -> 3.5, log gf=0.91): it is 0.105 A
from the rounded AGSS21 centroid, just outside that helper's window. Both have
zero damping fields. Do not scale one component to the whole feature or widen
all matching tolerances blindly. Full component completeness, allocations and
broadening remain unresolved. Evidence: `mgii_10092_components.json`.

The published [CDS table J/A+A/688/A72](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/688/A72)
was downloaded and inspected: 2646 records and 12 feature labels. Its Mg II
feature is 448.1 nm; none of the six AGSS21 Mg II red/NIR lines is included.
The exact hash and labels are recorded in `matsuno_cds_coverage.json`.

To close full-3D replication we still need the AGSS21 solar Mg I/Mg II model
and transfer inputs (including collision rates and broadening), or the authors'
validated line-specific synthesis/correction data for our selection, with their
atmosphere and abundance conventions. A 448.1-nm correction or a 1D GALAH result
cannot be substituted for those missing data. No author contact was sent.

All named JSON evidence is under `data/output/rya1250/`; the compact summary is
`mgii_validation_summary.json`. Culls and existing uncertainty-input holds
remain open. Publication and preparation apply have not run.

### Controlled turnover diagnosis

At fixed A(Mg)=7.55, perturbing only the 9218 line log gf by 0, +0.4 and
+0.8 dex gives LTE EWs 51.566, 74.539 and 103.361 mA, respectively. Unlike
the abundance-varied curve, this control rises monotonically and crosses the
65.033-mA NLTE target. This implicates the abundance-dependent background/EOS
response in the turnover, rather than a simple inability of the line to grow
with opacity. It does not isolate which internal background contribution is
responsible, or validate the catalog broadening. No gf perturbation was adopted,
and this experiment is not a replacement correction method.
Evidence: `mgii_fixed_background_control.json`.
