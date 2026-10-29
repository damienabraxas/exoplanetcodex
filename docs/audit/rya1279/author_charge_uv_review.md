# Sourced charge coefficients, UV controls and the surface boundary

RYA-1281 / parent RYA-1278, continuing from 4f01ebde. Same stacked draft PR #598;
all calculations below are bounded serial Mac diagnostics. Scientific admission
remains **HOLD**. Ryan's element order remains O, then C/N (CNO milestone), then
Fe, then another star. No adopted abundance or default collision mapping changed.

## Source tracing now has a numerical result

The author-linked [Barklem 2018 data](https://github.com/barklem/public-data/tree/cd2f54ed46687659017b343763e95478ed93653e/inelastic-O%2BH)
provide forward ion-pair and reverse mutual-neutralization (MN) coefficients.
ReadMe.simple/ReadMe.cds establish cm^3/s and matrix orientation (final, initial);
states.dat identifies the LS terms and ionic states. The paper's
[author preprint, Section 2](https://arxiv.org/html/1712.01166v2) identifies the
O+ + H- ionic channel. The ten 1000-10000 K tables, state list and two ReadMes
were acquired externally and verified against pinned Git blobs and SHA256.
barklem2018_oh is registered in the bibliography; full matrices remain external.
The repository carries BSD-3-Clause (Paul Barklem2020); its notice is retained
in docs/licenses/barklem-public-data-BSD-3-Clause.txt. The derived comparison
includes selected published coefficients, not the complete matrices or raw atom.
This primary evidence replaces the previous access-limited paper search.

An explicit LS/J mapping compares every held CH0 table: 24 channels x 10
temperatures = **240 entries**. The held coefficients closely follow

    q_raw ~= q_MN(term) * [g(O II)/g(O I J)]
                       * exp[-(E(O II)-E(O I J))/(kT)].

The ratio q_raw/q_candidate spans **0.99610079 to 1.00477394**, maximum absolute
relative difference **0.4774%**, using the unchanged O41 energies and weights.
No scale factor or energy shift was fitted. Source tables are rounded; exact
byte-equivalence or the generating transformation is not proved.
By contrast q_raw / published forward ion-pair coefficient spans **0.0003174
to 20.8333** across the entries. Treating q_raw as the published forward table
therefore is not established. The formula is strong numerical evidence of an
encoded reverse/Boltzmann relationship. It does not identify the modified MULTI
producer or license a silent correction to the default converter.

scripts/audit_bride_charge_contract.py records all 240 comparisons and a copy
of the verified acquisition receipt. The snapshot writer is not represented as
the original one-off acquisition tool. The raw atom and author matrices remain
unchanged. TERM_MAP explicitly records every manual term association and g sum.

## Independently sourced rate-amplitude control

A predeclared experiment replaces only CH0 with published MN term coefficients.
CE/CI/CH retain the preceding native-spline/source-direction control. Each MN
term coefficient is split among the represented final J states using
**g_J/sum(g_term)**, an explicit statistical-branching assumption. The local
H-minus collider is computed by ground-H Saha attachment; the inverse per-O
rate follows local LTE detailed balance. All source coefficients convert once
from cm^3/s to m^3/s. Native log TAUTSP interpolation remains in-domain.

This is sourced absolute MN plus assumed fine-structure branching and an LTE
inverse, not a claim to reproduce an independently supplied fine-structure
forward table. Tests sum all branches back to the author term rates, check
full-atom collision-only LTE closure and reject a non-LTE hydrogen reservoir.
The initial guard test incorrectly modified an LTE array aliased to nStar;
that test was corrected to supply a separate non-LTE population array. The
physical implementation and experiment inputs were unchanged by that correction.

A second predeclared control adds 47 native background lines with centers in
100-300 nm, fixed LTE/CRD: H 2, C 5, N 4, Mg 16, Al 7, Si 10, S 1, Fe 2.
[Departure comparison](../../../data/output/rya1281/author_charge_v1/departures.png) and
[flux comparison](../../../data/output/rya1281/author_charge_v1/flux.png) show both controls.
Moving those species to detailed-static also includes their native continua and
changes spectral sampling. It is **not** a pure line-only intervention or
complete UV blanketing. The constructor test verifies unchanged oxygen LTE
populations and LTE backgrounds. No C/N/Fe NLTE result is claimed.

Both full O41 runs retain MARCS 5750/4.5, xi=1 km/s, A(O)=8.70, 3 rays,
qscale=1, joint population/J/rate tolerance1e-6 and conservation1e-10.

| Control | Iterations | Seconds | Triplet fixed +/-0.5 A areas (mA) |
|---|---:|---:|---|
| Sourced MN + LTE inverse | 568 | 43.0445 | 68.7776 / 59.2304 / 46.2804 |
| Above + native LTE UV background | 396 | 250.2854 | 69.1461 / 59.5468 / 46.5247 |

Terminal physical residuals are 2.54e-15 and 3.76e-15 respectively; no invalid Ng update was rejected. All owned
science jobs finished. The attempted nice priority adjustment on the first run
was denied by the sandbox; the calculation ran one-thread at ordinary priority.
No Sirius science jobs or new dependency installation.

## Population agreement and what is still wrong

| Control | Upper b max difference, declared logtau[-3,0] | Upper b max, all depths | Upper b max, excluding first point |
|---|---:|---:|---:|
| Sourced MN pair | 0.005523 | 1.249476 | 0.063440 |
| Above + native LTE UV background | 0.055948 | 1.215358 | 0.198266 |

The photospheric subset is a diagnostic subset, **not** a contribution function.
The sourced MN pair reproduces the earlier effective-pair photospheric result
without fixing the surface. At logtau=-1.2621, lower/upper/O II b are
1.02710/0.84929/0.91040 versus reference 1.03018/0.84559/0.87156.
At the first point, logtau=-4.9177, sourced MN lower/upper b are
0.02515/0.03192 versus reference 0.98848/1.28139. With UV background they become
0.04910/0.06603, still far away. At the second point, the no-UV pair already
has lower/upper b=0.90807/1.14719. The first-point discontinuity deserves an
explicit boundary test; excluding it does not erase the 0.06344 upper difference
or the O II mismatch (all-level/depth maximum still 7.77491).
The small UV background control worsens the deeper agreement while approaching
the reference-replay fluxes. Spectral agreement remains insufficient for admission.

## Independent UV integration narrows the audit

scripts/audit_bride_photoionization.py uses the previously verified **raw-paired,
no-H-line fixed-population rate replay** (native-paired-rates-v2), not these new
MN/UV runs. It integrates 4*pi/h * integral[sigma(lambda)*Jnu/lambda d(lambda)]
using linear interpolation of the unchanged raw cross-sections. This is separate
from backend WENO and backend rate quadrature; their differences are reported.

For the lower-triplet 277.2477 nm continuum at the top point:

- Independent photoionization: **16245.8141 s^-1**.
- Backend: **16241.3588 s^-1**; independent/backend=**1.00027431**.
- Local-Planck integral: **251.4464 s^-1**; actual J-weighted rate / Planck=**64.6095**.
- At logtau=-1.2621, the corresponding enhancement is5.9233, with quadrature
  ratio1.00026937.

This rules against a large unit/quadrature error in this particular rate
calculation. A nonlocal superthermal field does not by itself prove missing
opacity. These values must not be attributed to the new UV-background run.
All 22 continua and three fixed-population modes are retained in the CSV.

## Strong next lead: the held native top boundary differs

Inspection of [pinned pymulti4d6aa3c transp_23.f](https://github.com/MingjieJian/pymulti/blob/4d6aa3c39d634b38da772e469a5692435eba1354/pymulti/files/multi/mul23/source/transp_23.f) finds an explicit finite
optical depth above the first grid point. TRANSP lines60-68 sets the ray optical
depth T=tau(1)*X(1)/mu. TRANF lines137-143 documents incoming radiation from a
constant-source cap: I_in=S(1)*(1-exp(-T)); lines177-200 implement its limiting
behavior. The analogous factors occur in the other solver branches.
File SHA256: d458213711b468899090ef8bbd108ad6a8387eb1d9657cd6df546966e0314eba.

Our constructor explicitly uses ZeroRadiation at the first geometric point.
Those conditions agree only when the overlying monochromatic optical depth is
negligible. Small continuum tau500 does not establish that for strong UV lines.
This is a verified difference in the held source family, **not** proof of which
boundary the exact grid producer used or proof that it explains the populations.
No cap, atmosphere extension or incoming radiation has been added in this phase.

Next: preserve the native MARCS top column/optical-depth metadata and construct
an explicit finite-depth-cap diagnostic from actual opacity/source values. Test
its zero-depth, optically thick and angular limits, then compare at identical
rates/atmosphere. Do not tune a cap thickness, incident temperature or UV factor
to match the reference. Exact producer opacity/run inputs remain a parallel need.

## Evidence and reproduction

Original immutable stores under CODEX_WORK_ROOT/bride_mac_checks:

- rya1281-charge-contract-v1/v2 preserved; **v3** is the final receipt snapshot.
- rya1281-author-charge-v1 and rya1281-uv-background-v1: two full solved cases.
- rya1281-author-analysis-v1 preserved; **v2** includes first-point exclusions.
- rya1281-photoionization-v1 preserved; **v2** includes source-control identity.

Author input: experiments.bride_barklem_oxygen in the path register;13 files,
pinned receipt SHA256 bcf3d889a969309a9c083aceaee4b129b1892de73cd9375db1b56b1969f1bdc0.
The small review package is data/results/rya1281/author_charge_v1 with figures
under data/output/rya1281/author_charge_v1. Raw author tables remain external.
Reproduce with the committed scripts' --help and the two predeclared plans;
use new immutable output directories. No default-path overwrite is supported.

Validation:167 focused science tests and3 installed iSpec tests pass. The
experimental UV case is refused by the fixed-physics grid CLI. All18new
artifacts/flags validate; the global generator guard retains exactly2436prior
branch diagnostics (0added/0removed). Source-page consistency and register
freshness pass. Standing library gaps remain. Current main is9d3ebbea; this
branch has not merged/rebased it. Details are in progress.md and the Linear handoff. No full-suite/CI-green claim, scientific admission, production
emission, merge or publication. The exact producer contract stays OPEN.
