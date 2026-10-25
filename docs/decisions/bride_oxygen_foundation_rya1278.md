# Bride oxygen-first foundation — RYA-1278

**Engineering GO; physical first-light HOLD.** This is a review proposal and an
input-validation implementation, not certification of an oxygen atom or solver.
No RT/SE, abundance fitting, production run, publication or merge was performed.

Audit date: 2026-10-25. Base `0ade2170ebeda715f7b23ee8b0f57efdae5d0f3a`
(main, includes RYA-1274 PR #596); branch
`ryandamienschmitt/rya-1278-bride-oxygen-foundation`; isolated worktree
`codex/rya1278` below the workstation home. Other Mg/Na worktrees were untouched.
The issue descriptions and chronological comments were reconciled with code,
including all 118 RYA-1214 comments and RYA-1220's later corrections.

## What the earlier tickets establish

| Ticket / current state | Established evidence | Proof still missing |
|---|---|---|
| RYA-1008 / In Progress | Frontier investigation and competing engine routes; historical cost estimates | No certified full-3D Al or Bride O result |
| RYA-1012 / Done | `voronoirt_rung1_rya1012.md`: searchlight within 9 ULP; homogeneous tiled-FALC continuum comparison within 0.5% | Structured solar atmosphere / independent full-3D NLTE validation; source cube was unavailable |
| RYA-1118 / Done | Lightweaver Ca II FALC example converged; reported 1.70 s/column for that example. Corrected x-fastest STAGGER reader at `0a47f754`; exact non-square ramp tests | Planned STAGGER LTE Fe reproduction and Judge2017 comparison are not demonstrated by these artifacts; no measured O runtime |
| RYA-1119 / Done | `pipeline/multi_atom.py`, merged ancestry `1c73a879`: MULTI→native LW conversion, unit/format/charge-exchange tests | Restricted Fe atom completeness, independently validated populations/spectra, production 1.5D Fe reproduction. Mg 32.15 s/column was a proxy, not oxygen timing |
| RYA-1137 / Backlog | An explicit first-light atom decision was requested; no comments | Converter prose claiming this decision was established is not a completed scientific decision |
| RYA-1032 / Done; RYA-1099 / Backlog | Existing Gerber mean-3D route and paired estimator; unresolved Fe model/line discrepancy | No proof that a horizontal mean or independent columns reproduce full-3D transfer |

Existing reader tests detect the historical permutation that finite/range/byte-count
checks missed. Do not reinstate the original row-major fixture assumption. The
current reader returns SI `(z,y,x)` arrays, decreasing geometric height and already
processed T/ne from `staggtom3d`; raw hydro EOS conversion is not established as owed.
Its Saha split is not a multi-level hydrogen NLTE solution. It uses `nan_to_num`;
new boundary validation rejects invalid fields but does not certify that approximation.

## Ranked findings

1. **Critical — historical gravity guidance has a unit error.** RYA-1118/1137
   describe passing conventional solar cgs logg directly to LW. Installed 0.16.1
   `Atmosphere.make_1d` documents logG in **m/s²**. Canonical stellar logg is cgs
   (`config/constants.py`, `config/stars.yaml`). Passing it raw would multiply
   gravity by 100. The adapter resolves the canonical value and subtracts
   `log10(100)`; it never changes the stellar parameter. The library default is
   approximately solar in SI, but using any implicit default loses star identity.
   No affected committed Bride science product was found. See bug ledger.
2. **Critical — no admitted resolved-triplet atom or independent O reference.**
   Counts and class names alone do not establish rate completeness. See below.
3. **Warning — smoke/converter tests are not independent RT/SE validation.**
   Finite populations, convergence or agreement with a familiar abundance do not
   close this gap. Numerical benchmark gates remain open.
4. **Warning — full-3D capability is unproved.** LW advertises 1D/2D transfer;
   `Layout.make_3d` is not evidence of horizontally coupled 3D RT+SE. Model roster
   row 8 still names the earlier VoronoiRT route. Its source history is retained;
   metadata now identifies the adopted experimental LW/STAGGER foundation and
   states that full 3D is pending. No treatment token is introduced.
5. **Warning — snapshot identity and resources limit comparisons.** The held
   cube is not demonstrated to be the exact published oxygen benchmark cube.
   Sirius was inspected read-only: load 14.18/13.95/13.61, 20 users; no jobs started.
   Historical `/srv/codex/grids-overflow/atmospheres/stagger_cubes/t5777g44m00`
   exists; the tested `/mnt/codex-ext/...` alias does not. Resolve host paths before
   any later job; do not reinterpret a missing alias as missing science data.
6. **Suggestion — retain independent gate artifacts and diagnostics.** Integrate
   into ORCH only after validation; never infer admission from process completion.

## Held inputs and oxygen atom decision

[`held_assets.json`](../audit/rya1278/held_assets.json) is generated read-only by
`scripts/audit_bride_foundation.py`. It records snapshot and installed-source hashes,
API signature, dependencies, canonical gravity, atomic structure/rate-class counts
and canonical triplet records. `engines.bride_eval` discovers the existing holding;
its presence does not authorize compute. No package was installed or changed.

The held snapshot `t5777g44m0005_00091_80x80x240` has shape `(240,80,80)` and
finite arrays passing the new boundary checks. Its info header identifies the
March 2024 `staggtom3d` processing, EOS table, sampling and Rosseland cuts.
The horizontal coordinate extent is about 7.9 Mm. Original acquisition/reuse
receipt and exact published snapshot equivalence remain prerequisites, not inferred
from the stem. Do not assert a global monotonic temperature in convective columns.

Installed environment: Python 3.12.13, LW 0.16.1; full versions and hashes in JSON.
The held source clone is `d33058f8857acf28d9187e778b6f41e49d3aad95` (clean).
LW code is MIT, copyright Christopher M J Osborne; preserve attribution.
The atom's original scientific data rights/provenance need their own receipt.
CRTAF is absent; the existing converter produces LW-native objects, not CRTAF.

| Bundled native factory | Levels / lines / continua / collisions | Triplet-region audit / verdict |
|---|---|---|
| O_atom, O_no_777_atom | 23 / 8 / 11 / 66 | No transition in 776–779 vacuum nm: HOLD |
| OI_fine_atom, OI_fine_crd_atom | 14 / 18 / 13 / 100 | One collapsed transition: HOLD |
| OI_ord_atom | 14 / 16 / 13 / 100 | One collapsed transition: HOLD |

The OI variants have one line at 777.570512 vacuum nm, f=0.9781, statistical
weights 5→15, labels 3s5So→3p5Pe. This is not three resolved components.
They contain CE/CI, proton collisions and H/proton charge exchange, but no
demonstrated benchmark neutral-H excitation parity. Radiative, Unsold and Stark
broadening are inventoried in JSON; this is not a damping validation. Continuum
counts do not validate photoionization cross-sections or background opacities.

**Preferred reference physics:** the reduced atom in Amarsi et al. 2018,
[full text](https://arxiv.org/html/1803.10531v2), §§2.1–2.4, Figs 3–5.
It has 50 levels including three O II levels, resolved relevant O I fine
structure, TOPbase photoionization, and a documented radiative/broadening chain.
Its neutral-H excitation combines LCAO and Impulse rates; electron excitation uses
BSR calculations. The study uses balder full-3D transfer on seven resampled
80×80×220 solar snapshots and centre-to-limb intensity observations. The reported
0.02 dex agreement is an observational result, not a general solver tolerance.

No reusable artifact with complete atom/rate license and independent numerical
baseline has been admitted here. RYA-1279 must acquire/verify one; do not assemble
missing physics from scratch. The older held Amarsi2016 atom (22 O I levels plus
O II ground, 43 radiative transitions; §2.3, Table 1) uses older collision physics,
including Drawin with SH=1. It could define a separately labelled exact replication
if its artifact is obtained; do not mix it with 2018 rates and claim replication.
An optional synthetic O-like atom can test numerics but cannot close G1.

Canonical triplet air wavelengths are 7771.944/7774.166/7775.388 Å;
`gf_108794/5/6`, physical IDs and adopted values are preserved in JSON. They carry
NIST grade A with tier NIST-C+, **not primary laboratory gf**. RYA-1214 adopted
NIST and left its approximately 0.016 dex offset from Tachiev open to Ryan.
Do not change gf to make the atom agree. Air/vacuum conversion must use the
canonical wavelength utility and physical identities, never a free matching tolerance.
RYA-1129 CROSSMATCH_REVIEW is an earlier intake checkpoint, not today's entire
admission state: RYA-1214 adjudicated gf and approved species eligibility, while
product gates and the O source decision remain. RYA-1220's later N artifacts
remain HOLD for explicit covariance/model/identity terms; no new debt is cleared.

## Ordered executable gates

| Gate | Objective evidence and failure rule | Current state / owner |
|---|---|---|
| G0 input boundary | Exact non-square byte-layout ramp; units and axis ordering; positive finite T/ne and total H; detached indexed columns; explicit boundaries, vturb and canonical SI gravity. Any silent repair/default fails | Implemented input controls; physical atmosphere/EOS parity remains open |
| G1 atom/reference | All three component identities, complete rates/units/provenance/license, background opacity and exact independent reference setup. Missing or mismatched physics fails | HOLD; RYA-1279, related RYA-1137/1214 |
| G2 formal transfer | Actual backend: vacuum Iout=Iin; constant-source slab Iout=Iin exp(-tau)+S(1-exp(-tau)); optically thin limit; both ray directions and explicit boundaries. Compare against analytic arrays and independent published code data | Open; RYA-1280, synthetic implementation can proceed |
| G3 SE/LTE | Nonnegative populations, total abundance conservation, normalized net rate residuals; detailed balance in thermal radiation; collision-dominated LTE limit; alternate initial guesses; spectral and population convergence jointly | Open; RYA-1280 |
| G4 one column | G1–G3 pass, exact immutable manifest, Ryan approves real run. Reproduce independent matched 1D/column spectra/populations before any observational inference; record failure histories and measured resources | HOLD; RYA-1281 blocked by 1279/1280 |
| G5 1.5D | Homogeneous tiled atmosphere reproduces its 1D solution; heterogeneous columns, area/time weights, rays, depths/frequencies and representative subsets converge. Horizontal velocities projected for each ray, not silently discarded | HOLD; RYA-1282 blocked by 1281 |
| G6 full 3D | Separate genuinely coupled RT+SE implementation, periodic horizontal/open vertical boundaries as appropriate, homogeneous horizontal limit and independent structured-atmosphere comparison | HOLD; separate milestone within RYA-1282 |
| G7 indicators/admission | Matched geometry/line set; verified exact-window tellurics; full absolute/differential uncertainty with measured covariance; model-domain and publication gates | Existing RYA-1214/1220/587 own this work; no duplicate campaign |

For algebraic/exact copy controls, tests use exact equality; the gravity identity
uses floating-point approximation. For real solver G2/G3, predeclare a versioned
error budget from machine precision, analytic discretization error and independent
reference uncertainty **before** examining results. Require refinement to stay below
that allocated budget, rate/population residuals below their specified norms, and
no material dependence on initialization. No generic dex cutoff is invented here.
Missing reference errors or tolerance justification means HOLD, not an automatic
pass. Include deliberately wrong units, boundary, ordering and rate normalization
controls proving tests fail. Refinement alone is necessary but not independent proof.

## First-light plan and budget

The runnable deliverables today are tests and a read-only inventory, not an O solver.
Resolve the existing evaluation interpreter from the `engines.bride_eval` holding;
from the checkout, using that interpreter:

```sh
PYTHONPATH=. "$BRIDE_PYTHON" -m pytest -q -p no:cacheprovider tests/test_bride_foundation_rya1278.py tests/test_multi3d_read_rya1118.py tests/test_multi_atom_rya1119.py
PYTHONPATH=. "$BRIDE_PYTHON" scripts/audit_bride_foundation.py --snapshot-stem t5777g44m0005_00091_80x80x240 --output /tmp/rya1278-held-assets-new.json
```

`BRIDE_PYTHON` must point to the verified held interpreter; the CLI fails for missing
assets and refuses to overwrite an audit file. The committed audit records the base
commit and hashes of its generator/validator, since these were uncommitted at capture.

RYA-1281 must implement the future solver CLI and dry-run before seeking approval.
Its manifest must freeze: code/env/source hashes, atom/rates/background opacity,
snapshot and original provenance, star parameter source, exact ix/iy and all SI
arrays, EOS/H treatment, velocity projection, explicit vturb policy, boundaries,
fixed abundance as an input (not a fitted target), observable/mu/azimuth quadrature,
wavelength frame, frequency/depth grids, iteration/residual norms, thread count,
resource ceiling and immutable output destination. No implicit solver fallback.
The first experiment is a controlled vertical intensity column. It is not disk flux,
a solar abundance determination or a centre-to-limb validation. An independently
computed reference for that identical column must exist before comparison.

Versioned outputs must include manifest, spectra/continuum, populations, rates,
conservation/LTE diagnostics, full iteration/residual histories, exit disposition,
setup/iteration/synthesis wall and CPU times, memory/threads and host load.
Measured O timing is **null** today. After approval, estimate serial subset cost
from measured setup plus actual per-column costs over representative extremes;
include convergence failures, rays/frequencies, snapshot count and measured
parallel efficiency. Bound both wall time and memory before expansion. The old
Ca/Mg/Fe estimates cannot furnish this measurement or authorize a full cube.

## Ranked CNO UV/IR feasibility matrix

Ranking is feasibility for a validated Bride extension, not product acceptance.
Coverage below is registry/manifest evidence, **not a new pixel-level coverage or
SNR measurement**. All candidate windows require actual pixels, gaps, resolution,
continuum and residual QA. SNR is unmeasured in this audit for every row; no
instrument headline SNR substitutes for it. Instrument ranges/resolution come from
`data/catalog/instrument_catalog.csv`; holding conditioning from the manifest registry.

| Rank / indicator | Atomic/molecular data and physics | Holdings / instrumental feasibility | Dominant gates and uncertainty |
|---|---|---|---|
| 1 O I 777 triplet, 777.2–777.5 nm air | Canonical NIST-evaluated gf; preferred 2018 atom above, artifact HOLD; optical/red atomic NLTE | KPNO corrected flux atlases, IAG; KPNO R≈300k–500k, IAG≈700k. HARPS ends near 691 nm and cannot supply it. ESPRESSO route evidence exists in RYA-1220, exact holding receipt still needed | Geometry: disk-centre/CLV reference is not disk-integrated flux. Exact telluric windows, continuum/blends, collisions, gf decision, snapshot and numerical errors |
| 2 C I optical anchors 505.2/538.0 nm, then near-IR ~1.07 µm | Existing canonical/multiband intake and Amarsi C grids are prior work, not a validated Bride atom. Audit transition/rate reach separately | HARPS/KPNO/IAG for anchors; IAG and CRIRES+ Y/wide-Y candidates near 1.07 µm, detector edges/gaps must be checked | Optical control before wavelength expansion; gf/blends, model/parameter covariance; IR tellurics and background opacity |
| 3 N I 10108.90 Å (~1.01 µm) plus selected optical N I controls | Amarsi2020 Table 1 and RYA-1220 method. 10108.90 remains a noncanonical overlay: no matching canonical row in this audit; identity/gf/line-set admission first | IAG nominal reach and corrected CRIRES+ Y-wide (9800–10796 Å) are candidates; narrow Y starts 10280 Å and misses it. Complete registered line-window coverage not established | Weak high-excitation line, CN blends, continuum/EW, exact tellurics and N atom applicability. Do not confuse this with H/K |
| 4 Longer O I 844.6 nm and audited C I/O I J/H/K lines | RYA-1214 inventories candidate atomic lines. O8446 fluorescence/pumping and blends require explicit atom/background reach; no automatic transfer from O777 | KPNO/IAG for 844.6; CRIRES+ J/K R≈50k–100k manifests in `data/audit/rya1219_crires_jk/`; H refreshed by RYA-1232 | Settings are not continuous coverage. Individual line residual qualification, damping/blends, IR continuum/opacity, native atom rates and uncertainty |
| 5 CN A-X 1087.5–1320.8 nm; CO overtone ~2.3 µm; OH H-band | Molecular line lists/provenance already owned by RYA-1130/1214/1220. Molecular LTE/3D equilibrium is distinct from atomic NLTE; no molecular NLTE implementation demonstrated | IAG/KPNO CN and corrected CRIRES+ J CN/K CO exist; H OH candidates. Compare shared transitions only: IAG CN and CRIRES J need not overlap | C/N/O equilibrium covariance, isotopes, dissociation/partition data, strong tellurics especially K, temperature/3D effects, correlated pixels and holding repeatability |
| 6 Near-UV CH ~387–432, CN ~388, NH ~336, OH ~306–330 nm; atomic FUV C/N/O separately | Existing TS molecular route; NH336 is not the mid-IR NH rotational/fundamental diagnostic in Amarsi2021. FUV atomic resonance lines require chromosphere/interstellar/model-domain audit | UVES/KPNO mode/holding dependent; HST STIS nominal 114–1027 nm and R≈500–114k is stellar UV opportunity, not a held solar CLV benchmark | Missing UV opacity/pseudo-continuum, blends, chromosphere, flux calibration and molecular data. RYA-1220 UV molecular covariance HOLD and HST grid/domain gates remain |

Corrected holding entries are not blanket acceptance: new CRIRES H points to
`data/audit/rya1232_crires_h/corrected_products_manifest.csv`; J/K to
`data/audit/rya1219_crires_jk/corrected_products_manifest.csv`. Never use raw Vesta
IDPs as corrected data. Delbouille/Liège disk-centre intensity is held but has
unresolved conditioning here; it is not silently substituted for flux or CLV data.
No IR abundance without verified correction and residual evidence at each accepted
line/window. Molecular NLTE and unavailable UV model cells stay unavailable.

All future products must carry the RYA-587 full uncertainty contract: observational
and correlated-fit error, continuum, gf/rates, stellar-parameter Jacobians/covariance,
blend abundances, CNO molecular coupling, telluric residuals, model/dimensionality,
numerics and holding repeatability. For target stars, preserve paired Sun/target
Jacobians and shared gf/model covariance; do not quadrature two aggregate errors.
Solar work defines its zero point; no artificial differential to a literature number.
No averaging discordant families to hide model deficiencies.

## Verification and handoff

Before the final constructor integration check, project environment: 36 passed,
5 skipped (LW absent); held Bride environment: 41 passed, no skips. These include
existing reader/converter tests and new synthetic boundary tests; they are not RT/SE
validation. The final validation receipt in `docs/audit/rya1278/validation.md` records
updated counts and exact commands. No frozen gold, FITS, Fe/Al result or adopted
science parameter was modified.

Library resolver succeeded. Library audit remains nonzero for 40 unindexed documents
and 12 unresolved declarations with candidates present, already owned by RYA-1277;
this is not evidence of 12 missing papers. Held primary PDFs inspected include
Amarsi2016 (stv2608, §§2.2–2.4), AGSS21, Amarsi2020 nitrogen Table 1,
Amarsi2021 molecular methods, and Bergemann/Hoppe2025 p35. The latter's favourable
Mn 1.5D example is not a universal oxygen/full-3D equivalence result. New online
primary sources are indexed in the bibliography; no duplicate PDFs acquired.
Linear shared-agent-context endpoint was advertised but returned Tool not found;
canonical repo/issue evidence was used and this connector limitation is preserved.

Next independent work: RYA-1279 (atom/reference receipts) and RYA-1280 (synthetic
actual-backend tests). Their completion does not authorize a real O run. Then
RYA-1281 prepares a concrete dry-run manifest for Ryan, and RYA-1282 owns the
later dimensionality gates. RYA-1137 retains its broader Fe atom decision; the O-first
foundation does not delete, retune or certify prior Fe/Al work.

Other inspected primary references:
[Lightweaver paper](https://arxiv.org/html/2107.00475v1),
[upstream implementation/license](https://github.com/Goobley/Lightweaver),
[API documentation](https://goobley.github.io/Lightweaver/lightweaver-api.html).
Current online API names differ from the inspected pinned camelCase API; the adapter
is explicitly version-scoped, with no silent compatibility/default fallback.
