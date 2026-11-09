# RYA-1281 physical direct-map acceleration qualification

Scientific HOLD. This follows the preserved [matched-column failures and operator
identity audit](solar_matched_columns_review.md). Its purpose is numerical
qualification before interpreting canonical-solar T/g effects, not selection of
an abundance, physical atom or source/reference model.

## Frozen-matrix conditioning

`audit_bride_solver_conditioning.py` reads the exact archived direct/MALI/direct
fixed-population audit. It changes no atmosphere, radiation, rate or population.
For each of56depths it inspects8predeclared mixtures of the physical R matrix and
approximate MALI P. Row-only and relative-population-coordinate systems use the
same largest-population conservation row. The frozen error operator is
I - P_c^-1 R_c, with conservation imposed in BOTH matrices, so the normalization
mode is not counted as a physical instability.

At pure MALI,52of56frozen systems have spectral radius>=1, maximum1194.558.
At0.75MALI,50are unstable, maximum2.98999. At0.5MALI, all are below1 but the
maximum0.998329is close to1. Physical-only systems have maximum2.77e-14,
as expected when the preconditioner equals the fixed physical operator.
This diagnoses the frozen-rate correction only. MALI is designed to precondition
radiation feedback; these eigenvalues are NOT derivatives or convergence proofs
for the full nonlinear coupled RT/SE map. No mixture is selected for production
or used to tune a spectrum.

Relative-population coordinates improve the worst MALI condition number from
1.32555e11to4885.46. Physical-only improves1.92425e10to1062.55. Conditioning
improves, but that coordinate change does not remove the frozen-rate unstable
modes. This separates scaling from the approximate-operator correction issue.

## Prospective direct-map Anderson control

Config `rya1281_solar_matched_anderson.json` retains the original common fractional
seed independently in all four cases, all physical inputs/owners/source receipts,
all56depths, original1e-6joint/physical threshold,1e-8nested radiation tolerance,
1500outer/100inner limits and600swhole-case budget. Native Ng is explicitly false,
replaced by log-fraction Anderson of the qualified physical direct population map:
memory5, delay10, least-squares rcond1e-10, maximum proposed log-population step0.05.
No MALI or prior failed seed is used. This changes a numerical algorithm only.

History stores actual direct-map input/output pairs. The unconstrained log proposal
is normalized per depth to the fixed oxygen total using log-sum-exp. Oversized,
nonfinite or nonpositive proposals are rejected in favor of the exact direct-map
output, with history reset. No population floor or clipping is used. Every
accepted population still faces the existing strict fresh physical C+transition-R,
conservation, radiation/cap and independent per-class reconstruction checks.
Rejection and acceleration metadata are recorded at each iteration.

Two discriminating tests show acceleration of a slowly contracting conservative
map while retaining its independent physical root, and exact fallback for an
oversized proposal. Together with matched-input and defect-correction controls,
6tests pass1.27s. These do not substitute for the full-oxygen endpoint checks.

All four dry Contexts pass. Actual native baseline passes3outer and reproduces the
archived direct baseline flux and continuum byte-for-value (maximum numeric
differences0/0), fresh independent residual3.956231536e-7. It stops before the
delay10acceleration threshold. Remaining full-column qualification is required.
The analysis endpoint contract allows only the explicitly declared Ng replacement;
all other source/physical and numerical controls remain exact, and both1e-6flux
parity screens remain unchanged. No failed or incomplete matrix enters analysis.

Current main remains7f023fd7afterfetch; the pinned ORCH source is unchanged.
Library audit remains212documents/127accounted/85unindexed/12offMac, nonzeroexit.
No new literature is cited, no source download or binary rebuild is performed.
All RT jobs remain Mac-only/serial/one-thread, no Sirius or CI workload.


## Minimal continuum-overlap control, prospective contract

The full source atom contains two O II levels: zero-based39is2p34S* and40is2p32D*.
Its4->41source continuum is the only continuum to the excited ion level; the
others use the ion ground. This motivates a discrimination test, not a causal
conclusion by itself. `audit_bride_overlapping_continua.py` uses a synthetic
three-level atom with one neutral state and two ion states, isothermal6000K,
thermal boundaries, fixed positive LTE populations, two simple hydrogenic
continua and ordinary CI/CE collisions. Controls enable only the ground-ion
continuum, only the excited-ion continuum, or both. It samples direct/MALI/direct
operators after fixed-population radiation convergence. There is no finite cap,
custom collision adapter, source atom, interpolation, Ng or population update.
The60swhole-audit guard and all attempted outputs remain. This tests whether the
operator mismatch persists in a minimal compiled-backend configuration; it does
not qualify the artificial atom as physical oxygen or validate the full solver.


## Actual Anderson and overlap results

Anderson baseline passes20.6097s/3outer with exact endpoint flux/continuum parity.
Temperature times out600.05098s/115complete outer iterations; last population
change0.00395699, old-radiation physical residual0.00191369, J0.00120828, O
conservation7.77e-16. The serial driver stops; no Anderson gravity/solar run is
started. This algorithm is unqualified for the matrix.93focused tests pass,
666existing warnings visible; tests alone did not establish full-O convergence.

The native-background synthetic-v1 fixture fails the100iteration radiation
check; no operator result is claimed. Its original stderr remains external.
A distinctv2fixture uses explicitly specified thermal absorption/emission and
Planck initialization, without weakening the radiation check or60sbudget.
All LTE single-ground/single-excited/overlap variants converge fixed-population
radiation in3iterations and preserve operator action at roundoff.

Prospectivev3uses a fixed non-LTE population pattern[0.9,1.4,0.6], normalized to
original O totals. Its physical residual near0.4is intentional: it is NOT a
claimed statistical equilibrium. Single continua retain action identity within
2.4e-16; combined continua show MALI mismatch1.57853e-6while direct remains
7.0e-17. No population is updated. Thus the mismatch can occur without the source
atom, custom collisions, finite cap or interpolation, under controlled non-LTE
state and radiative overlap. LTE-only controls miss this failure mode.

## Two isolated source-based numerical plugin controls

Both builds copy the held clean LW source d33058f8857acf28d9187e778b6f41e49d3aad95
into fresh external directories and compile small process-local iteration plugins.
The installed4,127,408byte backend and original source checkout remain untouched.
The builder records source/header hashes, patch/provider, backend, command and
plugin receipts. Dynamic loading checks the exact backend/plugin identity and
original5MiBguard; libraries are222648and222696bytes, no renaming/size bypass.
Default storage/formal-solver properties match the public LW plugin interface;
actual compiled toy and full-O controls, not ABI assumptions alone, qualify each
use. Neither plugin is installed or selected as a shared default.

The first guard disables the approximate operator only where continua to distinct
ion levels overlap. It restores the synthetic action identity(v4), but DOES NOT
change the full-source mismatch: source O continuum wavelength extents do not
match that simplified overlap condition. This unsuccessful intervention remains
preserved; the toy reproduction alone was not proof of the full-O cause/fix.

The second plugin uses per-transition preconditioning. For a transition i->j,
let K=PsiStar*Vij*Uji. Its offdiagonal update is
P_ij=Rji-n_i*K and P_ji=Rij-n_j*K, with the same quadrature weights; column
conservation supplies the diagonal. Multiplication by the old populations cancels
the two added terms exactly, preserving the physical transition action regardless
of other overlapping transitions. Formal transfer, opacity/emissivity, accumulated
physical Rij/Rji, all collisions and populations are unchanged by this operator
construction. This is a numerical control, not a new atomic model.

Syntheticv5non-LTE action identity is restored below2.9e-16. In the full-source
fixed-population auditv3, MALI action mismatch falls from0.00570313to9.74725e-15.
All three physical C+R arrays are EXACTLY equal to the original native audit;
physical residual stays3.95623e-7. The direct/MALI/direct return check and fixed
populations pass. This completes operator-identity qualification only; baseline
convergence, independent fresh physical rates and original endpoint parity remain
required before any solar-column or abundance interpretation.


The source-based per-transition control is running with the original native Ng
settings, original warm fractions and600sbudget. Its early unaccelerated steps
fall0.000668->0.000187->6.35e-5->2.16e-5->7.37e-6->2.53e-6, but larger updates
recur on the native acceleration cadence (e.g.iter9=0.000743). This is temporal
association, not proof that Ng alone causes the full convergence failure.
A future optional no-Ng control could isolate this cadence if convergence
remains problematic. It has not been run. The current source-based native-Ng
control continues under its original budget; no mid-run algorithm is changed.

The per-transition/native-Ng baseline PASSES197outer/398.33487s with fresh
independent physical residual6.810874767e-8 and conservation7.77e-16. Original
baseline endpoint spectrum differs1.45230e-7normalized and1.10134e-13relative
continuum, passing unchanged1e-6screens. Some rare-level populations differ
by up to0.0008769fractionally, so this remains a spectral numerical screen,
not a claim of all-level forward-error accuracy. The remaining T/g/solar cases
now use the same independently mapped common fractional start, not this result
as a sequential seed.
