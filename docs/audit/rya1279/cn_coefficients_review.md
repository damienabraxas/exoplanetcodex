# Native CN coefficient reporting check

2026-11-08; RYA-1281/1278; scientific HOLD. Branch
ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, worktree codex/rya1279,
stacked base4635c742, starting0df65b4d. Predeclared plan/code67418fef.

## What ran and what passed

An isolated copy of pinned native source21269948750e2177fe32c5df38e1d9bb81363374
gets a reporting-only WRITE block after native population/plez calculation. No
native physical equation, constant, source table, input or installed binary changes.
Original and reporting binaries each run off/CN/all, six serial Mac programs,
using accepted original inputs with only RESULTFILE paths changed. All six native
numeric spectrum files are byte-identical to their corresponding accepted spectra.
Thus reporting preserves the full20A spectrum in these three cases, not just the
three sampled lines. Native output precision is unchanged, so byte identity does
not bound changes below serialization precision or prove all possible inputs.

Deterministic source selection: nearest source row to AIR7771.944A in each of the
three delivered CN files, tie by original row order. No selection based on strengths
or agreement. All56depths exported for each,168records in CN and168in all.
Those two state CSVs are byte-identical; off exports zero CN records as expected.

| Isotope file | AIR A | source line | excitation eV | loggf |
|---|---:|---:|---:|---:|
|12C14N_750-800.bsyn|7771.957|3043|2.278|−3.648|
|13C14N_750-800.bsyn|7771.980|2544|2.088|−10.857|
|12C15N_750-800.bsyn|7772.751|102|1.011|−3.363|

The weak13C14Nrow remains included. These three are computational sampling, not
three observed CN identifications or a validated full line-list coefficient sweep.
The source-level record precedes native profile/weak-line inclusion decisions.

Native fields include T/Pg/Pe/rho, full-network CN pressure, Q, isotope fractions,
isotope mass, total/lower populations per gram, gf, constants, stimulated-emission
factor, Doppler width, standard opacity, plez, damping a and pre-Voigt amplitude.
Actual CN pressure spans8.54027e-7–3.93016e-4dyn/cm2 over the56depths; Q6430.37793–
51498.5625. Masses4.31729e-23–4.48388e-23g distinguish sampled isotopologues.
CN pressure is exported from the native equilibrium result; no independent EOS
solver validates it here, and no undepleted C-times-N approximation substitutes it.

## Separate arithmetic reconstruction

Python float64 reconstruction takes exported native physical state and source
coefficients, not native intermediate results, to compute:

- CN per gram = pressure/(k*T*rho) times one fraction per constituent;
- source partition polynomial, then Boltzmann lower populations with the native
 5040/0.43429 convention and g=1;
- stimulated emission via stable expm1, gf from loggf, sqrt(pi)*e^2/(m_e*c),
  plez and pre-Voigt amplitude.

Prospective bound128*float32epsilon=1.52587890625e-5 was committed before execution.
It screens accumulated default-REAL arithmetic; it is not a mathematically certified
error bound for every operation, a physical uncertainty or a spectrum/abundance budget.
No adjustment after seeing outcomes. Maximum relative differences (CN and all equal):

| Quantity | maximum relative difference |
|---|---:|
|isotope-weighted total/ntt|1.16744e-7|
|Q|0|
|lower population|7.31865e-7|
|stimulated emission|5.89566e-8|
|gf|3.22538e-8|
|constant|5.58342e-8|
|calf|7.77298e-8|
|plez|7.93670e-7|
|pre-Voigt amplitude / standard opacity|8.10020e-7|

Q now agrees with actual compiled native output at these56temperatures; the prior
source-only projection limitation is closed in this bounded runtime sense. Shared
source coefficients and exported native state make this same-source interoperability
evidence, not an independent primary partition/line-opacity/EOS benchmark.

The exported amplitude is calf*plez, normalized by standard opacity, before Voigt
profile evaluation. No claim of a complete monochromatic absorption coefficient,
profile normalization, SI conversion or emergent-flux reconstruction follows.
Native LTE source function, scattering/continuum and chemistry remain unchanged.
No hand-added nuclear-spin or statistical factor; such physical conventions are
still unqualified against primary molecular releases.

## Attempts and verification

Self-review caught a fixed-form continuation missing a newline before execution;
fix86631181 adds boundary regression coverage. First full attemptv1 then completed
the build and six spectrum checks but failed report serialization: NumPy float32
prospective bound was not JSON serializable. No accepted report was written.
All raw runs/state CSVs/logs remain intact. Fixf77853a0casts the bound to Python
float; freshv2reruns the unchanged plan. No physics, bound, inputs or masks changed.
Both attempts' CN/all state CSVs match byte for byte. Bug-ledger dispositions appended.

Acceptedv2exit0:139protected-input checks before/after (includes repeated authorities,
not139distinct files). Compile7.39654s; original/reporting off0.63926/0.78485s,
CN0.64112/0.65369s, all0.63378/0.62714s. Native timing excludes setup/preflight.
No new per-program peak memory measurement; earlier26MB measurements remain scoped
to their two original repeat programs. No Sirius, CI, installation or source download.

Preflight2tests pass0.09s. Final16tests pass1.34s, exit0, covering exact source
insertion boundaries, malformed/nonfinite reporting refusal and existing partition,
isotope, build/input and opacity-reader contracts. Runtime spectrum/state/139input
checks are additional evidence, not extra independent tests.

Reproduce after source config/hosts/mac.env with OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1:

```
/Users/ryanschmitt/codex/venv/bin/python scripts/probe_bride_cn_coefficients.py --output <fresh-output>
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q tests/test_bride_cn_coefficients_rya1281.py tests/test_bride_cn_interface_rya1281.py tests/test_bride_native_lte_build_rya1281.py tests/test_bride_native_molecules_input_rya1281.py tests/test_bride_ts_budget_rya1281.py
```

External registeredv1/v2directories retain private source/build/native spectra and
stdout/stderr. Committed10compact derived files (report, insertion block,2state
CSVs,6inputs) byte-matchv2; binaries/native source copies/full padded logs stay
external. Actual producer and invocation registered per file.

## What comes next

Extend the matched-state check through native Voigt/profile normalization and
explicit cgs-to-SI conversion into a diagnostic Bride opacity provider. Retain
native CN pressure as an explicitly external state input initially; this does not
establish Bride chemistry parity or authorize production use. Require independent
population/partition/isotope conventions and full-network EOS/abundance consistency
before a coherent molecular model is admitted. Broader line coverage, canonical
solar rather than5750/4.5state and fresh O NLTE response remain separate checks.
Primary provenance, independent model reference, solar abundance, adaptive3D and
all-band gates stay OPEN. No adopted product/default, production, publication or merge.

Final repository checks: MIRROR/version/sequence/whitespace pass. Generator
diagnostic set exactly2436historical, zero added/removed. Existing global ledger
RYA587blank science_impact and tracker/gold RYA777/669 debts remain; no global
green claim. Latest library audit in this continuing session remains196documents/
127accounted/69unindexed/12off-Mac, incomplete. No new literature used.
