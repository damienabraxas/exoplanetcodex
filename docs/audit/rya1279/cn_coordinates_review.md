# CN AIR/vacuum coordinate diagnostic and provider contract

2026-11-08, RYA-1281/1278, scientific HOLD. Starta2869fc9, worktree codex/rya1279,
branch ryandamienschmitt/rya-1279-1280-bride-oxygen-gates, stacked base4635c742.
Plan/codebc4bcee0; explicit AIR-domain guard/environment receiptfbb37c06.

## Controlled mapping checks

168unchanged native CN states: three sampled isotope lines x56depths. State report
and both CN/all CSVs byte-match committed accepted originals. Original source
newvoigt.f byte-matches pin21269948750e2177fe32c5df38e1d9bb81363374. Small
standalone caller reproduces native bsyn wavelength operations, including default
REAL assignment between steps and the rounded exported XL:

```
vt=(wave_air_A-centre_air_A)*1.d-8
vt=c*vt/xl**2
v=vt/doppler_width
```

2352native calls span14predeclared AIR coordinates per state. Python float64
native-map reconstruction agrees1.51402614e-7relative, below predeclared
128float32epsilon=1.52587890625e-5 arithmetic screen; zero offsets agree exactly.
This is the native-compatible linear AIR slope, not exact physical vacuum frequency.
Native red-positive v and conventional blue-positive frequency offsets give the
same symmetric Voigt profile; no velocity field/RV convention is validated here.

The kernel now accepts shared1D or explicit wavelength x depth2D offsets. New
`frequency_offsets` requires declared input medium and mode; incompatible frames,
unknown modes, invalid shapes/nonfinite states and undefined AIR wavelengths below
the canonical2000A boundary are refused. The latter guard was added after source
review, before any such input was used; no invalid UV data were measured.

## Four cases, no silent substitutions

All cases use identical populations/gf/native constants and the original4001point
AIR7763–7783Agrid at0.005A. Absolute extinctions are computed independently per
sampled line/depth, not summed into a spectrum.

1. Native-compatible linear AIR offsets, original width/a/stim.
2. Exact reciprocal AIR-frequency difference, original width/a/stim: isolates
   linearization plus native rounded-XL slope difference. AIR-frequency coordinates
   are a computational control, not physical vacuum frequencies.
3. Exact vacuum-frequency difference after canonical one-time AIR->vacuum conversion,
   original width/a/stim frozen: isolates the coordinate intervention.
4. Exact vacuum frequencies with width rescaled by native_lambda_cm/lambda_vac_cm,
   a inversely rescaled to preserve gamma, and stimulated emission recomputed at
   vacuum line-centre frequency with the same native constants. This coupled
   intervention is not an adopted correction or a physical model admission.

| Case | max absolute profile change / native peak | max pointwise relative change | max absolute relative finite AIR-area change |
|---|---:|---:|---:|
|Exact AIR frequency|7.50898e-6|0.00284408|6.91400e-8|
|Exact vacuum, frozen widths|0.000212397|0.00959677|0.000280483|
|Exact vacuum, width/a/stim rescaled|0.000262109|0.00339330|0.000267345|

Maxima span all168states. Pointwise maxima include weak profile wings; they do
not describe peak changes. Areas integrate extinction (m^-1) over sampled AIR
centre+/-0.5A, with units m^-1 A. These are **not observed EWs**, emergent flux,
abundance shifts, uncertainty bounds or a line-grid convergence demonstration.
Finite wavelength-area changes include coordinate Jacobian and window effects.
No reference-number tuning, favorable line/depth choice or post-result threshold.

## AIR/vacuum handling

Only canonical pipeline.wavelength_util is used, fingerprinted in the report.
AIR centres shift+2.13881553–2.13902989A for these three sampled lines. Original
AIR grid is retained for native controls; its vacuum counterpart is created once.
Canonical approximate inverse is evaluated separately as an audit, never used to
construct the native grid or to apply a second observation conversion. Its measured
roundtrip residual across20A is1.11673944e-5A. No iterative replacement or global
converter edit is introduced. Legacy TS-consumer AIR semantics still do not
independently establish primary-release molecular wavelengths.

No Kitt Peak/IAG observation, RV, telluric correction, continuum, source wavelength,
gf or abundance is changed or read as a target in this experiment.

## Concrete diagnostic provider input contract

The next implementation must make the following explicit before attaching opacity:

- Immutable source-state receipt, sampled transition/source identity, units and
  statistical/isotope convention. External native CN populations are external
  state, not a solution of the Bride chemistry. Never multiply isotope factors twice.
- Explicit vacuum-nm LW grid and one declared mapping policy; canonical conversion
  at the boundary. Legacy-compatible and exact-vacuum policies remain distinct.
  Depth-dependent widths require2Doffsets; no implicit broadcasting across lines.
- Matched atmosphere/depth ordering, T, density, electron/H normalization,
  composition and microturbulence receipts. CN pressure derives from the full
  native network; the current export does not by itself prove Bride EOS parity.
  Native/default-REAL versus model precision must be reported, never silently
  coerced into an arbitrary tolerance or used to substitute a solar column.
- Baseline background inventory before adding CN, to prevent counting an existing
  opacity component twice. With CN disabled, native background arrays must be
  unchanged. Enabled extinction/emissivity/scattering routing must be separately
  checked, including an explicit LTE Kirchhoff/Planck convention and constants.
- Initial provider/background-array tests with synthetic/disposable state and no
  RT/SE. A real O NLTE diagnostic follows only after matched state/identity checks;
  adding opacity to a canonical solar column using native5750/4.5populations does
  not satisfy those checks. Primary provenance and independent chemistry remain
  distinct scientific admission gates even if array plumbing works.

No provider class or new LTE emissivity is implemented in this phase. The mapping
and extinction kernel supply numerical pieces of that contract, not missing chemistry.

## Attempts, verification and reproduction

Bothv1andv2exit0. v1is preserved; source review added the undefined-AIR-input
refusal and explicit Python/NumPy/SciPy environment receipts, then freshv2reran.
All numeric output CSVs and native input/stdout are byte-identical between attempts;
only code/environment report metadata differs. No bound/physics adjustment.
Run bug ledger records the API validation gap as fixed, without an affected spectrum.

Project Mac interpreter Python3.12.13/NumPy2.2.6/SciPy1.18.0; OMP/OPENBLAS threads1.
Preflight10tests passed0.29s (guard repeat0.31s). Final33tests pass0.64s, exit0:
new frame/zero-centre/linear-limit/2Dshape tests, existing CN kernel/reporting/source
checks and canonical wavelength-scale tests. No actual LW provider/RT execution
claimed. Standalone native probe compile/runtime reported separately; no new peak
memory measurement or end-to-end performance claim. No Sirius/CI/install/download.

After source config/hosts/mac.env:

```
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /Users/ryanschmitt/codex/venv/bin/python scripts/check_bride_cn_coordinates.py --output <fresh-output>
/Users/ryanschmitt/codex/venv/bin/python -m pytest -q tests/test_bride_cn_coordinates_rya1281.py tests/test_bride_cn_diagnostic_rya1281.py tests/test_bride_cn_coefficients_rya1281.py tests/test_bride_cn_interface_rya1281.py tests/test_wavelength_scale_rya264.py
```

Registered external first/accepted directories retain unchanged native source copy,
private binary, logs and reports. Six compact outputs (report,672row summary,
2352row probe comparison, wrapper,input,stdout) are archived byte for byte fromv2.
Actual producer/invocation registered; no native source/binary committed or renamed.

Next: implement and test the explicit external-state provider contract above, first
at the background-array level. Primary line-list/frame/statistical provenance,
EOS/composition consistency, independent reference, canonical solar abundance and
adaptive3D/all-band goals remain OPEN. No adopted model/product/default, production,
publication or merge. Scientific HOLD remains.

Final guards: MIRROR/version/sequence/whitespace pass. Generator diagnostic set
exactly2436historical, zero added/removed. Acceptedv2compile0.50026s/caller0.26837s
exclude setup/grid evaluation. Existing global ledger RYA587blank science_impact
and tracker/gold debts remain; no global-green claim. Last library audit in this
continuing session196documents/127accounted/69unindexed/12off-Mac, incomplete.
No new literature used.
