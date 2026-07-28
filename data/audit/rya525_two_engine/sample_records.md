# RYA-525 — sample per-element two-engine records (constructed demo)

_Thresholds from config.TWO_ENGINE: saturation_knee=100.0 mA · cross_engine_mix_gate=0.1 dex · synth_chi2_gate=10.0._


## Si I (clean-weak → Engine-A)
- **reported value = 7.508 ± 0.031**  (n=2, engines used: engineA_1dnlte)
- diagnostic: Engine-A=7.510 · Engine-B=7.555 · mean cross-engine Δ(B−A)=+0.045 dex
- cross-engine mix=False · **mix_flagged=False**
    - 5701 Si I: **engineA_1dnlte** (7.500) [clean-weak] — clean-weak line → 1D-NLTE (cleanest for weak lines); rejected engineB_synth=7.560
    - 5772 Si I: **engineA_1dnlte** (7.520) [clean-weak] — clean-weak line → 1D-NLTE (cleanest for weak lines); rejected engineB_synth=7.550

## Mn I (HFS → Engine-B)
- **reported value = 5.424 ± 0.031**  (n=2, engines used: engineB_synth)
- diagnostic: Engine-A=5.290 · Engine-B=5.425 · mean cross-engine Δ(B−A)=+0.135 dex
- cross-engine mix=False · **mix_flagged=False**
    - 6013 Mn I: **engineB_synth** (5.420) [hard] — hard line (blend/saturation/HFS) → synthesis; rejected engineA_1dnlte=5.300
    - 6021 Mn I: **engineB_synth** (5.430) [hard] — hard line (blend/saturation/HFS) → synthesis; rejected engineA_1dnlte=5.280

## Ca I (mixed, agree → combined)
- **reported value = 6.315 ± 0.035**  (n=2, engines used: engineA_1dnlte, engineB_synth)
- diagnostic: Engine-A=6.305 · Engine-B=6.325 · mean cross-engine Δ(B−A)=+0.020 dex
- cross-engine mix=True · **mix_flagged=False**
    - 6122 Ca I: **engineA_1dnlte** (6.300) [clean-weak] — clean-weak line → 1D-NLTE (cleanest for weak lines); rejected engineB_synth=6.320
    - 6162 Ca I: **engineB_synth** (6.330) [hard] — hard line (blend/saturation/HFS) → synthesis; rejected engineA_1dnlte=6.310

## Ti I (mixed, disagree → FLAGGED)
- **reported value = 4.980 ± 0.035**  (n=2, engines used: engineA_1dnlte, engineB_synth)
- diagnostic: Engine-A=4.905 · Engine-B=5.055 · mean cross-engine Δ(B−A)=+0.150 dex
- cross-engine mix=True · **mix_flagged=True**  ⚠️ ADJUDICATE (disagreement > gate; not a silent mean)
    - 5689 Ti I: **engineA_1dnlte** (4.900) [clean-weak] — clean-weak line → 1D-NLTE (cleanest for weak lines); rejected engineB_synth=5.050
    - 5648 Ti I: **engineB_synth** (5.060) [hard] — hard line (blend/saturation/HFS) → synthesis; rejected engineA_1dnlte=4.910
