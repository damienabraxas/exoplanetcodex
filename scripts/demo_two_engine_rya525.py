"""Demo: sample per-element two-engine records (RYA-525).

Constructed lines only (encode-don't-tune — NOT real solar data). Shows each
selection branch + the element record the verdict will carry. Run:
    python3 scripts/demo_two_engine_rya525.py
"""
from pipeline import _runtime as _rt  # noqa: F401  force-fork + BLAS pins (RYA-514)
from pipeline.engine_selection import LineEngines, select_element, ENGINE_A
from config.constants import TWO_ENGINE

KNEE = TWO_ENGINE['saturation_knee_mA']
G = TWO_ENGINE['cross_engine_mix_gate']


def L(wl, sp, av, ae, bv, be, ew, blend=False, pc=False, chi2=1.0, hull=True):
    return LineEngines(wl, sp, av, ae, hull, bv, be, chi2, ew, blend, pc)


CASES = {
    # clean-weak metal → Engine-A wins every line
    'Si I (clean-weak → Engine-A)': [
        L(5701, 'Si I', 7.50, 0.04, 7.56, 0.05, ew=22),
        L(5772, 'Si I', 7.52, 0.05, 7.55, 0.06, ew=30),
    ],
    # HFS element → every line HARD → Engine-B (synthesis) wins
    'Mn I (HFS → Engine-B)': [
        L(6013, 'Mn I', 5.30, 0.09, 5.42, 0.04, ew=45),
        L(6021, 'Mn I', 5.28, 0.10, 5.43, 0.05, ew=52),
    ],
    # mixed winners, engines AGREE within the gate → combined, not flagged
    'Ca I (mixed, agree → combined)': [
        L(6122, 'Ca I', 6.30, 0.05, 6.32, 0.05, ew=28),                 # clean-weak → A
        L(6162, 'Ca I', 6.31, 0.06, 6.33, 0.05, ew=KNEE + 20, blend=True),  # hard → B
    ],
    # mixed winners, engines DISAGREE beyond the gate (Ti lesson) → FLAGGED for adjudication
    'Ti I (mixed, disagree → FLAGGED)': [
        L(5689, 'Ti I', 4.90, 0.05, 4.90 + G + 0.05, 0.05, ew=25),               # clean-weak → A
        L(5648, 'Ti I', 4.91, 0.05, 4.91 + G + 0.05, 0.05, ew=KNEE + 15, blend=True),  # hard → B
    ],
}

lines_out = ['# RYA-525 — sample per-element two-engine records (constructed demo)\n',
             f'_Thresholds from config.TWO_ENGINE: saturation_knee={KNEE} mA · '
             f'cross_engine_mix_gate={G} dex · synth_chi2_gate={TWO_ENGINE["synth_chi2_gate"]}._\n']
for title, lines in CASES.items():
    r = select_element(title.split(' (')[0], lines)
    lines_out.append(f'\n## {title}')
    lines_out.append(f'- **reported value = {r.value:.3f} ± {r.err:.3f}**  (n={r.n_lines}, '
                     f'engines used: {", ".join(r.selected_engines)})')
    lines_out.append(f'- diagnostic: Engine-A={r.engineA_value:.3f} · Engine-B={r.engineB_value:.3f} · '
                     f'mean cross-engine Δ(B−A)={r.mean_cross_engine_delta:+.3f} dex')
    lines_out.append(f'- cross-engine mix={r.cross_engine_mix} · **mix_flagged={r.mix_flagged}**'
                     + ('  ⚠️ ADJUDICATE (disagreement > gate; not a silent mean)' if r.mix_flagged else ''))
    for w in r.per_line:
        lines_out.append(f'    - {w.wavelength:.0f} {w.species}: **{w.engine}** ({w.value:.3f}) '
                         f'[{w.regime}] — {w.reason}; rejected {w.rejected_engine}='
                         f'{w.rejected_value if w.rejected_value is None else f"{w.rejected_value:.3f}"}')

out = 'data/audit/rya525_two_engine/sample_records.md'
open(out, 'w').write('\n'.join(lines_out) + '\n')
print('\n'.join(lines_out))
print(f'\n[wrote {out}]')
