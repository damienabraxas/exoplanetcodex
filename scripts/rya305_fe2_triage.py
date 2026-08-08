"""
scripts/rya305_fe2_triage.py
RYA-305 — Fe II three-tier triage (DECISION 2026-06-14): clean(EW) / recover(synth) / drop.

Discriminant = synthetic Fe II-only EW at expected A(Fe) (SPECTRUM
calculate_theoretical_ew_and_depth) vs a measurability floor. Pre-flag anchored
to expected A(Fe)=solar+[Fe/H] (NOT this run's Fe I median). Reports the mixed
Fe II mean + Fe I-Fe II balance + the EW-vs-synthesis homogeneity guardrail.

Read-only analysis driving the pipeline fix. Reuses solar_per_line.csv (MOOG EW
a_1dlte) + solar_per_line_synth_v2.csv (synthesis a_synth).
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.simplefilter('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / 'ispec'))
import ispec
import pipeline.abundances_derive as ad
from pipeline.abundances_derive import (_load_solar_ews, _build_ispec_line_regions,
                                        _load_atmosphere, _load_synth_resources,
                                        _ISPEC_SOLAR_ABUND_FILE)
from config.constants import SOLAR_ASPLUND2021
from pipeline import data_namespace as ns
from pipeline import two_engine_inputs as tei

TEFF, LOGG, FEH, XI = 5772.0, 4.438, 0.0, 1.0
COG_THRESH = 0.5      # |log10(obs/theo)| below this = clean (obs consistent with theory)
EW_FLOOR_MA = 5.0     # synthetic Fe II-only EW floor (mA): below = no real line -> drop

# RYA-682: these two per-line tables are RYA-469-namespaced products of the solar
# run — data/outputs/solar/, not the pre-namespacing data/processed/ this script
# used to read. That path drifted at RYA-469 and was never migrated: nothing writes
# a per-line table to data/processed/ any more, so a clean checkout hit a bare
# pandas FileNotFoundError and stale worktrees quietly served June copies.
PER_LINE = ns.output_path('solar', 'per_line.csv', create=False)
PER_LINE_SYNTH_V2 = tei.engine_b_per_line_path('solar')


def main():
    ew = _load_solar_ews()
    lm = _build_ispec_line_regions(ew)
    notes = np.array([str(n) for n in lm['note']])
    fe2 = lm[notes == 'Fe 2']
    atm = _load_atmosphere(TEFF, LOGG, FEH, XI)
    _, isotopes, _ = _load_synth_resources()
    solar_ab = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    out = ispec.calculate_theoretical_ew_and_depth(atm, TEFF, LOGG, FEH, 0.0,
            fe2.copy(), isotopes, solar_ab, microturbulence_vel=XI, verbose=0)
    theo = np.asarray(out['theoretical_ew'], dtype=float)
    obs = np.asarray(fe2['ew'], dtype=float)
    wave = np.asarray(fe2['wave_A'], dtype=float)

    moog = pd.read_csv(PER_LINE)
    syn = pd.read_csv(PER_LINE_SYNTH_V2)
    moog = moog[(moog.element == 'Fe') & (moog.ion.astype(str) == 'II')]
    syn = syn[(syn.element == 'Fe') & (syn.ion.astype(str) == 'II')]
    fe1 = pd.read_csv(PER_LINE)
    fe1 = fe1[(fe1.element == 'Fe') & (fe1.ion.astype(str) == 'I')]
    a_fe1 = float(np.nanmedian(fe1['a_1dlte']))

    def near(df, w, col):
        i = (df['wavelength_air_A'] - w).abs().idxmin()
        return float(df.loc[i, col])

    rows = []
    for k in range(len(fe2)):
        w = wave[k]
        logr = np.log10(max(obs[k], 1e-6) / max(theo[k], 1e-6))
        a_moog = near(moog, w, 'a_1dlte')
        a_syn = near(syn, w, 'a_synth')
        if abs(logr) < COG_THRESH:
            tier, a_use = 'clean', a_moog
        elif theo[k] >= EW_FLOOR_MA:
            tier, a_use = 'recover', a_syn
        else:
            tier, a_use = 'drop', np.nan
        rows.append(dict(wave=round(w, 3), obs_ew=round(obs[k], 1), theo_ew=round(theo[k], 3),
                         logratio=round(logr, 2), a_moog=round(a_moog, 3),
                         a_synth=round(a_syn, 3), tier=tier, a_use=round(a_use, 3) if np.isfinite(a_use) else np.nan))
    df = pd.DataFrame(rows).sort_values('a_moog')
    print(df.to_string(index=False))

    clean = df[df.tier == 'clean']; rec = df[df.tier == 'recover']; drop = df[df.tier == 'drop']
    fe2_final = df[df.tier != 'drop']['a_use'].dropna()
    print(f"\nTiers: clean={len(clean)} recover={len(rec)} drop={len(drop)}")
    print(f"Dropped (quarantine): {sorted(drop.wave.tolist())}")
    print(f"Recovered (synth): {sorted(rec.wave.tolist())}")
    # guardrail: EW vs synth on the clean overlap
    do = (clean.a_moog - clean.a_synth)
    print(f"\nGUARDRAIL EW-vs-synth on clean Fe II: mean offset={do.mean():+.3f}  scatter={do.std():.3f}  (n={len(clean)})")
    print(f"  Fe I: EW median={a_fe1:.3f}  synth median={float(np.nanmedian(syn['a_synth'])):.3f} (Fe II synth)")
    print(f"\nFe II final (mixed) median={fe2_final.median():.3f} mean={fe2_final.mean():.3f} n={len(fe2_final)}")
    print(f"Fe I (clean EW) median={a_fe1:.3f}")
    print(f"==> Fe I - Fe II (EW Fe I vs mixed Fe II) = {a_fe1 - fe2_final.median():+.3f} dex  (target <0.05)")
    # all-synthesis reference
    v2 = pd.read_csv(ns.output_path('solar','abundances_synth_v2.csv', create=False))
    fe2_syn = df['a_synth'].dropna()
    fe1_syn = float(v2[(v2.element == 'Fe') & (v2.ion.astype(str) == 'I')]['A_X_abs'].iloc[0])
    print(f"All-synthesis ref: Fe I {fe1_syn:.3f} - Fe II {fe2_syn.median():.3f} = {fe1_syn - fe2_syn.median():+.3f}")


if __name__ == '__main__':
    main()
