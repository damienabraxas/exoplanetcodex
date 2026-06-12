"""
scripts/validate_fe_rya238.py
==============================
RYA-238: Fe I/II focused validation — solar + Procyon.

Runs abundances_derive.run() for both stars with vmic fixed at literature
values, then produces the four required diagnostic outputs:

  1. Per-line table (wavelength, EW, a_1dlte, aberr, a_nlte, Δ from median)
  2. Excitation plot: A(Fe I) vs EP
  3. REW plot:        A(Fe I) vs reduced EW (log10(EW/λ))
  4. Ionisation balance scatter: A(Fe I) vs A(Fe II) per star

Gate checks (solar and Procyon) are printed verbatim for the Linear comment.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from pipeline.abundances_derive import run

# ── Gate definitions (from issue spec) ───────────────────────────────────────
GATES = {
    'solar': {
        'A_Fe_I_lo'  : 7.41, 'A_Fe_I_hi'  : 7.51,
        'A_Fe_II_lo' : 7.41, 'A_Fe_II_hi' : 7.51,
        'dFe_max'    : 0.05,
        'scatter_max': 0.10,
        'n_Fe2_min'  : 8,
        'vmic'       : 1.00,
    },
    'procyon': {
        'A_Fe_I_lo'  : 7.38, 'A_Fe_I_hi'  : 7.54,
        'A_Fe_II_lo' : 7.38, 'A_Fe_II_hi' : 7.54,
        'dFe_max'    : 0.08,
        'scatter_max': 0.15,
        'n_Fe2_min'  : 3,
        'vmic'       : 1.66,
    },
}

PLOT_DIR = REPO_ROOT / 'data' / 'processed'


def _check_gates(star_id: str, abundances: pd.DataFrame, per_line: pd.DataFrame,
                 vmic_val: float) -> dict:
    g = GATES[star_id]
    fe = abundances[abundances['element'] == 'Fe']

    def _get(ion):
        row = fe[fe['ion'] == ion]
        if row.empty:
            return {}
        r = row.iloc[0]
        return {
            'a_nlte'   : float(r.get('A_X_nlte', r['A_X'])),
            'a_1dlte'  : float(r['A_X']),
            'scatter'  : float(r.get('A_X_std_nlte', r.get('A_X_std', np.nan))),
            'n_lines'  : int(r['n_lines']),
            'delta_nlte': float(r.get('delta_nlte_mean', np.nan)),
        }

    fe1 = _get('I')
    fe2 = _get('II')

    a1  = fe1.get('a_nlte', np.nan)
    a2  = fe2.get('a_nlte', np.nan)
    dfe = abs(a1 - a2) if (np.isfinite(a1) and np.isfinite(a2)) else np.nan
    sc  = fe1.get('scatter', np.nan)
    n2  = fe2.get('n_lines', 0)

    results = {
        'A_Fe_I_nlte' : a1,
        'A_Fe_I_pass' : g['A_Fe_I_lo'] <= a1 <= g['A_Fe_I_hi'] if np.isfinite(a1) else False,
        'A_Fe_II_nlte': a2,
        'A_Fe_II_pass': g['A_Fe_II_lo'] <= a2 <= g['A_Fe_II_hi'] if np.isfinite(a2) else False,
        'dFe'         : dfe,
        'dFe_pass'    : dfe <= g['dFe_max'] if np.isfinite(dfe) else False,
        'scatter'     : sc,
        'scatter_pass': sc < g['scatter_max'] if np.isfinite(sc) else False,
        'n_Fe2'       : n2,
        'n_Fe2_pass'  : n2 >= g['n_Fe2_min'],
        'vmic'        : vmic_val,
        'vmic_lit'    : g['vmic'],
        'fe1'         : fe1,
        'fe2'         : fe2,
    }
    return results


def _print_gate_table(star_id: str, res: dict):
    g = GATES[star_id]
    print(f"\n{'='*62}")
    print(f"  GATE TABLE — {star_id.upper()}")
    print(f"{'='*62}")
    print(f"  {'Metric':<28} {'Value':>10}  {'Gate':>16}  {'Status'}")
    print(f"  {'-'*60}")

    def row(name, val, gate_str, passed):
        val_s  = f"{val:.4f}" if np.isfinite(val) else "  N/A  "
        status = "PASS" if passed else "FAIL"
        print(f"  {name:<28} {val_s:>10}  {gate_str:>16}  {status}")

    row("A(Fe I) NLTE",
        res['A_Fe_I_nlte'],
        f"[{g['A_Fe_I_lo']}, {g['A_Fe_I_hi']}]",
        res['A_Fe_I_pass'])
    row("A(Fe II) NLTE",
        res['A_Fe_II_nlte'],
        f"[{g['A_Fe_II_lo']}, {g['A_Fe_II_hi']}]",
        res['A_Fe_II_pass'])
    row("ΔFe(I−II)",
        res['dFe'],
        f"< {g['dFe_max']}",
        res['dFe_pass'])
    row("Fe I scatter (σ NLTE)",
        res['scatter'],
        f"< {g['scatter_max']}",
        res['scatter_pass'])
    row(f"Fe II n_lines",
        float(res['n_Fe2']),
        f">= {g['n_Fe2_min']}",
        res['n_Fe2_pass'])
    row("vmic (fixed)",
        res['vmic'],
        f"{g['vmic']} km/s (lit)",
        abs(res['vmic'] - g['vmic']) < 0.01)

    all_pass = all([res['A_Fe_I_pass'], res['A_Fe_II_pass'], res['dFe_pass'],
                    res['scatter_pass'], res['n_Fe2_pass']])
    print(f"\n  Overall: {'ALL PASS ✓' if all_pass else 'GATES FAIL ✗'}")
    return all_pass


def _per_line_table(star_id: str, per_line: pd.DataFrame) -> pd.DataFrame:
    fe = per_line[per_line['element'] == 'Fe'].copy()
    if fe.empty:
        print(f"  No Fe lines in per_line output for {star_id}")
        return fe

    for ion in ('I', 'II'):
        sub = fe[fe['ion'] == ion]
        if sub.empty:
            continue
        if 'a_nlte' in sub.columns and sub['a_nlte'].notna().any():
            median_nlte = float(np.nanmedian(sub['a_nlte']))
            fe.loc[fe['ion'] == ion, 'delta_from_median'] = (
                fe.loc[fe['ion'] == ion, 'a_nlte'] - median_nlte
            ).round(4)
        else:
            median_1d = float(np.nanmedian(sub['a_1dlte']))
            fe.loc[fe['ion'] == ion, 'delta_from_median'] = (
                fe.loc[fe['ion'] == ion, 'a_1dlte'] - median_1d
            ).round(4)

    cols = ['ion', 'wavelength_air_A', 'ew_mA', 'a_1dlte']
    if 'aberr' in fe.columns:
        cols += ['aberr']
    if 'a_nlte' in fe.columns:
        cols += ['a_nlte']
    cols += ['delta_from_median']

    for ion in ('I', 'II'):
        sub = fe[fe['ion'] == ion].sort_values('wavelength_air_A')[cols]
        print(f"\n  Fe {ion} per-line ({len(sub)} lines) — {star_id}")
        print(sub.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    return fe


def _plot_ep(star_id: str, per_line: pd.DataFrame, ax):
    fe1 = per_line[(per_line['element'] == 'Fe') & (per_line['ion'] == 'I')].copy()
    if fe1.empty:
        return
    col = 'a_nlte' if ('a_nlte' in fe1.columns and fe1['a_nlte'].notna().any()) else 'a_1dlte'
    ep = fe1['excitation_potential_eV'].values
    ab = fe1[col].values
    mask = np.isfinite(ep) & np.isfinite(ab)
    ax.scatter(ep[mask], ab[mask], s=18, alpha=0.7, color='steelblue')
    if mask.sum() >= 2:
        m, b = np.polyfit(ep[mask], ab[mask], 1)
        xs = np.linspace(ep[mask].min(), ep[mask].max(), 50)
        ax.plot(xs, m * xs + b, 'r--', lw=1, label=f"slope={m:+.4f} dex/eV")
    ax.axhline(np.nanmedian(ab[mask]), color='gray', lw=0.8, ls=':')
    ax.set_xlabel("Excitation potential (eV)")
    ax.set_ylabel("A(Fe I)" + (" NLTE" if col == 'a_nlte' else " 1D LTE"))
    ax.set_title(f"{star_id} — A(Fe I) vs EP")
    ax.legend(fontsize=8)


def _plot_rew(star_id: str, per_line: pd.DataFrame, ax):
    fe1 = per_line[(per_line['element'] == 'Fe') & (per_line['ion'] == 'I')].copy()
    if fe1.empty:
        return
    col = 'a_nlte' if ('a_nlte' in fe1.columns and fe1['a_nlte'].notna().any()) else 'a_1dlte'
    ew  = fe1['ew_mA'].values
    wav = fe1['wavelength_air_A'].values
    ab  = fe1[col].values
    rew = np.log10(ew / wav)
    mask = np.isfinite(rew) & np.isfinite(ab) & (ew > 0)
    ax.scatter(rew[mask], ab[mask], s=18, alpha=0.7, color='darkorange')
    if mask.sum() >= 2:
        m, b = np.polyfit(rew[mask], ab[mask], 1)
        xs = np.linspace(rew[mask].min(), rew[mask].max(), 50)
        ax.plot(xs, m * xs + b, 'r--', lw=1, label=f"slope={m:+.3f}")
    ax.axhline(np.nanmedian(ab[mask]), color='gray', lw=0.8, ls=':')
    ax.set_xlabel("log₁₀(EW/λ)  [reduced EW]")
    ax.set_ylabel("A(Fe I)" + (" NLTE" if col == 'a_nlte' else " 1D LTE"))
    ax.set_title(f"{star_id} — A(Fe I) vs REW")
    ax.legend(fontsize=8)


def _plot_ion_balance(results_solar, results_procyon, ax):
    for star_id, ab_df, color in [
        ('Solar',   results_solar,   'steelblue'),
        ('Procyon', results_procyon, 'darkorange'),
    ]:
        fe = ab_df[ab_df['element'] == 'Fe']
        fe1 = fe[fe['ion'] == 'I']
        fe2 = fe[fe['ion'] == 'II']
        if fe1.empty or fe2.empty:
            continue
        a1 = float(fe1.iloc[0].get('A_X_nlte', fe1.iloc[0]['A_X']))
        a2 = float(fe2.iloc[0].get('A_X_nlte', fe2.iloc[0]['A_X']))
        ax.scatter([a1], [a2], s=80, color=color, label=f"{star_id}: ΔFe={a1-a2:+.3f}")
    ax.plot([7.2, 7.8], [7.2, 7.8], 'k--', lw=0.8, label="1:1")
    ax.set_xlabel("A(Fe I) NLTE")
    ax.set_ylabel("A(Fe II) NLTE")
    ax.set_title("Ionisation balance: Fe I vs Fe II")
    ax.legend(fontsize=8)


def main():
    print("\n" + "="*62)
    print("  RYA-238: Fe I/II validation — solar + Procyon")
    print("="*62 + "\n")

    # ── Solar run ────────────────────────────────────────────────────
    print("\n>>> Running: solar (vmic fixed = 1.00 km/s)")
    conv_solar, ab_solar = run('solar', vmic_fixed=True)

    solar_pl_path = REPO_ROOT / 'data' / 'processed' / 'solar_per_line.csv'
    solar_pl = pd.read_csv(solar_pl_path) if solar_pl_path.exists() else pd.DataFrame()

    # ── Procyon run ──────────────────────────────────────────────────
    print("\n>>> Running: procyon (vmic fixed = 1.66 km/s)")
    conv_proc, ab_proc = run('procyon', vmic_fixed=True)

    proc_pl_path = REPO_ROOT / 'data' / 'processed' / 'procyon_per_line.csv'
    proc_pl = pd.read_csv(proc_pl_path) if proc_pl_path.exists() else pd.DataFrame()

    # ── Per-line tables ──────────────────────────────────────────────
    print("\n" + "="*62)
    print("  PER-LINE TABLES")
    print("="*62)
    _per_line_table('solar',   solar_pl)
    _per_line_table('procyon', proc_pl)

    # ── NLTE outliers (|a_nlte - median| > 0.2 dex) ─────────────────
    print("\n" + "="*62)
    print("  NLTE OUTLIERS (|a_nlte − median| > 0.2 dex)")
    print("="*62)
    for star_id, pl in [('solar', solar_pl), ('procyon', proc_pl)]:
        fe = pl[pl['element'] == 'Fe'] if not pl.empty else pl
        if 'a_nlte' not in fe.columns or fe.empty:
            print(f"  {star_id}: no NLTE column in per-line output")
            continue
        for ion in ('I', 'II'):
            sub = fe[fe['ion'] == ion]
            if sub.empty:
                continue
            med = float(np.nanmedian(sub['a_nlte']))
            outliers = sub[np.abs(sub['a_nlte'] - med) > 0.2]
            if outliers.empty:
                print(f"  {star_id} Fe {ion}: no outliers")
            else:
                for _, r in outliers.iterrows():
                    ab_str = f"{r['aberr']:+.4f}" if 'aberr' in r and np.isfinite(r['aberr']) else "  n/a"
                    print(f"  {star_id} Fe {ion} {r['wavelength_air_A']:.3f} Å  "
                          f"a_1dlte={r['a_1dlte']:.4f}  aberr={ab_str}  "
                          f"a_nlte={r['a_nlte']:.4f}  Δ={r['a_nlte']-med:+.4f}")

    # ── Gate tables ──────────────────────────────────────────────────
    res_solar  = _check_gates('solar',   ab_solar,  solar_pl,  conv_solar.get('vturb_kms', 1.00))
    res_procyon = _check_gates('procyon', ab_proc,   proc_pl,   conv_proc.get('vturb_kms', 1.66))

    pass_solar   = _print_gate_table('solar',   res_solar)
    pass_procyon = _print_gate_table('procyon', res_procyon)

    # ── Diagnostic plots ─────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("RYA-238: Fe I/II Validation — Solar & Procyon", fontsize=13)

    _plot_ep('solar',   solar_pl,  axes[0, 0])
    _plot_ep('procyon', proc_pl,   axes[1, 0])
    _plot_rew('solar',   solar_pl, axes[0, 1])
    _plot_rew('procyon', proc_pl,  axes[1, 1])
    _plot_ion_balance(ab_solar, ab_proc, axes[0, 2])
    axes[1, 2].axis('off')  # spare panel

    plt.tight_layout()
    plot_path = PLOT_DIR / 'rya238_fe_validation.png'
    fig.savefig(plot_path, dpi=150)
    plt.close(fig)
    print(f"\n  Plots saved → {plot_path}")

    # ── Summary ──────────────────────────────────────────────────────
    print("\n" + "="*62)
    print("  SUMMARY")
    print("="*62)
    print(f"  Solar:   {'ALL PASS' if pass_solar   else 'GATES FAIL'}")
    print(f"  Procyon: {'ALL PASS' if pass_procyon else 'GATES FAIL'}")
    if pass_solar and pass_procyon:
        print("\n  → Proceed to RYA-239 (full 27-element solar run).")
    else:
        print("\n  → Fix failing gates before full run.")
    print()


if __name__ == '__main__':
    main()
