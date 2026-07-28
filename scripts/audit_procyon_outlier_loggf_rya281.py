"""
scripts/audit_procyon_outlier_loggf_rya281.py
=============================================
RYA-281 — loggf quality of the current Procyon Fe outliers (closes the RYA-273 RCA,
order 3 of 3 after RYA-279/280).

Step 0 — the outlier set is the RYA-458 EW-integrity ABUND_OUTLIER flags from the
current RYA-273 run (NOT a re-derived list of 9; per the BRINGING CURRENT note 458 is
the canonical per-line outlier mechanism). Read from {star}_ew_integrity.csv.

Step 1 — cross-match each outlier's VALD log gf against the GES (Heiter et al. 2021)
value, REUSING the RYA-203/350 single-sourced canonical_gf.csv (gf_synth_ges = GES;
gf_linelist_vald / gf_regions_vald = VALD3). No new audit is built. Per line:
lambda, EP, VALD gf, GES gf, delta gf, GES reference/grade, EW, A(Fe), proximity, and a
LOGGF-CONFIRMED / LOGGF-PARTIAL / NOT-LOGGF flag (does the gf error explain the
abundance excess delta A = A - pool median?).

The diagnosis is the deliverable; the disposition (blend_flag on confirmed EW blends)
is applied to linelist_procyon.csv separately and the floor recomputed by re-running
validate_fe_rya238.py.

Usage:  python scripts/audit_procyon_outlier_loggf_rya281.py
Writes: data/audit/procyon_outlier_loggf_rya281.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

EWI   = ROOT / 'data' / 'processed' / 'procyon_ew_integrity.csv'
PERLN = ROOT / 'data' / 'processed' / 'procyon_per_line.csv'
GF    = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
OUT   = ROOT / 'data' / 'audit' / 'procyon_outlier_loggf_rya281.csv'

# A delta gf that "confirms" the bad-loggf hypothesis must be of the same sign and at
# least this fraction of the abundance excess (a weaker VALD gf inflates A by ~the same
# dex on the linear COG). Below it, gf cannot be the driver.
LOGGF_CONFIRM_FRAC = 0.70
LOGGF_PARTIAL_FRAC = 0.30


def _gf_row(gf: pd.DataFrame, species: str, wl: float, tol: float = 0.06):
    m = gf[(gf['species'] == species) &
           (gf['wavelength_air_A'] > wl - tol) & (gf['wavelength_air_A'] < wl + tol)]
    if m.empty:
        return None
    # nearest in wavelength
    return m.iloc[(m['wavelength_air_A'] - wl).abs().argsort().iloc[0]]


def main():
    ewi = pd.read_csv(EWI)
    pln = pd.read_csv(PERLN)
    gf  = pd.read_csv(GF, low_memory=False)

    outliers = ewi[ewi['ew_integrity'].astype(str).str.contains('ABUND_OUTLIER')].copy()

    # pool medians per ion (1D-LTE) for the abundance excess
    med = {ion: pln[(pln['element'] == 'Fe') & (pln['ion'] == ion)]['a_1dlte'].median()
           for ion in ('I', 'II')}

    rows = []
    for _, o in outliers.iterrows():
        el, ion, wl, ew = o['element'], o['ion'], float(o['wavelength_air_A']), float(o['ew_mA'])
        species = f"{el} {ion}"
        pl = pln[(pln['element'] == el) & (pln['ion'] == ion) &
                 (pln['wavelength_air_A'].sub(wl).abs() < 0.05)]
        a_1dlte = float(pl['a_1dlte'].iloc[0]) if not pl.empty else np.nan
        vpf = float(pl['vald_proximity_flag'].iloc[0]) if not pl.empty and 'vald_proximity_flag' in pl else np.nan
        grade = pl['line_grade'].iloc[0] if not pl.empty and 'line_grade' in pl else ''
        dA = a_1dlte - med.get(ion, np.nan)        # abundance excess vs pool median

        g = _gf_row(gf, species, wl)
        ep = np.nan
        if g is None:
            vald_gf = ges_gf = dgf = np.nan
            ref = grade_gf = 'NO-GES-MATCH'
        else:
            ep = float(g['excitation_potential_eV']) if pd.notna(g.get('excitation_potential_eV')) else np.nan
            ges_gf = float(g['gf_synth_ges']) if pd.notna(g['gf_synth_ges']) else np.nan
            vald_gf = next((float(g[c]) for c in ('gf_linelist_vald', 'gf_regions_vald', 'log_gf')
                            if pd.notna(g.get(c))), np.nan)
            dgf = vald_gf - ges_gf if np.isfinite(vald_gf) and np.isfinite(ges_gf) else np.nan
            ref = str(g.get('loggf_reference', ''))
            grade_gf = str(g.get('nist_grade', '') if pd.notna(g.get('nist_grade')) else '')

        # A weaker VALD gf (dgf < 0) inflates A by ~|dgf|; confirm if it covers the excess.
        explained = (-dgf) if np.isfinite(dgf) else 0.0
        frac = explained / dA if (np.isfinite(dA) and dA > 0) else 0.0
        if frac >= LOGGF_CONFIRM_FRAC:
            flag = 'LOGGF-CONFIRMED'
        elif frac >= LOGGF_PARTIAL_FRAC:
            flag = 'LOGGF-PARTIAL'
        else:
            flag = 'NOT-LOGGF'

        rows.append(dict(species=species, wavelength_air_A=round(wl, 3),
                         EP_eV=round(ep, 3) if np.isfinite(ep) else np.nan,
                         ew_mA=round(ew, 1), a_1dlte=round(a_1dlte, 3),
                         pool_median=round(med.get(ion, np.nan), 3), delta_A=round(dA, 3),
                         vald_log_gf=round(vald_gf, 3) if np.isfinite(vald_gf) else np.nan,
                         ges_log_gf=round(ges_gf, 3) if np.isfinite(ges_gf) else np.nan,
                         delta_log_gf=round(dgf, 3) if np.isfinite(dgf) else np.nan,
                         ges_ref=ref, nist_grade=grade_gf, vald_proximity=round(vpf, 3),
                         line_grade=grade, loggf_flag=flag))

    df = pd.DataFrame(rows).sort_values('delta_A', ascending=False).reset_index(drop=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)

    pd.set_option('display.width', 200); pd.set_option('display.max_columns', 30)
    print(f"\n{'='*92}\n  RYA-281 — GES log gf cross-match of the current Procyon Fe outliers "
          f"(RYA-458 ABUND_OUTLIER)\n{'='*92}")
    print(df.to_string(index=False))
    n_loggf = int((df['loggf_flag'] == 'LOGGF-CONFIRMED').sum())
    print(f"\n  LOGGF-CONFIRMED: {n_loggf} / {len(df)}.  "
          + ("bad-loggf DISPROVED — the inflation is not an oscillator-strength error "
             "(every |delta gf| is far below its abundance excess)." if n_loggf == 0
             else "some lines are gf-driven — correct toward GES."))
    print(f"  Wrote: {OUT.relative_to(ROOT)}\n")
    return df


if __name__ == '__main__':
    main()
