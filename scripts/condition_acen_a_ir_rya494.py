# =============================================================================
# THE EXOPLANET CODEX  —  exoplanetcodex.org
# =============================================================================
# File:         condition_acen_a_ir_rya494.py
# Module:       scripts (caller #1 of the generalized IR telluric module)
# Description:  Applies pipeline.ir_telluric (RYA-494) to the α Cen A IR — the
#               first non-Vesta caller of the generalized per-star/per-instrument
#               telluric module. CRIRES+ → molecfit (stellar-BERV branch); NIRPS →
#               select FLUX_TELL_* (already telluric-corrected + BARYCENT). Verifies
#               the telluric-specific residual, RV-conditions to α Cen A's rest, and
#               writes telluric-clean IR + a conditioning summary.
#
#               HONEST SCOPE: the brief's headline atomic lines O I 844.6/926.6 nm
#               fall BELOW both the CRIRES-Y (949.6 nm) and NIRPS (966.1 nm) blue
#               edges → NOT covered by α Cen A IR. The ¹³C/CO molecular isotope work
#               stays STAGGER-walled. α Cen B has NO NIRPS (RYA-439).
#
# Author:       Ryan Schmitt  |  Contributors: Claude (Anthropic) via Claude Code
# Linear issue: RYA-494
# =============================================================================
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from pipeline import ir_telluric as irt
from pipeline.crires_telluric import NM_TO_A

OUT_DIR = Path("data/audit/acen_a_ir_rya494")
# atomic IR lines of interest (vacuum Å). The brief's O I 844/926 + a few in-range refs.
ATOMIC_IR_LINES = {"O I 8446": 8446.5, "O I 9266": 9265.9,
                   "Fe I 11593": 11593.6, "Mg I 11831": 11831.5,  # NIRPS-covered atomic IR
                   "Si I 10827": 10827.1}


def _coverage(line_A, crires, nirps):
    cr = [f.wlen_id for f in crires if f.covers_nm(line_A / NM_TO_A)]
    nr = any(f.covers_A(line_A) for f in nirps)
    return {"crires_settings": cr, "nirps": nr, "covered": bool(cr) or nr}


def main(run_molecfit: bool = True):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t = irt.TARGETS["alpha_cen_a"]
    crires = irt.crires_frames(t)
    nirps = irt.nirps_frames(t)

    summary = {
        "ticket": "RYA-494", "target": "alpha_cen_a",
        "velocity_mode": t.velocity_mode.value,
        "crires_frames": {f.wlen_id: {"nm": [round(f.wmin_nm, 1), round(f.wmax_nm, 1)],
                                      "specsys": f.specsys, "snr": round(f.snr, 1)}
                          for f in crires},
        "nirps_n_frames": len(nirps),
        "nirps_attribution": t.nirps_attribution,
        "nirps_telluric_column": nirps[0].telluric_column if nirps else None,
        "nirps_specsys": nirps[0].specsys if nirps else None,
        "atomic_ir_coverage": {name: _coverage(wl, crires, nirps)
                               for name, wl in ATOMIC_IR_LINES.items()},
        "flags": {
            "O_I_844_926": "NOT COVERED — both fall below CRIRES-Y (9496Å) and NIRPS "
                           "(9661Å) blue edges; the brief's headline atomic lines are "
                           "not reachable in α Cen A IR (coverage gap, not a telluric gap).",
            "13C_CO": "STAGGER-pending — molecular CO ¹³C/¹²C needs the Amarsi 3D model "
                      "(RYA-373 wall); deferred, not silently dropped.",
            "alpha_cen_B_IR": "NO NIRPS (RYA-439: the 'ALPHACENB' program is α Cen A); "
                              "B has only 1 CRIRES frame. B IR is a data-acquisition gap.",
        },
    }

    # NIRPS path (fast): select telluric column + stellar RV-condition a representative
    if nirps:
        n0 = nirps[0]
        irt.telluric_correct_nirps(n0)
        irt.velocity_condition(n0, t)
        summary["nirps_demo"] = {
            "file": n0.path.name, "telluric_corrected": n0.telluric_corrected,
            "rv_applied": n0._rv_applied, "rest_frame": n0.rest_frame,
            "atm_transm_min": round(n0._atm_transm_min, 3)}

    # CRIRES path (molecfit): telluric-correct the best frame + verify the gate + RV
    if run_molecfit and crires:
        y = max(crires, key=lambda f: f.snr)          # Y1029, SNR 302
        win = (1100.0, 1120.0)                          # telluric-rich Y-band H2O
        irt.telluric_correct_crires(y, win, Path(f"/tmp/rya494_acenA_{y.wlen_id}"))
        gate = irt.telluric_residual_gate(y)
        irt.velocity_condition(y, t)
        summary["crires_demo"] = {
            "frame": y.wlen_id, "window_nm": win, "molecfit_gdas": y._gdas,
            "telluric_corrected": y.telluric_corrected,
            "telluric_gate": {k: (round(v, 4) if isinstance(v, float) else v)
                              for k, v in gate.items()},
            "rv_applied": y._rv_applied, "rest_frame": y.rest_frame}
        # write the telluric-clean, RV-conditioned IR product (the deliverable)
        _write_clean_product(y, gate)
        summary["crires_product"] = str(OUT_DIR / f"acen_a_CRIRES_{y.wlen_id}_telluric_clean_rya494.fits")

    (OUT_DIR / "acen_a_ir_conditioning_rya494.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    return summary


def _write_clean_product(frame, gate):
    """Persist the telluric-clean, stellar-rest CRIRES IR (the conditioned deliverable)."""
    ph = fits.PrimaryHDU(); h = ph.header
    h["RYA"] = "RYA-494"
    h["TARGET"] = "alpha_cen_a"
    h["TELL"] = ("molecfit", "ESO molecfit (esorex), topocentric, H2O+CH4+CO2")
    h["TELLEXC"] = (round(gate.get("excess", float("nan")), 4), "telluric-specific excess residual")
    h["TELLPASS"] = (bool(gate.get("passed", False)), "telluric gate <2% excess")
    h["RESTFRM"] = (bool(frame.rest_frame), "shifted to alpha Cen A rest (BERV+systemic)")
    h["VELMODE"] = ("stellar", "direct-star BERV branch (NOT asteroid ephemeris)")
    h["WLEN"] = frame.wlen_id
    h["GDAS"] = getattr(frame, "_gdas", "?")
    hdus = [ph]
    for s in frame.segments:
        if not np.isfinite(getattr(s, "_mtrans", np.array([np.nan]))).any():
            continue
        hdus.append(fits.BinTableHDU.from_columns([
            fits.Column(name="wave_rest_A", format="1D", array=s.wave_A),
            fits.Column(name="flux_tellclean", format="1D", array=s.flux),
            fits.Column(name="mtrans", format="1D", array=s._mtrans)],
            name=f"ORD{s.order}_DET{s.detector}"))
    out = OUT_DIR / f"acen_a_CRIRES_{frame.wlen_id}_telluric_clean_rya494.fits"
    fits.HDUList(hdus).writeto(out, overwrite=True)
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-molecfit", action="store_true", help="skip the slow molecfit pass")
    a = ap.parse_args()
    main(run_molecfit=not a.no_molecfit)
