"""RYA-1208 — publish every Gerber cell this ticket produced.

One invocation per tier so the reason travels with the change, then the RYA-1178 schema
emitter, which is the required second half of any publish (it re-derives grade, the xi
terms, sigma_syst_complete and sigma_reported; publishing alone strips 23 fields).
"""
import subprocess, sys, glob, os, json

ROOT = os.path.dirname(os.path.abspath(__file__))
BP = "data/results/band_products"
REASON = ("RYA-1208: completing the Gerber matrix. `synth-1D-LTE-gerber` is Turbospectrum "
          "LTE on the Gerber deck's OWN atmosphere (MARCS.GES) with departures withheld, "
          "so it needs no NLTE labels and is runnable in every band the synthesis route "
          "reaches -- it was simply never run outside VIS and red-optical. Nothing is "
          "calibrated: the value is whatever the fit returns (RYA-161).")

def publish(srcs, tier, extra=()):
    if not srcs:
        return None
    cmd = [sys.executable, "scripts/publish_product.py", "--from", *srcs,
           "--tier", tier, "--selector", tier, "--route", "SYNTH", "--star", "solar",
           "--host", "Sirius", "--reason", REASON, *extra]
    print("  $", " ".join(x if " " not in x else "…" for x in cmd[:12]), "…")
    r = subprocess.run(cmd, capture_output=True, text=True)
    print(r.stdout.strip()[-600:] or r.stderr.strip()[-600:])
    return r.returncode

if __name__ == "__main__":
    # holdings differ per file, so publish one holding at a time (--holding is file-wide)
    groups = {}
    for f in sorted(glob.glob(f"{BP}/*_synth-1D-LTE-gerber_products.csv")):
        base = os.path.basename(f)
        for h in ("solar_kpno_kurucz2005_corrected", "solar_kpno_molecfit_corrected",
                  "solar_iag", "solar_crires_plus_y_wide_rya1054",
                  "solar_crires_plus_h_rya1094", "solar_harps_molecfit_corrected"):
            if h in base:
                tier = "DEEPGRADED" if "_DEEPGRADED_" in base else "GRADED"
                groups.setdefault((h, tier), []).append(f)
                break
    for (h, tier), srcs in sorted(groups.items()):
        print(f"\n=== {h}  tier={tier}  ({len(srcs)} file(s))")
        publish(srcs, tier, extra=["--holding", h,
                                   "--origin-path", f"/home/damienabraxas/rya1208/{BP}/"])
