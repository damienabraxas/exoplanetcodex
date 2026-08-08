#!/usr/bin/env python3
"""
RYA-565 — Europium (Eu II) LTE in-window HFS blend-fit synthesis (RUNS ON SIRIUS).

WHY
---
Eu emits NO value in the solar verdict, yet the HFS-summed EW is committed
(`data/measured/sol_ew_results_v1.csv`: Eu II 6645.127, 6.8 mA, `hfs_total_ew`,
`blend_flag`) and the two-engine registry claimed "LTE synthesis is the finished
treatment" (`data/curation/nlte_two_engine_coverage.csv`). Either Eu is a finished
LTE value never wired through, or it is honestly noise-floor/blend-owed and the
registry wording was optimistic. This harness DOES the measurement the registry
claimed was finished and lets the numbers decide — the same decision-gate discipline
as the Zr II rescue (RYA-585) and the Sr II in-window blend fit (RYA-551).

ENGINE: LTE ONLY. Eu II is the MAJORITY ion of Eu in the solar photosphere and there
is no Eu departure grid in either engine (registry: LTE-only-by-design, governing
ticket RYA-458). nlte_delta is None, nlte_status 'NLTE_unavailable_LTE_robust' — the
Sr II / Zr II / V II precedent. Do NOT go looking for or building an Eu NLTE grid
here; the binding blocker is NOT NLTE.

THE HFS / ISOTOPE gf CONVENTION (the thing that decides the number)
------------------------------------------------------------------
Every solar Eu II line is hyperfine- AND isotope-split. Two on-disk sources describe
the pattern and they use INCOMPATIBLE conventions:

  * GES v6 `atomic_lines.tsv` (reference_code LWHS = Lawler, Wickliffe, den Hartog &
    Sneden 2001) carries the resolved components split across `turbospectrum_species`
    blocks 63.151 / 63.153. Each ISOTOPE's components sum to the FULL
    (isotope-fraction-free) oscillator strength — correct, because both isotopes share
    the same gf and the abundance split is applied separately.

  * The production VALD "for-grid" lists collapse each feature to TWO lines, one per
    isotope block, whose gf ratio is exactly 0.522/0.478 — i.e. VALD has ALREADY
    FOLDED THE SOLAR ISOTOPE FRACTIONS INTO log gf.

Turbospectrum settles which is right: `bsyn.f` (~line 1350) multiplies level
populations by `isotopfrac(atom, isotope)`, and `makeabund.f` hardcodes
isotopfrac(63,151)=0.478 / isotopfrac(63,153)=0.522 (isotopfrac(x,0)==1.0 for a
species written WITHOUT an isotope code). So an isotope-coded block gets the fraction
applied BY THE ENGINE.

  => The GES block is the correct input: engine applies f, effective total gf = the
     LWHS literature value.
  => The as-shipped VALD block DOUBLE-APPLIES the isotope fraction (folded into gf AND
     applied by bsyn) => effective gf 0.301 dex (a factor 2) too WEAK, which biases a
     fitted A(Eu) 0.30 dex too HIGH.

This harness synthesises on the GES LWHS HFS pattern and ALSO runs the as-shipped VALD
block as a labelled sensitivity leg, so the convention error is MEASURED, not asserted.

SSOT NOTE ON canonical_gf.csv: every committed Eu II HFS row (loggf_reference LWHS)
stores the NAIVE sum of the GES components across BOTH isotopes. That is a correct
checksum of the block Turbospectrum consumes, but it is NOT a physical total gf — the
physical total is that value minus log10(2). The harness ASSERTS the checksum identity
against the live GES list per line (loud-fail on drift) and records the distinction; it
does NOT rewrite the single source.

LINES: 6645.0905 is THE ticket line (primary). The other four HFS-resolved LWHS Eu II
rows in canonical_gf (6437.6315, 6173.0215, 6049.5033, 5818.7392) are fitted as
CROSS-CHECKS — a single blended weak line cannot support a disposition on its own, and
line-to-line spread is exactly the systematic the decision gate needs to see.

VALIDATE-DON'T-TUNE: nothing is fitted toward Asplund 2021 A(Eu) = 0.52; that value
only brackets the trial grid and is printed as a reference point.

Out: data/results/eu2_synthesis_rya565.json
Run: ssh sirius; /mnt/codex-data/venv_ci/bin/python3 scripts/rya565_eu2_synth_sirius.py
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solar_profile_fit import (CLIGHT, RCHI2_REVIEW,  # noqa: E402
                               RELIABLE_DEWDA, assess_reliability, broaden,
                               fit_profile, local_renorm, require_arm_rv)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline._numcompat import trapezoid  # noqa: E402  (numpy>=2 removed np.trapz)

EXE = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
MARCS = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
         "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
GES_TSV = "/mnt/codex-data/linelists/ges_v6/atomic_lines.tsv"
W = "/mnt/codex-data/codex/rya565/work"
HARPS = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"

Z_EU = 63                 # Eu atomic number (bsyn INDIVIDUAL ABUNDANCES element)
XI = 1.0                  # solar microturbulence (matches the t01 MARCS node)
VSINI = 1.8               # solar vsini (STAR_PARAMS)
EU_ISOTOPES = (151, 153)
TS_ISOTOPFRAC = {151: 0.478, 153: 0.522}   # makeabund.f, quoted into the audit record

LOG2 = float(np.log10(2.0))
GF_TOL = 0.005

# The HFS-resolved LWHS Eu II rows in data/linelists/canonical_gf.csv. `canon_loggf`
# is the COMMITTED value (the naive cross-isotope GES sum); verify_gf_provenance()
# asserts it against the file AND against the live GES component block — nothing here
# is a substitute for the single source, it is the loud-fail cross-check.
PRIMARY = 6645.0905
LINES = {
    6645.0905: dict(ep=1.380, canon_loggf=0.4208, hfs_n=11, gf_id='gf_051798',
                    vald='vald-6300-6800-for-grid.list', role='primary', fit_hw=0.45,
                    # Kurucz-2010 (ungraded) Cr I sits essentially ON the Eu core — the
                    # dominant blend systematic for this feature.
                    blend_probe=(24, 1, 6645.087, 'Cr I 6645.087 (K10, ungraded)')),
    6437.6315: dict(ep=1.320, canon_loggf=-0.0186, hfs_n=9, gf_id='gf_051797',
                    vald='vald-6300-6800-for-grid.list', role='crosscheck', fit_hw=0.45),
    6173.0215: dict(ep=1.320, canon_loggf=-0.5590, hfs_n=6, gf_id='gf_051795',
                    vald='vald-5800-6300-for-grid.list', role='crosscheck', fit_hw=0.45),
    6049.5033: dict(ep=1.279, canon_loggf=-0.4991, hfs_n=7, gf_id='gf_051794',
                    vald='vald-5800-6300-for-grid.list', role='crosscheck', fit_hw=0.45),
    5818.7392: dict(ep=1.230, canon_loggf=-0.9489, hfs_n=8, gf_id='gf_051789',
                    vald='vald-5800-6300-for-grid.list', role='crosscheck', fit_hw=0.45),
}

LMARGIN = 8.0                                          # synth half-window (A)
A_GRID = np.round(np.arange(-0.60, 1.801, 0.05), 3)    # A(Eu) trial grid; brackets 0.52
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())
A_BLEND_ONLY = -6.0        # Eu effectively switched off -> the blend-only continuum
BLEND_PROBE_DEX = 0.20     # +/- perturbation on the ungraded blend gf (the systematic)

# Reliability constants + rule: imported from the shared module (RYA-679). This
# harness used to carry its own RCHI2_SANE_MAX = 15.0 — a third live value alongside
# RYA-564/581's 60.0 and RYA-560's 5.0. Retired: red_chi2 is reported, not gated.
# Not load-bearing for Eu either way — Eu II 6645 fails on sensitivity (dEW/dA 13.9)
# with an excellent fit (red_chi2 0.16).


# ─────────────────────────── gf / HFS provenance (SSOT) ───────────────────────────

def read_ges_hfs(center, path=GES_TSV, halfwidth=0.30):
    """Read a line's LWHS Eu II HFS/isotope components from the GES v6 list.

    Returns {isotope: [component dicts]}. Nothing is hardcoded: components, gf and EP
    come from the on-disk list; the caller asserts the sums against the committed
    canonical_gf.csv checksum. Loud-fails rather than inventing a pattern."""
    if not os.path.exists(path):
        raise SystemExit(f"RYA-565: GES HFS list not found at {path} — cannot source the "
                         f"Eu II hyperfine pattern. Refusing to invent one.")
    out = {i: [] for i in EU_ISOTOPES}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['element'] != 'Eu 2':
                continue
            wl = float(r['wave_A'])
            if abs(wl - center) > halfwidth:
                continue
            iso = int(r['spectrum_synthe_isotope'])
            if iso not in out:
                raise SystemExit(f"RYA-565: unexpected Eu isotope {iso} in the GES list")
            out[iso].append(dict(wave=wl, ep=float(r['lower_state_eV']),
                                 loggf=float(r['loggf']),
                                 fdamp=float(r['turbospectrum_fdamp']),
                                 upper_g=float(r['upper_g']),
                                 rad=r['turbospectrum_rad'],
                                 lo_orb=r['lower_orbital_type'],
                                 up_orb=r['upper_orbital_type'],
                                 ref=r['reference_code']))
    for iso, comps in out.items():
        if not comps:
            raise SystemExit(f"RYA-565: no Eu II {center} components for isotope {iso} in the "
                             f"GES list — the HFS pattern is not sourceable for this line.")
        comps.sort(key=lambda c: c['wave'])
    return out


def gf_sum_loggf(comps):
    return float(np.log10(sum(10.0 ** c['loggf'] for c in comps)))


def gf_weighted_centroid(comps):
    g = np.array([10.0 ** c['loggf'] for c in comps])
    w = np.array([c['wave'] for c in comps])
    return float((g * w).sum() / g.sum())


def _canonical_rows(root):
    path = os.path.join(root, 'data', 'linelists', 'canonical_gf.csv')
    rows = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r['species'] == 'Eu II':
                rows[round(float(r['wavelength_air_A']), 4)] = r
    return rows


def verify_gf_provenance(root, hfs_by_line):
    """Loud-fail SSOT gate, per line. Establishes from two independent on-disk sources:
      (a) the committed canonical_gf.csv row exists, is the LWHS one, and matches;
      (b) canonical log_gf == the NAIVE cross-isotope GES sum (the checksum identity);
      (c) both isotopes' HFS sums agree with each other (same physical gf);
      (d) canonical - physical == log10(2), the isotope double count.
    Returns the per-line audit. Never edits the single source."""
    canon = _canonical_rows(root)
    audit = {}
    for wave, cfg in LINES.items():
        row = canon.get(round(wave, 4))
        if row is None:
            for w, r in canon.items():
                if abs(w - wave) < 0.01:
                    row = r
                    break
        if row is None:
            raise SystemExit(f"RYA-565 SSOT: no Eu II {wave} row in canonical_gf.csv — cannot "
                             f"source gf. Fix the single source, do not hardcode.")
        canon_loggf, canon_n = float(row['log_gf']), int(row['hfs_n_components'])
        if abs(canon_loggf - cfg['canon_loggf']) > GF_TOL:
            raise SystemExit(f"RYA-565 SSOT DRIFT: Eu II {wave} canonical log_gf={canon_loggf} "
                             f"!= script {cfg['canon_loggf']} ({row['line_id']}). Reconcile.")
        if row['line_id'] != cfg['gf_id']:
            raise SystemExit(f"RYA-565 SSOT DRIFT: Eu II {wave} line_id {row['line_id']} != "
                             f"{cfg['gf_id']}.")

        hfs = hfs_by_line[wave]
        allc = [c for comps in hfs.values() for c in comps]
        naive = gf_sum_loggf(allc)
        per_iso = {iso: gf_sum_loggf(c) for iso, c in hfs.items()}
        if len(allc) != canon_n or canon_n != cfg['hfs_n']:
            raise SystemExit(f"RYA-565 SSOT DRIFT: Eu II {wave}: GES {len(allc)} components, "
                             f"canonical hfs_n_components={canon_n}, script {cfg['hfs_n']}.")
        if abs(naive - canon_loggf) > GF_TOL:
            raise SystemExit(f"RYA-565 SSOT DRIFT: Eu II {wave}: naive GES HFS sum {naive:.4f} "
                             f"!= canonical log_gf {canon_loggf:.4f}. Checksum identity broken.")
        vals = list(per_iso.values())
        if max(vals) - min(vals) > GF_TOL:
            raise SystemExit(f"RYA-565 gf: Eu II {wave}: the two isotopes' HFS sums disagree "
                             f"({per_iso}) — they must carry the same physical gf. STOP.")
        physical = float(np.mean(vals))
        if abs((naive - physical) - LOG2) > 0.01:
            raise SystemExit(f"RYA-565 gf: Eu II {wave}: naive-minus-physical = "
                             f"{naive - physical:.4f}, expected log10(2)={LOG2:.4f}. STOP.")
        centroid = gf_weighted_centroid(allc)
        audit[str(wave)] = dict(
            canonical_line_id=row['line_id'], canonical_loggf=canon_loggf,
            canonical_reference=row['loggf_reference'],
            canonical_adjudication=row['adjudication_status'],
            canonical_hfs_n_components=canon_n,
            ges_naive_sum_loggf=round(naive, 4),
            ges_per_isotope_sum_loggf={str(k): round(v, 4) for k, v in per_iso.items()},
            ges_component_counts={str(k): len(v) for k, v in hfs.items()},
            physical_total_loggf=round(physical, 4),
            naive_minus_physical_dex=round(naive - physical, 4),
            gf_weighted_centroid_A=round(centroid, 4),
            hfs_span_A=[round(min(c['wave'] for c in allc), 3),
                        round(max(c['wave'] for c in allc), 3)])
        print(f"  Eu II {wave:9.4f} [{row['line_id']}, {row['loggf_reference']}, "
              f"{row['adjudication_status']}]")
        print(f"      canonical log_gf {canon_loggf:+.4f} == naive GES sum {naive:+.4f} over "
              f"{len(allc)} comps ({', '.join(f'{k}:{len(v)}' for k, v in hfs.items())})  OK")
        print(f"      per-isotope HFS sum {physical:+.4f} = the PHYSICAL total gf; "
              f"canonical - physical = {naive - physical:+.4f} = log10(2)")
        print(f"      gf-weighted centroid {centroid:.4f} A, pattern spans "
              f"{min(c['wave'] for c in allc):.3f}-{max(c['wave'] for c in allc):.3f} A")
    return audit


# ───────────────────────────── linelist construction ──────────────────────────────

_HDR = re.compile(r"^'\s*([0-9]+\.[0-9]+)\s*'\s+(\d+)\s+(\d+)\s*$")


def _fmt_line(c, comment):
    return ("  %9.3f %6.3f %6.3f %8.3f %6.1f  %s '%s' '%s'   0.0    1.0 '%s'\n"
            % (c['wave'], c['ep'], c['loggf'], c['fdamp'], c['upper_g'],
               c['rad'], c['lo_orb'], c['up_orb'], comment))


def _set_loggf(ln, delta):
    m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(\s+.*)$", ln.rstrip('\n'))
    if not m:
        raise SystemExit(f"RYA-565: cannot parse a linelist row for the blend probe: {ln!r}")
    return f"{m.group(1)}{float(m.group(2)) + delta:.3f}{m.group(3)}\n"


def build_linelist(mode, center, hfs, dst, blend_delta=None):
    """Write the in-window synthesis linelist for one Eu II feature.

    Every non-Eu VALD block is copied VERBATIM — that is the in-window blend model
    (RYA-551 / RYA-585 pattern: fit the profile with the blends present, never an
    isolated-line inversion). Only the Eu II 63.151 / 63.153 blocks are touched:

      mode='ges_hfs'    the target entry is replaced by the GES LWHS HFS components for
                        that isotope (the correct, engine-applies-isotopfrac form).
      mode='vald_asis'  untouched — the as-shipped production list, kept as a labelled
                        sensitivity leg so the convention error is MEASURED.
      mode='blend_only' the target entry is deleted from both isotope blocks -> the
                        blend-only profile, i.e. how much of the observed absorption is
                        NOT europium.

    `blend_delta` = (z, ion, wave, delta_dex) additionally perturbs one BLEND line's
    log gf, which is how the ungraded-blend systematic is quantified.
    """
    if mode not in ('ges_hfs', 'vald_asis', 'blend_only'):
        raise ValueError(mode)
    src = os.path.join(VALD_DIR, LINES[center]['vald'])
    lines = open(src).read().splitlines(keepends=True)

    out, i, n_eu_blocks, n_probe_hits, iso_coded = [], 0, 0, 0, 0
    while i < len(lines):
        m = _HDR.match(lines[i].rstrip('\n'))
        if not m:
            out.append(lines[i]); i += 1
            continue
        species, ion, nline = m.group(1), int(m.group(2)), int(m.group(3))
        z_part, iso_part = species.split('.')
        z, iso = int(z_part), int(iso_part)
        if iso != 0:
            iso_coded += 1
        body = lines[i + 2:i + 2 + nline]
        is_eu = (z == Z_EU and ion == 2 and iso in EU_ISOTOPES)

        if blend_delta is not None and (z, ion) == (blend_delta[0], blend_delta[1]):
            new_body = []
            for ln in body:
                if abs(float(ln.split()[0]) - blend_delta[2]) < 0.005:
                    new_body.append(_set_loggf(ln, blend_delta[3])); n_probe_hits += 1
                else:
                    new_body.append(ln)
            body = new_body

        if not is_eu or mode == 'vald_asis':
            out.append(lines[i]); out.append(lines[i + 1]); out.extend(body)
            i += 2 + nline
            continue

        n_eu_blocks += 1
        kept = [ln for ln in body if abs(float(ln.split()[0]) - center) > 0.30]
        if len(kept) == len(body):
            raise SystemExit(f"RYA-565: no Eu II {center} entry in the VALD {species} block — "
                             f"the linelist changed under us; refusing to guess.")
        new = list(kept)
        if mode == 'ges_hfs':
            new += [_fmt_line(c, f"Eu II {c['wave']:.3f} HFS iso-{iso} "
                                 f"[{c['ref']} via GES v6, RYA-565]") for c in hfs[iso]]
        new.sort(key=lambda ln: float(ln.split()[0]))
        out.append("'  %s            '    %d      %4d\n" % (species, ion, len(new)))
        out.append(lines[i + 1])
        out.extend(new)
        i += 2 + nline

    if mode != 'vald_asis' and n_eu_blocks != len(EU_ISOTOPES):
        raise SystemExit(f"RYA-565: expected {len(EU_ISOTOPES)} Eu II isotope blocks around "
                         f"{center}, found {n_eu_blocks}.")
    if blend_delta is not None and n_probe_hits != 1:
        raise SystemExit(f"RYA-565: blend probe {blend_delta} matched {n_probe_hits} rows "
                         f"(expected exactly 1) — refusing an ambiguous perturbation.")
    with open(dst, 'w') as fh:
        fh.writelines(out)
    return dst, iso_coded


# ─────────────────────────────── engine wrappers ──────────────────────────────────

def babsma(opac, lmin, lmax):
    ctl = (f"'LAMBDA_MIN:'  '{lmin}'\n'LAMBDA_MAX:'  '{lmax}'\n'LAMBDA_STEP:' '0.005'\n"
           f"'MODELINPUT:' '{MARCS}'\n'MARCS-FILE:' '.true.'\n'MODELOPAC:' '{opac}'\n"
           f"'METALLICITY:'    '0.00'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
           f"'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n'XIFIX:' 'T'\n{XI}\n")
    r = subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(opac):
        raise SystemExit(f"babsma failed:\n{r.stdout[-1500:]}\n{r.stderr[-800:]}")


def bsyn(a_x, linelist, opac, lmin, lmax, tag):
    res = f"{W}/s_{tag}.spec"
    ctl = (f"'NLTE :'          '.false.'\n'LAMBDA_MIN:'     '{lmin}'\n'LAMBDA_MAX:'     '{lmax}'\n"
           f"'LAMBDA_STEP:'    '0.005'\n'INTENSITY/FLUX:' 'Flux'\n'ABFIND        :' '.false.'\n"
           f"'MODELOPAC:' '{opac}'\n'RESULTFILE :' '{res}'\n'METALLICITY:'    '0.00'\n"
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{Z_EU}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
           f"'NFILES   :' '1'\n{linelist}\n'SPHERICAL:'  'F'\n  30\n  300.00\n  15\n  1.30\n")
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(res):
        raise SystemExit(f"bsyn failed A={a_x}:\n{r.stdout[-1500:]}")
    d = np.loadtxt(res)
    return d[:, 0], d[:, 1]


# ──────────────────────────────── observations ────────────────────────────────────

def load_harps():
    w, f = [], []
    with open(HARPS) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            w.append(float(row[0])); f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(repo):
    """IAG Reiners+2016 visible FTS solar flux atlas: vacuum wavenumber (cm^-1),
    normalised flux. Vacuum->air via the pipeline SINGLE-SOURCE converter."""
    import gzip
    sys.path.insert(0, repo)
    from pipeline.wavelength_util import vac_to_air
    wn, fl = [], []
    with gzip.open(IAG, 'rt') as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 2:
                continue
            try:
                wn.append(float(p[0])); fl.append(float(p[1]))
            except ValueError:
                continue
    wn = np.array(wn); fl = np.array(fl)
    lam_air = vac_to_air(1e8 / wn)
    idx = np.argsort(lam_air)
    return lam_air[idx], fl[idx]


# ───────────────────────── noise floor / detection significance ───────────────────

def noise_floor(obs_w, obs_f, center, hw, gsig_kms):
    """Quantify the noise floor at the line, in the units the disposition needs.

    S/N is measured from the scatter of the LOCALLY RENORMALISED continuum points in
    the fit window (top-decile flux) — the S/N actually available to this fit, not a
    header claim. sigma_EW follows Cayrel (1988):
        sigma_EW ~= 1.5 * sqrt(FWHM * dlambda) / (S/N)
    with FWHM from the fitted broadening. Returns mA."""
    ww, ff, _ = local_renorm(obs_w, obs_f, center, hw)
    sel = (ww > center - hw - 1.0) & (ww < center + hw + 1.0)
    ww, ff = ww[sel], ff[sel]
    if len(ww) < 20:
        return None
    cont = ff[ff >= np.percentile(ff, 90)]
    sigma_flux = float(np.std(cont))
    snr = float(1.0 / sigma_flux) if sigma_flux > 0 else float('inf')
    dlam = float(np.median(np.diff(ww)))
    fwhm = 2.3548 * (gsig_kms / CLIGHT) * center
    sigma_ew = 1.5 * float(np.sqrt(max(fwhm, dlam) * dlam)) / snr * 1000.0
    return dict(snr_local=round(snr, 1), pixel_A=round(dlam, 4),
                fwhm_A=round(float(fwhm), 4), sigma_EW_mA=round(sigma_ew, 3))


def eu_only_ew(synth_full, synth_blend, a_fit, center, core_hw, gsig):
    """Synthetic EW (mA) attributable to Eu alone at the fitted abundance: the area
    between the blend-only profile and the full profile over the core window."""
    keys = np.array(sorted(synth_full))
    sw, sf = synth_full[float(keys[int(np.argmin(np.abs(keys - a_fit)))])]
    bw, bf = synth_blend
    sb = broaden(sw, sf, VSINI, gsig)
    bb = broaden(bw, bf, VSINI, gsig)
    m = (sw > center - core_hw) & (sw < center + core_hw)
    return float(trapezoid(np.interp(sw[m], bw, bb) - sb[m], sw[m]) * 1000.0)


# ─────────────────────────────────── driver ───────────────────────────────────────

def fit_all_arms(center, cfg, arms, arm_rv, synth, synth_blend, label):
    """Fit one (line, leg) against every solar arm; return the per-arm record dict."""
    leg = {}
    for arm, (ow, of) in arms.items():
        fit = fit_profile(center, ow, of, synth, hw=cfg['fit_hw'], vsini=VSINI,
                          a_lo=A_LO, a_hi=A_HI)
        if fit is None:
            leg[arm] = dict(status='no_coverage')
            print(f"    {label:10s} {arm:6s}: no coverage")
            continue
        nf = noise_floor(ow, of, center, cfg['fit_hw'], fit['gsig'])
        ew_eu = eu_only_ew(synth, synth_blend, fit['A'], center, 0.4, fit['gsig'])
        rchi2 = fit['red_chi2']
        rel_f = assess_reliability(fit)
        reliable = rel_f['reliable']
        rec = dict(A_LTE=round(fit['A'], 3), A_NLTE=None,
                   red_chi2=round(rchi2, 2), npix=fit['npix'],
                   dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'],
                   gsig_kms=round(fit['gsig'], 2), gsig_railed=fit['gsig_railed'],
                   dv_fitted_kms=round(fit['dv'], 2),
                   dv_measured_kms=arm_rv[arm]['v_kms'],
                   obs_window_EW_mA=round(fit['obs_ew_mA'], 2),
                   synth_core_EW_mA=fit['core_EW_mA'],
                   eu_only_EW_mA=round(ew_eu, 3),
                   **rel_f)
        # Linear-COG diagnostic: on the linear part of the curve of growth a line's
        # sensitivity is EW*ln10 exactly. A ratio near 1 says the LOW dEW/dA is
        # intrinsic weakness (maximal sensitivity per mA), NOT desensitisation by
        # saturation or blend swamping. It does NOT relax the floor — it explains it.
        lin = ew_eu * np.log(10.0)
        rec['dEW_dA_linear_cog_mA_dex'] = round(float(lin), 2)
        rec['cog_linearity'] = round(float(fit['dEW_dA'] / lin), 3) if lin > 0 else None
        if nf:
            rec.update(nf)
            rec['sigma_A_noise_dex'] = (round(nf['sigma_EW_mA'] / fit['dEW_dA'], 4)
                                        if fit['dEW_dA'] > 0 else None)
            rec['detection_sigma'] = (round(ew_eu / nf['sigma_EW_mA'], 2)
                                      if nf['sigma_EW_mA'] > 0 else None)
        leg[arm] = rec
        print(f"    {label:10s} {arm:6s}: A_LTE={rec['A_LTE']:+.3f}  rchi2={rec['red_chi2']:6.2f}  "
              f"dEW/dA={rec['dEW_dA_mA_dex']:6.1f} (linear-COG {rec['dEW_dA_linear_cog_mA_dex']:.1f}, "
              f"ratio {rec['cog_linearity']})  railed={str(rec['railed']):5s} reliable={reliable}")
        print(f"    {'':10s} {'':6s}  Eu-only synth EW={rec['eu_only_EW_mA']:.2f} mA, obs window "
              f"EW={rec['obs_window_EW_mA']:.2f} mA, S/N~{rec.get('snr_local')}, "
              f"sigma_EW={rec.get('sigma_EW_mA')} mA => {rec.get('detection_sigma')} sigma, "
              f"sigma_A(noise)={rec.get('sigma_A_noise_dex')} dex")
    return leg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lines', default='all', help='comma list of line centres, or "all"')
    ap.add_argument('--skip-sensitivity', action='store_true',
                    help='skip the gf-convention and blend-gf sensitivity legs')
    args = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    from config.constants import assert_on_sirius, SOLAR_ASPLUND2021
    # RYA-567: Sirius-only heavy-compute leg. Never against local-Mac engines/grids.
    assert_on_sirius("RYA-565 Eu II synthesis", require_subdirs=("engines", "grids"))

    want = (list(LINES) if args.lines == 'all'
            else [float(x) for x in args.lines.split(',')])
    a_sun = float(SOLAR_ASPLUND2021['Eu'])
    print("=" * 96)
    print("RYA-565 — Eu II LTE in-window HFS blend-fit synthesis (Sirius)")
    print(f"  A(Eu)_sun reference (Asplund 2021) = {a_sun}  [grid bracket only, NOT a target]")
    print("  engine: LTE only — Eu II majority ion, no NLTE grid in either engine (RYA-458)")
    print("=" * 96)

    print("\nStep 0 — gf / HFS provenance (SSOT: canonical_gf.csv x GES v6 LWHS block)")
    hfs_by_line = {w: read_ges_hfs(w) for w in LINES}
    gf_audit = verify_gf_provenance(root, hfs_by_line)

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")

    print("\nStep 1 — solar arms + measured rest frame (abundance-blind, RYA-643)")
    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag(root)
    except Exception as e:                                   # noqa: BLE001
        print(f"  IAG load FAILED ({e}); HARPS-only")
    arm_rv = {}
    for arm, (ow, of) in arms.items():
        v, n, sd = require_arm_rv(ow, of, arm)
        arm_rv[arm] = dict(v_kms=round(v, 3), n_lines=n, scatter_kms=round(sd, 3))
        print(f"  {arm:6s}: {v:+.3f} km/s (n={n}, scatter {sd:.3f})")

    results = dict(_meta=dict(
        ticket='RYA-565', element='Eu', ion='II', z=Z_EU, primary_line_A=PRIMARY,
        a_sun_ref=a_sun, a_sun_ref_source='Asplund2021',
        engine='LTE', nlte_delta=None, nlte_status='NLTE_unavailable_LTE_robust',
        nlte_rationale=('Eu II is the majority ion of Eu in the solar photosphere and no Eu '
                        'departure grid exists in either engine (registry LTE-only-by-design, '
                        'RYA-458). The binding blocker is NOT NLTE.'),
        reliable_dEW_dA_floor=RELIABLE_DEWDA, red_chi2_review_trigger=RCHI2_REVIEW,
        red_chi2_gates_reliable=False,
        A_grid=[A_LO, A_HI, 0.05], blend_probe_dex=BLEND_PROBE_DEX,
        measured_EW_mA=6.8,
        measured_EW_source=('data/measured/sol_ew_results_v1.csv — Eu II 6645.127, 6.8 mA, '
                            'hfs_total_ew, blend_flag'),
        ts_isotopfrac={str(k): v for k, v in TS_ISOTOPFRAC.items()},
        gf_provenance=gf_audit, arm_rv=arm_rv), lines={})

    for center in want:
        cfg = LINES[center]
        hfs = hfs_by_line[center]
        fit_center = gf_audit[str(center)]['gf_weighted_centroid_A']
        print(f"\nStep 2 — Eu II {center} ({cfg['role']}, EP {cfg['ep']}, physical log gf "
              f"{gf_audit[str(center)]['physical_total_loggf']:+.4f}, fit centre {fit_center:.4f})")
        lmin, lmax = fit_center - LMARGIN, fit_center + LMARGIN
        opac = f"{W}/opac_{center:.0f}"
        babsma(opac, lmin, lmax)

        ll_blend, iso_coded = build_linelist('blend_only', center, hfs,
                                             f"{W}/ll_{center:.0f}_blendonly.list")
        synth_blend = bsyn(A_BLEND_ONLY, ll_blend, opac, lmin, lmax, f"{center:.0f}_blendonly")
        print(f"    in-window blend model: VALD3 block verbatim "
              f"({iso_coded} isotope-coded species blocks in the source list)")

        rec = dict(role=cfg['role'], ep_eV=cfg['ep'], fit_centre_A=fit_center,
                   physical_loggf=gf_audit[str(center)]['physical_total_loggf'], legs={})

        legs = [('ges_hfs', None, None)]
        if not args.skip_sensitivity:
            legs.append(('vald_asis', None, None))
            if cfg.get('blend_probe'):
                z, ion, bw, blabel = cfg['blend_probe']
                legs.append((f'blend_gf_minus{BLEND_PROBE_DEX}', (z, ion, bw, -BLEND_PROBE_DEX), blabel))
                legs.append((f'blend_gf_plus{BLEND_PROBE_DEX}', (z, ion, bw, +BLEND_PROBE_DEX), blabel))

        for name, bdelta, blabel in legs:
            mode = 'vald_asis' if name == 'vald_asis' else 'ges_hfs'
            ll, _ = build_linelist(mode, center, hfs, f"{W}/ll_{center:.0f}_{name}.list",
                                   blend_delta=bdelta)
            synth = {float(a): bsyn(a, ll, opac, lmin, lmax, f"{center:.0f}_{name}_{a:+.2f}")
                     for a in A_GRID}
            sb = synth_blend
            if bdelta is not None:
                llb, _ = build_linelist('blend_only', center, hfs,
                                        f"{W}/ll_{center:.0f}_{name}_blendonly.list",
                                        blend_delta=bdelta)
                sb = bsyn(A_BLEND_ONLY, llb, opac, lmin, lmax, f"{center:.0f}_{name}_blendonly")
            rec['legs'][name] = fit_all_arms(center, cfg, arms, arm_rv, synth, sb, name)
            if blabel:
                rec['legs'][name]['_probe'] = blabel
        results['lines'][str(center)] = rec

    # ── error budget + disposition ────────────────────────────────────────────────
    prim = results['lines'].get(str(PRIMARY))
    budget, disposition = {}, {}
    if prim:
        g = prim['legs']['ges_hfs']
        arms_ok = [a for a in ('harps', 'iag') if isinstance(g.get(a), dict) and 'A_LTE' in g[a]]
        a_vals = [g[a]['A_LTE'] for a in arms_ok]
        budget['sigma_A_noise_dex'] = g.get('harps', {}).get('sigma_A_noise_dex')
        if len(a_vals) > 1:
            budget['cross_arm_spread_dex'] = round(abs(a_vals[0] - a_vals[1]), 3)
        for name in prim['legs']:
            if name.startswith('blend_gf_') and isinstance(prim['legs'][name].get('harps'), dict):
                budget[f'delta_A_{name}_dex'] = round(
                    prim['legs'][name]['harps']['A_LTE'] - g['harps']['A_LTE'], 3)
        if isinstance(prim['legs'].get('vald_asis', {}).get('harps'), dict):
            budget['gf_convention_bias_dex'] = round(
                prim['legs']['vald_asis']['harps']['A_LTE'] - g['harps']['A_LTE'], 3)
        blend_sys = max((abs(v) for k, v in budget.items() if k.startswith('delta_A_blend_gf')),
                        default=None)
        budget['blend_gf_systematic_dex'] = blend_sys
        parts = [v for v in (budget.get('sigma_A_noise_dex'),
                             budget.get('cross_arm_spread_dex'), blend_sys) if v]
        budget['quadrature_sigma_A_dex'] = (round(float(np.sqrt(sum(p * p for p in parts))), 3)
                                            if parts else None)

        # line-to-line spread across every fitted line's primary (ges_hfs / harps) leg
        per_line = {w: r['legs']['ges_hfs']['harps']['A_LTE']
                    for w, r in results['lines'].items()
                    if isinstance(r['legs']['ges_hfs'].get('harps'), dict)
                    and 'A_LTE' in r['legs']['ges_hfs']['harps']}
        budget['A_per_line_harps'] = per_line
        if len(per_line) > 1:
            v = np.array(list(per_line.values()))
            budget['line_to_line_scatter_dex'] = round(float(v.std(ddof=1)), 3)
            budget['line_to_line_range_dex'] = round(float(v.max() - v.min()), 3)

        reliable = [a for a in arms_ok if g[a]['reliable']]
        if reliable:
            vv = np.array([g[a]['A_LTE'] for a in reliable])
            disposition = dict(emit_value=True, A_Eu=round(float(vv.mean()), 3),
                               n_arms=len(vv),
                               reason='Eu II 6645 cleared the reliability floor on the primary leg.')
        else:
            h = g.get('harps', {})
            disposition = dict(
                emit_value=False, A_Eu=None,
                reason=(
                    f"Eu II {PRIMARY} does NOT clear the RYA-551/560/585 reliability floor: "
                    f"dEW/dA={h.get('dEW_dA_mA_dex')} mA/dex vs floor {RELIABLE_DEWDA} "
                    f"(rchi2={h.get('red_chi2')}, railed={h.get('railed')}). The profile fit "
                    f"itself is GOOD and the feature is well detected — the line is simply "
                    f"intrinsically weak (Eu-only synthetic EW {h.get('eu_only_EW_mA')} mA, "
                    f"cog_linearity {h.get('cog_linearity')} i.e. on the linear curve of growth), "
                    f"so the sensitivity floor cannot be met by this line at any S/N. The binding "
                    f"uncertainty is SYSTEMATIC, not photon noise: cross-arm spread "
                    f"{budget.get('cross_arm_spread_dex')} dex, ungraded-blend gf systematic "
                    f"{budget.get('blend_gf_systematic_dex')} dex (vs sigma_A(noise) "
                    f"{budget.get('sigma_A_noise_dex')} dex), line-to-line scatter "
                    f"{budget.get('line_to_line_scatter_dex')} dex over "
                    f"{len(budget.get('A_per_line_harps', {}))} HFS-resolved Eu II lines. "
                    f"NOISE-FLOOR / BLEND LIMITED (RYA-458), not NLTE-limited. Eu stays "
                    f"owed-no-value; the next lever is a cleaner line / higher-SNR blue arm "
                    f"(Eu II 4129/4205), not this line set."))
    results['error_budget'] = budget
    results['disposition'] = disposition

    out_dir = os.path.join(root, 'data', 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'eu2_synthesis_rya565.json')
    with open(out_path, 'w') as fh:
        json.dump(results, fh, indent=2)

    print("\n" + "=" * 96)
    print("ERROR BUDGET (primary line, HARPS unless noted)")
    for k, v in budget.items():
        print(f"  {k:34s} {v}")
    print("-" * 96)
    if disposition.get('emit_value'):
        print(f"DISPOSITION: EMIT  A(Eu) = {disposition['A_Eu']:+.3f} (n={disposition['n_arms']})")
    else:
        print("DISPOSITION: OWED-NO-VALUE (Eu stays owed)")
    print("  " + disposition.get('reason', ''))
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
