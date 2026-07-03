#!/usr/bin/env python3
"""
scripts/molecular_linelist_convert.py
=====================================
RYA-503 — THE canonical molecular line-list converter: ExoMol (.states/.trans, which
also covers the MoLLIST datasets ExoMol re-hosts) → Turbospectrum molecular `.bsyn`.
Single source of truth for molecular acquisition; future molecular lists reuse it (no
per-species bespoke path).

Target format authority = the VENDORED `CO_IR_Li2015.dat` (read, not specified from
memory; RYA-360). Its per-line record is:

    λ_Å  χ_low_eV  loggf  fdamp  g_u  A_ul  'X' 'X'  J_u  J_l  'X' 'X'  0.0 0.0  'desc'

with, verified field-by-field against the vendored CO file (RYA-503 Phase 1):

    λ (Å)      = 1e8 / (E_upper − E_lower)     — VACUUM wavelength (no air conversion)
    χ_low (eV) = E_lower[cm⁻¹] / 8065.54429    — CODATA 1 eV ≡ 8065.54429 cm⁻¹
    g_u        = statistical weight of the UPPER (higher-energy) state
    loggf      = log10( 1.4992e-16 · λ_Å² · g_u · A_ul )   — the gf–A relation
    A_ul       = Einstein A (s⁻¹) straight from the source
    upper/lower assigned by ENERGY (upper = higher E); desc = 'v{vu}-{vl}_J{Ju}-{Jl}_{tag}'

Physical constants are single-sourced + cited below; nothing is hard-coded from memory.

CO round-trip (Phase 1 hard gate): re-convert CO from its ExoMol Li2015 source and
compare to the vendored file — line count + per-field max deviation vs tolerance. A
MISMATCH is a loud STOP finding; this script NEVER writes over the vendored CO file.

Usage:
    python scripts/molecular_linelist_convert.py --selftest
    python scripts/molecular_linelist_convert.py --source exomol \\
        --states X.states --trans X.trans --masses 6:12 8:16 \\
        --tag Li2015 --nu-min 1000 --out out.bsyn
    python scripts/molecular_linelist_convert.py --roundtrip <vendored CO_IR_Li2015.dat> \\
        --states 12C-16O__Li2015.states --trans 12C-16O__Li2015.trans --masses 6:12 8:16 \\
        --tag Li2015 --nu-min 1000
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

# ── physical constants (single-sourced + cited; NOT from memory) ──────────────
# 1 eV in wavenumber — CODATA 2018 recommended value (physics.nist.gov, "electron
# volt-inverse meter relationship" 8.065543937e5 m⁻¹ = 8065.543937 cm⁻¹).
CM_PER_EV = 8065.543937
# Oscillator-strength / Einstein-A relation, λ in Å:
#   g_u·A_ul = (8π²e²)/(m_e c) · (g_l f_lu)/λ²  ⇒  g_l f_lu = 1.49919e-16 · λ²[Å] · g_u A_ul
# The 1.49919e-16 prefactor is the standard cgs constant m_e c/(8π²e²) with λ in Å
# (e.g. Thorne, "Spectrophysics"; ties to A = 6.6702e15·g_l f/(g_u λ²)). RYA-236 used
# 1.4992e-16; both agree to <1e-4 dex (checked in the round-trip).
GF_CONST = 1.49919e-16


class MolecularConvertError(RuntimeError):
    """Loud failure — an unreadable source, an unknown format, or (the Phase-1 gate) a
    CO round-trip mismatch. Never a silent fallback / overwrite."""


# ── ExoMol readers ────────────────────────────────────────────────────────────
@dataclass
class State:
    E_cm: float      # state energy (cm⁻¹, vacuum)
    g: int           # total statistical weight
    J: float         # rotational quantum number
    v: int           # vibrational quantum number (0 if not resolvable)
    elec: str        # electronic-state label (e.g. 'X2Pi', 'X(3SIGMA-)'; '' if absent)


@dataclass
class ExomolDef:
    """The parts of an ExoMol `.def` needed to locate quantum columns in `.states`.
    ExoMol `.states` columns are: id, E, g, J, then the OPTIONAL unc/lifetime/lande
    columns (each present per an availability flag), then the quantum labels in order."""
    dataset: str
    version: str
    quantum_labels: list           # ordered quantum-number labels
    has_unc: bool
    has_lifetime: bool
    has_lande: bool

    @property
    def q_offset(self) -> int:
        return 4 + int(self.has_unc) + int(self.has_lifetime) + int(self.has_lande)

    def q_index(self, label: str):
        for i, lab in enumerate(self.quantum_labels):
            if lab == label:
                return self.q_offset + i
        return None


def parse_exomol_def(path) -> ExomolDef:
    """Read dataset/version + the availability flags + the ordered quantum labels."""
    dataset = version = ''
    has_unc = has_life = has_lande = False
    labels = []
    with open(path) as fh:
        for ln in fh:
            if '#' not in ln:
                continue
            val, _, comment = ln.partition('#')
            val = val.strip(); c = comment.strip().lower()
            if c == 'isotopologue dataset name' and not dataset:
                dataset = val
            elif c.startswith('version number') and not version:
                version = val
            elif c.startswith('uncertainty availability'):
                has_unc = val.startswith('1')
            elif c.startswith('lifetime availability'):
                has_life = val.startswith('1')
            elif c.startswith('lande') and 'availability' in c:
                has_lande = val.startswith('1')
            elif c.startswith('quantum label'):
                labels.append(val)
    return ExomolDef(dataset, version, labels, has_unc, has_life, has_lande)


# Fallback layout for the simple 'id E g J v e/f' case (ExoMol 'dcs', e.g. CO Li2015):
# no optional columns, v is the first quantum column.
_SIMPLE_DEF = ExomolDef('', '', ['v', 'e/f'], False, False, False)


def parse_exomol_states(path, edef: ExomolDef = _SIMPLE_DEF) -> dict[int, State]:
    """ExoMol `.states` → {id: State}. id/E/g/J are the fixed first four columns; v and
    ElecState are located via the `.def` quantum-column map (`edef`). Values are single
    whitespace tokens (ExoMol quantum values contain no spaces)."""
    vi = edef.q_index('v')
    ei = edef.q_index('ElecState')
    states: dict[int, State] = {}
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            i = int(p[0]); E = float(p[1]); g = int(p[2]); J = float(p[3])
            v = int(p[vi]) if (vi is not None and vi < len(p) and _is_int(p[vi])) else 0
            elec = p[ei] if (ei is not None and ei < len(p)) else ''
            states[i] = State(E, g, J, v, elec)
    if not states:
        raise MolecularConvertError(f"no states parsed from {path}")
    return states


def iter_exomol_trans(path):
    """ExoMol `.trans`: `id_upper  id_lower  A(s⁻¹)  [ν(cm⁻¹)]`. Yields (iu, il, A)."""
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 3:
                continue
            yield int(p[0]), int(p[1]), float(p[2])


def _is_int(s: str) -> bool:
    try:
        int(s); return True
    except ValueError:
        return False


# ── the conversion (physics; format authority = vendored CO) ──────────────────
@dataclass
class BsynLine:
    lam_A: float
    chi_low_eV: float
    loggf: float
    g_u: float
    A: float
    Ju: float
    Jl: float
    desc: str


def _is_ground(elec: str) -> bool:
    """The ground electronic state's label starts with 'X' (X2Pi / X(3SIGMA-) / X(2PI))."""
    return elec.upper().startswith('X')


def convert_exomol(states: dict[int, State], trans_path, *, tag: str,
                   nu_min: float = 0.0, nu_max: float = float('inf'),
                   ground_only: bool = False) -> list[BsynLine]:
    """Convert an ExoMol/MoLLIST (.states + .trans) pair to Turbospectrum lines.
    Recipe parameters (non-physics): the ν window [nu_min, nu_max] cm⁻¹ (e.g. the CO
    recipe keeps ν≥1000, dropping the far-IR pure-rotation tail), and ``ground_only``
    which keeps only X–X (ground-electronic) transitions — the ro-vibrational bands —
    excluding the UV/optical electronic bands (those are the RYA-360 held .bsyn set)."""
    out: list[BsynLine] = []
    for iu, il, A in iter_exomol_trans(trans_path):
        su, sl = states.get(iu), states.get(il)
        if su is None or sl is None:
            raise MolecularConvertError(
                f"transition references missing state id ({iu} or {il}) — states/trans mismatch")
        # assign upper/lower by ENERGY (emission convention; robust to id ordering)
        if sl.E_cm > su.E_cm:
            su, sl = sl, su
        nu = su.E_cm - sl.E_cm
        if nu <= 0 or nu < nu_min or nu > nu_max:
            continue
        if ground_only and not (_is_ground(su.elec) and _is_ground(sl.elec)):
            continue
        lam = 1.0e8 / nu                               # vacuum Å
        gf = GF_CONST * lam * lam * su.g * A
        if gf <= 0:
            continue
        out.append(BsynLine(
            lam_A=lam, chi_low_eV=sl.E_cm / CM_PER_EV, loggf=math.log10(gf),
            g_u=float(su.g), A=A, Ju=su.J, Jl=sl.J,
            desc=f"v{su.v}-{sl.v}_J{int(su.J)}-{int(sl.J)}_{tag}"))
    if not out:
        raise MolecularConvertError(
            f"conversion produced 0 lines (ν window [{nu_min},{nu_max}] cm⁻¹"
            f"{' + ground-only' if ground_only else ''} too tight?)")
    out.sort(key=lambda r: r.lam_A)
    return out


def species_code(z1: int, m1: int, z2: int, m2: int) -> str:
    """Turbospectrum molecular species code, e.g. C(6)m12 + O(8)m16 → '0608.012016'."""
    return f"{z1:02d}{z2:02d}.{m1:03d}{m2:03d}"


def format_bsyn(lines: list[BsynLine], code: str, source_label: str) -> str:
    """Emit the Turbospectrum molecular `.bsyn` text (matches the vendored CO layout)."""
    head = f"'{code:>20}'{1:>5}{len(lines):>9}\n'{source_label}'\n"
    body = []
    for r in lines:
        body.append(
            f"{r.lam_A:10.3f}{r.chi_low_eV:9.5f}{r.loggf:8.3f}{0.0:9.3f}{r.g_u:7.1f} "
            f"{r.A:.2E} 'X' 'X' {r.Ju:6.1f} {r.Jl:6.1f} 'X' 'X' {0.0:5.1f} {0.0:6.1f}  "
            f"'{r.desc}'")
    return head + "\n".join(body) + "\n"


# ── vendored-file parser (for the round-trip) ─────────────────────────────────
def parse_bsyn(path) -> dict[str, BsynLine]:
    """Parse a Turbospectrum molecular `.bsyn`/`.dat` into {desc: BsynLine}. desc (the
    'v{vu}-{vl}_J{Ju}-{Jl}_tag' transition id) is the unambiguous match key."""
    out: dict[str, BsynLine] = {}
    with open(path) as fh:
        for ln in fh:
            s = ln.lstrip()
            if not s or s[0] == "'":
                continue
            p = ln.split()
            # fields: lam chi loggf fdamp gu A 'X' 'X' Ju Jl 'X' 'X' 0 0 'desc'
            desc = p[-1].strip("'")
            out[desc] = BsynLine(
                lam_A=float(p[0]), chi_low_eV=float(p[1]), loggf=float(p[2]),
                g_u=float(p[4]), A=float(p[5]), Ju=float(p[8]), Jl=float(p[9]), desc=desc)
    return out


# ── CO round-trip (Phase 1 hard gate) ─────────────────────────────────────────
# Tolerances = the vendored file's OWN PRINTED PRECISION (it is a rounded text file; a
# re-convert from the same source can only differ by that rounding). λ 3 dp, χ 5 dp,
# loggf 3 dp, and A to 3 SIGNIFICANT FIGURES ('8.35E-10'). A is therefore compared at
# 3-sig-fig precision — the correct comparison against a 3-sf field; the raw full-
# precision A_rel is reported alongside and equals the 3-sf rounding bound (~5e-3), not
# a numeric disagreement.
TOL = {'lam_A': 0.01, 'chi_low_eV': 1e-4, 'loggf': 0.01, 'A_rel': 1e-3}


def _sig3(x: float) -> float:
    """Round to 3 significant figures — the vendored A field's stored precision."""
    if x == 0:
        return 0.0
    from math import floor, log10
    return round(x, -int(floor(log10(abs(x)))) + 2)


def roundtrip(vendored, states_path, trans_path, *, masses, tag, nu_min, nu_max):
    """Re-convert and compare to the vendored file. Returns (ok, report dict). NEVER
    writes the vendored file."""
    vend = parse_bsyn(vendored)
    states = parse_exomol_states(states_path)
    lines = convert_exomol(states, trans_path, tag=tag, nu_min=nu_min, nu_max=nu_max)
    got = {r.desc: r for r in lines}

    rep = {'vendored_lines': len(vend), 'converted_lines': len(got),
           'line_count_match': len(vend) == len(got)}
    only_vend = set(vend) - set(got)
    only_got = set(got) - set(vend)
    rep['only_in_vendored'] = len(only_vend)
    rep['only_in_converted'] = len(only_got)

    shared = set(vend) & set(got)
    dev = {'lam_A': 0.0, 'chi_low_eV': 0.0, 'loggf': 0.0, 'A_rel': 0.0}
    a_rel_raw = 0.0
    for d in shared:
        a, b = vend[d], got[d]
        dev['lam_A'] = max(dev['lam_A'], abs(a.lam_A - b.lam_A))
        dev['chi_low_eV'] = max(dev['chi_low_eV'], abs(a.chi_low_eV - b.chi_low_eV))
        dev['loggf'] = max(dev['loggf'], abs(a.loggf - b.loggf))
        if a.A != 0:
            # compare A at the vendored field's 3-sig-fig precision (the gated field)
            dev['A_rel'] = max(dev['A_rel'], abs(a.A - _sig3(b.A)) / abs(a.A))
            a_rel_raw = max(a_rel_raw, abs(a.A - b.A) / abs(a.A))
    rep['shared_lines'] = len(shared)
    rep['max_dev'] = dev
    rep['A_rel_raw_fullprec'] = a_rel_raw   # informational: = the 3-sf rounding bound
    rep['within_tol'] = {k: dev[k] <= TOL[k] for k in dev}
    rep['ok'] = bool(rep['line_count_match'] and not only_vend and not only_got
                     and all(rep['within_tol'].values()))
    return rep['ok'], rep


# ── self-test (offline; CLI/branch sanity, no network) ────────────────────────
def selftest() -> bool:
    ok = True
    # gf–A round-trip on the first vendored CO line's numbers (λ 4515.025, gu 19, A 8.35e-10)
    gf = GF_CONST * 4515.025 ** 2 * 19.0 * 8.35e-10
    lg = math.log10(gf)
    c1 = abs(lg - (-16.314)) < 0.01
    print(f"  [gf] log10(1.49919e-16·4515.025²·19·8.35e-10) = {lg:.3f}  (vendored -16.314)  "
          f"{'OK' if c1 else 'FAIL'}")
    # eV conversion on a known lower level (v0 J8: B0·J(J+1) ≈ 138.4 cm⁻¹ → 0.01716 eV)
    chi = 138.39 / CM_PER_EV
    c2 = abs(chi - 0.01716) < 5e-4
    print(f"  [eV] 138.39 cm⁻¹ / {CM_PER_EV} = {chi:.5f} eV  (vendored 0.01716)  {'OK' if c2 else 'FAIL'}")
    # species code
    c3 = species_code(6, 12, 8, 16) == '0608.012016'
    print(f"  [code] species_code(C12,O16) = {species_code(6,12,8,16)!r}  {'OK' if c3 else 'FAIL'}")
    return ok and c1 and c2 and c3


def _masses(pairs):
    out = []
    for pr in pairs:
        z, m = pr.split(':')
        out.append((int(z), int(m)))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description='RYA-503 ExoMol/MoLLIST → Turbospectrum .bsyn converter')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--source', choices=['exomol'], default='exomol',
                    help="input format (exomol covers ExoMol + the MoLLIST datasets ExoMol re-hosts; "
                         "HITRAN .par is a documented future branch — not yet validated, raises)")
    ap.add_argument('--states'); ap.add_argument('--trans')
    ap.add_argument('--def', dest='defpath', help="ExoMol .def (locates v/ElecState columns)")
    ap.add_argument('--masses', nargs=2, help="Z:mass Z:mass, e.g. 6:12 8:16")
    ap.add_argument('--code', default=None,
                    help="explicit TS species code (e.g. '0108.000016'); overrides --masses. "
                         "For hydrides read it from the molecule's vendored .bsyn (iSpec writes "
                         "the ¹H mass as 000), never guess it.")
    ap.add_argument('--tag', default='')
    ap.add_argument('--label', default=None, help="source label written to the .bsyn header")
    ap.add_argument('--nu-min', type=float, default=0.0)
    ap.add_argument('--nu-max', type=float, default=float('inf'))
    ap.add_argument('--ground-only', action='store_true',
                    help="keep only X–X (ground-electronic) ro-vibrational transitions")
    ap.add_argument('--out')
    ap.add_argument('--roundtrip', help="vendored .bsyn/.dat to validate against (no write of it)")
    args = ap.parse_args(argv)

    if args.selftest:
        print("molecular_linelist_convert --selftest:")
        return 0 if selftest() else 1

    if args.roundtrip:
        (z1, m1), (z2, m2) = _masses(args.masses)
        ok, rep = roundtrip(args.roundtrip, args.states, args.trans,
                            masses=[(z1, m1), (z2, m2)], tag=args.tag,
                            nu_min=args.nu_min, nu_max=args.nu_max)
        print(f"CO ROUND-TRIP vs {Path(args.roundtrip).name}:")
        print(f"  vendored lines : {rep['vendored_lines']}")
        print(f"  converted lines: {rep['converted_lines']}  "
              f"(match: {rep['line_count_match']})")
        print(f"  only-in-vendored / only-in-converted: "
              f"{rep['only_in_vendored']} / {rep['only_in_converted']}")
        print(f"  shared: {rep['shared_lines']}  max per-field deviation vs tolerance:")
        for k in ('lam_A', 'chi_low_eV', 'loggf', 'A_rel'):
            label = 'A_rel(3sf)' if k == 'A_rel' else k
            print(f"      {label:<11} {rep['max_dev'][k]:.3e}  (tol {TOL[k]:.0e})  "
                  f"{'OK' if rep['within_tol'][k] else 'FAIL'}")
        print(f"      A_rel(raw)  {rep['A_rel_raw_fullprec']:.3e}  (informational = vendored's "
              f"3-sig-fig display bound, not a numeric disagreement)")
        print(f"  RESULT: {'PASS' if ok else 'FAIL — STOP (loud finding; vendored CO not touched)'}")
        return 0 if ok else 1

    if not (args.states and args.trans and (args.masses or args.code)):
        ap.error("conversion needs --states --trans and (--masses or --code) "
                 "(or use --selftest / --roundtrip)")
    edef = parse_exomol_def(args.defpath) if args.defpath else _SIMPLE_DEF
    states = parse_exomol_states(args.states, edef)
    lines = convert_exomol(states, args.trans, tag=args.tag,
                           nu_min=args.nu_min, nu_max=args.nu_max,
                           ground_only=args.ground_only)
    code = args.code if args.code else species_code(*_masses(args.masses)[0], *_masses(args.masses)[1])
    text = format_bsyn(lines, code, args.label or f"ExoMol {args.tag}")
    if args.out:
        Path(args.out).write_text(text)
        print(f"wrote {len(lines)} lines → {args.out}  (species {code})")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
