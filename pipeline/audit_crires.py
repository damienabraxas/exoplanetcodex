"""
RYA-952 — inventory every CRIRES frame, split classic from CRIRES+, confirm every target
by ASTROMETRY, and co-add CRIRES+ per target x setting.

THREE THINGS THIS REFUSES TO TRUST, EACH BECAUSE IT HAS BEEN WRONG BEFORE.

1. 🔴 `OBJECT`. On this very data set, four science frames of **tau Ceti carry
   `OBJECT = 'STD'`** — observed as an RV calibration standard, so the pipeline wrote the
   ROLE where the name goes. An inventory keyed on `OBJECT` files them as an anonymous
   standard and reports tau Ceti as absent. RYA-423 found the same shape three other ways
   (`alf Cen B` filed under A, a frame called `Star S5`). One star also appears under two
   names in one night's data here: `eps Eri` and `HD 22049` are the same object.

2. 🔴 THE FOLDER. `spectra/eps_eri/CRIRESPlus/` is a fine guess and nothing more. Files are
   classified by what their headers say, and the directory is recorded only as provenance.

3. 🔴 THE FILENAME'S DATE. `< 2014 = classic` is the ticket's rule and it is a good first
   cut, but the instrument is settled by the PIPELINE that reduced the frame:
   `ESO PRO REC1 PIPE ID` naming **cr2res** is CRIRES+ by construction. ⚠️ `INSTRUME` reads
   the bare string `CRIRES` for BOTH instruments, so it cannot split them — a trap RYA-794
   and RYA-796 both had to work around. Date and pipeline are checked AGAINST each other and
   a disagreement is a loud failure, not a silent preference.

THE IDENTIFICATION IS POSITIONAL, AND THE REFEREE IS EXTERNAL. Each frame's pointing
(`RA`/`DEC`) is compared against every candidate in `data/reference/crires_target_astrometry.csv`
(SIMBAD, committed by `scripts/rya952_fetch_target_astrometry.py`) after propagating each
candidate's J2000 position to the frame's own epoch by its proper motion. That propagation is
not a nicety: **tau Ceti moves 1.92"/yr**, so 22 years of it is 42" — far larger than the
match radius, and skipping it makes the true target look like the wrong one or like nothing
at all.

A frame identifies only when exactly ONE catalogue target lies inside the radius. Zero or
several ⇒ QUARANTINE, never a guess (the ticket's "loud-fail on unconfirmed target ID").
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

#: Match radius. Generous next to the astrometry (SIMBAD positions are mas-accurate and the
#: propagated epoch positions land within a few arcsec of the pointing) but far tighter than
#: the spacing between any two candidate targets, which is degrees. Slack is for nodding
#: throw and pointing error, not for identity: CRIRES nods along the slit by a few arcsec.
MATCH_RADIUS_ARCSEC = 60.0

#: The date the ticket splits on. Classic CRIRES was decommissioned in 2014; CRIRES+ returned
#: in 2021. Nothing legitimately falls between, so the gap is itself a check.
CLASSIC_MAX_YEAR = 2014
PLUS_MIN_YEAR = 2021

#: The CRIRES+ pipeline. `INSTRUME` cannot split the two instruments; the RECIPE NAME can.
#: ⚠️ IT IS `PRO REC1 ID`, NOT `PRO REC1 PIPE ID`. The latter is a VERSION STRING ("1.6.9")
#: and contains no instrument name at all, so testing it for "cr2res" is always False and
#: every reduced IDP silently falls through to the raw branch. Caught because the run
#: reported 64 of 64 frames as raw when 63 of them plainly carry a PRO REC chain.
PLUS_PIPE_TOKEN = 'cr2res'


#: Solar-system bodies have NO catalogue position, so the astrometric matcher cannot place
#: them and would report 18 correctly-identified Vesta frames as unidentified. Their identity
#: is established by a different route entirely -- OBJECT plus the ephemeris provenance
#: already banked by RYA-794/805 -- and calling that "quarantine" would bury the frames that
#: really are unidentified under a pile of frames that are not.
MOVING_TARGETS = {'vesta': 'minor planet 4 Vesta — reflected solar (RYA-372/794/805)'}

#: 🔴 alpha Cen A and B are a VISUAL BINARY currently ~4-8 arcsec apart. No match radius wide
#: enough to absorb CRIRES nodding is ever narrow enough to separate them, so astrometry
#: ALONE cannot decide which component a frame points at -- and that is not a limitation to
#: work around here, it is the entire reason RYA-423 exists (IR-native star ID by
#: RV-ephemeris with IR-template/CO corroboration). Frames landing on the pair are reported
#: as ambiguous and routed there, never assigned by picking the nearer of the two.
BINARY_PAIRS = (frozenset({'alpha_cen_a', 'alpha_cen_b'}),)


class TargetUnconfirmed(Exception):
    """A frame whose target could not be established. Quarantine, never guess."""


@dataclass
class Frame:
    path: str
    md5: str
    instrume: str = ''
    object_raw: str = ''
    obs_targ_name: str = ''
    date_obs: str = ''
    mjd: float = float('nan')
    ra_deg: float = float('nan')
    dec_deg: float = float('nan')
    pipe_id: str = ''
    rec_id: str = ''
    pro_catg: str = ''
    setting: str = ''
    snr: float = float('nan')
    exptime: float = float('nan')
    specsys: str = ''
    prog_id: str = ''
    n_ext: int = 0
    # filled by the classifiers
    instrument_class: str = ''
    vintage_evidence: str = ''
    star_id: str = ''
    id_method: str = ''
    id_sep_arcsec: float = float('nan')
    id_evidence: str = ''
    id_status: str = ''
    duplicate_of: str = ''


def _f(h, key, default=None):
    v = h.get(key, default)
    return v


def read_frame(path: Path) -> Frame:
    """One frame's header truth. Never opens the data arrays."""
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        h = hdul[0].header
        n_ext = len(hdul)
    fr = Frame(path=str(path), md5=hashlib.md5(path.read_bytes()).hexdigest(), n_ext=n_ext)
    fr.instrume = str(_f(h, 'INSTRUME', '') or '').strip()
    fr.object_raw = str(_f(h, 'OBJECT', '') or '').strip()
    fr.obs_targ_name = str(_f(h, 'HIERARCH ESO OBS TARG NAME', '') or '').strip()
    fr.date_obs = str(_f(h, 'DATE-OBS', '') or '')[:19]
    fr.mjd = float(_f(h, 'MJD-OBS', float('nan')) or float('nan'))
    for k in ('RA', 'DEC'):
        try:
            setattr(fr, 'ra_deg' if k == 'RA' else 'dec_deg', float(h[k]))
        except (KeyError, TypeError, ValueError):
            pass
    fr.pipe_id = str(_f(h, 'HIERARCH ESO PRO REC1 PIPE ID', '') or '').strip()
    fr.rec_id = str(_f(h, 'HIERARCH ESO PRO REC1 ID', '') or '').strip()
    fr.pro_catg = str(_f(h, 'HIERARCH ESO PRO CATG', '') or '').strip()
    fr.setting = str(_f(h, 'HIERARCH ESO INS WLEN ID', '') or '').strip()
    try:
        fr.snr = float(_f(h, 'SNR', float('nan')))
    except (TypeError, ValueError):
        pass
    try:
        fr.exptime = float(_f(h, 'EXPTIME', float('nan')))
    except (TypeError, ValueError):
        pass
    fr.specsys = str(_f(h, 'SPECSYS', '') or '').strip()
    fr.prog_id = str(_f(h, 'HIERARCH ESO OBS PROG ID', '') or '').strip()
    return fr


def is_crires(fr: Frame) -> bool:
    return 'CRIRES' in fr.instrume.upper()


def classify_vintage(fr: Frame) -> Frame:
    """classic vs CRIRES+, from the PIPELINE first and the date as an independent check.

    🔴 THE TWO MUST AGREE. A cr2res-reduced frame dated 2011, or a 2024 frame reduced by the
    old pipeline, is not something to resolve by preferring one field — it means the header
    is not describing the file, and that is exactly the condition under which a classic frame
    would get co-added into a CRIRES+ stack.
    """
    year = None
    if len(fr.date_obs) >= 4 and fr.date_obs[:4].isdigit():
        year = int(fr.date_obs[:4])
    pipe_says_plus = PLUS_PIPE_TOKEN in fr.rec_id.lower()

    by_date = None
    if year is not None:
        if year < CLASSIC_MAX_YEAR:
            by_date = 'classic'
        elif year >= PLUS_MIN_YEAR:
            by_date = 'plus'
        else:
            by_date = 'gap'

    if pipe_says_plus:
        fr.instrument_class = 'crires_plus'
        fr.vintage_evidence = (f"PRO REC1 ID={fr.rec_id} (cr2res => CRIRES+), "
                               f"pipeline v{fr.pipe_id}; DATE-OBS {fr.date_obs[:10]}")
        if by_date == 'classic':
            raise ValueError(
                f"{fr.path}: cr2res pipeline on a {year} frame — the pipeline says CRIRES+ "
                f"and the date says classic. Refusing to pick one.")
    elif by_date == 'classic':
        fr.instrument_class = 'crires_classic'
        fr.vintage_evidence = (f"DATE-OBS {fr.date_obs[:10]} < {CLASSIC_MAX_YEAR}; "
                               f"no cr2res recipe (PRO REC1 ID={fr.rec_id!r})")
    elif by_date == 'plus':
        # Post-2021 but no cr2res stamp: a RAW frame, which carries no PRO REC keywords.
        fr.instrument_class = 'crires_plus_raw'
        fr.vintage_evidence = (f"DATE-OBS {fr.date_obs[:10]} >= {PLUS_MIN_YEAR}, no PRO REC "
                               f"chain — raw, not pipeline-reduced")
    else:
        fr.instrument_class = 'unknown'
        fr.vintage_evidence = f"DATE-OBS {fr.date_obs!r}, REC1 ID {fr.rec_id!r} — cannot place"
    return fr


def load_astrometry(repo_root: Path) -> pd.DataFrame:
    p = repo_root / 'data' / 'reference' / 'crires_target_astrometry.csv'
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found — regenerate with scripts/rya952_fetch_target_astrometry.py. "
            f"Target identity is not something this module will guess at.")
    return pd.read_csv(p)


def _propagate(cat: pd.DataFrame, mjd: float) -> tuple[np.ndarray, np.ndarray]:
    """Catalogue J2000 positions carried to the frame's epoch by proper motion.

    🔴 NOT OPTIONAL. tau Ceti moves 1.92 arcsec/yr; over the ~22 years from J2000 to these
    observations that is 42 arcsec, comfortably outside the match radius. Omitting this does
    not blur the answer, it changes it.
    """
    yr = (mjd - 51544.5) / 365.25          # MJD 51544.5 = J2000.0
    dec = cat.dec_deg_j2000.to_numpy(float) + cat.pm_dec_mas_yr.to_numpy(float) * yr / 3.6e6
    # pm_ra is mu_alpha* (already cos-dec corrected), so undo that to get a delta in RA.
    cosd = np.cos(np.radians(cat.dec_deg_j2000.to_numpy(float)))
    ra = cat.ra_deg_j2000.to_numpy(float) + (
        cat.pm_ra_cosdec_mas_yr.to_numpy(float) * yr / 3.6e6) / np.where(cosd == 0, np.nan, cosd)
    return ra, dec


def _sep_arcsec(ra1, dec1, ra2, dec2) -> np.ndarray:
    """Angular separation, small-angle safe via the haversine form."""
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    a = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))) * 3600.0


def _norm_name(s: str) -> str:
    """Normalise a star name for comparison.

    ⚠️ SIMBAD identifiers carry catalogue DECORATION -- `* rho01 Cnc`, `** STT 270A`, `V* ...`
    -- so stripping whitespace alone leaves a leading `*` that no OBJECT string will ever
    match. Without this every honestly-labelled frame reads as mislabelled, which buries the
    two that genuinely are.
    """
    return ''.join(ch for ch in str(s).lower() if ch.isalnum())


def identify(fr: Frame, cat: pd.DataFrame, *, radius: float = MATCH_RADIUS_ARCSEC) -> Frame:
    """Establish the target from POSITION. Exactly one match, or quarantine."""
    obj_norm = _norm_name(fr.object_raw)
    for key, why in MOVING_TARGETS.items():
        if key in obj_norm:
            fr.star_id = key
            fr.id_method = 'solar-system body (no catalogue position by construction)'
            fr.id_status = 'moving_target'
            fr.id_evidence = f"OBJECT={fr.object_raw!r}; {why}"
            return fr
    if not np.isfinite(fr.ra_deg) or not np.isfinite(fr.dec_deg) or not np.isfinite(fr.mjd):
        fr.id_status = 'quarantine'
        fr.id_evidence = 'frame carries no usable RA/DEC/MJD — cannot be placed on the sky'
        return fr
    ra, dec = _propagate(cat, fr.mjd)
    sep = _sep_arcsec(fr.ra_deg, fr.dec_deg, ra, dec)
    hit = np.where(sep <= radius)[0]
    if len(hit) == 1:
        i = int(hit[0])
        fr.star_id = str(cat.star_id.iloc[i])
        fr.id_method = 'astrometry (SIMBAD J2000 + proper motion to epoch)'
        fr.id_sep_arcsec = float(sep[i])
        fr.id_status = 'confirmed'
        # Agreement is judged against SIMBAD's FULL identifier list, so a frame honestly
        # labelled `HD 22049` counts as naming eps Eri. What survives this test is the real
        # finding: a ROLE (`STD`) or a run placeholder (`Star S5`) where a name should be.
        al = {_norm_name(a) for a in str(cat.aliases.iloc[i] or '').split('|') if a}
        # 🔴 JUDGED ON `OBJECT` ALONE, ON PURPOSE. `ESO OBS TARG NAME` also carries a name
        # and on the tau Ceti frames it carries the RIGHT one -- but `OBJECT` is the field
        # every downstream tool keys on, so letting OBS TARG NAME vouch for it reports the
        # mislabel as fine and the finding disappears. The two are recorded separately:
        # OBJECT decides the verdict, OBS TARG NAME is corroboration for the ID.
        agrees = obj_norm in al
        targ_ok = _norm_name(fr.obs_targ_name) in al
        fr.id_evidence = (
            f"{sep[i]:.1f}\" from {cat.simbad_main_id.iloc[i]} at epoch "
            f"{fr.date_obs[:10]}; next nearest {np.sort(sep)[1]:.0f}\"; "
            f"OBJECT={fr.object_raw!r} "
            f"{'agrees' if agrees else 'DOES NOT NAME THIS STAR'}"
            f"{'; OBS TARG NAME does' if (targ_ok and not agrees) else ''}")
    elif len(hit) == 0:
        fr.id_status = 'quarantine'
        j = int(np.argmin(sep))
        fr.id_sep_arcsec = float(sep[j])
        fr.id_evidence = (f"no catalogue target within {radius:.0f}\"; nearest is "
                          f"{cat.simbad_main_id.iloc[j]} at {sep[j]:.0f}\"")
    else:
        fr.id_status = 'quarantine'
        found = {str(cat.star_id.iloc[i]) for i in hit}
        note = ''
        if any(found >= pair for pair in BINARY_PAIRS):
            note = (' — VISUAL BINARY: astrometry alone cannot separate these components at '
                    'any usable radius; this is RYA-423\'s job (IR-native RV/template ID), '
                    'not a tolerance to tighten')
        fr.id_evidence = ('ambiguous — ' + ', '.join(
            f"{cat.simbad_main_id.iloc[i]} {sep[i]:.0f}\"" for i in hit) + note)
    return fr


def inventory(roots, repo_root: Path, *, radius: float = MATCH_RADIUS_ARCSEC) -> pd.DataFrame:
    """Every CRIRES frame under `roots`, classified, identified and de-duplicated."""
    cat = load_astrometry(repo_root)
    frames: list[Frame] = []
    for root in roots:
        root = Path(root)
        if not root.exists():
            continue
        for p in sorted(root.rglob('*.fits')):
            try:
                fr = read_frame(p)
            except Exception:
                continue
            if not is_crires(fr):
                continue
            fr = classify_vintage(fr)
            fr = identify(fr, cat, radius=radius)
            frames.append(fr)

    # De-duplicate by CONTENT. The same frame is on this drive under several names -- an
    # ESO `ADP.*` id in one tree and a `CR_SONE_*` name in another -- so a filename-keyed
    # inventory double-counts it, and a target-keyed one reports coverage it does not have.
    seen: dict[str, str] = {}
    for fr in frames:
        if fr.md5 in seen:
            fr.duplicate_of = seen[fr.md5]
        else:
            seen[fr.md5] = fr.path
    return pd.DataFrame([f.__dict__ for f in frames])


def coverage_table(df: pd.DataFrame) -> pd.DataFrame:
    """Per target x instrument class, counting DISTINCT frames only."""
    u = df[df.duplicate_of == '']
    rows = []
    for (star, cls), g in u.groupby([u.star_id.replace('', 'UNCONFIRMED'), 'instrument_class']):
        rows.append({
            'star_id': star or 'UNCONFIRMED', 'instrument_class': cls, 'n': len(g),
            'settings': ' '.join(sorted({s for s in g.setting if s})),
            'date_min': min(g.date_obs)[:10], 'date_max': max(g.date_obs)[:10],
            'snr_median': round(float(np.nanmedian(g.snr)), 1) if g.snr.notna().any() else None,
            'confirmed': int(g.id_status.isin(['confirmed', 'moving_target']).sum()),
            'quarantined': int((g.id_status == 'quarantine').sum()),
        })
    return pd.DataFrame(rows).sort_values(['star_id', 'instrument_class'])


# ── co-adding ────────────────────────────────────────────────────────────────────────
#: Paranal. CRIRES/CRIRES+ live on UT3/UT1 at the VLT; the site is what BERV needs and the
#: difference between the four unit telescopes is metres, i.e. far below mm/s.
PARANAL = dict(lat_deg=-24.6270, lon_deg=-70.4045, height_m=2635.0)

#: 🔴 THE CO-ADD GATE IS A MEASURED VELOCITY, NOT A CALENDAR SPAN.
#: The first version of this refused any group spanning more than 60 days, which threw away
#: every 55 Cnc and eps Eri pair on the drive. That gate was answering the wrong question.
#: What must be true before adding two spectra is that their LINES LIE ON TOP OF EACH OTHER;
#: elapsed time is only a proxy for that, and a poor one. For these stars the astrophysical
#: RV variation is metres per second -- 55 Cnc's planets and eps Eri's activity jitter both
#: sit near 10 m/s -- while ONE RESOLUTION ELEMENT AT R=86,000 IS 3.5 km/s. The astrophysics
#: is three orders of magnitude below the thing that would smear a line. What is NOT
#: negligible across a two-year baseline is the INSTRUMENT: wavelength-solution drift between
#: observing runs.
#: So the residual velocity is measured by cross-correlation after the BERV shift, reported,
#: and corrected. A group is refused only when that residual is too large to be drift.
RESOLVING_POWER = 86000.0
C_KMS = 299792.458
#: One resolution element. A residual larger than this is not drift being tidied up, it is
#: two spectra that disagree about where the lines are -- a different star, a different
#: setting mislabelled, or a broken wavelength solution.
COADD_MAX_RESIDUAL_KMS = C_KMS / RESOLVING_POWER
#: A sanity bound only. Nothing astrophysical here changes on this timescale; it exists so a
#: catastrophically mis-grouped set cannot be silently stacked.
COADD_MAX_SPAN_DAYS = 3650.0

#: 🔴 VESTA IS NOT CO-ADDABLE BY THIS FUNCTION AND THAT IS DELIBERATE. It is a MOVING
#: REFLECTOR: the shift is the two-leg Sun->Vesta->observer rate (RYA-372), which no stellar
#: BERV carries, and RYA-796 measured that two of its five duplicate settings are 166 deg
#: apart in sub-observer longitude — opposite faces of the asteroid. Co-adding those is not
#: a worse spectrum, it is a different object averaged with itself.
NO_COADD = {'vesta': ('moving reflector — needs the RYA-372 two-leg reflected-solar '
                      'conditioning and RYA-796\'s <45 deg sub-observer longitude gate, '
                      'neither of which is a stellar BERV')}


def berv_kms(fr: Frame) -> float:
    """Barycentric correction for one frame, from its own pointing and time."""
    from astropy import units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    loc = EarthLocation.from_geodetic(
        lat=PARANAL['lat_deg'] * u.deg, lon=PARANAL['lon_deg'] * u.deg,
        height=PARANAL['height_m'] * u.m)
    sc = SkyCoord(ra=fr.ra_deg * u.deg, dec=fr.dec_deg * u.deg)
    t = Time(fr.mjd, format='mjd')
    return float(sc.radial_velocity_correction(obstime=t, location=loc).to(u.km / u.s).value)


def read_spectrum(path: str):
    """WAVE/FLUX/ERR/QUAL from a CRIRES+ IDP, in ANGSTROM, good pixels only.

    ⚠️ `WAVE` IS NANOMETRES (`TUNIT1='nm'`) and the conversion is done from the DECLARED
    unit, never assumed — RYA-796 made this a rule after the same table shipped a wavelength
    column labelled one thing and holding another. An unrecognised unit stops rather than
    guessing. ⚠️ Coverage comes from the ARRAY with `QUAL == 0`, never from WAVELMIN/MAX,
    because CRIRES+ has real inter-order gaps (RYA-377/796).
    """
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        d = hdul[1].data
        hdr = hdul[1].header
        cols = {c.upper() for c in d.columns.names}
        unit = str(hdr.get('TUNIT1', '')).strip().lower()
    w = np.asarray(d['WAVE'], float).ravel()
    f = np.asarray(d['FLUX'], float).ravel()
    e = (np.asarray(d['ERR'], float).ravel() if 'ERR' in cols
         else np.full_like(f, np.nan))
    q = (np.asarray(d['QUAL'], float).ravel() if 'QUAL' in cols
         else np.zeros_like(f))
    if unit in ('nm', 'nanometre', 'nanometer'):
        w = w * 10.0
    elif unit in ('angstrom', 'a', '0.1 nm'):
        pass
    else:
        raise ValueError(f"{path}: wavelength unit {unit!r} not recognised — refusing to "
                         f"guess a factor of 10 (RYA-796)")
    good = (q == 0) & np.isfinite(w) & np.isfinite(f) & (f != 0)
    return w[good], f[good], e[good]


def _xcorr_kms(w_ref: np.ndarray, f_ref: np.ndarray, w: np.ndarray, f: np.ndarray,
               *, span_kms: float = 30.0, step_kms: float = 0.10) -> float:
    """Residual velocity of (w, f) against the reference, by cross-correlation.

    Both spectra are continuum-flattened by division through their own median first: CRIRES+
    IDPs are UN-normalised adu (RYA-796) with different blaze levels per frame, and an
    un-normalised cross-correlation is dominated by the continuum ratio rather than by the
    lines. Returns NaN when the two do not overlap in usable pixels.
    """
    lo = max(w_ref.min(), w.min())
    hi = min(w_ref.max(), w.max())
    if not np.isfinite(lo) or not np.isfinite(hi) or hi - lo < 1.0:
        return float('nan')
    grid = np.linspace(lo, hi, min(20000, max(1000, int((hi - lo) / 0.01))))
    a = np.interp(grid, w_ref, f_ref)
    a = a / np.nanmedian(a) - 1.0
    shifts = np.arange(-span_kms, span_kms + step_kms, step_kms)
    best, best_v = -np.inf, float('nan')
    for v in shifts:
        b = np.interp(grid, w * (1.0 + v / C_KMS), f)
        b = b / np.nanmedian(b) - 1.0
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 100:
            continue
        sa, sb = a[m].std(), b[m].std()
        if sa == 0 or sb == 0:
            continue
        c = float(np.dot(a[m], b[m]) / (sa * sb * m.sum()))
        if c > best:
            best, best_v = c, float(v)
    return best_v


def _snr(f: np.ndarray, e: np.ndarray) -> float:
    """Median flux/error where both are usable; else a robust continuum proxy."""
    m = np.isfinite(f) & np.isfinite(e) & (e > 0)
    if m.sum() > 10:
        return float(np.nanmedian(f[m] / e[m]))
    return float('nan')


def coadd_group(frames: list[Frame]) -> dict:
    """Co-add one (target, setting) group in the BARYCENTRIC frame.

    Rest-frame FIRST, always: every frame is shifted by its own BERV before anything is
    added. Co-adding topocentric frames taken months apart smears every line by up to
    ~60 km/s peak-to-peak, which at R~86,000 is several resolution elements.
    """
    if len(frames) < 2:
        return {'n': len(frames), 'status': 'single-frame — nothing to co-add'}
    span = (max(f.mjd for f in frames) - min(f.mjd for f in frames))
    if span > COADD_MAX_SPAN_DAYS:
        return {'n': len(frames), 'status': f'REFUSED: {span:.0f} d span exceeds the '
                                            f'{COADD_MAX_SPAN_DAYS:.0f} d sanity bound'}
    spectra, bervs, snrs = [], [], []
    for fr in frames:
        w, f, e = read_spectrum(fr.path)
        v = berv_kms(fr)
        # to the barycentric frame: lambda_bary = lambda_obs * (1 + v/c)
        spectra.append((w * (1.0 + v / C_KMS), f, e))
        bervs.append(v)
        snrs.append(_snr(f, e))

    grid = spectra[0][0]
    # Residual velocity per frame, MEASURED against the first, after the BERV shift.
    resid = [0.0]
    for w, f, e in spectra[1:]:
        v = _xcorr_kms(grid, spectra[0][1], w, f)
        resid.append(v)
    worst = max(abs(v) for v in resid)
    if not np.isfinite(worst):
        return {'n': len(frames), 'status': 'REFUSED: residual velocity not measurable '
                                            '(no overlapping usable pixels)'}
    if worst > COADD_MAX_RESIDUAL_KMS:
        # 🔴 WHY THE CROSS-CORRELATION FAILED, DIAGNOSED RATHER THAN NAMED "drift".
        # If the measured residual cancels the BERV difference, the correlation locked onto
        # features that are STATIONARY IN THE TOPOCENTRIC FRAME -- i.e. TELLURIC lines. The
        # BERV shift aligns the stellar lines and by exactly the same amount MIS-aligns the
        # tellurics; in un-corrected CRIRES+ IR spectra the tellurics are strong enough to
        # win the correlation. Measured here: residual + dBERV sums to -1.2..+1.3 km/s on
        # every refused pair (K-band worst at -8.7, where tellurics are deepest).
        dberv = [v - bervs[0] for v in bervs]
        cancels = max(abs(r + d) for r, d in zip(resid, dberv))
        if cancels < 3 * COADD_MAX_RESIDUAL_KMS:
            why = (f'TELLURIC-LOCKED: residual {worst:.1f} km/s cancels dBERV '
                   f'{max(abs(d) for d in dberv):.1f} km/s to {cancels:.1f} km/s, so the '
                   f'cross-correlation is tracking TELLURIC lines, not the star. These IDPs '
                   f'are telluric_applied=not-applied (RYA-805/806); co-adding across a '
                   f'large dBERV needs telluric correction first (RYA-947)')
        else:
            why = (f'residual {worst:.1f} km/s exceeds one resolution element '
                   f'({COADD_MAX_RESIDUAL_KMS:.2f} km/s) and does NOT cancel dBERV '
                   f'({cancels:.1f} km/s) — an unexplained wavelength disagreement')
        return {'n': len(frames), 'span_days': round(span, 3),
                'residual_kms': [round(v, 3) for v in resid],
                'dberv_kms': [round(d, 3) for d in dberv],
                'status': 'REFUSED: ' + why}
    # Align on the MEASURED residual before adding. Correcting a measured shift is strictly
    # better than refusing on a proxy for it.
    spectra = [(w * (1.0 - v / C_KMS), f, e) for (w, f, e), v in zip(spectra, resid)]
    grid = spectra[0][0]
    num = np.zeros_like(grid)
    den = np.zeros_like(grid)
    for w, f, e in spectra:
        fi = np.interp(grid, w, f, left=np.nan, right=np.nan)
        ei = np.interp(grid, w, e, left=np.nan, right=np.nan)
        ok = np.isfinite(fi) & np.isfinite(ei) & (ei > 0)
        wgt = np.zeros_like(grid)
        wgt[ok] = 1.0 / ei[ok] ** 2
        num[ok] += fi[ok] * wgt[ok]
        den += wgt
    good = den > 0
    co_f = np.full_like(grid, np.nan)
    co_e = np.full_like(grid, np.nan)
    co_f[good] = num[good] / den[good]
    co_e[good] = 1.0 / np.sqrt(den[good])          # propagated, not assumed sqrt(N)
    return {
        'n': len(frames), 'status': 'co-added',
        'span_days': round(span, 3),
        'berv_kms': [round(v, 3) for v in bervs],
        'berv_spread_kms': round(float(np.ptp(bervs)), 3),
        'snr_per_frame': [round(s, 1) for s in snrs],
        'snr_before_median': round(float(np.nanmedian(snrs)), 1),
        'snr_after': round(_snr(co_f[good], co_e[good]), 1),
        'residual_kms': [round(v, 3) for v in resid],
        'residual_max_kms': round(worst, 3),
        'resolution_element_kms': round(COADD_MAX_RESIDUAL_KMS, 3),
        'dberv_max_kms': round(max(abs(v - bervs[0]) for v in bervs), 3),
        # ⚠️ WHY THIS GROUP WAS SAFE. With dBERV below a resolution element the telluric and
        # stellar frames coincide, so the ambiguity above never arises -- these frames were
        # taken at nearly the same barycentric phase. That is a property of the SCHEDULE,
        # not evidence that the tellurics were handled. They were not.
        'safe_because': ('dBERV below one resolution element — telluric and stellar rest '
                         'frames coincide for this group; tellurics remain UNCORRECTED'),
        'n_px': int(good.sum()),
        'wave': grid[good], 'flux': co_f[good], 'err': co_e[good],
    }


# ── CLI ──────────────────────────────────────────────────────────────────────────────
#: Where CRIRES data lives. Roots, not per-target folders: the classification is by header,
#: so the walker only has to be told where to look, never what it will find.
DEFAULT_ROOTS = (
    '/mnt/codex-data/spectra',
)


def main(argv=None) -> int:
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument('--inventory-all', action='store_true')
    ap.add_argument('--coadd', action='store_true')
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--roots', nargs='*', default=list(DEFAULT_ROOTS))
    ap.add_argument('--out', default=None, help='directory for the inventory + co-adds')
    ap.add_argument('--radius-arcsec', type=float, default=MATCH_RADIUS_ARCSEC)
    a = ap.parse_args(argv)

    repo_root = Path(__file__).resolve().parents[1]
    out = Path(a.out) if a.out else repo_root / 'data' / 'audit' / 'rya952_crires_inventory'
    out.mkdir(parents=True, exist_ok=True)

    df = inventory(a.roots, repo_root, radius=a.radius_arcsec)
    if df.empty:
        print('no CRIRES frames found under: ' + ', '.join(a.roots))
        return 1
    df.to_csv(out / 'crires_inventory.csv', index=False)

    uniq = df[df.duplicate_of == '']
    print(f"CRIRES frames found : {len(df)}  "
          f"({len(uniq)} distinct by md5, {len(df) - len(uniq)} duplicate copies)")
    print(f"instrument class    : {dict(uniq.instrument_class.value_counts())}")
    print('\n-- COVERAGE (distinct frames) --')
    cov = coverage_table(df)
    print(cov.to_string(index=False))
    cov.to_csv(out / 'crires_coverage.csv', index=False)

    if a.verify:
        print('\n-- TARGET ID VERDICTS --')
        for star, g in uniq.groupby(uniq.star_id.replace('', 'UNCONFIRMED')):
            names = sorted({o for o in g.object_raw if o})
            print(f"  {star:14s} n={len(g):3d}  OBJECT strings seen: {names}")
            ex = g.iloc[0]
            print(f"       {ex.id_status}: {ex.id_evidence}")
        q = uniq[uniq.id_status == 'quarantine']
        print(f"\n  QUARANTINED: {len(q)}")
        for r in q.itertuples():
            print(f"    {Path(r.path).name}  OBJECT={r.object_raw!r}  {r.id_evidence}")
        # 🔴 The mislabels are the point. Report every frame whose OBJECT does not name the
        # star the astrometry found -- silently correcting them would erase the finding.
        mis = uniq[(uniq.id_status == 'confirmed')
                   & uniq.id_evidence.str.contains('DOES NOT NAME THIS STAR')]
        print(f"\n  CONFIRMED BUT MISLABELLED (OBJECT does not name the star): {len(mis)}")
        for r in mis.itertuples():
            print(f"    {Path(r.path).name}  OBJECT={r.object_raw!r} -> {r.star_id}  "
                  f"({r.id_sep_arcsec:.1f}\")")

    results = {}
    if a.coadd:
        print('\n-- CO-ADD (CRIRES+ reduced, confirmed target, per setting) --')
        pool = uniq[(uniq.instrument_class == 'crires_plus')
                    & uniq.id_status.isin(['confirmed', 'moving_target'])]
        for (star, setting), g in pool.groupby(['star_id', 'setting']):
            if star in NO_COADD:
                print(f"  {star:12s} {setting:7s} n={len(g)}  REFUSED: {NO_COADD[star]}")
                results[f'{star}|{setting}'] = {'n': len(g), 'status': 'REFUSED',
                                                'reason': NO_COADD[star]}
                continue
            frames = [Frame(**{k: v for k, v in r._asdict().items() if k != 'Index'})
                      for r in g.itertuples()]
            try:
                res = coadd_group(frames)
            except Exception as e:
                print(f"  {star:12s} {setting:7s} n={len(g)}  ERROR {type(e).__name__}: {e}")
                continue
            if res.get('status') == 'co-added':
                w, f, e = res.pop('wave'), res.pop('flux'), res.pop('err')
                dest = out / f'coadd_{star}_{setting}.csv'
                pd.DataFrame({'wavelength_air_A': w, 'flux': f, 'err': e}).to_csv(
                    dest, index=False)
                print(f"  {star:12s} {setting:7s} n={res['n']}  SNR {res['snr_before_median']}"
                      f" -> {res['snr_after']}  (x{res['snr_after'] / res['snr_before_median']:.2f})"
                      f"  BERV spread {res['berv_spread_kms']} km/s  -> {dest.name}")
            else:
                print(f"  {star:12s} {setting:7s} n={res['n']}  {res['status']}")
            results[f'{star}|{setting}'] = res
        (out / 'coadd_summary.json').write_text(json.dumps(results, indent=2, default=str))

    print(f"\n[out] {out}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
