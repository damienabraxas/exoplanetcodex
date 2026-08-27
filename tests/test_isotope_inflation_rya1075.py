"""
tests/test_isotope_inflation_rya1075.py
=======================================
RYA-1075 guards for the HFS-per-isotope INFLATION — the consumer-side half of the
isotope class that RYA-684 left open.

RYA-684 closed the ENGINE side (``isotopfrac`` applied on top of already-folded gf,
offset ``-log10(sum f_i^2)``) and was scoped to "sweep every HFS/isotope species". It
was closed with 54 live instances of the CONSUMER side still sitting in
``canonical_gf.csv``, and it left no guard — so nothing failed for two months, until
RYA-1070's linemake cross-reference happened to trip over them.

That is the failure this file exists to prevent, so the tests are built the way that
failure would have been caught:

  A. The ARITHMETIC must hold anywhere — synthetic form-(A) blocks with the answer
     built in, so these run on the Mac with no engines mounted. A guard whose only
     tests are grid-gated is a guard that is silently skipped where it is authored.
  B. The CORRECTION must stay applied — the 54 rows are pinned against the committed
     provenance sidecar, arithmetically, not by remembering the numbers.
  C. The DEFECT must be re-detectable from the SOURCE, and REINTRODUCING it must fail.
  D. The SELECTOR must be specific — Li I 6707 is a real multi-isotope cluster whose
     isotopes do NOT each carry the full gf. It must never be corrected.

⚠️ The offset here is ``log10(n_isotopes)`` — a COUNT — not RYA-684's
``-log10(sum f_i^2)``. La II settles it: La is 99.911% La-139, so RYA-684's term is
+0.0008 while a count gives +0.3010, and +0.3010 is what was measured.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_resolver import physical_total, ISOTOPE_FORM_A_SPREAD_DEX  # noqa: E402
from pipeline.isotope_gf_convention import isotope_inflated_rows            # noqa: E402

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
SIDECAR = ROOT / "data" / "linelists" / "canonical_gf_isotope_corrections_rya1075.csv"


def _ges_path():
    """The GES v6 synth list, or None where the engines are not mounted."""
    try:
        import pipeline.abundances_derive as ad
        p = Path(ad._SYNTH_LINELIST_FILE)
        return p if p.exists() else None
    except Exception:
        return None


needs_ges = pytest.mark.skipif(_ges_path() is None,
                               reason="GES v6 synth line list not mounted (engines absent)")


# ── A. the arithmetic, anywhere ──────────────────────────────────────────────
def test_physical_total_divides_out_the_isotope_count():
    """Two isotopes each carrying the full gf must total the gf, not twice it."""
    # one transition, gf = 1.0 (log 0.0), split into 3 HFS components per isotope
    comp = np.log10(np.array([0.5, 0.3, 0.2]))
    gf = np.concatenate([comp, comp])
    iso = np.array([151] * 3 + [153] * 3)
    assert physical_total(gf, iso) == pytest.approx(0.0, abs=1e-12)
    # the isotope-blind sum is the defect, and it is exactly log10(2) high
    naive = float(np.log10(np.sum(10.0 ** gf)))
    assert naive - physical_total(gf, iso) == pytest.approx(math.log10(2), abs=1e-12)


@pytest.mark.parametrize("n", [2, 5, 7])
def test_offset_is_log10_of_the_count_not_of_the_abundances(n):
    """The correction is a COUNT. It must not move when the abundances do.

    This is the RYA-684 confusion pinned: that defect's term depends on the isotope
    fractions, this one does not. Same components, wildly different natural abundances,
    identical answer.
    """
    comp = np.log10(np.array([0.6, 0.4]))
    gf = np.concatenate([comp] * n)
    iso = np.repeat(np.arange(100, 100 + n), 2)
    naive = float(np.log10(np.sum(10.0 ** gf)))
    assert naive - physical_total(gf, iso) == pytest.approx(math.log10(n), abs=1e-12)


def test_not_form_A_is_left_alone():
    """If the isotopes do NOT each carry the full gf, the count rule is undefined."""
    gf = np.array([-0.002, -0.303])          # Li I 6707's shape: one component each
    iso = np.array([7, 6])
    naive = float(np.log10(np.sum(10.0 ** gf)))
    assert physical_total(gf, iso) == pytest.approx(naive, abs=1e-12)


def test_single_isotope_and_uncoded_are_untouched():
    gf = np.array([-0.5, -0.7, -1.2])
    assert physical_total(gf, np.array([151, 151, 151])) == pytest.approx(
        float(np.log10(np.sum(10.0 ** gf))), abs=1e-12)
    assert physical_total(gf, np.zeros(3, int)) == pytest.approx(
        float(np.log10(np.sum(10.0 ** gf))), abs=1e-12)
    assert physical_total(gf, None) == pytest.approx(
        float(np.log10(np.sum(10.0 ** gf))), abs=1e-12)


def test_form_A_spread_threshold_separates_the_real_populations():
    """The threshold is nowhere near a boundary: measured form-(A) clusters sit at
    <=0.0097 dex and the one real exception at 0.301."""
    assert 0.0097 < ISOTOPE_FORM_A_SPREAD_DEX < 0.301


# ── B. the correction stays applied ──────────────────────────────────────────
def test_sidecar_and_table_agree_and_the_arithmetic_holds():
    assert SIDECAR.exists(), "the RYA-1075 provenance sidecar must ship with the correction"
    side = pd.read_csv(SIDECAR)
    canon = pd.read_csv(CANONICAL, low_memory=False).set_index("line_id")
    assert len(side) == 54, f"expected 54 corrected rows, sidecar carries {len(side)}"
    for r in side.itertuples():
        # the correction term IS -log10(n), derived, not transcribed
        assert r.correction_term == pytest.approx(-math.log10(r.n_isotopes), abs=1e-6), \
            f"{r.line_id}: correction term is not -log10({r.n_isotopes})"
        assert r.corrected_log_gf == pytest.approx(
            r.published_log_gf + r.correction_term, abs=6e-4), \
            f"{r.line_id}: corrected != published - log10(n)"
        assert float(canon.at[r.line_id, "log_gf"]) == pytest.approx(
            r.corrected_log_gf, abs=1e-9), f"{r.line_id}: table no longer carries it"
        assert canon.at[r.line_id, "adjudication_status"] == "isotope_rya1075"


def test_every_corrected_row_matches_its_vald_sibling_where_one_exists():
    """The independent corroboration. One row is a known, reported exception."""
    side = pd.read_csv(SIDECAR)
    have = side[side.sibling_gf_linelist_vald.notna()]
    assert len(have) >= 45, "lost the sibling corroboration"
    off = (have.corrected_log_gf - have.sibling_gf_linelist_vald).abs()
    # Cu I 5782.122 disagrees by 0.116 dex — GES vs VALD, reported in RYA-1075, and
    # RYA-684 already established Cu I's VALD surface is folded and so is not a
    # trustworthy referee for this species.
    assert (off > 0.01).sum() == 1, "the set of sibling disagreements changed"
    assert have.loc[off.idxmax(), "line_id"] == "gf_048042"
    # Everything else agrees to better than 0.01 dex; the largest of those is Ba II
    # 4934 at 0.0072, a plain GES-vs-VALD catalogue difference.
    assert off.drop(off.idxmax()).max() < 0.01


# ── C. re-detection from the source, and reintroduction ──────────────────────
@needs_ges
def test_corrected_table_is_clean():
    canon = pd.read_csv(CANONICAL, low_memory=False)
    assert isotope_inflated_rows(canon, _ges_path()) == []


@needs_ges
def test_reintroducing_one_row_fails_loud():
    """The guard RYA-684 lacked. A single re-inflated row must be named."""
    canon = pd.read_csv(CANONICAL, low_memory=False)
    i = canon.index[canon.line_id == "gf_051798"][0]      # Eu II 6645.09
    canon.at[i, "log_gf"] = 0.4208                        # the isotope-blind sum
    hits = isotope_inflated_rows(canon, _ges_path())
    assert len(hits) == 1
    assert hits[0]["line_id"] == "gf_051798"
    assert hits[0]["n_isotopes"] == 2
    assert hits[0]["inflation_dex"] == pytest.approx(math.log10(2), abs=1e-3)


@needs_ges
def test_reintroducing_all_54_fails_loud():
    canon = pd.read_csv(CANONICAL, low_memory=False).set_index("line_id")
    for r in pd.read_csv(SIDECAR).itertuples():
        canon.at[r.line_id, "log_gf"] = r.published_log_gf
    hits = isotope_inflated_rows(canon.reset_index(), _ges_path())
    assert len(hits) == 54


# ── D. specificity ───────────────────────────────────────────────────────────
@needs_ges
def test_li_6707_is_never_corrected():
    """Positive control. Li I 6707 spans two isotopes whose sets do NOT agree
    (spread 0.301), so it is a genuine multi-isotope cluster that this rule must
    leave alone. If the selector ever widens enough to swallow it, that is a blanket
    correction and the whole fix is unsafe."""
    canon = pd.read_csv(CANONICAL, low_memory=False).set_index("line_id")
    for r in pd.read_csv(SIDECAR).itertuples():
        canon.at[r.line_id, "log_gf"] = r.published_log_gf
    hits = isotope_inflated_rows(canon.reset_index(), _ges_path())
    assert "gf_087098" not in {h["line_id"] for h in hits}
    assert "Li I" not in set(pd.read_csv(SIDECAR).species)


def test_no_uncorrected_species_leaked_into_the_sidecar():
    """La II 5971 is NOT this defect — it is a duplicated row in the GES source
    (two byte-identical entries, iso=0), which also sums to +log10(2). Correcting it
    under an isotope rationale would attach the wrong provenance to the right number."""
    side = pd.read_csv(SIDECAR)
    assert set(side.species) == {"Ba II", "Cu I", "Eu II", "Nd II", "Sm II"}
    assert "La II" not in set(side.species)
