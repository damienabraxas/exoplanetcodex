"""RYA-1229 — the per-line product is a projection OF THE FEED, and says so out loud.

🔴 WHAT WENT WRONG, SO THAT THESE TESTS SAY WHY THEY EXIST. `generate_perline_product.py`
discovered its inputs by globbing two directories named in a module constant:

    DEFAULT_BAND_PRODUCTS = [ROOT/"data/results/rya847/gated", ROOT/"data/results/rya877"]

The measurements moved to `data/results/band_products/` and the constant did not. Nothing
failed. The generator read 13 of 178 per-line files, so the artifact the Sun page links to
covered ONE of four instruments, carried zero rows for the engine the published headline
rests on, and labelled 932 of its 1039 rows `1D-LTE (ts-lte)` — the parent-directory name
from the pre-RYA-906 layout, a vocabulary no published product has ever used. It stayed
that way for a month.

⚠️ THE OBVIOUS FIX IS NOT A FIX, and this was measured rather than assumed: repointing the
constant at `data/results/band_products/` refuses with "6505 rows share a
(line x instrument x engine) key", because that directory holds every tier, selector and
route of every holding and the emitted key could not tell them apart.

So the tests below hold two lines. The first is that the stale constant cannot come back.
The second, and the one with teeth, is `test_the_committed_artifact_covers_the_published
_set`: the artifact's identity set must equal what the feed publishes and can be reached.
That test fails the day a new product is published, which is the alarm the old design had
no way to sound.
"""
from __future__ import annotations

import csv
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import perline_sources as ps            # noqa: E402
from pipeline.perline_product import (                # noqa: E402
    ROW_COLUMNS, PerLineProductError, _assert_inputs_committed,
)

PRODUCT = ROOT / "data" / "products" / "solar" / "Fe_perline.csv"
FEED = ROOT / "data" / "products" / "solar" / "Fe.json"
CLI = ROOT / "scripts" / "generate_perline_product.py"


def _rows(path: Path):
    body = "".join(l for l in path.read_text().splitlines(keepends=True)
                   if not l.startswith("#"))
    return list(csv.DictReader(io.StringIO(body)))


def _header(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if ": " in line:
            k, v = line[1:].split(": ", 1)
            out[k.strip()] = v.strip()
    return out


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def resolved():
    return ps.resolve_published("solar", "Fe")


# --------------------------------------------------------------------------------------
# 1. the constant that caused it cannot return
# --------------------------------------------------------------------------------------

def test_the_snapshot_default_is_gone_from_the_cli():
    """🔴 The defect was a PATH CONSTANT, so the regression test is about the constant.

    A future edit that reintroduces `DEFAULT_BAND_PRODUCTS` reintroduces a default that
    can silently outlive the directories it names.
    """
    src = CLI.read_text()
    assert "DEFAULT_BAND_PRODUCTS" not in src.split('"""')[-1], (
        "DEFAULT_BAND_PRODUCTS is back in executable code — discovery must come from the "
        "feed, not from a directory list a constant remembers")
    assert 'add_argument("--band-products"' not in src, (
        "--band-products is back: an argument that selects the input directory is an "
        "argument that can select a superseded snapshot")


def test_the_retired_roots_are_named_nowhere_as_an_input():
    for stale in ("rya847/gated", "results/rya877"):
        assert stale not in CLI.read_text().split('"""')[-1], (
            f"{stale} is referenced in executable CLI code again")


# --------------------------------------------------------------------------------------
# 2. resolution: the feed points at the evidence, and the pointer is gated
# --------------------------------------------------------------------------------------

def test_candidate_paths_strip_a_trailing_treatment_before_appending_one():
    """⚠️ `provenance.copied_to` names the treatment inconsistently across families.

    Appending without stripping produced `..._ENGINE-A_ENGINE-A_lines.csv` and lost 90 of
    160 products — a silent 56% coverage loss that looked like missing artifacts.
    """
    with_t = {"treatment": "ENGINE-A", "provenance": {
        "copied_to": "data/results/band_products/FeII_1_2_h_hold_SYNTH_T_ENGINE-A_products.csv"}}
    without = {"treatment": "ENGINE-A", "provenance": {
        "copied_to": "data/results/band_products/FeII_1_2_h_hold_SYNTH_T_products.csv"}}
    a = [p.name for p in ps.candidate_paths(with_t)]
    b = [p.name for p in ps.candidate_paths(without)]
    assert a == b, f"the two conventions must resolve to the same name: {a} != {b}"
    assert a[0] == "FeII_1_2_h_hold_SYNTH_T_ENGINE-A_lines.csv"
    assert a[0].count("ENGINE-A") == 1


def test_a_product_with_no_provenance_is_unresolved_not_guessed():
    path, frame, _, why = ps.resolve({"treatment": "1D-LTE", "n_lines": 3})
    assert path is None and frame is None
    assert "provenance.copied_to" in why


def test_resolution_is_gated_on_n_lines_value_equality(feed):
    """🔴 A PATH THAT EXISTS IS NOT A MATCH. RYA-1112 found two live products differing
    only in `selector` and publishing 0.196 dex apart, so 'the file is there' can mean
    'the neighbour's file is there'."""
    p = dict(next(q for q in feed["products"]
                  if (q.get("provenance") or {}).get("copied_to")))
    p["n_lines"] = int(p["n_lines"]) + 1
    path, frame, _, why = ps.resolve(p)
    assert path is None and frame is None
    assert "n_lines MISMATCH" in why


def test_every_3dnlte_product_reproduces_its_published_n_lines(feed):
    """The 3D-NLTE families keep every line and blank `a_3dnlte` for the ones the network
    refused, so membership is a DERIVED conjunction. It is only legitimate to derive it
    because it reproduces the published count on all twelve products, exactly."""
    ps3 = [q for q in feed["products"] if q["treatment"] == "ENGINE-A-3DNLTE"]
    assert len(ps3) == 12, f"expected 12 ENGINE-A-3DNLTE products, found {len(ps3)}"
    for q in ps3:
        path, frame, family, why = ps.resolve(q)
        assert path is not None, f"{q['tier']}/{q['instrument']}: {why}"
        assert family == "3dnlte"
        kept = int(frame["in_aggregate"].astype(str).eq("True").sum())
        assert kept == int(q["n_lines"]), (
            f"{q['tier']}/{q['instrument']}: derived {kept}, published {q['n_lines']}")


def test_the_3dnlte_adapter_publishes_the_3d_value_not_the_1d_base(feed):
    """🔴 The same file carries `a_1dlte` — the base the correction was applied to.
    Emitting it under the 3D-NLTE treatment would report different physics under this
    product's name."""
    import pandas as pd
    q = next(x for x in feed["products"] if x["treatment"] == "ENGINE-A-3DNLTE"
             and x["tier"] == "REFERENCE")
    path, frame, _, _ = ps.resolve(q)
    raw = pd.read_csv(path)
    keep = frame["in_aggregate"].astype(str).eq("True")
    got = frame.loc[keep, "abundance"].astype(float).reset_index(drop=True)
    want3d = raw.loc[keep.values, "a_3dnlte"].astype(float).reset_index(drop=True)
    want1d = raw.loc[keep.values, "a_1dlte"].astype(float).reset_index(drop=True)
    assert got.equals(want3d)
    assert not got.equals(want1d), "the adapter is publishing the 1D-LTE base"


def test_the_ew_column_of_the_other_3dnlte_schema_is_not_silently_blanked():
    """⚠️ RYA-1106 writes `ew_mA_agss21`, RYA-1095/1213 write `ew_mA`. Reading one name
    only would blank the measured EW on a third of the 3D-NLTE products while looking
    entirely healthy."""
    import pandas as pd
    feed = json.loads(FEED.read_text(encoding="utf-8"))
    q = next(x for x in feed["products"] if x["treatment"] == "ENGINE-A-3DNLTE"
             and x["selector"] == "ASPLUND_AGSS21")
    path, frame, _, _ = ps.resolve(q)
    raw = pd.read_csv(path)
    assert "ew_mA_agss21" in raw.columns and "ew_mA" not in raw.columns
    assert frame["ew_mA"].notna().any(), "ew_mA came through all-blank for this family"


# --------------------------------------------------------------------------------------
# 3. the emitted schema: every row says which published number it is evidence for
# --------------------------------------------------------------------------------------

def test_the_row_schema_carries_the_full_published_identity():
    for f in ps.IDENTITY_FIELDS:
        assert f in ROW_COLUMNS, f"{f} is not an emitted column"


def test_an_untracked_input_is_refused_and_not_merely_unmodified(tmp_path):
    """🔴 "NOT DIRTY" IS NOT "TRACKED". band_products is gitignored, and
    `git status --porcelain` says nothing about an ignored file — so the one check standing
    between this product and a commit_sha that describes nothing was blind to the exact
    case that directory's gitignore guarantees."""
    probe = ps.BAND_PRODUCTS / "_rya1229_untracked_probe.csv"
    probe.parent.mkdir(parents=True, exist_ok=True)
    probe.write_text("element,ion\nFe,I\n")
    try:
        assert not subprocess.run(
            ["git", "-C", str(ROOT), "status", "--porcelain", "--", str(probe)],
            capture_output=True, text=True).stdout.strip(), (
            "precondition: the probe must be invisible to `git status`, which is what "
            "made the old check useless")
        with pytest.raises(PerLineProductError) as e:
            _assert_inputs_committed([probe])
        assert "not tracked" in str(e.value)
    finally:
        probe.unlink(missing_ok=True)


# --------------------------------------------------------------------------------------
# 4. the alarm: the committed artifact must not drift from the feed
# --------------------------------------------------------------------------------------

@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_no_emitted_engine_label_is_outside_the_feeds_vocabulary(feed):
    """932 of 1039 rows of the old artifact failed this. A deck-suffixed label such as
    `1D-LTE (ts-lte)` cannot be joined to any published product."""
    vocab = {q["treatment"] for q in feed["products"]}
    bad = sorted({r["engine"] for r in _rows(PRODUCT)} - vocab)
    assert not bad, f"engine labels no published product uses: {bad}"


@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_engine_is_exactly_the_treatment():
    off = [r for r in _rows(PRODUCT) if r["engine"] != r["treatment"]]
    assert not off, f"{len(off)} rows where engine != treatment, e.g. {off[0]}"


@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_the_committed_artifact_covers_the_published_set(resolved):
    """🔴 THE ANTI-STALENESS TRIPWIRE. The identity set the artifact carries must be
    exactly the set of published products whose evidence is reachable. This fails the day
    a product is published that the artifact does not carry — which is the alarm the
    directory-glob design had no way to sound."""
    sources, _ = resolved
    want = {tuple(str(s.identity[k]) for k in ps.IDENTITY_FIELDS) for s in sources}
    have = {tuple(str(r[k]) for k in ps.IDENTITY_FIELDS) for r in _rows(PRODUCT)}
    missing = sorted(want - have)
    extra = sorted(have - want)
    assert not missing, (
        f"{len(missing)} published products have reachable per-line evidence that the "
        f"committed artifact does not carry — regenerate it. First: {missing[0]}")
    assert not extra, (
        f"{len(extra)} identities in the artifact are not published products with "
        f"reachable evidence. First: {extra[0]}")


@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_each_products_kept_rows_equal_its_published_n_lines(feed, resolved):
    """The projection must agree with the number it is evidence for, per product."""
    sources, _ = resolved
    by_id = {}
    for r in _rows(PRODUCT):
        if r["status"] == "in_aggregate":
            by_id[tuple(str(r[k]) for k in ps.IDENTITY_FIELDS)] = \
                by_id.get(tuple(str(r[k]) for k in ps.IDENTITY_FIELDS), 0) + 1
    bad = []
    for s in sources:
        k = tuple(str(s.identity[f]) for f in ps.IDENTITY_FIELDS)
        if k in by_id and by_id[k] != int(s.product["n_lines"]):
            bad.append((k, by_id[k], int(s.product["n_lines"])))
    assert not bad, f"kept-count != published n_lines for {len(bad)} products: {bad[:5]}"


@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_the_header_reports_reachability_rather_than_implying_completeness():
    h = _header(PRODUCT)
    for k in ("discovery", "products_published", "products_with_perline_evidence",
              "products_unresolved"):
        assert k in h, f"the header does not state {k}"
    assert "FEED-DRIVEN" in h["discovery"]
    assert (int(h["products_with_perline_evidence"]) + int(h["products_unresolved"])
            == int(h["products_published"]))


@pytest.mark.skipif(not PRODUCT.exists(), reason="product not committed")
def test_the_published_headline_engine_has_per_line_evidence():
    """The Fe I anchor the Sun page reports rests on ENGINE-A-3DNLTE. The old artifact had
    zero rows for it, so the headline number was the one number with no downloadable
    per-line evidence at all."""
    n = sum(1 for r in _rows(PRODUCT) if r["treatment"] == "ENGINE-A-3DNLTE")
    assert n > 0, "no ENGINE-A-3DNLTE rows: the headline still has no per-line evidence"
