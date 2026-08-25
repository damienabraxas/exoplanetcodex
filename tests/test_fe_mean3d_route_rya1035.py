"""RYA-1035: the Fe ⟨3D⟩ route exists, and the vendor aux carries a defect that would
have selected the wrong star.

Step 0 of RYA-1035 asked whether Fe 3D-NLTE already exists as public-fetchable before
building anything. It does: `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` has been on the
same MPG Keeper share we already fetch every Gerber deck from since 2021. What Step 0 then
found is that consuming it is NOT a registry line — two things in the vendor's own files
break the RYA-821 reader, and one of them breaks it silently.

These tests run against the REAL aux tables (19 KB and 31 KB, committed under
`data/nlte_grids/gerber_ts/fe_mean3d_aux/`), not synthetic fixtures, because the whole
finding is a property of the vendor's actual bytes. The 92.9 MB deck binary stays
Sirius-only; the numbers quoted from it in the comments were measured against it directly.
"""
from pathlib import Path

import pytest

from pipeline import gerber_nlte as G

AUX_DIR = Path(__file__).resolve().parents[1] / "data" / "nlte_grids" / "gerber_ts" / "fe_mean3d_aux"
PLAIN = AUX_DIR / "auxData_Fe_STAGGERmean3D_May-21-2021.txt"
MARCS_NAMES = AUX_DIR / "auxData_Fe_STAGGERmean3D_May-21-2021_marcs_names.txt"

#: The solar star as the pipeline knows it (IAU 5772/4.438) against the deck's STAGGER
#: solar member (5777/4.44) — a convention difference for the same star, inside the
#: RYA-821 node tolerance.
SOLAR = dict(teff=5772.0, logg=4.438, feh=0.0)


@pytest.fixture(scope="module")
def plain_rows():
    return G._parse_aux_text(PLAIN.read_text())


@pytest.fixture(scope="module")
def marcs_rows():
    return G._parse_aux_text(MARCS_NAMES.read_text())


# ── the name parser: three conventions, and an `m` that means two different things ──

def test_the_short_teff_form_parses():
    """🔴 182 of the Fe ⟨3D⟩ deck's 189 rows are `p50g25m40`-shaped. A parser that knew
    only the 4-digit form returned None for all of them, and `read_deck_node` refuses a
    record it cannot identify — so the deck was unusable at every node but the solar one."""
    assert G._node_from_model_name("p50g25m40") == (5000.0, 2.5, -4.0)
    assert G._node_from_model_name("p60g45m10") == (6000.0, 4.5, -1.0)
    assert G._node_from_model_name("p50g50p05") == (5000.0, 5.0, 0.5)


def test_the_long_teff_form_still_parses():
    assert G._node_from_model_name("t5777g44m00") == (5777.0, 4.4, -0.0)
    assert G._node_from_model_name("p5777g44m10") == (5777.0, 4.4, -1.0)


def test_in_the_marcs_style_name_the_metallicity_is_z_and_m_is_the_MASS():
    """🔴 THE TRAP THAT ONLY SHOWS UP OFF-SOLAR. The two conventions put an `m` field in
    the same place and mean different things by it: STAGGER's `m00` is [Fe/H], MARCS's
    `m1.0` is the stellar mass, and MARCS's metallicity is `z-4.00` further along.

    Reading the MARCS `m` as metallicity is right at the solar node BY COINCIDENCE — mass
    0.0 and z+0.00 — and wrong at every other node. That coincidence is why the previous
    parser passed a test that only ever asked it about the Sun."""
    assert G._node_from_model_name(
        "s5000_g+2.5_m1.0_t02_st_z-4.00_a+0.00_c+0.00") == (5000.0, 2.5, -4.0)
    assert G._node_from_model_name(
        "p6000_g+4.5_m0.0_t02_st_z-1.00_a+0.00_c+0.00") == (6000.0, 4.5, -1.0)
    # the solar case, where the two readings agree and the bug is invisible
    assert G._node_from_model_name(
        "p5777_g+4.4_m0.0_t02_st_z+0.00_a+0.00") == (5777.0, 4.4, 0.0)


def test_the_giant_model_type_letter_is_accepted():
    """MARCS spherical models (logg < 3) are named `s...`, not `p...`. 24 of the Fe deck's
    rows are giants."""
    assert G._node_from_model_name("s5500_g+2.5_m1.0_t02_st_z-3.00_a+0.00") is not None


def test_an_unparseable_name_is_still_not_silently_accepted():
    assert G._node_from_model_name("garbage") is None
    assert G._node_from_model_name("") is None


# ── the vendor defect, on the vendor's real bytes ────────────────────────────

def test_the_plain_aux_zeroes_metallicity_on_exactly_the_seven_solar_teff_rows(plain_rows):
    """🔴 The deck's seven Teff=5777 members span the FULL metallicity axis — m00 / m05 /
    m10 / m20 / m30 / m40 / p05 — and the file's [Fe/H] column reads +0.00 for all seven."""
    assert len(plain_rows) == 189
    over = [r for r in plain_rows if r["feh_from_name"]]
    assert len(over) == 6, "6 of the 7 disagree; m00 is the one the column happens to get right"
    assert all(r["feh_aux"] == 0.0 for r in over)
    assert {r["id"] for r in over} == {
        "p5777g44m10", "p5777g44m20", "p5777g44m30", "p5777g44m40",
        "p5777g44p05", "p5777g44m05"}
    assert sorted(r["feh"] for r in over) == [-4.0, -3.0, -2.0, -1.0, -0.5, 0.5]


def test_the_other_182_rows_are_the_positive_control(plain_rows):
    """The column is wrong and the name is right, not the other way round — and that is
    MEASURED, not assumed. Every row outside the solar-Teff sequence agrees with its own
    name exactly, so the naming convention decodes correctly and the disagreement on the
    seven is a defect in the COLUMN."""
    non_solar = [r for r in plain_rows if r["teff"] != 5777.0]
    assert len(non_solar) == 182
    assert not any(r["feh_from_name"] for r in non_solar)


def test_the_referee_makes_the_solar_node_unambiguous(plain_rows):
    """🔴 WITHOUT IT THE SOLAR LOOKUP IS A 7-WAY TIE BROKEN BY FILE ORDER. All seven carry
    A(X)=7.50, so the A(X) tie-break is a no-op and the first row wins: `p5777g44m10`,
    [Fe/H] = −1.0. That is a well-formed departure block for a metal-poor star."""
    hits = [r for r in plain_rows
            if abs(r["teff"] - SOLAR["teff"]) <= G._NODE_TOL_TEFF
            and abs(r["logg"] - SOLAR["logg"]) <= G._NODE_TOL_LOGG
            and abs(r["feh"] - SOLAR["feh"]) <= G._NODE_TOL_FEH]
    assert [r["id"] for r in hits] == ["p5777g44m00"]

    uncorrected = [r for r in plain_rows
                   if abs(r["teff"] - SOLAR["teff"]) <= G._NODE_TOL_TEFF
                   and abs(r["logg"] - SOLAR["logg"]) <= G._NODE_TOL_LOGG
                   and abs(r["feh_aux"] - SOLAR["feh"]) <= G._NODE_TOL_FEH]
    assert len(uncorrected) == 7
    assert uncorrected[0]["id"] == "p5777g44m10", "file order picks the [Fe/H]=-1.0 deck"


def test_the_correction_is_recorded_never_silent(plain_rows):
    """A silently corrected input is the same class of defect as the one being corrected."""
    over = [r for r in plain_rows if r["feh_from_name"]]
    assert all(r["feh_aux"] != r["feh"] for r in over)
    clean = [r for r in plain_rows if not r["feh_from_name"]]
    assert all(r["feh_aux"] == r["feh"] for r in clean)


def test_the_ABUNDANCE_column_is_wrong_on_exactly_the_same_six_rows(plain_rows):
    """🔴 [Fe/H] IS NOT THE ONLY METALLICITY-DEPENDENT FIELD IN THE ROW. The deck's own
    relation A(X) = 7.50 + [Fe/H] is EXACT on all 183 clean rows and violated on exactly
    the six: they ship A(X) = 7.50 at [Fe/H] = -4.0 .. +0.5. Both fields were written as
    if the row were solar-metallicity."""
    clean = [r for r in plain_rows if not r["feh_from_name"]]
    assert {round(r["abundance"] - r["feh"], 4) for r in clean} == {7.5}, \
        "the relation that makes the six detectable"
    over = [r for r in plain_rows if r["feh_from_name"]]
    assert all(r["abundance"] == 7.5 for r in over)
    assert all(round(r["abundance"] - r["feh"], 4) != 7.5 for r in over)


def test_a_suspect_row_is_REFUSED_not_repaired():
    """🔴 THE NAME CAN REFEREE [Fe/H] BECAUSE IT ENCODES IT. NOTHING REFEREES A(X).

    So the override exists to get the bad rows OUT of the solar candidate set, not to
    repair them. Reconstructing A(X) from the deck's relation would be inventing vendor
    data, so a caller that lands on one of the six is refused."""
    rows = G._parse_aux_text(PLAIN.read_text())
    suspect = next(r for r in rows if r["id"] == "p5777g44m10")
    assert suspect["feh_from_name"] and suspect["feh"] == -1.0

    G._AUX_ROW_CACHE["FeTest@mean3D"] = rows
    G.DECKS["FeTest@mean3D"] = dict(Z=26, atom="atom.fe607a", aux="x", grid="y",
                                    read_via="direct")
    try:
        with pytest.raises(G.GerberDeckError, match="SUSPECT"):
            # ask for the node the SUSPECT row occupies -- Teff 5777, [Fe/H] -1.0
            G.read_deck_node("FeTest@mean3D", teff=5777.0, logg=4.44, feh=-1.0,
                             abundance=7.5)
    finally:
        G._AUX_ROW_CACHE.pop("FeTest@mean3D", None)
        G.DECKS.pop("FeTest@mean3D", None)


# ── why the canonical aux is the one we must NOT use ─────────────────────────

def test_the_marcs_names_aux_is_UNRECOVERABLE(marcs_rows):
    """🔴 THE VENDOR'S CONVERTER DESTROYED THE ONLY SURVIVING EVIDENCE, and it is the file
    TSFitPy's own downloader points at ([Fe] 3d_aux_link).

    `convert_3d_grid_to_marcs_names.py` builds the new name FROM the [Fe/H] column, so the
    zeroing propagated into the name: all seven Teff=5777 rows come out as one
    byte-identical string. Name and column now agree — and both are wrong for six of them —
    so the referee finds nothing to correct and the 7-way tie returns."""
    assert len(marcs_rows) == 189
    assert not any(r["feh_from_name"] for r in marcs_rows), \
        "nothing left to referee with: the converter overwrote the name too"

    solar_names = {r["id"] for r in marcs_rows if r["teff"] == 5777.0}
    assert len(solar_names) == 1, "seven distinct models collapsed to ONE name"

    hits = [r for r in marcs_rows
            if abs(r["teff"] - SOLAR["teff"]) <= G._NODE_TOL_TEFF
            and abs(r["logg"] - SOLAR["logg"]) <= G._NODE_TOL_LOGG
            and abs(r["feh"] - SOLAR["feh"]) <= G._NODE_TOL_FEH]
    assert len(hits) == 7, "still a 7-way tie -- this aux cannot address the solar node"


def test_the_two_aux_files_describe_the_same_deck(plain_rows, marcs_rows):
    """Same 189 records, same pointers — they differ only in how the models are NAMED. So
    choosing the plain one costs nothing but the name convention."""
    assert ([r["pointer"] for r in plain_rows] == [r["pointer"] for r in marcs_rows])
    assert ([r["abundance"] for r in plain_rows] == [r["abundance"] for r in marcs_rows])


# ── what the deck is, so a later wiring ticket does not re-derive it ─────────

def test_fe_mean3d_has_no_abundance_axis(plain_rows):
    """🔴 Fe IS NOT SHAPED LIKE Al, AND HERE THAT MAKES IT SIMPLER. A(X) = 7.50 + [Fe/H]
    for every one of the 189 rows: one abundance per atmosphere node, no free axis. So the
    v118 abundance-axis machinery Al needed (31 nodes, departures differing by up to 10.28
    between adjacent ones) is not required for Fe ⟨3D⟩ — the pre-v118 hoisting, which was
    always correct for Fe, stays correct."""
    solar_abundances = {round(r["abundance"], 4) for r in plain_rows if abs(r["feh"]) < 1e-9}
    assert solar_abundances == {7.5}
    clean = [r for r in plain_rows if not r["feh_from_name"]]
    assert {round(r["abundance"] - r["feh"], 3) for r in clean} == {7.5}


def test_the_deck_covers_the_sun_and_a_wide_box(plain_rows):
    """Teff 4000-7000, logg 1.5-5.0, [Fe/H] -4..+0.5, with a DEDICATED solar member
    dropped into an otherwise 500 K Teff grid. Far wider than the vendored Amarsi+2022 Fe
    MLP (Teff 5000-6500, logg 4.0-4.5, [Fe/H] -3..0), and unlike that MLP it is a
    DEPARTURE deck — not tied to the Jofre golden-line list or to 4787-6810 A."""
    assert sorted({r["teff"] for r in plain_rows}) == [
        4000.0, 4500.0, 5000.0, 5500.0, 5777.0, 6000.0, 6500.0, 7000.0]
    assert min(r["logg"] for r in plain_rows) == 1.5
    assert max(r["logg"] for r in plain_rows) == 5.0
    assert sorted({r["feh"] for r in plain_rows}) == [-4.0, -3.0, -2.0, -1.0, -0.5, 0.0, 0.5]


def test_fe_mean3d_is_registered_against_the_PLAIN_aux():
    """🔴 THE MATRIX DERIVES 'WIRED' FROM `DECKS` (RYA-821/v117), so this entry may only
    exist while the bytes really are on Sirius — otherwise the status surface claims a
    capability we cannot run, which is the exact lie v117 removed. RYA-1035 staged and
    md5-verified the deck; RYA-710 added the registry line, which is its charter.

    🔴 AND IT MUST BE THE PLAIN AUX, THE OPPOSITE OF Al. The vendor's converter builds the
    MARCS-style name FROM the [Fe/H] column it zeroed, collapsing seven distinct
    atmospheres to one string and making the solar node unaddressable. Measured both ways:
    `_marcs_names` refuses the solar node, plain returns `p5777g44m00`. Pinning the FILE
    here is the point — this is not a stylistic choice between two equivalent inputs."""
    assert "Fe@mean3D" in G.DECKS
    cfg = G.DECKS["Fe@mean3D"]
    assert cfg["aux"] == "auxData_Fe_STAGGERmean3D_May-21-2021.txt"
    assert "_marcs_names" not in cfg["aux"], "the converted aux is unrecoverable for Fe"
    assert cfg["read_via"] == "direct"
    assert cfg["atom"] == "atom.fe607a"
    # Al's, by contrast, IS the _marcs_names one -- and correctly so (0/6345 defects)
    assert "_marcs_names" in G.DECKS["Al@mean3D"]["aux"]


def test_fe_mean3d_carries_its_own_3D_atmosphere():
    """⟨3D⟩ departures on a 1D MARCS structure produce a perfectly well-formed file, which
    is why this is a guard and not a comment. The deck is keyed at STAGGER's 5777/4.44, a
    node MARCS_SOLAR (5750/4.5) does not contain."""
    assert G.deck_atmosphere("Fe@mean3D") != G.MARCS_SOLAR
    assert "stagger" in G.deck_atmosphere("Fe@mean3D").lower()
    assert G.deck_atmosphere("Fe") == G.MARCS_SOLAR, "the 1D deck is untouched"


def test_registering_fe_mean3d_did_not_give_it_an_abundance_axis(plain_rows):
    """🔴 THE DIFFERENCE THAT IS INVISIBLE FROM THE REGISTRY LINE. Al's ⟨3D⟩ deck resolves
    31 abundances and needs the v118 per-trial interpolation; Fe's resolves ONE, so its
    departures are a property of the atmosphere node and may be hoisted out of the χ² loop
    — which is what the pre-v118 code always did, correctly, for Fe."""
    G._AUX_ROW_CACHE["Fe@mean3D"] = plain_rows
    G._AXIS_CACHE.pop("Fe@mean3D", None)
    try:
        assert G.abundance_axis("Fe@mean3D") == (7.5,)
        assert not G.has_abundance_axis("Fe@mean3D")
    finally:
        G._AUX_ROW_CACHE.pop("Fe@mean3D", None)
        G._AXIS_CACHE.pop("Fe@mean3D", None)


# ── (B) the deck abundance: resolved by provenance, and guarded against re-drift ──

def test_the_model_atom_declares_7_50_and_the_aux_agrees(plain_rows):
    """🔴 THE GRID'S OWN DECLARATION. `atom.fe607a` line 2 is `7.50  55.85` — A(Fe) and the
    atomic mass — and both Fe aux tables encode A(X) = 7.50 + [Fe/H] exactly. The atom file
    is Sirius-only, so what is pinned here is the aux half plus the md5 that ties the atom
    we read to the one staged (`d08dc8232ed68eec65f9bb6631e82ea8`, in Fe_gerber2023.prov.json).
    """
    clean = [r for r in plain_rows if not r["feh_from_name"]]
    assert {round(r["abundance"] - r["feh"], 4) for r in clean} == {7.5}
    solar = [r for r in clean if abs(r["feh"]) < 1e-9]
    assert {r["abundance"] for r in solar} == {7.5}


def test_the_provenance_record_says_7_50_not_7_46():
    """🔴 7.46 WAS OUR OWN INPUT COMING BACK. `abu_ref` is read from stdin
    (interpol_modeles_nlte.f:206), written verbatim into the departure file (:761), loaded
    as `abundance_nlte` (read_departure.f) and printed back by bsyn (:988) — and
    `gerber_nlte` fed that stdin from `deck_abundance()` itself. The record was its own
    echo. ⚠️ 7.46 is the Asplund solar A(Fe): the number that looks right on arrival."""
    import json
    p = Path(__file__).resolve().parents[1] / "data" / "nlte_grids" / "gerber_ts" / \
        "Fe_gerber2023.prov.json"
    d = json.loads(p.read_text())
    assert float(d["deck_abundance"]["a_sun"]) == 7.50
    assert d["deck_abundance"]["corrected_from"] == "7.46"
    assert "atom.fe607a" in d["deck_abundance"]["why"]


def test_the_abundance_we_pass_does_not_move_the_departures(plain_rows):
    """The measurement that makes `abu_ref` a LABEL rather than a physical input. Fe's deck
    has no abundance axis, so every request returns the same node — and the reported A(X)
    is the AUX's own value, never the caller's. This is why the ⟨3D⟩ direct-read path
    always said 7.50 while the interpolator path echoed 7.46: one reports the grid, the
    other reports us."""
    G.DECKS["FeAbTest@mean3D"] = dict(Z=26, atom="atom.fe607a", aux="x", grid="y",
                                      read_via="direct")
    G._AUX_ROW_CACHE["FeAbTest@mean3D"] = plain_rows
    try:
        assert not G.has_abundance_axis("FeAbTest@mean3D")
        assert G.abundance_axis("FeAbTest@mean3D") == (7.5,)
        # the node lookup is indifferent to the abundance asked for
        hits = {a: min([r for r in plain_rows if r["teff"] == 5777.0
                        and abs(r["feh"]) < 1e-9],
                       key=lambda r: abs(r["abundance"] - a))["id"]
                for a in (7.46, 7.50, 7.60, 6.90)}
        assert set(hits.values()) == {"p5777g44m00"}
    finally:
        G._AUX_ROW_CACHE.pop("FeAbTest@mean3D", None)
        G.DECKS.pop("FeAbTest@mean3D", None)
        G._AXIS_CACHE.pop("FeAbTest@mean3D", None)


def test_a_provenance_abundance_that_contradicts_the_deck_is_REFUSED(plain_rows, tmp_path):
    """🔴 THE GUARD THAT MAKES THE LOOP UNCLOSEABLE. Nothing ever compared the record
    against the deck, so a number that came from us was handed BACK to the vendor binary as
    though it had come from the vendor. Now the aux cross-examines the record."""
    import json
    G.DECKS["FeProv@mean3D"] = dict(Z=26, atom="atom.fe607a", aux="x", grid="y",
                                    read_via="direct")
    G._AUX_ROW_CACHE["FeProv@mean3D"] = plain_rows
    old_prov = G.PROV_DIR
    try:
        G.PROV_DIR = tmp_path
        (tmp_path / "FeProv_gerber2023.prov.json").write_text(
            json.dumps({"deck_abundance": {"a_sun": 7.46}}))
        with pytest.raises(G.GerberDeckError, match="the deck's own aux declares"):
            G.deck_abundance("FeProv@mean3D")
        # ...and the corrected value passes
        (tmp_path / "FeProv_gerber2023.prov.json").write_text(
            json.dumps({"deck_abundance": {"a_sun": 7.50}}))
        assert G.deck_abundance("FeProv@mean3D") == 7.50
    finally:
        G.PROV_DIR = old_prov
        G._AUX_ROW_CACHE.pop("FeProv@mean3D", None)
        G._AXIS_CACHE.pop("FeProv@mean3D", None)
        G.DECKS.pop("FeProv@mean3D", None)


def test_the_provenance_record_is_found_through_an_atmosphere_suffix():
    """`Fe@mean3D` must resolve to `Fe_gerber2023.prov.json`. Both Fe decks were solved
    with the same model atom at the same A(Fe), and without the split the deck would look
    unregistered the moment it is staged."""
    assert G.deck_abundance("Fe") == 7.50
