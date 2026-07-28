"""
tests/test_hfs_synthesis_rya411.py
==================================
RYA-411 — the HFS-resolved synthesis row-expansion (pipeline.pysme_nlte._linelist_rows).

The actual PySME synthesis needs the gitignored departure grid (manual/offline tool), so
here we unit-test only the pure row-building logic: an HFS feature must expand to one PySME
line row PER hyperfine component, each carrying the feature's SHARED lower/upper NLTE level
labels (so all components map to the same departure coefficients) but its OWN wavelength +
gflog. This is the fix for the RYA-410 Mn over-saturation (gf-summed single line).
"""
from pipeline.pysme_nlte import _linelist_rows


# a non-HFS feature (single line) and an HFS feature (3 components), shared labels.
_SINGLE = (6304.906, -2.0, 4.889, 4.0, 6.85, 4.0, 'e8S', '8P*', 0.0)
_HFS = (6016.67, -0.5, 3.073, 2.0, 5.13, 2.0, 'z6P*', 'e6S', 0.0,
        [(6016.586, -1.54), (6016.635, -0.86), (6016.682, -1.03)])


def test_single_line_emits_one_row():
    rows = _linelist_rows('Mn', [_SINGLE])
    assert len(rows) == 1
    assert rows[0]['wlcent'] == 6304.906
    assert rows[0]['gflog'] == -2.0


def test_hfs_feature_expands_to_one_row_per_component():
    rows = _linelist_rows('Mn', [_HFS])
    assert len(rows) == 3                                   # NOT collapsed to one gf-summed line
    assert [r['wlcent'] for r in rows] == [6016.586, 6016.635, 6016.682]
    assert [r['gflog'] for r in rows] == [-1.54, -0.86, -1.03]


def test_hfs_components_share_the_feature_nlte_levels():
    rows = _linelist_rows('Mn', [_HFS])
    # every component maps to the SAME lower/upper level + EP/J (HFS splits the level by
    # ~ueV -> identical departure coefficients); only wl + gf differ.
    for r in rows:
        assert r['term_lower'] == 'z6P*' and r['term_upper'] == 'e6S'
        assert r['excit'] == 3.073 and r['e_upp'] == 5.13
        assert r['species'] == 'Mn 1' and r['atom_number'] == 25


def test_mixed_single_and_hfs():
    rows = _linelist_rows('Mn', [_SINGLE, _HFS])
    assert len(rows) == 1 + 3
