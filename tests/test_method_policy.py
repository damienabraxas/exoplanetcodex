"""
tests/test_method_policy.py
===========================
Tests for pipeline/method_policy.py — the per (star × species) method-selection
resolver (RYA-315 / RYA-306).

Mirrors the five smoke cases from the RYA-315 spec, rewritten to the adopted
species-key convention "<Element> <RomanIon>" (e.g. "Fe I"). The canonical
element set is sourced from constants.TARGET_ELEMENTS (the 27-element list,
RYA-109) rather than hardcoded — the resolver stays star-agnostic and the test
exercises the real element set.
"""
import pytest

from config.constants import TARGET_ELEMENTS
from pipeline.method_policy import get_method

# 26 unique atomic symbols (Fe covers Fe I and Fe II via the ion key).
ELEMS = set(TARGET_ELEMENTS)


class TestGetMethod:
    def test_fe_ii_override_synthesis(self):
        """Fe II -> synthesis, evidence-confirmed (seeded override, RYA-305)."""
        cell = get_method('sun', 'Fe II', canonical_elements=ELEMS)
        assert cell['method'] == 'synthesis'
        assert cell['confidence'] == 'evidence-confirmed'

    def test_fe_i_override_ew(self):
        """Fe I -> ew, evidence-confirmed (seeded override, RYA-305)."""
        cell = get_method('sun', 'Fe I', canonical_elements=ELEMS)
        assert cell['method'] == 'ew'
        assert cell['confidence'] == 'evidence-confirmed'

    def test_unoverridden_atomic_is_ew_default_not_a_gap(self):
        """Ti I -> ew, default-pending-data. MUST NOT raise.

        The key behavior of the fail-loud narrowing: an un-overridden *known*
        atomic species is the EW default, not an unclassifiable gap.
        """
        cell = get_method('procyon', 'Ti I', canonical_elements=ELEMS)
        assert cell['method'] == 'ew'
        assert cell['confidence'] == 'default-pending-data'

    def test_molecular_is_synthesis_a_priori(self):
        """CN -> synthesis, a-priori (no EW on a band)."""
        cell = get_method('sun', 'CN', canonical_elements=ELEMS)
        assert cell['method'] == 'synthesis'
        assert cell['confidence'] == 'a-priori'

    def test_unclassifiable_species_fails_loud(self):
        """Xy I -> ValueError (not a molecule, not a known element)."""
        with pytest.raises(ValueError):
            get_method('sun', 'Xy I', canonical_elements=ELEMS)
