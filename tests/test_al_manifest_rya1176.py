import pytest
import sys
from pipeline.al_manifest import load_manifest, require_product_manifest

def test_manifest_has_explicit_axes_and_frozen_denominator():
    d = load_manifest()
    assert len(d) == 505 and d.canonical_line_id.is_unique
    assert d.line_set.value_counts().to_dict() == {"our-all": 493, "our-deep-graded": 7, "our-graded": 5}
    assert d.telluric_applied.eq("unknown").all()
    assert d.normalization_state.eq("unknown").all()
    assert d.observed_conditioning.eq("unknown").all()

def test_product_guard_requires_holding_and_conditioning():
    d = load_manifest()
    with pytest.raises(ValueError, match="no holding supplied"):
        require_product_manifest(d)
    with pytest.raises(ValueError, match="normalization_state=unknown"):
        require_product_manifest(d, holding={"telluric_applied": "not-applied", "normalization_state": "unknown", "observed_conditioning": "solar"})
    out = require_product_manifest(d, holding={"telluric_applied": "not-applied", "normalization_state": "normalised", "observed_conditioning": "solar_kpno_v1"})
    assert out.telluric_applied.eq("not-applied").all() and out.normalization_state.eq("normalised").all()

def test_publisher_refuses_Al_without_observed_conditioning(monkeypatch):
    from scripts import publish_product
    monkeypatch.setattr(sys, "argv", ["publish_product", "--from", "data/results/band_products/AlI_4200_6910_harps_solar_harps_SYNTH_products.csv", "--holding", "solar_harps", "--tier", "ALL", "--dry-run"])
    assert publish_product.main() == 10
