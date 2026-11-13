"""An inherited operator cannot silently fall back to the backend default."""
import pytest
from pipeline.bride_coupled_contract import validate_operator_contract

def test_operator_receipt_must_match_inherited_report():
    receipt=dict(name='receipt.json',bytes=123,sha256='source-qualified-receipt')
    source=dict(iteration_plugin=dict(receipt=receipt))
    correct=dict(resource='experiments.bride_full_cap_plugin',build_receipt=receipt)
    assert validate_operator_contract(dict(iteration_plugin=correct),source)==correct
    for plan in [{},dict(iteration_plugin={}),dict(iteration_plugin=dict(correct,build_receipt=dict(receipt,sha256='other-build'))),dict(iteration_plugin=dict(correct,resource='default'))]:
        with pytest.raises(ValueError,match='rate-preserving operator'):validate_operator_contract(plan,source)
