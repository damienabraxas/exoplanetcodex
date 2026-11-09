import pytest
from scripts.analyze_bride_cn_binding_closure import read_binding_table


def test_binding_csv_preserves_boolean_without_spurious_nan(tmp_path):
    p=tmp_path/'state.csv';p.write_text('temperature_K,arithmetic_passed\n5000,True\n6000,True\n')
    data=read_binding_table(p)
    assert data.dtype['arithmetic_passed'].kind=='b'
    assert data['temperature_K'].tolist()==[5000,6000]


def test_invalid_numeric_or_false_arithmetic_flag_refused(tmp_path):
    p=tmp_path/'state.csv'
    for value,flag in [('NaN','True'),('5000','False'),('5000','maybe')]:
        p.write_text(f'temperature_K,arithmetic_passed\n{value},{flag}\n')
        with pytest.raises(ValueError):read_binding_table(p)
