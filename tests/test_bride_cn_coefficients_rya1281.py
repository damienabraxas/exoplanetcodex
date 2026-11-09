import pytest
from scripts.probe_bride_cn_coefficients import reporting_patch,parse_states,ANCHOR,FIELDS


def test_reporting_patch_preserves_source_outside_single_insertion():
    source='prefix\n'+ANCHOR+'\nsuffix\n';rows=[dict(C=12,N=14,wave=7771.944)]
    patched,block=reporting_patch(source,rows)
    assert patched.replace(block,'')==source
    assert 'BRIDE_CN_STATE'in block
    assert all('&'not in line[6:]for line in block.splitlines()if not line.startswith('*'))
    assert all(len(x)<=72 for x in block.splitlines()if not x.startswith('*'))
    for bad in [source+ANCHOR,'no anchor']:
        with pytest.raises(ValueError):reporting_patch(bad,rows)


def test_state_parser_requires_complete_finite_numeric_rows():
    line='BRIDE_CN_STATE 12 14 1 1 '+' '.join(['1.0']*len(FIELDS))
    assert parse_states('other log\n'+line)[0]['depth_index']==1
    for bad in [line+' 1',line.replace('1.0','NaN',1)]:
        with pytest.raises(ValueError):parse_states(bad)
