"""Native rate formula, writer accounting and refusal disposition."""
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_rate_audit import summarize,review,constants
from rya1354_native_parts import load_parts


def test_all_four_native_cases_have_complete_model_accounting():
    report=summarize(ROOT/'data/audit/rya1354')
    assert report['accepted_lifetime'] is None and report['scientific_admission']=='HOLD'
    for row in report['cases'].values():
        assert row['printed_rows_verified']==10 and row['raw_nonzero_channels']==11
        assert len(row['explicit_model_angular_zeros'])==2
        assert 2<row['exploratory_model_lifetime_ns']<3
        assert 0<row['writer_suppressed_nonzero'][0]['A_s1']<1e-18


@pytest.mark.parametrize('mutation',['table','refusal','basis','occupation','missing'])
def test_corruption_and_unexplained_omission_block_totals(mutation):
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_rates_complete_v21')
    root='runs/fac_model_rates_v21/';p=root+'coherent100_B/'
    raw=files[p+'raw.tsv'].decode();table=files[p+'rates.tr'].decode();basis=files[p+'o_export_000.c'].decode()
    c=constants(files[root+'source_documents/faclib/consts.h'].decode())
    if mutation=='table':table=table.replace('2.302055E+08','2.302056E+08')
    if mutation=='refusal':raw=raw.replace('-3 2 0 -1 0 0','-3 2 0 -2 0 0')
    if mutation=='basis':basis=basis.replace('3s','3d')
    if mutation=='occupation':basis=basis.replace('2p ( 2)','2p ( 1)',1)
    if mutation=='missing':raw='\n'.join(raw.splitlines()[:-1])
    with pytest.raises(ValueError):review(raw,table,basis,c,2)
