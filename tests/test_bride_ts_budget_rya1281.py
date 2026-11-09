"""Native TS formats, quantization and the dormant >36-component print failure."""
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

from pipeline.bride_ts_budget import read_ts_opacity,read_ts_components,native_real32_print


def test_signed_geometry_and_native_opacity_units(tmp_path):
    path=tmp_path/'opacity'
    path.write_text(" 'MRXF' 2 5000\n2\n2700 2800\n"
        '-100 0.1 4000 2 3 4e-7 1 5\n1 2 3 4\n'
        '100 0.2 5000 3 4 5e-7 1 10\n5 6 7 8\n')
    result=read_ts_opacity(path)
    np.testing.assert_array_equal(result['absorption_m2_kg'],[[.5,5],[1.5,7]])
    np.testing.assert_array_equal(result['scattering_m2_kg'],[[1,6],[2,8]])
    np.testing.assert_allclose(result['electron_pressure_Pa'],[.2,.3])
    path.write_text(path.read_text()+' 12')
    with pytest.raises(ValueError):read_ts_opacity(path)


def component_fixture():
    labels=[f'a{i}'for i in range(40)]+['e-sc','H-sc','H2sc','He-sc']
    text=' BRIDE_FE_STATE 1 4000 2 4e-7 25 .01 .02\n'
    for wave in (5000,2770):
        text+=f'0WAVEL.={wave:7.0f}    ABS       SCAT  '+''.join(f'{x:>10}'for x in labels)+'\n'
        text+=' T= 4000.0 .1 .2 '+ ' '.join(['.025']*40+['.25']*4)+'\n'
    return text


def test_complete_ordered_components_and_real_precision(tmp_path):
    path=tmp_path/'detail';path.write_text(component_fixture())
    d=read_ts_components(path,1,1)
    assert len(d['labels'])==44
    assert d['records'].shape==(1,2,48)
    assert d['states'][0,3]==25
    np.testing.assert_allclose(native_real32_print([4045.9]),[4045.8999],rtol=0,atol=0)


def test_post_detabs_additions_are_separate_nonnegative_records(tmp_path):
    path=tmp_path/'detail';text=component_fixture()
    text+='BRIDE_POST_STATE 1 5000 .01 .02 .13 .2\nBRIDE_POST_STATE 1 2770 .01 .02 .13 .2\n'
    path.write_text(text);d=read_ts_components(path,1,1)
    assert d['post'].shape==(1,2,5)
    np.testing.assert_array_equal(d['post'][0,0],[5000.,.01,.02,.13,.2])
    path.write_text(text.replace('5000 .01','5000 -.01'))
    with pytest.raises(ValueError):read_ts_components(path,1,1)


@pytest.mark.parametrize('change',[
    lambda x:x.replace('BRIDE_FE_STATE','MISSING_STATE'),
    lambda x:x.replace('.025','-.025',1),
    lambda x:x.replace('.025','.125',1),
    lambda x:x[:x.rfind(' T=')],
])
def test_reject_incomplete_or_nonconserving_reporting(tmp_path,change):
    path=tmp_path/'detail';path.write_text(change(component_fixture()))
    with pytest.raises((ValueError,AssertionError)):read_ts_components(path,1,1)


def test_native_36_field_format_reversion_reproduces_failure(tmp_path):
    compiler=shutil.which('gfortran')
    if compiler is None:pytest.skip('Native compiler absent')
    for width,expect_success in ((36,False),(100,True)):
        program=tmp_path/f'format{width}.f90'
        program.write_text("program diagnostic\ninteger i\ncharacter(len=8) names(44)\n"
            "names='component'\nwrite(*,10) 2770.,names\n"
            f"10 format('0WAVEL.=',F7.0,' ABS SCAT ',{width}A10)\nend program\n")
        executable=tmp_path/f'format{width}'
        subprocess.run([compiler,str(program),'-o',str(executable)],check=True,capture_output=True,timeout=30)
        result=subprocess.run([str(executable)],capture_output=True,text=True,timeout=5)
        assert (result.returncode==0)==expect_success
        if not expect_success:assert 'Expected REAL' in result.stderr
