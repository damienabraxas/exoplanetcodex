import sys
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rya1354_mcdhf_s0_controls as C  # noqa: E402


def test_variants_and_method3_lists():
    assert C.VARIANTS==('default','method3_f','method3_all')
    assert C.method3_list('method3_f','even')=='4f*,5f*' and C.method3_list('method3_f','odd')==''
    assert C.method3_list('method3_all','odd')=='*'


def test_default_stdin_is_v53_with_more_cycles():
    assert C.rmcdhf_stdin('default','even','1\n2\n')=='y\n1\n2\n5\n*\n*\n500\n'


def test_non_default_stdin_follows_getscd_prompt_order():
    s=C.rmcdhf_stdin('method3_f','even','1\n').split('\n')
    # default?, debug?, light/grid?, ACCY?, levels, weights, varied, spectroscopic, NSCF
    assert s[:9]==['n','n','n','n','1','5','*','*','500']
    # modify? y; threshold? n; methods? y; method 1,2 none; method 3 f; method 4 none
    assert s[9:16]==['y','n','y','','','4f*,5f*','']
    # first-oscillation, subshell accel, eigenvector accel, NSIC, NSOLV, Schmidt: defaults; SCF order 1
    assert s[16:22]==['n']*6 and s[22:]==['1','']
    assert C.rmcdhf_stdin('method3_f','odd','1\n') is None


def test_prompt_check(tmp_path):
    good=('Modify other defaults?  (y/n) Select a different integration method for any\n Method            3 :\n'
          ' Subshell radial wavefunctions will be Schmidt orthogonalised immediately\n')
    (tmp_path/'rmcdhf.stderr').write_text(good)
    assert C.prompt_check(tmp_path,'method3_all')=='checked'
    (tmp_path/'rmcdhf.stderr').write_text(good+' GETRSL: Unable to interpret string\n redo ...\n')
    with pytest.raises(ValueError,match='mismatch'):C.prompt_check(tmp_path,'method3_all')
    (tmp_path/'rmcdhf.stderr').write_text('Default settings?')
    assert C.prompt_check(tmp_path,'default')=='checked'
    with pytest.raises(ValueError,match='not reached'):C.prompt_check(tmp_path,'method3_f')


def test_method3_all_matches_the_proven_ladder_sequence_without_damping():
    """rya1354_lifetime_ladder (reference_method3, n=0, full grid) ran this prompt path natively; v55
    differs only in declining subshell acceleration (ladder: y, n, -0.9) and NSIC kept default."""
    roots='1,2\n1\n'
    ladder=('n\nn\nn\nn\n'+roots+'5\n'+'*\n*\n500\n'+'y\nn\n'+'y\n\n\n*\n\n'+'n\ny\nn\n-0.9\nn\n'+'n\n'+'n\nn\n1\n')
    ours=C.rmcdhf_stdin('method3_all','even',roots)
    assert ladder.replace('n\ny\nn\n-0.9\nn\n','n\nn\nn\n')==ours


def test_v54_archive_records_the_missing_debug_prompt():
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_controls_failed_v54',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_controls_v54/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and 'even_method3_f/rmcdhf: rc=2' in r['error']
    assert r['runs']['even_default']['outcome']=='NOT_CONVERGED'
    err=files['runs/mcdhf_s0_controls_v54/even_method3_f/rmcdhf.stderr'].decode()
    assert 'Generate debug output?' in err and 'Bad real number' in err


def test_v55_archive_improv_stop_is_classified_as_measured_non_convergence():
    """v55 FAILED_RETAINED: method 3 solved 4f/5f, then 6p failed methods 1 and 2 -> native IMPROV
    ERROR STOP. v56 classifies exactly that as a measured variant outcome; other refusals still raise."""
    import json
    from rya1354_native_parts import load_parts
    files,_=load_parts(ROOT/'data/audit/rya1354/mcdhf_s0_controls_failed_v55',repo_root=ROOT)
    r=json.loads(files['runs/mcdhf_s0_controls_v55/receipt.json'])
    assert r['status']=='FAILED_RETAINED' and r['error'].endswith("rc=1 ['ERROR STOP']")
    out=files['runs/mcdhf_s0_controls_v55/even_method3_f/rmcdhf.stdout'].decode()
    msg=r['error'].split(': ',1)[1]
    assert C.classify_stop(msg,out)=='improv'
    assert 'could not be solved using method 2' in out and '6p' in out.rsplit('Failure;',1)[1][:40]
    assert C.classify_stop("Native refusal in x/rmcdhf: rc=0 ['MAXIMUM ITERATIONS IN SCF EXCEEDED']",'')=='iteration_limit'
    assert C.classify_stop("Native refusal in x/rmcdhf: rc=1 ['ERROR STOP']",'some other stop')is None
    assert C.classify_stop("Native refusal in x/rmcdhf: rc=2 []",out) is None
    # v55 consumed every non-default prompt in order
    err=files['runs/mcdhf_s0_controls_v55/even_method3_f/rmcdhf.stderr'].decode()
    assert 'Orthonomalization order?' in err and C.prompt_check.__name__
