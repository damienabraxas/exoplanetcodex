"""Joint outer-layer controls and actual qualified native n5 initializer."""
from pathlib import Path
from types import SimpleNamespace
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_joint_outer import controls,qualified_source
from rya1354_frozen_extend import controls as n5_controls


def test_joint_scope_changes_only_declared_varied_shells():
    for state in ('ground','quintet','triplet'):
        for verify in (False,True):
            text=n5_controls(state,verify)
            assert controls(state,verify)==text.replace('\n5*\n','\n4*,5*\n')
            assert '2*,3*' not in controls(state,verify)


def test_real_n5_native_source_passes_and_changed_mixing_is_refused(tmp_path):
    base=ROOT/'data/audit/rya1354'
    args=SimpleNamespace(source=tmp_path/'n5',n4_source=tmp_path/'n4',parent=tmp_path/'parent',
                         n4_core=tmp_path/'core4',core=tmp_path/'core5')
    for bundle,prefix,dest in [('n5_complete_v15','runs/lifetime_frozen_prior_n5_v1/',args.source),
            ('frozen_complete_v12','runs/lifetime_frozen_prior_n4_v1/',args.n4_source),
            ('from_n3_terminal_v12','runs/lifetime_relax_n4_from_relaxed_n3_v7/',args.parent),
            ('lifetime_aux_v4b','runs/lifetime_core_n4_v1/',args.n4_core),
            ('lifetime_core_n5_v4','runs/lifetime_core_n5_v1/',args.core)]:
        files,_=load_parts(base/bundle)
        for name,value in files.items():
            if name.startswith(prefix):
                p=dest/name.removeprefix(prefix);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(value)
    result=qualified_source(args)
    assert len(result['states'])==6 and len(result['decay']['channels'])==13
    assert not result['adoption_allowed'] and result['accepted_lifetime'] is None
    p=args.source/'outer_n5_triplet_verify/state.m';p.write_bytes(p.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='source artifact differs'):
        qualified_source(args)
