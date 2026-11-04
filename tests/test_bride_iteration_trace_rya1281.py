import json
from types import SimpleNamespace
import numpy as np
import pytest


def test_native_call_return_state_and_invalid_retry_are_not_replaced(tmp_path):
    from pipeline.bride_iteration_trace import StatEquilTrace
    from pipeline.bride_parameter_checks import positive_stat_equil
    history=tmp_path/'history.jsonl'; history.touch()
    active=SimpleNamespace(n=np.ones((2, 1)), atomicModel=SimpleNamespace(element='O'))
    class Context:
        activeAtoms=[active]
        def __init__(self): self.calls=0; self.clears=0; self.last=None
        def stat_equil(self):
            self.calls+=1; active.n[:]=(-1 if self.calls==1 else 2)
            self.last=SimpleNamespace(updatedPops=True, ngAccelerated=[self.calls==1])
            return self.last
        def clear_ng(self): self.clears+=1
    context=Context()
    with StatEquilTrace(context,tmp_path/'updates.jsonl',history) as trace:
        assert positive_stat_equil(trace,active,active.n.copy(),reject_invalid_ng=True)
        assert context.calls==2 and context.clears==1
        np.testing.assert_array_equal(active.n,[[2],[2]])
        with history.open('a') as f: f.write('{}\n')
        assert trace.stat_equil() is context.last
    rows=[json.loads(x) for x in (tmp_path/'updates.jsonl').read_text().splitlines()]
    assert [(r['outer_iteration'],r['attempt_in_outer']) for r in rows]==[(1,1),(1,2),(2,1)]
    assert [r['ngAccelerated'] for r in rows]==[[True],[False],[False]]
    assert [r['populations_positive'] for r in rows]==[[False],[True],[True]]


def test_missing_or_wrong_species_flags_refuse_instead_of_defaulting_false(tmp_path):
    from pipeline.bride_iteration_trace import StatEquilTrace
    history=tmp_path/'history.jsonl';history.touch()
    a=SimpleNamespace(n=np.ones((2,1)),atomicModel=SimpleNamespace(element='O'))
    ctx=SimpleNamespace(activeAtoms=[a],stat_equil=lambda:SimpleNamespace(updatedPops=True,ngAccelerated=[]))
    with StatEquilTrace(ctx,tmp_path/'updates.jsonl',history) as trace:
        with pytest.raises(ValueError,match='exact'):trace.stat_equil()


@pytest.mark.parametrize('ng',[False,True])
def test_recording_real_native_backend_preserves_analytic_fixture_arrays(monkeypatch,tmp_path,ng):
    import lightweaver as lw
    from pipeline.bride_iteration_trace import StatEquilTrace
    from pipeline.bride_numerical_validation import collision_lte_control
    options=dict(temperature=6000.,initial_weights=[100.,1.,.01],radiative=True,lambda_iterate=True,ng=ng)
    original=collision_lte_control(**options)
    history=tmp_path/'history.jsonl';history.touch()
    constructor=lw.Context;traces=[]
    def observed(*args,**kwargs):
        native=constructor(*args,**kwargs)
        trace=StatEquilTrace(native,tmp_path/'updates.jsonl',history);traces.append(trace)
        return trace
    monkeypatch.setattr(lw,'Context',observed)
    try: recorded=collision_lte_control(**options)
    finally:
        for trace in traces:trace.__exit__(None,None,None)
    for key in ('actual_fraction','conservation','expected_fraction'):
        np.testing.assert_array_equal(original[key],recorded[key])
    assert original['history']==recorded['history']
    rows=[json.loads(x) for x in (tmp_path/'updates.jsonl').read_text().splitlines()]
    assert rows and len(rows)==len(recorded['history'])
    assert all(len(r['ngAccelerated'])==1 and r['updatedPops'] for r in rows)
    assert any(any(r['ngAccelerated']) for r in rows)==ng
