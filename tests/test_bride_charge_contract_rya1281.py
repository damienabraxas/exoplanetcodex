"""Physical limits and fine-structure accounting for the sourced MN experiment."""
import numpy as np
import pytest


def test_sourced_pair_full_atom_balance_and_complete_branching():
    pytest.importorskip('lightweaver')
    from config.constants import require_codex_path
    from pipeline.bride_reference import read_marcs_structure,make_reference_context
    from pipeline.multi_atom import read_multi_atom
    from pipeline.bride_charge_contract import TERM_MAP,read_author_rates,install_author_mn_pair
    from pipeline.bride_native_spline import NativeTautsp,install_native_collision_control
    from pipeline.bride_rate_audit import process_flow
    import json
    refpath=require_codex_path('experiments.bride_oxygen_atom').parent/'reference-5750-4.5-A8.70-v1.json'
    if not refpath.exists():pytest.skip('Held benchmark unavailable')
    ref=json.loads(refpath.read_text())
    source=read_multi_atom(require_codex_path('experiments.bride_oxygen_atom'))
    structure=read_marcs_structure(require_codex_path('grids.marcs_standard')/ref['atmosphere_file']['name'])
    ctx=make_reference_context(structure,source,rays=3,qscale=1)
    oracle=NativeTautsp(require_codex_path('experiments.bride_multi_spline_oracle'))
    install_native_collision_control(ctx,source,oracle,source_directions=True)
    directory=require_codex_path('experiments.bride_barklem_oxygen')
    receipt=install_author_mn_pair(ctx,source,directory,oracle)
    assert receipt['CH0_replaced']==24
    atom=ctx.activeAtoms[0]
    rates={rate.i+1:rate for rate in atom.atomicModel.collisions if rate.source_kind=='CH0'}
    grid,tables,_=read_author_rates(directory)
    # The sum of fine-structure MN branches must equal the independently tabulated term rate.
    for levels,state,total in TERM_MAP:
        expected=oracle.interpolate(grid,tables[:,state-1,18],structure['temperature_K'],logarithmic=True)*1e-6
        np.testing.assert_allclose(sum(rates[j].coefficient for j in levels),expected,rtol=1e-14)
    atom.compute_collisions(fillDiagonal=False)
    net,gross=process_flow(np.array(atom.C),np.array(atom.nStar))
    np.testing.assert_array_equal(net[gross==0],0.)
    assert np.max(abs(net[gross>0])/gross[gross>0])<1024*np.finfo(float).eps
    # Conditional inverse requires the declared reservoir, even after installation.
    h=ctx.kwargs['eqPops'].atomicPops['H'];saved=h.n.copy();h.n[:]*=1.01
    try:
        with pytest.raises(ValueError,match='requires LTE hydrogen'):
            atom.compute_collisions(fillDiagonal=False)
    finally:h.n[:]=saved
