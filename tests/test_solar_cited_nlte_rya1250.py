import json
from types import SimpleNamespace
import pytest
from pipeline import solar_nlte_tables as tab
from pipeline.run_matrix import preparation_hold_for_cell


def test_cited_mgii_differences_use_marcs_lte_baseline():
    waves=[7877.054,7896.366,9218.25,9244.265,10914.244]
    assert list(tab.deltas('Mg','II',waves,'1D',star='solar').values()) == [-.016,-.029,-.078,-.058,-.063]
    assert list(tab.deltas('Mg','II',waves,'3D',star='solar').values()) == [-.002,-.034,-.028,-.016,-.024]
    assert tab.deltas('Mg','I',waves,'3D',star='solar') == {}
    assert tab.deltas('Mg','II',waves,'3D',star='procyon') == {}
    assert tab.deltas('Mg','II',[7877.5],'3D',star='solar') == {}


def report():
    return dict(ready=False,scoped_readiness=True,nlte_missing=[],adopt=[],cull_candidates=[],
                missing_lines=[dict(species='Mg II',wavelength_A=10092.2)],holdings=[],
                literature={'litscan':True,'missing_bib_keys':[]},
                blocking={'missing_solar_literature_lines':1})


def test_ir_or_other_ion_hold_does_not_block_prepared_optical_cell():
    d=SimpleNamespace(element='Mg',ion='I',lo_A=4000,hi_A=6900,holding='h')
    assert preparation_hold_for_cell(report(),d) == ''
    d.ion='II';d.lo_A=9000;d.hi_A=13000
    assert 'missing_lines' in preparation_hold_for_cell(report(),d)
    d.ion='I'
    assert preparation_hold_for_cell(report(),d) == ''


def test_pending_adoption_and_unknown_reports_still_block():
    d=SimpleNamespace(element='Mg',ion='II',lo_A=7000,hi_A=9500,holding='h')
    r=report();r['adopt']=[{'species':'Mg II','wavelength_A':7877.}]
    assert 'adopt' in preparation_hold_for_cell(r,d)
    assert preparation_hold_for_cell({'ready':False},d)
    r=report();r['blocking']['unrecognized']=1
    assert 'unknown' in preparation_hold_for_cell(r,d)


def test_adding_unknown_sigma_cannot_mislabel_a_line_as_graded():
    from pipeline.published_line_add import plan,LineAddError
    with pytest.raises(LineAddError,match='uncertainty'):
        plan({'authors_sigma_dex':None})


def test_registered_atom_reach_refuses_ionized_nlte_but_not_lte():
    from pipeline.gerber_nlte import nlte_ion_capability
    from pipeline.run_matrix import _deck_hold
    assert nlte_ion_capability('Mg','II')[0] is False
    assert nlte_ion_capability('Mg','I')[0] is True
    d=SimpleNamespace(element='Mg',ion='II',engine_deck='gerber-nlte')
    assert 'no bound-bound' in _deck_hold(d)


def test_published_recipe_can_request_mean_without_changing_default_median():
    from pipeline.band_products import LineMeasurement,build_product
    from pipeline.line_sets import aggregation_for
    lines=[LineMeasurement('Mg','II',7800+i,'test',10,'synthesis',abundance=a,
                           treatment='1D-LTE') for i,a in enumerate([7.,7.,10.])]
    args=('Mg','II','test','red-optical','1D-LTE',lines)
    assert build_product(*args,handler='SynthesisHandler').value == 7.
    p=build_product(*args,handler='SynthesisHandler',aggregation=aggregation_for('MG_AGSS21=x.csv'))
    assert p.value == 8.
    assert 'unweighted mean' in p.provenance
    assert aggregation_for(None) == 'median'
    with pytest.raises(ValueError,match='Unsupported'):
        build_product(*args,handler='SynthesisHandler',aggregation='bogus')


def test_adoption_preserves_unrelated_high_precision_fields(tmp_path,monkeypatch):
    from pipeline import element_prepare as prep,line_sets
    path=tmp_path/'canonical.csv'
    path.write_text('line_id,species,wavelength_air_A,log_gf,gf_tier,loggf_reference\n'
                    'target,Mg II,7877.054,0.39,OTHER,old\n'
                    'untouched,Fe I,15879.557896921768,-1.9638878507483917,OTHER,old\n')
    monkeypatch.setattr(prep,'CANONICAL',path)
    monkeypatch.setattr(prep,'ROOT',tmp_path)
    monkeypatch.setattr(line_sets,'load_sources',lambda:[{'name':'test','citation':'ref'}])
    monkeypatch.setattr(line_sets,'build',lambda *a,**k:None)
    prep.apply('solar','Mg',dict(adopt=[dict(line_id='target',species='Mg II',
               canonical_wavelength_A=7877.054,authors_loggf=.389,authors_sigma_dex=.0086,
               set='test')],missing_lines=[],cull_candidates=[]))
    import csv
    row=list(csv.DictReader(path.open()))[1]
    assert row['wavelength_air_A']=='15879.557896921768'
    assert row['log_gf']=='-1.9638878507483917'


def test_grid_loader_separates_ions_and_keeps_declared_legacy_tables(tmp_path,monkeypatch):
    import pandas as pd
    from pipeline import nlte_corrections as nc
    from config.constants import NLTE_CORRECTION_ELEMENTS
    rows=[]
    for ion,wave,delta in [(1,5000.,.1),(2,6000.,.9)]:
        for teff,logg,feh in [(5700.,4.,0.),(5800.,4.,0.),(5700.,4.5,0.),(5700.,4.,.2)]:
            rows.append(dict(element='Mg',ion=ion,wave_A=wave,teff_K=teff,logg=logg,
                             feh=feh,delta_nlte=delta))
    path=tmp_path/NLTE_CORRECTION_ELEMENTS['Mg']['grid']
    df=pd.DataFrame(rows);df.to_csv(path,index=False)
    monkeypatch.setattr(nc,'_MPIA_ELEMENT_DIR',tmp_path)
    monkeypatch.setattr(nc,'_mpia_element_cache',{})
    assert list(nc._load_mpia_element_grid('Mg')['waves']) == [5000.]
    df[df.ion==1].drop(columns='ion').to_csv(path,index=False)
    nc._mpia_element_cache.clear()
    assert list(nc._load_mpia_element_grid('Mg')['waves']) == [5000.]


def test_provenance_names_only_sources_that_supplied_the_pool():
    kw=dict(star='solar',native_source='native grid')
    assert tab.applied_source('Mg','I',{5711.088:-.03},**kw) == 'native grid'
    cited=tab.applied_source('Mg','II',{7877.054:-.016},**kw)
    assert 'CITED solar' in cited and 'native grid' not in cited
    missing=tab.applied_source('Mg','II',{7877.054:-.016,9000:float('nan')},**kw)
    assert missing == cited
    mixed=tab.applied_source('Mg','I',{6318.717:0.,5711.088:-.03},**kw)
    assert 'native grid' in mixed and 'CITED solar' in mixed


def test_declared_but_ungraded_reference_selection_remains_visible():
    from pipeline import run_matrix as rm
    cells=rm.expand('solar','Mg',ions=['I'],engines=['ts-lte'],methods=['synthesis'],pools=['set'])
    vis=[d for d in cells if d.lo_A < 6318 < d.hi_A]
    assert vis, 'The eight declared Mg I features must not disappear from coverage'
    # Derivation is allowed from the declared literature selection.  The missing
    # per-line gf uncertainty remains a final uncertainty/publication hold.
    assert all(rm.cell_process_hold(d) == '' for d in vis)


def test_incomplete_asplund_ion_set_derives_from_all_declared_lines():
    from pipeline.run_descriptor import line_set, set_derivation_csv
    row = line_set("MG_AGSS21", "Mg", "II")
    assert set_derivation_csv(row, "solar").endswith(
        "data/reference/line_sets/mg_agss21_MgII.csv")


def test_budget_pairs_use_the_declared_statistic_on_paired_lines():
    import pandas as pd
    from scripts.rya1230_cno_budget import _paired
    nom=pd.DataFrame(dict(wavelength_air_A=[5000.,5100.,5200.],abundance=[0.,0.,0.],
                         in_aggregate=[True]*3))
    leg=nom.copy();leg['abundance']=[1.,2.,9.]
    assert _paired(leg,nom,3)['effect'] == 2.
    mean=_paired(leg,nom,3,aggregation='mean')
    assert mean['effect'] == 4.
    assert mean['n_paired'] == 3 and not mean['moved']
