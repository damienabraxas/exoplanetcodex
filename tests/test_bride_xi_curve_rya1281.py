import copy,json
from pathlib import Path
import numpy as np
import pytest


def reports():
    root=Path(__file__).resolve().parents[1]/'data/results/rya1281/xi_slice_v1/runs'
    return [json.loads(p.read_text())for p in sorted(root.glob('*/report.json'))]


def test_comparison_accepts_declared_xi_variation_but_rejects_other_changes():
    from scripts.compare_bride_xi_curves import check_shared_inputs
    r=reports();check_shared_inputs(r)
    bad=copy.deepcopy(r);bad[1]['pressure_source']['sha256']='altered'
    with pytest.raises(ValueError,match='physics'):check_shared_inputs(bad)
    bad=copy.deepcopy(r);fields=bad[1]['local_binding']['local_binding']['fields'];fields['temperature_K'][0]=np.nextafter(fields['temperature_K'][0],np.inf)
    with pytest.raises(ValueError,match='fixed atmosphere'):check_shared_inputs(bad)
    bad=copy.deepcopy(r);bad[1]['coupled']['controls']['rtol']=1e-4
    with pytest.raises(ValueError,match='coupled control'):check_shared_inputs(bad)
    bad=copy.deepcopy(r);bad[1]['population_mode']='fixed'
    with pytest.raises(ValueError,match='physics mode'):check_shared_inputs(bad)
