"""A convergence-only comparison must refuse changes to physics or solver map."""
import copy,json,shutil,sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]


def test_ng_comparison_rejects_changed_background_map_and_budget(tmp_path,monkeypatch):
    from scripts.audit_bride_ng_baseline import main
    on=ROOT/'data/results/rya1281/solar_transition_v1/rya1281-solar-transition-v1/baseline'
    source=ROOT/'data/results/rya1281/solar_transition_no_ng_v1/rya1281-solar-transition-no-ng-v1/baseline'
    off=tmp_path/'off';shutil.copytree(source,off)
    original=json.loads((off/'report.json').read_text());output=tmp_path/'should_not_exist'
    for change in ['background','helper','method','budget']:
        bad=copy.deepcopy(original)
        if change=='background':bad['background']='Different opacity provider'
        elif change=='helper':bad['coupled']['helper']['sha256']='wrong-map'
        elif change=='method':bad['numerical_algorithm']['method']='different-map'
        else:bad['controls']['seconds_per_case']+=1
        (off/'report.json').write_text(json.dumps(bad))
        monkeypatch.setattr(sys,'argv',['audit','--on',str(on),'--off',str(off),'--output',str(output)])
        with pytest.raises(ValueError,match='Ng control|Expected Ng'):main()
        assert not output.exists()
