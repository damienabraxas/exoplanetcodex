import pytest
from scripts.build_bride_native_lte import source_list


def test_recipe_preserves_hyphen_and_module_dependency_order():
    names=['clock-gf.f','model.f','strings_module.f90']+[f'file{i}.f'for i in range(30)]
    text='bsyn_luobj = '+' '.join('../$(SOURCE)/'+n for n in names)+'\n#clean :'
    result=source_list(text,'bsyn_lu',None)
    assert result[0]=='strings_module.f90'
    assert 'clock-gf.f'in result and 'clock'not in result
    assert set(result)==set(names) and len(result)==len(names)


def test_truncated_or_duplicate_recipe_refuses():
    with pytest.raises(ValueError):source_list('bsyn_luobj = ../$(SOURCE)/clock-gf.f\n#clean','bsyn_lu',None)
    names=['same.f']*31
    with pytest.raises(ValueError):source_list('bsyn_luobj = '+' '.join('../$(SOURCE)/'+n for n in names)+'\n#clean','bsyn_lu',None)
