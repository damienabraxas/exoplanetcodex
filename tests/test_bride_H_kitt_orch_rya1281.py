import numpy as np
import pytest
from scripts.analyze_bride_H_kitt_orch import validate_processed_flux


def test_padding_is_preserved_but_cannot_enter_continuum_support():
 w=np.arange(-8.,9.);f=np.ones(len(w));f[[0,-1]]=0;before=f.copy()
 assert validate_processed_flux(w,f,[0.])==[0,len(w)-1]
 np.testing.assert_array_equal(f,before)
 for offset in [-5.,0.,5.]:
  invalid=f.copy();invalid[w==offset]=0
  with pytest.raises(ValueError):validate_processed_flux(w,invalid,[0.])


def test_nonfinite_and_negative_padding_are_never_accepted():
 w=np.arange(-8.,9.)
 for bad in [np.nan,np.inf,-1.]:
  f=np.ones(len(w));f[0]=bad
  with pytest.raises(ValueError):validate_processed_flux(w,f,[0.])
