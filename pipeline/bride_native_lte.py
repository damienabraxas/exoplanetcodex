"""Serialized O-only calls to the pinned native LTEPOP common-block oracle."""
import ctypes as ct
import json
from pathlib import Path
import threading

import numpy as np

_LOCK=threading.Lock()


class NativeLte:
    def __init__(self,directory):
        from pipeline.bride_reference import fingerprint
        directory=Path(directory);manifest=json.loads((directory/'manifest.json').read_text())
        if manifest['source_commit']!='4d6aa3c39d634b38da772e469a5692435eba1354':
            raise ValueError('Wrong native LTE source')
        library=directory/'libmulti_lte.dylib'
        if fingerprint(library)!=manifest['library']:raise ValueError('Changed native LTE library')
        self.receipt=fingerprint(directory/'manifest.json');self.library=ct.CDLL(str(library))
        d=ct.POINTER(ct.c_double);i=ct.POINTER(ct.c_int)
        self.library.bride_lte_.argtypes=[i,i,d,d,d,d,i,d,d]
        self.library.bride_lte_.restype=None

    def fractions(self,temperature,ne_m3,energies_cm,weights,stages):
        T=np.ascontiguousarray(temperature,dtype=np.float64)
        ne=np.ascontiguousarray(ne_m3,dtype=np.float64)*1e-6
        E=np.ascontiguousarray(energies_cm,dtype=np.float64)
        g=np.ascontiguousarray(weights,dtype=np.float64)
        original_stages=np.asarray(stages)
        stage=np.ascontiguousarray(stages,dtype=np.int32)
        if (T.ndim!=1 or ne.shape!=T.shape or E.ndim!=1 or g.shape!=E.shape
                or stage.shape!=E.shape or not np.array_equal(original_stages,stage)
                or not 1<=len(E)<=64 or not 1<=len(T)<=400
                or not all(np.isfinite(x).all()for x in (T,ne,E,g))
                or np.any(T<=0) or np.any(ne<=0) or E[0]!=0 or np.any(E<0)
                or np.any(g<=0) or stage[0]!=1 or np.any((stage<1)|(stage>2))):
            raise ValueError('Native LTE requires finite O-only atomic arrays within declared capacity')
        n=ct.c_int(len(E));nd=ct.c_int(len(T));output=np.empty((len(E),len(T)),order='F')
        constants=np.empty(8)
        ptr=lambda x:x.ctypes.data_as(ct.POINTER(ct.c_double))
        with _LOCK:
            self.library.bride_lte_(ct.byref(n),ct.byref(nd),ptr(T),ptr(ne),ptr(E),ptr(g),
                stage.ctypes.data_as(ct.POINTER(ct.c_int)),ptr(output),ptr(constants))
        if not np.isfinite(output).all() or np.any(output<=0):raise ValueError('Invalid native LTE fractions')
        return output,dict(zip(('eV_erg','h_erg_s','c_cm_s','k_erg_K','electron_g','amu_g','pi','eV_over_k_K'),constants.tolist()))
