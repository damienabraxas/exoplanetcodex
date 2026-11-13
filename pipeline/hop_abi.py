"""Reflect the pinned HOP header types; unsupported ABI declarations refuse."""
import ctypes as ct
import re


def bind_header(library,header):
    types={'int':ct.c_int,'float':ct.c_float,'double':ct.c_double,'char *':ct.c_char_p,'const char *':ct.c_char_p}
    for name,base in [('int',ct.c_int),('float',ct.c_float),('double',ct.c_double)]:
        types[name+' *']=ct.POINTER(base);types['const '+name+' *']=ct.POINTER(base)
    bound={}
    for result,name,args in re.findall(r'(int|double)\s+(hop_\w+)\s*\(([^)]*)\);',header):
        params=[]
        if args.strip():
            for arg in args.split(','):
                kind=re.sub(r'\b[A-Za-z_]\w*\s*$','',arg).strip();kind=re.sub(r'\s+',' ',kind)
                if kind not in types:raise ValueError('Unsupported native ABI type '+kind)
                params.append(types[kind])
        fun=getattr(library,name);fun.argtypes=params;fun.restype=types[result];bound[name]=fun
    return bound
