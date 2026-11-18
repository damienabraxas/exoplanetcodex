"""Pinned process-local numerical operator for the native linear ghost response."""
import ctypes,json
from pathlib import Path

def load_ghost_operator(directory,expected):
    from lightweaver import LwCompiled
    from pipeline.reference_engines import fingerprint
    directory=Path(directory);receipt=directory/'receipt.json'
    if fingerprint(receipt)!=expected:raise ValueError('Changed predeclared ghost operator receipt')
    r=json.loads(receipt.read_text());library=directory/'bride_overlap_guard.dylib';name='bride_full_coupling_ghost_scalar'
    if r['schema']!='codex.bride_native_ghost_plugin_build/1'or not r['completed']or r['scheme']!=name or r['method']!='full_coupling_ghost'or r['backend']!=fingerprint(LwCompiled.__file__)or r['plugin']!=fingerprint(library)or library.stat().st_size>5*1024**2:raise ValueError('Unqualified ghost operator build')
    if fingerprint(directory/'Source/SimdFullIterationTemplates.hpp')!=r['patched_template']or fingerprint(directory/'overlap_guard.cpp')!=r['provider']:raise ValueError('Changed ghost numerical recipe')
    ctypes.CDLL(LwCompiled.__file__,mode=ctypes.RTLD_GLOBAL)
    if name not in LwCompiled.FsIterationSchemes.names:LwCompiled.FsIterationSchemes.load_fns_from_path(str(library))
    properties=LwCompiled.FsIterationSchemes.scheme_properties(name)
    if not all(properties[k]for k in ['respectsFormalSolver','defaultPerAtomStorage','defaultWlaGijStorage']):raise ValueError('Invalid ghost operator storage')
    return name,dict(receipt=fingerprint(receipt),plugin=r['plugin'],backend=r['backend'],properties=properties)
