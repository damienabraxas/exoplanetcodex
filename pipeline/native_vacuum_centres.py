"""Read a source-bound native photon-centre packet; never fit or reconvert wavelengths."""
import json
from pathlib import Path
import numpy as np
from pipeline.reference_engines import fingerprint,load_engine


def load_centres(item,table,atom):
    from config.constants import require_codex_path
    from pipeline.reference_engines import ROOT
    root=require_codex_path(item['resource']);path=root/'output.txt';report_file=root/'report.json'
    if fingerprint(path)!=item['receipt']or fingerprint(report_file)!=item['report_receipt']:raise ValueError('Changed native vacuum-centre packet')
    report=json.loads(report_file.read_text());_,_,_,engine=load_engine('openmulti3d')
    if report['engine']!=engine or report['output']!=item['receipt']or report['generator']!=fingerprint(ROOT/'scripts/probe_native_wavelength_roundtrip.py'):raise ValueError('Unqualified native wavelength producer')
    checks={str(path):item['receipt'],str(report_file):item['report_receipt']}
    for name,expected in report['inputs'].items():
        if fingerprint(name)!=expected:raise ValueError('Changed native wavelength producer input')
        checks[name]=expected
    if report['inputs'].get(str(table))!=fingerprint(table):raise ValueError('Packet refers to another consumed atom table')
    d=np.loadtxt(path)
    if d.shape!=(len(atom.lines),5)or not np.isfinite(d).all()or np.any(d<=0):raise ValueError('Invalid native photon-centre array')
    expected=np.array([1e8/(atom.levels[x.j-1].E_cm-atom.levels[x.i-1].E_cm)for x in atom.lines]);np.testing.assert_allclose(d[:,0],expected,rtol=1e-14,atol=1e-12)
    if np.max(abs(d[:,2]-expected))>.001:raise ValueError('Native centre packet exceeds declared consumed-table precision')
    result=d[:,2].copy();result.flags.writeable=False
    return result,checks
