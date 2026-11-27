"""RYA-1395 v72 declared controls for rya1395_hydrogenic_rwfn (Sirius):
  C1 our Z_eff=1 4f j=5/2 vs the v60 converged 4f- (runs/mcdhf_s0_sequential_v60/even_5_4f/rwfn.out):
     overlap >= 0.999 and |eps - eps_v60| <= 1e-4.
  C2 round trip: our file read by GRASP rwfnestimate (option 1) in a copy of the v71 4f-step inputs reproduces our
     P/Q (overlap >= 0.99999 per orbital).
"""
import json, math, shutil, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_hydrogenic_rwfn as H

RNT, HH, N = math.exp(-65 / 16) / 8, 0.03125, 439
root = Path(sys.argv[1]); grasp = Path(sys.argv[2]); out = Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=False)
res = {}
ours = {H.label(o): o for o in (H.write_rwfn(out / 'hyd.w', H.build(['4f'], RNT, HH, N)) or H.read_rwfn(out / 'hyd.w'))}
v60 = {H.label(o): o for o in H.read_rwfn(root / 'runs/mcdhf_s0_sequential_v60/even_5_4f/rwfn.out')}
ov = H.overlap(ours['4f-'], v60['4f-']); de = ours['4f-']['eps'] - v60['4f-']['eps']
res['C1'] = dict(overlap=ov, d_eps=de, eps_v60=v60['4f-']['eps'], pass_=bool(abs(ov) >= 0.999 and abs(de) <= 1e-4))
src = root / 'runs/mcdhf_s0_v71_even/even_5_4f'
for f in ('rcsf.inp', 'isodata'):
    shutil.copyfile(src / f, out / f)
for i in range(5):
    shutil.copyfile(src / f'prev{i}.w', out / f'prev{i}.w')
grid = f'{RNT:.17g}\n{HH:.17g}\n0\n{N}\n'
stdin = 'n\nn\n\ny\nn\ny\n' + grid + ''.join(f'1\nprev{i}.w\n*\n' for i in range(5)) + '1\nhyd.w\n*\nn\n'
p = subprocess.run([str(grasp / 'rwfnestimate')], input=stdin, text=True, cwd=out, capture_output=True)
(out / 'rwfnestimate.stdout').write_text(p.stdout); (out / 'rwfnestimate.stderr').write_text(p.stderr)
back = {H.label(o): o for o in H.read_rwfn(out / 'rwfn.inp')}
c2 = {lab: H.overlap(ours[lab], back[lab]) for lab in ('4f-', '4f')}
res['C2'] = dict(rc=p.returncode, overlaps=c2, eps_back={lab: back[lab]['eps'] for lab in c2},
                 pass_=bool(p.returncode == 0 and all(abs(v) >= 0.99999 for v in c2.values())))
res['pass_'] = res['C1']['pass_'] and res['C2']['pass_']
(out / 'control.json').write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
sys.exit(0 if res['pass_'] else 1)
