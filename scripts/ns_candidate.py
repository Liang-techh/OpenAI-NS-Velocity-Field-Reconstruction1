"""Stable repository entry: status, integrity, evaluation, and independent validation."""
from __future__ import annotations
import argparse,json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from research_baseline import load_best, verify_integrity


def main():
    p=argparse.ArgumentParser(description='ST006 research baseline; not a validated NS solution')
    sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('status');sub.add_parser('verify')
    e=sub.add_parser('evaluate');e.add_argument('--point',type=float,nargs=3,default=[.1,0.,.1]);e.add_argument('--time',type=float,default=.5)
    v=sub.add_parser('validate');v.add_argument('--seed',type=int,default=9172801);v.add_argument('--out',type=Path,default=ROOT/'outputs/ST006_validation.json')
    a=p.parse_args()
    if a.command=='status':
        print((ROOT/'artifacts/research/ST006/manifest.json').read_text());return 0
    if a.command=='verify':
        print(json.dumps(verify_integrity(),indent=2));return 0
    if a.command=='evaluate':
        f=load_best();u,pr=f.fields(a.point,a.time)
        print(json.dumps({'candidate':'ST006','point':a.point,'time':a.time,'velocity':u.tolist(),'pressure':float(pr),'force':f.forcing(a.point,a.time).tolist(),'pde_validated':False},indent=2));return 0
    verify_integrity();out=a.out.resolve()
    if out.exists():p.error('Refusing to overwrite an existing report; choose a new --out')
    out.parent.mkdir(parents=True,exist_ok=True)
    source=ROOT/'experiments/root_st030'
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    subprocess.run([sys.executable,str(source/'validate.py'),str(ROOT/'artifacts/research/ST006/candidate.json'),'--seed',str(a.seed),'--out',str(out)],cwd=source,env=env,check=True)
    return subprocess.run([sys.executable,str(source/'check_acceptance.py'),str(out)],cwd=source,env=env).returncode

if __name__=='__main__':raise SystemExit(main())
