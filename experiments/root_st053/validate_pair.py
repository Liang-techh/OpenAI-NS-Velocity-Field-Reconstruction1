"""Freeze first, then run original independent validators on both fields/seeds."""
import argparse,os,sys,json,hashlib,subprocess,concurrent.futures
from pathlib import Path
from datetime import datetime,timezone
from moment_step import ROOT
from pressure_morph import atomic_json

def worker(item):
    ident,seed,path,folder=item;folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    out=folder/f'{ident}_{seed}.json';log=folder/f'{ident}_{seed}.log'
    code='import sys;sys.path.insert(0,sys.argv[1]);import aligned_continuation;from validate import validate;validate(sys.argv[2],sys.argv[3],int(sys.argv[4]))'
    cmd=[sys.executable,'-c',code,str(ROOT/'experiments/root_st051'),path,str(out),str(seed)]
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    with log.open('w') as fp:ret=subprocess.run(cmd,stdout=fp,stderr=subprocess.STDOUT,env=env)
    if ret.returncode!=0:return dict(id=ident,seed=seed,error='validator process failed',returncode=ret.returncode)
    check=[sys.executable,str(ROOT/'experiments/root_st030/check_acceptance.py'),str(out)]
    accept=subprocess.run(check,capture_output=True,text=True,env=env)
    r=json.loads(out.read_text());fine=[v for v in r['spatial_refinement'] if v['space_step']==.005]
    return dict(id=ident,seed=seed,max=max(v['momentum_max'] for v in fine),volume_L2=max(v['momentum_L2'] for v in fine),gates=r['gates'],acceptance_exit=accept.returncode,acceptance_stdout=accept.stdout)

def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    if (out/'FREEZE.json').exists():raise ValueError('A freeze already exists')
    parent=ROOT/'artifacts/research/ST052-M/candidate.json';paths={'ST052-M':parent,'ST053-Q':Path(args.candidate).resolve()}
    freeze=dict(utc=datetime.now(timezone.utc).isoformat(),fields={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in paths.items()},seeds=[9175391,9175392],no_further_fitting_or_selection=True)
    atomic_json(out/'FREEZE.json',freeze)
    tasks=[(k,s,str(v),str(out/'validation')) for s in freeze['seeds'] for k,v in paths.items()]
    rows=[]
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(worker,tasks):
            rows.append(row);atomic_json(out/'partial_results.json',rows);print(json.dumps(row),flush=True)
    for k,v in paths.items():assert hashlib.sha256(v.read_bytes()).hexdigest()==freeze['fields'][k]['sha256']
    atomic_json(out/'paired_results.json',rows)
if __name__=='__main__':main()
