"""Replay the recorded bounded experiments without overwriting frozen outputs."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import sys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path,default=Path('rerun'))
    args=parser.parse_args();root=Path(__file__).resolve().parent;out=args.out.resolve()
    if out.exists() and any(out.iterdir()):
        parser.error('Output directory must be empty to preserve prior evidence')
    out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    def run(script,*options):
        subprocess.run([sys.executable,str(root/script),*map(str,options)],cwd=root,env=env,check=True)
    warm=None
    for name,nr,nz,nt,moments in [('run_3x3x4',3,3,4,False),('run_4x4x5',4,4,5,False),('run_5x5x6',5,5,6,False),('moment_5x5x6',5,5,6,True),('moment_7x7x8',7,7,8,True)]:
        options=['--out',out/name,'--nr',nr,'--nz',nz,'--nt',nt,'--maxfun',2000]
        if warm:options+=['--warm',warm]
        if moments:options+=['--moments']
        run('spacetime.py',*options);warm=out/name/'candidate.json'
    run('gauss_newton.py','--out',out/'gauss_7x7x8','--warm',warm,'--maxfun',100)
    for name,seed in [('run_5x5x6',9172622),('moment_7x7x8',9172624),('gauss_7x7x8',9172626)]:
        run('validate.py',out/name/'candidate.json','--out',out/name/'validation.json','--seed',seed)
    run('harmonic_moments.py',out/'gauss_7x7x8/candidate.json','--out',out/'gauss_7x7x8/harmonic_moments.json')
    print('Replay complete. Read validation gates: command completion is NOT PDE acceptance.')
if __name__=='__main__':main()
