"""Sequential replay; preserve evidence and never interpret CLI completion as acceptance."""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import subprocess
import sys


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen',type=Path,default=Path(__file__).resolve().parents[1]/'frozen')
    p.add_argument('--fit',action='store_true',help='Recompute bounded fits instead of validating supplied frozen candidates')
    p.add_argument('--diagnostics',action='store_true',help='Also replay the ineligible nonaxisymmetric capacity screen at three orders')
    a=p.parse_args();src=Path(__file__).resolve().parent;out=a.out.resolve();frozen=a.frozen.resolve()
    if out.exists() and any(out.iterdir()):p.error('Output directory must be empty; previous results are never overwritten')
    names=('ST006',) if a.fit else ('ST006','ST030','ST032','ST033')
    for n in names:
        if not (frozen/n/'candidate.json').is_file():p.error(f'Missing actual frozen input: {frozen/n/"candidate.json"}')
    out.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(script,*args):
        subprocess.run([sys.executable,str(src/script),*map(str,args)],cwd=src,env=env,check=True)
    base=frozen/'ST006/candidate.json';paths={n:frozen/n/'candidate.json' for n in ('ST030','ST032','ST033')}
    if a.fit:
        # Replay the actual segmented budget without intentionally causing an OOM.
        run('mixed_exact_newton.py','--warm',base,'--out',out/'ST030_first42','--iterations',42,'--seed',9172901,'--count',6144)
        run('mixed_exact_newton.py','--warm',out/'ST030_first42/candidate.json','--out',out/'ST030_recovery8','--iterations',8,'--seed',9172901,'--count',6144)
        paths['ST030']=out/'ST030_recovery8/candidate.json'
        run('asymmetric_exact_run.py','--warm',base,'--out',out/'ST032','--iterations',12)
        paths['ST032']=out/'ST032/fit/candidate.json'
        run('peak_refinement.py','--warm',base,'--out',out/'ST033','--rounds',2,'--iterations',20)
        paths['ST033']=out/'ST033/candidate.json'
    val=out/'validation';val.mkdir(exist_ok=True)
    for n,seed in [('ST030',9172910),('ST032',9172911),('ST033',9172912)]:
        run('validate.py',base,'--seed',seed,'--out',val/f'ST006_seed{seed}.json')
        run('validate.py',paths[n],'--seed',seed,'--out',val/f'{n}_seed{seed}.json')
    if a.diagnostics:
        for order in (24,32,48):
            run('annular_refine.py',base,'--out',out/f'ST031_order{order}','--order',order)
    print('Replay finished. Read every validation gate: completion is not PDE acceptance.')
if __name__=='__main__':main()
