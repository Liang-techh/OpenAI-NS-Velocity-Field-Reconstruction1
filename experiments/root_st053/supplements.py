"""Supplemental finite-grid structure, global moments, and edge maxima after freezing."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from moment_step import ROOT
from aligned_continuation import Family
from pressure_morph import atomic_json
from audit_st051 import audit
from weak_audit import report
from spacetime import force
from validate import cartesian_residual

def edge_report(path,out):
    f,raw=Family.load(path);r=np.linspace(0,1.45,42);z=np.r_[-np.linspace(1.82,1.985,69)[::-1],np.linspace(1.82,1.985,69)]
    R,Z=np.meshgrid(r,z,indexing='ij');x=np.c_[R.ravel(),np.zeros(R.size),Z.ravel()];rows=[]
    field=lambda x,t:f.fields(raw,x,t);forcing=lambda x,t:force(x,t,*raw[-2:])
    for t in np.linspace(.25,.75,13):
        res=np.concatenate([f.analytic_residual(raw,x[i:i+256],t) for i in range(0,len(x),256)])
        norm=np.linalg.norm(res,axis=1);k=int(norm.argmax());fd,_=cartesian_residual(field,forcing,x[k:k+1],t,.00125,.000625)
        rows.append(dict(time=float(t),point=x[k].tolist(),analytic_max=float(norm[k]),FD_max=float(np.linalg.norm(fd[0])),FD_vector_difference=float(np.linalg.norm(res[k]-fd[0]))))
    result=dict(rows=rows,maximum=max(r['analytic_max'] for r in rows),scope='New42x138 axial-edge grid,13times; no radial-edge search or continuum bound; FD at each sampled peak',pde_validated=False)
    atomic_json(out,result);return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',required=True);p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out);freeze=json.loads((out/'FREEZE.json').read_text());parent=ROOT/'artifacts/research/ST052-M/candidate.json';ref=ROOT/'artifacts/research/ST048-S/candidate.json'
    for ident,path in [('ST052-M',parent),('ST053-Q',Path(args.candidate))]:
        assert hashlib.sha256(path.read_bytes()).hexdigest()==freeze['fields'][ident]['sha256']
        d=out/'supplementary';d.mkdir(parents=True,exist_ok=True)
        audit(path,parent,ref,d/f'{ident}_structure.json',seed=9175393)
        report(path,d/f'{ident}_weak.json');print(ident,'WEAK_DONE',flush=True)
        e=edge_report(path,d/f'{ident}_edge.json');print(ident,'EDGE',e['maximum'],flush=True)
if __name__=='__main__':main()
