"""Instantaneous axial-velocity and axial-pressure-gradient zeros, NOT material surfaces."""
from __future__ import annotations
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import aligned_continuation
from spacetime import Family
from structure_audit import jets
from pressure_morph import atomic_json

def compute(path,out):
    f,raw=Family.load(path);R,T=np.meshgrid(np.linspace(.035,.215,13),[.25,.5,.75],indexing='ij')
    r=(R*np.sqrt(1-T)).ravel();times=T.ravel();scale=(1-times)**.495
    Z=np.linspace(-.10,.10,101);vals=jets(f,raw,r[:,None],scale[:,None]*Z,times[:,None])
    records=[]
    for label,column in [('axial_velocity',2),('axial_pressure_gradient',7)]:
        brackets=[];counts=[]
        for v in vals[:,:,column]:
            ids=np.flatnonzero(v[:-1]*v[1:]<=0);counts.append(len(ids))
            brackets.append(int(ids[np.argmin(abs((Z[ids]+Z[ids+1])/2))]) if len(ids) else -1)
        inds=np.flatnonzero(np.array(brackets)>=0);i=np.array(brackets)[inds];a=Z[i].copy();b=Z[i+1].copy();fa=vals[inds,i,column].copy()
        for _ in range(26):
            mid=(a+b)/2;v=jets(f,raw,r[inds],scale[inds]*mid,times[inds])[:,column]
            same=fa*v>0;a=np.where(same,mid,a);fa=np.where(same,v,fa);b=np.where(same,b,mid)
        root=(a+b)/2;value=jets(f,raw,r[inds],scale[inds]*root,times[inds])[:,column]
        rows=[dict(R=float(R.ravel()[j]),time=float(times[j]),r=float(r[j]),scaled_z_root=float(z),physical_z_root=float(scale[j]*z),scan_sign_change_count=int(counts[j]),bracket_width=float(b[k]-a[k]),function_value=float(value[k])) for k,(j,z) in enumerate(zip(inds,root))]
        records.append(dict(quantity=label,rows=rows,missing_roots=int(len(r)-len(rows))))
    result=dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),scaled_z_search=[-.1,.1],records=records,scope='Finite instantaneous zeros. An uz=0 sheet is NOT automatically a Lagrangian/material dividing surface; upward bias implies some nearby negative-z points also move upward by continuity. No source-identification or NS acceptance.')
    atomic_json(out,result);print(json.dumps({r['quantity']:{'roots':len(r['rows']),'scaled_range':[min(v['scaled_z_root'] for v in r['rows']),max(v['scaled_z_root'] for v in r['rows'])]} for r in records}),flush=True);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();compute(a.candidate,a.out)
