"""Fresh structural/dynamical probes; no image matching or continuum certification."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import acceleration_fit
from spacetime import Family,force,NU
from structure_audit import jets

def terms(f,raw,r,z,t):
    r,z,t=np.broadcast_arrays(r,z,t);shape=r.shape;r,z,t=r.ravel(),z.ravel(),t.ravel()
    aa,bb,pp,fc,_=f.coefficients(raw);out=[]
    for start in range(0,len(r),256):
        rr,zz,tt=r[start:start+256],z[start:start+256],t[start:start+256];s=rr*rr;D=f.bundle(s,zz,tt)
        V={k:D[k]@(bb if k.startswith('B') else pp if k.startswith('Q') else aa) for k in ('A','B','C','Ct','Cs','Cz','CL','Qz')}
        xyz=np.c_[rr,np.zeros(len(rr)),zz];fz=force(xyz,tt,*fc)[:,2]
        transport=2*s*V['A']*V['Cs']+V['C']*V['Cz'];visc=-NU*V['CL'];material=V['Ct']+transport+visc-fz
        out.append(np.c_[V['Ct'],transport,visc,-fz,V['Qz'],material,material+V['Qz']])
    return np.concatenate(out).reshape(shape+(7,))

def audit(path):
    f,raw=Family.load(path);cases=[]
    times=np.r_[.25,np.sort(np.random.default_rng(9174691).uniform(.25,.75,17)),.75]
    definitions=[('comparable',np.linspace(.04,.2,11),np.linspace(.04,.2,11),np.array([.25,.3125,.4375,.5625,.6875,.75])),('fresh_offgrid',np.linspace(.035,.205,17),np.linspace(.031,.201,17),times)]
    for name,rg,zg,ts in definitions:
        R,Z=np.meshgrid(rg,np.r_[-zg[::-1],zg],indexing='ij');initial=None;rows=[]
        for t in ts:
            tau=1-t;r=R*tau**.5;z=Z*tau**.495;v=jets(f,raw,r,z,t);sc=np.array([tau**.5,tau**.505,tau**.505]);u=v[...,:3]*sc
            if initial is None:initial=u.copy()
            a=terms(f,raw,r,z,t);sg=np.sign(Z);signed=sg[...,None]*a
            row=dict(time=float(t),profile_drift=float(np.linalg.norm(u-initial)/np.linalg.norm(initial)),radial_inward_fraction=float(np.mean(v[...,0]<0)),positive_swirl_fraction=float(np.mean(v[...,1]>0)),bipolar_fraction=float(np.mean(Z*v[...,2]>0)),radial_pressure_inward_fraction=float(np.mean(v[...,6]>0)),axial_pressure_toward_fraction=float(np.mean(Z*v[...,7]>0)),adverse_pressure_rms=float(np.sqrt(np.mean(np.minimum(signed[...,4],0)**2))),positive_nonpressure_momentum_rms=float(np.sqrt(np.mean(np.maximum(signed[...,5],0)**2))),positive_nonpressure_max=float(np.max(np.maximum(signed[...,5],0))),axial_residual_rms=float(np.sqrt(np.mean(a[...,6]**2))),signed_term_means=dict(zip(['time_derivative','transport','negative_viscosity','negative_force','pressure_gradient','nonpressure_momentum','residual'],np.mean(signed,axis=(0,1)).tolist())))
            rows.append(row)
        cases.append(dict(grid=name,points_per_time=R.size,rows=rows))
    mid=[]
    for t in [.25,.371,.5,.633,.75]:
        v=jets(f,raw,np.linspace(.04,.6,65),0,t);mid.append(dict(time=t,uz_min=float(v[:,2].min()),uz_max=float(v[:,2].max()),max_abs_radial_axial_shear=float(abs(v[:,8]).max())))
    return dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),cases=cases,midplane=mid,scope='Fresh autonomous off-grid structure/time probes, source-alignment diagnostics only; no source identity or continuum residual certificate.',pde_validated=False,source_correspondence_verified=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();d=audit(a.candidate);Path(a.out).write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({v['grid']:v['rows'][-1] for v in d['cases']},indent=2))
