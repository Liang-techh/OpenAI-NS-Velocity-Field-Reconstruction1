"""Independent finite-grid pressure audit; exact raw velocity invariance is separate."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from pressure_completion import Family,pressure_basis,residual,save
from structure_audit import jets

def audit(parent,child,out,seed=9175493):
    f,b=Family.load(parent);fc,c=Family.load(child)
    if (f.nr,f.nz,f.nt,f.basis_kind)!=(fc.nr,fc.nz,fc.nt,fc.basis_kind):raise ValueError('Changed field family')
    vel_equal=bool(np.array_equal(b[:2*f.n],c[:2*f.n]));force_equal=bool(np.array_equal(b[-2:],c[-2:]))
    rng=np.random.default_rng(seed);R,Z=np.meshgrid(np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],indexing='ij')
    from spacetime import TIMES
    specs=[('prior_comparable_grid',R.ravel(),Z.ravel(),TIMES),('new_offgrid',rng.uniform(.035,.215,800),rng.choice([-1.,1.],800)*rng.uniform(.025,.215,800),np.r_[.25,np.sort(rng.uniform(.25,.75,15)),.75])]
    cases=[]
    for name,R,Z,times in specs:
        rows=[];initial=None
        for t in times:
            tau=1-t;j=jets(f,c,R*np.sqrt(tau),Z*tau**.495,t);base=jets(f,b,R*np.sqrt(tau),Z*tau**.495,t)
            scaled=j[:,:3]*[tau**.5,tau**.505,tau**.505]
            if initial is None:initial=scaled.copy()
            psign=Z*j[:,7];vsign=Z*j[:,2]
            rows.append(dict(time=float(t),velocity_difference_max=float(abs(j[:,:3]-base[:,:3]).max()),inward_fraction=float(np.mean(j[:,0]<0)),swirl_fraction=float(np.mean(j[:,1]>0)),bipolar_fraction=float(np.mean(vsign>0)),radial_pressure_inward_fraction=float(np.mean(j[:,6]>0)),axial_pressure_toward_fraction=float(np.mean(psign>0)),minimum_signed_axial_pressure=float(np.min(np.sign(Z)*j[:,7])),profile_drift=float(np.linalg.norm(scaled-initial)/np.linalg.norm(initial))))
        cases.append(dict(name=name,points=len(R),rows=rows))
    r,t=np.meshgrid(np.sort(rng.uniform(.035,.62,41)),np.r_[.25,np.sort(rng.uniform(.25,.75,23)),.75],indexing='ij')
    j=jets(f,c,r,0,t).reshape(-1,9);base=jets(f,b,r,0,t).reshape(-1,9)
    pts=rng.uniform(-1.9,1.9,(64,3));ts=rng.uniform(.25,.75,64)
    r0=residual(f,b,pts,ts);r1=residual(f,c,pts,ts);xy=np.linalg.norm(pts[:,:2],axis=1)
    theta=np.c_[-pts[:,1]/xy,pts[:,0]/xy,np.zeros(len(xy))]
    theta_error=float(abs(np.sum(theta*(r1-r0),axis=1)).max())
    result=dict(parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),child_sha256=hashlib.sha256(Path(child).read_bytes()).hexdigest(),velocity_coefficients_identical=vel_equal,force_coefficients_identical=force_equal,shear_probe_count=int(r.size),shear_max_difference=float(abs(j[:,8]-base[:,8]).max()),midplane_velocity_max_difference=float(abs(j[:,2]-base[:,2]).max()),azimuthal_residual_difference=theta_error,cases=cases,scope='Finite pressure direction audit. Identical velocity coefficients in identical basis prove unchanged velocity, vorticity, strain, normalized core profiles and exact-flow trajectories. Pressure is a different function. No continuum residual certificate or source identity.',pde_validated=False)
    save(out,result);print('AUDIT',str(child),json.dumps({case['name']:{key:min(row[key] for row in case['rows']) for key in ['axial_pressure_toward_fraction','radial_pressure_inward_fraction','inward_fraction','swirl_fraction','bipolar_fraction']} for case in cases}),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--child',required=True);p.add_argument('--out',required=True);a=p.parse_args();audit(a.parent,a.child,a.out)
