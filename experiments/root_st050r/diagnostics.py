"""Frozen-field audits on separate grids; no coefficient updates or image matching."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
import pressure_morph as pm
from spacetime import Family, TIMES
from structure_audit import jets

def pressure_quantities(f,raw,s,z,t):
    s,z,t=np.broadcast_arrays(s,z,t);s,z,t=s.ravel(),z.ravel(),t.ravel();out=[]
    a,b,q,fc,_=f.coefficients(raw)
    for i in range(0,len(s),256):
        sl=slice(i,i+256);S=s[sl];D=f.bundle(S,z[sl],t[sl]);A=D['A']@a;As=D['As']@a;Az=D['Az']@a;Cs=D['Cs']@a;Cz=D['Cz']@a;B=D['B']@b;Bs=D['Bs']@b
        Q=f.basis(S,z[sl],f.nr,f.nz,False);T=np.column_stack([Chebyshev.basis(j)(4*t[sl]-2) for j in range(f.nt)])
        co=q.reshape(f.ns,f.nt)@T.T
        ev=lambda key:np.einsum('ij,ji->i',Q[key]@f.Tq,co)
        ps=ev((1,0));pzz=ev((0,2));lap=4*ps+4*S*ev((2,0))+pzz
        tr=(A+2*S*As)**2+A*A+Cz*Cz-2*B*(B+2*S*Bs)+4*S*Az*Cs
        out.append(np.c_[tr+lap,tr,lap,A,B,4*ps,pzz])
    return np.concatenate(out)

def report(candidate,parent,out):
    f,raw=Family.load(candidate);fb,rb=Family.load(parent);rng=np.random.default_rng(9175093)
    R,Z=np.meshgrid(np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],indexing='ij')
    cases=[]
    for name,R,Z,times in [('comparable_core_grid',R.ravel(),Z.ravel(),TIMES),('new_offgrid_core',rng.uniform(.035,.215,480),rng.choice([-1.,1.],480)*rng.uniform(.025,.215,480),np.r_[.25,np.sort(rng.uniform(.25,.75,9)),.75])]:
        rows=[];initial=None
        for t in times:
            tau=1-t;v=jets(f,raw,R*np.sqrt(tau),Z*tau**.495,t);sv=v[:,:3]*[np.sqrt(tau),tau**.505,tau**.505];pg=np.sign(Z)*v[:,7]
            if initial is None:initial=sv.copy()
            rows.append(dict(time=float(t),profile_drift=float(np.linalg.norm(sv-initial)/np.linalg.norm(initial)),inward_fraction=float(np.mean(v[:,0]<0)),positive_swirl_fraction=float(np.mean(v[:,1]>0)),bipolar_fraction=float(np.mean(Z*v[:,2]>0)),radial_pressure_inward_fraction=float(np.mean(v[:,6]>0)),axial_pressure_toward_fraction=float(np.mean(pg>0)),adverse_axial_pressure_rms=float(np.sqrt(np.mean(np.minimum(pg,0)**2)))))
        cases.append(dict(name=name,points=len(R),rows=rows))
    rr=np.sort(rng.uniform(.035,.62,31));tt=np.r_[.25,np.sort(rng.uniform(.25,.75,19)),.75];R,T=np.meshgrid(rr,tt,indexing='ij');v=jets(f,raw,R,0,T);vp=jets(fb,rb,R,0,T)
    ratio=v[...,8]*np.sign(vp[...,8])/np.maximum(abs(vp[...,8]),1e-30)
    shear=dict(points=int(R.size),minimum_signed_ratio_to_parent=float(ratio.min()),sign_reversal_count=int(np.sum(ratio<0)),minimum_midplane_uz=float(v[...,2].min()))
    mid=jets(f,raw,np.linspace(.04,.6,65),0,.5);shear['midplane_shear_max_t05']=float(abs(mid[:,8]).max())
    morphology=[];pprows=[]
    for nr,nz in [(48,72),(64,96)]:
        x,wx=leggauss(nr);z,wz=leggauss(nz);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij');w=(4*np.pi*np.outer(wx,wz)).ravel();s,z=S.ravel(),Z.ravel()
        for t in [.25,.4375,.75]:
            j=jets(f,raw,np.sqrt(s),z,t);omega2=np.sum(j[:,3:6]**2,axis=1);en=w@omega2
            morphology.append(dict(order=[nr,nz],time=t,enstrophy=float(en),omega_z_rms_extent=float(np.sqrt((w*z*z)@omega2/en)),omega_z_fourth_moment=float((w*z**4)@omega2/en),omega_r_rms_extent=float(np.sqrt((w*s)@omega2/en))))
            Q=pressure_quantities(f,raw,s,z,t)
            pprows.append(dict(order=[nr,nz],time=t,div_residual_volume_L2=float(np.sqrt(w@(Q[:,0]**2))),div_residual_sampled_max=float(abs(Q[:,0]).max())))
    axis=[]
    for t in TIMES:
        z=np.array([-.15,-.075,0,.075,.15])*(1-t)**.495;Q=pressure_quantities(f,raw,np.zeros(5),z,t)
        axis.append(dict(time=float(t),z=z.tolist(),swirl_strain_ratio=(abs(Q[:,4])/np.maximum(abs(Q[:,3]),1e-30)).tolist(),residual_divergence=Q[:,0].tolist(),radial_pressure_curvature=Q[:,5].tolist(),axial_pressure_curvature=Q[:,6].tolist()))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),cases=cases,shear=shear,morphology=morphology,pressure_poisson=pprows,axis=axis,pde_validated=False,source_correspondence_verified=False,scope='Separate finite-grid diagnostics. Pressure-Poisson residual is div(R), not the acceptance residual itself; no continuum certificate.')
    pm.atomic_json(out,result);print(json.dumps(dict(file=str(out),core_final=cases[0]['rows'][-1],shear=shear,poisson_fine=pprows[-3:]),indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--parent',required=True);p.add_argument('--out',required=True);a=p.parse_args();report(a.candidate,a.parent,a.out)
