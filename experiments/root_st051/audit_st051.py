"""Frozen actual fields on independent core/shear grids and volume quadratures."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import aligned_continuation
from pressure_morph import atomic_json
from diagnostics import pressure_quantities
from structure_audit import jets
from spacetime import Family,TIMES

def audit(candidate,parent,shear_reference,out,seed=9175193):
    f,raw=Family.load(candidate);fp,rp=Family.load(parent);fs,rs=Family.load(shear_reference)
    rng=np.random.default_rng(seed)
    R,Z=np.meshgrid(np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],indexing='ij')
    cases=[]
    newR=rng.uniform(.035,.215,800);newZ=rng.choice([-1.,1.],800)*rng.uniform(.025,.215,800)
    newtimes=np.r_[.25,np.sort(rng.uniform(.25,.75,15)),.75]
    for name,R,Z,times in [('comparable_core_grid',R.ravel(),Z.ravel(),TIMES),('fresh_offgrid_core',newR,newZ,newtimes)]:
        rows=[];initial=None
        for t in times:
            tau=1-t;j=jets(f,raw,R*np.sqrt(tau),Z*tau**.495,t)
            scaled=j[:,:3]*[np.sqrt(tau),tau**.505,tau**.505]
            if initial is None:initial=scaled.copy()
            pg=np.sign(Z)*j[:,7];cz=np.sign(Z)*j[:,2]
            rows.append(dict(time=float(t),profile_drift=float(np.linalg.norm(scaled-initial)/np.linalg.norm(initial)),inward_fraction=float(np.mean(j[:,0]<0)),positive_swirl_fraction=float(np.mean(j[:,1]>0)),bipolar_fraction=float(np.mean(cz>0)),radial_pressure_inward_fraction=float(np.mean(j[:,6]>0)),axial_pressure_toward_fraction=float(np.mean(pg>0)),min_signed_pressure_gradient=float(pg.min()),min_signed_axial_velocity=float(cz.min()),adverse_axial_pressure_rms=float(np.sqrt(np.mean(np.minimum(pg,0)**2)))))
        cases.append(dict(name=name,points=len(R),rows=rows))
    r,t=np.meshgrid(np.sort(rng.uniform(.035,.62,41)),np.r_[.25,np.sort(rng.uniform(.25,.75,23)),.75],indexing='ij')
    j=jets(f,raw,r,0,t);jp=jets(fp,rp,r,0,t);js=jets(fs,rs,r,0,t)
    def ratio(reference):return j[...,8]*np.sign(reference[...,8])/np.maximum(abs(reference[...,8]),1e-30)
    pr,sr=ratio(jp),ratio(js)
    sh=dict(points=r.size,minimum_signed_ratio_to_P=float(pr.min()),minimum_signed_ratio_to_S=float(sr.min()),violations_P_0999=int(np.sum(pr<.999)),violations_S_0995=int(np.sum(sr<.995)),sign_reversals=int(np.sum(sr<0)),minimum_midplane_uz=float(j[...,2].min()),max_abs_shear_t05=float(abs(jets(f,raw,np.linspace(.04,.6,65),0,.5)[:,8]).max()))
    morphology=[];poisson=[]
    for nr,nz in [(48,72),(64,96)]:
        x,wx=leggauss(nr);z,wz=leggauss(nz);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij');s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(wx,wz)).ravel()
        for t in [.25,.4375,.75]:
            v=jets(f,raw,np.sqrt(s),z,t);om=np.sum(v[:,3:6]**2,axis=1);en=w@om
            morphology.append(dict(order=[nr,nz],time=t,enstrophy=float(en),z_second_moment=float((w*z*z)@om/en),z_fourth_moment=float((w*z**4)@om/en),r_second_moment=float((w*s)@om/en)))
            pp=pressure_quantities(f,raw,s,z,t)
            poisson.append(dict(order=[nr,nz],time=t,divR_L2=float(np.sqrt(w@(pp[:,0]**2))),divR_max=float(abs(pp[:,0]).max())))
    axis=[]
    for t in TIMES:
        z=np.array([-.15,-.075,0,.075,.15])*(1-t)**.495;Q=pressure_quantities(f,raw,np.zeros(5),z,t)
        axis.append(dict(time=float(t),swirl_strain=(abs(Q[:,4])/np.maximum(abs(Q[:,3]),1e-30)).tolist(),divR=Q[:,0].tolist(),radial_pressure_curvature=Q[:,5].tolist(),axial_pressure_curvature=Q[:,6].tolist()))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),shear_reference_sha256=hashlib.sha256(Path(shear_reference).read_bytes()).hexdigest(),seed=seed,cases=cases,shear=sh,morphology=morphology,pressure_poisson=poisson,axis=axis,pde_validated=False,source_correspondence_verified=False,scope='Frozen, finite grid and quadrature diagnostics; no continuum residual or source-identity certificate. divR is not the original residual gate.')
    atomic_json(out,result)
    print(json.dumps(dict(file=str(out),core_final=cases[0]['rows'][-1],new_pressure_min_fraction=min(r['axial_pressure_toward_fraction'] for r in cases[1]['rows']),new_bipolar_min_fraction=min(r['bipolar_fraction'] for r in cases[1]['rows']),shear=sh),indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--parent',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',required=True);p.add_argument('--seed',type=int,default=9175193);a=p.parse_args();audit(a.candidate,a.parent,a.reference,a.out,a.seed)
