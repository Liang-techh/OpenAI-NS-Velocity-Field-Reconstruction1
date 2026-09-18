"""Frozen-field diagnostics, not an exact source-field correspondence certificate."""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parents[2]
for name in ('root_st030','root_st040','root_st043'):
    sys.path.insert(0,str(ROOT/'experiments'/name))
from spacetime import Family
from structure_audit import jets

def integrate(field,initial,rtol=1e-8,max_step=.02):
    initial=np.asarray(initial,float)
    def rhs(t,y):return np.asarray(field(y.reshape(-1,3),t)).ravel()
    ret=solve_ivp(rhs,(.25,.75),initial.ravel(),method='DOP853',t_eval=np.linspace(.25,.75,51),rtol=rtol,atol=rtol*.01,max_step=max_step)
    if not ret.success:raise ValueError(ret.message)
    return ret.t,ret.y.T.reshape(-1,len(initial),3),ret.nfev

def diagnostics(candidate,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    f,raw=Family.load(candidate);rng=np.random.default_rng(9174791)
    specs=[('prior_comparable_grid',np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],np.array([.25,.3125,.4375,.5625,.6875,.75])),
           ('fresh_offgrid',rng.uniform(.031,.219,19),np.r_[-rng.uniform(.029,.223,19),rng.uniform(.029,.223,19)],np.r_[.25,np.sort(rng.uniform(.25,.75,17)),.75])]
    cases=[]
    for name,rr,zz,times in specs:
        R,Z=np.meshgrid(rr,zz,indexing='ij');initial=None;rows=[]
        for t in times:
            tau=1-t;r=R*np.sqrt(tau);z=Z*tau**.495;j=jets(f,raw,r,z,t);u=j[...,:3];su=u*np.array([tau**.5,tau**.505,tau**.505])
            if initial is None:initial=su.copy()
            pgrad=np.sign(Z)*j[...,7]
            # Independent diagnostic grid, never an optimization point selector.
            points=np.c_[r.ravel(),np.zeros(R.size),z.ravel()]
            res=np.concatenate([f.analytic_residual(raw,points[k:k+256],t) for k in range(0,len(points),256)])
            mz=np.sign(Z.ravel())*(res[:,2]-j[...,7].ravel())
            rows.append(dict(time=float(t),profile_drift=float(np.linalg.norm(su-initial)/np.linalg.norm(initial)),inward_fraction=float(np.mean(u[...,0]<0)),positive_swirl_fraction=float(np.mean(u[...,1]>0)),bipolar_fraction=float(np.mean(Z*u[...,2]>0)),radial_pressure_inward_fraction=float(np.mean(j[...,6]>0)),axial_pressure_toward_fraction=float(np.mean(pgrad>0)),adverse_pressure_rms=float(np.sqrt(np.mean(np.minimum(pgrad,0)**2))),positive_nonpressure_Mz_rms=float(np.sqrt(np.mean(np.maximum(mz,0)**2))),core_residual_sampled_max=float(np.linalg.norm(res,axis=1).max())))
        cases.append(dict(name=name,points_per_time=R.size,rows=rows))
    mid=[]
    for t in [.25,.371,.5,.633,.75]:
        r=np.linspace(.04,.6,65);j=jets(f,raw,r,0,t)
        mid.append(dict(time=t,min_uz=float(j[:,2].min()),max_uz=float(j[:,2].max()),max_abs_radial_axial_shear=float(np.max(abs(j[:,8])))))
    R,Z=np.meshgrid([.05,.10,.16,.20],[-.20,-.16,-.10,-.04,0.,.04,.10,.16,.20],indexing='ij')
    init=np.c_[R.ravel()*np.sqrt(.75),np.zeros(R.size),Z.ravel()*.75**.495]
    field=lambda p,t:f.fields(raw,p,t)[0]
    t,path,nf=integrate(field,init);_,fine,nff=integrate(field,init,rtol=1e-10,max_step=.01)
    r0=np.linalg.norm(init[:,:2],axis=1);rf=np.linalg.norm(fine[-1,:,:2],axis=1)
    angles=np.unwrap(np.arctan2(fine[:,:,1],fine[:,:,0]),axis=0);angle=angles[-1]-angles[0]
    nz=init[:,2]!=0;outward=np.sign(init[nz,2])*(fine[-1,nz,2]-init[nz,2])
    trajectories=dict(seed_count=len(init),inward_final_fraction=float(np.mean(rf<r0)),positive_turn_fraction=float(np.mean(angle>0)),bipolar_outward_displacement_fraction=float(np.mean(outward>0)),radial_final_initial_ratio_range=[float(np.min(rf/r0)),float(np.max(rf/r0))],rotation_radians_range=[float(angle.min()),float(angle.max())],midplane_positive_displacement_fraction=float(np.mean(fine[-1,~nz,2]>0)),integrator_comparison_max=float(np.max(abs(fine-path))),function_evaluations=[nf,nff],scope='Actual nonautonomous particle trajectories on36autonomous seeds; not steady streamlines, source matching or a blow-up proof')
    np.savez_compressed(out/'particle_paths.npz',times=t,initial=init,positions=fine)
    report=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),cases=cases,midplane=mid,trajectories=trajectories,pde_validated=False,source_correspondence_verified=False,scope='Finite sampled structure diagnostics with separately varied integration tolerances. No universal source-prescribed thresholds.')
    (out/'structure.json').write_text(json.dumps(report,indent=2));print(json.dumps(dict(last=cases[0]['rows'][-1],midplane_t05=mid[2],trajectories=trajectories),indent=2),flush=True)
    return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();diagnostics(a.candidate,a.out)
