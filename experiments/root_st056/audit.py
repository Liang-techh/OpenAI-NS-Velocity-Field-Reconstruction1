"""Frozen independent diagnostics; no fitting or parameter changes."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss
import joint_step
from joint_step import Family,save,force
from validate import cartesian_residual

def jets(f,raw,r,z,t):
    r,z,t=np.broadcast_arrays(r,z,t);r,z,t=r.ravel(),z.ravel(),t.ravel();a,b,p,fc,_=f.coefficients(raw);rows=[]
    for k in range(0,len(r),256):
        sl=slice(k,k+256);rr=r[sl];s=rr*rr;D=f.bundle(s,z[sl],t[sl])
        A=D['A']@a;B=D['B']@b;C=D['C']@a
        rows.append(np.c_[rr*A,rr*B,C,-rr*(D['Bz']@b),rr*(D['Az']@a-2*D['Cs']@a),2*(B+s*(D['Bs']@b)),2*rr*(D['Qs']@p),D['Qz']@p,2*rr*(D['Cs']@a)])
    return np.concatenate(rows)

def exact_residual(f,raw,x,t):
    t=np.broadcast_to(t,(len(x),));return np.concatenate([f.analytic_residual(raw,x[i:i+256],t[i:i+256]) for i in range(0,len(x),256)])

def moment_audit(candidate,out):
    f,raw=Family.load(candidate);rows=[]
    for n in [48,72,96]:
        x,w=leggauss(n);S,Z=np.meshgrid(2*(x+1),2*x,indexing='ij');s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(w,w)).ravel();pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
        for t in [.25,.5,.75]:
            u,_=f.fields(raw,pts,t);vals=w@np.square(u);D=vals[2]-.5*(vals[0]+vals[1]);row=dict(order=n,time=t,component_integrals=vals.tolist(),anisotropy=float(D),volume_L2_lower_bound_estimate=float(abs(D)/np.sqrt(88*np.pi/3)))
            if n==96:
                R=exact_residual(f,raw,pts,t);lhs=w@(np.sqrt(s)/2*R[:,0]-z*R[:,2]);row.update(independent_weak_residual_integral=float(lhs),identity_mismatch=float(lhs-D),volume_L2_quadrature=float(np.sqrt(w@np.sum(R*R,axis=1))))
            rows.append(row)
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),rows=rows,test_function='(x/2,y/2,-z)',test_function_L2_on_support_cylinder=float(np.sqrt(88*np.pi/3)),scope='Exact weak identity under smooth compact u,p and divergence-free prescribed force. Numerical integrals are converged quadrature estimates, not interval-certified lower bounds. Original sampled validation is separate.',pde_validated=False)
    save(out,result);return result

def structure(candidate,parent,out):
    f,raw=Family.load(candidate);fb,rb=Family.load(parent);rng=np.random.default_rng(9205693)
    R=rng.uniform(.035,.215,800);Z=rng.choice([-1.,1.],800)*rng.uniform(.025,.215,800);times=np.r_[.25,np.sort(rng.uniform(.25,.75,15)),.75];rows=[];v0=None
    for t in times:
        tau=1-t;v=jets(f,raw,R*np.sqrt(tau),Z*tau**.495,t);vp=jets(fb,rb,R*np.sqrt(tau),Z*tau**.495,t);scaled=v[:,:3]*[np.sqrt(tau),tau**.505,tau**.505]
        if v0 is None:v0=scaled.copy()
        rows.append(dict(time=float(t),inward_fraction=float(np.mean(v[:,0]<0)),positive_swirl_fraction=float(np.mean(v[:,1]>0)),bipolar_fraction=float(np.mean(Z*v[:,2]>0)),radial_pressure_fraction=float(np.mean(v[:,6]>0)),axial_pressure_fraction=float(np.mean(Z*v[:,7]>0)),min_signed_axial_pressure_gradient=float(np.min(np.sign(Z)*v[:,7])),profile_drift=float(np.linalg.norm(scaled-v0)/np.linalg.norm(v0)),core_swirl_relative_max=float(np.max(abs(v[:,1]-vp[:,1])/np.maximum(abs(vp[:,1]),1e-8))),core_poloidal_relative_L2=float(np.linalg.norm(v[:,[0,2]]-vp[:,[0,2]])/np.linalg.norm(vp[:,[0,2]]))))
    rr,tt=np.meshgrid(rng.uniform(.035,.62,41),np.r_[.25,np.sort(rng.uniform(.25,.75,23)),.75],indexing='ij')
    v=jets(f,raw,rr,0,tt);vp=jets(fb,rb,rr,0,tt)
    sh=dict(points=len(v),maximum_axial_velocity_difference=float(abs(v[:,2]-vp[:,2]).max()),maximum_shear_difference=float(abs(v[:,8]-vp[:,8]).max()),minimum_signed_shear_ratio=float(np.min(v[:,8]*np.sign(vp[:,8])/np.maximum(abs(vp[:,8]),1e-12))),minimum_midplane_uz=float(v[:,2].min()))
    morph=[]
    for nr,nz in [(48,72),(64,96)]:
        x,w=leggauss(nr);z,wz=leggauss(nz);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij');s,z=S.ravel(),Z.ravel();ww=(4*np.pi*np.outer(w,wz)).ravel()
        for t in [.25,.5,.75]:
            v=jets(f,raw,np.sqrt(s),z,t);en=np.sum(v[:,3:6]**2,axis=1);total=ww@en
            morph.append(dict(order=[nr,nz],time=t,enstrophy=float(total),radial_rms=float(np.sqrt((ww*s)@en/total)),axial_rms=float(np.sqrt((ww*z*z)@en/total)),axial_fourth=float((ww*z**4)@en/total)))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),seed=9205693,core_rows=rows,shear=sh,morphology=morph,pde_validated=False,source_correspondence_verified=False,scope='Fresh finite core/shear grids and separate enstrophy quadratures; not full-domain structural certification. Core R[.035,.215],|Z|[.025,.215].')
    save(out,result);return result

def dense(candidate,out):
    f,raw=Family.load(candidate)
    r=np.unique(np.r_[np.linspace(0,1.995,41),np.linspace(.03,.55,14)])
    z=np.unique(np.r_[np.linspace(-1.995,1.995,67),np.linspace(-1.985,-1.71,11),np.linspace(1.71,1.985,11)])
    R,Z=np.meshgrid(r,z,indexing='ij');pts=np.c_[R.ravel(),np.zeros(R.size),Z.ravel()];rows=[]
    for t in np.linspace(.25,.75,13):
        res=exact_residual(f,raw,pts,t);n=np.linalg.norm(res,axis=1);idx=int(n.argmax());p=pts[idx:idx+1]
        fd,_=cartesian_residual(lambda x,t:f.fields(raw,x,t),lambda x,t:force(x,t,*raw[-2:]),p,float(t),.00125,.000625)
        rows.append(dict(time=float(t),sampled_max=float(n[idx]),point=p[0].tolist(),residual=res[idx].tolist(),independent_FD_vector_error=float(np.linalg.norm(fd[0]-res[idx]))))
    save(out,dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),grid=[len(r),len(z),13],rows=rows,scope='New axis-inclusive physical grid plus independent Cartesian FD at each sampled peak. Not a continuum supremum bound.',pde_validated=False));return rows

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--kind',choices=['moment','structure','dense'],required=True);p.add_argument('--out',required=True);p.add_argument('--parent');a=p.parse_args()
    if a.kind=='moment':moment_audit(a.candidate,a.out)
    elif a.kind=='dense':dense(a.candidate,a.out)
    else:structure(a.candidate,a.parent,a.out)
