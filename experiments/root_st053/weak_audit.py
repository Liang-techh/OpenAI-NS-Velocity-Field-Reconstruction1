"""Global weak-moment audit on separately refined quadrature, not interval bounds."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import moment_step
from spacetime import Family,force,TIMES
from pressure_morph import atomic_json
def integrate(f,raw,nr,nz,t):
    from numpy.polynomial.legendre import leggauss
    x,wx=leggauss(nr);z,wz=leggauss(nz);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij')
    s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(wx,wz)).ravel();pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    U=[];R=[]
    for i in range(0,len(s),256):
        U.append(f.fields(raw,pts[i:i+256],t)[0]);R.append(f.analytic_residual(raw,pts[i:i+256],t))
    u=np.concatenate(U);r=np.concatenate(R);phi=pts*np.array([.5,.5,-1.])
    D=w@(u[:,2]**2-.5*(u[:,0]**2+u[:,1]**2))
    lhs=w@np.sum(phi*r,axis=1);F=force(pts,t,*raw[-2:]);force_moment=w@np.sum(phi*F,axis=1)
    E=.5*w@np.sum(u*u,axis=1);alpha=(w@(u[:,2]**2))/(2*E)
    return np.array([D,lhs,force_moment,E,alpha,np.sqrt(w@np.sum(r*r,axis=1))])

def report(path,out):
    f,r=Family.load(path);rows=[]
    for order in [(48,72),(64,96),(80,120)]:
        for t in TIMES:
            d,lhs,fm,E,alpha,L2=integrate(f,r,*order,float(t))
            rows.append(dict(order=list(order),time=float(t),D=float(d),residual_pairing=float(lhs),force_pairing=float(fm),energy=float(E),axial_fraction=float(alpha),quadrature_residual_L2=float(L2),weak_L2_lower_estimate=float(abs(d)/np.sqrt(moment_step.NORM2))))
    result=dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),rows=rows,scope='Analytical identity, floating quadrature estimates only. Not interval certified, not the independent original Cartesian acceptance gate.',pde_validated=False)
    atomic_json(out,result);return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();report(a.candidate,a.out)
