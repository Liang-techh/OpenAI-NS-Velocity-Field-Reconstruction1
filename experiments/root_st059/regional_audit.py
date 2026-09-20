"""Independent post-freeze mapped-slab residual and fixed-time pressure budgets.
No candidate changes, no resampling based on validation, no alternative gate.
"""
from pathlib import Path
import argparse,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import qr
from frozen_checks import Family,exact_residual,save

def run(candidate,out):
 f,raw=Family.load(candidate);rows=[]
 times=np.r_[.25,np.sort(np.random.default_rng(9205994).uniform(.25,.75,5)),.75]
 for nr,nz in [(40,32),(64,48)]:
  xr,wr=leggauss(nr);xz,wz=leggauss(nz)
  for a,b,label in [(-2,-1.9,'outer_collar'),(-1.9,-1.75,'adjacent_collar'),(-1.75,0,'bulk'),(0,1.75,'bulk'),(1.75,1.9,'adjacent_collar'),(1.9,2,'outer_collar')]:
   S,Z=np.meshgrid(2*(xr+1),(a+b)/2+(b-a)*xz/2,indexing='ij');w=(np.pi*(b-a)*np.outer(wr,wz)).ravel();s,z=S.ravel(),Z.ravel();pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
   for t in times:
    R=exact_residual(f,raw,pts,t);comp=w@(R*R)
    rows.append(dict(time=float(t),order=[nr,nz],slab=[a,b],region=label,component_squared_integrals=comp.tolist(),residual_squared_integral=float(comp.sum()),sampled_max=float(np.linalg.norm(R,axis=1).max())))
 result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),times=times.tolist(),seed=9205994,rows=rows,pde_validated=False,scope='NEW interior times and finer mapped-slab quadrature after freeze; disjoint positive/negative slabs expose compensation. Numerical estimates, not continuous bounds.')
 save(out,result);return result

def pressure_budget(candidate,out):
 f,raw=Family.load(candidate);rows=[]
 for n in (64,96):
  x,w1=leggauss(n);S,Z=np.meshgrid(2*(x+1),2*x,indexing='ij');s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(w1,w1)).ravel();pts=np.c_[np.sqrt(s),np.zeros(len(s)),z]
  basis=f.basis(s,z,f.nr,f.nz,False);br=2*np.sqrt(s)[:,None]*(basis[1,0]@f.Tq);bz=basis[0,1]@f.Tq
  P=np.vstack((np.sqrt(w)[:,None]*br,np.sqrt(w)[:,None]*bz));Q,T=qr(P,mode='economic')
  for t in (.25,.5,.75):
   R=exact_residual(f,raw,pts,t);mer=np.r_[np.sqrt(w)*R[:,0],np.sqrt(w)*R[:,2]];ang=float(w@(R[:,1]**2))
   projection=Q.T@mer;removed=float(projection@projection);remaining=mer-Q@projection
   after=float(remaining@remaining)+ang;before=float(w@np.sum(R*R,axis=1))
   rows.append(dict(order=n,time=t,L2=float(np.sqrt(before)),azimuthal_L2=float(np.sqrt(ang)),unconstrained_fixed_time_pressure_floor_estimate=float(np.sqrt(after)),pressure_removable_fraction=removed/before,meridional_unremoved_L2=float(np.linalg.norm(remaining)),qr_smallest_diagonal=float(np.min(abs(np.diag(T)))),projection_orthogonality=float(np.linalg.norm(Q.T@remaining))))
 save(out,dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),rows=rows,pde_validated=False,scope='Weighted QR least squares over108 existing spatial pressure directions at each fixed time. Velocity and force fixed; pressure/core/coefficient/time-compatibility restrictions discarded ONLY for this diagnostic. Not a proposed candidate or certified continuum lower bound.'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--pressure',action='store_true');a=p.parse_args()
 if a.pressure:pressure_budget(a.candidate,a.out)
 else:run(a.candidate,a.out)
