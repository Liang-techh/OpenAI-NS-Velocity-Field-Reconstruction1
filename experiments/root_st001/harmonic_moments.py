"""Pressure/force-independent weak residual certificates (quadrature estimates)."""
import argparse
import json
from pathlib import Path
import numpy as np
import sympy as sp
from spacetime import Family,quad,TIMES


def audit(candidate,out):
    f,raw=Family.load(candidate);s,z=sp.symbols('s z',real=True)
    degrees=(2,4,6,8);derivatives=[]
    for degree in degrees:
        H=sp.expand(sum((-1)**k*sp.factorial(2*degree-2*k)*z**(degree-2*k)*(s+z*z)**k/(2**degree*sp.factorial(k)*sp.factorial(degree-k)*sp.factorial(degree-2*k)) for k in range(degree//2+1)))
        hs,hz=sp.diff(H,s),sp.diff(H,z)
        hss,hsz,hzz=sp.diff(H,s,2),sp.diff(H,s,z),sp.diff(H,z,2)
        assert sp.expand(4*hs+4*s*hss+hzz)==0
        derivatives.append([sp.lambdify((s,z),v,'numpy') for v in (hs,hz,hss,hsz,hzz)])
    rows=[]
    for order in (48,96):
        sq,zq,w=quad(order);r=np.sqrt(sq);x=np.column_stack((r,np.zeros(len(r)),zq))
        terms=[tuple(np.broadcast_to(fun(sq,zq),sq.shape) for fun in group) for group in derivatives]
        radial=np.array([2*r*v[0] for v in terms]).T;axial=np.array([v[1] for v in terms]).T
        G=radial.T@(w[:,None]*radial)+axial.T@(w[:,None]*axial)
        scales=np.sqrt(np.diag(G));cor=G/scales[:,None]/scales[None,:]
        for time in TIMES:
            u,_=f.fields(raw,x,time);ur,ut,uz=u.T
            D=np.array([-w@(ur*ur*(2*hs+4*sq*hss)+ut*ut*(2*hs)+4*r*ur*uz*hsz+uz*uz*hzz) for hs,hz,hss,hsz,hzz in terms])
            scaled=D/scales;bound=np.sqrt(max(0.,scaled@np.linalg.solve(cor,scaled)))
            rows.append(dict(order=order,time=float(time),degrees=list(degrees),individual_L2_lower_bound_estimates=np.abs(scaled).tolist(),combined_L2_lower_bound_estimate=float(bound),normalized_gram_condition=float(np.linalg.cond(cor))))
    report=dict(identity='integral grad(H).R = -integral u^T Hess(H) u for every harmonic H',hypotheses='smooth compact u,p, div u=0, compact divergence-free prescribed force',combined_bound='sqrt(D^T G^{-1}D), G_ij=integral_support grad(Hi).grad(Hj)',rows=rows,scope='Exact inequality under hypotheses; numerical values are quadrature estimates, not interval-certified lower bounds or proofs that all velocity families fail.',pde_validated=False)
    Path(out).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps([r for r in rows if r['order']==96],indent=2))
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',default='harmonic_moments.json');a=p.parse_args();audit(a.candidate,a.out)
