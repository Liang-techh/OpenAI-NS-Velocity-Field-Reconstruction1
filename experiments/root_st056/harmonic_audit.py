"""Frozen harmonic-gradient weak witnesses. No candidate fitting or optimization.
Exact polynomial Gram matrices; velocity integrals remain floating quadrature.
"""
from __future__ import annotations
import argparse,hashlib
from pathlib import Path
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from joint_step import Family,save

def basis(max_degree=8):
    s,z,q=sp.symbols('s z q',real=True);Hs=[];H=[]
    for ell in range(2,max_degree+1):
        p=sp.Poly(sp.legendre(ell,q),q)
        h=sp.expand(sum(coef*z**k[0]*(s+z*z)**((ell-k[0])//2) for k,coef in p.terms()))
        assert sp.expand(4*sp.diff(h,s)+4*s*sp.diff(h,s,2)+sp.diff(h,z,2))==0
        H.append(h);Hs.append(sp.diff(h,s))
    def cylinder_integral(expr):
        total=0
        for (a,b),coef in sp.Poly(sp.expand(expr),s,z).terms():
            if b%2==0:total+=coef*sp.Rational(4**(a+1),a+1)*sp.Rational(2*2**(b+1),b+1)
        return total*sp.pi
    G=sp.Matrix(len(H),len(H),lambda i,j:cylinder_integral(4*s*Hs[i]*Hs[j]+sp.diff(H[i],z)*sp.diff(H[j],z)))
    fun=sp.lambdify((s,z),[(2*sp.diff(h,s)+4*s*sp.diff(h,s,2),2*sp.diff(h,s),sp.diff(h,s,z),sp.diff(h,z,2)) for h in H],'numpy')
    return H,np.array(G,float),fun

def run(candidate,out):
    H,G,fun=basis();f,raw=Family.load(candidate);rows=[];scale=np.sqrt(np.diag(G));Gn=G/scale[:,None]/scale[None,:]
    for n in (48,72,96):
        x,wx=leggauss(n);S,Z=np.meshgrid(2*(x+1),2*x,indexing='ij');s,z=S.ravel(),Z.ravel();r=np.sqrt(s);w=(4*np.pi*np.outer(wx,wx)).ravel();pts=np.c_[r,0*r,z];hess=fun(s,z)
        for t in (.25,.5,.75):
            u=f.fields(raw,pts,t)[0];a,b,c=u.T
            d=np.array([-w@(np.broadcast_to(rr,s.shape)*a*a+np.broadcast_to(tt,s.shape)*b*b+4*r*np.broadcast_to(rz,s.shape)*a*c+np.broadcast_to(zz,s.shape)*c*c) for rr,tt,rz,zz in hess])
            levels=[]
            for degree in (2,4,6,8):
                k=degree-1;v=d[:k]/scale[:k];M=Gn[:k,:k];coeff=np.linalg.solve(M,v)
                bound=abs(v@coeff)/np.sqrt(coeff@M@coeff) if np.linalg.norm(coeff)>0 else 0.
                levels.append(dict(max_degree=degree,L2_lower_bound_estimate=float(bound),gram_condition=float(np.linalg.cond(M))))
            rows.append(dict(order=n,time=t,weak_moments=d.tolist(),bounds=levels))
    result=dict(candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),harmonic_polynomials=[str(h) for h in H],exact_polynomial_gram_evaluated_float=G.tolist_value() if False else G.tolist(),rows=rows,scope='Exact integration-by-parts identity for compact solenoidal u/f and compact p. Gram matrix computed from exact polynomial integrals; velocity moments use floating quadrature, NOT interval-certified lower bounds. Post-freeze diagnosis only; no candidate retuning.',pde_validated=False)
    save(out,result);print(candidate,[(x['time'],x['bounds'][-1]['L2_lower_bound_estimate']) for x in rows[-3:]],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);a=p.parse_args();run(a.candidate,a.out)
