"""Finite annular matching: positive swirl and smooth compact profile edits.
Moment feasibility does NOT imply admissible stress or full nonlinear NS balance.
"""
from __future__ import annotations
from pathlib import Path
import sys,json
import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import chebyshev as ch
sys.path.insert(0,str(Path(__file__).parent/'upstream'))
from core_series import Core,evaluate,deriv
from handoff import moment_polynomials
from heat_exterior import profile as heat_profile,tail_moments
ROOT=Path(__file__).parent

def step(t):
    t=np.asarray(t);out=np.zeros_like(t,dtype=float);out[t>=1]=1
    m=(t>0)&(t<1);a=t[m];out[m]=1/(1+np.exp(np.clip(1/a-1/(1-a),-700,700)))
    return out

def bump(t):
    t=np.asarray(t);o=np.zeros_like(t,dtype=float);m=(t>0)&(t<1)
    o[m]=np.exp(4-1/(t[m]*(1-t[m])))
    return o

class Annulus:
    def __init__(self,Xb=4.,c=.25,core=None,n=28):
        self.core=core or Core.load(ROOT/'data/ST068-I.npz');self.Xc=.25;self.Xb=float(Xb);self.c=float(c);self.h=self.core.h;self.n=n
        self.centers=np.linspace(.15,.85,6);self.width=.29
        # Boundaries capture all smooth bump supports, as well as blend endpoints.
        self.breaks=np.unique(np.clip(np.r_[.25,.65,self.Xb,.25+(self.Xb-.25)*(self.centers-self.width),.25+(self.Xb-.25)*(self.centers+self.width)],.25,self.Xb))
        x,w=leggauss(n);self.x=np.concatenate([a+(x+1)*(b-a)/2 for a,b in zip(self.breaks[:-1],self.breaks[1:])]);self.w=np.concatenate([w*(b-a)/2 for a,b in zip(self.breaks[:-1],self.breaks[1:])])
        self.B=self.basis(self.x);self.mpol=moment_polynomials(self.core)
    def basis(self,X):
        y=(np.asarray(X)-self.Xc)/(self.Xb-self.Xc)
        return bump((y[...,None]-self.centers+self.width)/(2*self.width))*bump(y)[...,None]
    def base(self,X,eta):
        X,e=np.broadcast_arrays(np.asarray(X),np.asarray(eta));dx=X-self.Xc
        # The exact finite ST068 polynomial is used as an autonomous seam extension.
        # Flat C-infinity blending preserves every inner radial/eta derivative.
        # No assertion that the leading equations hold outside X<=.25.
        xc=np.minimum(X,.65)
        fp=self.core.dval('F',xc,e);up=self.core.dval('U',xc,e)
        if np.any(fp<=0):raise ValueError('Nonpositive polynomial extension')
        t=step((X-self.Xc)/.4);heat=heat_profile(X,e,self.c,self.h,n=64)['F']
        f=np.exp((1-t)*np.log(fp)+t*np.log(heat));u=(1-t)*up
        return f,u
    def setup_eta(self,eta):
        f,u=self.base(self.x,eta);inner={k:float(evaluate(v,self.Xc,eta)) for k,v in self.mpol.items()}
        cp,st,di=tail_moments(self.Xb,eta,self.c,self.h,n=96)
        pi0=float(self.core.dval('P',0,eta))
        target=np.array([0.,np.sqrt(2)*self.c*self.Xb**(1-self.h)/(1-self.h)-di,0.,st,-pi0-cp])
        old=np.array([inner[k] for k in ['M','I','J','S','Cp']]);scale=np.maximum(abs(target),[1,1,1,1,1])
        return f,u,old,target,scale
    def moments(self,a,eta,setup=None,jac=False):
        f0,u0,inner,target,scale=self.setup_eta(eta) if setup is None else setup
        f=f0*np.exp(self.B@a[:6]);u=u0+self.B@a[6:];x=self.x;w=self.w;H=2*x*f
        values=inner+np.array([w@u,w@H,w@(u*H),w@(u*u-x*f*f),w@(f*f)])
        r=(values-target)/scale
        if not jac:return r,values,dict(Fmin=float(f.min()),Umin=float(u.min()),Umax=float(u.max()))
        df=f[:,None]*self.B;du=self.B
        J=np.zeros((5,12));J[0,6:]=w@du;J[1,:6]=(w*2*x)@df
        J[2,:6]=(w*2*x*u)@df;J[2,6:]=(w*H)@du
        J[3,:6]=(-2*w*x*f)@df;J[3,6:]=(2*w*u)@du;J[4,:6]=(2*w*f)@df
        return r,J/scale[:,None]
    def fit(self,etas,start=None,max_nfev=500):
        records=[];sol=[];a=np.zeros(12) if start is None else np.array(start)
        for eta in etas:
            setup=self.setup_eta(float(eta));fun=lambda a:self.moments(a,eta,setup,jac=True)
            ret=least_squares(lambda a:fun(a)[0],a,jac=lambda a:fun(a)[1],bounds=(np.r_[np.full(6,-6.),np.full(6,-40.)],np.r_[np.full(6,6.),np.full(6,40.)]),ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=max_nfev)
            a=ret.x;r,values,diag=self.moments(a,eta,setup)
            records.append(dict(eta=float(eta),success=bool(ret.success),message=str(ret.message),nfev=ret.nfev,max_scaled_error=float(abs(r).max()),moments=values.tolist(),targets=setup[3].tolist(),**diag));sol.append(a.copy())
        return np.array(sol),records
    def fields(self,X,eta,coeff):
        X,eta=np.broadcast_arrays(X,eta);f0,u0=self.base(X,eta);B=self.basis(X)
        return f0*np.exp(np.sum(B*coeff[...,:6],axis=-1)),u0+np.sum(B*coeff[...,6:],axis=-1)

if __name__=='__main__':
    import time
    etas=[0,.25,.5,-.25,-.5]
    rows=[]
    for xb in [4.,8.,16.]:
      for c in [.25,.5,1.]:
        t=time.monotonic();s=Annulus(xb,c);co,re=s.fit(etas)
        row=dict(Xb=xb,c=c,elapsed=time.monotonic()-t,records=re,coefficients=co.tolist());rows.append(row)
        print(xb,c,[round(q['max_scaled_error'],8) for q in re],time.monotonic()-t,flush=True)
    p=ROOT/'evidence/feasibility_scan.json'
    if p.exists():raise FileExistsError(p)
    p.write_text(json.dumps(rows,indent=2)+'\n')
