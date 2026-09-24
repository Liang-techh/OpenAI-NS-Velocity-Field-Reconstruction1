"""Equation-preserving radial recurrence in Y=Lambda X and Phi=F/g(eta).

Finite local prototype based on the independently inspected source B.1-B.3.
Axis jets are expanded about each requested eta, NOT fitted by a global
power polynomial. Long-double arithmetic limits coefficient cancellation.
A radial truncation is still only an approximate leading solution. No global
field, heat matching, stress waves, energy normalization or PDE pass implied.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from numpy.polynomial import polynomial as pp

DT=np.longdouble
@dataclass(frozen=True)
class Parameters:
    lam:float=64.
    sigma:float=.025
    pressure:float=4.
    g_peak:float=.05
    bias:float=.05
    h:float=.005
    order:int=26
    jet_degree:int=32
    Y_max:float=4.
    eta_max:float=.5
    def check(self):
        if not (0<self.h<.01 and self.lam>=1 and self.sigma>0 and self.pressure>0 and self.g_peak>0):raise ValueError('Invalid positive parameters')
        if not (0<self.bias<=.05 and 2<=self.order<=64 and self.jet_degree>=self.order+3 and 0<self.Y_max<=4.1 and 0<self.eta_max<1):raise ValueError('Invalid numerical/domain parameters')

def mul(a,b,n):return np.convolve(a,b)[:n]
def pad(a,n):
    out=np.zeros(n,dtype=DT);out[:min(len(a),n)]=a[:n];return out

def divide(a,b,n):
    a=pad(a,n);b=pad(b,n);out=np.zeros(n,dtype=DT)
    if b[0]==0:raise ValueError('Singular series division')
    for k in range(n):out[k]=(a[k]-np.dot(b[1:k+1],out[k-1::-1]))/b[0] if k else a[0]/b[0]
    return out

def der(a,w):
    out=np.zeros_like(a);out[:-1]=np.arange(1,len(a),dtype=DT)*a[1:]/w;return out

def exponential_jet(a,n):
    a=pad(a,n);out=np.zeros(n,dtype=DT);out[0]=np.exp(a[0])
    for k in range(1,n):out[k]=np.dot(np.arange(1,k+1,dtype=DT)*a[1:k+1],out[k-1::-1])/k
    return out

class NormalizedCore:
    def __init__(self,params:Parameters):
        params.check();self.params=params;self.h=params.h;self.A=.5+self.h;self.D=.5-self.h
        self.eta0=brentq(lambda e:self.Hstar(e),-.1,.1,xtol=1e-15)
        self.meta=dict(id='ST072-N',schema='source_normalized_local_recurrence_v1',h=self.h,nu=.01,X_max=params.Y_max/params.lam,eta_max=params.eta_max,tau_range=[.5/64,.5],parameters=asdict(params),pde_validated=False,leading_validated=False,global_field_ready=False)
    def Hstar(self,e):return self.D*e+(1-e*e)*(4*e+self.params.bias)
    def zeta(self,e):
        H=self.Hstar(e);return -(1-2*self.h*e*e)*H/(H*H+self.params.sigma**2)
    @lru_cache(maxsize=12000)
    def log_g(self,e):
        val,err=quad(self.zeta,self.eta0,float(e),epsabs=2e-13,epsrel=2e-13,limit=160)
        return np.log(DT(self.params.g_peak))+DT(self.params.lam)*DT(val)
    @lru_cache(maxsize=3000)
    def coefficients(self,eta):
        p=self.params;n=p.jet_degree+1;N=p.order;la=DT(p.lam);A=DT(self.A);D=DT(self.D)
        e0=DT(eta)
        if abs(e0)>p.eta_max+1e-12:raise ValueError('Outside declared axial strip')
        # Pure conditioning scale: output is always evaluated at the jet center.
        w=DT(p.sigma/(abs(self.D+4-12*eta*eta-2*p.bias*eta)+1)*.25)
        one=pad([1],n);e=pad([e0,w],n);ee=mul(e,e,n);d=one-ee;L=one-2*DT(p.h)*ee;Li=divide(one,L,n)
        V0=4*e+DT(p.bias)*one;H0=D*e+mul(d,V0,n);ze=divide(-mul(L,H0,n),mul(H0,H0,n)+DT(p.sigma)**2*one,n)
        logG2=pad([2*self.log_g(float(eta))],n)
        logG2[1:]=2*la*w*ze[:-1]/np.arange(1,n,dtype=DT)
        g2=exponential_jet(logG2,n)
        phi=np.zeros((N+1,n),dtype=DT);V=np.zeros_like(phi);P=np.zeros((N+2,n),dtype=DT)
        phi[0]=one;V[0]=V0;P[0]=-DT(p.pressure)*one+DT(.5*p.pressure)*ee
        for k in range(N):
            W=np.zeros((k+1,n),dtype=DT);H=np.zeros_like(W)
            for i in range(k+1):
                W[i]=-(2*D*mul(e,V[i],n)+mul(d,der(V[i],w),n))/(i+1)
                H[i]=mul(d,V[i],n)
            W[0]+=one;H[0]+=D*e
            rf=np.zeros(n,dtype=DT);ru=np.zeros(n,dtype=DT);sphi=np.zeros(n,dtype=DT);hp=np.zeros(n,dtype=DT)
            for i in range(k+1):
                j=k-i
                rf+=(j+1)*mul(W[i],phi[j],n)+mul(H[i],der(phi[j],w),n)-2*DT(p.h)*mul(e,mul(V[i],phi[j],n),n)
                hp+=mul(H[i],phi[j],n)
                ru+=j*mul(W[i],V[j],n)+mul(H[i],der(V[j],w),n)-2*A*mul(e,mul(V[i],V[j],n),n)
                sphi+=mul(phi[i],phi[j],n)
            rf+=DT(p.h)*phi[k]+la*mul(ze,hp,n)
            ru+=A*V[k]+mul(d,der(P[k],w),n)-(4*A+2*k)*mul(e,P[k],n)
            phi[k+1]=mul(Li,rf,n)/(2*la*(k+1)*(k+2))
            V[k+1]=mul(Li,ru,n)/(2*la*(k+1)**2)
            P[k+1]=mul(g2,sphi,n)/(la*(k+1))
        if not np.isfinite(phi).all() or not np.isfinite(V).all():raise FloatingPointError('Nonfinite recurrence')
        avg=V/np.arange(1,N+2,dtype=DT)[:,None];vn=np.empty_like(avg)
        for i in range(N+1):vn[i]=2*mul(e,V[i],n)-2*D*mul(e,avg[i],n)-mul(d,der(avg[i],w),n)
        return phi,V,vn,w,ze
    def point_profiles(self,X,eta):
        p=self.params;Y=DT(p.lam)*DT(X);e=DT(eta);la=DT(p.lam)
        phi,V,vn,w,ze=self.coefficients(float(eta));g=np.exp(self.log_g(float(eta)));xi=la*ze[0];xie=la*ze[1]/w
        def js(c):
            return {(i,j):pp.polyval(Y,pp.polyder(c[:,j],m=i))*la**i*DT([1,1,2,6][j])/w**j for i,j in [(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]}
        phiJ=js(phi);U=js(V);nv=js(vn)
        F={(0,0):g*phiJ[0,0],(1,0):g*phiJ[1,0],(2,0):g*phiJ[2,0],(0,1):g*(phiJ[0,1]+xi*phiJ[0,0]),(1,1):g*(phiJ[1,1]+xi*phiJ[1,0]),(0,2):g*(phiJ[0,2]+2*xi*phiJ[0,1]+(xi*xi+xie)*phiJ[0,0])}
        c0=phi[:,0];c1=phi[:,1]/w;c2=2*phi[:,2]/w**2
        def inte(c):return pp.polyval(Y,pp.polyint(c))
        s0=np.convolve(c0,c0);s1=2*np.convolve(c0,c1);s2=2*np.convolve(c1,c1)+2*np.convolve(c0,c2)
        i0,i1,i2=[inte(c) for c in (s0,s1,s2)];g2=g*g
        pressure=DT(p.pressure);Pi={(0,0):-pressure+pressure*e*e/2+g2/la*i0,(1,0):F[0,0]**2,(2,0):2*F[0,0]*F[1,0],(0,1):pressure*e+g2/la*(i1+2*xi*i0),(1,1):2*F[0,0]*F[0,1],(0,2):pressure+g2/la*(i2+4*xi*i1+(4*xi*xi+2*xie)*i0)}
        L=1-2*DT(self.h)*e*e;inv=1/L;iv1=4*DT(self.h)*e/L**2;iv2=4*DT(self.h)/L**2+32*DT(self.h)**2*e*e/L**3
        v={(0,0):nv[0,0]*inv,(1,0):nv[1,0]*inv,(2,0):nv[2,0]*inv,(0,1):nv[0,1]*inv+nv[0,0]*iv1,(1,1):nv[1,1]*inv+nv[1,0]*iv1,(0,2):nv[0,2]*inv+2*nv[0,1]*iv1+nv[0,0]*iv2}
        av=V/np.arange(1,len(V)+1,dtype=DT)[:,None];av0=pp.polyval(Y,av[:,0]);ave=pp.polyval(Y,av[:,1])/w;d=1-e*e
        W=1-2*DT(self.D)*e*av0-d*ave;H=DT(self.D)*e+d*U[0,0]
        ef=2*L*(X*F[2,0]+2*F[1,0])-(W+self.h*(1-2*e*U[0,0]))*F[0,0]-W*X*F[1,0]-H*F[0,1]
        eu=2*L*(X*U[2,0]+U[1,0])-W*X*U[1,0]-self.A*(1-2*e*U[0,0])*U[0,0]-H*U[0,1]-d*Pi[0,1]+4*self.A*e*Pi[0,0]+2*e*X*Pi[1,0]
        # Conditioning metric guards against suppressing F while faking its equation.
        efscaled=ef/(g*la) if g>0 else DT('nan')
        return dict(F=F,U=U,P=Pi,v=v,leading_F_defect=ef,leading_U_defect=eu,normalized_F_defect=efscaled,log_g=self.log_g(float(eta)),Phi=phiJ[0,0],pressure_defect=Pi[1,0]-F[0,0]**2,divergence_defect=v[0,0]+X*v[1,0]-(2*self.A*e*U[0,0]-d*U[0,1]+2*e*X*U[1,0])/L)
    def profiles(self,X,eta,check=True):
        X,eta=np.broadcast_arrays(np.asarray(X,float),np.asarray(eta,float))
        if not np.isfinite(X).all() or not np.isfinite(eta).all():raise ValueError('Nonfinite coordinates')
        if check and (np.any(X<0) or np.any(X>self.meta['X_max']+1e-13) or np.any(abs(eta)>self.params.eta_max+1e-13)):raise ValueError('Outside local profile')
        rows=[self.point_profiles(x,e) for x,e in zip(X.ravel(),eta.ravel())]
        out={}
        for key in ['F','U','P','v']:
            out[key]={ij:np.array([r[key][ij] for r in rows],dtype=float).reshape(X.shape) for ij in rows[0][key]} if rows else {}
        for key in ['leading_F_defect','leading_U_defect','normalized_F_defect','pressure_defect','divergence_defect','log_g','Phi']:
            out[key]=np.array([r[key] for r in rows],dtype=float).reshape(X.shape)
        return out
    def save(self,path):
        path=Path(path)
        if path.exists():raise FileExistsError(path)
        path.write_text(json.dumps(self.meta,indent=2)+'\n')
    @classmethod
    def load(cls,path):return cls(Parameters(**json.loads(Path(path).read_text())['parameters']))
