"""CR-ROOT-ST001 compact space-time NS research; NOT an accepted NS solution.
F and B define psi=r^2 F and u_theta=r B; s=r^2. Evaluation: numpy/scipy.
PyTorch is used only for optimizer coefficient gradients, not spatial derivatives.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import time as clock
import numpy as np
from numpy.polynomial import Legendre, Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky, solve_triangular
NU=.01
TIMES=np.array([.25,.3125,.4375,.5625,.6875,.75])

def bump_derivatives(q):
    q=np.asarray(q,float)
    if not np.all(np.isfinite(q)): raise ValueError('Non-finite bump input')
    d=np.where(q<1,1-q,1.)
    b=np.where(q<1,np.exp(1-1/d),0.)
    return b,-b/d**2,b*(d**-4-2*d**-3),b*(-d**-6+6*d**-5-6*d**-4)

def bumped_basis(s,z,nr,nz,odd):
    from math import comb
    s,z=np.broadcast_arrays(s,z);s,z=s.ravel(),z.ravel()
    rb=[v/4.**k for k,v in enumerate(bump_derivatives(s/4))]
    b,b1,b2,b3=bump_derivatives(z*z/4)
    zb=[b,b1*z/2,b2*z*z/4+b1/2,b3*z**3/8+3*b2*z/4]
    R=[];Z=[]
    for d in range(4):
        R.append(np.column_stack([sum(comb(d,k)*rb[k]*Legendre.basis(i).deriv(d-k)(s/2-1)/2.**(d-k) for k in range(d+1)) for i in range(nr)]))
        Z.append(np.column_stack([sum(comb(d,k)*zb[k]*Legendre.basis(2*j+int(odd)).deriv(d-k)(z/2)/2.**(d-k) for k in range(d+1)) for j in range(nz)]))
    return {(i,j):(R[i][:,:,None]*Z[j][:,None,:]).reshape(len(s),nr*nz) for i in range(4) for j in range(4-i)}

def quad(order):
    x,w=leggauss(order);s,z=np.meshgrid(2*(x+1),2*x,indexing='ij')
    return s.ravel(),z.ravel(),(np.pi*np.outer(2*w,2*w)).ravel()

def force(points,t,a,c):
    if not (np.isfinite(a) and np.isfinite(c) and 0<=a<=10 and 0<=c<=10): raise ValueError('Force bound violation')
    x,y,z=np.moveaxis(np.asarray(points,float),-1,0);r2=x*x+y*y
    br,br1,*_=bump_derivatives(r2/4);bz,bz1,*_=bump_derivatives(z*z/4)
    g=bump_derivatives((2*np.asarray(t)-1)**2)[0];b=br*bz
    br_over_r=br1*bz/2;b_z=br*bz1*z/2
    A=-a*(b+z*b_z);B=c*(b+r2*br_over_r/2);C=a*z*(2*b+r2*br_over_r)
    return g[...,None]*np.stack((x*A-y*B,y*A+x*B,C),axis=-1)

@dataclass
class Family:
    basis = staticmethod(bumped_basis)
    basis_kind = "legendre_v1"
    nr:int=3
    nz:int=3
    nt:int=4
    def __post_init__(self):
        if min(self.nr,self.nz,self.nt)<1 or max(self.nr,self.nz,self.nt)>12: raise ValueError('Unsupported order')
        self.ns=self.nr*self.nz;self.n=self.ns*self.nt
        s,z,w=quad(80);P=self.basis(s,z,self.nr,self.nz,True);W=self.basis(s,z,self.nr,self.nz,False)
        A=-P[0,1];C=2*P[0,0]+2*s[:,None]*P[1,0];B=W[0,0]
        Gp=A.T@(w[:,None]*s[:,None]*A)+C.T@(w[:,None]*C);Gw=B.T@(w[:,None]*s[:,None]*B)
        def orth(G): return solve_triangular(cholesky(G,lower=True).T,np.eye(self.ns),lower=False)
        self.Tp=orth(Gp);self.Tw=orth(Gw);self.Tq=self.Tw.copy()
        self.Gz=self.Tp.T@(C.T@(w[:,None]*C))@self.Tp;self.Jb=(w*s)@B@self.Tw
        self.q0=np.array([(-1.)**k for k in range(self.nt)])
    def coefficients(self,raw):
        raw=np.asarray(raw,float)
        if raw.shape!=(3*self.n+2,) or not np.all(np.isfinite(raw)): raise ValueError('Malformed parameters')
        al=raw[:self.n].reshape(self.ns,self.nt);be=raw[self.n:2*self.n].reshape(self.ns,self.nt)
        norm=np.sum((al@self.q0)**2)+np.sum((be@self.q0)**2)
        if norm<=1e-20: raise ValueError('Collapsed initial velocity')
        amp=np.sqrt(2/norm)
        return amp*al.ravel(),amp*be.ravel(),raw[2*self.n:3*self.n],raw[-2:],amp
    def bundle(self,s,z,t):
        s,z,t=np.broadcast_arrays(s,z,t);s,z,t=s.ravel(),z.ravel(),t.ravel()
        P={k:v@self.Tp for k,v in self.basis(s,z,self.nr,self.nz,True).items()}
        B={k:v@self.Tw for k,v in self.basis(s,z,self.nr,self.nz,False).items()}
        Q={k:v@self.Tq for k,v in self.basis(s,z,self.nr,self.nz,False).items()}
        T=np.column_stack([Chebyshev.basis(k)(4*t-2) for k in range(self.nt)])
        Td=np.column_stack([4*Chebyshev.basis(k).deriv()(4*t-2) for k in range(self.nt)])
        def mat(X,dt=False): return (X[:,:,None]*(Td if dt else T)[:,None,:]).reshape(len(s),self.n)
        S=s[:,None];C=2*P[0,0]+2*S*P[1,0]
        return dict(s=s,z=z,t=t,A=mat(-P[0,1]),As=mat(-P[1,1]),Az=mat(-P[0,2]),At=mat(-P[0,1],True),AL=mat(-8*P[1,1]-4*S*P[2,1]-P[0,3]),C=mat(C),Cs=mat(4*P[1,0]+2*S*P[2,0]),Cz=mat(2*P[0,1]+2*S*P[1,1]),Ct=mat(C,True),CL=mat(16*P[1,0]+32*S*P[2,0]+8*S*S*P[3,0]+2*P[0,2]+2*S*P[1,2]),B=mat(B[0,0]),Bs=mat(B[1,0]),Bz=mat(B[0,1]),Bt=mat(B[0,0],True),BL=mat(8*B[1,0]+4*S*B[2,0]+B[0,2]),Q=mat(Q[0,0]),Qs=mat(Q[1,0]),Qz=mat(Q[0,1]))
    def initial(self):
        a=np.zeros((self.ns,self.nt));b=a.copy();v=np.zeros(self.ns);v[0]=.5
        a[:,0]=np.linalg.solve(self.Tp,v);v[0]=1.;b[:,0]=np.linalg.solve(self.Tw,v)
        x=np.r_[a.ravel(),b.ravel(),np.zeros(self.n),0.,0.];x[:2*self.n]/=max(1.,np.max(np.abs(x[:2*self.n])))
        return x
    def fields(self,raw,points,t):
        points=np.asarray(points,float)
        if points.ndim!=2 or points.shape[1]!=3 or not np.all(np.isfinite(points)): raise ValueError('Expected finite (n,3) points')
        t=np.broadcast_to(np.asarray(t,float),(len(points),))
        if np.any((t<.25-1e-13)|(t>.75+1e-13)) or not np.all(np.isfinite(t)): raise ValueError('Time outside window')
        x,y,z=points.T;s=x*x+y*y
        P=self.basis(s,z,self.nr,self.nz,True);W=self.basis(s,z,self.nr,self.nz,False)
        T=np.column_stack([Chebyshev.basis(k)(4*t-2) for k in range(self.nt)])
        a,b,q,fc,_=self.coefficients(raw)
        a=a.reshape(self.ns,self.nt)@T.T;b=b.reshape(self.ns,self.nt)@T.T;q=q.reshape(self.ns,self.nt)@T.T
        A=np.einsum('ij,ji->i',-P[0,1]@self.Tp,a);C=np.einsum('ij,ji->i',(2*P[0,0]+2*s[:,None]*P[1,0])@self.Tp,a)
        B=np.einsum('ij,ji->i',W[0,0]@self.Tw,b);pressure=np.einsum('ij,ji->i',W[0,0]@self.Tq,q)
        return np.column_stack((x*A-y*B,y*A+x*B,C)),pressure
    def analytic_residual(self,raw,points,t):
        x,y,z=np.asarray(points).T;s=x*x+y*y;D=self.bundle(s,z,t);aa,bb,qq,fc,_=self.coefficients(raw)
        v={k:D[k]@(bb if k.startswith('B') else qq if k.startswith('Q') else aa) for k in D if k not in ['s','z','t']}
        A,B,C=v['A'],v['B'],v['C']
        rr=v['At']+A*A+2*s*A*v['As']+C*v['Az']-B*B+2*v['Qs']-NU*v['AL']
        rt=v['Bt']+2*A*B+2*s*A*v['Bs']+C*v['Bz']-NU*v['BL']
        rz=v['Ct']+2*s*A*v['Cs']+C*v['Cz']+v['Qz']-NU*v['CL']
        return np.column_stack((x*rr-y*rt,y*rr+x*rt,rz))-force(points,t,*fc)
    def save(self,raw,path,metadata=None):
        payload=dict(schema='root_st001_compact_spacetime_v1',nr=self.nr,nz=self.nz,nt=self.nt,coefficients=np.asarray(raw).tolist(),metadata=metadata or {},pde_validated=False,paper_exact=False,blowup_proved=False,field_identity_claim=False)
        if self.basis_kind != "legendre_v1": payload["basis_kind"]=self.basis_kind
        Path(path).write_text(json.dumps(payload,indent=2)+'\n')
    @classmethod
    def load(cls,path):
        p=json.loads(Path(path).read_text())
        if p.get('schema')!='root_st001_compact_spacetime_v1' or p.get('pde_validated') is not False: raise ValueError('Unknown schema or unsupported acceptance claim')
        kind=p.get('basis_kind','legendre_v1')
        if kind=='hybrid_legendre7_inverse_even_v1':
            from hybrid_basis import HybridFamily
            cls=HybridFamily
        elif kind=='hybrid9_axial_opposite3_v1':
            from asymmetric_basis import AsymmetricFamily
            cls=AsymmetricFamily
        elif kind!='legendre_v1': raise ValueError('Unknown spatial basis')
        obj=cls(p['nr'],p['nz'],p['nt']);raw=np.array(p['coefficients']);obj.coefficients(raw)
        if np.max(np.abs(raw[:2*obj.n]))>4+1e-12 or np.max(np.abs(raw[2*obj.n:-2]))>100+1e-12 or np.any((raw[-2:]<0)|(raw[-2:]>10)): raise ValueError('Bound violation')
        return obj,raw
