"""Derivative-enriched compact basis; original force and support are unchanged.
First 7 one-dimensional columns are the original Legendre-bump columns.
Extra columns b/(1-q)^2, b/(1-q)^4, ... resolve the bump's derivatives.
Every such factor has a smooth zero extension at q=1.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from math import comb
import numpy as np
from numpy.polynomial import Legendre
from scipy.linalg import solve_triangular,qr
from spacetime import Family,bump_derivatives,quad


def rational_derivatives(q,k):
    q=np.asarray(q,float)
    if not np.all(np.isfinite(q)): raise ValueError('Nonfinite basis input')
    d=np.where(q<1,1-q,1.)
    y=1/d;b=np.where(q<1,np.exp(1-y),0.)
    polys={k:1.};out=[]
    for derivative in range(4):
        # Here k <=10 and y on finite float arguments; truncate far below
        # float range before multiplying to avoid 0*overflow near the edge.
        yc=np.minimum(y,1000.)
        out.append(np.where(y<=1000,b*sum(c*yc**power for power,c in polys.items()),0.))
        nxt={}
        for power,c in polys.items():
            nxt[power+1]=nxt.get(power+1,0.)+power*c
            nxt[power+2]=nxt.get(power+2,0.)-c
        polys=nxt
    return out


def hybrid_basis(s,z,nr,nz,odd):
    s,z=np.broadcast_arrays(s,z);s,z=s.ravel(),z.ravel()
    rb=[v/4.**d for d,v in enumerate(bump_derivatives(s/4))]
    b,b1,b2,b3=bump_derivatives(z*z/4)
    zb=[b,b1*z/2,b2*z*z/4+b1/2,b3*z**3/8+3*b2*z/4]
    R=[];Z=[]
    r_extra=[[v/4.**d for d,v in enumerate(rational_derivatives(s/4,2*(i-6)))] for i in range(7,nr)]
    z_extra=[]
    for j in range(7,nz):
        b,b1,b2,b3=rational_derivatives(z*z/4,2*(j-6))
        ez=[b,b1*z/2,b2*z*z/4+b1/2,b3*z**3/8+3*b2*z/4]
        z_extra.append([z*ez[d]/2+(d*ez[d-1]/2 if d else 0.) for d in range(4)] if odd else ez)
    for d in range(4):
        rvals=[sum(comb(d,k)*rb[k]*Legendre.basis(i).deriv(d-k)(s/2-1)/2.**(d-k) for k in range(d+1)) for i in range(min(nr,7))]
        zvals=[sum(comb(d,k)*zb[k]*Legendre.basis(2*j+int(odd)).deriv(d-k)(z/2)/2.**(d-k) for k in range(d+1)) for j in range(min(nz,7))]
        R.append(np.column_stack(rvals+[row[d] for row in r_extra]));Z.append(np.column_stack(zvals+[row[d] for row in z_extra]))
    return {(i,j):(R[i][:,:,None]*Z[j][:,None,:]).reshape(len(s),nr*nz) for i in range(4) for j in range(4-i)}


class HybridFamily(Family):
    basis=staticmethod(hybrid_basis)
    basis_kind='hybrid_legendre7_inverse_even_v1'
    def __post_init__(self):
        if min(self.nr,self.nz,self.nt)<1 or max(self.nr,self.nz,self.nt)>12:raise ValueError('Unsupported order')
        self.ns=self.nr*self.nz;self.n=self.ns*self.nt
        s,z,w=quad(128);P=self.basis(s,z,self.nr,self.nz,True);W=self.basis(s,z,self.nr,self.nz,False)
        A=-P[0,1];C=2*P[0,0]+2*s[:,None]*P[1,0];B=W[0,0]
        MP=np.vstack((np.sqrt(w*s)[:,None]*A,np.sqrt(w)[:,None]*C));MB=np.sqrt(w*s)[:,None]*B
        def orth(M):
            R=qr(M,mode='economic',check_finite=False)[1]
            return solve_triangular(R,np.eye(self.ns),lower=False)
        self.Tp=orth(MP);self.Tw=orth(MB);self.Tq=self.Tw.copy()
        self.Gz=self.Tp.T@(C.T@(w[:,None]*C))@self.Tp;self.Jb=(w*s)@B@self.Tw
        self.q0=np.array([(-1.)**k for k in range(self.nt)])


def embed(warm,out,nr=9,nz=9,nt=8):
    old,raw=Family.load(warm)
    if old.basis_kind=='legendre_v1' and (old.nr>7 or old.nz>7):raise ValueError('Polynomial source larger than preserved block')
    if nr<old.nr or nz<old.nz or nt<old.nt:raise ValueError('Embedding would truncate source')
    f=HybridFamily(nr,nz,nt);aa,bb,qq,fc,_=old.coefficients(raw);coeffs=[]
    for co,To,Tn in [(aa,old.Tp,f.Tp),(bb,old.Tw,f.Tw),(qq,old.Tq,f.Tq)]:
        oldpol=(To@co.reshape(old.ns,old.nt)).reshape(old.nr,old.nz,old.nt)
        padded=np.zeros((nr,nz,nt));padded[:old.nr,:old.nz,:old.nt]=oldpol
        coeffs.append(np.linalg.solve(Tn,padded.reshape(f.ns,nt)).ravel())
    x=np.r_[*coeffs,fc];x[:2*f.n]/=max(1.,np.max(np.abs(x[:2*f.n]))/3.9)
    if np.max(np.abs(x[2*f.n:-2]))>100:raise ValueError('Pressure bounds violated by embedding')
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    f.save(x,out,dict(stage='exact embedding before fitting, modulo independent mass quadrature',source=str(warm)))
    rng=np.random.default_rng(9172740);pts=rng.uniform(-2,2,(512,3));times=rng.uniform(.25,.75,512)
    us,ps=old.fields(raw,pts,times);ut,pt=f.fields(x,pts,times)
    info=dict(velocity_max_difference=float(np.max(np.abs(us-ut))),pressure_max_difference=float(np.max(np.abs(ps-pt))),source_basis=old.basis_kind,target_basis=f.basis_kind,source_parameters=len(raw),target_parameters=len(x),pde_validated=False)
    Path(out).with_name('embedding.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
    return f,x

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('warm');p.add_argument('out');p.add_argument('--nr',type=int,default=9);p.add_argument('--nz',type=int,default=9);p.add_argument('--nt',type=int,default=8);a=p.parse_args();embed(a.warm,a.out,a.nr,a.nz,a.nt)
