"""Reduced-model derivatives reused by ST041/ST042.
Full rejected ST040 optimizer retained in the user research archive.
"""
from __future__ import annotations

import bootstrap

import argparse,json,time,hashlib

from pathlib import Path

import numpy as np

from scipy.optimize import least_squares

from numpy.polynomial.legendre import leggauss

from spacetime import Family,NU,force

from asymmetric_basis import AsymmetricFamily

VK=('A','As','Az','At','AL','C','Cs','Cz','Ct','CL','B','Bs','Bz','Bt','BL')

def embed(path):
 old,x=Family.load(path);f=AsymmetricFamily(9,12,8);a,b,p,fc,_=old.coefficients(x);part=[]
 for c,O,T in [(a,old.Tp,f.Tp),(b,old.Tw,f.Tw),(p,old.Tq,f.Tq)]:
  p0=(O@c.reshape(old.ns,old.nt)).reshape(old.nr,old.nz,old.nt);p1=np.zeros((9,12,8));p1[:old.nr,:old.nz,:old.nt]=p0
  part.append(np.linalg.solve(T,p1.reshape(f.ns,f.nt)).ravel())
 return f,np.r_[*part,fc]

class Model:
 def norm(self,c):
   e=self.e0+self.EM@c[:self.nv];d=e@e
   if d<1e-20:raise ValueError('Collapsed energy')
   lam=np.sqrt(2/d);dl=-lam*(e@self.EM)/d
   return lam,dl
 def candidate(self,c):
   a=self.a+self.Ma@c[:self.na];b=self.b+self.Mb@c[self.na:self.nv];p=self.p+self.Mp@c[self.nv:self.nv+self.np];fc=self.fc+c[-2:]
   raw=np.r_[a,b,p,fc];f=self.f
   if np.max(abs(raw[:2*f.n]))>4 or np.max(abs(raw[2*f.n:-2]))>100 or np.any((fc<0)|(fc>10)):raise ValueError('Original parameter bound violation')
   return raw
 def cache(self,s,z,t):
   s,z,t=np.broadcast_arrays(s,z,t);s,z,t=s.ravel(),z.ravel(),t.ravel();out={k:[] for k in VK+('Qs','Qz')};outb={k:[] for k in out}
   for i in range(0,len(s),256):
    sl=slice(i,i+256);D=self.f.bundle(s[sl],z[sl],t[sl])
    for k in out:
     base,M=(self.b,self.Mb) if k.startswith('B') else (self.p,self.Mp) if k.startswith('Q') else (self.a,self.Ma)
     out[k].append(D[k]@M);outb[k].append(D[k]@base)
   out={k:np.concatenate(v) for k,v in out.items()};outb={k:np.concatenate(v) for k,v in outb.items()}
   return dict(s=s,z=z,t=t,M=out,v=outb)
 def velocity(self,c,D):
   lam,dl=self.norm(c);n=len(D['s']);r=np.sqrt(D['s']);values=[];mat=[]
   for k,scale in [('A',r),('B',r),('C',np.ones(n))]:
    E=np.zeros((n,self.nv));sl=slice(self.na,self.nv) if k=='B' else slice(0,self.na);E[:,sl]=D['M'][k];bar=D['v'][k]+E@c[:self.nv];V=scale*bar;values.append(lam*V);mat.append(lam*scale[:,None]*E+V[:,None]*dl)
   return np.stack(values,axis=1),np.stack(mat,axis=1)
 def momentum(self,c,D):
   n=len(D['s']);s=D['s'];r=np.sqrt(s);a=c[:self.na];b=c[self.na:self.nv];p=c[self.nv:self.nv+self.np];lam,dl=self.norm(c);v={k:D['v'][k]+D['M'][k]@(b if k.startswith('B') else p if k.startswith('Q') else a) for k in D['M']};M=D['M'];A,B,C=v['A'],v['B'],v['C']
   L=np.column_stack((r*(v['At']-NU*v['AL']),r*(v['Bt']-NU*v['BL']),v['Ct']-NU*v['CL']))
   Cn=np.column_stack((r*(A*A+2*s*A*v['As']+C*v['Az']-B*B),r*(2*A*B+2*s*A*v['Bs']+C*v['Bz']),2*s*A*v['Cs']+C*v['Cz']))
   P=np.column_stack((2*r*v['Qs'],np.zeros(n),v['Qz']));pts=np.column_stack((r,np.zeros(n),D['z']));FA=force(pts,D['t'],1,0);FC=force(pts,D['t'],0,1);fc=self.fc+c[-2:];R=lam*L+lam**2*Cn+P-fc[0]*FA-fc[1]*FC
   J=np.zeros((n,3,self.dim));qa=slice(0,self.na);qb=slice(self.na,self.nv);qp=slice(self.nv,self.nv+self.np)
   J[:,0,qa]=r[:,None]*(lam*(M['At']-NU*M['AL'])+lam**2*((2*A+2*s*v['As'])[:,None]*M['A']+(2*s*A)[:,None]*M['As']+v['Az'][:,None]*M['C']+C[:,None]*M['Az']))
   J[:,0,qb]=-2*lam**2*(r*B)[:,None]*M['B']
   J[:,1,qa]=lam**2*r[:,None]*((2*B+2*s*v['Bs'])[:,None]*M['A']+v['Bz'][:,None]*M['C'])
   J[:,1,qb]=r[:,None]*(lam*(M['Bt']-NU*M['BL'])+lam**2*(2*A[:,None]*M['B']+(2*s*A)[:,None]*M['Bs']+C[:,None]*M['Bz']))
   J[:,2,qa]=lam*(M['Ct']-NU*M['CL'])+lam**2*((2*s*v['Cs'])[:,None]*M['A']+(2*s*A)[:,None]*M['Cs']+v['Cz'][:,None]*M['C']+C[:,None]*M['Cz'])
   J[:,:,:self.nv]+=(L+2*lam*Cn)[:,:,None]*dl[None,None,:]
   J[:,0,qp]=2*r[:,None]*M['Qs'];J[:,2,qp]=M['Qz'];J[:,:,-2]=-FA;J[:,:,-1]=-FC
   return R,J

