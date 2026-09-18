"""Low-rank Gaussian-shaped corrections projected into the EXISTING ST042 basis.
No support, field schema, force family or PDE gate is changed. Projection is a
basis-design operation; the actual residual is evaluated on the projected field.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import lstsq
from coupled_fit import CoupledModel
from spacetime import quad,bump_derivatives

class LocalizedModel(CoupledModel):
 def __init__(self,path,radial=2,axial=2,time_degree=4):
  super().__init__(path,radial,axial,time_degree);f=self.f
  s,z,w=quad(96);r=np.sqrt(s);P={k:v@f.Tp for k,v in f.basis(s,z,f.nr,f.nz,True).items()};Q={k:v@f.Tw for k,v in f.basis(s,z,f.nr,f.nz,False).items()}
  Vp=np.vstack((-np.sqrt(w*s)[:,None]*P[0,1],np.sqrt(w)[:,None]*(2*P[0,0]+2*s[:,None]*P[1,0])))
  Vb=np.sqrt(w*s)[:,None]*Q[0,0]
  Vq=np.vstack((2*np.sqrt(w*s)[:,None]*Q[1,0],np.sqrt(w)[:,None]*Q[0,1]))
  rb,rs,*_=bump_derivatives(s/4);rs=rs/4;zb,zd,*_=bump_derivatives(z*z/4);zd=zd*z/2
  new=[];projection=[]
  for block,M,Tm in [('A',self.Ma,Vp),('B',self.Mb,Vb),('Q',self.Mp,Vq)]:
   cols=[]
   for width in (.45,.8):
    for degree in (0,1):
     ir=degree if block!='Q' else 0;iz=1 if block=='A' else 2*degree if block=='Q' else 0
     e=np.exp(-(s+z*z)/width**2);pol=(s/width**2)**ir*(z/width)**iz
     ps=(ir/width**2*(s/width**2)**(ir-1)*(z/width)**iz if ir else np.zeros(len(s)))
     pz=(iz/width*(s/width**2)**ir*(z/width)**(iz-1) if iz else np.zeros(len(s)))
     value=rb*zb*e*pol;ds=rs*zb*e*pol+rb*zb*e*(ps-pol/width**2);dz=rb*zd*e*pol+rb*zb*e*(pz-2*z*pol/width**2)
     target=np.r_[-np.sqrt(w*s)*dz,np.sqrt(w)*(2*value+2*s*ds)] if block=='A' else np.sqrt(w*s)*value if block=='B' else np.r_[2*np.sqrt(w*s)*ds,np.sqrt(w)*dz]
     vec,_,rank,_=lstsq(Tm,target,cond=1e-12,lapack_driver='gelsy');res=float(np.linalg.norm(Tm@vec-target)/np.linalg.norm(target));vec/=np.linalg.norm(vec)
     projection.append(dict(block=block,width=width,radial_degree=ir,axial_degree=iz,projection_relative_L2_error=res,rank=int(rank)))
     for k in range(time_degree):
      mat=np.zeros((f.ns,f.nt));mat[:,k]=vec;cols.append(mat.ravel())
   new.append(np.column_stack((M,np.column_stack(cols))))
  na0,nb0,np0=self.na,self.nb,self.np;oldbounds=self.bounds
  boundgroups=[oldbounds[:na0],oldbounds[na0:na0+nb0],oldbounds[na0+nb0:-2]]
  self.Ma,self.Mb,self.Mp=new;self.na,self.nb,self.np=[M.shape[1] for M in new];self.nv=self.na+self.nb;self.dim=self.nv+self.np+2
  self.bounds=sum([b+[(-.05,.05) if i<2 else (-.15,.15)]*(4*time_degree) for i,b in enumerate(boundgroups)],[])+oldbounds[-2:]
  self.EM=np.zeros((2*f.ns,self.nv));self.EM[:f.ns,:self.na]=np.einsum('stk,t->sk',self.Ma.reshape(f.ns,f.nt,-1),f.q0);self.EM[f.ns:,self.na:]=np.einsum('stk,t->sk',self.Mb.reshape(f.ns,f.nt,-1),f.q0)
  self.projection=projection
