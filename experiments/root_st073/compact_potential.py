"""Compact physical-space solenoidal localization via an axisymmetric vector potential.

Defined only on the registered ST073 time slab. The velocity is curl(chi*A),
with A_theta=psi/r and A_z=-integral_0^r u_theta(s,z) ds.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,Legendre
from numpy.polynomial import Polynomial as Poly
from wide_modes import CachedWidth,WideJointModes
from poloidal_bridge import coefficients,septic
from joined_field import ROOT,coordinates
from compact_control import cutoff

class CompactPotentialField:
 def __init__(self,quadrature_order=32,experimental_time_extension=False):
  self.base=CachedWidth(6.);self.nu=self.base.nu;self.ratio=6.
  self.experimental_time_extension=bool(experimental_time_extension)
  self.base.experimental_time_extension=self.experimental_time_extension
  self.a=np.asarray(json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes'])
  self.field=WideJointModes(self.base,self.a)
  self.nodes,self.weights=leggauss(quadrature_order)
 def support(self,tau):
  sn=np.sqrt(self.nu);D=self.base.inner.D;X=3/64
  qmax=tau/(1-.43**2);ri_max=sn*np.sqrt(2*qmax*X)
  rflat,rsupp=6.5*ri_max,8*ri_max
  zflat=sn*(tau/(1-.25**2))**D*.25
  zsupp=sn*qmax**D*.43
  return rflat,rsupp,zflat,zsupp
 def streamfunction(self,r,z,tau):
  sn=np.sqrt(self.nu);co=coordinates(r/sn,z/sn,tau,self.base.inner.h);q=float(co['q']);eta=float(co['eta']);ri=sn*np.sqrt(2*q*3/64);ro=6*ri;X=3/64
  if r<=ri:
   c=self.base.inner.coefficients(eta,q)[2,:,0];s=r*r/(2*self.nu);x=s/q
   return float(self.nu**1.5*q*np.polynomial.polynomial.polyval(x,np.r_[0,c/np.arange(1,len(c)+1)]))
  if r>=ro:return 0.
  _,_,left,_=coefficients(self.base.inner,eta,tau);width=5*ri;y=(r-ri)/width;psi=float(np.polynomial.polynomial.polyval(y,septic(left,width)))
  for j in range(3):
   P=Legendre.basis(j).convert(kind=Poly)(Poly([-1,2]));B=256*Poly([0,1])**4*Poly([1,-1])**4*P
   for m in range(2):psi+=self.a[:12].reshape(2,3,2)[0,j,m]*self.nu**1.5*(3/32)*q**(1-self.base.inner.A)*B(y)*(eta/.3)**m
  return psi
 def az(self,r,z,tau):
  if r==0:return 0.
  radii=r*(self.nodes+1)/2;pts=np.column_stack((radii,np.zeros(len(radii)),np.full(len(radii),z)))
  swirl=self.field.fields(pts,tau)[0][:,1]
  return float(-r/2*np.dot(self.weights,swirl))
 def fields(self,points,tau):
  pts=np.asarray(points,float);ts=np.broadcast_to(tau,(len(pts),));out=np.zeros((len(pts),3));pressure=np.zeros(len(pts));sn=np.sqrt(self.nu)
  for i,(x,t) in enumerate(zip(pts,ts)):
   if t<=0 or (not self.experimental_time_extension and not .5/64<=t<=.5):
    raise ValueError('Registered finite time slab only')
   r=np.hypot(x[0],x[1]);z=x[2];rflat,rsupp,zflat,zsupp=self.support(t)
   if r>=rsupp or abs(z)>=zsupp:continue
   fr,frp=cutoff((r-rflat)/(rsupp-rflat));fz,fzp=cutoff((abs(z)-zflat)/(zsupp-zflat))
   fr=float(fr);fz=float(fz);frp=float(frp)/(rsupp-rflat);fzp=float(fzp)*np.sign(z)/(zsupp-zflat)
   chi=fr*fz;chir=frp*fz;chiz=fr*fzp
   uv,pv=self.field.fields(x[None,:],t);ca,sa=(x[0]/r,x[1]/r) if r>0 else (1.,0.)
   ur=uv[0,0]*ca+uv[0,1]*sa;uth=-uv[0,0]*sa+uv[0,1]*ca;uz=uv[0,2]
   atheta=self.streamfunction(r,z,t)/r if r>0 and (chiz!=0 or chir!=0) else 0.
   az=self.az(r,z,t) if chir!=0 else 0.
   vr=chi*ur-chiz*atheta;vt=chi*uth-chir*az;vz=chi*uz+chir*atheta
   out[i]=[vr*ca-vt*sa,vr*sa+vt*ca,vz];pressure[i]=chi*pv[0]
  return out,pressure

if __name__=='__main__':
 f=CompactPotentialField();print('support k5.5',f.support(.5*2**(-5.5)))
