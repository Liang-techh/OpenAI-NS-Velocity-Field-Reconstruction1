"""Add eta-squared streamfunction and pressure modes to the width6 field."""
import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import Legendre
from wide_modes import CachedWidth,WideAngularModes
from joined_field import coordinates

class EvenPoloidalModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.a=np.asarray(amplitudes).reshape(2,3,3)
 def fields(self,points,tau):
  pts=np.asarray(points,float);u,p=self.base.fields(pts,tau)
  r=np.hypot(pts[:,0],pts[:,1]);sn=np.sqrt(self.nu);safe=np.where(r>0,r,1)
  co=coordinates(r/sn,pts[:,2]/sn,np.broadcast_to(tau,r.shape),self.inner.h)
  q=np.asarray(co['q']);eta=np.asarray(co['eta']);qz=np.asarray(co['q_z'])/sn;ez=np.asarray(co['eta_z'])/sn
  ri=np.sqrt(2*self.nu*q*3/64);width=(self.ratio-1)*ri;raw_y=(r-ri)/width;active=(raw_y>0)&(raw_y<1);y=np.clip(raw_y,0,1);yz=-(1+(self.ratio-1)*y)*qz/(2*q*(self.ratio-1))
  psi_scale=self.nu**1.5*(3/32)*q**(1-self.inner.A)
  psi_r=np.zeros_like(r);psi_z=np.zeros_like(r);dp=np.zeros_like(r)
  for j in range(3):
   leg=Legendre.basis(j).convert(kind=Poly)(Poly([-1,2]))
   bubble=256*Poly([0,1])**4*Poly([1,-1])**4*leg
   pressure=64*Poly([0,1])**3*Poly([1,-1])**3*leg
   B=bubble(y);By=bubble.deriv()(y)
   for m in range(3):
    axial=(eta/.3)**m;axial_z=np.zeros_like(eta) if m==0 else m*(eta/.3)**(m-1)*ez/.3
    psi_r+=self.a[0,j,m]*psi_scale*By*axial/width
    psi_z+=self.a[0,j,m]*psi_scale*((1-self.inner.A)*qz/q*B*axial+By*yz*axial+B*axial_z)
    dp+=self.a[1,j,m]*self.nu*q**(-2*self.inner.A)*pressure(y)*axial
  ur=np.where(active,-psi_z/safe,0);uz=np.where(active,psi_r/safe,0);dp=np.where(active,dp,0)
  u[:,0]+=ur*pts[:,0]/safe;u[:,1]+=ur*pts[:,1]/safe;u[:,2]+=uz
  return u,p+dp
class EvenJointModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.a=np.asarray(amplitudes)
  self.field=WideAngularModes(EvenPoloidalModes(base,self.a[:18]),self.a[18:])
 def fields(self,points,tau):return self.field.fields(points,tau)

def expand_old(old):
 old=np.asarray(old);result=np.zeros(30)
 result[:18].reshape(2,3,3)[:,:,:2]=old[:12].reshape(2,3,2)
 result[18:]=old[12:]
 return result
