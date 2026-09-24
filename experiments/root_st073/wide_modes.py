"""Endpoint-preserving swirl/streamfunction/pressure modes for arbitrary radial width."""
import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import Legendre,leggauss,legval
from width_field import WidthField,coordinates,independent_fd

class CachedWidth(WidthField):
 def __init__(self,ratio=6.):super().__init__(ratio);self.cache={}
 def fields(self,points,tau):
  pts=np.asarray(points,float);ts=np.broadcast_to(tau,(len(pts),));key=(pts.shape,pts.tobytes(),ts.tobytes())
  if key not in self.cache:self.cache[key]=super().fields(pts,ts)
  u,p=self.cache[key];return u.copy(),p.copy()

class WidePoloidalModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.a=np.asarray(amplitudes).reshape(2,3,2)
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
   for m in range(2):
    axial=(eta/.3)**m;axial_z=np.zeros_like(eta) if m==0 else ez/.3
    psi_r+=self.a[0,j,m]*psi_scale*By*axial/width
    psi_z+=self.a[0,j,m]*psi_scale*((1-self.inner.A)*qz/q*B*axial+By*yz*axial+B*axial_z)
    dp+=self.a[1,j,m]*self.nu*q**(-2*self.inner.A)*pressure(y)*axial
  ur=np.where(active,-psi_z/safe,0);uz=np.where(active,psi_r/safe,0);dp=np.where(active,dp,0)
  u[:,0]+=ur*pts[:,0]/safe;u[:,1]+=ur*pts[:,1]/safe;u[:,2]+=uz
  return u,p+dp
class WideAngularModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.amplitudes=np.asarray(amplitudes).reshape(4,3)
 def fields(self,points,tau):
  pts=np.asarray(points,float);u,p=self.base.fields(pts,tau)
  r=np.hypot(pts[:,0],pts[:,1]);sn=np.sqrt(self.nu)
  co=coordinates(r/sn,pts[:,2]/sn,np.broadcast_to(tau,r.shape),self.inner.h)
  q=np.asarray(co['q']);eta=np.asarray(co['eta']);ri=np.sqrt(2*self.nu*q*3/64)
  y=np.clip((r/ri-1)/(self.ratio-1),0,1);amp=np.zeros_like(r)
  for j in range(4):
   radial=legval(2*y-1,[0]*j+[1])
   for m in range(3):amp+=self.amplitudes[j,m]*radial*(eta/.3)**m
  delta=sn*q**(-self.inner.A)*64*y**3*(1-y)**3*amp
  safe=np.where(r>0,r,1);u[:,0]-=delta*pts[:,1]/safe;u[:,1]+=delta*pts[:,0]/safe
  return u,p

class WideJointModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.ratio=base.ratio;self.a=np.asarray(amplitudes)
  self.field=WideAngularModes(WidePoloidalModes(base,self.a[:12]),self.a[12:])
 def fields(self,points,tau):return self.field.fields(points,tau)

def sample(f,k,eta,order,step=.0005):
 tau=.5*2**(-k);g,w=leggauss(order);ip=f.inner.from_similarity([3/64],[eta],tau)[0];ri=ip[0];width=(f.ratio-1)*ri;r=ri+width*(g+1)/2
 pts=np.column_stack((r,np.zeros(order),np.full(order,ip[2])))
 res,div=independent_fd(f,pts,tau,step*np.sqrt(f.nu*tau),.00025*tau)
 mw=w*width/2*r*r/(f.ratio*ri)**2
 return res,div,mw
