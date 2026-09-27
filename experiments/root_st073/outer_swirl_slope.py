"""Axisymmetric swirl time-slope directions supported outside the inner cone."""
import numpy as np
from joined_field import coordinates
from separated_moment_modes import flat_bump
class OuterSwirlSlope:
    def __init__(self,base,coefficients,k0,windows):
        self.base=base;self.windows=tuple(windows);self.a=np.asarray(coefficients).reshape(len(windows),3);self.k0=k0
        for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(base,name))
    def fields(self,points,tau):
        points=np.asarray(points);u,p=self.base.fields(points,tau);r=np.hypot(points[:,0],points[:,1]);safe=np.where(r>0,r,1.)
        co=coordinates(r/np.sqrt(self.nu),points[:,2]/np.sqrt(self.nu),tau,self.inner.h);q,eta=np.asarray(co['q']),np.asarray(co['eta'])
        ri=np.sqrt(2*self.nu*q*self.join_X);y=(r-ri)/((self.ratio-1)*ri);k=-np.log2(2*float(np.asarray(tau).ravel()[0]));swirl=np.zeros(len(points))
        for j,(lo,hi) in enumerate(self.windows):
            bump,_=flat_bump(y,lo,hi)
            for power in range(3):swirl+=self.a[j,power]*np.sqrt(self.nu)*q**(-self.inner.A)*bump*(eta/.3)**power*(k-self.k0)
        u=u.copy();u[:,0]-=swirl*points[:,1]/safe;u[:,1]+=swirl*points[:,0]/safe
        return u,p
