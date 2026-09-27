"""Smooth axisymmetric swirl direction with monotone angular momentum.

Inside its plateau u_theta is proportional to r**(-power), power > 1.
Inner and outer radial cutoffs preserve the existing core and exterior.
This is a candidate mean correction, not a balanced NS solution.
"""
import numpy as np
from joined_field import coordinates
from separated_moment_modes import flat_transition


class BroadAnnularShear:
    def __init__(self, base, amplitude, power=2., onset=(.01,.10),
                 offset=(.88,.99), reference_y=.325):
        self.base, self.amplitude = base, float(amplitude)
        self.power = float(power)
        self.onset, self.offset = tuple(onset), tuple(offset)
        self.reference_y = float(reference_y)
        for name in ('inner','nu','join_X','ratio'):
            setattr(self, name, getattr(base,name))
        if not (self.power>1 and 0<=self.onset[0]<self.onset[1]
                <self.offset[0]<self.offset[1]<=1):
            raise ValueError('Require power>1 and ordered interior radial cutoffs.')

    def fields(self, points, tau):
        points=np.asarray(points,float)
        u,p=self.base.fields(points,tau)
        r=np.hypot(points[:,0],points[:,1])
        co=coordinates(r/np.sqrt(self.nu),points[:,2]/np.sqrt(self.nu),tau,self.inner.h)
        q=np.asarray(co['q'])
        ri=np.sqrt(2*self.nu*q*self.join_X)
        y=(r/ri-1)/(self.ratio-1)
        cutoff=np.array([flat_transition((v-self.onset[0])/(self.onset[1]-self.onset[0]))
                         *flat_transition((self.offset[1]-v)/(self.offset[1]-self.offset[0]))
                         for v in y])
        added=np.zeros(len(points))
        active=cutoff>0
        reference=ri*(1+(self.ratio-1)*self.reference_y)
        added[active]=(self.amplitude*np.sqrt(self.nu)*q[active]**(-self.inner.A)
                       *(reference[active]/r[active])**self.power*cutoff[active])
        safe=np.where(r>0,r,1.)
        u=u.copy()
        u[:,0]-=added*points[:,1]/safe
        u[:,1]+=added*points[:,0]/safe
        return u,p

    def parameters(self):
        return dict(amplitude=self.amplitude,power=self.power,onset=list(self.onset),
                    offset=list(self.offset),reference_y=self.reference_y)
