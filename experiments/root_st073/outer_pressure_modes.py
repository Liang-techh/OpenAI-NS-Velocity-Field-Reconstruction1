"""Compact radial pressure directions and a general outer stress cache."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from joined_field import coordinates
from separated_moment_modes import flat_bump
from affine_momentum import jets
WINDOWS=((.12,.38),(.40,.72),(.58,.92))
BREAKS=sorted({v for lo,hi in WINDOWS for v in (lo,(lo+hi)/2,hi)})

class OuterPressure:
    def __init__(self,base,coefficients,windows=None):
        self.base=base;self.windows=tuple(WINDOWS[1:] if windows is None else windows)
        self.coefficients=np.asarray(coefficients).reshape(len(self.windows),3)
        for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(base,name))
    def fields(self,points,tau):
        points=np.asarray(points);u,p=self.base.fields(points,tau)
        r=np.hypot(points[:,0],points[:,1]);co=coordinates(r/np.sqrt(self.nu),points[:,2]/np.sqrt(self.nu),tau,self.inner.h)
        q,eta=np.asarray(co['q']),np.asarray(co['eta']);ri=np.sqrt(2*self.nu*q*self.join_X);y=(r-ri)/((self.ratio-1)*ri)
        p=p.copy()
        for j,(lo,hi) in enumerate(self.windows):
            bump,_=flat_bump(y,lo,hi)
            for power in range(3):p+=self.coefficients[j,power]*self.nu*q**(-2*self.inner.A)*bump*(eta/.3)**power
        return u,p


def outer_cache(field,units,k=11,order=48):
    inner=field.inner;tau=.5*2.**-k;centers=[]
    for eta in (-.2,0.,.2):
        for y in (.5,.75):centers.append(inner.from_similarity([inner.p.X_max*(1+15*y)**2],[eta],tau)[0])
    points=list(centers);panels=[];g,w=leggauss(order)
    for R,_,z in centers:
        q=float(coordinates(0.,z/np.sqrt(inner.nu),tau,inner.h)['q']);ri=np.sqrt(2*inner.nu*q*inner.p.X_max)
        edges=sorted(set(np.clip([0.,ri,R,*[ri*(1+15*b) for b in BREAKS]],0.,R)))
        rr=np.concatenate([(lo+hi)/2+(hi-lo)*g/2 for lo,hi in zip(edges[:-1],edges[1:])]);ww=np.concatenate([(hi-lo)*w/2 for lo,hi in zip(edges[:-1],edges[1:])])
        start=len(points);points.extend(np.column_stack((rr,np.zeros_like(rr),np.full_like(rr,z))));panels.append((slice(start,len(points)),rr,ww,R))
    args=(np.asarray(points),tau,.0005*np.sqrt(inner.nu*tau),.0001*tau)
    baseline=jets(field,*args);uj=[jets(unit,*args) for unit in units]
    return dict(baseline=baseline,modes=tuple(np.stack([j[i] for j in uj]) for i in range(3)),panels=panels,indices=np.arange(len(units)),centers=centers)
