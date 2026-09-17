"""Compact full-space candidate from a nonlinear inner streamfunction.

The smooth exterior is an independent construction, not the paper's matched
heat exterior. Multiplying the streamfunction preserves incompressibility;
pressure and full Navier-Stokes acceptance remain separate questions.
"""
from functools import lru_cache
import math
import numpy as np
from .coordinates import similarity_coordinates
from .paper_core_reference import PaperCoreReference
from .paper_core_series import PaperCoreSeries


def cutoff(value,inner,outer):
    """C-infinity plateau and derivative with respect to value."""
    if value<=inner:
        return 1.,0.
    if value>=outer:
        return 0.,0.
    s=(value-inner)/(outer-inner)
    logratio=1/(1-s)-1/s
    if logratio>700:
        return 0.,0.
    if logratio < -700:
        return 1.,0.
    c=1/(1+math.exp(logratio))
    return c,-c*(1-c)*(1/s**2+1/(1-s)**2)/(outer-inner)


class PaperCompactField:
    def __init__(self):
        self.core=PaperCoreSeries(PaperCoreReference(sigma=.5),maxdegree=14,eta_nodes=257)

    def _terms(self,x,y,z,t):
        if not all(math.isfinite(v) for v in (x,y,z,t)) or not 0<=t<1:
            raise ValueError('finite coordinates and 0<=t<1 required')
        r=math.hypot(x,y)
        if r>=2 or abs(z)>=2:
            return None
        c=similarity_coordinates(r,z,t,self.core.reference.h)
        if c.X>=.4:
            return None
        bx,dx=cutoff(c.X,.25,.4)
        br,dr=cutoff(r*r,1.,4.)
        bz,dz=cutoff(z*z,.25,4.)
        B=bx*br*bz
        # Here s=r^2; X_s=1/(2q), while X_z follows Eq4.2.
        Bs=bz*(dx*br/(2*c.q)+bx*dr)
        Xz=-2*c.eta*c.X/(c.q**c.D*c.L)
        Bz=br*(dx*Xz*bz+bx*dz*2*z)
        return c,B,Bs,Bz

    def velocity(self,x,y,z,t):
        terms=self._terms(x,y,z,t)
        if terms is None:
            return np.zeros(3)
        c,B,Bs,Bz=terms
        profile=self.core.profile
        avg=profile.radial_average_U(c.X,c.eta)
        # psi=r^2 H, H=(1/2)q^(-A)avgU; multiply psi, not velocity.
        H=.5*c.q**(-c.A)*avg
        radial=B*profile.radial_flux_factor(c.X,c.eta,c.h,d=c.d,L=c.L)/(2*c.q)-Bz*H
        swirl=B*c.q**(-c.A-.5)*profile.smooth_swirl_factor(c.X,c.eta)
        axial=B*c.q**(-c.A)*profile.U(c.X,c.eta)+2*(x*x+y*y)*Bs*H
        result=np.array([x*radial-y*swirl,y*radial+x*swirl,axial])
        if not np.all(np.isfinite(result)):
            raise OverflowError('velocity exceeds floating-point range')
        return result

    def pressure(self,x,y,z,t):
        terms=self._terms(x,y,z,t)
        if terms is None:
            return 0.
        c,B,_,_=terms
        return B*c.q**(-2*c.A)*float(self.core.Pi(c.X,c.eta))

    def at_points(self,points,t):
        points=np.asarray(points,dtype=float)
        if points.ndim<1 or points.shape[-1]!=3:
            raise ValueError('points must have shape (...,3)')
        times=np.broadcast_to(t,points.shape[:-1])
        values=[self.velocity(*p,float(time)) for p,time in zip(points.reshape(-1,3),times.ravel())]
        return np.array(values).reshape(points.shape)

    def metadata(self):
        return dict(core=self.core.metadata(),support='r<2, |z|<2, X<.4',
                    inner_unchanged='X<=.25, r<=1, |z|<=.5',
                    time_domain='0<=t<1; nonzero initial candidate',
                    construction='psi=q^D*X*avgU; localized as B*psi; swirl localized as B*u_theta',
                    acceptance='independent compact candidate; full momentum and forcing not accepted',
                    correspondence='qualitative structural candidate; source image numerical field unverified')


@lru_cache(maxsize=1)
def default_field():
    return PaperCompactField()


def velocity(x,y,z,t):
    """Return Cartesian (u,v,w) for any finite spatial point, 0<=t<1."""
    return tuple(default_field().velocity(x,y,z,t))
