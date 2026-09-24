"""Source Appendix A.6 radial swirl heat exterior, not a global finite-energy field.
Unit-viscosity K(r,t)=c*s**(-A)*H(2*tau/s); physical u=sqrt(nu)*u1.
Positive generalized Gauss-Laguerre quadrature independently checked against quad.
"""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from scipy.special import roots_genlaguerre,gamma,poch
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

@lru_cache(None)
def laguerre(n,h):
    x,w=roots_genlaguerre(n,h)
    return x,w/gamma(1+h)

def heat_factor(Z,h=.005,order=0,n=96):
    Z=np.asarray(Z,float)
    if np.any(~np.isfinite(Z)) or np.any(Z<0):raise ValueError('Z must be nonnegative and finite')
    v,w=laguerre(n,float(h))
    return (-1.)**order*poch(h,order)*np.sum(w*v**order*np.exp((-h-order)*np.log1p(Z[...,None]*v)),axis=-1)

def factor_minus_one(Z,h=.005,n=96):
    v,w=laguerre(n,float(h))
    return np.sum(w*np.expm1(-h*np.log1p(np.asarray(Z)[...,None]*v)),axis=-1)

def heat_factor_quad(Z,h=.005,order=0):
    ans,err=quad(lambda v: np.exp(-v)*v**(h+order)*(1+Z*v)**(-h-order),0,np.inf,epsabs=1e-12,epsrel=1e-12)
    return (-1.)**order*poch(h,order)*ans/gamma(1+h)

def tail_moments(Xb,eta,c=1.,h=.005,n=96):
    """Infinite-tail moments using Jacobi quadrature to resolve weak power tails.
    Returns Cp_tail, S_tail=int E^2/2, deltaI=int(Hangular-Hpower).
    The deltaI integrand avoids H-1 cancellation using expm1.
    """
    from scipy.special import roots_jacobi
    eta=np.asarray(eta);d=1-eta**2;A=.5+h
    if Xb<=0 or np.any(d<0):raise ValueError('Invalid exterior domain')
    def rule(beta):
        t,w=roots_jacobi(n,0,beta);return (t+1)/2,w/2**(beta+1)
    v,w=rule(2*h)
    H=heat_factor(2*d[...,None]*v/Xb,h,n=n)
    cp=c*c/(2*Xb**(1+2*h))*np.sum(w*H*H,axis=-1)
    v,w=rule(2*h-1)
    H=heat_factor(2*d[...,None]*v/Xb,h,n=n)
    sp=c*c/(2*Xb**(2*h))*np.sum(w*H*H,axis=-1)
    v,w=rule(h-1)
    delta=factor_minus_one(2*d[...,None]*v/Xb,h,n=n)/v
    di=np.sqrt(2.)*c*Xb**(1-h)*np.sum(w*delta,axis=-1)
    return cp,sp,di

def profile(X,eta,c=1.,h=.005,n=96):
    X,eta=np.broadcast_arrays(np.asarray(X,float),np.asarray(eta,float));A=.5+h
    if np.any(X<=0) or np.any(abs(eta)>1):raise ValueError('Exterior requires X>0, |eta|<=1')
    Z=2*(1-eta*eta)/X;H=heat_factor(Z,h,n=n);Hp=heat_factor(Z,h,1,n)
    return dict(E=c*X**(-A)*H,F=c/np.sqrt(2)*X**(-1-h)*H,
      F_X=c/np.sqrt(2)*X**(-2-h)*((-1-h)*H-Z*Hp),
      F_eta=c/np.sqrt(2)*X**(-1-h)*Hp*(-4*eta/X))

def physical(points,tau,c=1.,h=.005,nu=.01,n=96):
    """Pure exterior swirl, radial pressure normalized to vanish at radial infinity.
    Singular/unsupported near axis; do not extend to a whole-space solution.
    """
    p=np.asarray(points,float);tau=np.broadcast_to(np.asarray(tau,float),(len(p),));r=np.hypot(p[:,0],p[:,1]);rho=r/np.sqrt(nu)
    if np.any(r<=0) or np.any(tau<=0):raise ValueError('Use exterior points r>0 and tau>0')
    s=rho*rho/2;A=.5+h;Z=2*tau/s
    H=heat_factor(Z,h,n=n);Hp=heat_factor(Z,h,1,n);Hpp=heat_factor(Z,h,2,n)
    K=c*s**(-A)*H;Kt=-2*c*s**(-A-1)*Hp
    Lk=2*c*s**(-A-1)*(Z*Z*Hpp+(2*A+1)*Z*Hp+(A*A-.25)*H)
    er=np.c_[p[:,0]/r,p[:,1]/r,np.zeros(len(p))];et=np.c_[-p[:,1]/r,p[:,0]/r,np.zeros(len(p))]
    vel=np.sqrt(nu)*K[:,None]*et;grad=(np.sqrt(nu)*K*K/rho)[:,None]*er
    # Pressure integral in source radial s; evaluate q=s,d=tau for the tail factor.
    from scipy.special import roots_jacobi
    x,w=roots_jacobi(n,0,2*h);v=(x+1)/2;w=w/2**(1+2*h)
    f=heat_factor(Z[:,None]*v,h,n=n)
    pres=-nu*c*c/(2*s**(1+2*h))*np.sum(w*f*f,axis=1)
    return dict(velocity=vel,pressure=pres,pressure_gradient=grad,
      residual=np.sqrt(nu)*(Kt-Lk)[:,None]*et,
      time_derivative=np.sqrt(nu)*Kt[:,None]*et,viscous_term=np.sqrt(nu)*Lk[:,None]*et,
      radial_advection=-grad,scope='Pure swirl exterior only; no global axial localization or finite total energy')
