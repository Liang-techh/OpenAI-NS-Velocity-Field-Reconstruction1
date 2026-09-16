"""Independent finite-difference check of leading-profile Eqs 4.7/4.13.

These are not the full physical Navier-Stokes residual or force acceptance.
"""
import numpy as np


def residual(profile, X, eta, h=.005, step=1e-5):
    if X <= 2*step or abs(eta)+2*step >= 1 or step <= 0:
        raise ValueError('centered stencil requires an interior point')
    if profile.Pi is None:
        raise ValueError('pressure profile required')
    def first(fn,x):
        return (fn(x-2*step)-8*fn(x-step)+8*fn(x+step)-fn(x+2*step))/(12*step)
    def second(fn,x):
        return (-fn(x+2*step)+16*fn(x+step)-30*fn(x)
                +16*fn(x-step)-fn(x-2*step))/(12*step*step)
    F=profile.smooth_swirl_factor(X,eta)
    if F <= 0:
        raise ValueError('positive swirl profile required for logarithmic equation')
    U=profile.U(X,eta); P=profile.Pi(X,eta)
    Fx=first(lambda s:profile.smooth_swirl_factor(s,eta),X)
    Fxx=second(lambda s:profile.smooth_swirl_factor(s,eta),X)
    Fe=first(lambda e:profile.smooth_swirl_factor(X,e),eta)
    Ux=first(lambda s:profile.U(s,eta),X)
    Uxx=second(lambda s:profile.U(s,eta),X)
    Ue=first(lambda e:profile.U(X,e),eta)
    Pe=first(lambda e:profile.Pi(X,e),eta)
    Px=first(lambda s:profile.Pi(s,eta),X)
    avg=profile.radial_average_U(X,eta)
    avge=first(lambda e:profile.radial_average_U(X,e),eta)
    A,D=.5+h,.5-h; d,L=1-eta*eta,1-2*h*eta*eta
    W=1-2*D*eta*avg-d*avge
    Hc=D*eta+d*U
    Sq=-W*(1+X*Fx/F)-h*(1-2*eta*U)-Hc*Fe/F
    Sn=-W*X*Ux-A*(1-2*eta*U)*U-Hc*Ue-d*Pe+4*A*eta*P+2*eta*X*Px
    return np.array([-2*L*(X*Fxx+2*Fx)/F-Sq,
                     -2*L*(X*Uxx+Ux)-Sn,Px-F*F])
