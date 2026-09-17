"""Pole-checked Pade evaluation of a nonlinear core Taylor prefix.

Exploratory continuation only. Denominator poles on the integration path are
rejected rather than hidden; eta derivatives remain numerical.
"""
from functools import lru_cache
import math
import numpy as np
from scipy.interpolate import pade
from scipy.integrate import quad
from .profiles import LeadingProfile
from .velocity import leading_velocity_cartesian


class PaperCorePade:
    def __init__(self,series,denominator=None):
        self.series=series
        self.denominator=series.maxdegree//2 if denominator is None else denominator
        if not 0<=self.denominator<=series.maxdegree:
            raise ValueError('denominator degree outside series degree')

    @lru_cache(maxsize=4096)
    def _pair(self,eta,component):
        rows=self.series.phi if component==0 else self.series.u
        c=self.series.grid.interpolate(rows,eta)
        p,q=pade(c,self.denominator)
        roots=q.r
        poles=[float(r.real) for r in roots if abs(r.imag)<1e-8 and r.real>=0]
        return p,q,min(poles,default=float('inf'))

    def _value(self,X,eta,component):
        self.series._broadcast_domain(X,eta)
        p,q,pole=self._pair(float(eta),component)
        if pole<=X:
            raise ValueError('Pade pole on radial integration path')
        value=float(p(X)/q(X))
        if not math.isfinite(value):
            raise ArithmeticError('nonfinite Pade value')
        return value

    def F(self,X,eta):
        return self.series.reference.amplitude(eta)*self._value(X,eta,0)

    def U(self,X,eta):
        return self._value(X,eta,1)

    def Pi(self,X,eta):
        # Check the whole interval before numerical quadrature.
        self._value(X,eta,0)
        return (-self.series.reference.pressure_scale**2/(1+eta*eta)**2
                +quad(lambda s:self.F(s,eta)**2,0,X,epsabs=1e-11,epsrel=1e-11)[0])

    def average(self,X,eta):
        self._value(X,eta,1)
        return quad(lambda s:self.U(s*X,eta),0,1,epsabs=1e-11,epsrel=1e-11)[0]

    @staticmethod
    def eta_derivative(fn,X,eta):
        delta=1e-5
        if abs(eta)+2*delta>1:
            raise ValueError('numerical eta derivative needs interior eta')
        return (fn(X,eta-2*delta)-8*fn(X,eta-delta)
                +8*fn(X,eta+delta)-fn(X,eta+2*delta))/(12*delta)

    @property
    def profile(self):
        return LeadingProfile(E=lambda X,e:math.sqrt(2*X)*self.F(X,e),F=self.F,
            U=self.U,dU_deta=lambda X,e:self.eta_derivative(self.U,X,e),Pi=self.Pi,
            average_U=self.average,
            average_dU_deta=lambda X,e:self.eta_derivative(self.average,X,e),
            name='Pade-core-experiment',paper_exact=False,
            provenance='Pole-checked rational continuation of finite nonlinear core series; numerical eta derivatives')

    def velocity(self,x,y,z,t):
        if not math.isfinite(t) or not 0<=t<1:
            raise ValueError('require 0<=t<1')
        return leading_velocity_cartesian(x,y,z,t,self.profile,h=self.series.reference.h)
