"""Numerical radial continuation of leading profiles, not a full NS solution.

Integrates Eqs 4.7/4.13 on an eta collocation grid from a finite-series seed.
No artificial forcing, damping, clipping or imposed exterior is used.
"""
from dataclasses import dataclass, field
import math
import numpy as np
from scipy.integrate import solve_ivp
from .paper_core_reference import PaperCoreReference
from .paper_core_series import PaperCoreSeries
from .profiles import LeadingProfile
from .velocity import leading_velocity_cartesian


@dataclass
class PaperCoreContinuation:
    eta_nodes: int = 257
    degree: int = 12
    sigma: float = .3
    start: float = .03
    end: float = .41
    rtol: float = 1e-10
    seed: PaperCoreSeries = field(init=False,repr=False)
    solution: object = field(init=False,repr=False)

    def __post_init__(self):
        if not (0 < self.start < self.end <= .41 and 0 < self.rtol < 1):
            raise ValueError('require 0<start<end<=.41 and 0<rtol<1')
        self.seed=PaperCoreSeries(PaperCoreReference(sigma=self.sigma),
                                 maxdegree=self.degree,eta_nodes=self.eta_nodes)
        e=self.seed.grid.eta; ref=self.seed.reference
        self._d=1-e*e; self._L=1-2*ref.h*e*e
        self._a=self.seed.a_values; self._g2=self.seed.g_values**2
        initial=np.array([self.seed.phi_value(self.start,e),
                          self.seed.phi_radial_derivative(self.start,e),
                          self.seed.U(self.start,e),
                          self.seed.U_radial_derivative(self.start,e),
                          self.seed.radial_average_U(self.start,e),
                          self.seed.Pi(self.start,e)])
        self.solution=solve_ivp(self.rhs,(self.start,self.end),initial.ravel(),
            method='DOP853',rtol=self.rtol,atol=self.rtol*.01,
            max_step=.01,dense_output=True)
        if not self.solution.success:
            raise RuntimeError(self.solution.message)

    def rhs(self,X,flat):
        phi,phix,U,Ux,avg,P=flat.reshape(6,self.eta_nodes)
        e=self.seed.grid.eta; ref=self.seed.reference
        A,D=.5+ref.h,.5-ref.h
        phie,Ue,avge,Pe=self.seed.grid.differentiate(np.array([phi,U,avg,P]))
        W=1-2*D*e*avg-self._d*avge
        H=D*e+self._d*U
        px=self._g2*phi*phi
        angular=(W*(phi+X*phix)+ref.h*(1-2*e*U)*phi
                 +H*(phie+self._a*phi))
        axial=(W*X*Ux+A*(1-2*e*U)*U+H*Ue+self._d*Pe
               -4*A*e*P-2*e*X*px)
        return np.array([phix,(angular-4*self._L*phix)/(2*self._L*X),
                         Ux,(axial-2*self._L*Ux)/(2*self._L*X),
                         (U-avg)/X,px]).ravel()

    def _state(self,X):
        if not math.isfinite(X) or not self.start <= X <= self.end:
            raise ValueError('X outside integrated interval')
        return self.solution.sol(X).reshape(6,self.eta_nodes)

    def value(self,component,X,eta,eta_derivative=False):
        if not (math.isfinite(X) and math.isfinite(eta) and 0<=X<=self.end and abs(eta)<=1):
            raise ValueError('outside continuation profile domain')
        if X<self.start:
            methods={0:self.seed.phi_eta if eta_derivative else self.seed.phi_value,
                     2:self.seed.dU_deta if eta_derivative else self.seed.U,
                     4:self.seed.radial_average_dU_deta if eta_derivative else self.seed.radial_average_U,
                     5:self.seed.Pi_eta if eta_derivative else self.seed.Pi}
            return float(methods[component](X,eta))
        row=self._state(X)[component]
        if eta_derivative:
            row=self.seed.grid.differentiate(row)
        return float(self.seed.grid.interpolate(row,eta))

    @property
    def profile(self):
        def F(X,e):
            return self.seed.reference.amplitude(e)*self.value(0,X,e)
        return LeadingProfile(
            E=lambda X,e: math.sqrt(2*X)*F(X,e), F=F,
            U=lambda X,e:self.value(2,X,e),
            dU_deta=lambda X,e:self.value(2,X,e,True),
            Pi=lambda X,e:self.value(5,X,e),
            average_U=lambda X,e:self.value(4,X,e),
            average_dU_deta=lambda X,e:self.value(4,X,e,True),
            name='numerical-radial-continuation',paper_exact=False,
            provenance='Independent finite parameters; spectral eta collocation; DOP853 radial continuation')

    def velocity(self,x,y,z,t):
        if not math.isfinite(t) or not 0<=t<1:
            raise ValueError('require 0<=t<1')
        return leading_velocity_cartesian(x,y,z,t,self.profile,h=self.seed.reference.h)

    def metadata(self):
        return dict(eta_nodes=self.eta_nodes,degree=self.degree,sigma=self.sigma,
                    start=self.start,end=self.end,rtol=self.rtol,
                    nfev=self.solution.nfev,steps=len(self.solution.t),
                    status='numerical leading profile only; exterior and full NS missing',
                    seed=self.seed.metadata())
