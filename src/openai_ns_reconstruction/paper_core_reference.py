"""Computable B.13 reference approximation, not the nonlinear paper solution.

Independent finite parameters and analytic pressure datum are explicit below.
Only the near-axis region Lambda*X<=4.1 is exposed; no exterior is invented.
"""
from dataclasses import dataclass, asdict
from functools import cached_property, lru_cache
import math
from scipy.integrate import quad
from scipy.special import hyp0f1
from .natural_axis import real_gradient
from .profiles import LeadingProfile
from .velocity import leading_velocity_cartesian


@dataclass(frozen=True)
class PaperCoreReference:
    h: float = .005
    j: float = .02
    sigma: float = .1
    Lambda: float = 10.
    C: float = 2.
    pressure_scale: float = 1.

    def __post_init__(self):
        if not all(math.isfinite(v) for v in asdict(self).values()):
            raise ValueError('parameters must be finite')
        if not (0 < self.h < .01 and 0 < self.j <= .05 and self.sigma > 0
                and self.Lambda >= 1 and self.C > 1 and self.pressure_scale > 0):
            raise ValueError('parameters outside reference domain')

    def _check(self, X, eta):
        if not (math.isfinite(X) and math.isfinite(eta)
                and 0 <= self.Lambda*X <= 4.1 and abs(eta) <= 1):
            raise ValueError('reference only defined for 0<=Lambda*X<=4.1, |eta|<=1')

    @lru_cache(maxsize=512)
    def amplitude(self, eta):
        self._check(0, eta)
        phase = quad(lambda e: real_gradient(self.h, self.j, self.sigma, e),
                     0, eta, epsabs=1e-12, epsrel=1e-12)[0]
        return math.exp(self.Lambda*phase-math.log(self.C))

    def axial_coefficient(self, e):
        """B=-Z*/(2L) and its analytic eta derivative, B.1/B.13."""
        A, D = .5+self.h, .5-self.h
        d, L, Lp = 1-e*e, 1-2*self.h*e*e, -4*self.h*e
        Us = 4*e+self.j
        H, Hp = D*e+d*Us, D-2*e*Us+4*d
        p2 = self.pressure_scale**2
        p = -p2/(1+e*e)**2
        pp = 4*p2*e/(1+e*e)**3
        ppp = 4*p2*(1-5*e*e)/(1+e*e)**4
        Z = -A*(1-2*e*Us)*Us-4*H-d*pp+4*A*e*p
        Zp = -A*((-2*Us-8*e)*Us+4*(1-2*e*Us))-4*Hp+2*e*pp-d*ppp+4*A*(p+e*pp)
        return -Z/(2*L), -Zp/(2*L)+Z*Lp/(2*L*L)

    def F(self, X, eta):
        self._check(X, eta)
        H = (.5-self.h)*eta+(1-eta*eta)*(4*eta+self.j)
        chi = H*H/(H*H+self.sigma**2)
        # B.11: sum (-s/2)^n/[n!(n+1)!] = 0F1(;2;-s/2).
        return self.amplitude(eta)*float(hyp0f1(2, -.5*self.Lambda*X*chi))

    def U(self, X, eta):
        self._check(X, eta)
        return 4*eta+self.j+X*self.axial_coefficient(eta)[0]

    def dU(self, X, eta):
        self._check(X, eta)
        return 4+X*self.axial_coefficient(eta)[1]

    @cached_property
    def profile(self):
        return LeadingProfile(
            E=lambda X,e: math.sqrt(2*X)*self.F(X,e), U=self.U,
            dU_deta=self.dU, F=self.F,
            average_U=lambda X,e: self.U(X/2,e),
            average_dU_deta=lambda X,e: self.dU(X/2,e),
            name='paper-B13-core-reference-independent-parameters',
            provenance='Eqs B.1/B.3/B.11/B.13; autonomous pressure datum and finite parameters; no nonlinear correction')

    def velocity(self, x, y, z, t):
        if not math.isfinite(t) or not 0 <= t < 1:
            raise ValueError('time must satisfy 0<=t<1')
        return leading_velocity_cartesian(x,y,z,t,self.profile,h=self.h)

    def metadata(self):
        return dict(parameters=asdict(self), paper_exact=False,
                    status='near-axis reference approximation; nonlinear correction not solved',
                    domain='0<=t<1; Lambda*(x*x+y*y)/(2*q)<=4.1',
                    pressure_datum='-pressure_scale^2/(1+eta^2)^2, autonomous; not matched to exterior',
                    missing=['nonlinear profile solve','outer matching','higher-order and oscillatory corrections'])
