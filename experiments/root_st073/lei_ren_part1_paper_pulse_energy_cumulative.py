"""Actual cumulative axial quadratic moment on the source pulse interval."""
from decimal import Decimal
import mpmath as mp
import numpy as np
from lei_ren_part1_paper_axial_correction import from_signed_log
from lei_ren_part1_paper_outer_closure import _paper_raw_bump
from lei_ren_part1_paper_axial_primitive import source_pulse_value_mp


class PulseEnergyCumulative:
    def __init__(self,profile):
        self.profile=profile; self.precision=profile.precision

    def moment(self,logR,Z):
        p=self.profile; s=p.schedule
        start=p.offset(logR,s.logR_p); end=p.offset(logR,s.logR_v)
        if start<0 or end>0:
            raise ValueError('Axial quadratic moment currently supports Rp<=R<=Rv')
        with mp.workdps(self.precision):
            data=p.coefficients(float(Z)); axial=data['axial']; incoming=data['incoming']
            mu=mp.mpf(str(s.mu)); xi=mu*mp.mpf(str(start)); order=axial['quadrature_order']
            nodes,weights=np.polynomial.legendre.leggauss(order)
            logRmp=mp.mpf(str(logR)); logE=mp.mpf(str(s.at_log_radius(logR,Z)['log_angular_amplitude']))
            # The incoming axial energy is constant after the first turnoff.
            initial=mp.exp(mp.mpf(str(s.logRref)))*mp.mpf(incoming['dimensionless_integrals']['I_uz2'])
            integral=mp.mpf(0); upper=min(xi,mp.mpf(11))
            cuts=[mp.mpf(0)]+[v for v in (mp.mpf('.02'),mp.mpf(10)) if v<upper]+[upper]
            for left,right in zip(cuts,cuts[1:]):
                half=(right-left)/2; midpoint=(right+left)/2
                for n,w in zip(nodes,weights):
                    v=midpoint+half*mp.mpf(str(float(n)))
                    gp=source_pulse_value_mp(mp.nstr(v,self.precision),precision=self.precision)
                    integral+=half*mp.mpf(str(float(w)))*mp.exp(-2*v)*gp*gp
            norm_pulse=mp.mpf(axial['a_p'])**2/mu*mp.exp(2*xi)*integral
            # End bumps have disjoint support from gp, so no cross term.
            norm_bumps=mp.mpf(0); ell=mp.mpf('.15'); endmp=mp.mpf(str(end))
            raw=[mp.mpf(str(float(_paper_raw_bump(float(n))))) for n in nodes]
            normalization=sum(mp.mpf(str(float(w)))*v for w,v in zip(weights,raw))
            for center,coefficient in zip((-3,-1),axial['c']):
                left=mp.mpf(center)-ell; right=min(endmp,mp.mpf(center)+ell)
                if right<=left:continue
                half=(right-left)/2; midpoint=(right+left)/2
                c=from_signed_log(coefficient)
                for index,(n,w) in enumerate(zip(nodes,weights)):
                    t=midpoint+half*mp.mpf(str(float(n)))
                    bump=raw[index] if right==mp.mpf(center)+ell else mp.mpf(
                        str(float(_paper_raw_bump(float((t-center)/ell)))))
                    beta=bump/(ell*normalization)
                    norm_bumps+=half*mp.mpf(str(float(w)))*mp.exp(-2*mu*(t-endmp))*c*c*beta*beta
            scale=mp.exp(logRmp+2*logE)
            return initial+scale*(norm_pulse+norm_bumps)
