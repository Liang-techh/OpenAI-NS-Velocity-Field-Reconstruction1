"""Continuous MP angular amplitude on the incoming cutoff support.

This provider has the exact source primitive definition rather than the
legacy binary64 Gauss primitive. Point quadrature is nominal; independent
interval bounds are provided by continuous_incoming_angular_enclosure.
It is not yet a replacement for the complete exterior angular field.
"""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse


class ContinuousIncomingAngular:
    def __init__(self,logPstar='14',Md='.5',*,precision=100):
        if precision<50:raise ValueError('At least 50 working digits required')
        self.precision=int(precision)
        with mp.workdps(self.precision):
            self.logPstar=mp.mpf(logPstar);self.Md=mp.mpf(Md)
            if self.Md<=0:raise ValueError('Positive Md required')
            self.cutoff_y=mp.exp(self.Md)

    @lru_cache(maxsize=512)
    def _primitive(self,text):
        with mp.workdps(self.precision):
            y=mp.mpf(text)
            if y<=0:return mp.mpf(0)
            if y>=1:return y-mp.mpf('.5')
            # Symmetry gives J(y)=y-.5+J(1-y), avoiding full-minus-tail loss.
            if y>mp.mpf('.5'):return y-mp.mpf('.5')+self._primitive(mp.nstr(1-y,self.precision))
            return mp.quad(lambda t:ContinuousAxialPulse.sigma_pair(t)[0],[0,y])

    def J(self,y):
        with mp.workdps(self.precision):return self._primitive(mp.nstr(mp.mpf(y),self.precision))

    def log_amplitude(self,y):
        with mp.workdps(self.precision):
            y=mp.mpf(y)
            if y>self.cutoff_y:raise ValueError('Incoming-only angular provider used beyond cutoff support')
            return self.logPstar+mp.mpf('.1')*y-mp.mpf('.6')*self.J(y)

    def values(self,y,Z):
        with mp.workdps(self.precision):
            y=mp.mpf(y);z=mp.mpf(Z)
            switch,derivative=ContinuousAxialPulse.sigma_pair(y)
            logu=self.log_amplitude(y)-mp.log(1+z*z)
            return dict(Utheta=mp.exp(logu),log_Utheta=logu,
                Utheta_Z=-2*z/(1+z*z)*mp.exp(logu),
                log_amplitude_y=mp.mpf('.1')-mp.mpf('.6')*switch,
                log_amplitude_yy=-mp.mpf('.6')*derivative)

    def mixed_factor(self,*,order=48):
        """Nominal coefficient of Z/(1+Z^2) in the full I_theta_z row."""
        from lei_ren_part1_paper_continuous_incoming import ContinuousIncomingAxial
        with mp.workdps(self.precision):
            nodes,weights=mp.gauss_quadrature(order,'legendre')
            unit=sum(weight*mp.exp(mp.mpf('1.6')*(node+1)/2-
                mp.mpf('.6')*self.J((node+1)/2)) for node,weight in zip(nodes,weights))/2
            incoming=ContinuousIncomingAxial(self.Md,precision=self.precision)
            transition=incoming.weighted_integral()-mp.expm1(1)
            return 4*mp.exp(self.logPstar)*(mp.mpf(1)/mp.mpf('1.6')+unit+mp.exp(mp.mpf('.3'))*transition)
