"""Prescribed Section 9.25 exit ODE with pressure and collar-width atoms.

State is integrated in s=y/h_b. Both the physical width and pressure tail
remain separate formal parameters; no materialized whole field is used.
This supplies the base joint ODE. Z tangents and downstream reshaping remain
explicitly uninstalled, so this is not a completed global matching claim.
"""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_axial_primitive import _sigma_mp


class PressureWidthExitBridge:
    def __init__(self,comparison,*,epsilon=None,steps=8):
        self.comparison=comparison;self.precision=comparison.precision
        self.pressure_order=comparison.pressure_order;self.width_order=comparison.width_order
        self.width=comparison.W;self.r=comparison.Ra
        self.epsilon_width_tied=epsilon is None
        with mp.workdps(self.precision):
            self.epsilon=self.width if epsilon is None else self.jet(mp.mpf(str(epsilon)))
        if steps<8:raise ValueError('At least 8 RK steps required')
        if epsilon is not None and not 0<mp.mpf(str(epsilon))<1:raise ValueError('Require 0<epsilon<1')
        self.steps=int(steps)
    def jet(self,value=0):
        return PressureWidthJet(value,pressure_order=self.pressure_order,width_order=self.width_order)
    def multiplier(self,s):
        x=mp.mpf(str(s))
        if x<=0:return self.jet(1)
        if x>=1:return self.epsilon
        return self.jet(_sigma_mp(1-x))+self.epsilon*_sigma_mp(x)
    @lru_cache(maxsize=32)
    def initial(self,Z):
        with mp.workdps(self.precision):
            z=mp.mpf(str(Z));a=self.comparison.core_snapshot(0,z)
            coeff=self.comparison.core_coefficients(z);r=self.r;Fa=a['F']
            axial=self.jet(0);swirl=self.jet(0)
            for i,ui in enumerate(coeff['Uz']):
                for j,uj in enumerate(coeff['Uz']):axial+=ui[0]*uj[0]*r**(i+j+1)/(i+j+1)
            for i,fi in enumerate(coeff['F']):
                for j,fj in enumerate(coeff['F']):swirl+=fi[0]*fj[0]*r**(i+j+2)/(i+j+2)
            m=a['moments']
            state=[self.jet(0),a['Uz'],m['theta']/(Fa*r*r),m['z']/r,
                m['theta_z']/(Fa*r*r),axial/r,swirl/(Fa*Fa*r*r),m['p']/(Fa*Fa*r)]
            # P0 comes from the original axis jet, avoiding cancellation of a
            # pressure integral against its reconstructed total.
            return dict(Fa=Fa,P0=coeff['P'][0][0],state=state,core=a)
    def rhs_y(self,s,state,Z):
        with mp.workdps(self.precision):
            initial=self.initial(Z);Fa=initial['Fa'];g,u=state[:2]
            bar=self.comparison.evaluate(s,Z);y=mp.mpf(str(s))*self.width
            chi=self.multiplier(s);relative_bar=(bar['F']/Fa).log()
            du=-chi*(bar['R']/2).sqrt()*(g-relative_bar).exp()*bar['I_z']
            eg=g.exp();ey=y.exp();e2y=(2*y).exp()
            return [-chi*bar['D']/2,du,2*e2y*eg,ey*u,2*e2y*eg*u,
                ey*u*u,e2y*eg*eg,ey*eg*eg]
    def rhs_s(self,s,state,Z):return [self.width*v for v in self.rhs_y(s,state,Z)]
    @lru_cache(maxsize=64)
    def evaluate(self,s,Z):
        with mp.workdps(self.precision):
            s=mp.mpf(str(s));z=mp.mpf(str(Z))
            if not 0<=s<=2 or abs(z)>=1:raise ValueError('Require 0<=s<=2, |Z|<1')
            init=self.initial(z);state=list(init['state']);h=s/self.steps
            def shifted(a,k,factor):return [v+factor*h*w for v,w in zip(a,k)]
            if s:
                for n in range(self.steps):
                    t=n*h;k1=self.rhs_s(t,state,z)
                    k2=self.rhs_s(t+h/2,shifted(state,k1,mp.mpf('.5')),z)
                    k3=self.rhs_s(t+h/2,shifted(state,k2,mp.mpf('.5')),z)
                    k4=self.rhs_s(t+h,shifted(state,k3,1),z)
                    state=[v+h*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(state,k1,k2,k3,k4)]
            g,u,theta,mz,mixed,axial,swirl,p=state
            Fa=init['Fa'];r=self.r;y=s*self.width;R=r*y.exp()
            moments=dict(theta=Fa*r*r*theta,z=r*mz,theta_z=Fa*r*r*mixed,
                z_theta=r*axial-Fa*Fa*r*r*swirl,p=Fa*Fa*r*p)
            slopes=self.rhs_y(s,state,z)
            return dict(s=s,y=y,R=R,F=Fa*g.exp(),Uz=u,log_F_over_Fa=g,
                F_R_over_F=slopes[0]/R,Uz_R=slopes[1]/R,moments=moments,
                P=init['P0']+moments['p'],chi=self.multiplier(s),steps=self.steps,
                pressure_order=self.pressure_order,width_order=self.width_order,
                Z_moment_derivatives_available=False,cone_certified=False,
                source_collar_constants_certified=False,outer_matching_complete=False)
    def metadata(self):
        return dict(source_equations='Section 9.25-9.26 joint base exit ODE',
            independent_variable='s=y/h_b',epsilon_width_tied_as_current_source=self.epsilon_width_tied,
            pressure_parameter_order=self.pressure_order,width_parameter_order=self.width_order,
            pressure_parameter_remainder_enclosed=False,width_parameter_remainder_enclosed=False,
            ODE_error_enclosed=False,Z_tangents_installed=False,reshape_installed=False,
            finite_energy_certified=False,cone_certified=False,global_field_installed=False)
