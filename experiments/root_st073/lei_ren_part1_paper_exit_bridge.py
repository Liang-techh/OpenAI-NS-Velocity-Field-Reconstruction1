"""Actual initial Section 9.25 exit, integrating five cumulative moments.

One joint ODE carries the velocity and its moment increments. No moment
boundary targets are imposed. Collar constants are numerical inputs,
not certified Section 9 contraction constants.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_primitive import _sigma_mp


class ExitBridge:
    def __init__(self,comparison,*,epsilon=None,steps=32):
        self.comparison=comparison; self.core=comparison.core
        self.precision=comparison.precision
        with mp.workdps(self.precision):
            self.hb=mp.mpf(str(comparison.h_b))
            if epsilon is None:
                axis=comparison.axis
                Gmax=max(axis.G(-1),axis.G(1))
                self.log_epsilon=-axis.Lambda*(2*Gmax+1)
                self.epsilon=mp.exp(self.log_epsilon)
                self.epsilon_policy='Numerical real-domain amplitude margin; Section 9 constants not certified'
            else:
                self.epsilon=mp.mpf(str(epsilon)); self.log_epsilon=mp.log(self.epsilon)
                self.epsilon_policy='Caller supplied'
            self.r=4/self.core.Lambda
        if not 0<self.epsilon<1 or steps<8:
            raise ValueError('Require 0<epsilon<1 and at least 8 integration steps')
        self.steps=int(steps)

    def multiplier(self,y):
        x=y/self.hb
        if x<=0:return mp.mpf(1)
        if x>=1:return self.epsilon
        return _sigma_mp(1-x)+self.epsilon*_sigma_mp(x)

    @lru_cache(maxsize=16)
    def initial(self,Z):
        with mp.workdps(self.precision):
            a=self.core.evaluate(self.r,Z); coeff=self.core.coefficients(Z)
            Fa=a['F']; r=self.r; f=coeff['F']; u=coeff['Uz']
            axial=sum(ui[0]*uj[0]*r**(i+j+1)/(i+j+1)
                      for i,ui in enumerate(u) for j,uj in enumerate(u))
            swirl=sum(fi[0]*fj[0]*r**(i+j+2)/(i+j+2)
                      for i,fi in enumerate(f) for j,fj in enumerate(f))
            m=a['moments']
            state=[mp.mpf(0),a['Uz'],m['theta']/(Fa*r*r),m['z']/r,
                   m['theta_z']/(Fa*r*r),axial/r,swirl/(Fa*Fa*r*r),
                   m['p']/(Fa*Fa*r)]
            return {'Fa':Fa,'P0':a['P']-m['p'],'state':state,'core':a}

    def rhs(self,y,state,Z):
        with mp.workdps(self.precision):
            initial=self.initial(Z); Fa=initial['Fa']; g,u=state[:2]
            bar=self.comparison.evaluate(y,Z); R=self.r*mp.exp(y)
            chi=self.multiplier(y)
            relative_bar=mp.log(bar['F']/Fa)
            Iz=bar.get('I_z',bar['E']*bar['F'])
            du=-chi*mp.sqrt(R/2)*mp.exp(g-relative_bar)*Iz
            expg=mp.exp(g)
            return [-chi*bar['D']/2,du,2*mp.exp(2*y)*expg,
                    mp.exp(y)*u,2*mp.exp(2*y)*expg*u,mp.exp(y)*u*u,
                    mp.exp(2*y)*expg*expg,mp.exp(y)*expg*expg]

    def evaluate(self,y,Z):
        with mp.workdps(self.precision):
            y=mp.mpf(str(y))
            if not 0<=y<=2*self.hb:
                raise ValueError('Initial exit currently supports 0<=y<=2hb')
            init=self.initial(Z); state=list(init['state'])
            h=y/self.steps
            def shifted(a,k,factor):return [v+factor*h*w for v,w in zip(a,k)]
            for n in range(self.steps):
                t=n*h; k1=self.rhs(t,state,Z)
                k2=self.rhs(t+h/2,shifted(state,k1,mp.mpf('.5')),Z)
                k3=self.rhs(t+h/2,shifted(state,k2,mp.mpf('.5')),Z)
                k4=self.rhs(t+h,shifted(state,k3,1),Z)
                state=[v+h*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(state,k1,k2,k3,k4)]
            g,u,theta,z,mixed,axial,swirl,p=state
            Fa=init['Fa']; r=self.r; R=r*mp.exp(y)
            moments={'theta':Fa*r*r*theta,'z':r*z,'theta_z':Fa*r*r*mixed,
                     'z_theta':r*axial-Fa*Fa*r*r*swirl,'p':Fa*Fa*r*p}
            slopes=self.rhs(y,state,Z)
            return {'R':R,'y':y,'F':Fa*mp.exp(g),'Uz':u,'log_F_over_Fa':g,
                    'F_R_over_F':slopes[0]/R,'Uz_R':slopes[1]/R,
                    'moments':moments,'P':init['P0']+moments['p'],
                    'chi':self.multiplier(y),'steps':self.steps,
                    'Z_moment_derivatives_available':False,'cone_certified':False,
                    'source_collar_constants_certified':False,'outer_matching_complete':False}


def run():
    from lei_ren_part1_paper_exit_comparison import build_comparison
    comparison=build_comparison()
    bridges=[ExitBridge(comparison,steps=n) for n in (16,32,64)]
    rows=[]
    with mp.workdps(comparison.precision):
        def number(x):return mp.nstr(x,35)
        def signed_log(x):
            return {'sign':int(mp.sign(x)), 'log_abs':number(mp.log(abs(x))) if x else None}
        for Z in (comparison.axis.Z0,mp.mpf('.3')):
            initial=bridges[-1].initial(Z)
            start=bridges[-1].evaluate(0,Z)
            endpoints=[bridge.evaluate(2*comparison.hb,Z) for bridge in bridges]
            scales={'theta':initial['Fa']*bridges[-1].r**2,
                    'z':bridges[-1].r,
                    'theta_z':initial['Fa']*bridges[-1].r**2,
                    'z_theta':bridges[-1].r,
                    'p':initial['Fa']**2*bridges[-1].r}
            refinements=[]
            for coarse,fine in zip(endpoints,endpoints[1:]):
                refinements.append({
                    'steps':[coarse['steps'],fine['steps']],
                    'log_F_over_Fa_difference':number(abs(coarse['log_F_over_Fa']-fine['log_F_over_Fa'])),
                    'Uz_difference':number(abs(coarse['Uz']-fine['Uz'])),
                    'normalized_moment_differences':{k:number(abs(coarse['moments'][k]-fine['moments'][k])/scales[k]) for k in scales}})
            end=endpoints[-1]
            rows.append({'Z':number(Z),
                'start_field_relative_error':number(abs(start['F']/initial['Fa']-1)),
                'start_Uz_error':number(abs(start['Uz']-initial['core']['Uz'])),
                'start_normalized_moment_errors':{k:number(abs(start['moments'][k]-initial['core']['moments'][k])/scales[k]) for k in scales},
                'start_logarithmic_slope_defect':number(abs(bridges[-1].r*(start['F_R_over_F']-initial['core']['F_R']/initial['Fa']))),
                'terminal_log_F_over_Fa':number(end['log_F_over_Fa']),
                'terminal_Uz':number(end['Uz']),
                'moment_increments':{k:signed_log(end['moments'][k]-start['moments'][k]) for k in scales},
                'normalized_moment_increments':{k:number((end['moments'][k]-start['moments'][k])/scales[k]) for k in scales},
                'refinement':refinements})
    receipt={'source':'https://arxiv.org/html/2609.35406v1','equations':'9.25-9.26',
        'precision':comparison.precision,'radial_degree':18,'Lambda':'1e36',
        'h_b':number(comparison.hb),'log_epsilon':number(bridges[-1].log_epsilon),
        'epsilon_policy':bridges[-1].epsilon_policy,
        'comparison_transition_steps':comparison.transition_steps,
        'refinement_scope':'Actual joint ODE only; comparison discretization held fixed',
        'Z_moment_derivatives_available':False,'cone_certified':False,
        'source_collar_constants_certified':False,'outer_matching_complete':False,
        'rows':rows}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'rows':len(rows),'refinement':[row['refinement'] for row in rows]},indent=2))
    return receipt


if __name__=='__main__':run()
