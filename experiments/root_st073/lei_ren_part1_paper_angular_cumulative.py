"""Actual forward angular/energy moments through the source pulse.

Normalized forward equations avoid exponentiating the huge absolute radius.
They use the candidate's actual unflattened swirl, not terminal moment targets.
"""
from decimal import Decimal
import json
from pathlib import Path
import mpmath as mp
import numpy as np


class AngularCumulative:
    def __init__(self,profile,*,order=128):
        self.profile=profile; self.schedule=profile.schedule
        self.precision=profile.precision; self.order=order

    def normalized(self,logR):
        p=self.profile; s=self.schedule
        y=p.offset(logR,s.logRref)
        if y>s.y_v:
            raise ValueError('Angular cumulative implementation stops at Rv')
        with mp.workdps(self.precision):
            # reference Utheta=constant*R^.1
            state=[mp.mpf(5)/8,mp.mpf(5)/6]
            if y<=0:return dict(theta=state[0],swirl_energy=state[1])
            nodes,weights=np.polynomial.legendre.leggauss(self.order)
            for name in ('slope_transition_ref','axial_turnoff',
                         'slope_transition_mu','power_buffer','pulse_reserved'):
                left,right=s._make_stage_bounds()[name]
                if y<=left:break
                end=min(y,right)
                L=mp.mpf(str(p.offset(end,left)))
                startlog=p.log_at(s.logRref,left); endlog=p.log_at(s.logRref,end)
                if name not in ('slope_transition_ref','slope_transition_mu'):
                    slope=mp.mpf('-.5')
                    if name in ('power_buffer','pulse_reserved'):
                        slope-=mp.mpf(str(s.mu))
                    for i,(a,b) in enumerate(((mp.mpf('1.5'),1),(1,2))):
                        rate=a+b*slope
                        state[i]=state[i]*mp.exp(-rate*L)+(-mp.expm1(-rate*L)/rate if rate else L)
                else:
                    E0=mp.mpf(str(s.at_log_radius(startlog,0)['log_angular_amplitude']))
                    E1=mp.mpf(str(s.at_log_radius(endlog,0)['log_angular_amplitude']))
                    for i,(a,b) in enumerate(((mp.mpf('1.5'),1),(1,2))):
                        value=state[i]*mp.exp(-a*L-b*(E1-E0))
                        for node,weight in zip(nodes,weights):
                            t=L*(mp.mpf(str(float(node)))+1)/2
                            Et=mp.mpf(str(s.at_log_radius(p.log_at(startlog,t),0)['log_angular_amplitude']))
                            value+=L/2*mp.mpf(str(float(weight)))*mp.exp(-a*(L-t)-b*(E1-Et))
                        state[i]=value
                if end==y:break
            return dict(theta=state[0],swirl_energy=state[1])

    def moments(self,logR,Z):
        with mp.workdps(self.precision):
            ratios=self.normalized(logR)
            R=mp.exp(mp.mpf(str(logR)))
            swirl=self.profile.values(logR,Z)['Utheta']
            return dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*swirl*ratios['theta'],
                        swirl_energy=R*swirl**2*ratios['swirl_energy'])


def run():
    from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
    from lei_ren_part1_paper_reference_moments import PaperReferenceMoments
    p=CorrectedSourceProfile(); engine=AngularCumulative(p)
    reference=PaperReferenceMoments(p,quadrature_order=128)
    with mp.workdps(p.precision):
        logR=p.log_at(p.schedule.logRref,'.5')
        actual=engine.moments(logR,.3); direct=reference(logR,.3)
        R=mp.exp(mp.mpf(str(logR)))
        # Uz=4Z at this point, so recover swirl energy from actual mixed moment.
        direct_energy=2*((mp.mpf('1.2')**2)*R-direct['z_theta'])
        agreement={'theta':mp.nstr(abs(actual['theta']/direct['theta']-1),30),
                   'swirl_energy':mp.nstr(abs(actual['swirl_energy']/direct_energy-1),30)}
        rows=[]; equations=[]
        for checkpoint,offset in ((p.schedule.logR_p,'0'),(p.schedule.logR_p,
             mp.nstr(5/mp.mpf(str(p.schedule.mu)),p.precision)),(p.schedule.logR_v,'0')):
            radius=p.log_at(checkpoint,offset)
            ratios=engine.normalized(radius)
            rows.append({'logR':str(radius),'ratios':{k:mp.nstr(v,p.precision) for k,v in ratios.items()}})
            h=mp.mpf('.001')
            neighbors={i:engine.normalized(p.log_at(radius,i*h)) for i in (-2,-1,1,2)
                       if p.offset(p.log_at(radius,i*h),p.schedule.logR_v)<=0}
            if len(neighbors)==4:
                slope=mp.mpf(str(p.schedule.at_log_radius(radius,0)['logarithmic_slope']))
                defects={}
                for key,a,b in (('theta',mp.mpf('1.5'),1),('swirl_energy',1,2)):
                    dy=(neighbors[-2][key]-8*neighbors[-1][key]+8*neighbors[1][key]-neighbors[2][key])/(12*h)
                    defects[key]=mp.nstr(abs(dy+(a+b*slope)*ratios[key]-1),30)
                equations.append({'logR':str(radius),'normalized_ODE_absolute_defects':defects})
        assert max(mp.mpf(v) for v in agreement.values())<mp.mpf('1e-10')
        assert max(mp.mpf(v) for row in equations for v in row['normalized_ODE_absolute_defects'].values())<mp.mpf('1e-7')
    return {'independent_initial_quadrature_agreement':agreement,'samples':rows,
            'normalized_forward_ODE_checks':equations,
            'scope':'Actual source swirl through Rv; no terminal moment substitution.',
            'later_outer_supported':False,'global_closure_certified':False}


if __name__=='__main__':
    receipt=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt['independent_initial_quadrature_agreement']))
