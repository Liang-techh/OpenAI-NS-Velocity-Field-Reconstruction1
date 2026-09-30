"""Source angular bump solve with arbitrary-exponent coefficients.

The inputs are candidate-derived numerical r_pre and bounded Taylor heat
defects. The finite-dimensional equations are solved on their small branch;
this does not certify exact heat moments, derivatives, cone or global core.
"""
from decimal import Decimal, localcontext
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from lei_ren_part1_paper_outer import PaperOuterSchedule
from lei_ren_part1_paper_waiting_length import (
    incoming_angular_ratio, solve_waiting_length, apply_waiting_root)
from lei_ren_part1_paper_heat_defects import heat_defects, preheat_angular_defect
from lei_ren_part1_paper_outer_closure import paper_bump_weights, _paper_raw_bump


class AngularCorrection:
    def __init__(self, schedule, incoming, *, precision=120, order=128):
        if precision<64:
            raise ValueError('At least 64 decimal digits are required')
        self.schedule=schedule; self.incoming=incoming
        self.precision=precision; self.order=order
        if float(schedule.mu)==0:
            raise ValueError('mu underflows the finite-quadrature bump weights; arbitrary precision weights are required')
        self.weights=paper_bump_weights(mu=float(schedule.mu),quadrature_order=order)

    def coefficients(self,Z):
        heat=heat_defects(self.schedule,Z,order=self.order)
        pre=preheat_angular_defect(self.schedule,self.incoming,Z)
        with mp.workdps(self.precision):
            def value(row,key):
                log=row.get(key)
                return mp.mpf(0) if log is None else mp.exp(mp.mpf(log))
            r_pre=value(pre,'log_r_pre'); r_H=value(heat,'log_r_H')
            r=r_pre+r_H; s=value(heat,'log_s_H')
            if not r>0:
                raise ArithmeticError('Actual angular target must remain positive')
            mu=mp.mpf(str(self.schedule.mu)); A=mp.mpf(str(self.weights.A_mu))
            B=mp.mpf(str(self.weights.B_mu)); D=mp.mpf(str(self.weights.D_mu))
            alpha=mp.exp(-2*(1-mu)); beta=mp.exp(-2*(1+2*mu))
            rho=mp.exp(1-mu)/A
            sigma=mp.exp(-3*(1+2*mu))*s/(B*r); q=D/(2*B)
            # Divide d_j by r before solving. This leaves finite unknowns
            # even when log(r) is -2.7e28. Use the cancellation-free root.
            c2=q*r*(1+beta*alpha**2)
            c1=1-beta*alpha-2*q*beta*alpha*rho*r
            c0=beta*rho+q*beta*rho**2*r-sigma
            discriminant=c1*c1-4*c2*c0
            if c1<=0 or discriminant<=0:
                raise ArithmeticError('Source small-branch assumptions failed')
            e1=-2*c0/(c1+mp.sqrt(discriminant)); e2=rho-alpha*e1
            d1=r*e1; d2=r*e2
            fmt=lambda x:mp.nstr(x,self.precision)
            signed=lambda x:{'sign':int(mp.sign(x)),
                'log_abs':fmt(mp.log(abs(x))) if x else None,
                'arbitrary_exponent_value':fmt(x)}
            angular=alpha*e1+e2-rho
            pressure=e1+beta*e2+q*r*(e1*e1+beta*e2*e2)-sigma
            negative=max(-min(d1,mp.mpf(0)),-min(d2,mp.mpf(0)))
            if negative*mp.mpf(str(self.weights.beta_inf))>=1:
                raise ArithmeticError('Bump multiplier lost positivity')
            return {'Z':float(Z),'r_pre':pre,'heat':heat,'log_r':fmt(mp.log(r)),
                'log_s':fmt(mp.log(s)) if s else None,
                'd1':signed(d1),'d2':signed(d2),
                'd1_over_r':fmt(e1),'d2_over_r':fmt(e2),
                'scaled_algebraic_residuals':[fmt(angular),fmt(pressure)],
                'bump_multiplier_positive':True,
                'log_maximum_negative_multiplier_deficit':fmt(mp.log(negative*self.weights.beta_inf)) if negative else None,
                'decimal_precision':self.precision,'bump_quadrature_order':self.order,
                'input_status':'Numerical preheat difference and bounded first-Taylor heat inputs; not exact heat closure',
                'full_outer_closed':False}

    def corrected_at_log_radius(self,log_radius,Z,*,coefficients=None):
        """Attach the actual multiplicative bump to the source profile.

        Keep the log1p term separately: summing it into an O(1e28) log
        amplitude would discard it even with arbitrary-exponent arithmetic.
        """
        base=self.schedule.at_log_radius(log_radius,Z)
        with localcontext() as ctx:
            ctx.prec=self.schedule.decimal_precision
            offset=Decimal(str(log_radius))-self.schedule.logR_rel
        if offset<Decimal('-3.15') or offset>Decimal('-.85'):
            correction={'sign':0,'log_abs':None,'log1p_relative_correction':'0.0'}
        else:
            correction=self.relative_bump(float(offset),Z,coefficients=coefficients)
        return {'base_profile':base,'relative_angular_bump':correction,
                'log_angular_amplitude_terms':[str(base['log_angular_amplitude']),
                    correction['log1p_relative_correction']],
                'scope':'Corrected angular profile representation only; corrected pressure and axial pulse still pending.'}

    def relative_bump(self,t,Z,*,coefficients=None):
        """h(t,Z) for t=log(R/Rrel), in signed-log representation."""
        row=self.coefficients(Z) if coefficients is None else coefficients
        if row['Z']!=float(Z):
            raise ValueError('Coefficient receipt uses a different Z')
        with mp.workdps(self.precision):
            result=mp.mpf(0)
            for name,center in (('d1',-3.),('d2',-1.)):
                beta=float(_paper_raw_bump((float(t)-center)/self.weights.ell))\
                     /(self.weights.ell*self.weights.normalization_integral)
                data=row[name]
                if beta and data['sign']:
                    result+=data['sign']*mp.exp(mp.mpf(data['log_abs']))*beta
            return {'sign':int(mp.sign(result)),
                'log_abs':mp.nstr(mp.log(abs(result)),self.precision) if result else None,
                'log1p_relative_correction':mp.nstr(mp.log1p(result),self.precision)}

    def independent_bump_replay(self,row,*,order=192):
        """Actual bump quadrature after dividing by r; no huge radii."""
        nodes,weights=leggauss(order)
        with mp.workdps(self.precision):
            mu=float(self.schedule.mu)
            e1=float(row['d1_over_r']); e2=float(row['d2_over_r'])
            angular=0.; pressure=0.
            for center,e in ((-3.,e1),(-1.,e2)):
                for node,weight in zip(nodes,weights):
                    s=self.weights.ell*float(node); t=center+s
                    beta=float(_paper_raw_bump(float(node)))\
                        /(self.weights.ell*self.weights.normalization_integral)
                    measure=self.weights.ell*float(weight)
                    angular+=measure*np.exp((1-mu)*t)*e*beta
                    pressure+=measure*np.exp(-(1+2*mu)*t)*e*beta
            # The omitted quadratic term divided by r is O(r), reported
            # explicitly. Float replay cannot resolve the tiny pressure target.
            r=mp.exp(mp.mpf(row['log_r']))
            target_ratio=mp.exp(mp.mpf(row['log_s'])-mp.mpf(row['log_r']))\
                if row['log_s'] is not None else mp.mpf(0)
            return {'angular_increment_over_r':angular,
                'angular_scaled_difference_from_one':abs(angular-1),
                'pressure_linear_increment_over_r':pressure,
                'log_pressure_target_over_r':mp.nstr(mp.log(target_ratio),30) if target_ratio else None,
                'log_quadratic_term_scale_over_r':row['log_r'],
                'quadrature_order':order,
                'scope':'Independent finite-precision bump replay; does not resolve cancellation at the tiny pressure-target scale.'}


def run():
    base=PaperOuterSchedule(logPstar=14,logRref=10,delta='1e-32',Md='.5',
                           c_mu='.001',c_delta='.001',c_epsilon='.01')
    incoming=incoming_angular_ratio(base)
    matched=apply_waiting_root(base,solve_waiting_length(base,incoming['incoming_X']))
    correction=AngularCorrection(matched,incoming)
    rows=[correction.coefficients(z) for z in (0.,.5,-.5,1.)]
    replay=[correction.independent_bump_replay(row) for row in rows]
    for row,check in zip(rows,replay):
        if check['angular_scaled_difference_from_one']>1e-8:
            raise ArithmeticError('Actual bump replay failed angular identity')
        row['bump_center_values']=[correction.relative_bump(t,row['Z'],coefficients=row)
                                   for t in (-3.,-1.,0.)]
        with localcontext() as ctx:
            ctx.prec=matched.decimal_precision
            log_radius=matched.logR_rel-Decimal(3)
        represented=correction.corrected_at_log_radius(log_radius,row['Z'],coefficients=row)
        row['corrected_center_log_amplitude_terms']=represented['log_angular_amplitude_terms']
        if represented['relative_angular_bump']['sign']!=-1:
            raise ArithmeticError('Corrected profile lost its first signed bump')
    report={'source':'https://arxiv.org/html/2609.35406v1','equations':'7.20-7.23',
            'schedule':matched.metadata(),'coefficients':rows,'independent_replay':replay,
            'scope':'Small-branch coefficients for candidate-derived approximate inputs, preserving arbitrary exponents. Exact heat/pressure cancellation, derivatives, cone and regular core remain open.',
            'full_outer_closed':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'angular_replay_max':max(r['angular_scaled_difference_from_one'] for r in replay),
                      'central_d1':rows[0]['d1'],'off_axis_d1':rows[1]['d1']}))
    return report


if __name__=='__main__':run()
