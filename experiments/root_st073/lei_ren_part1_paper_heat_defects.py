"""Logarithmic heat replacement defects (7.13),(7.17).

Returns the first Taylor term and an analytic relative remainder bound,
not a claim that the exact heat defects vanish. Quadrature error is separate.
H(x)=E[(1+x V)^(-h)], V~Gamma(1+h), gives
0 <= a1*x-(1-H(x)) <= a2*x*x/2.
Astronomical source radii never need to be exponentiated.
"""
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.integrate import solve_ivp
from scipy.special import gamma
from lei_ren_part1_paper_outer import PaperOuterSchedule, _sigma, _sigma_prime, _flat_edge
from lei_ren_part1_paper_waiting_length import (
    incoming_angular_ratio, solve_waiting_length, apply_waiting_root)


def preheat_angular_defect(schedule, incoming, Z, *, rtol=1e-11, order=256):
    """r_pre by a difference ODE, not subtraction of rounded tiny moments.

    Before flattening the normalized angular moment is independent of Z.
    Afterwards its difference contracts exactly by mu**(30*(1-mu)).
    Any discarded early-stage contraction is inherited from incoming.
    """
    z=float(Z)
    if not math.isfinite(z) or abs(z)>1:
        raise ValueError('Z must lie in [-1,1]')
    if z==0:
        return {'Z':z,'exact_zero':True,'log_r_pre':None}
    a=float(1-schedule.mu); duration=float(schedule.Tf)
    base_log=-math.log(2.); difference_log=math.log1p(z*z)
    def rhs(t,state):
        derivative=_sigma_prime(t/duration)/duration
        base_rate=a+base_log*derivative
        return [1-base_rate*state[0]]
    result=solve_ivp(rhs,(0.,duration),
                     [float(incoming['incoming_flattening_X'])],dense_output=True,
                     rtol=rtol,atol=rtol/100,max_step=duration/200)
    if not result.success:
        raise ArithmeticError('Flattening difference ODE failed positive defect')
    # Positive variation-of-constants integral preserves the very small
    # endpoint difference; an absolute ODE tolerance on Delta would erase it.
    nodes,weights=leggauss(order)
    deficit=0.
    z_log=base_log+difference_log
    for node,weight in zip(nodes,weights):
        t=duration*(float(node)+1)/2
        exponent=-a*(duration-t)-z_log*(1-_sigma(t/duration))
        deficit+=float(weight)/2*math.exp(exponent)*difference_log\
            *_sigma_prime(t/duration)*float(result.sol(t)[0])
    if deficit<=0:
        raise ArithmeticError('Flattening integral lost positive defect')
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        log_defect=(Decimal(str(deficit)).ln()
                    +30*(1-schedule.mu)*schedule.mu.ln())
    return {'Z':z,'exact_zero':False,'log_r_pre':str(log_defect),
            'flattening_X_difference':-deficit,'ode_rtol':rtol,
            'positive_difference_integral_order':order,
            'inherited_discarded_contraction_log_factors':
                incoming['discarded_constant_stage_contraction_log_factors'],
            'scope':'Pre-heat defect relative to the matched Z=0 target; waiting-root numerical error remains separate.'}


def heat_defects(schedule, Z, *, order=96):
    z=float(Z)
    if not math.isfinite(z) or abs(z)>1:
        raise ValueError('Z must lie in [-1,1]')
    if order<16:
        raise ValueError('quadrature order must be at least 16')
    if abs(z)==1:
        return {'Z':z,'exact_zero':True,'log_r_H':None,'log_s_H':None}
    nodes, weights=leggauss(order)
    h=float(schedule.delta/2); eps=float(schedule.epsilon)
    # g=sigma*f is exactly one beyond t=3; integrate that infinite tail
    # analytically, including the 1/h factor in the angular moment.
    angular_short=0.; pressure_short=0.
    for left in range(3):
        for node,weight in zip(nodes,weights):
            t=left+(float(node)+1)/2
            sigma=_sigma(t); g=sigma*(1-eps*_flat_edge((3-t)/2))
            pre=(1-sigma)*(1-eps)+g
            angular_short+=float(weight)/2*math.exp(-h*t)*g
            pressure_short+=float(weight)/2*math.exp(-(2+2*h)*t)*pre*g
    pressure_integral=pressure_short+math.exp(-3*(2+2*h))/(2+2*h)
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        hd=schedule.delta/2; d=Decimal(str(1-z*z))
        angular_scaled=(-3*hd).exp()+hd*Decimal(str(angular_short))
        logE=schedule._log_A(schedule.y_rel)-Decimal(2).ln()
        log_r=(schedule._log_c_inf+(2*d*(1+hd)).ln()
               -hd*schedule.logR_tail-Decimal('1.5')*schedule.logR_rel
               -logE+angular_scaled.ln())
        a1=hd*(1+hd)
        log_s=(2*schedule._log_c_inf+(2*d*a1).ln()
               -(2+2*hd)*schedule.logR_tail-2*logE
               +Decimal(str(pressure_integral)).ln())
        # The bounds below concern truncation of H, not finite quadrature.
        # For pressure use g/pre <= 1/(1-epsilon) and 0<=1-H<=a1*x.
        angular_bound=(1+hd)*(2+hd)/2
        pressure_bound=angular_bound+a1/(2*(1-schedule.epsilon))
        log_xi=(2*d).ln()-schedule.logR_tail
        return {'Z':z,'exact_zero':False,
                'log_r_H':str(log_r),'log_s_H':str(log_s),
                'log_relative_Taylor_bound_r_H':str(log_xi+angular_bound.ln()),
                'log_relative_Taylor_bound_s_H':str(log_xi+pressure_bound.ln()),
                'scaled_angular_integral':str(angular_scaled),
                'pressure_integral':pressure_integral,'quadrature_order':order,
                'scope':'Positive leading heat defects with analytic relative truncation bounds; r_pre and coefficient solve are not included.'}


def run():
    base=PaperOuterSchedule(logPstar=14,logRref=10,delta='1e-32',Md='.5',
                           c_mu='.001',c_delta='.001',c_epsilon='.01')
    incoming=incoming_angular_ratio(base)
    matched=apply_waiting_root(base,solve_waiting_length(base,incoming['incoming_X']))
    rows=[heat_defects(matched,z) for z in (0.,.5,-.5,1.,-1.)]
    pre_rows=[preheat_angular_defect(matched,incoming,z) for z in (0.,.5,-.5,1.,-1.)]
    pre_coarse=preheat_angular_defect(matched,incoming,.5,rtol=1e-9,order=128)
    coarse=heat_defects(matched,0.,order=48)
    with localcontext() as ctx:
        ctx.prec=matched.decimal_precision
        differences={key:str(abs(Decimal(str(rows[0][key]))-Decimal(str(coarse[key]))))
                     for key in ('scaled_angular_integral','pressure_integral')}
    # Independently integrate the Gamma expectation at representable xi.
    # expm1 avoids the very cancellation that the production path removes.
    heat_checks=[]
    for h in (.001,.01):
        a1=h*(1+h); a2=h*(1+h)**2*(2+h)
        for xi in (.001,.01,.1):
            loss,error=quad(lambda v: math.exp(-v)*v**h
                *(-math.expm1(-h*math.log1p(xi*v)))/gamma(1+h),
                0.,np.inf,epsabs=1e-13,epsrel=1e-11)
            deficit=a1*xi-loss; bound=a2*xi*xi/2
            heat_checks.append({'h':h,'xi':xi,'actual_1_minus_H':loss,
                'linear_overestimate':deficit,'Taylor_bound':bound,
                'quadrature_error':error,
                'bound_passed':-10*error<=deficit<=bound+10*error})
    if not all(row['bound_passed'] for row in heat_checks):
        raise ArithmeticError('Independent heat expectation failed Taylor bound')
    report={'source':'https://arxiv.org/html/2609.35406v1',
            'equations':'7.13,7.17','schedule':matched.metadata(),
            'heat_defects':rows,'quadrature_refinement':differences,
            'preheat_angular_defects':pre_rows,
            'preheat_difference_ODE_refinement':abs(pre_rows[1]['flattening_X_difference']
                -pre_coarse['flattening_X_difference']),
            'independent_Gamma_expectation_checks':heat_checks,
            'parity_exact_in_computed_logs':rows[1]['log_r_H']==rows[2]['log_r_H']
                and rows[1]['log_s_H']==rows[2]['log_s_H'],
            'actual_angular_correction_solved':False,'full_outer_closed':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'quadrature_refinement':differences,'central_defects':rows[0]}))
    return report


if __name__=='__main__': run()
