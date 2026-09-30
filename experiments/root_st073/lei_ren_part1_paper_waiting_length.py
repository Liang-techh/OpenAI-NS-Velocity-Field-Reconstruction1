"""Normalized pre-heat waiting-length equation (7.9)--(7.10).

Incoming angular moments are obtained from the actual candidate ODE
X'+(3/2+logarithmic_slope)X=1, X(0)=5/8. Long constant stages are
propagated analytically. Decimal arithmetic retains the tiny heat deficit
that the waiting segment amplifies. This does not restore actual heat or
other terminal moments; the Section 7 angular/axial repairs remain required.
"""
from decimal import Decimal, localcontext
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from lei_ren_part1_paper_outer import PaperOuterSchedule, _sigma, _flat_edge


def incoming_angular_ratio(schedule,*,rtol=1e-11):
    omitted=[]
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        ctx.Emax=999999999; ctx.Emin=-999999999
        x=Decimal('0.625')
        def propagate(value,length,a):
            equilibrium=1/a; exponent=-a*length
            if exponent < Decimal('-1000000'):
                # A controlled contraction, not an attempted float exp.
                omitted.append(str(exponent))
                return equilibrium
            return equilibrium+(value-equilibrium)*exponent.exp()
        def transition(value,start,length):
            def equation(t,state):
                with localcontext() as inner:
                    inner.prec=schedule.decimal_precision
                    row=schedule.at_log_radius(schedule.logRref+start+Decimal(str(t)),0.)
                return [1-(1.5+float(row['logarithmic_slope']))*state[0]]
            solution=solve_ivp(equation,(0.,float(length)),[float(value)],
                               rtol=rtol,atol=rtol/100,max_step=float(length)/16)
            if not solution.success: raise ArithmeticError(solution.message)
            return Decimal(str(solution.y[0,-1]))
        x=transition(x,Decimal(0),Decimal(1))
        x=propagate(x,schedule.y_d-1,Decimal(1))
        x=transition(x,schedule.y_d,Decimal(1))
        x=propagate(x,schedule.y_v-schedule.y_w,1-schedule.mu)
        incoming_flattening_X=x
        x=transition(x,schedule.y_v,schedule.Tf)
        x=propagate(x,schedule.y_rel-schedule.y_f,1-schedule.mu)
        return {'incoming_X':x,'incoming_flattening_X':incoming_flattening_X,'ode_rtol':rtol,
                'discarded_constant_stage_contraction_log_factors':omitted,
                'discarded_stage_note':'Each discarded correction is (incoming-equilibrium)*exp(log_factor); no physical radius is truncated.'}


def target_coefficients(schedule,*,order=96):
    nodes,weights=leggauss(order)
    def integrate(a,b,function):
        return (b-a)/2*sum(float(w)*function(a+(b-a)*(float(t)+1)/2)
                           for t,w in zip(nodes,weights))
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        k=1-schedule.delta/2; a=1-schedule.mu
        # J(t)=integral sigma, evaluated on a short interval only.
        def J(t): return integrate(0.,t,lambda v:_sigma(v))
        Lrestore=integrate(0.,1.,lambda t:np.exp(float(k)*J(t)))
        Ldrop=integrate(0.,1.,lambda t:np.exp(float(a)*(t-J(t))))
        K=sum(integrate(float(left),float(left+1),lambda t:
              np.exp(float(k)*t)*(1-_sigma(t)+_sigma(t)*_flat_edge((3-t)/2)))
              for left in range(3))
        eq=1/k; er=(k/2).exp(); ed=(a/2).exp()
        constant=ed*(er*eq-Decimal(str(Lrestore))-schedule.Ts)-Decimal(str(Ldrop))
        exponential=ed*er*schedule.epsilon/(1-schedule.epsilon)*(eq+Decimal(str(K)))
        return {'constant':constant,'exponential_coefficient':exponential,
                'growth_rate':k,'collar_deficit_integral_K':K,
                'restore_transition_integral':Lrestore,'drop_transition_integral':Ldrop,
                'quadrature_order':order}


def solve_waiting_length(schedule,incoming_X,*,order=96):
    if not (0<schedule.mu<=Decimal(1)/60 and 0<schedule.delta<=schedule.mu/4):
        raise ValueError('Source waiting equation requires 0<mu<=1/60 and delta<=mu/4')
    coefficients=target_coefficients(schedule,order=order)
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        incoming=Decimal(str(incoming_X))
        ratio=(incoming-coefficients['constant'])/coefficients['exponential_coefficient']
        if ratio<1: raise ValueError('No nonnegative waiting-length root for these incoming data')
        waiting=ratio.ln()/coefficients['growth_rate']
        replay=coefficients['constant']+coefficients['exponential_coefficient']*(coefficients['growth_rate']*waiting).exp()
        return {**coefficients,'incoming_X':incoming,'waiting_length':waiting,
                'preheat_target_X_at_root':replay,'normalized_replay_defect':replay-incoming,
                'actual_heat_moment_closed':False,'other_four_moments_closed':False}


def apply_waiting_root(schedule,solution):
    return PaperOuterSchedule(logPstar=schedule.logPstar,logRref=schedule.logRref,
        delta=schedule.delta,Md=schedule.Md,c_mu=schedule.c_mu,
        c_delta=schedule.c_delta,c_epsilon=schedule.c_epsilon,
        waiting_length=solution['waiting_length'],
        sigma_quadrature_order=schedule.sigma_quadrature_order)


def independent_backward_replay(schedule,solution):
    """Reverse the two actual short source stages, independently of J formulas."""
    with localcontext() as ctx:
        ctx.prec=schedule.decimal_precision
        k=solution['growth_rate']; eq=1/k
        tail_deficit=schedule.epsilon/(1-schedule.epsilon)*(eq+Decimal(str(solution['collar_deficit_integral_K'])))
        x=float(eq+tail_deficit*(k*solution['waiting_length']).exp())
        def reverse(value,start):
            def equation(t,state):
                with localcontext() as inner:
                    inner.prec=schedule.decimal_precision
                    row=schedule.at_log_radius(schedule.logRref+start+Decimal(str(t)),0.)
                return [1-(1.5+float(row['logarithmic_slope']))*state[0]]
            result=solve_ivp(equation,(1.,0.),[value],rtol=1e-12,atol=1e-12,max_step=.025)
            if not result.success: raise ArithmeticError(result.message)
            return float(result.y[0,-1])
        x=reverse(x,schedule.y_q)
        x-=float(schedule.Ts)
        x=reverse(x,schedule.y_rel)
        return {'actual_source_backward_ODE_X':x,
                'absolute_difference_from_incoming_X':abs(x-float(solution['incoming_X'])),
                'interpretation':'Independent finite-accuracy source-stage replay; algebraic Decimal replay precision is not numerical accuracy.'}


def run():
    schedule=PaperOuterSchedule(logPstar=14,logRref=10,delta='1e-32',Md='.5',
                               c_mu='.001',c_delta='.001',c_epsilon='.01')
    incoming=incoming_angular_ratio(schedule)
    coarse=solve_waiting_length(schedule,incoming['incoming_X'],order=48)
    fine=solve_waiting_length(schedule,incoming['incoming_X'],order=96)
    matched=apply_waiting_root(schedule,fine)
    report={'source':'https://arxiv.org/html/2609.35406v1','equations':'7.9-7.10',
            'schedule':schedule.metadata(),'incoming_actual_candidate_ODE':incoming,
            'coarse':coarse,'fine':fine,
            'preheat_waiting_matched_schedule':matched.metadata(),
            'independent_backward_source_stage_replay':independent_backward_replay(schedule,fine),
            'waiting_quadrature_difference':str(abs(coarse['waiting_length']-fine['waiting_length'])),
            'full_outer_closed':False,'source_absolute_constants_certified':False,
            'scope':'Numerical root of the pre-heat scalar equation; actual heat, pressure and axial/mixed/quadratic closure still required.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2,default=str)+'\n',encoding='utf-8')
    print(json.dumps({'waiting':str(fine['waiting_length']),
                      'replay_defect':str(fine['normalized_replay_defect']),
                      'quadrature_difference':report['waiting_quadrature_difference']}),flush=True)
    return report


if __name__=='__main__': run()
