"""Logarithmic paper collar markers and conditional inlet cone scales.

K and c_star are the existing numerical policy, not certified C3 norms or
theorem constants. Preserve log(Ran/Ra) without subtracting nearby radii.
"""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_collar_interval_receipt_reader import read_inlet
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def log_sigma(ctx,s):
    lo,hi=endpoints(s)
    if not 0<lo<=hi<mp.mpf('.5'):
        raise ValueError('Logarithmic inverse uses 0<s<1/2')
    logit=1/(1-s)**2-1/s**2
    # log(1+exp(logit)) is between zero and exp(logit). This avoids
    # erasing the correction by forming 1+an exponentially small number.
    return logit-ctx.mpf([0,endpoints(ctx.exp(logit))[1]])


def inverse_log_sigma(ctx,target,iterations=160):
    low_target,high_target=endpoints(target)
    if high_target>=-10:raise ValueError('Require a small positive target')
    center=1/mp.sqrt(-low_target);left=center/2;right=center*2
    if not endpoints(log_sigma(ctx,ctx.mpf(left)))[1]<low_target:
        raise AssertionError('Left inverse bracket')
    if not endpoints(log_sigma(ctx,ctx.mpf(right)))[0]>high_target:
        raise AssertionError('Right inverse bracket')
    for _ in range(iterations):
        middle=(left+right)/2;lo,hi=endpoints(log_sigma(ctx,ctx.mpf(middle)))
        if hi<low_target:left=middle
        elif lo>high_target:right=middle
        else:break
    return ctx.mpf([left,right])


def factored_margins(ctx,D,E,chi,rho):
    """Cone lower bounds assuming |e|<=rho*(1-chi), I/F=q+e.

    Divide out positive edge factors before interval evaluation. At chi=0
    or 1 these are limiting coefficients, not strict physical margins.
    """
    if endpoints(D)[0]<=0 or endpoints(rho)[0]<0:
        raise ValueError('Require D>0 and rho>=0')
    cl,ch=endpoints(chi)
    if not 0<=cl<=ch<=1:raise ValueError('Require chi in [0,1]')
    Q2=D*D+E*E;Q=ctx.sqrt(Q2);kappa=chi*Q2/D
    dot_lower=Q2-rho*Q
    # Angular cone bound requires the dot lower bound to be positive
    # before squaring it. Otherwise no lower bound is claimed.
    angular=None
    if endpoints(dot_lower)[0]>0 and endpoints(kappa)[0]>2:
        angular=2*dot_lower**2-(kappa-2)*rho*rho*Q2
    relaxed=Q2-2*D-(1-chi)*rho*Q
    return dict(kappa=kappa,negative_dot_factored_lower=dot_lower,
                strict_angular_factored_lower=angular,
                relaxed_small_kappa_factored_lower=relaxed)


def run():
    inlet,source=read_inlet();ctx=inlet['ctx']
    with mp.workdps(ctx.dps+40):
        logK=ctx.mpf('1e152');logc=ctx.mpf(-100);gamma=ctx.mpf('.01')
        logh=logc-100*logK;target=ctx.log(gamma/10)-10*logK
        marker=inverse_log_sigma(ctx,target)
        # Work with a logarithmic displacement: Ran/Ra-1 cannot be
        # resolved by adding exp(logh)*marker to 1 at ordinary precision.
        log_displacement=logh+ctx.log(marker)
        F=inlet['F'].value.evaluate(width=0).value
        It=inlet['I_theta'].value.evaluate(width=0).value
        Iz=inlet['I_z'].value.evaluate(width=0).value
        if endpoints(F)[0]<=0 or endpoints(It)[0]<=0:
            raise AssertionError('Conditional inlet cone prerequisites')
        log_ratio=ctx.log(It*It+Iz*Iz)-ctx.log(F)-ctx.log(It)
        log_exit_kappa=logh+log_ratio
        if endpoints(log_exit_kappa)[1]>=0:
            raise AssertionError('Numerical inlet direction lacks exit kappa<1')
        lower_K_logs=dict(inverse_F=-ctx.log(F),absolute_E=ctx.log(abs(Iz)/F),
                          inverse_R=-ctx.log(inlet['R'].value.evaluate(width=0).value))
        for value in lower_K_logs.values():
            if endpoints(value)[1]>=endpoints(logK)[0]:
                raise AssertionError('Assumed K below a retained inlet norm contribution')
        # Independent inverse at a resolvable target checks sigma convention
        # and directed bracketing, independently of the huge actual logK.
        fixture=inverse_log_sigma(ctx,ctx.log(ctx.mpf('1e-20')))
        with mp.workdps(100):
            root=mp.findroot(lambda x:1/(1-x)**2-1/x**2-mp.log(1+mp.exp(1/(1-x)**2-1/x**2))-mp.log(mp.mpf('1e-20')),
                             (mp.mpf('.13'),mp.mpf('.16')))
        lo,hi=endpoints(fixture)
        if not lo<=root<=hi:raise AssertionError('Independent inverse fixture')
        # Resolvable independent fixture of the edge-factored inequality.
        margins=factored_margins(ctx,ctx.mpf(4),ctx.mpf(3),ctx.mpf('.9'),ctx.mpf('.1'))
        # Q=5, kappa=45/8, dot lower=49/2, angular=38387/32,
        # relaxed=339/20. Check exact rationals in directed enclosures.
        for name,value in (('kappa',mp.mpf(45)/8),('negative_dot_factored_lower',mp.mpf(49)/2),
                           ('strict_angular_factored_lower',mp.mpf(38387)/32),
                           ('relaxed_small_kappa_factored_lower',mp.mpf(339)/20)):
            lo,hi=endpoints(margins[name])
            if not lo<=value<=hi:raise AssertionError('Factored cone fixture: '+name)
        result=dict(source=source,gamma=gamma,assumed_log_K=logK,assumed_log_c_star=logc,
            log_h_b=logh,log_epsilon_b=logh,s_an=marker,
            log_sigma_s_an_target=target,log_of_log_Ran_over_Ra=log_displacement,
            strict_region='0<s<=s_an, excluding the zero-stress inlet',
            relaxed_region='s>s_an through the later connecting interval',
            inlet_log_normalized_direction_ratio=log_ratio,
            inlet_direction_log_kappa_after_switch=log_exit_kappa,
            pointwise_log_K_lower_contributions=lower_K_logs,
            pointwise_necessary_K_tests_passed=True,independent_inverse_fixture_passed=True,
            factored_margin_fixture_passed=True,factored_margin_fixture=margins,
            global_C3_K_bound_certified=False,c_star_theorem_restrictions_certified=False,
            actual_finite_width_cone_certified=False,omitted_orders_enclosed=False,
            marker_scope='Conditional marker for existing logK/logc policy; not a certified collar',
            temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n')
        print('conditional s_an',mp.nstr(endpoints(marker)[0],18),flush=True)
        print('inlet-direction log kappa after switch',mp.nstr(endpoints(log_exit_kappa)[1],18),flush=True)


if __name__=='__main__':run()
