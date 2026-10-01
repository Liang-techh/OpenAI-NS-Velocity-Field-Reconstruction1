"""Whole-axis positive Bessel reference; actual nonlinear error stays open."""
import json,math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def run():
    ctx=MPIntervalContext();ctx.dps=260
    with mp.workdps(300):
        z=ctx.mpf([-1,1]);dt=ctx.mpf('1e-200');j=ctx.mpf('1e-14');lam=ctx.mpf('1e36')
        L=1-dt*z**2;beta=(3-dt/2+j*z)/L
        # chi=H0²/(H0²+sigma²) in [0,1] identically, including H0=0.
        q=ctx.mpf([0,endpoints(ctx.mpf('4.1')*(1+beta/lam)/2)[1]])
        qmax=ctx.mpf(endpoints(q)[1])
        if endpoints(q)[1]>=mp.mpf('2.1'):raise AssertionError('Bessel argument bound')
        # P3 is decreasing: P3'= -((q-4)^2+8)/48. The alternating
        # tail after degree3 is nonnegative because q/20<1 and decreases.
        lower=1-qmax/2+qmax*qmax/12-qmax**3/144
        margin=lower-ctx.mpf(1)/8
        if endpoints(margin)[0]<=0:raise AssertionError('Reference positivity margin')
        degree18_tail=qmax**19/(math.factorial(19)*math.factorial(20))
        # Independent resolved endpoint value checks the exact B series.
        value=mp.besselj(1,2*mp.sqrt(endpoints(qmax)[1]))/mp.sqrt(endpoints(qmax)[1])
        if value<endpoints(lower)[0] or value>1:raise AssertionError('Independent Bessel endpoint')
        report=dict(axial_domain=['-1','1'],scaled_radial_domain=['0','4.1'],beta=beta,q=q,
            Bessel_reference_lower=lower,Bessel_reference_upper=1,
            degree18_Bessel_series_tail_upper=degree18_tail,
            normalized_error_budget_to_retain_one_eighth_lower=margin,
            sufficient_K8_upper_if_error_le_K8_over_Lambda=lam*margin,
            independent_Bessel_endpoint=value,reference_positivity_certified=True,
            core_model='B(q)=sum (-q)^n/(n!(n+1)!), q=scaled_R*(chi+beta/Lambda)/2',
            zero_H0_case='Retain beta/Lambda: q does not identically vanish for scaled_R>0',
            constant_scope='K8 is the normalized core comparison constant, not the later global data-size K',
            actual_nonlinear_comparison_error_enclosed=False,actual_core_positivity_certified=False,
            infinite_radial_remainder_enclosed=False,full_K_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('whole-axis Bessel lower',mp.nstr(endpoints(lower)[0],18),'error budget',mp.nstr(endpoints(margin)[0],18),flush=True)


if __name__=='__main__':run()
