"""Explicit factorial linear inverse bound from paper (8.32)-(8.35).

Uses complex-axis modulus enclosures and Cauchy estimates. This conservative
operator bound does not supply the nonlinear Lipschitz constant Kstar.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def factorial_series_bound(ctx,K,relative_tail='1e-100'):
    """Directed sum of K^n/[n!(n+1)!], including its infinite tail."""
    K=ctx.mpf(K);threshold=ctx.mpf(relative_tail)
    if endpoints(K)[0]<0:raise ValueError('Nonnegative factorial constant required')
    term=ctx.mpf(1);total=term;n=0
    while n<10000:
        next_ratio=K/((n+1)*(n+2))
        next_term=term*next_ratio
        if endpoints(next_ratio)[1]<1 and endpoints(next_term)[1]<=endpoints(total*threshold)[0]:
            tail=next_term/(1-next_ratio)
            return dict(partial_sum=total,infinite_tail_upper=tail,
                        resolvent_norm_upper=total+tail,last_summed_order=n,
                        unsummed_term_ratio_upper=next_ratio)
        n+=1;term=next_term;total+=term
    raise RuntimeError('Factorial inverse summation did not close')


def bessel_model_norm_bound(ctx,coefficient_modulus,weight):
    """Bound the model norm directly, avoiding an operator inverse bound.

For each n Cauchy bounds its coefficient by chi_mod^n/(2^n n!(n+1)!).
The axial weight divided by binomial(n+m,m) is no larger than weight.
The positive radial sequence decreases forever after its ratio falls below1.
"""
    K=10*coefficient_modulus;term=ctx.mpf(1);upper=endpoints(term)[1];n=0
    while n<10000:
        ratio=K*(n+2)/(n+1)**3
        if endpoints(ratio)[1]<=1:
            return dict(model_Xh_norm_upper=ctx.mpf([0,upper])*weight,
                        last_checked_radial_order=n,subsequent_weighted_coefficient_ratio_upper=ratio)
        n+=1;term*=ratio;upper=max(upper,endpoints(term)[1])
    raise RuntimeError('Model coefficient maximum did not close')


def run():
    path=Path(__file__).with_name('lei_ren_part1_paper_shared_analytic_tube.json')
    source=json.loads(path.read_text());ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def restore(name):
            row=source[name]
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        if not source['common_complex_axis_poles_excluded']:
            raise AssertionError('Analytic input guard missing')
        eta=restore('complex_tube_radius');h=restore('Xh_parameter')
        ratio=h/(eta/2)
        if endpoints(ratio)[1]>mp.mpf('0.4'):
            raise AssertionError('Cauchy square-weight bound unavailable')
        # For radial-independent analytic data Cauchy on eta/2 gives
        # |d^m f|/m! <= modulus/(eta/2)^m. Since (m+1)^2/4^m<=1,
        # the entire X_h norm is bounded by the complex modulus.
        # Use max(1,4r), with outward rounding. For m>=1 the successive
        # weight ratio is <=(3/2)^2 r<1, so the maximum is at m=0 or1.
        weight=ctx.mpf([1,max(mp.mpf(1),endpoints(4*ratio)[1])])
        chi_norm=restore('chi_modulus_upper')*weight;beta_norm=restore('beta_modulus_upper')*weight
        multiplier=256*(chi_norm+beta_norm/500)
        factorial_constant=40*multiplier
        bound=factorial_series_bound(ctx,factorial_constant)
        direct_model=bessel_model_norm_bound(ctx,restore('chi_modulus_upper'),weight)
        report=dict(input_receipt=path.name,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             precision=160,Xh_parameter=h,Cauchy_disk_radius=eta/2,
             Cauchy_weight_supremum_upper=weight,product_constant=256,
             chi_Xh_norm_upper=chi_norm,beta_Xh_norm_upper=beta_norm,
             uniform_epsilon_range=['0','1/500'],multiplier_norm_upper=multiplier,
             factorial_inverse_constant=factorial_constant,**bound,
             resolvent_norm_log_upper=ctx.log(bound['resolvent_norm_upper']),
             Phi_model_Xh_norm_upper=direct_model['model_Xh_norm_upper'],
             direct_Bessel_model_norm=direct_model,
             Phi_model_inverse_based_alternative_upper=bound['resolvent_norm_upper'],
             analytic_linear_inverse_bound_certified=True,
             bound_is_conservative_not_actual_operator_norm=True,
             analytic_pressure_modulus_certified=False,Psi_model_Xh_norm_certified=False,
             nonlinear_Kstar_certified=False,actual_Lambda_contraction_certified=False,
             infinite_nonlinear_core_remainder_enclosed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('factorial inverse log upper',mp.nstr(endpoints(report['resolvent_norm_log_upper'])[1],22),
              'summed orders',bound['last_summed_order'],flush=True)
        return report


if __name__=='__main__':run()
