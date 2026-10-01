"""Tighter X_h inverse bound using commuting radial/axial operators.

A=(chi/2)J2, hence A^n=(chi^n/2^n)J2^n. Bound whole powers
directly instead of repeatedly applying a generic product constant.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def weighted_factorial_tail(ctx,K,relative_tail='1e-100'):
    """Directed sum over n>=1 of K^n(n+1)^2/[n!(n+1)!]."""
    term=ctx.mpf(1);total=ctx.mpf(0);K=ctx.mpf(K)
    if endpoints(K)[0]<=0:raise ValueError('Positive K required')
    for n in range(1,10001):
        term*=K*(n+1)/n**3
        total+=term
        ratio=K*(n+2)/(n+1)**3
        next_term=term*ratio
        if endpoints(ratio)[1]<1 and endpoints(next_term)[1]<=endpoints(total*ctx.mpf(relative_tail))[0]:
            remainder=next_term/(1-ratio)
            return dict(positive_series_upper=total+remainder,partial_sum=total,
                        tail_upper=remainder,last_summed_order=n,next_ratio_upper=ratio)
    raise RuntimeError('Weighted factorial tail did not close')


def run():
    path=Path(__file__).with_name('lei_ren_part1_paper_analytic_radial_tail.json')
    source=json.loads(path.read_text());ctx=MPIntervalContext();ctx.dps=160
    with mp.workdps(200):
        def read(name):
            row=source[name]
            return ctx.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                            mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
        eta=read('complex_tube_radius');h=read('Xh_parameter')
        r=h/(eta/2)
        if endpoints(r)[1]>=1:raise AssertionError('Cauchy convolution not summable')
        S=(1+r)/(1-r)**3
        C=read('chi_modulus_upper');K=10*C
        series=weighted_factorial_tail(ctx,K)
        upper=1+S*series['positive_series_upper']
        report=dict(input_receipt=path.name,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
             precision=160,Xh_parameter=h,Cauchy_radius=eta/2,axial_ratio=r,
             axial_convolution_factor_upper=S,chi_complex_modulus_upper=C,
             radial_factorial_constant=K,**series,resolvent_norm_upper=upper,
             resolvent_norm_log_upper=ctx.log(upper),
             operator='R=(I+(chi/2)J2)^(-1)',
             exact_power_identity='A^n=(chi^n/2^n)J2^n',
             per_power_norm_bound='S(r)*(10C)^n*(n+1)^2/[n!(n+1)!], n>=1',
             all_radial_and_axial_orders_enclosed=True,
             no_generic_product_constant_per_power=True,
             nonlinear_contraction_certified=False,infinite_core_remainder_enclosed=False,
             temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('commuting resolvent log upper',mp.nstr(endpoints(ctx.log(upper))[1],22),
              'summed orders',series['last_summed_order'],flush=True)
        return report


if __name__=='__main__':run()
