"""Uniform pressure approximation derivative bounds, using analytic q powers.

The comparison is the stored finite mathematical datum, not runtime roundoff.
All bounds concern |Z|<=a and stored schedule parameters only.
"""
import json, math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval


def derivative_envelope(ctx, radius, n):
    """For beta in [0,2], bound |d^n (1+Z^2)^(-beta)|.

    Expand q(Z+h)=q(Z)+2Zh+h^2. For j quadratic choices,
    k=n-j and m=n-2j. The absolute Taylor coefficient is at most
    (k+1)! (2a)^m/(m! j!), since q>=1 and (beta)_k<=(2)_k.
    Multiplication by n! converts coefficients to derivatives.
    """
    if not isinstance(n,int) or isinstance(n,bool) or n<0:raise ValueError('Nonnegative integer order required')
    a=ctx.mpf(radius)
    if endpoints(a)[0]<0 or endpoints(a)[1]>=1:raise ValueError('Require 0<=radius<1')
    total=ctx.mpf(0)
    for j in range(n//2+1):
        k=n-j;m=n-2*j
        total+=ctx.mpf(math.factorial(k+1))*(2*a)**m/(math.factorial(m)*math.factorial(j))
    return math.factorial(n)*total


def run(degree=24):
    iv=MPIntervalContext();iv.dps=300
    fixed=json.loads(Path(__file__).with_name('lei_ren_part1_paper_uniform_fixed_beta_error.json').read_text())
    full=json.loads(Path(__file__).with_name('lei_ren_part1_paper_uniform_preheat_error_budget.json').read_text())
    with mp.workdps(330):
        varying=iv.mpf(0);constant=iv.mpf(0)
        for stage,row in fixed['stages'].items():
            err=read_interval(iv,row['mass_error_interval']);lo,hi=endpoints(err);size=iv.mpf(max(abs(lo),abs(hi)))
            if row['beta']==2:varying+=size
            else:constant+=size
        true_mass=iv.mpf(mp.make_mpf(tuple(full['flatten_true_mass_upper']['exact_mpf_tuple'])))
        finite_mass=read_interval(iv,full['flatten_finite_mass_interval'])
        flatten= true_mass+finite_mass
        bounds=[];rows=[]
        for n in range(degree+1):
            factor=derivative_envelope(iv,full['radius'],n)
            error=(varying+flatten)*factor+(constant if n==0 else 0)
            bounds.append(endpoints(error)[1])
            rows.append(dict(order=n,q_power_derivative_envelope=factor,
                normalized_derivative_error_upper=endpoints(error)[1],
                normalized_Taylor_coefficient_error_upper=endpoints(error/math.factorial(n))[1]))
        report=dict(radius=full['radius'],maximum_derivative_order=degree,stage_count=14,
            derivatives=rows,normalized_derivative_error_upper_bounds=bounds,
            uniform_Z_error_enclosed=True,all_pressure_stages_included=True,
            flatten_relative_error_enclosed=False,adapter_evaluation_roundoff_enclosed=False,
            original_parameter_errors_enclosed=False,core_RK_error_enclosed=False,
            five_defect_interval_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        print('uniform pressure derivatives through',degree,'order3 error',mp.nstr(bounds[3],16),'order24 coefficient error',mp.nstr(bounds[-1]/math.factorial(degree),16))


if __name__=='__main__':run()
