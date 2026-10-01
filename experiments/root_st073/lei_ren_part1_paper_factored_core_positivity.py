"""Bernstein lower bounds for the amplitude-factored finite radial core.

Uniform intervals that cross zero are inconclusive, never counterexamples.
The radial polynomial is finite: omitted radial orders remain open.
"""
import json,math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_difference import IntervalDifference
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_global_finite_core_bounds import accepted_pressure_axis_rows
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def bernstein_enclosure(ctx,coefficients,left=0,right=1):
    """Power-to-Bernstein conversion followed by the convex-hull bound."""
    left=ctx.mpf(left);right=ctx.mpf(right)
    if endpoints(left)[1]>=endpoints(right)[0]:raise ValueError('Increasing radial interval required')
    values=[v.value if isinstance(v,IntervalDifference) else ctx.mpf(v) for v in coefficients]
    degree=len(values)-1;shifted=[]
    for k in range(degree+1):
        shifted.append((right-left)**k*sum((values[n]*math.comb(n,k)*left**(n-k)
                                           for n in range(k,degree+1)),ctx.mpf(0)))
    b=[sum((shifted[i]*ctx.mpf(math.comb(k,i))/math.comb(degree,i) for i in range(k+1)),ctx.mpf(0))
       for k in range(degree+1)]
    lower=min(endpoints(v)[0] for v in b);upper=max(endpoints(v)[1] for v in b)
    return dict(coefficients=b,range=ctx.mpf([lower,upper]),strictly_positive=lower>0,
                whole_radial_interval_enclosed=True)


def squared_axis_rows(ctx,gradient,lam,F0,length):
    """Taylor rows of F0² from its logarithmic derivative, without division."""
    ratio=[ctx.mpf(1)]
    for n in range(length-1):
        ratio.append(-2*lam*sum((gradient[k]*ratio[n-k] for k in range(n+1)),ctx.mpf(0))/(n+1))
    return [F0**2*v for v in ratio]


def run():
    from lei_ren_part1_paper_amplitude_factored_core import factored_core_coefficients
    ctx=MPIntervalContext();ctx.dps=473;base=Path(__file__).parent;degree=18;length=degree+4
    with mp.workdps(ctx.dps+40):
        lam=ctx.mpf('1e36');r=4/lam;globalz=ctx.mpf([-1,1])
        p0,alignment=accepted_pressure_axis_rows(ctx,globalz,length,base/'lei_ren_part1_paper_global_pressure_high_derivatives.json')
        convert=IntervalDifference.converter(ctx);rows=[]
        for label,z in (('whole_axis',globalz),('point_Z03',ctx.mpf('.3'))):
            axis=uniform_axis_jets(ctx,radius=1,j='1e-14',Lambda=lam,logC='5e151',delta='1e-200',length=length,axial_interval=z)
            gradient=axis['gradient_coefficients'];ell=[-lam*v for v in gradient]
            squared=squared_axis_rows(ctx,gradient,lam,axis['F0_interval'],length)
            core=factored_core_coefficients(z,'1e-200',ell_Z_taylor=ell,S_Z_taylor=squared,
                U0_Z_taylor=axis['U0'],P0_Z_taylor=p0,radial_degree=degree,precision=ctx.dps,scalar_converter=convert)
            coefficients=[row[0]*r**n for n,row in enumerate(core['A'])]
            bound=bernstein_enclosure(ctx,coefficients)
            data=dict(label=label,Z_interval=z,normalized_core_F_over_F0=bound)
            if bound['strictly_positive']:
                # C3 reciprocal bounds at this axial slice over every radius.
                # The axis factor remains outside all normalized arithmetic.
                A=[bound['range']]
                for k in range(1,4):
                    box=bernstein_enclosure(ctx,[row[k]*r**n for n,row in enumerate(core['A'])])
                    A.append(box['range'])
                reciprocal=[1/A[0]]
                for k in range(1,4):
                    reciprocal.append(-sum((A[i]*reciprocal[k-i] for i in range(1,k+1)),ctx.mpf(0))/A[0])
                invaxis=[ctx.mpf(1)]
                for n in range(3):
                    invaxis.append(lam*sum((gradient[k]*invaxis[n-k] for k in range(n+1)),ctx.mpf(0))/(n+1))
                inverse_ratios=[math.factorial(k)*sum((invaxis[i]*reciprocal[k-i] for i in range(k+1)),ctx.mpf(0)) for k in range(4)]
                ratio_sum=sum((ctx.mpf(max(abs(v) for v in endpoints(value))) for value in inverse_ratios),ctx.mpf(0))
                data['reciprocal_F_axial_C3_log_upper']=endpoints(-ctx.log(axis['F0_interval'])+ctx.log(ratio_sum))[1]
                data['reciprocal_F_normalized_axial_derivative_bounds']=inverse_ratios
                data['reciprocal_scope']='At the stated Z slice, all R in [0,Ra]; not global-Z or mixed C3'
            rows.append(data)
            print(label,'positive',bound['strictly_positive'],'lower',mp.nstr(endpoints(bound['range'])[0],15),flush=True)
        # Independent scalar polynomial catches basis/index conventions.
        fixture=bernstein_enclosure(ctx,[ctx.mpf('1.1'),ctx.mpf(-2),ctx.mpf(1)])
        lo,hi=endpoints(fixture['range'])
        if not fixture['strictly_positive'] or not lo<=mp.mpf('.1')<=hi:
            raise AssertionError('Bernstein independent polynomial fixture')
        result=dict(accepted_schedule_sha256=alignment['accepted_schedule']['sha256'],
            radial_degree=degree,normalized_radial_domain=['0','1'],axis_domain=['-1','1'],enclosures=rows,
            independent_polynomial_fixture_passed=True,axis_amplitude_factored_exactly=True,
            nonlinear_swirl_coupling_retained=True,pressure_datum_preserved=True,
            source_pressure_error_included=True,whole_axis_finite_core_positivity_certified=rows[0]['normalized_core_F_over_F0']['strictly_positive'],
            crossed_zero_interpretation='inconclusive interval dependence, not an observed zero',
            infinite_radial_remainder_enclosed=False,full_K_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n')


if __name__=='__main__':run()
