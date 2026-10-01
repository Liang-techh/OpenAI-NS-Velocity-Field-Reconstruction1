"""Analytic first-width collar integrals with directed primitive errors.

For epsilon tied to W, chi0(s)=sigma(1-s) on [0,1] and zero afterwards.
At W=0 the comparison drivers are constant, so the first width coefficient
uses Jchi(s)=s-Jsigma(s) for s<=1, and Jchi(s)=1/2 for s>=1.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_high_order_preheat_integrals import HighOrderPreheatIntegrals


def chi0(enclosures,s):
    iv=enclosures.iv;lo,hi=endpoints(s)
    if hi<=0:return iv.mpf(1)
    if lo>=1:return iv.mpf(0)
    return enclosures.sigma_interval(1-s)


def rk_first_width_quadrature(enclosures,s,steps=8):
    iv=enclosures.iv;h=s/steps;result=iv.mpf(0)
    for n in range(steps):
        t=n*h
        result+=h*(chi0(enclosures,t)+4*chi0(enclosures,t+h/2)+chi0(enclosures,t+h))/6
    return result


def run(panels=128,order=12):
    profile,alignment=accepted_profile();e=ScheduleEndpointEnclosures(profile.schedule);iv=e.iv
    with mp.workdps(e.precision+40):
        grid=HighOrderPreheatIntegrals(e).primitive_grid(panels,order);rows={}
        for text in ('.25','.5','.75','1','2'):
            s=iv.mpf(text);sv=mp.mpf(text)
            exact=s-grid[int(sv*panels)] if sv<1 else iv.mpf('.5')
            quadrature=rk_first_width_quadrature(e,s)
            exact_RK_symmetry=(sv in (1,2))
            if exact_RK_symmetry:
                # The node/weight pairs on the unit interval are symmetric;
                # sigma(x)+sigma(1-x)=1 makes the mathematical RK sum 1/2.
                assert endpoints(quadrature)[0]<=mp.mpf('.5')<=endpoints(quadrature)[1]
                quadrature=iv.mpf('.5')
            diff=exact-quadrature
            lo,hi=endpoints(diff)
            rows[text]=dict(true_integral_interval=exact,mathematical_RK_quadrature_interval=quadrature,
                true_minus_RK_interval=diff,absolute_quadrature_error_upper=max(abs(lo),abs(hi)),
                endpoint_symmetry_proves_exact_RK=exact_RK_symmetry)
        report=dict(accepted_schedule_sha256=alignment['accepted_schedule']['sha256'],
            primitive_panels=panels,Taylor_order=order,RK_steps=8,points=rows,
            first_width_ODE_coefficient_integral_enclosed=True,endpoint_first_width_quadrature_exact=True,
            second_width_ODE_error_enclosed=False,full_ODE_error_enclosed=False,
            source_driver_errors_enclosed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(report),indent=2)+'\n')
        for text,row in rows.items():print('s',text,'first-width switch quadrature error upper',mp.nstr(row['absolute_quadrature_error_upper'],16),flush=True)


if __name__=='__main__':run()
