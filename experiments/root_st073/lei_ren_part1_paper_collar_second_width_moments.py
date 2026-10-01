"""Analytic second-width endpoint coefficients for all collar states.

The dynamic g/u coefficients use core-inlet radial inertial derivatives;
the six integral coefficients require only the first-width solutions.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_coherent_pressure_error_transfer import accepted_profile
from lei_ren_part1_paper_schedule_endpoint_enclosures import ScheduleEndpointEnclosures,endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_high_order_preheat_integrals import switch_taylor
from lei_ren_part1_paper_interval_taylor import IntervalTaylor,integrate_symmetric


def endpoint_integrated_switch(e,panels=128,order=12):
    """Enclose integral_0^2 Jchi(s) ds = 1/2 + integral_0^1 x sigma(x) dx."""
    if not isinstance(panels,int) or panels<4 or not isinstance(order,int) or order<2:
        raise ValueError('Require panels>=4 and order>=2')
    iv=e.iv;h=iv.mpf(1)/panels;radius=h/2;edge_switch=e.sigma_interval(h)
    first_upper=endpoints(h*h*edge_switch/2)[1]
    total=iv.mpf([0,first_upper])
    for i in range(1,panels-1):
        center=(iv.mpf(i)+iv.mpf('.5'))/panels
        cell=iv.mpf([endpoints(iv.mpf(i)/panels)[0],endpoints(iv.mpf(i+1)/panels)[1]])
        def density(x):return IntervalTaylor.variable(iv,x,order)*switch_taylor(iv,x,order)
        total+=integrate_symmetric(density(center),density(cell),radius)
    # Flat-end symmetry: sigma(x)=1-sigma(1-x), so the last cell's
    # deficit is positive and bounded by its length times sigma(h).
    last_exact=h-h*h/2;last_deficit=endpoints(h*edge_switch)[1]
    total+=last_exact-iv.mpf([0,last_deficit])
    return dict(weighted_sigma_integral=total,integrated_switch_primitive=iv.mpf('.5')+total,
                panels=panels,Taylor_order=order)


def second_width_integral_coefficients(inlet,integrated_switch):
    """Return endpoint width^2 coefficients divided by h_b^2.

inlet quantities are directed second-Z jets and integrated_switch encloses
K=integral_0^2 Jchi. No second-width comparison driver is needed here.
"""
    u=inlet['Uz'];D=inlet['D'];B=(inlet['R']/2).sqrt()*inlet['I_z'];K=integrated_switch
    # Right multiplication dispatches through the interval axial ring.
    DK=D*K;BK=B*K
    return dict(theta=8-DK,mz=2*u-BK,mixed=8*u-u*DK-2*BK,
                axial=2*u*u-2*u*BK,swirl=4-DK,p=2-DK)


def second_width_coefficients(inlet,integrated_switch):
    """All endpoint coefficients / h_b^2 for width-tied epsilon.

chi0 is supported in s<=1, where the comparison equals the exact finite
core. Thus D1=s Ra D_R and Iz1=s Ra Iz_R; no transition RK replay is needed.
The supplied derivative jets must be from this same inlet and pressure.
"""
    K=integrated_switch;M=1-K;R=inlet['R'];D=inlet['D'];Iz=inlet['I_z']
    radial_log_slope=R*inlet['F_R']/inlet['F']
    result=second_width_integral_coefficients(inlet,K)
    result['g']=-(R*inlet['D_R']*M+D*3/2)/2
    result['u']=-(R/2).sqrt()*(R*inlet['I_z_R']*M+
        Iz*((1/2-radial_log_slope)*M-D/16+3/2))
    return result


def run():
    profile,_=accepted_profile();e=ScheduleEndpointEnclosures(profile.schedule)
    with mp.workdps(e.precision+40):
        result=endpoint_integrated_switch(e)
        # Independent resolved quadrature checks this directed scalar integral.
        from lei_ren_part1_paper_axial_primitive import _sigma_mp
        with mp.workdps(90):
            resolved=mp.quad(lambda x:x*_sigma_mp(x),[0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
        lo,hi=endpoints(result['weighted_sigma_integral'])
        if not lo<=resolved<=hi:raise AssertionError('Weighted switch enclosure misses independent quadrature')
        result.update(independent_weighted_integral=resolved,
            independent_quadrature_contained=True,six_second_width_integral_formulas_available=True,
            eight_second_width_endpoint_formulas_available=True,
            actual_inlet_second_width_coefficients_evaluated=False,full_ODE_error_enclosed=False,
            source_driver_error_enclosed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(result),indent=2)+'\n')
        print('integrated switch',mp.nstr(endpoints(result['integrated_switch_primitive'])[0],18),
              'enclosure width',mp.nstr(hi-lo,10),flush=True)


if __name__=='__main__':run()
