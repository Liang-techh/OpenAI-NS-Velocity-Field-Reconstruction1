"""Independent derivative checks and nonlinear pressure-error containment."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_core_recursion import core_coefficients
from lei_ren_part1_paper_component_pressure_core import fixture as component_fixture,evaluate_component_core_coefficients
from lei_ren_part1_paper_uniform_pressure_high_derivatives import derivative_envelope
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    iv=MPIntervalContext();iv.dps=70
    checks=0
    with mp.workdps(110):
        for beta in (mp.mpf('.3'),mp.mpf('1.4'),mp.mpf(2)):
            for z in (mp.mpf('-.8'),mp.mpf(0),mp.mpf('.4'),mp.mpf('.8')):
                for n in range(9):
                    value=mp.diff(lambda w:(1+w*w)**(-beta),z,n)
                    upper=endpoints(derivative_envelope(iv,'.8',n))[1]
                    assert abs(value)<=upper;checks+=1
        z=mp.mpf('.3');delta=mp.mpf('.1');length=6
        f0=[mp.mpf(3),mp.mpf('.2'),mp.mpf('.1')]+[mp.mpf(0)]*3
        u0=[mp.mpf(1),mp.mpf('.4'),mp.mpf('.2')]+[mp.mpf(0)]*3
        base=[mp.mpf('1.1'),mp.mpf('-.2'),mp.mpf('.4')]+[mp.mpf(0)]*3
        error=[mp.mpf('.02')/(k+1) for k in range(length)]
        uncertain=[iv.mpf(p)+iv.mpf([-e,e]) for p,e in zip(base,error)]
        kwargs=dict(F0_Z_taylor=f0,U0_Z_taylor=u0,radial_degree=3,precision=110)
        enclosure=core_coefficients(z,delta,P0_Z_taylor=uncertain,scalar_converter=iv.mpf,**kwargs)
        containment=0;moment_containment=0
        interval_field=evaluate_component_core_coefficients(enclosure,iv.mpf(".2"),iv.mpf(z),iv.mpf(delta),square_root=iv.sqrt)
        for direction in (-1,0,1):
            sample=[p+direction*(-1)**k*e for k,(p,e) in enumerate(zip(base,error))]
            result=core_coefficients(z,delta,P0_Z_taylor=sample,**kwargs)
            for key in ('F','Uz','P'):
                for row,intervals in zip(result[key],enclosure[key]):
                    for value,interval in zip(row,intervals):
                        lo,hi=endpoints(interval)
                        assert lo<=value<=hi,(key,value,lo,hi)
                        containment+=1
            sample_field=evaluate_component_core_coefficients(result,mp.mpf('.2'),z,delta)
            for key in ('moments','moments_Z'):
                for name,value in sample_field[key].items():
                    lo,hi=endpoints(interval_field[key][name]);assert lo<=value<=hi
                    moment_containment+=1
        component=component_fixture()
    report=dict(all_checks_passed=True,independent_q_power_derivative_checks=checks,
        nonlinear_interval_coefficient_containment_checks=containment,
        five_moment_interval_containment_checks=moment_containment,
        pressure_polynomial_arithmetic_regression=component,
        point_samples_are_regression_only=True,uniform_pressure_proof='analytic coefficient envelope',
        core_fixture_scope='resolved finite radial recurrence; no exact infinite profile or RK claim')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Derivative envelope checks',checks,'nonlinear coefficient containment checks',containment,'component fixture passed')


if __name__=='__main__':run()
