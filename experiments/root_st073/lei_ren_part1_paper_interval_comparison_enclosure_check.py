"""Closed-form frozen-field and nonzero-slope numerical comparisons."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets
from lei_ren_part1_paper_candidate_shared_inlet import radial_product
from lei_ren_part1_paper_interval_comparison_enclosure import integrate, frozen_extension
from lei_ren_part1_paper_interval_taylor import IntervalTaylor, constant
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def fixture(nonzero_slope=False):
    calc = CandidateComparisonJets.__new__(CandidateComparisonJets)
    c = MPIntervalContext()
    c.dps = 80
    calc.ctx,calc.precision = c,80
    calc.center = mp.mpf('.5')
    calc.center_family = '.5'
    calc.z = IntervalTaylor.variable(c,c.mpf('.5'),3)
    calc.lam,calc.eps,calc.delta = c.mpf(10),c.mpf('.1'),c.mpf('.01')
    calc.ell,calc.S,calc.p0 = constant(c,0,3),constant(c,'.5',3),constant(c,-3,3)
    calc.phi = [1+calc.z*calc.z/10]
    calc.u = [2+calc.z/5]
    if nonzero_slope:
        calc.phi.append(constant(c,'.1',3))
        calc.u.append(constant(c,'-.05',3))
    calc.phi_squared = radial_product(calc.phi,calc.phi)
    calc.phi_u = radial_product(calc.phi,calc.u)
    calc.u_squared = radial_product(calc.u,calc.u)
    calc.initial_phi = calc.radial(calc.phi,c.mpf(4))
    calc.hb = mp.mpf('.005')
    calc.steps = 256
    calc.cache,calc.radial_cache = {},{}
    calc.degree,calc.state_hash = len(calc.phi)-1,'polynomial-fixture'
    calc.tail = dict(target_met=True,rows=[dict(scaled_radial_order=i,axial_order=k,
        Phi_tail=c.mpf(0),Psi_tail=c.mpf(0)) for i in range(2) for k in range(4-i)])
    return calc


def containment(actual, expected):
    count = 0
    for key,value in expected.items():
        for k in range(value.order+1):
            lo,hi = endpoints(actual[key][k])
            el,eh = endpoints(value[k])
            if not lo <= el <= eh <= hi:
                raise AssertionError((key,k,'reference outside directed comparison enclosure'))
            count += 1
    return count


def run():
    with mp.workdps(120):
        frozen = fixture()
        bound = integrate(frozen,16)
        c = frozen.ctx
        s = 4*c.exp(c.mpf('.01'))
        f,u = frozen.phi[0].truncate(2),frozen.u[0].truncate(2)
        expected = dict(phi=f,U=u,theta=f*s*s,z=u*s,theta_z=f*u*s*s,
            p=f*f*s,u_squared=u*u*s,weighted_phi_squared=f*f*s*s/2)
        exact_count = containment(bound['endpoint_state'],expected)
        extended = frozen_extension(frozen,bound,c.mpf(100))
        r = c.mpf(100)
        expected_extended = dict(phi=f,U=u,theta=f*r*r,z=u*r,theta_z=f*u*r*r,
            p=f*f*r,u_squared=u*u*r,weighted_phi_squared=f*f*r*r/2)
        extended_count = containment(extended['state'],expected_extended)
        varying = fixture(True)
        varied_bound = integrate(varying,32)
        reference = varying.evaluate('.01','.5')
        numerical_count = containment(varied_bound['endpoint_state'],
            {key:value.truncate(2) for key,value in reference['state'].items()})
        driver_count = containment(varied_bound,
            {key:reference[key].truncate(1) for key in ('D','I_z','pressure')})
        report = dict(passed=True,fixture_only=True,closed_form_field_coefficient_checks=exact_count,
            nonzero_slope_RK256_coefficient_checks=numerical_count,
            nonzero_slope_driver_coefficient_checks=driver_count,
            closed_form_frozen_extension_coefficient_checks=extended_count,
            directed_cells_closed_form=16,directed_cells_nonzero_slope=32,
            axial_field_orders=[0,1,2],driver_axial_orders=[0,1],
            numerical_reference_is_accuracy_diagnostic_not_proof=True,
            construction_parameter_errors_certified=False,
            source_hashes={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                for name in ('lei_ren_part1_paper_interval_comparison_enclosure.py',
                             'lei_ren_part1_paper_interval_comparison_enclosure_check.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Directed comparison fixtures:',exact_count,'closed-form,',numerical_count,
              'nonzero-slope field,',driver_count,'driver coefficients',flush=True)
        return report


if __name__ == '__main__':
    run()
