"""Independent closed-form polynomial check of core moment error bounds."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets
from lei_ren_part1_paper_candidate_shared_inlet import radial_product
from lei_ren_part1_paper_interval_taylor import IntervalTaylor, constant
from lei_ren_part1_paper_interval_core_inlet_bound import enclose
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def run():
    c = MPIntervalContext()
    c.dps = 80
    with mp.workdps(120):
        z = IntervalTaylor.variable(c,c.mpf('.5'),3)
        p0 = 1+z*z/c.mpf(10)
        u0 = 2+z/c.mpf(5)
        p1, u1 = constant(c,c.mpf('.1'),3), constant(c,c.mpf('-.05'),3)
        phi, u = [p0,p1], [u0,u1]
        calc = SimpleNamespace(ctx=c,precision=80,phi=phi,u=u,eps=c.mpf('.1'),lam=c.mpf(10),
            z=z,delta=c.mpf('.01'),ell=constant(c,0,3),S=constant(c,'.5',3),
            p0=constant(c,-3,3),degree=1,state_hash='polynomial-fixture',center_family='.5')
        calc.phi_squared = radial_product(phi,phi)
        calc.phi_u = radial_product(phi,u)
        calc.u_squared = radial_product(u,u)
        calc.core_state = lambda s: CandidateComparisonJets.core_state(calc,s)
        calc.radial = CandidateComparisonJets.radial
        calc.tail = dict(target_met=False,rows=[dict(scaled_radial_order=0,axial_order=k,
            Phi_tail=c.mpf(['.011','.002','0','0'][k]),
            Psi_tail=c.mpf(['.2','.03','0','0'][k])) for k in range(4)])
        out = enclose(calc,'4')
        # The exact full fields add (.01+.002 Z) and (.02-.003 Z),
        # independent of radius. Integrate their linear radial polynomials
        # in closed form rather than using the production integral helper.
        a = p0+c.mpf('.01')+z*c.mpf('.002')
        b = u0+c.mpf('.02')-z*c.mpf('.003')
        r = c.mpf(4)
        expected = dict(phi=a+p1*r,U=b+u1*r,
            theta=a*r*r+p1*(2*r**3/3),z=b*r+u1*r*r/2,
            theta_z=a*b*r*r+(a*u1+p1*b)*(2*r**3/3)+p1*u1*r**4/2,
            p=a*a*r+a*p1*r*r+p1*p1*r**3/3,
            u_squared=b*b*r+b*u1*r*r+u1*u1*r**3/3,
            weighted_phi_squared=a*a*r*r/2+a*p1*(2*r**3/3)+p1*p1*r**4/4)
        count = 0
        for name, value in expected.items():
            for k in range(4):
                lo,hi = endpoints(out['full_core_enclosures'][name][k])
                el,eh = endpoints(value[k])
                if not lo <= el <= eh <= hi:
                    raise AssertionError((name,k,'closed-form field outside tail enclosure'))
                count += 1
        report = dict(closed_form_polynomial_coefficients_contained=count,
            axial_orders=[0,1,2,3],nonzero_axial_tail_derivatives=True,
            Psi_to_physical_U_factor='1/Lambda',passed=True,
            fixture_only=True,production_source_errors_certified=False,
            input_hashes={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                for name in ('lei_ren_part1_paper_interval_core_inlet_bound.py',
                             'lei_ren_part1_paper_interval_core_inlet_bound_check.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Independent closed-form moment enclosure checks:',count,flush=True)
        return report


if __name__ == '__main__':
    run()
