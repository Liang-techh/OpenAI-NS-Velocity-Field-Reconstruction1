"""Independent RK diagnostic for actual exit and analytic axial tangents."""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_interval_comparison_enclosure_check import fixture
from lei_ren_part1_paper_interval_comparison_cells import build_cells
from lei_ren_part1_paper_interval_exit_bridge_enclosure import integrate_exit,multiplier
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet,radial_product
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def source(center):
    calc = fixture()
    c = calc.ctx
    calc.z = IntervalTaylor.variable(c,c.mpf(center),3)
    calc.center_family = center
    calc.phi = [1+calc.z*calc.z/10]
    calc.u = [2+calc.z/5]
    calc.phi_squared = radial_product(calc.phi,calc.phi)
    calc.phi_u = radial_product(calc.phi,calc.u)
    calc.u_squared = radial_product(calc.u,calc.u)
    calc.initial_phi = calc.phi[0]
    calc.F0 = c.sqrt(c.mpf('.5'))
    # This known polynomial has zero coefficients and zero tail beyond
    # its constant radial row, including at the required truncation124.
    calc.degree = 124
    return calc


def reference(calc,epsilon,steps=512):
    c = calc.ctx
    hb = c.mpf('.005')
    initial = calc.core_state(c.mpf(4))
    phi_a = initial['phi'].truncate(1)
    initial_state = {key:value.truncate(1) for key,value in initial.items() if key not in ('phi','U')}
    initial_state.update(g=phi_a*0,U=initial['U'].truncate(1))
    def rhs(y,state):
        s = 4*c.exp(y)
        bar = calc.core_state(s)
        zero = bar['phi']*0
        transfer = normalized_inlet(calc.phi,calc.u,bar,calc.S,calc.ell,calc.p0,
            calc.lam,s,calc.z,calc.delta,endpoints_override=dict(
                phi_exit=bar['phi'],u_exit=bar['U'],phi_s=zero,u_s=zero))
        chi = multiplier(c,y,hb,epsilon)
        A = -transfer['ratio'].truncate(1)*chi/2
        B = -(phi_a/bar['phi'].truncate(1))*transfer['iz'].truncate(1)*c.sqrt(calc.eps)*c.sqrt(s*calc.eps/2)*chi
        eg = state['g'].exp()
        f,u = phi_a*eg,state['U']
        return dict(g=A,U=eg*B,theta=f*(2*s*s),z=u*s,theta_z=f*u*(2*s*s),
            p=f*f*s,u_squared=u*u*s,weighted_phi_squared=f*f*s*s)
    def shifted(state,rhs,weight):
        return {key:value+rhs[key]*weight for key,value in state.items()}
    state = initial_state
    h = 2*hb/steps
    for n in range(steps):
        y = h*n
        k1 = rhs(y,state)
        k2 = rhs(y+h/2,shifted(state,k1,h/2))
        k3 = rhs(y+h/2,shifted(state,k2,h/2))
        k4 = rhs(y+h,shifted(state,k3,h))
        state = {key:value+(k1[key]+k2[key]*2+k3[key]*2+k4[key])*h/6 for key,value in state.items()}
    return dict(phi=phi_a*state['g'].exp(),**state)


def run():
    with mp.workdps(120):
        family = source(['.49','.51'])
        comparison = build_cells(family,16)
        bound = integrate_exit(family,comparison,16,exit_epsilon='.01')
        scalar = source('.5')
        ref = reference(scalar,scalar.ctx.mpf('.01'))
        full_bound = dict(g=bound['log_F_over_Fa'],**bound['normalized_exit_state'])
        count = 0
        for name,jet in ref.items():
            for k in range(2):
                lo,hi = endpoints(full_bound[name][k])
                el,eh = endpoints(jet[k])
                if not lo <= el <= eh <= hi:
                    raise AssertionError((name,k,'independent exit RK diagnostic outside enclosure'))
                count += 1
        if not bound['physical_F_strictly_positive'] or bound['midpoint_projection_used']:
            raise AssertionError('Exit field positivity or interval preservation failed')
        report = dict(passed=True,fixture_only=True,scalar_center='.5',family=['.49','.51'],
            caller_exit_epsilon='.01',comparison_cells=32,independent_RK_steps=512,
            exit_and_first_axial_derivative_coefficient_checks=count,
            numerical_reference_is_accuracy_diagnostic_not_proof=True,
            source_hashes={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()
                for name in ('lei_ren_part1_paper_interval_exit_bridge_enclosure.py',
                    'lei_ren_part1_paper_interval_exit_bridge_enclosure_check.py',
                    'lei_ren_part1_paper_interval_comparison_cells.py')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Independent actual exit RK/tangent diagnostic:',count,'coefficients contained',flush=True)
        return report


if __name__ == '__main__':
    run()
