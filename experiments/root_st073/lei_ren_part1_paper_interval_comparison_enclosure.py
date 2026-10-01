"""Directed integral enclosure of the switched comparison ODE.

Phi' = g(y,Z) Phi and U' = j(y,Z) have explicit integrating-factor
solutions. Interval cell integrals bound these solutions and all moment
integrals, including radial core tails. Axial order is two: the mixed C3
core budget supplies one radial plus two axial derivatives, not C4.
Bounds remain relative to accepted construction parameters.
"""
import hashlib
import json
import math
import operator
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_core_inlet_bound import absolute_jet_envelope
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def tail_jet(calc, name, radial_order=0, scale=1):
    c = calc.ctx
    rows = {row['axial_order']:row[name] for row in calc.tail['rows']
            if row['scaled_radial_order'] == radial_order}
    values = []
    for k in range(3):
        bound = rows[k]*c.mpf(scale)/math.factorial(k)
        upper = endpoints(bound)[1]
        values.append(c.mpf([-upper,upper]))
    return IntervalTaylor(c,values)


def core_box(calc, s):
    """Full core/moment jets over a radius interval through scaled 4.1."""
    c = calc.ctx
    lo,hi = endpoints(s)
    if lo <= 0 or hi > mp.mpf('4.1'):
        raise ValueError('Core box outside (0,4.1]')
    finite = {name:jet.truncate(2) for name,jet in calc.core_state(s).items()}
    ep,eu = tail_jet(calc,'Phi_tail'),tail_jet(calc,'Psi_tail',scale=calc.eps)
    bp = absolute_jet_envelope(calc.phi,c.mpf(hi)).truncate(2)
    bu = absolute_jet_envelope(calc.u,c.mpf(hi)).truncate(2)
    pp = bp*ep*2+ep*ep
    pu = bp*eu+bu*ep+ep*eu
    uu = bu*eu*2+eu*eu
    return dict(phi=finite['phi']+ep,U=finite['U']+eu,
        theta=finite['theta']+ep*s*s,z=finite['z']+eu*s,
        theta_z=finite['theta_z']+pu*s*s,p=finite['p']+pp*s,
        u_squared=finite['u_squared']+uu*s,
        weighted_phi_squared=finite['weighted_phi_squared']+pp*s*s/2)


def alpha_box(c, y, hb):
    """Monotone decreasing smooth cutoff, with exact endpoint limits."""
    qlo,qhi = endpoints((y-hb)/hb)
    def point(q):
        if q <= 0:
            return c.mpf(1)
        if q >= 1:
            return c.mpf(0)
        x = c.mpf(q)
        left,right = c.exp(-1/(x*x)),c.exp(-1/((1-x)*(1-x)))
        return right/(left+right)
    return c.mpf([endpoints(point(qhi))[0],endpoints(point(qlo))[1]])


def integrate(calc, cells=32):
    cells = operator.index(cells)
    if cells < 1:
        raise ValueError('At least one integral cell required')
    if not calc.tail['target_met']:
        raise ValueError('Controlled comparison requires passing analytic core tail')
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        hb = c.mpf('.005')
        step = hb/cells
        start = core_box(calc,4*c.exp(hb))
        zero = start['phi']*0
        G,J = zero,zero
        moments = {key:value for key,value in start.items() if key not in ('phi','U')}
        ep = tail_jet(calc,'Phi_tail')
        dp = tail_jet(calc,'Phi_tail',radial_order=1)
        du = tail_jet(calc,'Psi_tail',radial_order=1,scale=calc.eps)
        minimum_source_phi = None
        for n in range(cells):
            left = hb*(1+c.mpf(n)/cells)
            right = hb*(1+c.mpf(n+1)/cells)
            y = c.mpf([endpoints(left)[0],endpoints(right)[1]])
            s = 4*c.exp(y)
            source_phi = calc.radial(calc.phi,s).truncate(2)+ep
            lower = endpoints(source_phi[0])[0]
            if lower <= 0:
                raise ValueError('Source slope denominator crosses zero; refine axial family')
            minimum_source_phi = lower if minimum_source_phi is None else min(lower,minimum_source_phi)
            alpha = alpha_box(c,y,hb)
            g = (calc.radial(calc.phi,s,1).truncate(2)+dp)*s*alpha/source_phi
            j = (calc.radial(calc.u,s,1).truncate(2)+du)*s*alpha
            # A partial-cell integral is enclosed by [0, step_upper] times
            # the complete cell range, for each ordinary Taylor coefficient.
            partial = c.mpf([0,endpoints(step)[1]])
            f = start['phi']*(G+g*partial).exp()
            u = start['U']+J+j*partial
            rhs = dict(theta=f*(2*s*s),z=u*s,theta_z=f*u*(2*s*s),
                p=f*f*s,u_squared=u*u*s,weighted_phi_squared=f*f*s*s)
            moments = {key:value+rhs[key]*step for key,value in moments.items()}
            G,J = G+g*step,J+j*step
        end = dict(phi=start['phi']*G.exp(),U=start['U']+J,**moments)
        if endpoints(end['phi'][0])[0] <= 0:
            raise ValueError('Controlled comparison endpoint loses positivity')
        s = 4*c.exp(2*hb)
        transfer = normalized_inlet(calc.phi,calc.u,end,calc.S.truncate(2),calc.ell.truncate(2),
            calc.p0.truncate(2),calc.lam,s,calc.z.truncate(2),calc.delta,
            endpoints_override=dict(phi_exit=end['phi'],u_exit=end['U'],phi_s=zero,u_s=zero))
        iz = transfer['iz']*c.sqrt(calc.eps)
        initial_phi = calc.initial_phi.truncate(2)+ep
        A = -transfer['ratio']/2
        B = -(initial_phi/end['phi'])*iz*c.sqrt(s*calc.eps/2)
        return dict(center_family=calc.center_family,state_sha256=calc.state_hash,
            completed_degree=calc.degree,cells=cells,y_interval=['.005','.01'],
            retained_axial_order=2,driver_axial_order=transfer['ratio'].order,
            endpoint_state=end,endpoint_scaled_radius=s,
            D=transfer['ratio'],I_z=iz,pressure=transfer['pressure'],
            initial_phi=initial_phi,endpoint_driver_A=A,endpoint_driver_B=B,
            endpoint_driver_chi='1',
            minimum_source_phi_lower=minimum_source_phi,
            Phi_log_increment=G,U_increment=J,
            cell_integrals_enclosed=True,analytic_core_tail_propagated=True,
            comparison_discretization_enclosed=True,
            method='integrating factor and directed cell range integrals',
            relative_to_accepted_construction_data=True,
            construction_parameter_errors_enclosed=False,
            exit_ODE_error_enclosed=False,whole_axis_transition_generated=False,
            terminal_five_moment_closure=False,stress_derivatives_supplied=False,
            temporal_recursion=False)


def pack(value):
    if isinstance(value,IntervalTaylor):
        return list(value.coefficients)
    if isinstance(value,dict):
        return {key:pack(item) for key,item in value.items()}
    return value


def frozen_extension(calc, result, scaled_radius):
    """Exact moment continuation of the comparison field after alpha=0.

    This is the frozen comparison, not the physical exit continuation ODE.
    It uses its already-enclosed switched endpoint and no core extrapolation.
    """
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        s = c.mpf(scaled_radius)
        sb = result['endpoint_scaled_radius']
        if endpoints(s)[0] < endpoints(sb)[1]:
            raise ValueError('Frozen extension requires radius beyond switched endpoint')
        state = dict(result['endpoint_state'])
        f,u = state['phi'],state['U']
        ds,ds2 = s-sb,s*s-sb*sb
        for name,increment in dict(theta=f*ds2,z=u*ds,theta_z=f*u*ds2,
            p=f*f*ds,u_squared=u*u*ds,weighted_phi_squared=f*f*ds2/2).items():
            state[name] = state[name]+increment
        zero = f*0
        transfer = normalized_inlet([f],[u],state,calc.S.truncate(2),calc.ell.truncate(2),
            calc.p0.truncate(2),calc.lam,s,calc.z.truncate(2),calc.delta,
            endpoints_override=dict(phi_exit=f,u_exit=u,phi_s=zero,u_s=zero))
        iz = transfer['iz']*c.sqrt(calc.eps)
        return dict(scaled_radius=s,state=state,D=transfer['ratio'],I_z=iz,
            pressure=transfer['pressure'],
            endpoint_driver_A=-transfer['ratio']/2,
            endpoint_driver_B=-(result['initial_phi']/f)*iz*c.sqrt(s*calc.eps/2),
            endpoint_driver_chi='1',frozen_comparison_extension=True,
            old_center_tensor_extrapolated=False,
            analytic_core_radial_domain_extrapolated=False,
            physical_exit_ODE_solved=False,stress_cone_certified=False)


def run(cells=32):
    calc = IntervalComparisonJets(4)
    result = integrate(calc,cells)
    result['frozen_comparison_at_physical_R110'] = frozen_extension(calc,result,calc.lam*110)
    here = Path(__file__).parent
    result['input_hashes'] = {name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
        ('lei_ren_part1_paper_interval_comparison_enclosure.py',
         'lei_ren_part1_paper_interval_comparison_jets.py',
         'lei_ren_part1_paper_interval_core_inlet_bound.py',
         'lei_ren_part1_paper_candidate_shared_inlet.py',
         'lei_ren_part1_paper_interval_taylor.py')}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Directed comparison integral enclosure:',cells,'cells; axial field order2, driver order',
          result['driver_axial_order'],flush=True)
    return result


if __name__ == '__main__':
    run()
