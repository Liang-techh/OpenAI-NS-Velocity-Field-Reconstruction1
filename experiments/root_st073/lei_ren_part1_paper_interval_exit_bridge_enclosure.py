"""Initial physical exit ODE with directed comparison/solution errors.

The actual exit uses g'=A and U'=exp(g) B, not the frozen comparison field.
Cell ranges integrate g first, then U and the six cumulative moments.
All fields carry one analytic axial derivative over the saved center family.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_comparison_enclosure import core_box,alpha_box,pack
from lei_ren_part1_paper_interval_comparison_cells import build_cells
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_candidate_general_center_factory import PARAMETERS
from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def shear_epsilon(calc):
    c = calc.ctx
    axis = uniform_axis_jets(c,radius=1,j=PARAMETERS['j'],Lambda=PARAMETERS['Lambda'],
        logC=PARAMETERS['logC'],delta=PARAMETERS['delta'],length=4)
    upper = endpoints(axis['G_absolute_upper'])[1]
    logeps = -calc.lam*(2*c.mpf(upper)+1)
    epsilon = c.exp(logeps)
    if not 0 < endpoints(epsilon)[0] <= endpoints(epsilon)[1] < 1:
        raise ValueError('Exit shear epsilon must be strictly between zero and one')
    return logeps,epsilon


def multiplier(c,y,hb,epsilon):
    # sigma(x)+sigma(1-x)=1, so this is exactly the source multiplier.
    return epsilon+(1-epsilon)*alpha_box(c,y+hb,hb)


def integrate_exit(calc,comparison,cells_per_half=16,exit_epsilon=None):
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        if exit_epsilon is None:
            logeps,epsilon = shear_epsilon(calc)
            policy = 'same conservative uniform analytic G upper bound as nominal exit'
        else:
            epsilon = c.mpf(exit_epsilon)
            if not 0 < endpoints(epsilon)[0] <= endpoints(epsilon)[1] < 1:
                raise ValueError('Supplied exit epsilon must be strictly between zero and one')
            logeps = c.ln(epsilon)
            policy = 'caller supplied; no Section9 contraction certificate'
        if len(comparison['cells']) != 2*cells_per_half:
            raise ValueError('Comparison cell partition does not match requested half count')
        hb = c.mpf('.005')
        initial = {key:jet.truncate(1) for key,jet in core_box(calc,c.mpf(4)).items()}
        phi_a,u_a = initial['phi'],initial['U']
        zero = phi_a*0
        g,J = zero,zero
        late_g = zero
        moments = {key:value for key,value in initial.items() if key not in ('phi','U')}
        for index,packet in enumerate(comparison['cells']):
            y,step,s = packet['y'],packet['step'],packet['scaled_radius']
            bar_phi = packet['state']['phi'].truncate(1)
            if endpoints(bar_phi[0])[0] <= 0:
                raise ValueError('Comparison angular denominator crosses zero')
            chi = multiplier(c,y,hb,epsilon)
            A = -packet['D'].truncate(1)*chi/2
            B = -(phi_a/bar_phi)*packet['I_z'].truncate(1)*c.sqrt(s*calc.eps/2)*chi
            partial = c.mpf([0,endpoints(step)[1]])
            g_range = g+A*partial
            eg = g_range.exp()
            u_rhs = eg*B
            u_range = u_a+J+u_rhs*partial
            phi_range = phi_a*eg
            rhs = dict(theta=phi_range*(2*s*s),z=u_range*s,
                theta_z=phi_range*u_range*(2*s*s),p=phi_range*phi_range*s,
                u_squared=u_range*u_range*s,weighted_phi_squared=phi_range*phi_range*s*s)
            moments = {key:value+rhs[key]*step for key,value in moments.items()}
            g,J = g+A*step,J+u_rhs*step
            if index >= cells_per_half:
                late_g += A*step
        phi,u = phi_a*g.exp(),u_a+J
        state = dict(phi=phi,U=u,**moments)
        endpoint = comparison['endpoint']
        s = endpoint['scaled_radius']
        R = s*calc.eps
        A_end = -endpoint['D'].truncate(1)*epsilon/2
        B_end = -(phi_a/endpoint['state']['phi'].truncate(1))*endpoint['I_z'].truncate(1)*c.sqrt(s*calc.eps/2)*epsilon
        U_y = g.exp()*B_end
        amplitude = IntervalTaylor(c,[calc.F0,calc.F0*calc.ell[0]])
        F = amplitude*phi
        S = calc.S.truncate(1)
        physical_moments = dict(theta=amplitude*moments['theta']*calc.eps**2,
            z=moments['z']*calc.eps,theta_z=amplitude*moments['theta_z']*calc.eps**2,
            z_theta=moments['u_squared']*calc.eps-S*moments['weighted_phi_squared']*calc.eps**2,
            p=S*moments['p']*calc.eps)
        P = calc.p0.truncate(1)+physical_moments['p']
        z = calc.z[0]
        mz = physical_moments['z']
        Ur = (2*z*R*u[0]-(1-calc.delta)*z*mz[0]-(1-z*z)*mz[1])/((1-calc.delta*z*z)*c.sqrt(2*R))
        if endpoints(F[0])[0] <= 0:
            raise ValueError('Actual exit angular field loses strict positivity')
        return dict(center_family=calc.center_family,state_sha256=calc.state_hash,
            accepted_schedule_sha256=ACCEPTED_SHA,completed_core_degree=calc.degree,
            cells_per_half=cells_per_half,endpoint_y='.01',physical_R=R,
            radial_scale_epsilon=calc.eps,exit_shear_epsilon=epsilon,log_exit_shear_epsilon=logeps,
            epsilon_policy=policy,
            log_F_over_Fa=g,U_increment=J,late_shear_log_increment=late_g,
            normalized_exit_state=state,F=F,Uz=u,Ur_value=Ur,P=P,
            physical_moments=physical_moments,
            actual_terminal_g_y=A_end,actual_terminal_U_y=U_y,
            actual_F_y=F*A_end,actual_F_R=F*A_end/R,actual_Uz_R=U_y/R,
            retained_axial_order=1,Ur_axial_derivative_enclosed=False,
            physical_exit_ODE_error_enclosed=True,comparison_error_propagated=True,
            analytic_core_tail_propagated=True,midpoint_projection_used=False,
            prescribed_terminal_shear_preserved=bool(endpoints(epsilon)[0]>0),
            physical_F_strictly_positive=True,
            original_construction_parameter_errors_enclosed=False,
            exit_shear_policy_section9_contraction_constants_certified=False,
            stress_cone_certified=False,post_collar_exit_R100_enclosed=False,
            terminal_five_moment_closure=False,whole_axis_transition=False,
            temporal_recursion=False)


def run(cells_per_half=16):
    calc = IntervalComparisonJets(4)
    comparison = build_cells(calc,cells_per_half=cells_per_half)
    result = integrate_exit(calc,comparison,cells_per_half)
    here = Path(__file__).parent
    result['input_hashes'] = {name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
        ('lei_ren_part1_paper_interval_exit_bridge_enclosure.py',
         'lei_ren_part1_paper_interval_comparison_cells.py',
         'lei_ren_part1_paper_interval_comparison_enclosure.py',
         'lei_ren_part1_paper_interval_comparison_jets.py',
         'lei_ren_part1_paper_uniform_axis_jets.py')}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual physical exit interval bridge complete; axial order1; midpoint projection unused',flush=True)
    return result


if __name__ == '__main__':
    run()
