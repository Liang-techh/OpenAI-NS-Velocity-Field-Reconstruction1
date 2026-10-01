"""Whole-interval actual exit continuation to physical R=100.

The comparison field is frozen, but the actual exit retains its nonzero
implicit shear. Analytic radial weights avoid enormous log-radius grids.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_comparison import LogarithmicComparison
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_state,restore_jet,_pack
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_candidate_shared_inlet import diff

HERE=Path(__file__).parent

def run(target_R='100'):
    calc=LogarithmicComparison();c=calc.ctx
    names=('lei_ren_part1_paper_logarithmic_comparison.json','lei_ren_part1_paper_logarithmic_exit_bridge.json')
    comp,bridge=[json.loads((HERE/n).read_bytes()) for n in names]
    for record in (comp,bridge):
        if record['identity']!=calc.identity:raise ValueError('continuation source changed')
        for n,d in record['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('dependency changed: '+n)
    if not bridge['actual_exit_ODE_integral_enclosed'] or not bridge['shear_zero_lower_is_bound_not_selected_value'] or bridge['y_interval']!=['0','.01']:
        raise ValueError('bridge admission or coverage missing')
    if not comp['cell_integral_discretization_enclosed'] or not comp['core_analytic_tail_propagated'] or comp['y_interval']!=['.005','.01']:
        raise ValueError('comparison admission missing')
    with mp.workdps(c.dps+40):
        R=c.mpf(target_R)
        if endpoints(R)[0]<=0 or endpoints(R)[1]>100:raise ValueError('continuation target must be in (0,100]')
        R=c.mpf(target_R);s0=read_interval(c,comp['endpoint_scaled_radius']);st=R/calc.eps
        if endpoints(st)[0]<=endpoints(s0)[1]:raise ValueError('target precedes bridge')
        s=c.mpf([endpoints(s0)[0],endpoints(st)[1]])
        length=c.ln(st/s0);partial=c.mpf([0,endpoints(length)[1]])
        frozen=restore_state(c,comp['endpoint_state'],2)
        f,u=frozen['phi'],frozen['U']
        def positive_weight(v):
            lo,hi=endpoints(v);return c.mpf([max(mp.mpf(0),lo),max(mp.mpf(0),hi)])
        ds,ds2=positive_weight(s-s0),positive_weight(s*s-s0*s0)
        increments=dict(theta=f*ds2,z=u*ds,theta_z=f*u*ds2,p=f*f*ds,
            u_squared=u*u*ds,weighted_phi_squared=f*f*ds2/2)
        comparison=dict(phi=f,U=u,**{n:frozen[n]+value for n,value in increments.items()})
        base=calc.transfer(comparison,s)
        upper=read_interval(c,bridge['shear_upper_envelope']);gamma=c.mpf([0,endpoints(upper)[1]])
        A=base['unmodulated_driver_A'].truncate(1)*gamma
        B=base['unmodulated_driver_B'].truncate(1)*gamma
        initial=restore_state(c,bridge['normalized_actual_exit_state'],1)
        g0=restore_jet(c,bridge['log_F_over_Fa'],1)
        gbox=g0+A*partial
        eg=gbox.exp();U_rhs=eg*B
        phi_a=calc.initial_phi.truncate(1)
        phi=phi_a*eg;uv=initial['U']+U_rhs*partial
        if endpoints(phi[0])[0]<=0:raise ValueError('long exit normalized swirl positivity lost')
        def actual_packet(radius,y_increment):
            ds=positive_weight(radius-s0);ds2=positive_weight(radius*radius-s0*s0)
            # Every actual field/axial coefficient is bounded by phi/uv
            # over the full interval. Integrate weights exactly: integral
            # r dy=Delta r; integral 2r^2 dy=Delta r^2.
            weights=dict(theta=phi*ds2,z=uv*ds,theta_z=phi*uv*ds2,
                p=phi*phi*ds,u_squared=uv*uv*ds,weighted_phi_squared=phi*phi*ds2/2)
            moments={n:initial[n]+value for n,value in weights.items()}
            ge=g0+A*y_increment;phie=phi_a*ge.exp();ue=initial['U']+U_rhs*y_increment
            pressure=calc.p0_scaled.truncate(1)+calc.S_scaled.truncate(1)*moments['p']
            z=calc.z.truncate(1);L=1-z*z*calc.delta;d=1-z*z
            Ur=(z*ue*(2*radius)-z*moments['z']*(1-calc.delta)-d*diff(moments['z']))*c.sqrt(calc.eps)/(L*c.sqrt(2*radius))
            return dict(scaled_radius=radius,normalized_state=dict(phi=phie,U=ue,**moments),
                log_F_over_Fa=ge,pressure_scaled=pressure,physical_Ur=Ur,
                physical_moments_without_F0=dict(theta_over_F0=moments['theta']*calc.eps**2,
                    z=moments['z']*calc.eps,theta_z_over_F0=moments['theta_z']*calc.eps**2,
                    z_theta=moments['u_squared']*calc.eps-calc.S_scaled.truncate(1)*moments['weighted_phi_squared']))
        target=actual_packet(st,length)
        # The range packet uses independent positive radial integration
        # weights. Preserve their nonnegativity despite s-s0 dependency.
        radius_range=c.mpf([endpoints(s0)[0],endpoints(st)[1]])
        hashes={**calc.input_hashes,**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
            Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        result=dict(identity=calc.identity,input_hashes=hashes,
            implicit_source_sha256=calc.identity['implicit_source_sha256'],
            datum_enclosure_sha256=calc.identity['datum_enclosure_sha256'],
            target_physical_R=target_R,start_physical_R=s0*calc.eps,
            scaled_radius_range=radius_range,log_radius_length=length,
            actual_log_shape_range=gbox,actual_Phi_range=phi,actual_Uz_range=uv,
            actual_driver_A_range=A,actual_driver_B_range=B,
            continuation_log_shape_increment_bound=A*partial,
            continuation_axial_increment_bound=U_rhs*partial,
            target=target,shear_upper_envelope=upper,
            true_positive_shear_kept_implicit=True,whole_interval_actual_exit_field_enclosed=True,
            paper_hb_cstar_K_inverse100_relation_verified=False,full_Section9_parameter_admission=False,
            comparison_frozen_but_actual_exit_not_frozen=True,
            exact_radial_weight_integrals_used=True,retained_axial_order=1,
            finite_log_grid_not_used=True,physical_F0_not_materialized=True,
            terminal_reference_inlet_built=False,terminal_five_moment_repair_completed=False,
            whole_axis=False,admissible_stress_lift_constructed=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Md40 actual exit continuation enclosed to physical R='+target_R+'; C1 axial family')
    return result

if __name__=='__main__':run()
