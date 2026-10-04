"""Recompute the whole original power cone and its correlated history proof."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_outer_power_cone import (
    CompliantOuterPowerStressC3,CompliantGlobalPhysicalAssembly,DOMAIN,TAIL_DOMAIN,
    outer_power_cone_identities,outer_power_cone_source_bridge,whole_outer_power_bounds)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_correlated_flatten_history_fixture():
    """Moderate direct native integrals; samples do not establish the cone.

    Xv=.02 is below equilibrium. The actual continuous original-domain proof
    is the source-bound analytic inequality checked separately below.
    """
    with mp.workdps(75):
        a,mu,xv=map(mp.mpf,('.1','.2','.02'));r=1-mu;k=1-a;b=(1-2*a)/2
        floor=mp.exp(-100*r)*(b*mp.expm1(50*r)/2-k)/r
        if floor<=0:raise ArithmeticError('Correlated fixture parameter condition failed')
        def sigma(v):
            x=v/100
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            odds=1/(1-x)**2-1/x**2
            return 1/(1+mp.exp(-odds))
        errors=[];gaps=[]
        for z in map(mp.mpf,('0','.5','1')):
            f=(1+z*z)/2;j=2*z*z/(1+z*z);rho=mp.log(f)
            weight=lambda v:mp.exp(-r*(100-v))
            I=mp.quad(lambda v:weight(v)*mp.exp(rho*sigma(v)),[0,25,50,75,100])
            Iz=mp.quad(lambda v:weight(v)*mp.exp(rho*sigma(v))*sigma(v)*2*z/(1+z*z),[0,25,50,75,100])
            numerator=I+xv*mp.exp(-100*r)
            Xf=numerator/f;Xfz=(Iz*f-numerator*z)/(f*f)
            direct=k*(Xf-1/r)-b*z*Xfz
            correlated=mp.exp(-100*r)*((k+b*j)*xv/f-k/r)
            correlated+=mp.quad(lambda v:weight(v)*(k*(mp.exp(-(1-sigma(v))*rho)-1)
                +b*j*(1-sigma(v))*mp.exp(-(1-sigma(v))*rho)),[0,25,50,75,100])
            error=abs(direct-correlated)
            if error>mp.mpf('1e-60') or direct<floor:raise ArithmeticError('Independent correlated flatten history check failed')
            errors.append(str(error));gaps.append(str(direct-floor))
        return dict(passed=True,direct_native_integral_identity_checks=3,positive_lemma_comparisons=3,
            Xv_below_equilibrium=True,identity_errors=errors,positive_gaps=gaps,
            comparison_tolerance='1e-60',sampled_values_used_as_actual_cone_proof=False,
            actual_cone_proof_is_source_bound_continuous_inequality=True,fixture_is_global_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'outer_power_cone.json';record=json.loads((HERE/name).read_bytes());hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power cone source changed: '+source)
        stress=CompliantOuterPowerStressC3();assembly=CompliantGlobalPhysicalAssembly();c=stress.ctx
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Power cone actual family differs')
        if outer_power_cone_identities()!=record['exact_identities']:raise ValueError('Power cone source identities changed')
        if encode(pack(outer_power_cone_source_bridge(stress)))!=record['source_bridge']:raise ValueError('Power cone actual history/log bindings changed')
        for flag in ('actual_F_endpoint_equals_original_q_over_two_AST_verified',
            'actual_Fc_over_endpoint_F_equals_f_power_minus_omega_AST_verified',
            'original_source_Xint_equals_f_times_same_normalized_integral_bound'):
            if not record['source_bridge']['identities'][flag]:raise ValueError('Actual correlated-history normalization missing: '+flag)
        bounds=whole_outer_power_bounds(stress,assembly)
        if encode(pack(bounds))!=record['bounds']:raise ValueError('Current continuous whole-power bounds changed')
        for label,value in bounds['positive_margins'].items():
            if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Power cone margin failed: '+label)
        if endpoints(bounds['log_directional_term_upper'])[1]>=endpoints(c.ln(2))[0]:raise ArithmeticError('Power directional inequality failed')
        if bounds['raw_whole_theta_box_used_as_positivity_proof'] or bounds['phase_grid_sampling_used_as_proof']:
            raise ValueError('Power cone proof must use continuous correlated history')
        for flag in ('continuous_native_cumulative_history_bound_used','correlated_flatten_history_uniform_lower_used',
            'whole_original_phase_Z_domain_covered','interval_S_lower_zero_not_used_for_strict_source_sign'):
            if not bounds[flag]:raise ValueError('Actual whole-power proof fact missing: '+flag)
        if endpoints(bounds['actual_kappa_minus2_enclosure'])[0]<=0 or endpoints(bounds['actual_kappa_minus2_upper'])[1]>=2:
            raise ArithmeticError('Whole original constant shear range failed')
        if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:raise ValueError('Power cone original domain changed')
        for flag in ('outer_power_cone_certified','actual_outer_power_theta_stress_positive',
            'actual_outer_power_source_shear_strictly_negative','actual_outer_power_directional_cone_uniform_margin',
            'outer_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
            'same_source_power_angular_similarity_and_physical_joins_consumed',
            'fixed_positive_viscosity_two_vector_cone_transfer_verified'):
            if not record[flag]:raise ValueError('Power cone admission missing: '+flag)
        for flag in ('outer_power_regional_remainder_exact_zero','left_flatten_physical_stress_join_verified',
            'completed_full_tensor_cone_certified','whole_outer_cone_certified','global_admissible_stress_lift_constructed',
            'independently_bounded_global_flat_remainder','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Power cone scope overclaimed: '+flag)
        fixture=independent_correlated_flatten_history_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
            outer_power_cone_certified=True,domain=DOMAIN,tail_domain=TAIL_DOMAIN,
            positive_source_margins_checked=len(bounds['positive_margins']),whole_original_phase_Z_source_bounds_recomputed=True,
            continuous_native_W_history_uniform_bound_verified=True,correlated_original_flatten_Hf_positive_lemma_verified=True,
            phase_grid_sampling_used_as_proof=False,actual_constant_power_kappa_minus2_range_verified=True,
            actual_source_shear_sign_uses_exact_positive_S=True,actual_Bmax_is_original_source_inlet_verified=True,
            outer_power_angular_entry_exit_waiting_collar_two_vector_cone_certified=True,
            same_source_power_angular_similarity_and_physical_joins_consumed=True,
            fixed_positive_viscosity_two_vector_cone_transfer_verified=True,independent_correlated_flatten_history_fixture=fixture,
            outer_power_regional_remainder_exact_zero=False,left_flatten_physical_stress_join_verified=False,
            completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,global_admissible_stress_lift_constructed=False,
            independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS whole original preceding-power cone and angular tail composition; flatten/global/recursion pending',flush=True)
        return result


if __name__=='__main__':run()
