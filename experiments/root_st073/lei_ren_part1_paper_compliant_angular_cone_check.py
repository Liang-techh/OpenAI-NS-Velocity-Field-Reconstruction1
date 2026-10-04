"""Recompute the original angular cone from stable native continuous history bounds."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_angular_cone import (
    CompliantAngularStressC3,CompliantGlobalPhysicalAssembly,DOMAIN,TAIL_DOMAIN,
    angular_cone_identities,angular_cone_source_bridge,whole_angular_bounds)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def run():
    with mp.workdps(300):
        name=PREFIX+'angular_cone.json'; record=json.loads((HERE/name).read_bytes()); hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Angular cone source changed: '+source)
        stress=CompliantAngularStressC3(); assembly=CompliantGlobalPhysicalAssembly(); c=stress.ctx
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Angular cone actual family differs')
        if angular_cone_identities()!=record['exact_identities']:raise ValueError('Angular cone source identities changed')
        if angular_cone_source_bridge(stress)!=record['source_bridge']:raise ValueError('Angular cone native W/log bindings changed')
        bounds=whole_angular_bounds(stress,assembly)
        if encode(pack(bounds))!=record['bounds']:raise ValueError('Current continuous whole-entry bounds changed')
        for label,value in bounds['positive_margins'].items():
            if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Angular margin failed: '+label)
        if endpoints(bounds['log_directional_term_upper'])[1]>=endpoints(c.ln(2))[0]:raise ArithmeticError('Angular directional inequality failed')
        if bounds['raw_whole_theta_box_used_as_positivity_proof'] or bounds['phase_grid_sampling_used_as_proof']:
            raise ValueError('Angular proof must use continuous source history bounds')
        if not bounds['continuous_native_cumulative_history_bound_used']:raise ValueError('Continuous native W history bound missing')
        if not bounds['whole_s_Z_domain_and_original_support_crossings_covered']:raise ValueError('Original angular support/crossing domain missing')
        if endpoints(bounds['actual_kappa_minus2_enclosure'])[0]<=0 or endpoints(bounds['actual_kappa_minus2_upper'])[1]>=2:
            raise ArithmeticError('Actual whole angular shear range failed')
        if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:raise ValueError('Angular cone domain changed')
        for flag in ('angular_cone_certified','actual_angular_theta_stress_positive',
                     'actual_angular_source_shear_strictly_negative','actual_angular_directional_cone_uniform_margin',
                     'angular_entry_power_exit_waiting_collar_two_vector_cone_certified',
                     'same_source_angular_entry_similarity_and_physical_joins_consumed',
                     'fixed_positive_viscosity_two_vector_cone_transfer_verified'):
            if not record[flag]:raise ValueError('Angular cone fact missing: '+flag)
        for flag in ('angular_regional_remainder_exact_zero','preceding_power_physical_stress_join_verified',
                     'completed_full_tensor_cone_certified','whole_outer_cone_certified','global_admissible_stress_lift_constructed',
                     'independently_bounded_global_flat_remainder','physical_energy_integral_certified','temporal_recursion'):
            if record[flag]:raise ValueError('Angular cone scope overclaimed: '+flag)
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    angular_cone_certified=True,domain=DOMAIN,tail_domain=TAIL_DOMAIN,
                    positive_source_margins_checked=len(bounds['positive_margins']),
                    whole_s_Z_source_bounds_recomputed=True,phase_grid_sampling_used_as_proof=False,
                    continuous_native_W_history_uniform_bound_verified=True,
                    actual_variable_angular_kappa_minus2_range_verified=True,
                    actual_source_shear_sign_uses_exact_positive_S=True,
                    actual_Bmax_is_original_source_inlet_verified=True,
                    angular_entry_power_exit_waiting_collar_two_vector_cone_certified=True,
                    same_source_angular_entry_similarity_and_physical_joins_consumed=True,
                    fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                    angular_regional_remainder_exact_zero=False,preceding_power_physical_stress_join_verified=False,
                    completed_full_tensor_cone_certified=False,whole_outer_cone_certified=False,
                    global_admissible_stress_lift_constructed=False,independently_bounded_global_flat_remainder=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS whole original angular cone and entry/power/exit/waiting/collar composition; upstream/global pending',flush=True)
        return result


if __name__=='__main__':run()
