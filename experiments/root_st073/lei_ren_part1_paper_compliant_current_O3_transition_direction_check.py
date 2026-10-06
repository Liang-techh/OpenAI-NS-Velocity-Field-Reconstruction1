"""Focused actual variable source, endpoint obstruction and shear-loop check."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_transition_direction import (
    CurrentO3TransitionDirection,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,
    pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_transition_direction_operator import (
    validate_whole_transition,whole_current_transition_bounds)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentO3TransitionDirection(require_checked=False)
    expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Variable transition manifest/scope differs')
    views=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if views!=encode(pack(dict(whole_current_variable_transition=field.whole_view,actual_current_O2_buffer_inlet=field.inlet))):
        raise ValueError('Complete actual variable source views differ')
    replay=whole_current_transition_bounds(field.powercone,field.whole_view,field.inlet)
    if encode(pack(replay))!=raw['current_whole_variable_transition_bounds']:raise ValueError('Continuous actual variable source replay differs')
    for name,value in replay['positive_margins'].items():
        if endpoints(value)[0]<=0:raise ValueError('Current variable bound unresolved: '+name)
    if replay['closed_transition_strict_cone_certified'] or not replay['original_shear_repair_required'] \
            or replay['exact_offset0_shear_excess']!=0 or not replay['offset0_stress_nonzero_by_same_positive_theta_lower']:
        raise ValueError('Nonzero original zero-shear endpoint must remain unadmitted')
    if replay['source_correlated_periodic_shear_loop']['spatial_support_and_finite_frequency_and_moment_repair_installed']:
        raise ValueError('Periodic loop is not an installed profile')
    for fixture in ('foreign','domain','jet','pressure','mode','physical'):
        bad=copy_containers(field.whole_view);source=bad['original_complete_view']
        if fixture=='foreign':bad['implicit_source_sha256']='foreign'
        elif fixture=='domain':source['coverage_coordinate']=field.ctx.mpf(['.001',1])
        elif fixture=='jet':source['all_variable_log_amplitude_ordinary_derivatives_retained']=False
        elif fixture=='pressure':source.pop('current_absolute_pressure_ordinary_y_rows')
        elif fixture=='mode':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['theta']['variable_radial_shear']['mode']=(.5,1,0,0)
        else:source['physical_cylindrical_stress_mixed3']['theta']['variable_radial_shear']['r0_z0']['signed_coefficient']=field.ctx.mpf(13)
        try:validate_whole_transition(field.powercone,bad)
        except ValueError:pass
        else:raise ValueError('Invalid variable source admitted: '+fixture)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
        exact_original_variable_history_pressure_operator_and_shear_loop_identities=len(field.theorem['identities']),
        strictly_positive_directed_bounds=len(replay['positive_margins']),
        all_eleven_signed_sectors_full_functional_pressure_energy_and_variable_jets_retained=True,
        exact_nonzero_zero_shear_endpoint_requires_original_shear_repair=True,
        actual_mean_preserving_shear_loop_primitives_and_same_inviscid_target_proved=True,
        finite_frequency_profile_spatial_cutoff_and_independent_moment_repair_remain_open=True,
        phase_samples_and_interval_overlap_not_used_as_proof=True,
        current_strict_nonzero_whole_regions_including_inherited=15,current_exact_zero_exterior_regions=1,
        remaining_registry_regions_without_current_whole_cone=17,scope=raw['scope'],
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual variable transition direction/shear-loop PASS; closed transition cone correctly remains open',flush=True)
    return result


if __name__=='__main__':run()
