"""Focused continuous O3 power source correlation and full-history check."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_O3_power_cone import (
    CurrentO3PowerCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,
    pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_O3_power_cone_operator import (
    validate_whole_O3_power_view,whole_current_O3_power_bounds)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentO3PowerCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current O3 power manifest/scope differs')
    stored=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if stored!=encode(pack(dict(whole_current_O3_power=field.whole_view,
        actual_O2_slope_endpoint=field.slope,actual_O3_power_phase0=field.initial))):
        raise ValueError('Complete current O3 power and canonical source views differ')
    replay=whole_current_O3_power_bounds(field.entrancecone,field.whole_view,field.slope,field.initial)
    if encode(pack(replay))!=raw['current_whole_O3_power_correlated_bounds']:
        raise ValueError('Current continuous O3 power source-bound replay differs')
    for name,value in replay['positive_margins'].items():
        if endpoints(value)[0]<=0:raise ValueError('Unresolved current O3 power bound: '+name)
    for fixture in ('foreign','domain','theta_memory','energy_pressure','history','full_pressure','mode','physical'):
        bad=copy_containers(field.whole_view);source=bad['original_complete_view']
        if fixture=='foreign':bad['implicit_source_sha256']='foreign'
        elif fixture=='domain':source['coverage_coordinate']=field.ctx.mpf(['.001',1])
        elif fixture=='theta_memory':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['theta'].pop('signed_original_memory')
        elif fixture=='energy_pressure':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['axial'].pop('full_energy_and_pressure')
        elif fixture=='history':source.pop('actual_upstream_original_pre_power_source')
        elif fixture=='full_pressure':source['actual_pre_pressure_function_used_in_tensor']=False
        elif fixture=='mode':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['axial']['linear_axial_moment']['extra_source']='incoming1'
        else:source['physical_cylindrical_stress_mixed3']['theta']['radial_shear']['r0_z0']['signed_coefficient']=field.ctx.mpf(13)
        try:validate_whole_O3_power_view(field.entrancecone,bad)
        except ValueError:pass
        else:raise ValueError('Invalid O3 source admitted: '+fixture)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,exact_original_O2_O3_and_axial_stress_identities=len(field.theorem['identities']),
        exact_original_full_theta_identities=len(field.theorem['current_full_theta_theorem']['identities']),
        exact_whole_Z_canonical_source_identities=len(field.bindings['exact_whole_Z_canonical_source_identities']),
        exact_actual_composed_phase0_correlation_and_pressure_endpoint_identities=len(field.bindings['actual_composed_source_phase0_M_K_X_correlations']),
        consumed_current_incoming_production_function_identities=len(field.bindings['current_complete_production_function_proof']['identities']),
        strictly_positive_directed_whole_domain_bounds=len(replay['positive_margins']),
        actual_O2_shared_mass_and_positive_transition_integral_correlated=True,
        all_ten_original_signed_O3_sectors_and_full_nonzero_energy_pressure_retained=True,
        original_absolute_pressure_backward_ODE_and_same_phase1_inlet_consumed=True,
        canonical_axial_stress_zero_sectors_proved_from_actual_operator_not_silenced=True,
        energy_direction_bound_uses_correlated_time_and_Z_AMGM=True,
        both_current_transition_power_and_power_entrance_function_joins_consumed=True,
        foreign_partial_missing_history_changed_mode_and_changed_physical_row_rejected=True,
        whole_domain_proof_not_phase_samples_or_interval_overlap=True,historical_cone_admission_not_promoted=True,
        current_strict_nonzero_regions_including_inherited=15,current_exact_zero_exterior_regions=1,
        remaining_registry_regions_without_current_cone=17,scope=raw['scope'],
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole O3 power cone PASS; variable O3 transition/O2/inner and actual recursion/waves open',flush=True)
    return result


if __name__=='__main__':run()
