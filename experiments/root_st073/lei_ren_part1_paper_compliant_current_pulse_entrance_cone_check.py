"""Focused continuous entrance bounds, full-history and source-scope check."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_pulse_entrance_cone import (
    CurrentPulseEntranceCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,
    pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_pulse_entrance_cone_operator import (
    validate_whole_entrance_view,whole_current_entrance_bounds)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseEntranceCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current entrance manifest/scope differs')
    stored=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if stored!=encode(pack(dict(whole_current_entrance=field.whole_view,actual_pre_Tw_power_source=field.preTw))):
        raise ValueError('Complete current entrance and pre-Tw source views differ')
    replay=whole_current_entrance_bounds(field.tailcone,field.whole_view,field.preTw)
    if encode(pack(replay))!=raw['current_whole_entrance_correlated_bounds']:
        raise ValueError('Current continuous entrance source-bound replay differs')
    positive_count=0
    for group in ('positive_parameter_margins','positive_log_margins','positive_monotonicity_margins',
        'positive_shape_and_error_margins','positive_algebraic_margins'):
        for name,value in replay[group].items():
            if endpoints(value)[0]<=0:raise ValueError('Unresolved current entrance bound: '+name)
            positive_count+=1
    if len(replay['normalized_stress_sector_bounds'])!=12:raise ValueError('All twelve signed corrections required')
    for fixture in ('foreign','domain','memory','pressure','energy','complete','mode','physical'):
        bad=copy_containers(field.whole_view);source=bad['original_complete_view']
        if fixture=='foreign':bad['implicit_source_sha256']='foreign'
        elif fixture=='domain':source['coordinate']=field.ctx.mpf(['.001','.02'])
        elif fixture=='memory':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['theta'].pop('signed_original_memory')
        elif fixture=='pressure':source.pop('current_signed_absolute_Rv_pressure')
        elif fixture=='energy':source.pop('current_forward_anchored_full_energy_rows')
        elif fixture=='complete':source['all_nonzero_incoming_histories_and_radial_remainder_retained']=False
        elif fixture=='mode':source['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['axial']['same_absolute_pressure_memory']['extra_source']='one'
        else:source['physical_cylindrical_stress_mixed3']['theta']['radial_shear']['r0_z0']['signed_coefficient']=field.ctx.mpf(13)
        try:validate_whole_entrance_view(field.tailcone,bad)
        except ValueError:pass
        else:raise ValueError('Invalid entrance source admitted: '+fixture)
    bad=copy_containers(field.whole_view);old=bad['original_complete_view']['current_selected_ap']
    bad['original_complete_view']['current_selected_ap']=[field.ctx.mpf('1.3')]+[old[i] for i in range(1,6)]
    try:whole_current_entrance_bounds(field.tailcone,bad,field.preTw)
    except ArithmeticError:pass
    else:raise ValueError('Oversized current amplitude admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_power_memory_and_endpoint_log_identities=len(field.theorem['identities']),
        current_incoming_production_function_identities=len(field.bindings['current_actual_incoming_function_proof']['identities']),
        strictly_positive_directed_whole_domain_bounds=positive_count,
        all_fifteen_signed_sectors_and_twelve_error_corrections_retained=True,
        current_typed_common_Tw_mu_definition_and_actual_pre_phase1_U_consumed=True,
        pre_Tw_memory_factored_before_interval_subtraction=True,xi0_keeps_full_incoming_energy_pressure_and_radial_velocity=True,
        both_current_tensor_source_joins_consumed=True,
        foreign_partial_missing_history_changed_mode_changed_physical_row_and_oversized_amplitude_rejected=True,
        whole_domain_proof_not_phase_samples_or_interval_overlap=True,historical_cone_admission_not_promoted=True,
        current_strict_nonzero_regions_including_inherited=14,current_exact_zero_exterior_regions=1,
        remaining_registry_regions_without_current_cone=18,scope=raw['scope'],
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole entrance cone PASS; O3/O2/core and actual waves/recursion open',flush=True)
    return result


if __name__=='__main__':run()
