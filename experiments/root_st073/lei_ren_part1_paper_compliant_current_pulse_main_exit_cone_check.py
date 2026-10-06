"""Current-domain source replay, strict margins and wrong-source rejection."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone import (
    CurrentPulseMainExitCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,
    sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_operator import (
    validate_whole_views,whole_current_main_exit_bounds)


def copy_containers(value):
    # Interval/Taylor source objects carry contexts that cannot be deepcopied.
    # Copy only the fixture containers; immutable source intervals stay shared.
    if isinstance(value,dict):return {key:copy_containers(v) for key,v in value.items()}
    if isinstance(value,list):return [copy_containers(v) for v in value]
    if isinstance(value,tuple):return tuple(copy_containers(v) for v in value)
    return value


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseMainExitCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()))
    expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current whole main/exit manifest/scope differs')
    stored=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))
    if stored!=encode(pack(field.whole_views)):raise ValueError('Complete current signed whole views differ')
    replay=whole_current_main_exit_bounds(field.cone,field.whole_views)
    if encode(pack(replay))!=raw['current_whole_main_exit_correlated_bounds']:
        raise ValueError('Current exact whole-source bounds replay differs')
    positive_count=0
    for group in ('positive_parameter_margins','positive_log_margins',
        'positive_monotonicity_margins','positive_shape_and_error_margins','positive_algebraic_margins'):
        for name,value in replay[group].items():
            if endpoints(value)[0]<=0:raise ValueError('Unresolved whole current cone margin: '+name)
            positive_count+=1
    if len(replay['normalized_stress_sector_bounds'])!=12:raise ValueError('All twelve current signed corrections required')
    bad=copy_containers(field.whole_views)
    rawmain=bad['pulse_main']['original_complete_signed_tensor_view']['original_complete_view']
    rawmain['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']['theta'].pop('radial_shear')
    try:validate_whole_views(field.cone,bad)
    except ValueError:pass
    else:raise ValueError('Omitted signed correction admitted')
    bad=copy_containers(field.whole_views)
    bad['pulse_exit']['original_complete_signed_tensor_view']['implicit_source_sha256']='foreign source'
    try:validate_whole_views(field.cone,bad)
    except ValueError:pass
    else:raise ValueError('Foreign source admitted')
    bad=copy_containers(field.whole_views)
    for region in bad:
        actual=bad[region]['original_complete_signed_tensor_view']['original_complete_view']
        actual['current_selected_ap']=[field.ctx.mpf('1.3')]+[actual['current_selected_ap'][i] for i in range(1,6)]
    try:whole_current_main_exit_bounds(field.cone,bad)
    except ArithmeticError:pass
    else:raise ValueError('Amplitude outside analytic theorem admitted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_signed_kernel_shear_algebra_identities=len(field.theorem['identities']),
        current_exact_grouped_source_log_identities=len(field.bindings['current_grouped_log_identities']),
        strictly_positive_directed_whole_domain_bounds=positive_count,
        both_whole_current_signed_view_records_retained=True,all_fifteen_signed_sectors_per_region_retained=True,
        all_twelve_corrections_bounded_only_as_errors_to_signed_baseline=True,
        original_full_kernel_two_IBP_remainder_and_axial_amplitude_derivative_retained=True,
        missing_signed_sector_foreign_source_and_oversized_amplitude_rejected=True,
        whole_domain_proof_not_phase_samples_or_direct_interval_overlap=True,
        historical_source_admission_not_promoted_to_current_field=True,
        actual_wave_uniform_edge_and_completed_tensor_admissions=False,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),
            Path(__file__).name:sha(Path(__file__).name)},**dict.fromkeys(GATES,True),
        **dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole current main/exit signed cone PASS; actual covariance waves and recursion remain open',flush=True)
    return result


if __name__=='__main__':run()
