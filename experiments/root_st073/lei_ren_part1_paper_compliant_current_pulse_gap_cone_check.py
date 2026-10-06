"""Whole current gap source bounds and signed-factor/shear rejection checks."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_pulse_gap_cone import (
    CurrentPulseGapCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,pack,encode,
    endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_pulse_gap_cone_operator import validate_whole_gap_views,whole_current_gap_bounds
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPulseGapCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current whole gap manifest/scope differs')
    if json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))!=encode(pack(field.whole_views)):
        raise ValueError('Complete signed gap views differ')
    replay=whole_current_gap_bounds(field.maincone,field.whole_views)
    if encode(pack(replay))!=raw['current_whole_gap_correlated_bounds']:raise ValueError('Current gap bound replay differs')
    count=0
    for group in ('positive_parameter_margins','positive_log_margins','positive_monotonicity_margins','positive_algebraic_margins'):
        for value in replay[group].values():
            if endpoints(value)[0]<=0:raise ValueError('Inconclusive current whole gap bound')
            count+=1
    if len(replay['normalized_signed_stress_sector_bounds'])!=8:raise ValueError('Eight signed gap error sectors required')
    for view in field.whole_views.values():
        if not view['original_complete_view']['full_nonzero_meridional_velocity_and_radial_remainder_retained']:
            raise ValueError('Nonzero radial history/remainder discarded')
    for fixture in ('foreign','D2','axial','missing'):
        bad=copy_containers(field.whole_views);view=bad['pulse_gap'];native=view['original_complete_view']
        sectors=native['current_actual_source_stress_packet']['full_meridional_stress_log_sectors']
        if fixture=='foreign':view['implicit_source_sha256']='foreign'
        elif fixture=='D2':sectors['theta']['mixed_angular_axial_history']['selected_D_recipe']='D1'
        elif fixture=='axial':
            rows=native['current_source_three_component_velocity_rows']['axial']
            rows[0]=[field.ctx.mpf(1)]+[rows[0][n] for n in range(1,6)]
        else:sectors['axial'].pop('same_absolute_pressure_memory')
        try:validate_whole_gap_views(field.maincone,bad)
        except ValueError:pass
        else:raise ValueError('Incorrect current gap source admitted: '+fixture)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_signed_gap_algebra_identities=len(field.theorem['identities']),
        current_grouped_log_and_exact_endpoint_identities=len(field.bindings['current_grouped_log_and_endpoint_identities']),
        strictly_positive_directed_whole_gap_bounds=count,both_whole_current_signed_gap_views_retained=True,
        all_nine_signed_sectors_per_region_retained=True,nonzero_radial_velocity_and_remainder_not_discarded=True,
        foreign_source_wrong_D2_nonzero_axial_shear_and_missing_pressure_sector_rejected=True,
        exact_correlated_endpoints_not_phase_derivatives_or_interval_log_subtraction=True,
        historical_admissions_not_promoted_to_current_graph=True,actual_global_wave_and_completed_tensor_admission=False,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole gap/gap-end strict signed cone PASS; global/wave/recursion remain open',flush=True)
    return result


if __name__=='__main__':run()
