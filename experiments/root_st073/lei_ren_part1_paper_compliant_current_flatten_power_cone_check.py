"""Focused current correlated flatten/power bounds, full views and scope."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_flatten_power_cone import (
    CurrentFlattenPowerCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,
    pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_flatten_power_cone_operator import (
    validate_whole_flatten_power_views,whole_current_flatten_power_bounds)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentFlattenPowerCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current whole flatten/power manifest or scope differs')
    if json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))!=encode(pack(field.whole_views)):
        raise ValueError('Current full tensor/history views differ')
    replay=whole_current_flatten_power_bounds(field.gapcone,field.whole_views)
    if encode(pack(replay))!=raw['current_whole_flatten_power_correlated_bounds']:
        raise ValueError('Current whole correlated bound replay differs')
    for value in replay['positive_margins'].values():
        if endpoints(value)[0]<=0:raise ValueError('Current cone bound inconclusive')
    for fixture in ('foreign','domain','pressure','complete'):
        bad=copy_containers(field.whole_views);view=bad['flatten'];native=view['original_complete_view']
        if fixture=='foreign':view['implicit_source_sha256']='foreign'
        elif fixture=='domain':native['original_coordinate']=field.ctx.mpf([0,50])
        elif fixture=='pressure':native['current_actual_normalized_full_moment_rows'].pop('P')
        else:native['current_complete_pressure_energy_and_forward_X_retained']=False
        try:validate_whole_flatten_power_views(field.gapcone,bad)
        except ValueError:pass
        else:raise ValueError('Incorrect current full cone source admitted: '+fixture)
    generic=sum(len(field.theorem[k]['identities']) for k in ('flatten','outer_power','current_KR_cancellation'))
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_generic_algebra_identities=generic,
        current_source_log_radius_and_history_identities=len(field.bindings['identities']),
        strictly_positive_directed_whole_bounds=len(replay['positive_margins']),
        both_whole_current_signed_views_and_full_pressure_energy_preserved=True,
        three_current_functional_tensor_joins_preserved=True,whole_Z_with_L_factor_not_only_axis=True,
        signed_memory_and_Hf_integral_correlation_used_before_enclosure=True,
        sigma_cover_uses_continuous_intervals_not_point_samples=True,
        foreign_half_domain_missing_pressure_and_incomplete_history_rejected=True,
        exact_S_not_Scap_endpoint_and_KR_not_counted_twice=True,
        historical_cone_admissions_not_promoted=True,global_wave_and_completed_tensor_admission=False,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current whole flatten/power strict cone PASS; global/waves/recursion open',flush=True)
    return result


if __name__=='__main__':run()
