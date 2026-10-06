"""Focused whole-domain current angular/tail source, cone and scope check."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_angular_tail_cone import (
    CurrentAngularTailCone,NAME,RECEIPT,VIEWS_NAME,GATES,OPEN,HERE,PREFIX,sha,
    pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_angular_tail_cone_operator import (
    validate_whole_angular_tail_views,whole_current_angular_tail_bounds)
from lei_ren_part1_paper_compliant_current_pulse_main_exit_cone_check import copy_containers


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentAngularTailCone(require_checked=False)
    field.assert_graph();expected=encode(pack(field.manifest()));expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if expected!=raw or any(raw[k] for k in GATES+OPEN):raise ValueError('Current cone manifest or regional scope differs')
    if json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes()))!=encode(pack(field.whole_views)):
        raise ValueError('Current complete angular/tail signed views differ')
    replay=whole_current_angular_tail_bounds(field.flattencone,field.whole_views)
    if encode(pack(replay))!=raw['current_whole_angular_tail_correlated_bounds']:
        raise ValueError('Current whole continuous bound replay differs')
    for value in replay['positive_margins'].values():
        if endpoints(value)[0]<=0:raise ValueError('Current cone bound inconclusive')
    for fixture in ('foreign','domain','pressure','complete','zero'):
        bad=copy_containers(field.whole_views);v=bad['outer_angular'];native=v['original_complete_view']
        if fixture=='foreign':v['implicit_source_sha256']='foreign'
        elif fixture=='domain':native['angular_offset']=field.ctx.mpf([-4,-2])
        elif fixture=='pressure':native['current_actual_normalized_full_moment_rows'].pop('P')
        elif fixture=='complete':native['actual_full_stress_not_local_difference']=False
        else:bad['heat_exterior']['original_complete_view']['source_exact_Gamma_tensor_and_remainder_zero']=False
        try:validate_whole_angular_tail_views(field.flattencone,bad)
        except ValueError:pass
        else:raise ValueError('Incorrect current cone source admitted: '+fixture)
    count=sum(len(v['identities']) if 'identities' in v else sum(k.endswith('_verified') and val is True for k,val in v.items())
        for k,v in field.theorem.items() if k not in ('input_hashes','historical_cone_receipts_loaded_or_promoted'))
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_generic_algebra_identities=count,
        current_actual_AST_source_identities=len(field.bindings['identities']),
        strictly_positive_directed_whole_bounds=len(replay['positive_margins']),
        original_continuous_partial_integral_FTC_positivity_used=True,
        same_current_waiting_heat_endpoint_full_moments_and_actual_O7_operator_bound=True,
        current_complete_pressure_energy_X_and_K_Z_preserved=True,
        current_seven_adjacent_and_four_internal_tensor_joins_preserved=True,
        strict_heat_cone_excludes_zero_endpoint_and_unbounded_zero_exterior=True,
        uniform_terminal_direction_finite_exact_source_constant_not_materialized=True,
        foreign_partial_domain_missing_pressure_incomplete_and_unproved_zero_rejected=True,
        exact_S_not_cap_endpoint_and_KR_not_counted_twice=True,
        historical_cone_admissions_not_promoted=True,
        current_strict_nonzero_regions_including_inherited=13,current_exact_zero_exterior_regions=1,
        remaining_registry_regions_without_current_cone=19,
        scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        **dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current angular/tail cone PASS; global waves/recursion open',flush=True)
    return result


if __name__=='__main__':run()
