"""Accept actual modified native dispatch/continuation, not global cones/NS."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_modified_tensor_dispatch import (
    CurrentModifiedTensorDispatch,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,LOCAL_GATES,GATES,OPEN,
    QUERIES,query,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_tensor_dispatch_operator import exact_modified_dispatch_theorem
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_check import compare_exit_physical_groups

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedTensorDispatch(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[key] for key in LOCAL_GATES+GATES+OPEN):
        raise ValueError('Actual modified native dispatch manifest/scope changed')
    if encode(pack(exact_modified_dispatch_theorem(field)))!=raw['exact_actual_modified_dispatch_source_theorem']:
        raise ValueError('Same actual repaired source/continuation/operator lineage changed')
    views={key:query(field,kind,region,value) for key,kind,region,value in QUERIES}
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Complete changed/native power continuation views differ')
    if len(field.registry.routes)!=33 or raw['actual_modified_33_chart_native_dispatch']!=list(field.registry.routes):
        raise ValueError('Actual modified native route table incomplete')
    if raw['affected_native_routes']!=['O2_buffer','O3_slope_mu','O3_power']:
        raise ValueError('Modified affected native source chart inventory differs')
    for view in views.values():
        if any(view[key] for key in OPEN) or not view['pieces_are_alternative_source_charts_not_summed_fields']:
            raise ValueError('Directed dispatch must preserve alternative source charts and open cone/NS scope')
        for piece in view['source_pieces']:
            inner=piece['complete_tensor_source']['original_complete_view']
            if endpoints(inner['exact_physical_divergence'])!=(0,0) or endpoints(inner['exact_completed_tensor_radial_divergence'])!=(0,0):
                raise ValueError('Original/modified physical divergence completion missing')
    # A whole normalized phase[0,1] query must include true phase1, while
    # its local part stays inside the parent's accepted domain.
    whole=views['whole_modified_native_power']['source_pieces']
    if len(whole)!=2 or endpoints(whole[0]['source_coordinate_box'])[0]!=0 or endpoints(whole[1]['source_coordinate_box'])[1]!=1:
        raise ValueError('Full native power phase coverage missing')
    if endpoints(whole[0]['source_coordinate_box'])[1]!=endpoints(whole[1]['source_coordinate_box'])[0]:
        raise ValueError('Power phase source partition has a gap')
    if endpoints(whole[0]['same_actual_q_offset_enclosure'])[1]>=2:
        raise ValueError('Native local query escaped the checked parent domain')
    phase_end=views['exact_original_power_endpoint']['source_pieces']
    if len(phase_end)!=1 or endpoints(phase_end[0]['source_coordinate_box'])!=(1,1):
        raise ValueError('Actual correlated power endpoint lost')
    endpoint=phase_end[0]['complete_tensor_source']['original_complete_view']
    if not endpoint['actual_q_equals_same_original_Tw_times_phase']:
        raise ValueError('Tw enclosure cap was used in place of the actual phase1 source')
    # Keep the nonzero last compact repair part in a box crossing q=2.
    seam=views['straddling_q2_with_nonzero_repair']['source_pieces']
    if len(seam)!=2 or endpoints(seam[0]['source_coordinate_box'])[1]!=2 or endpoints(seam[1]['source_coordinate_box'])[0]!=2:
        raise ValueError('Straddling q2 box must use two source branches')
    local_source=seam[0]['complete_tensor_source']['original_complete_view']['current_actual_modified_signed_source']['actual_source']
    if local_source['exact_functional_exit_reduction_applied']:
        raise ValueError('A straddling repair box was incorrectly reset to zero')
    if not any(endpoints(row)!=(0,0) for row in local_source['current_independent_bump_logR4_rows'][2]):
        raise ValueError('Actual final support source lost in the straddling box')
    # Same complete physical endpoint source, independently assembled
    # from original retained raw histories/pressure/velocity programs.
    actual=views['exact_q2']['source_pieces'][0]['complete_tensor_source']['original_complete_view']
    reference=field.continuation((-1,1),2,'log_radius_offset',('-3','-1'),None,'1')
    compare_exit_physical_groups(actual,reference)
    # Fresh continuation source compared with the accepted original native
    # incoming owner. Exact source identity is established before bounds;
    # differently factored interval tensor rows need not be identical.
    fresh=field.continuation('.537','2.413337','log_radius_offset','-2.337','.337','.8')
    original=field.registry.native('O3_power','.537',field.ctx.mpf('2.413337')/field.Tw,'-2.337','.337','.8')
    if encode(pack(fresh['actual_original_continuation_signed_source']['actual_source']['current_original_pre_source']))!=encode(pack(original['original_complete_view']['actual_upstream_original_pre_power_source'])):
        raise ValueError('Fresh continuation source differs from the actual accepted incoming provider')
    # One unchanged chart exercises the complementary30-route delegation.
    unchanged=field.native('pulse_entrance','.537','.013337','-2.337','.337','.8')
    direct=field.registry.native('pulse_entrance','.537','.013337','-2.337','.337','.8')
    if encode(pack(unchanged['source_pieces'][0]['complete_tensor_source']))!=encode(pack(direct)):
        raise ValueError('Unchanged native route did not preserve the original complete source')
    try:field.power_offset(0,field.ctx.mpf(endpoints(field.Tw)[1])+1)
    except ValueError:pass
    else:raise ValueError('Independent q beyond actual source bounds was admitted')
    result=dict(actual_five_defect_family_sha256=field.registry.family,implicit_source_sha256=field.registry.source,
      datum_enclosure_sha256=field.registry.datum_sha,modified_dispatch_definition_sha256=field.modified_dispatch_definition_sha256,
      finite_integer_N=field.N,new_exact_modified_power_transport_dispatch_source_identities=len(field.theorem['identities']),
      consumed_same_actual_repaired_history_identities=len(field.theorem['consumed_same_implicit_repaired_history_source_theorem']['identities']),
      modified_native_route_count=33,actually_changed_native_route_count=3,unchanged_same_source_native_route_count=30,
      actual_original_power_full_phase_domain_and_true_Tw_endpoint_covered=True,
      q2_straddling_box_retains_actual_nonzero_repair_source=True,
      q2_complete_physical_exit_source_equals_independent_original_source=True,
      post_support_source_closure_uses_same_implicit_root_no_pressure_reset=True,
      ordinary_logR_derivatives_not_phase_rescaled=True,
      own_energy_modified_interfaces_common_cone_N_global_admissibility_and_full_NS_remain_open=True,
      current_strict_original_nonzero_whole_regions=15,remaining_original_whole_regions_without_cone=17,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(LOCAL_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual modified33 tensor/source dispatcher and original full power continuation PASS; interfaces/common N/cones/full NS open',flush=True)
    return result

if __name__=='__main__':run()
