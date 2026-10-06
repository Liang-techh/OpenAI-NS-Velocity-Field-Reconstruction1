"""Admit exact analytic pressure/heat inheritance, keep global NS open."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_modified_heat_inheritance import (
    CurrentModifiedHeatInheritance,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,ANCESTOR_GATES,GATES,OPEN,
    DOWNSTREAM,DOWNSTREAM_SEAMS,QUERIES,TRACE_QUERIES,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_heat_inheritance_operator import exact_modified_heat_inheritance_theorem
from lei_ren_part1_paper_compliant_current_heat_background_tensor_check import check_view,check_physical,check_interface

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedHeatInheritance(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in ANCESTOR_GATES+GATES+OPEN):
        raise ValueError('Inherited heat source manifest/source/scope changed')
    if encode(pack(exact_modified_heat_inheritance_theorem(field)))!=raw['exact_actual_modified_downstream_heat_inheritance_theorem']:
        raise ValueError('Actual pressure/compact repair/heat source-function composition changed')
    views={key:field.downstream(region,(-1,1),value,('-3','-1'),None,'1') for key,region,value in QUERIES}
    views['pressure_datum']=field.pressure_datum((-1,1))
    for name in TRACE_QUERIES:views['trace_'+name]=field.trace(name)
    views['unbounded_Gamma_exterior']=field.unbounded_exterior()
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Actual inherited complete downstream/heat source views differ')
    if raw['inherited_downstream_regions']!=list(DOWNSTREAM) or raw['inherited_downstream_interfaces']!=list(DOWNSTREAM_SEAMS):
        raise ValueError('Original fifteen-region/source-trace inheritance inventory changed')
    for key,region,value in QUERIES:
        view=views[key]
        if any(view[k] for k in OPEN) or not view['actual_implicit_source_closure_and_whole_Z_datum_function_link_consumed']:
            raise ValueError('Regional inheritance changed global admission scope')
        piece=view['complete_same_source_view']['source_pieces']
        if len(piece)!=1 or piece[0]['source_branch']!='unchanged_same_checked_original_route':
            raise ValueError('An original downstream source was replaced')
    # Function identities above prove congruence on full domains. Fresh
    # queries below ensure runtime delegation keeps actual source packets.
    fresh=field.downstream('heat_exterior','.537','4.173337','-2.337','.337','.8')
    direct=field.registry.native('heat_exterior','.537','4.173337','-2.337','.337','.8')
    inherited=fresh['complete_same_source_view']['source_pieces'][0]['complete_tensor_source']
    if encode(pack(inherited))!=encode(pack(direct)):raise ValueError('Fresh heat source differed from accepted original provider')
    if not inherited['original_complete_view']['actual_Gamma_regional_physical_NS_identity']:
        raise ValueError('Actual regional heat identity was not carried through the consumer')
    for key in ('collar','exterior_far'):
        inner=views[key]['complete_same_source_view']['source_pieces'][0]['complete_tensor_source']['original_complete_view']
        check_view(inner)
    for name in ('waiting_collar','collar_exterior'):
        check_interface(views['trace_'+name]['complete_same_source_view']['original_complete_trace'])
    exterior=views['unbounded_Gamma_exterior']['complete_same_source_view']['original_complete_view']
    if not exterior['original_unbounded_exterior_covered'] or not exterior['finite_anchor_only_supplies_tensor_component_layout']:
        raise ValueError('An exterior cutoff replaced the full unbounded Gamma source')
    if not exterior['actual_stress_tensor_divergence_remainder_and_momentum_exactly_zero']:
        raise ValueError('Full original Gamma source identity missing')
    rows,_=check_physical(exterior,exact_zero=True)
    datum=views['pressure_datum']['complete_same_source_view']
    if encode(pack(datum))!=encode(pack(field.heat_datum.normalized_jets((-1,1),5))):
        raise ValueError('Original analytic P0 datum differs across the restored source graph')
    if len(datum['normalized_pressure_coefficients'])!=6:
        raise ValueError('Absolute analytic P0 axial5 coefficient source lost')
    try:field.downstream('O2_buffer',0,0)
    except ValueError:pass
    else:raise ValueError('Inherited downstream API admitted a changed upstream source without its own wrapper')
    result=dict(actual_five_defect_family_sha256=field.registry.family,implicit_source_sha256=field.registry.source,
      datum_enclosure_sha256=field.registry.datum_sha,finite_integer_N=field.N,
      modified_heat_inheritance_definition_sha256=field.modified_heat_inheritance_definition_sha256,
      new_exact_modified_downstream_heat_source_identities=len(field.theorem['identities']),
      consumed_actual_compact_repair_power_transport_identities=len(field.theorem['consumed_actual_compact_repair_original_power_transport']['identities']),
      consumed_actual_modified_interface_identities=len(field.interfaces.theorem['identities']),
      inherited_downstream_native_region_count=15,inherited_downstream_function_trace_count=15,
      inherited_unbounded_exterior_exact_zero_physical_rows=rows,
      original_analytic_P0_defining_function_and_axial5_projection_preserved=True,
      same_native_raw_waiting_root_and_complete_forward_Gamma_pressure_FTC_retained=True,
      full_exterior_identity_is_regional_and_not_radial_truncation_or_global_NS=True,
      own_modified_total_kinetic_energy_common_N_cones_higher_interfaces_and_recursion_remain_open=True,
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(ANCESTOR_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Modified downstream/analytic preheat/full unbounded Gamma inheritance PASS; higher interfaces/cones/energy/recursion/global NS open',flush=True)
    return result

if __name__=='__main__':run()
