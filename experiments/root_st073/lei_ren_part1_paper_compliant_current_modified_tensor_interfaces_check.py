"""Accept affected source-function tensor traces at their retained orders."""
import gzip
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces import (
    CurrentModifiedTensorInterfaces,HERE,PREFIX,NAME,RECEIPT,VIEWS_NAME,LOCAL_GATES,DISPATCH_GATES,
    GATES,OPEN,SEAMS,REPAIR_EDGES,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces_operator import exact_modified_interfaces_theorem
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_check import compare_exit_physical_groups

@source_precision
def run(field=None):
    field=field if field is not None else CurrentModifiedTensorInterfaces(require_checked=False)
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=field.manifest();expected['input_hashes'][VIEWS_NAME]=sha(VIEWS_NAME)
    if encode(pack(expected))!=raw or any(raw[k] for k in LOCAL_GATES+DISPATCH_GATES+GATES+OPEN):
        raise ValueError('Actual affected interface manifest or acceptance scope changed')
    if encode(pack(exact_modified_interfaces_theorem(field)))!=raw['exact_actual_affected_interface_source_theorem']:
        raise ValueError('Actual callable profile/history/endpoint/operator bindings changed')
    views={name:field.interface(name) for name in SEAMS}
    if encode(pack(views))!=json.loads(gzip.decompress((HERE/VIEWS_NAME).read_bytes())):
        raise ValueError('Actual complete left/right tensor traces differ from stored views')
    if len(views)!=12 or raw['affected_named_interfaces']!=list(SEAMS):raise ValueError('Affected interface inventory incomplete')
    for name,view in views.items():
        if view['current_common_tensor_contribution_count']!=71 or len(view['common_actual_tensor_rows'])!=71:
            raise ValueError('Full signed completed tensor/divergence/remainder trace omitted')
        if any(view[k] for k in OPEN) or not view['interval_overlap_not_used_as_source_identity']:
            raise ValueError('Retained-order source identity cannot promote global/higher-order/cone admission')
        for side in ('left','right'):
            inner=view[side+'_complete_tensor_source']
            if endpoints(inner['exact_physical_divergence'])!=(0,0) or endpoints(inner['exact_completed_tensor_radial_divergence'])!=(0,0):
                raise ValueError('Own moment incompressibility and physical tensor completion missing')
    # Shared t=0 means the actual same cumulative integral, not endpoints
    # selected independently or a zero replacement at the old seam.
    source=lambda view:view['current_actual_modified_signed_source']['actual_source']
    join=views['buffer_transition']
    left=source(join['left_complete_tensor_source']);right=source(join['right_complete_tensor_source'])
    if encode(pack(left['signed_actual_cumulative_scalar_enclosures']))!=encode(pack(right['signed_actual_cumulative_scalar_enclosures'])):
        raise ValueError('Buffer/O3 seam did not retain the same actual cumulative integral')
    # The support-end fields have zero local profiles but live own history.
    modulation=views['transition_modulation_end']
    a=modulation['left_complete_tensor_source']['actual_named_endpoint_source']
    b=modulation['right_complete_tensor_source']['actual_named_endpoint_source']
    if encode(pack(a['signed_actual_cumulative_scalar_enclosures']))!=encode(pack(b['signed_actual_cumulative_scalar_enclosures'])):
        raise ValueError('Flat modulation endpoint reset the cumulative source')
    if endpoints(a['signed_actual_cumulative_scalar_enclosures']['e_kinetic'])[0]<=0:
        raise ValueError('Positive modulation kinetic history was discarded')
    # Every named repair edge must preserve source-correlated zero/full
    # integrals. No partial support is rounded into another endpoint.
    for name in REPAIR_EDGES:
        a=field.endpoint_source(name,(-1,1),'left');b=field.endpoint_source(name,(-1,1),'right')
        for row in a['current_independent_bump_logR4_rows']:
            if any(endpoints(v)!=(0,0) for v in row):raise ValueError('Named flat repair profile has a nonzero jet')
        for key in ('current_same_source_partial_bump_weights','current_same_source_partial_repair_primitive_changes',
          'modified_five_histories_in_original_normalized_units_Pstar_sectors',
          'modified_absolute_pressure_over_Pstar2_ordinary_logR_rows'):
            if encode(pack(a[key]))!=encode(pack(b[key])):raise ValueError('Named repair one-sided FTC source mismatch: '+key)
        if a['exact_functional_exit_reduction_applied']:raise ValueError('Named repair edge used the q2 interval reset branch')
    compare_exit_physical_groups(views['repair_exit']['left_complete_tensor_source'],views['repair_exit']['right_complete_tensor_source'])
    # Fresh axial/time/angle/viscosity arguments exercise the named API.
    fresh=field.interface('repair1_out','.537','-2.337','.337','.8')
    if len(fresh['common_actual_tensor_rows'])!=71:raise ValueError('Fresh named interface missing components')
    try:field.endpoint_source('repair1_out','.537','interval')
    except ValueError:pass
    else:raise ValueError('Named source germ accepted a non-source side')
    result=dict(actual_five_defect_family_sha256=field.dispatch.registry.family,
      implicit_source_sha256=field.dispatch.registry.source,datum_enclosure_sha256=field.dispatch.registry.datum_sha,
      finite_integer_N=field.N,modified_interfaces_definition_sha256=field.modified_interfaces_definition_sha256,
      new_exact_affected_source_interface_identities=len(field.theorem['identities']),
      consumed_actual_repaired_history_identities=len(field.theorem['consumed_actual_repaired_history_theorem']['identities']),
      affected_named_interface_count=12,complete_common_tensor_groups_per_interface=71,
      named_edge_correlations_and_zero_full_integrals_are_source_level_reductions=True,
      flat_local_profiles_do_not_zero_own_moments_pressure_or_kinetic_history=True,
      original_native_inventory_and_unaffected_sources_preserved=True,
      source_orders=views['repair_exit']['source_orders'],
      scope=raw['scope'],input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
      **dict.fromkeys(LOCAL_GATES+DISPATCH_GATES+GATES,True),**dict.fromkeys(OPEN,False),all_passed=True)
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Twelve affected full tensor source traces PASS; higher orders/heat/cones/energy/recursion/NS open',flush=True)
    return result

if __name__=='__main__':run()
