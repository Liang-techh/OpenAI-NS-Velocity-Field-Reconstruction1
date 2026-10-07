"""Check original graph dispatch, explicit provider frontier and context bridge."""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_signed_u_factored_oracle as current
import lei_ren_part1_paper_compliant_current_native_Rc_factored_source_oracle_check as inherited

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
cover=current.cover;packets=cover.packets


def require(ok,message):
    if not ok:raise AssertionError(message)


@cover.native.inlet.source_precision
def run(role_owner,report=None):
    began=time.monotonic();built=role_owner.build();oracle=current.SignedUCoverFactoredOracle(role_owner,built)
    report=json.loads((HERE/current.NAME).read_bytes()) if report is None else report
    require(report[current.GATE] and report['source_family']==oracle.source_family,'Same original all17 producer required')
    for name,digest in report['input_hashes'].items():require(sha(name)==digest,'Changed original range prerequisite: '+name)
    records=report['original_chart_query_records']
    require(set(records)==set(current.POINTS) and len(records)==17,'All original providers must be attempted')
    enclosed=[name for name,row in records.items() if row['status']=='enclosed']
    unresolved=[name for name,row in records.items() if row['status']!='enclosed']
    require(len(enclosed)==report['enclosed_original_chart_queries']==14,'Declared original source count differs')
    require(set(unresolved)=={'bridge_first','actual_patch','O2_axial'},'Explicit unresolved source frontier differs')
    require(report['actual_original_role_dispatches']==142,'Declared original role dispatch count differs')
    require(not report['all17_whole_chart_or_continuous_route_ranges_admitted'] and not any(report.get(key) for key in packets.OPEN),
        'Provider queries cannot complete whole-chart/global gates')
    reasons={}
    for name in unresolved:
        source=records[name]['actual_original_spatial_source']
        branches=[row for cell in source['spatial_signed_density_branch_cells'] for row in cell['original_conditional_branch_first_jets']]
        require(any(row['status']!='enclosed' for row in branches),'Unresolved provider must retain its failed branch')
        reasons[name]=dict(original_q_status=source['original_q_slow_jet_source']['status'],
            branches=[dict(branch=row['signed_u_branch'],status=row['status'],obstruction=row.get('arithmetic_obstruction')) for row in branches])
    c=oracle.ctx;box=c.mpf((cover.ep(c.mpf('.12'))[0],cover.ep(c.mpf('.15'))[1]))
    broad=oracle.density_frame(chart='O2_slope',Z=(-1,1),coordinate=box,N=2048)
    require(broad.values is not None,'Actual broad signed-crossing O2 cannot regress')
    rows={(row.get('chart'),row.get('function_role')):row for row in built['graph'].nodes if row['operation']=='original_function_graph'}
    dispatched=0
    for key in cover.density.RATES:
        for order in ('C0','Z'):
            value=oracle.dispatch_function_range(rows['O2_slope','density_'+key+'_'+order],broad)
            require(value is broad.values[order][key] and value.ctx is c,'Dispatch must return the same original factored range')
            dispatched+=1
    quiet=oracle.density_frame(chart='bridge_first',Z=(-1,1),coordinate={'selected_sc_multiple':'1/2'},N=2048)
    # The existing analytic collar transfer has an exact-flat theorem, but
    # this generic q-jet range provider still cannot detect its cutoff branch.
    # Preserve that obstruction rather than inventing a zero source frame.
    require(quiet.values is None and quiet.record['status']=='requires_source_or_phase_refinement',
        'Unresolved native collar provider must remain explicit')
    require('separate_original_P0' in quiet.record and 'separate_original_P0_Z' in quiet.record,'Original P0/P0_Z must remain separate')
    require(quiet.record['derived_original_global_phase']['phase_independent_of_Z'],'Original collar phase must remain Z independent')
    # bridge_second previously failed the global dispatch context guard.
    bridge=oracle.density_frame(chart='bridge_second',Z=(-1,1),coordinate='1.831',N=2048)
    require(bridge.values is not None and all(value.ctx is c for group in bridge.values.values() for value in group.values()),
        'Actual original bridge ranges must bind to global arithmetic')
    bridge_q=bridge.record['actual_original_spatial_source']['original_q_slow_jet_source']
    require(bridge_q.get('original_source_ranges_directed_context_bridge') and bridge_q['original_source_ledger_unchanged'],
        'Actual native bridge context conversion must be explicit')
    require(bridge_q['conservative_range_conversion_not_added_source_correlations'],'Range conversion cannot invent correlations')
    result=dict(all_passed=True,source_family=oracle.source_family,**{current.GATE:True},
        original_chart_queries_attempted=17,enclosed_declared_original_queries=14,unresolved_declared_original_queries=3,
        explicit_remaining_original_source_frontier=reasons,
        inherited_role_namespace_root_checks=inherited.role_checks(role_owner,built,oracle),
        inherited_fail_closed_dispatch_checks=inherited.negative_checks(oracle,built,broad),
        additional_broad_original_role_dispatches=dispatched,
        original_selected_initial_collar_frame=quiet.record,original_native_bridge_global_context_frame=bridge.record,
        unresolved_native_collar_provider_not_replaced_by_analytic_flat_shortcut=True,
        actual_original_context_bridge_and_separate_P0_Z_checked=True,
        all17_whole_chart_or_continuous_route_ranges_admitted=False,full_factored_function_graph_evaluator_installed=False,
        actual_original_full_route_numerical_integrals_evaluated=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,
        input_hashes={**oracle.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        scope='Original role/context/namespace dispatch, all seventeen explicit provider outcomes, broad O2 frame, unresolved native initial-collar provider and directed native bridge context conversion. Three remaining interior-provider obstructions are retained; no whole-route numerical targets, controls, terminal closure/global N or recursion admission.')
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original signed-u factored dispatch focused check PASS;14 enclosed /3 unresolved',flush=True)
    return result
