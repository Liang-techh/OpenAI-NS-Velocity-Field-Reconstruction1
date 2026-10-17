"""Build new all-order radial rules and actual C5 signed-history IBP ports."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/'experiments/root_st073'
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_finite_N_signed_source import load,sha,require
from openai_ns_reconstruction.original_n1_radial_conjugacy import prove_all_order_radial_conjugacy
from openai_ns_reconstruction.original_C5_mean_IBP_ports import build as build_IBP

PARENT='lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json'
PARENT_SHA256='2ac1269c3da29466067d1dfd50a92be321e85291b2634ffa6d2eafd679e9941e'
NAME='lei_ren_part1_paper_compliant_current_original_radial_IBP.json'
CODE=('src/openai_ns_reconstruction/original_n1_radial_conjugacy.py',
    'src/openai_ns_reconstruction/original_C5_mean_IBP_ports.py','scripts/run_original_radial_IBP.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted inlet/source parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():
            require(sha(base/name)==h,'Changed source pin: '+name);count+=1
    return count


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256 and
        r['source_family']==p['source_family'],'Different parent or construction family')
    for key in ('Rc178_discrepancy_views','quiet179_discrepancy_views','inherited_m15_packet',
                'O3_Nscaled_inlet_functional','n1_source_vector_identity'):
        require(r[key]==p[key],'Accepted archive/operator changed: '+key)
    radial=r['n1_all_order_radial_conjugacy'];ibp=r['actual_C5_signed_history_IBP']
    require(radial['passed'] is True and radial['source_executed'] is False and
        radial['accepted_receipts']['passed'] is True and
        radial['point_values_selected'] is False and radial['finite_N_numerical_calculation'] is False,
        'Radial formal identity or scope failed')
    require(all(v['passed'] is True for v in radial['source_pins'].values()),'Radial source pins failed')
    rules=radial['symbolic']['radial_rules']
    require(all(rules['checks'].values()) and rules['coefficient_uniqueness']['formal_solution_unique'] is True and
        radial['symbolic']['extra_order_budget']['passed'] is True,'All-order rules/budget failed')
    require(radial['old_pure_definitions_executed'] is False and radial['archive_source_executed'] is False and
        radial['symbolic']['inherited_matrix_conjugacy']['replayed'] is False and
        radial['symbolic']['inherited_source_map']['replayed'] is False and
        radial['symbolic']['triangularity']['checks']['passed'] is True,
        'Accepted algebra replay or triangularity scope changed')
    for key in ('full_function_equality','radial_Taylor_remainder','numeric_radial_tail',
                'temporal_n_ge_2_recursion','global_n1_solution','all_error_views_mutated'):
        require(radial['limitations'][key] is False,'Unsupported radial completion: '+key)
    require(ibp['passed'] is True,'Actual C5 IBP binding failed')
    require(all(v['passed'] is True for v in ibp['source_pins'].values()) and
        all(v['passed'] is True for v in ibp['AST'].values()) and
        ibp['source_anchor_checks']['passed'] is True and ibp['accepted_zero_mean_receipt']['passed'] is True,
        'Actual C5/slow source provenance failed')
    require(all(ibp[key] is False for key in ('numeric_value_available','signed_numeric_H_values_installed',
        'cap_contraction_or_improvement','physical_scale_materialized',
        'archived_owner_builder_checker_fixture_or_quadrature_executed','full_n1_or_terminal_closure',
        'original_sources_numerically_executed')),'IBP function was promoted to a numeric value/completion')
    require(ibp['scaled_history_update_control']['passed'] is True and
        all(ibp['scaled_history_update_control']['checks'].values()),'New transformed transport control failed')
    ports=ibp['source_ports'];nodes=ibp['exact_local_expression_nodes']
    for chart,original_upper in (('O2_buffer',4),('O3_transition',1)):
        c=ports[chart];require(len(c['rows'])==4 and c['endpoint_is_separate_from_moving_query'] is True,
            'Missing ordinary-Z endpoint rows')
        require(nodes[c['domain_upper_node']]['original_node']==original_upper,'Wrong actual chart endpoint')
        for k,(moving,endpoint) in enumerate(zip(c['H_out_moving_query_rows'],c['H_out_endpoint_rows'])):
            node=nodes[endpoint]
            require(node['operation']=='function_substitution' and node['expression']==moving and
                node['variable']==c['coordinate'] and node['value']==c['domain_upper_node'],
                'Moving query was not substituted at the actual outgoing endpoint')
            row=c['rows'][k];slow=nodes[row['C_y_slow']]
            require(slow['operation']=='current_original_C5_phase_integral_of_slow_B_y' and
                slow['body_role']=='actual_slow_B_y_integrand' and slow['fixed_phase'] is True and
                row['physical_Jacobian_factor_count']==1 and row['selected_N_division_count']==1,
                'Wrong slow-primitive schema or scale/measure factor')
    require([row['incoming_H'] for row in ports['O3_transition']['rows']]==ports['O2_buffer']['H_out_endpoint_rows'] and
        ibp['H']['O3_transition_outgoing_Z0_to_Z3']==ports['O3_transition']['H_out_endpoint_rows'],
        'History composition retained a moving query instead of an endpoint')
    for node in ibp['exact_local_expression_nodes']:
        if node['operation']=='original_function_reference':
            require(node['report_sha256']==r['input_hashes'][node['report']],'Unpinned original IBP source')
    require(all(v is False for v in r['limits'].values()),'Global construction completion introduced')
    return dict(passed=True,active_discrepancy_views_unchanged=True,active_terms=178,
        all_order_radial_formal_conjugacy=True,radial_Taylor_remainder=False,
        actual_C5_history_IBP_source_ports=True,actual_signed_inlet_numeric_values=False,
        global_five_moment_closure=False,accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    radial=prove_all_order_radial_conjugacy(ROOT);ibp=build_IBP(ROOT,p)
    r=copy.deepcopy(p);r.update(stage='all_order_n1_radial_rules_and_actual_C5_signed_history_IBP',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        n1_all_order_radial_conjugacy=radial,actual_C5_signed_history_IBP=ibp,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{name:sha(ROOT/name) for name in CODE}))
    r['contract']=dict(actual_C5_B_and_fixed_phase_slow_B_y_bound=True,
        incoming_history_and_both_endpoints_retained=True,exact_N_division_once_after_scaling=True,
        source_bound_functions_are_not_numeric_H=True,radial_rules_are_not_temporal_n_ge_2=True,
        radial_analytic_remainder_unproved=True,Rc178_and_all_errors_unchanged=True,
        target_quiet_or_viscous_packets_not_consumed_again=True,accepted_ancestors_replayed=False)
    r['checks']=dict(passed=True,packet_invariants=packet_check(r,p),accepted_ancestors_replayed=False)
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True,'Saved passed report required')
    return dict(passed=True,pins=pins(r),packet_invariants=packet_check(r,parent()),serialized_only=True,
        owner_fixture_recurrence_or_quadrature_replayed=False)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-report',action='store_true');group.add_argument('--inspect-report',action='store_true')
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n');tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=r['execution_seconds'],bytes=path.stat().st_size,sha256=sha(path))))
    else:print(json.dumps(inspect(path),indent=2))
