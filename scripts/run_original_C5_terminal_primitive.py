"""Publish terminal C5 flat contraction and new directed primitive/source kernels."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/'experiments/root_st073'
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_finite_N_signed_source import load,sha,require,digest
from openai_ns_reconstruction.original_C5_directed_primitive import build as build_primitive,PARENT,PARENT_SHA256
from openai_ns_reconstruction.original_C5_terminal_flat import build as build_flat
from openai_ns_reconstruction.original_C5_terminal_history import build as build_terminal
from openai_ns_reconstruction.original_n1_complex_source_bounds import build as build_source_bounds

NAME='lei_ren_part1_paper_compliant_current_original_C5_terminal_primitive.json'
CODE=('src/openai_ns_reconstruction/original_C5_directed_primitive.py',
      'src/openai_ns_reconstruction/original_C5_terminal_flat.py',
      'src/openai_ns_reconstruction/original_C5_terminal_history.py',
      'src/openai_ns_reconstruction/original_n1_complex_source_bounds.py',
      'scripts/run_original_C5_terminal_primitive.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted joined-history parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():
            require(sha(base/name)==h,'Changed source pin: '+name);count+=1
    return count


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,
        'Different accepted parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes',
        'contract','checks','execution_seconds'}
    for key,value in p.items():
        if key not in mutable:require(r[key]==value,'Inherited accepted packet changed: '+key)
    flat=r['actual_C5_terminal_flat'];terminal=r['terminal_C5_signed_history']
    require(flat['passed'] is True and flat['outgoing_C_Z0_Z3_exact_zero'] is True,
        'Actual outgoing flat-loop function identity missing')
    require(terminal['passed'] is True and terminal['actual_flat_certificate_sha256']==digest(flat) and
        terminal['exact_controls']['passed'] is True and all(terminal['exact_controls']['checks'].values()),
        'Terminal source contraction failed')
    require(terminal['ordinary_Z_orders']==[0,1,2,3] and len(terminal['rows'])==4 and
        terminal['predecessor_initial_primitive_and_both_slow_integrals_retained'] is True and
        terminal['contraction_only_at_actual_transition_endpoint1'] is True and
        terminal['moving_history_unchanged'] is True,'A nonterminal history/source was discarded')
    for key in ('old_positive_caps_subtracted','old_graphs_or_error_views_modified',
        'numeric_signed_H_values_installed','actual_phase_chosen_or_materialized',
        'accepted_ancestor_or_quadrature_replayed'):
        require(terminal[key] is False,'Unsupported terminal completion: '+key)
    for node in terminal['exact_local_expression_nodes']:
        if node['operation']=='accepted_signed_expression_reference':
            require(node['report']==PARENT and node['report_sha256']==PARENT_SHA256,
                'Changed nested function source')
            graph=p
            for key in node['graph_path']:graph=graph[key]
            require(graph[node['expression_node']]['operation']==node['expression_operation'] and
                node['reference_is_function_not_numeric_enclosure'] is True,'Missing nested function expression')
    q=terminal['O3_quiet_functional']
    require(q['H_input_rows']==terminal['H_terminal_Z0_Z3'] and
        q['ordinary_Z1_and_H_Z3_retained'] is True and
        q['signed_H_or_pressure_increment_numerically_evaluated'] is False,'Quiet source scale/input lost')
    kernel=r['directed_C5_angle_primitive']
    require(kernel['passed'] is True and all(kernel['source_AST_guards'].values()) and
        kernel['numeric_kernel_controls']['passed'] is True and
        all(kernel['numeric_kernel_controls']['checks'].values()),'Directed angle primitive failed')
    for key in ('actual_root_mixed_enclosures_installed','actual_N_phase_enclosed',
        'actual_signed_C_or_slow_integral_values_installed','actual_signed_H_installed',
        'old_errors_or_caps_modified','accepted_owner_fixture_recurrence_or_quadrature_replayed'):
        require(kernel[key] is False,'Conditional primitive was promoted to an actual input value')
    bounds=r['n1_complex_actual_source_bounds']
    require(bounds['passed'] is True and bounds['actual_source_bounds_proved'] is True and
        bounds['actual_weighted_operator_bounds_proved'] is True and bounds['diagnostic_only'] is False and
        bounds['parent_report_sha256']==PARENT_SHA256 and bounds['parent_mapping_matches_file'] is True and
        bounds['missing_input_contract']==[], 'Actual complex source/operator bounds failed')
    require(all(v is True for k,v in bounds['gates'].items()
        if k not in ('solution_norm','radial_Taylor_remainder')),
        'Actual source/domain/delta/Cauchy-loss hypothesis missing')
    require(bounds['actual_solution_norm_proved'] is False and
        bounds['actual_radial_Taylor_remainder_certified'] is False and
        bounds['accepted_owner_checker_fixture_recurrence_quadrature_replayed'] is False and
        bounds['domain']['actual_collar_used_as_core_extension'] is False and
        bounds['green_input']['physical_Lambda_materialized'] is False,
        'Source/operator bound promoted to an unbound solution or physical tail')
    require(all(v is False for v in r['limits'].values()),'Global construction completion introduced')
    return dict(passed=True,all_inherited_source_and_error_packets_unchanged=True,
        outgoing_terminal_C_exact_zero=True,ordinary_Z0_Z3=True,
        actual_complex_source_operator_bounds=True,canonical_solution_norm_or_radial_tail=False,
        predecessor_Cin_and_both_signed_slow_integrals_retained=True,
        actual_signed_H_installed=False,global_five_moment_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    flat=build_flat(ROOT,p);terminal=build_terminal(ROOT,p,flat)
    primitive=build_primitive(ROOT,p);bounds=build_source_bounds(ROOT,p)
    r=copy.deepcopy(p);r.update(stage='terminal_C5_flat_contraction_and_directed_primitive',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_C5_terminal_flat=flat,terminal_C5_signed_history=terminal,
        directed_C5_angle_primitive=primitive,n1_complex_actual_source_bounds=bounds,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{name:sha(ROOT/name) for name in CODE}))
    r['contract']=dict(actual_outgoing_phase_primitive_zero_by_flat_source=True,
        only_terminal_C_removed_from_new_history=True,
        directed_angle_kernel_is_conditional_on_actual_mixed_root_enclosures=True,
        actual_signed_history_numeric_values_still_uninstalled=True,
        all_inherited_graphs_errors_and_target_consumption_unchanged=True,
        accepted_ancestors_replayed=False)
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
