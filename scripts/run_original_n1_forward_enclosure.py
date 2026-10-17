"""Bind the infinite forward n1 response and enclose its raw coefficient fields."""
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
from openai_ns_reconstruction.original_n1_forward_enclosure import build as build_enclosure,collar_dependency_scope,PARENT,PARENT_SHA256
from openai_ns_reconstruction.original_n1_forward_analytic import build as build_forward
from openai_ns_reconstruction.original_n1_forward_scale_jets import build as build_scales

NAME='lei_ren_part1_paper_compliant_current_original_n1_forward_enclosure.json'
CODE=('src/openai_ns_reconstruction/original_n1_forward_analytic.py',
    'src/openai_ns_reconstruction/original_n1_forward_scale_jets.py',
    'src/openai_ns_reconstruction/original_n1_forward_enclosure.py',
    'scripts/run_original_n1_forward_enclosure.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted canonical norm/tail parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():
            require(sha(base/name)==h,'Changed source pin: '+name);count+=1
    return count


def scale_source_pins(p,scales):
    name=scales['source_report'];path=DIRECTORY/name
    require(sha(path)==scales['source_report_sha256'],'Scale log-amplitude source changed')
    report=load(path)
    inputs=dict(p['input_hashes']);implementations=dict(p['implementation_hashes'])
    def merge(target,items):
        for n,h in items.items():
            require(target.get(n,h)==h,'Conflicting transitive scale source pin: '+n)
            target[n]=h
    merge(inputs,report['input_hashes']);merge(inputs,{name:sha(path),PARENT:PARENT_SHA256})
    merge(implementations,report['implementation_hashes'])
    merge(implementations,{row['path']:row['sha256'] for row in scales['source_pins'].values()})
    pins(dict(input_hashes=inputs,implementation_hashes=implementations))
    return inputs,implementations


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,
        'Different source parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes',
        'contract','checks','execution_seconds'}
    for key,value in p.items():
        if key not in mutable:require(r[key]==value,'Inherited source/error packet changed: '+key)
    forward=r['n1_infinite_forward_analytic_response']
    scales=r['n1_forward_actual_scale_jets']
    for packet,label in ((forward,'forward'),(scales,'scale')):
        require(packet['passed'] is True and packet['parent_report_sha256']==PARENT_SHA256,
            'Actual '+label+' source binding failed')
        require(packet.get('checks') and all(v is True for v in packet['checks'].values()),
            'An actual '+label+' binding condition failed')
    require(all(forward.get(key) is True for key in ('actual_infinite_forward_analytic_response',
        'actual_forward_analytic_definition_bound','actual_forward_regular_uniqueness',
        'actual_forward_prefix_identity','actual_source_scale_conjugacy')),
        'Infinite forward construction was not source-bound')
    enclosure=r['n1_full_raw_coefficient_field_enclosures']
    require(enclosure['passed'] is True and enclosure['parent_report_sha256']==PARENT_SHA256 and
        enclosure['forward_binding_sha256']==digest(forward) and
        enclosure['scale_binding_sha256']==digest(scales) and
        enclosure['accepted_coefficient_packet_sha256']==digest(p['n1_actual_raw_prefix_binding']) and
        enclosure['accepted_raw_tail_packet_sha256']==digest(p['n1_actual_canonical_norm_tail']),
        'Enclosure detached from actual analytic function, scale or coefficient receipts')
    require(enclosure['exact_coefficient_polynomial_and_analytic_tail_combined'] is True and
        enclosure['full_raw_function_enclosures'] is True,'Full raw inclusion not established')
    require([v['center'] for v in enclosure['whole_real_radial_core_at_archive_centers']]==['-0.3','0','0.3'],
        'Three archived exact centers required')
    require(enclosure['endpoint_interval_in_whole_radial_box'] is True,
        'Core exit is not enclosed by whole radial core')
    for ev in (enclosure['whole_real_radial_core_at_archive_centers']+
            enclosure['core_exit_at_archive_centers']+enclosure['exact_axis_checks']):
        require(ev['passed'] is True and ev['common_raw_degree']==14 and len(ev['rows'])==24 and
            ev['ordinary_Z_orders']==[0,1,2,3] and ev['exact_interval_records_roundtripped']>0,
            'Whole radial core/ordinary-Z enclosure failed')
        require(ev['coefficient_midpoint_selected'] is False and
            ev['coefficient_decimal_approximation_selected'] is False and
            ev['enclosure_is_function_value'] is False and
            ev['interval_encloses_actual_infinite_raw_function'] is True,
            'An enclosure or archived approximation was promoted to a function value')
    require(all(enclosure[k] is False for k in ('arbitrary_real_Z_query_or_global_matching',
        'physical_R_derivative_enclosure','physical_NS_velocity_enclosure',
        'collar_six_source_equality_or_global_moment_closure','old_error_caps_modified',
        'temporal_n_ge_2','accepted_owner_builder_fixture_recurrence_or_quadrature_replayed')),
        'Unsupported physical/global completion')
    scope=r['actual_collar_dependency_scope_audit']
    require(scope['passed'] is True and
        scope['actual_leading_profile_dependencies']==['Phi','Psi','mean_Psi'] and
        scope['raw_response_nodes_absent_from_lower_graph'] is True and
        scope['pressure_rho_nodes_absent_from_lower_graph'] is True and
        scope['radial_pressure_rho_binding_and_shape_checked'] is True and
        scope['raw_n1_response_is_not_a_leading_Phi_Psi_input'] is True and
        scope['saved_graph_node_counts']==dict(lower=36594,volterra=49024,radial=52487),
        'Actual collar consumer dependency type not established')
    require(all(scope[k] is False for k in ('actual_raw_response_to_collar_inlet_identity',
        'actual_q_dependent_leading_collar_cross_producer_identity',
        'geometry_alone_proves_consumer_identity','six_signed_source_equality','source_or_owner_executed')),
        'A geometry/owner coincidence was promoted to a collar source identity')
    require(all(v is False for v in r['limits'].values()),'Global completion introduced')
    return dict(passed=True,all_inherited_source_and_error_packets_unchanged=True,
        new_infinite_forward_definition_equation_and_prefix_bound=True,
        actual_scale_ordinary_Z0_Z3_and_raw_function_enclosures=True,
        six_signed_core_lower_terms_retained=True,
        full_old_forward_collar_equality=False,global_five_moment_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    forward=build_forward(ROOT,p);scales=build_scales(ROOT,p)
    input_hashes,implementation_hashes=scale_source_pins(p,scales)
    enclosure=build_enclosure(ROOT,p,forward,scales)
    collar_scope=collar_dependency_scope(ROOT,p)
    r=copy.deepcopy(p);r.update(stage='infinite_forward_n1_response_and_directed_full_raw_enclosures',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        n1_infinite_forward_analytic_response=forward,n1_forward_actual_scale_jets=scales,
        n1_full_raw_coefficient_field_enclosures=enclosure,
        actual_collar_dependency_scope_audit=collar_scope,
        input_hashes=input_hashes,
        implementation_hashes=dict(implementation_hashes,**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(new_infinite_forward_function_from_actual_canonical_solution=True,
        ordinary_scale_jets_bound_to_same_source=True,
        exact_coefficient_intervals_and_analytic_radial_tail_combined=True,
        F0_exponential_kept_factored_with_proved_logarithmic_dominance=True,
        whole_real_radial_core_only_at_three_archived_Z_centers=True,
        no_coefficient_midpoint_or_error_cap_promoted_to_value=True,
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
