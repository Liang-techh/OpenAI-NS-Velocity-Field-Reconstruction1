"""Certify the actual O2/O3 phase seam and join its signed mean history."""
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
from openai_ns_reconstruction.original_C5_mean_seam import build as build_seam
from openai_ns_reconstruction.original_C5_joined_history import build as build_joined,PARENT,PARENT_SHA256
from openai_ns_reconstruction.original_n1_complex_radial_majorant import build as build_green

NAME='lei_ren_part1_paper_compliant_current_original_C5_seam_joined_history.json'
CODE=('src/openai_ns_reconstruction/original_C5_mean_seam.py',
    'src/openai_ns_reconstruction/original_C5_joined_history.py',
    'src/openai_ns_reconstruction/original_n1_complex_radial_majorant.py',
    'scripts/run_original_C5_seam_joined_history.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted radial/IBP parent changed')
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
                'O3_Nscaled_inlet_functional','n1_source_vector_identity',
                'n1_all_order_radial_conjugacy','actual_C5_signed_history_IBP'):
        require(r[key]==p[key],'Accepted archive/operator changed: '+key)
    seam=r['actual_C5_O2_O3_seam'];joined=r['joined_C5_signed_history']
    require(seam['passed'] is True and seam['source_only'] is True and
        seam['accepted_parent_sha256']==PARENT_SHA256 and seam['ordinary_Z_orders']==[0,1,2,3],
        'Actual seam source or ordinary derivative scope failed')
    require(all(seam[k] is True for k in ('arbitrary_issued_phase_B_identity',
        'anchored_primitive_identity','common_unwrapped_phase_identity')),
        'Independent-phase and common-phase seam identities required')
    require(all(seam['checks'].values()) and seam['source_pins_ok'] is True and
        all(v['passed'] is True for v in seam['AST'].values()) and
        seam['numeric_C_H_evaluated'] is False,'Actual seam source or numeric scope failed')
    inputs=seam['root_input_congruence'];ep=inputs['actual_endpoint_equalities']
    gr=inputs['actual_generic_and_root_equalities']
    fields=[row for rows in ep['field_rows'].values() for row in rows]
    radial=[row for rows in ep['radial_rows'].values() for row in rows]
    generic=[row for rows in gr['generic_rows'].values() for row in rows]
    roots=[row for rows in gr['root_rows'].values() for row in rows]
    require(inputs['passed'] is True and ep['passed'] is True and gr['passed'] is True and
        (len(fields),len(radial),len(generic),len(roots))==(20,8,20,32) and
        all(row['equal'] and row['equal_to_actual_incoming'] for row in fields) and
        all(row['equal'] for row in radial+generic+roots) and ep['radius']['equal'] is True,
        'Complete actual endpoint/generic/root input identities are missing')
    require(joined['passed'] is True and joined['source_only'] is True and
        joined['accepted_seam_certificate_sha256']==digest(seam) and
        joined['exact_controls']['passed'] is True and all(joined['exact_controls']['checks'].values()),
        'Joined signed expression source/controls failed')
    require(joined['ordinary_Z_orders']==[0,1,2,3] and len(joined['rows'])==4 and
        joined['original_graph_and_all_error_views_unchanged'] is True and
        joined['predecessor_history_and_outer_endpoints_retained'] is True and
        joined['both_signed_slow_integrals_retained'] is True,'History or endpoint was discarded')
    for key in ('numeric_signed_H_values_installed','physical_N_or_scale_materialized',
                'old_positive_caps_subtracted','global_five_moment_closure',
                'accepted_owner_fixture_recurrence_or_quadrature_replayed'):
        require(joined[key] is False,'Unsupported joined-history completion: '+key)
    for node in joined['exact_local_expression_nodes']:
        if node['operation']=='accepted_signed_expression_reference':
            require(node['report']==PARENT and node['report_sha256']==PARENT_SHA256,
                'Changed nested source namespace')
            graph=p
            for key in node['graph_path']:graph=graph[key]
            require(graph[node['expression_node']]['operation']==node['expression_operation'] and
                node['reference_is_function_not_numeric_enclosure'] is True,'Missing actual function reference')
    quiet=joined['O3_quiet_functional']
    require(quiet['H_input_rows']==joined['H_endpoint_Z0_Z3'] and
        quiet['H_already_N_scaled'] is True and quiet['ordinary_Z1_and_H_Z3_retained'] is True and
        quiet['signed_H_or_pressure_increment_numerically_evaluated'] is False,
        'Quiet functional received a different inlet or scale')
    require(all(v is False for v in r['limits'].values()),'Global construction completion introduced')
    green=r['n1_weighted_complex_Green_kernel']
    require(green['passed'] is True and green['actual_weighted_matrix_entries']['passed'] is True and
        len(green['actual_weighted_matrix_entries']['entries'])==16 and green['exact_controls']['passed'] is True,
        'Actual operator parity or exact weighted Green kernel failed')
    require(all(green[k] is False for k in ('actual_complex_source_sup_bounds_installed',
        'actual_directed_solution_norm_installed','actual_radial_Taylor_remainder_certified',
        'actual_collar_used_as_analytic_core_extension','accepted_owner_fixture_recurrence_or_quadrature_replayed',
        'physical_scale_or_N_materialized')),'Conditional Green kernel promoted to an actual norm/remainder')
    return dict(passed=True,active_discrepancy_views_unchanged=True,active_terms=178,
        same_source_internal_phase_pair_contracted=True,ordinary_Z0_Z3=True,
        predecessor_and_two_slow_integrals_retained=True,quiet_pressure_scales_separate=True,
        actual_signed_inlet_numeric_values=False,global_five_moment_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    seam=build_seam(ROOT,p);joined=build_joined(ROOT,p,seam);green=build_green(ROOT,p)
    r=copy.deepcopy(p);r.update(stage='same_source_C5_phase_seam_and_joined_signed_history',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_C5_O2_O3_seam=seam,joined_C5_signed_history=joined,
        n1_weighted_complex_Green_kernel=green,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{name:sha(ROOT/name) for name in CODE}))
    r['contract']=dict(seam_proved_as_arbitrary_phase_function=True,
        internal_endpoint_contraction_is_source_identity_not_cap_subtraction=True,
        original_IBP_graph_and_all_discrepancy_errors_unchanged=True,
        quiet_kernel_applied_to_new_joined_source=True,signed_inlet_numeric_value_still_uninstalled=True,
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
