"""Bind the actual canonical n1 solution and certify its raw radial Taylor tail."""
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
from openai_ns_reconstruction.original_n1_canonical_norm_tail import build as build_tail,PARENT,PARENT_SHA256
from openai_ns_reconstruction.original_n1_canonical_solution_binding import build as build_canonical
from openai_ns_reconstruction.original_n1_raw_prefix_binding import build as build_prefix
from openai_ns_reconstruction.original_C5_buffer_swirl_jet import build as build_buffer

NAME='lei_ren_part1_paper_compliant_current_original_n1_canonical_norm_tail.json'
CODE=('src/openai_ns_reconstruction/original_n1_canonical_solution_binding.py',
    'src/openai_ns_reconstruction/original_n1_raw_prefix_binding.py',
    'src/openai_ns_reconstruction/original_n1_canonical_norm_tail.py',
    'src/openai_ns_reconstruction/original_C5_buffer_swirl_jet.py',
    'scripts/run_original_n1_canonical_norm_tail.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted terminal/source parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():
            require(sha(base/name)==h,'Changed source pin: '+name);count+=1
    return count


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,
        'Different source parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes',
        'contract','checks','execution_seconds'}
    for key,value in p.items():
        if key not in mutable:require(r[key]==value,'Inherited accepted packet changed: '+key)
    canonical=r['n1_actual_canonical_solution_binding'];prefix=r['n1_actual_raw_prefix_binding']
    for packet,name in ((canonical,'canonical'),(prefix,'prefix')):
        require(packet['passed'] is True and packet['parent_report_sha256']==PARENT_SHA256 and
            all(v is True for v in packet['checks'].values()),'Actual '+name+' source binding failed')
    tail=r['n1_actual_canonical_norm_tail']
    require(tail['passed'] is True and tail['parent_report_sha256']==PARENT_SHA256 and
        tail['actual_canonical_binding_sha256']==digest(canonical) and
        tail['actual_raw_prefix_binding_sha256']==digest(prefix) and
        tail['accepted_source_bound_packet_sha256']==digest(p['n1_complex_actual_source_bounds']) and
        tail['accepted_Green_packet_sha256']==digest(p['n1_weighted_complex_Green_kernel']),
        'Solution norm/tail is not bound to the actual source receipts')
    require(tail['actual_canonical_regular_solution_norm_proved'] is True and
        tail['actual_canonical_raw_radial_Taylor_remainder_certified'] is True and
        tail['raw_base_prefix_degree']==14 and tail['weighted_prefix_degree']==13 and
        tail['raw_odd_divided_by_x_prefix_degree']==13 and
        tail['maximal_value_prefix_degrees']==dict(raw_base=16,weighted=15,raw_odd_divided_by_x=15) and
        [v['ordinary_Z_order'] for v in tail['ordinary_Z0_Z3_tail_rows']]==[0,1,2,3],
        'Actual solution/tail or derivative degree failed')
    require(tail['exact_interval_records_roundtripped']>0 and
        tail['axial_domain']['real_center_interval_open']==['-1','1'] and
        tail['axial_domain']['target_in_solution_disk'] is True,
        'Ordinary-Z Cauchy domain or exact interval serialization failed')
    for key in ('computed_word_partial_sum_or_source_point_values','approximate_coefficient_errors_combined',
        'full_forward_response_function_identity','q_dependent_collar_bridge',
        'six_core_lower_terms_removed','global_five_moment_closure','temporal_n_ge_2',
        'actual_N_or_physical_scale_materialized','old_error_caps_modified',
        'accepted_owner_fixture_recurrence_or_quadrature_replayed'):
        require(tail[key] is False,'Unsupported canonical/physical completion: '+key)
    buffer=r['actual_O2_buffer_swirl_mixed_jet']
    require(buffer['passed'] is True and buffer['accepted_parent_sha256']==PARENT_SHA256 and
        all(v is True for v in buffer['checks'].values()) and buffer['errors']==[] and
        buffer['ordinary_Z_orders']==[0,1,2,3] and buffer['mixed_y_orders']==[0,1] and
        len(buffer['evaluation']['E_rows'])==8 and buffer['exact_interval_records_roundtripped']>0 and
        buffer['leading_background_swirl_not_modified_E_N'] is True,
        'Actual buffer leading swirl mixed-jet source failed')
    require(all(buffer[k] is False for k in ('numeric_C_H_evaluated','numeric_p2_mixed_evaluated',
        'numeric_q_dstar_evaluated','cap_contraction_or_improvement',
        'archived_owner_builder_checker_fixture_or_quadrature_executed')),
        'One swirl input promoted to a full signed inlet solution')
    require(all(v is False for v in r['limits'].values()),'Global construction completion introduced')
    return dict(passed=True,all_inherited_source_and_error_packets_unchanged=True,
        actual_canonical_equation_axis_uniqueness_and_prefix_bound=True,
        canonical_raw_solution_norm_and_radial_tail=True,ordinary_Z0_Z3=True,
        actual_buffer_leading_E_Z0_Z3_y0_y1_enclosures=True,
        full_forward_response_or_collar_identity=False,numeric_signed_H=False,
        global_five_moment_closure=False,accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    canonical=build_canonical(ROOT,p);prefix=build_prefix(ROOT,p)
    tail=build_tail(ROOT,p,canonical,prefix)
    buffer=build_buffer(ROOT,p)
    r=copy.deepcopy(p);r.update(stage='actual_n1_canonical_regular_solution_norm_and_raw_radial_tail',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        n1_actual_canonical_solution_binding=canonical,n1_actual_raw_prefix_binding=prefix,
        n1_actual_canonical_norm_tail=tail,
        actual_O2_buffer_swirl_mixed_jet=buffer,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(canonical_regular_solution_and_axis_bound=True,
        actual_complex_bounds_consumed_without_recomputation=True,
        raw_value_degree16_and_common_Z0_Z3_degree14_prefix_bound=True,
        ordinary_Z_remainder_uses_fresh_local_half_tube_gap=True,
        actual_buffer_swirl_mixed_source_component_evaluated=True,
        no_coefficient_center_or_error_cap_promoted_to_value=True,
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
