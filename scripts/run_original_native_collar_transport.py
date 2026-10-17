"""Prove the actual native first-collar leading transport source identity."""
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
from openai_ns_reconstruction.original_native_core_angular_direction import build as build_angular
from openai_ns_reconstruction.original_native_core_axial_direction import build as build_axial
from openai_ns_reconstruction.original_native_collar_transport import build as build_transport,PARENT,PARENT_SHA256

NAME='lei_ren_part1_paper_compliant_current_original_native_collar_transport.json'
CODE=('src/openai_ns_reconstruction/original_native_core_angular_direction.py',
    'src/openai_ns_reconstruction/original_native_core_axial_direction.py',
    'src/openai_ns_reconstruction/original_native_collar_transport.py',
    'scripts/run_original_native_collar_transport.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted core/collar pressure parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for n,h in r[key].items():
            require(sha(base/n)==h,'Changed source pin: '+n);count+=1
    return count


def merge_source_pins(p,*packets):
    inputs=dict(p['input_hashes']);impl=dict(p['implementation_hashes']);inputs[PARENT]=PARENT_SHA256
    for packet in packets:
        for row in packet['source_pins'].values():
            require(isinstance(row,dict) and row.get('passed') is True,'Invalid source pin row')
            n=row['path'];h=row.get('expected_sha256',row.get('sha256'))
            require(isinstance(h,str) and sha(ROOT/n)==h,'Invalid source pin: '+n)
            if n.startswith('experiments/root_st073/'):
                target=inputs;n=n.split('experiments/root_st073/',1)[1]
            else:target=impl
            require(target.get(n,h)==h,'Conflicting source pin: '+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Different accepted source parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes','contract','checks','execution_seconds'}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,'Inherited source or error packet changed: '+k)
    a=r['actual_native_core_angular_direction'];x=r['actual_native_core_axial_direction']
    t=r['actual_native_collar_transport']
    for packet,name in ((a,'angular direction'),(x,'axial drive'),(t,'collar transport')):
        require(packet['passed'] is True and packet['parent_report_sha256']==PARENT_SHA256 and
            packet.get('checks') and all(v is True for v in packet['checks'].values()),'Failed '+name)
    require(t['angular_packet_sha256']==digest(a) and t['axial_packet_sha256']==digest(x) and
        t['accepted_leading_packet_sha256']==digest(p['actual_leading_collar_source_binding']),
        'Transport is not bound to the actual direction packets')
    require(t['native_defining_leading_function_equals_canonical_collar'] is True and
        t['actual_native_D_over_R_and_total_drive_bound'] is True and
        t['actual_quotient_equals_ratio'] is True and t['ordinary_Z_orders']==[0,1,2,3] and
        t['exact_algebra']['passed'] is True and all(t['exact_algebra']['checks'].values()),
        'Actual source-level first-collar identity incomplete')
    require(t['domain']['first_window_only'] is True and t['domain']['alpha_identically_one'] is True and
        t['domain']['later_modification_identically_zero'] is True and
        t['quotient_conventions']['extra_exponential_after_evolved_field'] is False,
        'Wrong source domain, quotient or first/second-window scope')
    require(t['comparison_core_primitive_identity'] is True and
        t['native_actual_six_histories_equal_comparison'] is False and
        t['actual_own_mean_and_pressure_equal_canonical_actual_histories'] is True and
        t['history_distinction']['no_actual_history_replaced_by_comparison_history'] is True,
        'Actual histories must use the evolved fields, independently of the comparison histories')
    require(all(t[k] is False for k in ('old_finite_numeric_producer_is_exact_function',
        'raw_n1_cut_response_identity','full_C_inner_m15_identity','numerical_signed_A_minus_A0_value',
        'saved_six_terms_or_error_ledger_contracted','six_signed_source_equality','global_five_moment_closure',
        'accepted_ancestor_or_quadrature_replayed','later_second_window_or_global_native_identity','full_NS_validation')),
        'Unsupported source, numerical or global identity')
    require(t['all_old_error_caps_retained'] is True and all(v is False for v in r['limits'].values()),
        'An inherited error or global scope changed')
    return dict(passed=True,actual_same_core_native_angular_and_three_scale_axial_directions=True,
        actual_first_window_comparison_equals_core=True,native_actual_quotient_equals_ratio=True,
        leading_F_U_mean_and_pressure_same_functions_ordinary_Z0_Z3=True,
        inherited_source_graphs_six_signed_terms_and_errors_unchanged=True,
        raw_n1_cut_response_or_full_C_inner_identity=False,full_six_source_or_global_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    a=build_angular(ROOT,p);x=build_axial(ROOT,p);t=build_transport(ROOT,p,a,x)
    inputs,impl=merge_source_pins(p,a,x,t)
    r=copy.deepcopy(p);r.update(stage='native_same_core_directions_and_actual_first_collar_transport_identity',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_native_core_angular_direction=a,actual_native_core_axial_direction=x,
        actual_native_collar_transport=t,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(same_core_native_radial_directions_bound=True,
        native_defining_first_collar_leading_functions_equal_canonical=True,
        alpha_one_and_comparison_core_moments_bound=True,pointwise_quotient_equals_ratio=True,
        nonzero_mean_and_pressure_history_retained=True,ordinary_Z0_Z3=True,
        old_finite_rows_not_promoted_to_exact_function_values=True,
        old_six_signed_terms_and_all_error_caps_retained=True,global_five_moment_closure=False,
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
