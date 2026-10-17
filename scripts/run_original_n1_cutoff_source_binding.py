"""Bind the actual raw first-chart n1 function, cutoff histories and C_inner."""
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
from openai_ns_reconstruction.original_n1_first_raw_function import build as build_raw
from openai_ns_reconstruction.original_n1_first_volterra_remainder import build as build_tail
from openai_ns_reconstruction.original_n1_cutoff_history_identity import build as build_cutoff
from openai_ns_reconstruction.original_n1_m15_tail_transport import build as build_error,PARENT,PARENT_SHA256

NAME='lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json'
CODE=('src/openai_ns_reconstruction/original_n1_first_raw_function.py',
    'src/openai_ns_reconstruction/original_n1_first_volterra_remainder.py',
    'src/openai_ns_reconstruction/original_n1_cutoff_history_identity.py',
    'src/openai_ns_reconstruction/original_n1_m15_tail_transport.py',
    'scripts/run_original_n1_cutoff_source_binding.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted native collar parent changed')
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
            if not isinstance(row,dict):
                require(row is True,'Failed source pin summary');continue
            require(isinstance(row,dict) and isinstance(row.get('path'),str),'Invalid source pin row')
            n=row['path'];h=row.get('expected_sha256',row.get('sha256'))
            require(isinstance(h,str) and sha(ROOT/n)==h,'Invalid source pin: '+n)
            if n.startswith('experiments/root_st073/'):
                target=inputs;n=n.split('experiments/root_st073/',1)[1]
            else:target=impl
            require(target.get(n,h)==h,'Conflicting source pin: '+n);target[n]=h
    return inputs,impl


def source_binding(p,a,c,t,e):
    require(all(q['passed'] is True for q in (a,c,t,e)),'Passed raw/cutoff/tail/error source packets required')
    require(a['actual_first_raw_function_bound'] is True and
        a['actual_first_raw_all_word_convergence_bound'] is True and
        a['actual_first_raw_core_entry_bound'] is True and
        c['raw_input_gate']['function_source_gate']['passed'] is True and
        c['N1p_units']['passed'] is True and t['certified_raw_preconditioned_tail_majorant'] is True,
        'Actual first-raw source or physical cutoff pressure admission incomplete')
    h=c['history_identity']
    require(h['core_m15_C0_C3_nodes']==[52601,52602,52603,52604] and
        h['rho4_port_nodes']==[53125,53126] and h['first_increment_C0_C1_nodes']==[53137,53138] and
        h['combined_C0_C1_nodes']==[53141,53142] and
        h['C_inner_endpoint_C0_C3_nodes']==[53161,53162,53163,53164],
        'Actual signed C_inner moment handles changed')
    require(p['actual_core_pressure_m15_primitive']['actual_canonical_core_endpoint_handle_bound'] is True and
        p['actual_native_collar_transport']['native_defining_leading_function_equals_canonical_collar'] is True,
        'Accepted actual core and native leading source identities required')
    return dict(passed=True,parent_report_sha256=PARENT_SHA256,
        raw_packet_sha256=digest(a),cutoff_packet_sha256=digest(c),tail_packet_sha256=digest(t),
        error_transport_packet_sha256=digest(e),
        actual_C_inner_m15_source_function_bound=True,
        formula='C_inner15(Z)=M15_core(4,Z)+int_0^(sc/4) hb*rho*[2*epsilon*U_lead*eta*Uz1_raw-epsilon^2*rho*(2*F_lead*eta*F1_raw+N1p_physical)] ds',
        eta='1-sigma(4*(s/(sc/4))-2)',rho='4*exp(hb*s)',R='epsilon*rho',
        canonical_C0_C3_handles=h['C_inner_endpoint_C0_C3_nodes'],
        core_inlet_m15_C0_C1=h['rho4_port_nodes'],first_increment_C0_C1=h['first_increment_C0_C1_nodes'],
        ordinary_Z_orders=[0,1,2,3],nonzero_core_mean_and_pressure_history_retained=True,
        leading_N1p_survives_after_velocity_cutoff=True,P1_m13_is_not_m15=True,
        actual_finite_word_omitted_majorant_source_bound=True,
        actual_finite_word_values_evaluated=False,old_numeric_prefix_equals_word_prefix=False,
        actual_signed_C_inner_m15_value_evaluated=False,actual_m15_numeric_error_budget_instantiated=False,
        old_six_signed_terms_or_error_caps_removed=False,full_six_source_equality=False,
        full_five_moment_terminal_closure=False,accepted_ancestor_or_quadrature_replayed=False)


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Different accepted source parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes','contract','checks','execution_seconds'}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,'Inherited source or error packet changed: '+k)
    a=r['actual_first_raw_n1_function'];c=r['actual_n1_cutoff_history_identity']
    t=r['actual_first_raw_word_remainder'];e=r['actual_cutoff_m15_tail_transport']
    for packet,name in ((a,'raw function'),(c,'cutoff history'),(t,'raw tail'),(e,'m15 error transport')):
        require(packet['passed'] is True and packet['parent_report_sha256']==PARENT_SHA256,'Failed '+name)
    require(r['actual_C_inner_m15_source_binding']==source_binding(p,a,c,t,e),
        'Saved C_inner binding is not the actual source composition')
    require(c['raw_input_gate']['numerical_prefix_gate']['passed'] is False and
        e['actual_m15_numeric_error_budget_instantiated'] is False and
        e['actual_signed_m15_value_evaluated'] is False and
        e['old_pressure_or_moment_error_caps_removed'] is False and
        t['numeric_actual_first_chart_field_evaluated'] is False,
        'A finite approximation, norm cap or unresolved numerical moment was promoted')
    require(e['raw_packet_sha256']==digest(a) and e['cutoff_packet_sha256']==digest(c) and
        e['tail_packet_sha256']==digest(t) and all(v is False for v in r['limits'].values()),
        'Error transport inputs or inherited global scope changed')
    return dict(passed=True,actual_first_raw_all_word_function_and_core_entry_bound=True,
        raw_Akin_tail_provenance_bound_without_normalized_Kchi_substitution=True,
        actual_cutoff_mean_pressure_and_signed_C_inner_source_bound=True,ordinary_Z0_Z3=True,
        physical_N1p_and_Pstar_squared_once=True,conditional_raw_error_transport_formula=True,
        inherited_six_signed_terms_source_graphs_and_all_errors_unchanged=True,
        numerical_signed_C_inner_value=False,full_six_source_or_global_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    a=build_raw(ROOT,p);t=build_tail(ROOT,p);c=build_cutoff(ROOT,p,a);e=build_error(ROOT,p,a,c,t)
    binding=source_binding(p,a,c,t,e);inputs,impl=merge_source_pins(p,a,t,c,e)
    r=copy.deepcopy(p);r.update(stage='actual_first_raw_n1_cutoff_histories_and_C_inner_source_binding',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_first_raw_n1_function=a,actual_first_raw_word_remainder=t,
        actual_n1_cutoff_history_identity=c,actual_cutoff_m15_tail_transport=e,
        actual_C_inner_m15_source_binding=binding,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(actual_raw_first_function_and_own_cutoff_history=True,
        actual_signed_C_inner_m15_source_function_bound=True,raw_finite_word_omitted_majorant_provenance=True,
        common_inlet_pressure_error_cancellation_conditional_on_paired_reconstruction=True,
        actual_finite_word_values_or_signed_m15_numeric_values=False,
        old_six_signed_terms_and_all_error_caps_retained=True,
        global_five_moment_closure=False,accepted_ancestors_replayed=False)
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
