"""Bind core pressure FTC and the signed leading-collar pressure partition."""
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
from openai_ns_reconstruction.original_n1_core_pressure_primitive import build as build_core
from openai_ns_reconstruction.original_leading_collar_identity import build as build_leading
from openai_ns_reconstruction.original_collar_pressure_partition import build as build_partition,PARENT,PARENT_SHA256

NAME='lei_ren_part1_paper_compliant_current_original_core_collar_pressure_identity.json'
CODE=('src/openai_ns_reconstruction/original_n1_core_pressure_primitive.py',
    'src/openai_ns_reconstruction/original_leading_collar_identity.py',
    'src/openai_ns_reconstruction/original_collar_pressure_partition.py',
    'scripts/run_original_core_collar_pressure_identity.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted forward enclosure parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed parent required');return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for n,h in r[key].items():
            require(sha(base/n)==h,'Changed source pin: '+n);count+=1
    return count


def merge_source_pins(p,*packets):
    inputs=dict(p['input_hashes']);impl=dict(p['implementation_hashes'])
    inputs[PARENT]=PARENT_SHA256
    for packet in packets:
        for row in packet['source_pins'].values():
            if not isinstance(row,dict) or not isinstance(row.get('path'),str):continue
            n=row['path'];h=row.get('expected_sha256',row.get('sha256'))
            require(isinstance(h,str) and sha(ROOT/n)==h,'Invalid source pin: '+n)
            if n.startswith('experiments/root_st073/'):
                target=inputs;n=n.split('experiments/root_st073/',1)[1]
            else:target=impl
            require(target.get(n,h)==h,'Conflicting source pin: '+n);target[n]=h
    pins(dict(input_hashes=inputs,implementation_hashes=impl))
    return inputs,impl


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256,'Different source parent')
    mutable={'stage','parent_report','parent_report_sha256','input_hashes','implementation_hashes','contract','checks','execution_seconds'}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,'Inherited source or error packet changed: '+k)
    core=r['actual_core_pressure_m15_primitive'];leading=r['actual_leading_collar_source_binding']
    partition=r['actual_collar_pressure_partition']
    for packet,name in ((core,'core pressure'),(leading,'leading collar'),(partition,'pressure partition')):
        require(packet['passed'] is True and packet['parent_report_sha256']==PARENT_SHA256 and
            packet.get('checks') and all(v is True for v in packet['checks'].values()),'Failed '+name+' source binding')
    require(core['actual_core_pressure_m15_functional_identity_bound'] is True and
        core['actual_raw_pressure_rho_equation_bound'] is True and
        core['actual_canonical_core_cumulative_binding'] is True and
        core['actual_canonical_core_endpoint_handle_bound'] is True and
        core['canonical_actual_core_endpoint_handle']['requested']==[53125,53126] and
        core['symbolic_identity']['passed'] is True,'Actual core pressure primitive incomplete')
    pressure=core['actual_raw_pressure_rho_equation']
    require(all(pressure[k] is True for k in ('raw_matrix_row4_literal_bound',
        'raw_matrix_field_meaning_bound','raw_matrix_B1_row4_absent','symbolic_pressure_row_factor',
        'symbolic_normalized_pressure_row_factor','symbolic_physical_core_density')),
        'Raw pressure RHS or physical core density is not actually source-bound')
    cutoff=core['accepted_inner_cutoff_binding']
    require(cutoff['passed'] is True and cutoff['checks']['core_m15_exact_product_sign_terms'] is True and
        cutoff['checks']['complete_history_not_reset_at_core_seam'] is True and
        cutoff['complete_history_not_reset_at_core_seam'] is True,
        'Actual signed integrands or incoming core history is not source-bound')
    require(partition['leading_collar_binding_sha256']==digest(leading) and
        partition['full_pressure_partition_for_same_leading_fields'] is True and
        partition['physical_Pstar_squared_canceled_once'] is True and
        partition['ordinary_Z0_Z1'] is True and partition['exact_algebra']['passed'] is True,
        'Same-leading source pressure partition incomplete')
    require(partition['signed_primitive_binding']['packet_sha256']==digest(p['canonical_collar_signed_primitive']) and
        partition['signed_primitive_binding']['consumed_once'] is True and
        partition['local_same_source_singular_q1_pair_contracted'] is True,
        'Signed primitive or common internal endpoint lost')
    require(all(partition[k] is False for k in ('native_D_over_R_bridge_equivalence_claimed',
        'saved_canonical_inner_plus_lower_or_ledger_contracted',
        'numerical_signed_A_minus_A0_value','first_chart_raw_cut_response_identity',
        'complete_forward_regular_plus_cut_coupling_identity','six_signed_source_equality',
        'six_terms_removed','global_moment_closure','accepted_ancestor_or_quadrature_replayed')),
        'Unsupported signed numerical/source contraction')
    require(partition['all_old_error_caps_retained'] is True and all(v is False for v in r['limits'].values()),
        'An inherited error or global scope changed')
    progress=r['six_source_identity_progress']
    require(progress['core_source_formula_and_actual_rho4_port_identified'] is True and
        progress['old_six_numeric_rows_and_errors_unchanged'] is True and
        progress['old_numeric_rows_equal_exact_functions'] is False and
        progress['complete_six_source_identity'] is False,'Unsupported six-source progress view')
    return dict(passed=True,actual_core_pressure_m15_FTC_C0_C1=True,
        actual_same_leading_regular_singular_pressure_identity=True,
        signed_A_minus_A0_primitive_reused_once=True,local_singular_q1_source_pair_contracted=True,
        all_inherited_source_graphs_six_terms_and_errors_unchanged=True,
        first_chart_raw_cut_response_identity=False,full_six_source_or_global_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    core=build_core(ROOT,p);leading=build_leading(ROOT,p)
    partition=build_partition(ROOT,p,leading)
    inputs,impl=merge_source_pins(p,core,leading)
    r=copy.deepcopy(p);r.update(stage='core_m15_pressure_primitive_and_same_leading_collar_pressure_partition',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        actual_core_pressure_m15_primitive=core,actual_leading_collar_source_binding=leading,
        actual_collar_pressure_partition=partition,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r['contract']=dict(core_pressure_FTC_IBP_and_actual_equation_bound=True,
        same_leading_N1p_regular_singular_partition_bound=True,
        old_forward_finite_rows_not_promoted_to_function_identity=True,
        old_six_signed_terms_and_all_error_caps_retained=True,
        signed_primitive_not_its_norm_cap_consumed_once=True,
        only_local_common_singular_endpoint_source_pair_contracted=True,
        global_five_moment_closure=False,accepted_ancestors_replayed=False)
    r['six_source_identity_progress']=dict(
        core_source_formula_and_actual_rho4_port_identified=True,canonical_core_m15_port=[53125,53126],
        exact_core_integrand_products_and_signs_bound=True,complete_core_seam_history_retained=True,
        same_leading_pressure_regular_singular_identity=True,
        accepted_signed_A_minus_A0_function_reused=True,
        native_bridge_D_over_R_to_ratio_Psi_identity=False,
        first_chart_raw_cut_response_identity=False,
        old_six_numeric_rows_and_errors_unchanged=True,old_numeric_rows_equal_exact_functions=False,
        complete_six_source_identity=False,numerical_signed_A_minus_A0_value=False,
        target_quiet_consumed_once=True,old_six_terms_removed=False)
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
