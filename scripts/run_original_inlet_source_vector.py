"""Build new actual O3 inlet coefficients and n1 source-vector algebra."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/'experiments/root_st073'
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_finite_N_signed_source import load,sha,digest,require
from openai_ns_reconstruction.original_O3_inlet_functional import PARENT,PARENT_SHA256,build,exact_controls
from openai_ns_reconstruction.original_n1_source_vector_identity import prove_source_vector_identity

NAME='lei_ren_part1_paper_compliant_current_original_inlet_source_vector.json'
CODE=('src/openai_ns_reconstruction/original_O3_inlet_functional.py',
    'src/openai_ns_reconstruction/original_n1_source_vector_identity.py',
    'scripts/run_original_inlet_source_vector.py')


def parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,'Accepted Rc parent changed')
    p=load(DIRECTORY/PARENT);require(p['checks']['passed'] is True,'Passed Rc parent required');return p


def pins(r):
    n=0
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,h in r[key].items():require(sha(base/name)==h,'Changed source pin: '+name);n+=1
    return n


def packet_check(r,p):
    require(r['parent_report']==PARENT and r['parent_report_sha256']==PARENT_SHA256 and
        r['source_family']==p['source_family'],'Different parent or construction family')
    for key in ('Rc178_discrepancy_views','quiet179_discrepancy_views','inherited_m15_packet'):
        require(r[key]==p[key],'Accepted discrepancy or original archive changed')
    inlet=r['O3_Nscaled_inlet_functional'];n1=r['n1_source_vector_identity']
    require(inlet['passed'] is True and inlet['source_executed'] is False and
        all(v['passed'] is True for v in inlet['AST'].values()),'Actual inlet source provenance missing')
    require(inlet['unknown_H_not_set_to_zero'] is True and inlet['finite_N_error_logs_not_used_as_H_values'] is True and
        inlet['N_scaled_history_values_or_signed_centers_evaluated'] is False,'Unknown inlet was assigned a cap or value')
    bg=inlet['physical_background_Q0'];coeff=inlet['actual_directed_coefficient_enclosures']
    require(bg['strictly_negative'] is True and bg['ordinary_Z1_Z2_exact_zero'] is True and
        bg['diagnostic_agreement_is_not_function_identity_proof'] is True,'Background pressure provenance lost')
    require(coeff['passed'] is True and len(coeff['cells'])==2 and
        coeff['full_Z_domain_or_terminal_closure_claim'] is False,'Coefficient cell/domain scope changed')
    for cell in coeff['cells']:
        require(cell['actual_H_values_evaluated'] is False and cell['full_cell_directed_enclosures'] is True and
            len(cell['linear_C0_C1'])==2 and all(len(v)==4 for v in cell['linear_C0_C1']) and
            len(cell['quadratic_C0_C1'])==2 and
            all(len(v)==4 and all(len(row)==4 for row in v) for v in cell['quadratic_C0_C1']),
            'Incomplete ordinary-C1 coefficient tensors')
    for node in inlet['exact_local_expression_nodes']:
        if node['operation']=='original_function_reference':
            require(node['report_sha256']==r['input_hashes'][node['report']],'Unpinned actual function port')
    require(n1['passed'] is True and n1['source_executed'] is False and
        all(n1['symbolic']['checks'].values()) and all(v is False for v in n1['limits'].values()),
        'n1 source algebra failed or full-solution completion introduced')
    require(n1['accepted_common_fixed_point']['passed'] is True and
        n1['accepted_common_fixed_point']['common_input_hashes']['passed'] is True and
        n1['accepted_common_fixed_point']['scaled_S_jet_inclusion']['passed'] is True and
        all(n1['symbolic']['inherited_premises'].values()) and
        all(v['passed'] is True for v in n1['source_pins'].values()),'Actual n1 input functions are not bound')
    require(all(v is False for v in r['limits'].values()),'Global construction completion introduced')
    return dict(passed=True,active_discrepancy_views_unchanged=True,active_terms=178,
        signed_background_pressure_enclosure=True,source_bound_unknown_H=True,
        actual_coefficient_whole_cells=2,ordinary_linear_coefficients=16,ordinary_quadratic_coefficients=64,
        universal_n1_source_vector_algebra=True,complete_n1_solution_equality=False,
        actual_finite_N_signed_inlet_values=False,global_five_moment_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)


def run():
    start=time.perf_counter();p=parent();pins(p)
    inlet=build(ROOT,p);controls=exact_controls();n1=prove_source_vector_identity(ROOT)
    r=copy.deepcopy(p);r.update(stage='actual_O3_background_and_Nscaled_inlet_functional_with_n1_source_vector',
        parent_report=PARENT,parent_report_sha256=PARENT_SHA256,O3_Nscaled_inlet_functional=inlet,n1_source_vector_identity=n1,
        input_hashes=dict(p['input_hashes'],**{PARENT:PARENT_SHA256}),
        implementation_hashes=dict(p['implementation_hashes'],**{name:sha(ROOT/name) for name in CODE}))
    r['contract']=dict(actual_background_coefficient_enclosure_used=True,unknown_signed_H_not_numerically_installed=True,
        physical_Pstar_over_N_sectors_separate=True,ordinary_Z1_and_all_nonlinear_terms_retained=True,
        Rc178_and_all_error_terms_unchanged=True,target_or_quiet_not_subtracted_again=True,
        n1_algebra_does_not_admit_full_radial_solution_or_six_term_zero=True,accepted_ancestors_replayed=False)
    r['checks']=dict(passed=True,new_inlet_exact_controls=controls,packet_invariants=packet_check(r,p),
        accepted_ancestors_replayed=False)
    r['execution_seconds']=round(time.perf_counter()-start,3);return r


def inspect(path):
    r=load(path);require(r['checks']['passed'] is True and r['checks']['new_inlet_exact_controls']['passed'] is True,
        'Saved passed controls required')
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
