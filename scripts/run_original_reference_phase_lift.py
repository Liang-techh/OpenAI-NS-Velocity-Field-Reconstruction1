"""Produce or inspect the new source-bound reference periodic-lift packet."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_core_numeric import DIRECTORY,PREFIX
from openai_ns_reconstruction.original_micro_viscous_m15 import OriginalMicroViscousM15
from openai_ns_reconstruction.original_reference_phase_lift import build_phase_lift_certificate
from openai_ns_reconstruction.original_phase_graph_guards import build_phase_graph_guards

PARENT=PREFIX+'compliant_current_original_reference_physical_germ.json'
PARENT_SHA='a29bf10e7fa2fe6f17cc39505790596b6864f788061d8ee36fc68e6dfc18769e'
NAME=PREFIX+'compliant_current_original_reference_phase_lift.json'
CODE=('src/openai_ns_reconstruction/original_reference_phase_lift.py',
      'src/openai_ns_reconstruction/original_phase_graph_guards.py',
      'scripts/check_original_phase_lift.py','scripts/run_original_reference_phase_lift.py')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def require(ok,message):
    if not ok:raise ValueError(message)


def load_parent():
    path=DIRECTORY/PARENT
    require(sha(path)==PARENT_SHA,'Pinned accepted physical germ required')
    p=json.loads(path.read_text(encoding='utf-8'))
    require(p['checks']['passed'] and p['checks']['no_inlet_B_y_or_RQy_alias'],'Accepted unchanged ledger required')
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,digest in p[key].items():require(sha(base/name)==digest,'Changed source '+name)
    return p


def check_packet(data,parent):
    require(data['parent_report']==PARENT and data['parent_report_sha256']==PARENT_SHA,'Parent changed')
    require(data['source_family']==parent['source_family'],'Source family changed')
    cert=data['phase_lift_certificate']
    require(cert['passed'] and cert['global_periodic_lift_proved'] and
            cert['mixed_phase_parameter_jets_through_order2_proved'],'Periodic theorem missing')
    require(cert['graph_guards']['passed'] and cert['physical_seam_domain']['passed'],'Actual source binding missing')
    require(len(cert['physical_seam_domain']['rows'])==8 and
            all(r['passed'] for r in cert['physical_seam_domain']['rows']),'Exact root domain changed')
    require(cert['formal_identities']['passed'] and
            cert['formal_identities']['derivative_row_count']==108 and
            all(r['exact_zero'] for r in cert['formal_identities']['derivative_rows']),'Mixed derivative theorem missing')
    require(cert['source_only'] and not cert['numeric_caps_used'] and
            not cert['B_y_alias_authorized'],'Unproved native B_y alias')
    require(data['unchanged_m15_ledger']==parent['unchanged_m15_ledger'],'Ledger changed')
    require(data['focused_fixture']['passed'] and data['graph_negative_controls']['passed'],'Focused check failed')
    for key,value in data['limits'].items():require(value is False,'Open construction layer overclaimed: '+key)
    return dict(passed=True,source_only_smooth_periodic_lift=True,ordinary_Z1_physical_y_mixed_jets_admitted=True,
        raw_213_reduced_199_and_target_unchanged=True,all_old_error_budgets_unchanged=True,
        direct_modified_Uy_fingerprint=False,actual_native_B_y_naturality_proved=False,
        no_new_trace_alias=True,no_error_cap_subtraction=True)


def graph_controls(nodes,loop):
    original=build_phase_graph_guards(nodes,loop)
    require(original['passed'],'Unchanged graph rejected')
    changed=copy.deepcopy(loop)
    changed['charts']['reference']['actual_phase']=0
    wrong_phase=not build_phase_graph_guards(nodes,changed)['passed']
    require(wrong_phase,'Wrong actual phase accepted')
    inverse=loop['charts']['restoration']['loop']['inverse'][0]
    changed_nodes=list(nodes)
    changed_nodes[inverse]=dict(nodes[inverse],target_phase=0)
    wrong_inverse=not build_phase_graph_guards(changed_nodes,loop)['passed']
    require(wrong_inverse,'Wrong inverse target accepted')
    changed=copy.deepcopy(loop)
    changed['charts']['restoration']['loop']['inverse'][0]=0
    require(not build_phase_graph_guards(nodes,changed)['passed'],'Foreign inverse handle accepted')
    return dict(passed=True,wrong_actual_phase_rejected=True,
        wrong_inverse_target_rejected=True,foreign_inverse_handle_rejected=True)


def produce():
    start=time.perf_counter();parent=load_parent()
    # Archive-only hydration: no accepted owner/producer constructor is run.
    loader=object.__new__(OriginalMicroViscousM15)
    loader.directory=DIRECTORY;loader.hashes=dict(parent['input_hashes']);loader.family=parent['source_family']
    modified=loader.archive('modified_Q4_pressure_C3_functions')
    loop=loader.archive('prescribed_loop_C5_functions')
    nodes=modified['exact_graph_nodes'];prefix=loop['exact_graph_nodes']
    require(nodes[:len(prefix)]==prefix,'Original loop prefix changed')
    bindings=loop['prescribed_loop_C5_bindings']
    cert=build_phase_lift_certificate(nodes,bindings,parent['source_guards'],ROOT)
    require(cert['passed'],'New periodic-lift certificate failed')
    from check_original_phase_lift import run as fixture
    inputs=dict(loader.hashes);inputs[PARENT]=PARENT_SHA
    implementations=dict(parent['implementation_hashes'])
    implementations.update({name:sha(ROOT/name) for name in CODE})
    data=dict(source_family=parent['source_family'],parent_report=PARENT,parent_report_sha256=PARENT_SHA,
        phase_lift_certificate=cert,unchanged_m15_ledger=parent['unchanged_m15_ledger'],
        focused_fixture=fixture(),graph_negative_controls=graph_controls(nodes,bindings),
        input_hashes=inputs,implementation_hashes=implementations,
        contract=dict(accepted_ancestor_constructors_not_run=True,actual_fast_phase_not_assumed_integer=True,
            old_input_germ_and_m15_ledger_unchanged=True,native_chart_endpoint_is_not_fast_phase_endpoint=True),
        limits=dict(actual_modified_B_y_naturality=False,direct_modified_Uy_fingerprint=False,
            new_RQy_alias=False,new_budget_contraction=False,numeric_Rp_seam=False,
            global_five_moment_common_N_closure=False,exact_heat_and_energy=False,
            admissible_stress_flat_remainder=False,temporal_n_recursion_and_summation=False,
            mean_oscillatory_correction=False,independent_corrected_Cartesian_NS=False))
    data['checks']=check_packet(data,parent);data['execution_seconds']=round(time.perf_counter()-start,3)
    return data


def inspect(path):
    p=load_parent();data=json.loads(path.read_text(encoding='utf-8'))
    for base,key in ((DIRECTORY,'input_hashes'),(ROOT,'implementation_hashes')):
        for name,digest in data[key].items():require(sha(base/name)==digest,'Changed packet pin '+name)
    require(data['checks']==check_packet(data,p),'Saved checks changed')
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    data=inspect(path) if args.inspect_report else produce()
    if args.write_report:
        require(not args.inspect_report,'Inspection cannot rewrite accepted report')
        tmp=path.with_suffix('.json.tmp')
        tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
        tmp.replace(path)
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2))
