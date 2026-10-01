"""Ordered reproduction of the compliant pressure/core/outer/heat family.

Stops on the first failed prerequisite. Each module writes its separate,
source-bound receipt; no legacy .01 file is edited. The angular stage solves
the functional angular/pressure correction. Actual ap selection, complete
outer assembly and temporal recursion remain unfinished.
"""
import argparse
import importlib
import json
from pathlib import Path

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
SOURCE=('pressure_source','core_transfer','core_transfer_check','finite_core','finite_core_check')
INNER=('core_uniform_bounds','frozen_angular_bounds','frozen_H_bounds',
       'physical_norm_family','pressure_Kp','tolerance_exit_parameters',
       'core_tail_admission','seed_inclusion','Cstar_envelope_transfer',
       'K1_ledger','global_exit_certificate','reference_join_bounds',
       'reference_join_check','five_defect_admission','five_moment_repair',
       'five_moment_repair_check')
OUTER=('outer_initial','outer_initial_check','outer_buffer','outer_buffer_check',
       'outer_pulse_map','outer_pulse_map_check','outer_angular_candidate',
       'outer_angular_candidate_check','exact_heat_component','exact_heat_component_check')
ANGULAR=('outer_angular_repair','outer_angular_repair_check')


def stages(stage):
    return {'source':SOURCE,'inner':INNER,'outer':OUTER,'angular':ANGULAR,
            'all':SOURCE+INNER+OUTER+ANGULAR}[stage]


def run(stage='all',list_only=False):
    selected=stages(stage)
    if list_only:
        for name in selected:print(PREFIX+name+'.py')
        return None
    completed=[]
    for name in selected:
        print('Build compliant stage:',name,flush=True)
        importlib.import_module(PREFIX+name).run()
        completed.append(name)
    result=dict(stages_completed_this_run=completed,
        source_relation='epsilon_source=.001*delta; epsilon_core=1/Lambda',
        legacy_source_preserved=True,full_NS_background_completed=False,
        actual_angular_pressure_corrections_completed='outer_angular_repair_check' in completed,actual_ap_selected=False,
        temporal_recursion=False,
        next_dependency='Complete corrected future swirl energy, select ap, compose all five outer primitives and prove the matched outer cone')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('source','inner','outer','angular','all'),default='all')
    parser.add_argument('--list',action='store_true',help='Print the ordered modules without running them')
    args=parser.parse_args()
    run(args.stage,args.list)
