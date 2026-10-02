"""Ordered reproduction of the compliant pressure/core/outer/heat family.

Stops on the first failed prerequisite. Each module writes its separate,
source-bound receipt; no legacy .01 file is edited. The angular stage solves
the functional angular/pressure correction. The energy stage supplies the
complete future swirl energy and selects actual ap/c1/c2 functions. The pulse
stage installs the full O.4 interval in segmented coordinates. The postpulse
stage composes corrected outer fields and primitives through Gamma heat.
The closure stage composes absolute pressure/angular source identities and
admits all five terminal moments. C4/cone and temporal recursion remain unfinished.
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
ENERGY=('future_swirl_energy','future_swirl_energy_check',
        'axial_amplitude_selection','axial_amplitude_selection_check')
PULSE=('axial_pulse_field','axial_pulse_field_check')
POSTPULSE=('corrected_outer_field','corrected_outer_field_check')
CLOSURE=('absolute_moment_closure','absolute_moment_closure_check')


def stages(stage):
    return {'source':SOURCE,'inner':INNER,'outer':OUTER,'angular':ANGULAR,'energy':ENERGY,'pulse':PULSE,'postpulse':POSTPULSE,'closure':CLOSURE,
            'all':SOURCE+INNER+OUTER+ANGULAR+ENERGY+PULSE+POSTPULSE+CLOSURE}[stage]


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
        actual_angular_pressure_corrections_completed=any(name in completed for name in ('outer_angular_repair_check','axial_amplitude_selection_check','axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check')),
        actual_ap_selected=any(name in completed for name in ('axial_amplitude_selection_check','axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check')),
        actual_O4_partial_pulse_field_installed=any(name in completed for name in ('axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check')),
        corrected_post_Rv_fields_and_primitives_installed=any(name in completed for name in ('corrected_outer_field_check','absolute_moment_closure_check')),
        all_five_absolute_terminal_moments_admitted='absolute_moment_closure_check' in completed,
        temporal_recursion=False,
        next_dependency='Recover higher axial derivatives and certify full interface/C4/cone bounds; build admissible stress and flat remainder, then genuine n-dependent temporal recursion')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('source','inner','outer','angular','energy','pulse','postpulse','closure','all'),default='all')
    parser.add_argument('--list',action='store_true',help='Print the ordered modules without running them')
    args=parser.parse_args()
    run(args.stage,args.list)
