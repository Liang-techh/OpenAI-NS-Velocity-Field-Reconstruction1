"""Ordered reproduction of the compliant pressure/core/outer/heat family.

Stops on the first failed prerequisite. Each module writes its separate,
source-bound receipt; no legacy .01 file is edited. The angular stage solves
the functional angular/pressure correction. The energy stage supplies the
complete future swirl energy and selects actual ap/c1/c2 functions. The pulse
stage installs the full O.4 interval in segmented coordinates. The postpulse
stage composes corrected outer fields and primitives through Gamma heat.
The closure stage composes absolute pressure/angular source identities and
admits all five terminal moments. Angularjets adds the actual angular branch's
axial derivatives through4. Axialjets supplies complete future energy, exact
incoming functions and actual selected ap/c1/c2 through4. Full outer C4/cone
and temporal recursion remain unfinished.
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
ANGULARJETS=('angular_high_jets','angular_high_jets_check')
AXIALJETS=('future_energy_high_jets','future_energy_high_jets_check',
           'axial_high_jets','axial_high_jets_check')


def stages(stage):
    return {'source':SOURCE,'inner':INNER,'outer':OUTER,'angular':ANGULAR,'energy':ENERGY,'pulse':PULSE,'postpulse':POSTPULSE,'closure':CLOSURE,'angularjets':ANGULARJETS,'axialjets':AXIALJETS,
            'all':SOURCE+INNER+OUTER+ANGULAR+ENERGY+PULSE+POSTPULSE+CLOSURE+ANGULARJETS+AXIALJETS}[stage]


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
    inherited_angular=any(name in completed for name in ('future_energy_high_jets_check','axial_high_jets_check'))
    result=dict(stages_completed_this_run=completed,
        source_relation='epsilon_source=.001*delta; epsilon_core=1/Lambda',
        legacy_source_preserved=True,full_NS_background_completed=False,
        actual_angular_pressure_corrections_completed=inherited_angular or any(name in completed for name in ('outer_angular_repair_check','axial_amplitude_selection_check','axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check','angular_high_jets_check')),
        actual_ap_selected=inherited_angular or any(name in completed for name in ('axial_amplitude_selection_check','axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check','angular_high_jets_check')),
        actual_O4_partial_pulse_field_installed=inherited_angular or any(name in completed for name in ('axial_pulse_field_check','corrected_outer_field_check','absolute_moment_closure_check','angular_high_jets_check')),
        corrected_post_Rv_fields_and_primitives_installed=inherited_angular or any(name in completed for name in ('corrected_outer_field_check','absolute_moment_closure_check','angular_high_jets_check')),
        all_five_absolute_terminal_moments_admitted=inherited_angular or any(name in completed for name in ('absolute_moment_closure_check','angular_high_jets_check')),
        actual_angular_coefficient_C4_available=inherited_angular or 'angular_high_jets_check' in completed,
        complete_future_corrected_energy_C4_available=inherited_angular,
        actual_selected_ap_c1_c2_C4_available='axial_high_jets_check' in completed,
        temporal_recursion=False,
        next_dependency='Recover higher axial derivatives and certify full interface/C4/cone bounds; build admissible stress and flat remainder, then genuine n-dependent temporal recursion')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('source','inner','outer','angular','energy','pulse','postpulse','closure','angularjets','axialjets','all'),default='all')
    parser.add_argument('--list',action='store_true',help='Print the ordered modules without running them')
    args=parser.parse_args()
    run(args.stage,args.list)
