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
and temporal recursion remain unfinished. Pulsejets transports all five
partial primitives through4 and radial velocity through3 in every O.4 chart.
Radialjets admits separate fifth-order source data and recovers radial
velocity axial derivatives through4; higher radial/interface bounds remain.
Flatjets supplies original sigma/gp/beta derivatives and flat support
majorants. Mixedjets recovers every velocity/pressure multiindex in (y,Z)
with total order<=4 from the same primitive ODEs and physical prefactors.
Uniform whole-pulse interfaces and post-pulse C4 remain unfinished.
Physicaljets maps actual leading profiles to physical cylindrical r/z
derivatives and finite log-bound ledgers. Interfacejets identifies internal
pulse charts by exact uncapped source equations, not interval overlap.
Quantitative flat velocity and two-sided external joins remain pending.
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
PULSEJETS=('pulse_high_jets','pulse_high_jets_check')
RADIALJETS=('fifth_axial_jets','fifth_axial_jets_check','pulse_radial_C4','pulse_radial_C4_check')
FLATJETS=('flat_pulse_derivatives','flat_pulse_derivatives_check')
MIXEDJETS=('pulse_mixed_C4','pulse_mixed_C4_check')
PHYSICALJETS=('pulse_physical_bounds','pulse_physical_bounds_check')
INTERFACEJETS=('pulse_interface_certificate',)


def stages(stage):
    return {'source':SOURCE,'inner':INNER,'outer':OUTER,'angular':ANGULAR,'energy':ENERGY,'pulse':PULSE,'postpulse':POSTPULSE,'closure':CLOSURE,'angularjets':ANGULARJETS,'axialjets':AXIALJETS,'pulsejets':PULSEJETS,'radialjets':RADIALJETS,'flatjets':FLATJETS,'mixedjets':MIXEDJETS,'physicaljets':PHYSICALJETS,'interfacejets':INTERFACEJETS,
            'all':SOURCE+INNER+OUTER+ANGULAR+ENERGY+PULSE+POSTPULSE+CLOSURE+ANGULARJETS+AXIALJETS+PULSEJETS+RADIALJETS+FLATJETS+MIXEDJETS+PHYSICALJETS+INTERFACEJETS}[stage]


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
    inherited_physical=any(name in completed for name in ('pulse_physical_bounds_check','pulse_interface_certificate'))
    inherited_mixed=inherited_physical or 'pulse_mixed_C4_check' in completed
    inherited_angular=inherited_mixed or any(name in completed for name in ('future_energy_high_jets_check','axial_high_jets_check','pulse_high_jets_check','fifth_axial_jets_check','pulse_radial_C4_check','flat_pulse_derivatives_check'))
    inherited_fifth=inherited_mixed or any(name in completed for name in ('fifth_axial_jets_check','pulse_radial_C4_check','flat_pulse_derivatives_check'))
    inherited_radial=inherited_mixed or any(name in completed for name in ('pulse_radial_C4_check','flat_pulse_derivatives_check'))
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
        actual_selected_ap_c1_c2_C4_available=inherited_fifth or any(name in completed for name in ('axial_high_jets_check','pulse_high_jets_check')),
        actual_selected_ap_c1_c2_C5_available=inherited_fifth,
        pulse_primitives_axial_C4_installed=inherited_radial or 'pulse_high_jets_check' in completed,
        pulse_radial_velocity_axial_C3_installed=inherited_radial or 'pulse_high_jets_check' in completed,
        pulse_primitives_axial_C5_installed=inherited_radial,
        pulse_radial_velocity_axial_C4_installed=inherited_radial,
        original_radial_shape_derivatives_C4_available=inherited_mixed or 'flat_pulse_derivatives_check' in completed,
        pulse_all_mixed_derivatives_total_order_le4_available=inherited_mixed,
        quantitative_flat_support_majorants_available=inherited_mixed or 'flat_pulse_derivatives_check' in completed,
        all_cylindrical_r_z_derivatives_total_order_le4_mapped=inherited_physical,
        complete_pulse_physical_spatial_supremum_log_ledger_available=inherited_physical,
        exact_functional_main_gap_and_gap_end_identities_certified='pulse_interface_certificate' in completed,
        full_pulse_C4_installed=False,
        temporal_recursion=False,
        next_dependency='Quantitative flat velocity/moment bounds and two-sided external joins; compose post-pulse C4, Cartesian map, energy and stress cone; then genuine temporal recursion')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('source','inner','outer','angular','energy','pulse','postpulse','closure','angularjets','axialjets','pulsejets','radialjets','flatjets','mixedjets','physicaljets','interfacejets','all'),default='all')
    parser.add_argument('--list',action='store_true',help='Print the ordered modules without running them')
    args=parser.parse_args()
    run(args.stage,args.list)
