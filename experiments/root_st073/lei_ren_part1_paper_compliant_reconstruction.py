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
Subsequent stages admit the uniform local pulse interfaces; complete
post-pulse/whole-field C4 remains unfinished.
Physicaljets maps actual leading profiles to physical cylindrical r/z
derivatives and finite log-bound ledgers. Interfacejets identifies internal
pulse charts by exact uncapped source equations, not interval overlap.
Flatcomparison supplies quantitative moment/velocity differences against
continuations with the same exact interface histories. Externaljets rebuilds
O3 power high derivatives and the entire original O5 flatten, and admits
both local external joins. Leading O4 pulse C4 is certified in its stated
local spatial scope; whole outer/core/axis C4 and stress remain unfinished.
Postmixedjets adds the entire following power and both actual angular
supports. Correlated future-tail factorization keeps q^2 cancellation and
positive energy before enclosure, with source-identified O5/O6 joins.
Physicalfield maps the accepted local O3/O4 and complete postpulse chain to
Cartesian vector derivatives through4, including the moving cylindrical basis,
and the first fixed-x physical-time derivative. It also restores physical
volume/energy weights and source-bound local postpulse energy domains; the
original unlocalized Gamma field has infinite whole-space kinetic energy.
Core/axis physical assembly, stress, terminal-time energy and recursion remain.
Corephysical adds SAME selected-Cstar whole-Z analytic core jets through5,
nonsingular Cartesian spatial4 and first physical-time derivative enclosures,
exact axis data/slopes, primitive pressure and local physical energy. The
nonlinear fixed point is enclosed separately from its infinite Bessel model;
inner/matching annuli and full physical assembly remain unfinished.
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
FLATCOMPARISON=('pulse_flat_comparison','pulse_flat_comparison_check')
EXTERNALJETS=('power_inlet_C4','power_inlet_C4_check','flatten_mixed_C4','flatten_mixed_C4_check')
POSTMIXEDJETS=('power_angular_C4','power_angular_C4_check')
STEEPJETS=('steep_waiting_C4','steep_waiting_C4_check')
HEATJETS=('collar_Gamma_C4','collar_Gamma_C4_check')
PHYSICALFIELD=('cartesian_field','cartesian_field_check','physical_energy','physical_energy_check')
COREPHYSICAL=('core_physical_field','core_physical_field_check')
FROZENFIELD=('frozen_comparison_field','frozen_comparison_field_check')
BRIDGEPROFILES=('inner_bridge_profiles','inner_bridge_profiles_check')


def stages(stage):
    return {'source':SOURCE,'inner':INNER,'outer':OUTER,'angular':ANGULAR,'energy':ENERGY,'pulse':PULSE,'postpulse':POSTPULSE,'closure':CLOSURE,'angularjets':ANGULARJETS,'axialjets':AXIALJETS,'pulsejets':PULSEJETS,'radialjets':RADIALJETS,'flatjets':FLATJETS,'mixedjets':MIXEDJETS,'physicaljets':PHYSICALJETS,'interfacejets':INTERFACEJETS,'flatcomparison':FLATCOMPARISON,'externaljets':EXTERNALJETS,'postmixedjets':POSTMIXEDJETS,'steepjets':STEEPJETS,'heatjets':HEATJETS,'physicalfield':PHYSICALFIELD,'corephysical':COREPHYSICAL,
            'frozenfield':FROZENFIELD,
            'bridgeprofiles':BRIDGEPROFILES,
            'all':SOURCE+INNER+OUTER+ANGULAR+ENERGY+PULSE+POSTPULSE+CLOSURE+ANGULARJETS+AXIALJETS+PULSEJETS+RADIALJETS+FLATJETS+MIXEDJETS+PHYSICALJETS+INTERFACEJETS+FLATCOMPARISON+EXTERNALJETS+POSTMIXEDJETS+STEEPJETS+HEATJETS+PHYSICALFIELD+COREPHYSICAL+FROZENFIELD+BRIDGEPROFILES}[stage]


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
    inherited_bridge='inner_bridge_profiles_check' in completed
    inherited_frozen=inherited_bridge or 'frozen_comparison_field_check' in completed
    inherited_core=inherited_frozen or 'core_physical_field_check' in completed
    inherited_energy=inherited_core or 'physical_energy_check' in completed
    inherited_cartesian=inherited_energy or 'cartesian_field_check' in completed
    inherited_heat=inherited_cartesian or 'collar_Gamma_C4_check' in completed
    inherited_steep=inherited_heat or 'steep_waiting_C4_check' in completed
    inherited_post=inherited_steep or 'power_angular_C4_check' in completed
    inherited_external=inherited_post or 'flatten_mixed_C4_check' in completed
    inherited_power=inherited_external or 'power_inlet_C4_check' in completed
    inherited_flat=inherited_power or 'pulse_flat_comparison_check' in completed
    inherited_physical=inherited_flat or any(name in completed for name in ('pulse_physical_bounds_check','pulse_interface_certificate'))
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
        exact_functional_main_gap_and_gap_end_identities_certified=inherited_flat or 'pulse_interface_certificate' in completed,
        quantitative_flat_velocity_moment_comparison_ledger_available=inherited_flat,
        two_sided_O3_pulse_join_certified=inherited_power,
        two_sided_pulse_O5_join_certified=inherited_external,
        entire_original_100_unit_flatten_mixed_C4_available=inherited_external,
        entire_following_power_high_mixed_derivatives_available=inherited_post,
        both_actual_angular_supports_high_mixed_derivatives_available=inherited_post,
        flatten_power_and_power_angular_joins_certified=inherited_post,
        entire_steep_waiting_high_mixed_derivatives_available=inherited_steep,
        angular_steep_and_internal_joins_certified=inherited_steep,
        waiting_exit_Gamma_future_source_retained=inherited_steep,
        entire_collar_and_unbounded_Gamma_high_mixed_derivatives_available=inherited_heat,
        waiting_collar_and_collar_Gamma_joins_certified=inherited_heat,
        accepted_outer_chart_cartesian_spatial_C4_mapped=inherited_cartesian,
        accepted_outer_chart_first_physical_time_derivative_mapped=inherited_cartesian,
        physical_volume_and_kinetic_energy_functional_restored=inherited_energy,
        complete_postpulse_local_physical_energy_bounds_available=inherited_energy,
        whole_Z_analytic_core_profile_enclosures_through_order5_available=inherited_core,
        whole_core_and_axis_cartesian_spatial4_enclosures_available=inherited_core,
        whole_core_and_axis_first_physical_time_derivative_enclosures_available=inherited_core,
        core_local_physical_energy_bounds_available=inherited_core,
        actual_core_exit_five_primitive_axial5_enclosures_available=inherited_frozen,
        frozen_comparison_profile_mixed4_available=inherited_frozen,
        frozen_inertial_direction_axial4_enclosures_available=inherited_frozen,
        source_bound_smoothed_comparison_axial6_enclosures_available=inherited_bridge,
        actual_prescribed_shear_bridge_axial5_enclosures_available=inherited_bridge,
        actual_bridge_radial_recovery_axial4_enclosures_available=inherited_bridge,
        actual_bridge_radial_mixed4_certified=False,
        bridge_100_110_switches_installed=False,
        actual_smooth_comparison_installed=False,
        actual_prescribed_shear_bridge_installed=False,
        original_unlocalized_whole_space_kinetic_energy_is_infinite=True if inherited_energy else None,
        full_cartesian_vector_derivatives_certified=False,
        core_axis_interfaces_certified=False,
        physical_energy_integral_certified=False,
        global_physical_energy_integral_certified=False,
        full_pulse_C4_installed=inherited_external,
        full_pulse_C4_scope='Leading O4 and two local external interfaces, spatial/profile mixed total<=4 at fixed positive tau; whole outer/core/axis and time derivatives excluded',
        full_outer_C4_certified=False,
        temporal_recursion=False,
        next_dependency='Actual bridge phase/radial mixed4 with formal inverse-hb scales -> original100..110 switches ->110 inlet/reshape/restoration/moment patch and preceding O3/full physical assembly; energy/stress/flat remainder -> actual recursion/correction')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('source','inner','outer','angular','energy','pulse','postpulse','closure','angularjets','axialjets','pulsejets','radialjets','flatjets','mixedjets','physicaljets','interfacejets','flatcomparison','externaljets','postmixedjets','steepjets','heatjets','physicalfield','corephysical','frozenfield','bridgeprofiles','all'),default='all')
    parser.add_argument('--list',action='store_true',help='Print the ordered modules without running them')
    args=parser.parse_args()
    run(args.stage,args.list)
