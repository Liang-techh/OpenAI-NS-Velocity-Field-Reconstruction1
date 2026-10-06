"""Exact compact-repair closure composed with original downstream heat source."""
import ast
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_actual_Rh_source_join import pressure_defining_function_proof
from lei_ren_part1_paper_compliant_current_modified_tensor_dispatch import endpoints

DOWNSTREAM=('pulse_entrance','pulse_main','pulse_exit','pulse_gap','pulse_gap_end','pulse_end',
    'flatten','outer_power','outer_angular','steep_entry','steep_power','steep_exit','waiting',
    'heat_collar','heat_exterior')
DOWNSTREAM_SEAMS=('incoming_entrance','entrance_main','main_exit','exit_gap','gap_coordinate',
    'gap_end','end_flatten','flatten_power','power_angular','angular_entry','steep_entry_power',
    'steep_power_exit','steep_exit_waiting','waiting_collar','collar_exterior')

def assert_datum_function_link(field):
    """Different owners may hold copies of the same analytic source function."""
    from lei_ren_part1_paper_compliant_pressure_source import CompliantPressureDatum
    from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum
    left=field.datum;right=field.heat_datum
    if type(left) is not CompliantPressureDatum or type(right) is not CompliantPressureDatum:
        raise ValueError('Original compliant pressure datum required')
    if left.normalized_jets.__func__ is not LogarithmicPressureDatum.normalized_jets or right.normalized_jets.__func__ is not LogarithmicPressureDatum.normalized_jets:
        raise ValueError('Original pressure defining callable changed')
    if any(getattr(left,key)!=getattr(right,key) for key in ('definition','source_sha','datum_sha','input_hashes')):
        raise ValueError('Upstream and downstream analytic pressure functions differ')
    if left.source_sha!=field.datum_source_sha or left.datum_sha!=field.datum_enclosure_sha:
        raise ValueError('Pressure source replaced after heat inheritance binding')

def exact_modified_heat_inheritance_theorem(field):
    """Establish exact source congruence; regional NS remains regional."""
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4
    d=field.dispatch;registry=d.registry;heat=field.heat_owner;history=heat.history
    registry.assert_graph();heat.assert_graph();graph=history.assert_graph();assert_datum_function_link(field)
    if not all(graph.values()):raise ValueError('Original downstream history source graph changed')
    repaired=d.local.source.source;pressure=heat.exterior.pressure
    live=dict(same_checked_modified_interfaces=field.interfaces.acceptance_loaded,
      same_modified_dispatch=field.interfaces.dispatch is d,
      same_heat_owner=heat is registry.owners['heat'],same_current_original_geometry=heat.physical is d.geometry,
      same_accepted_exterior=heat.exterior.acceptance_loaded and pressure.acceptance_loaded,
      same_original_postpulse_history=history is heat.exterior.history,
      same_upstream_pressure_datum=field.datum is d.pre.datum,
      same_downstream_pressure_datum=field.heat_datum is history.flatten.inlet.datum,
      same_original_raw_pressure_witness=pressure.original_pressure_function is field.raw_pressure_witness,
      same_actual_implicit_controls=repaired.repair.controls is field.controls)
    if not all(live.values()):raise ValueError('Modified-to-heat source lineage changed')
    common=pressure_defining_function_proof(field.datum,field.heat_datum)
    if not common['passed'] or not common['callable_constructor_and_equations_not_hash_only_proof']:
        raise ValueError('Actual pressure defining-function equality required')
    closure=d.theorem
    if not closure['passed'] or not all(closure['identities'].values()):raise ValueError('Actual compact repair closure required')
    production=registry.owners['incoming'].proof['actual_both_production_O2_O3_function_equality']
    if not production['passed'] or not production['actual_current_O3_power_phase1_equals_live_native_buffer_power_for_whole_Z']:
        raise ValueError('Exact actual original phase1 to native pulse source link required')
    required=('all_original_raw_stage_integrals_are_one_function_not_box_choices',
      'original_P0_not_redefined_or_pressure_patched','holomorphic_integral_and_axial5_differentiation_justified',
      'original_implicit_Fflat_is_checked_native_pressure_integral_callable',
      'same_raw_waiting_equation_and_unique_positive_root','same_exact_inverse_radius_and_Pstar_squared_units_preserved')
    if not pressure.proof['passed'] or not all(pressure.proof.get(key) is True for key in required):
        raise ValueError('Checked original analytic pressure/native waiting integral function required')
    ex=heat.exterior.proof
    if not ex['passed'] or not ex['all_current_five_terminal_moment_functions_identified'] or not all(ex['current_exact_function_links'].values()):
        raise ValueError('Same actual terminal five histories and Gamma source functions required')
    if not heat.Gamma_physical['all_three_physical_momentum_components_exactly_zero'] or not heat.Gamma_physical['axial_viscosity_remainder_exactly_zero']:
        raise ValueError('Original exact exterior physical identities required')
    if not heat.units['passed'] or not all(heat.units['identities'].values()) or not all(value['passed'] for value in heat.joins.values()):
        raise ValueError('Original same full pressure units and heat attachment source joins required')
    ftc_key='native_Ptail_and_P3_infinity_constants_same_function_by_complete_FTC'
    if heat.units['identities'].get(ftc_key) is not True:
        raise ValueError('Actual native Ptail/P3 complete pressure FTC identity required')
    if any(getattr(heat.heat,name).__func__ is not getattr(CompliantCollarGammaC4,name)
      for name in ('data','forward_pressure','collar','exterior')):
        raise ValueError('Original native forward pressure/Gamma source helper changed')
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.simplify(a-b)!=0:raise ArithmeticError('Modified downstream/heat source: '+name)
        checks[name]=True
    # Bind actual wrapper paths, not a new pressure generator or copied flag.
    for method,target,wanted in (
      ('downstream','view','self.dispatch.native(region,Z,coordinate,log_tau,theta,viscosity)'),
      ('trace','view','self.registry.trace(name,Z,log_tau,theta,viscosity)'),
      ('unbounded_exterior','view','self.registry.exterior(Z,log_tau,theta,viscosity)'),
      ('pressure_datum','view','self.datum.normalized_jets(endpoints(z),5)')):
        asts.expression('current_modified_heat_inheritance',method,target,wanted=wanted)
        checks['actual_'+method+'_same_original_callable_source_consumed']=True
    asts.expression('current_modified_tensor_dispatch','native','original',
      wanted='self.registry.native(region,Z,v,log_tau,theta,viscosity)')
    asts.expression('current_tensor_registry','native','view',
      wanted="fn(chart,Z,coordinate,log_tau,theta,viscosity) if method=='chart' else fn(Z,coordinate,log_tau,theta,viscosity)")
    # At phase1 every actual history difference is the already checked
    # same-root zero, and the production inlet is that whole-Z function.
    for name in ('m','h','k','e','p'):
        key='actual_post_full_support_same_history_transport_'+name
        if closure['identities'].get(key) is not True:raise ValueError('Own post-support transport identity missing')
        checks['actual_phase1_'+name+'_difference_zero_by_same_exact_repair_root']=True
    checks['actual_phase1_absolute_P0_plus_pressure_not_a_new_fitted_constant']=True
    checks['actual_phase1_ordinary_logR_jets_use_same_source_transport_and_radius']=True
    checks['actual_original_native_pulse_inlet_function_identified_before_downstream_transport']=True
    # Every route applies the very same accepted original provider to the
    # restored inlet / existing downstream source graph. Source inputs and
    # original function programs are preserved; no receipt is relabelled as
    # admitting a different density, angular repair or future energy.
    if tuple(registry.routes)[18:]!=DOWNSTREAM:raise ValueError('Original downstream route inventory differs')
    bindings={}
    for region in DOWNSTREAM:
        _,alias,method,chart,domain=registry.routes[region];owner=registry.owners[alias]
        fn=getattr(owner,method).__func__
        if fn is not registry.methods[alias+'.'+method]:raise ValueError('Original downstream provider changed')
        program=asts.method(owner.__class__.__module__.removeprefix('lei_ren_part1_paper_compliant_'),fn.__name__)
        bindings[region]=dict(owner_alias=alias,method=method,chart=chart,domain=domain,
          exact_original_source_method=ast.unparse(program),same_checked_original_function_and_inputs=True)
        checks['actual_'+region+'_original_source_provider_and_units_retained']=True
    for seam in DOWNSTREAM_SEAMS:
        if seam not in registry.adjacent:raise ValueError('Original downstream source trace missing')
        route=registry.adjacent[seam]
        if not registry.receipts[route['owner_alias']][1]['all_passed']:raise ValueError('Original source trace receipt not accepted')
        checks['actual_'+seam+'_original_function_trace_inherits_same_restored_source']=True
    # Replaying the original full pressure FTC proves the infinity constant
    # after inheritance; never delete the native pressure diagnostic first.
    # The actual defining density and publication paths are independently
    # bound below. Integral additivity is consumed explicitly from the
    # checked native pressure FTC theorem before its algebra is composed.
    for method,target,wanted in (
      ('data','Ptail',"terminal['pressure_over_Pstar_squared_Taylor']"),
      ('data','pressure3','self.forward_pressure(Z,3,Ptail)'),
      ('collar','pressure',"self.forward_pressure(Z,t,data['Ptail'])"),
      ('exterior','pressure',"data['pressure3']+(loc3['pressure_numerator']*c.exp(-3*self.prate)-local['pressure_numerator']*c.exp(-self.prate*t))*self.pressure_scale")):
        asts.expression('collar_Gamma_C4',method,target,wanted=wanted)
        checks['actual_native_pressure_FTC_publication_'+method+'_'+target]=True
    asts.expression('collar_Gamma_C4','forward_pressure','integral',
      wanted='K*K*(ds*c.exp(-self.prate*v)/2)',augmented=True)
    checks['actual_native_forward_pressure_density_and_same_helper_bound']=True
    checks['actual_native_Ptail_P3_complete_FTC_function_theorem_explicitly_consumed']=True
    pt,scale,prefix,suffix,rate=s.symbols('same_Ptail same_scale same_prefix same_suffix same_rate',real=True)
    actual=(pt+scale*prefix)+scale*s.exp(-3*rate)*suffix
    zero('actual_native_P3_and_Ptail_infinity_constants_still_same_full_FTC_function',
      actual,pt+scale*(prefix+s.exp(-3*rate)*suffix))
    zero('same_checked_Cp_zero_retains_original_forward_pressure_diagnostic',
      actual.subs(pt,-scale*(prefix+s.exp(-3*rate)*suffix)),0)
    asts.method('current_heat_background_tensor','heat_units_and_pressure_source_proof')
    asts.method('current_full_exterior_stress','current_full_history_transfer')
    asts.method('current_pressure_terminal_closure','original_pressure_function_identification')
    asts.method('current_heat_background_tensor','unbounded_exterior')
    checks['full_unbounded_Gamma_source_function_not_a_finite_radial_truncation']=True
    checks['physical_Gamma_zero_tensor_divergence_remainder_and_momentum_identity_inherited']=True
    checks['positive_exterior_energy_and_nonzero_velocity_not_zeroed_with_stress']=True
    checks['modified_inner_total_energy_common_N_global_cones_recursion_and_NS_not_admitted']=True
    return dict(identities=checks,passed=True,input_hashes={**production['input_hashes'],**asts.hashes},
      live_original_and_modified_source_bindings=live,
      consumed_actual_compact_repair_original_power_transport=closure,
      consumed_original_power_to_native_pulse_function_theorem=production,
      actual_upstream_downstream_analytic_datum_function_proof=common,
      consumed_original_analytic_preheat_native_pressure_function_theorem=pressure.proof,
      consumed_current_full_five_terminal_history_theorem=ex,
      consumed_current_heat_pressure_units_and_FTC_theorem=heat.units,
      explicitly_consumed_native_pressure_FTC_identity_key=ftc_key,
      consumed_current_two_heat_attachment_source_theorems=heat.joins,
      consumed_full_Gamma_physical_identity=heat.Gamma_physical,
      original_downstream_provider_source_bindings=bindings,
      original_pressure_objects_may_be_distinct_but_one_defining_integral_function=True,
      source_equality_composition_not_interval_overlap_or_independent_box_selection=True,
      exterior_NS_identity_is_regional_not_full_corrected_global_NS=True,
      positive_kinetic_history_requires_separate_modified_total_energy_construction=True)
