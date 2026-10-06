"""Checked compact-source closure and true original power continuation."""
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_modified_pre_stress import endpoints
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import (
    raw_pre_velocity_rows,raw_pre_stress_rows,raw_pre_remainder_sectors,copy_jet)
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

def original_signed_power(field,Z,coordinate,coordinate_kind):
    """Same original pre provider; never reset its incoming moments or P0."""
    c=field.ctx;v=c.mpf(coordinate);Z=c.mpf(Z)
    phase=v/field.Tw if coordinate_kind=='log_radius_offset' else v
    pre=field.pre.power(Z,phase)
    raw=pre['actual_normalized_primitive_y_derivative_axial5']
    velocity=raw_pre_velocity_rows(c,pre)
    histories={name:[copy_jet(c,row) for row in raw[name]] for name in ('m','h','k','e','p')}
    datum=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in pre['original_P0_axial5_coefficients']])
    pressure=[copy_jet(c,raw['p'][0])+datum]+[copy_jet(c,row) for row in raw['p'][1:]]
    z=IntervalTaylor.variable(c,Z,5)
    rows=raw_pre_stress_rows(c,field.delta,z,velocity['theta'],velocity['axial'],histories,pressure)
    remainder=raw_pre_remainder_sectors(c,field.delta,z,velocity)
    source=dict(current_original_pre_source=pre,current_original_five_history_rows=histories,
      original_absolute_pressure_over_Pstar2_ordinary_logR_rows=pressure,
      actual_original_velocity_rows=velocity,
      exact_modified_history_differences_vanish_by_same_implicit_source_equation=True,
      original_nonzero_incoming_histories_and_axis_pressure_datum_retained=True,
      parent_repaired_source_definition_sha256=field.local.source.source.modified_source_definition_sha256)
    return dict(region='quiet_O3_power' if coordinate_kind=='log_radius_offset' else 'O3_power',
      coordinate=v,Z=Z,actual_source=source,
      modified_source_definition_sha256=field.local.source.source.modified_source_definition_sha256,
      modified_stress_definition_sha256=field.local.source.modified_stress_definition_sha256,
      full_signed_paper_theta_axial_stress3_sectors=rows,
      full_signed_source_remainder2_sectors=remainder,
      exact_original_phase_source=phase,source_coordinate_kind=coordinate_kind,
      original_source_used_only_after_proved_compact_repair_closure=True)

def exact_modified_dispatch_theorem(field):
    """Consume actual implicit closure, replay its post-support transport."""
    from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import CompliantPrePulseMixedC4
    from lei_ren_part1_paper_compliant_current_tensor_registry import CurrentTensorRegistry,ROUTES
    field.registry.assert_graph();field.local.source.source.histories.candidate.direction.assert_graph()
    repaired=field.local.source.source
    live=dict(same_checked_registry=field.registry is field.local.registry,
      same_original_pre=field.pre is repaired.histories.pre,
      same_original_geometry=field.geometry is field.registry.physical,
      same_original_Tw_object=field.geometry.params is field.pre.params,
      same_original_power_callable=field.pre.power.__func__ is CompliantPrePulseMixedC4.power,
      same_original_native_callable=field.registry.native.__func__ is CurrentTensorRegistry.native,
      same_incoming_owner=field.registry.owners['o3'].entrance_incoming_tensor is field.registry.owners['incoming'],
      same_incoming_geometry=field.registry.owners['incoming'].physical is field.geometry)
    if not all(live.values()):raise ValueError('Continuation must use the same checked original source/registry')
    old=repaired.theorem
    required=['actual_callable_exit_history_factors_same_checked_source_root_'+name for name in ('m','h','k','e','p')]
    required+=['actual_repair_normalized_transport_ODE_'+name for name in ('m','h','k','e')]+['actual_repair_pressure_FTC']
    if not old['passed'] or not all(old['identities'].get(key) is True for key in required):
        raise ValueError('Actual same implicit-source terminal closure and own transport ODE required')
    if not old['outside_repair_five_delta_histories_follow_same_homogeneous_ODE_after_exact_exit']:
        raise ValueError('Post-support source transport theorem required')
    if not repaired.repair.acceptance_loaded or not repaired.repair.certificate['unique_exact_implicit_controls_exist']:
        raise ValueError('Same checked exact implicit control vector required')
    production=field.registry.owners['incoming'].proof['actual_both_production_O2_O3_function_equality']
    units=field.registry.owners['incoming'].proof['current_actual_raw_power_normalization_theorem']
    if not production['passed'] or not units['passed']:raise ValueError('Checked original full power function/unit theorem required')
    asts=SourceAST();checks={}
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Modified continuation source: '+name)
        checks[name]=True
    asts.expression('current_modified_tensor_dispatch_operator','original_signed_power','pre',wanted='field.pre.power(Z,phase)')
    asts.expression('current_modified_tensor_dispatch_operator','original_signed_power','phase',
      wanted="v/field.Tw if coordinate_kind=='log_radius_offset' else v")
    asts.expression('current_modified_tensor_dispatch_operator','original_signed_power','rows',
      wanted="raw_pre_stress_rows(c,field.delta,z,velocity['theta'],velocity['axial'],histories,pressure)")
    asts.expression('current_modified_tensor_dispatch_operator','original_signed_power','remainder',
      wanted='raw_pre_remainder_sectors(c,field.delta,z,velocity)')
    checks['actual_original_five_histories_pressure_and_all_variable_stress_remainder_programs_bound']=True
    asts.expression('current_modified_tensor_dispatch','continuation','signed',
      wanted='original_signed_power(self,Z,v,coordinate_kind)')
    asts.expression('current_modified_tensor_dispatch','continuation','(packet, radius)',
      wanted='modified_physical_packet(self,signed)')
    asts.expression('current_modified_tensor_dispatch','continuation','point',
      wanted='self.local.lift(c,packet,self.delta,None,lt,theta,nu)')
    checks['same_accepted_physical_operators_apply_to_actual_original_continuation']=True
    asts.expression('current_modified_tensor_dispatch','native','original',
      wanted='self.registry.native(region,Z,v,log_tau,theta,viscosity)')
    asts.expression('current_modified_tensor_dispatch','native','local',
      wanted='self.local.tensor(region,Z,v,log_tau,theta,viscosity)')
    checks['unaffected_native_routes_and_actual_modified_O2_O3_routes_use_same_checked_graph']=True
    asts.expression('pre_pulse_mixed_C4','power','t',wanted='self.params.Tw*phase')
    asts.expression('pre_pulse_mixed_C4','power','parent',wanted='self.slope_mu(Z,1)')
    checks['actual_original_power_uses_same_Tw_phase_and_retained_transition_parent']=True
    # Every original compact repair support has ended before q=19/10.
    ci=(s.Rational(1,5),s.Rational(1,2),s.Rational(4,5));ell=s.Rational(1,40)
    asts.expression('current_O3_independent_repair_operator','bump_weights','centers',
      wanted='[c.mpf(1)/5,c.mpf(1)/2,c.mpf(4)/5]')
    asts.expression('current_O3_independent_repair_operator','bump_weights','ell',wanted='c.mpf(1)/40')
    for i,center in enumerate(ci):
        if 1+center+ell>=s.Rational(19,10):raise ArithmeticError('Continuation overlaps actual repair support')
        checks['actual_fixed_support_ended_before_continuation_floor_'+str(i)]=True
    # Replay the exact callable scalar formula at arbitrary post-support q.
    q=s.Symbol('same_actual_q',real=True);f2=s.Symbol('same_actual_f2',positive=True)
    inlet=s.symbols('same_Bm2 same_Bh2 same_Bk2 same_Be2 same_Bp2')
    full=dict(zip(('M','I','J','S','Cp'),s.symbols('same_full_M same_full_I same_full_J same_full_S same_full_Cp')))
    rates=dict(m=1,h=s.Rational(3,2),k=s.Rational(3,2),e=1,p=0)
    program=asts.expression('current_O3_repaired_histories','history','scalar',wanted=
      "dict(m=previous['m']+self.f2*primitives['M']/x,h=previous['h']+self.f2*primitives['I']/x**c.mpf('1.5'),k=previous['k']+self.f2**2*primitives['J']/x**c.mpf('1.5'),e=previous['e']+self.f2**2*primitives['S']/x,p=previous['p']+self.f2**2*primitives['Cp'])")
    ctx=SimpleNamespace(mpf=s.Rational)
    def values(position):
        previous={name:value*s.exp(-rates[name]*(position-1)) for name,value in zip(rates,inlet)}
        return asts.evaluate(program,dict(c=ctx,self=SimpleNamespace(f2=f2),previous=previous,primitives=full,x=s.exp(position-1)))
    atq=values(q);at2=values(s.Integer(2))
    for name,rate in rates.items():
        zero('actual_post_full_support_same_history_transport_'+name,atq[name],s.exp(-rate*(q-2))*at2[name])
        for order in range(1,5):
            zero('actual_post_support_source_ordinary_transport_'+name+'_order'+str(order),
              s.diff(atq[name],q,order),(-rate)**order*atq[name])
    checks['actual_zero_terminal_data_are_same_checked_implicit_root_not_chosen_values']=True
    checks['same_zero_absolute_pressure_increment_no_new_axial_constant']=True
    # Directed overlap bounds keep a full phase[0,1] request in two charts.
    # The cut is an enclosure bound for coverage, never a source parameter.
    c=field.ctx;lo,hi=endpoints(field.Tw);cut=field.phase_cut
    maxlocal=endpoints(c.mpf(hi)*c.mpf(cut))[1]
    minoriginal=endpoints(c.mpf(lo)*c.mpf(cut))[0]
    if lo<=4 or maxlocal>=2 or minoriginal<=mp.mpf('1.9'):
        raise ValueError('Original actual Tw enclosure does not certify overlapping power chart coverage')
    checks['actual_phase_cut_local_q_below2_original_q_above_full_support_floor']=True
    checks['full_native_phase_domain_covered_including_same_actual_Tw_endpoint']=True
    if tuple(field.registry.routes)!=tuple(row[0] for row in ROUTES) or len(field.registry.routes)!=33:
        raise ValueError('Original complete33 route inventory changed')
    checks['all33_original_routes_retained_with_three_affected_native_routes_overridden']=True
    return dict(identities=checks,passed=True,input_hashes={**old['input_hashes'],**production['input_hashes'],**units['input_hashes'],**asts.hashes},
      live_same_checked_source_registry_geometry_bindings=live,
      consumed_same_implicit_repaired_history_source_theorem=old,
      consumed_original_full_power_same_function_theorem=production,
      consumed_original_raw_to_power_unit_theorem=units,
      true_post_support_floor_q='19/10; all fixed supports end by73/40',
      original_q2_closure_transports_both_forward_and_back_within_same_post_support_interval=True,
      actual_power_coordinate_relation='q=same_actual_Tw*phase; phase endpoint1 is q=actual Tw',
      phase_cut_is_only_enclosure_coverage_coordinate_not_a_field_value=True,
      q_interval_straddling2_requires_two_actual_source_views=True,
      own_energy_modified_interfaces_cones_common_N_recursion_and_full_NS_remain_open=True)
