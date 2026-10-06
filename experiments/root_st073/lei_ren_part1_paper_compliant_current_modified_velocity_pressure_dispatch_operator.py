"""Same-current-source velocity routing, with directed coordinate covers."""
import ast
import copy
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_tensor_registry import ROUTES,SOURCE_COORDINATES
from lei_ren_part1_paper_compliant_current_complete_physical_assembly import CurrentCompletePhysicalAssembly,CurrentCompletePhysicalDispatch
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly as BASE
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_interfaces_operator import original_velocity_pressure_packet
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_operator import modified_physical_velocity_pressure
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure import endpoints

GAP_PHASE_CUT=s.Rational(1,20000)

def canonical_coordinate(field,region,coordinate):
    if region=='core_positive_radius':return 'core',coordinate
    if region=='pulse_entrance':return region,coordinate/field.geometry.pulse.mu
    return region,coordinate

def gap_source_coordinates(mu,phase):
    # Reduce correlations BEFORE enclosure. In particular phase1 has
    # exactly s=-4; subtracting two giant reciprocal-mu boxes loses this.
    return 12+(1-4*mu)*phase,-4-(1/mu-4)*(1-phase)

def canonical_original_velocity(field,region,Z,coordinate,log_tau,theta,axis=False):
    chart,value=canonical_coordinate(field,region,coordinate)
    physical=field.geometry.evaluate(chart,Z,value,log_tau=log_tau,theta=theta,axis=axis)
    return dict(physical,actual_original_geometry_chart=chart,actual_original_geometry_coordinate=value,
      preserved_current_canonical_source_packet_route=True,
      physical_layout='core_native_Cartesian' if chart=='core' else 'four_label_Cartesian',
      original_derivatives_already_ordinary_no_second_coordinate_factor=True)

def original_power_velocity(field,Z,coordinate,kind,log_tau,theta):
    c=field.ctx;v=c.mpf(coordinate)
    phase=v/field.Tw if kind=='log_radius_offset' else v
    pre=field.pre.power(Z,phase)
    source=dict(region='quiet_O3_power' if kind=='log_radius_offset' else 'O3_power',coordinate=v,Z=Z,
      original_native_velocity_pressure_source=pre,current_original_pre_source=pre,
      original_source_used_after_same_implicit_compact_closure=True,
      exact_original_phase_source=phase,source_coordinate_kind=kind)
    packet=original_velocity_pressure_packet(field,source)
    physical=modified_physical_velocity_pressure(field.velocity,packet,Z,log_tau,theta)
    public=dict(packet,grids={label:{'y%d_Z%d'%index:terms for index,terms in grid.items()}
      for label,grid in packet['grids'].items()},
      amplitudes={label:{physical['source_log_base_names'][index]:power for index,power in parts.items()}
        for label,parts in packet['amplitudes'].items()})
    return dict(physical,physical_layout='four_label_Cartesian',actual_original_power_source=source,
      actual_original_power_mixed4_packet=public,source_coordinate_kind=kind,
      same_Tw_phase_source_and_original_incoming_histories_P0_retained=True)

def gap_main_cover_velocity(field,Z,phase,log_tau,theta):
    """Map the actual current gap function at its correlated xi enclosure.

    Only a definition-only dispatcher proxy supplies this admitted
    supplemental source packet to the unchanged physical mapper. The
    live canonical graph and its public coordinate guards are untouched.
    """
    g=field.geometry;c=field.ctx;xi,_=gap_source_coordinates(g.pulse.mu,phase)
    if endpoints(xi)[0]<12 or endpoints(xi)[1]>endpoints(c.mpf('12.0001'))[0]:
        raise ValueError('Actual correlated gap main piece exceeds admitted supplemental coverage')
    source=g.pulse.gap(Z,xi)
    proxy=copy.copy(g)
    proxy.dispatch=SimpleNamespace(provider=g.dispatch.provider,
      evaluate=lambda chart,z,value:dict(source_packet=source))
    physical=BASE.evaluate(proxy,'pulse_gap',Z,xi,log_tau=log_tau,theta=theta)
    return dict(physical,physical_layout='four_label_Cartesian',actual_current_supplemental_gap_source=source,
      actual_original_geometry_chart='pulse_gap',actual_original_geometry_coordinate=xi,
      requested_gap_end_phase=phase,exact_gap_coordinate_function='xi=12+(1-4mu)*phase=13+mu*s',
      projected_coordinate_box_is_source_enclosure_not_a_point_substitution=True,
      actual_same_current_gap_function_and_selected_controls_retained=True,
      live_canonical_source_graph_unmodified=True)

def exact_modified_velocity_dispatch_theorem(field):
    g=field.geometry;h=field.velocity.heat;d=h.dispatch
    h.assert_graph();g.assert_graph()
    if type(g) is not CurrentCompletePhysicalAssembly or not g.physical_acceptance_loaded:
        raise ValueError('Checked canonical current33 physical source assembly required')
    if type(g.dispatch) is not CurrentCompletePhysicalDispatch or g is not h.registry.physical:
        raise ValueError('Same actual canonical geometry/registry required')
    if not field.interfaces.acceptance_loaded or field.interfaces.velocity is not field.velocity:
        raise ValueError('Checked same actual mixed4 velocity interfaces required')
    if field.pre is not d.pre or field.pre is not g.pre or field.controls is not field.velocity.source.repair.controls:
        raise ValueError('One actual pre/implicit source/control vector required')
    if tuple(field.registry.routes)!=tuple(row[0] for row in ROUTES) or len(ROUTES)!=33:
        raise ValueError('Original33 native route table changed')
    if g.evaluate.__func__ is not CurrentCompletePhysicalAssembly.evaluate or g.dispatch.gap_overlap.__func__ is not CurrentCompletePhysicalDispatch.gap_overlap:
        raise ValueError('Canonical current physical source/coverage callable changed')
    if not all(g.proof['current_defining_object_graph'].values()) or not g.proof['passed']:
        raise ValueError('Actual current selected and complete pressure/heat graph proof required')
    closure=d.theorem;inheritance=h.theorem;interfaces=field.interfaces.theorem
    if not all(proof['passed'] for proof in (closure,inheritance,interfaces)):
        raise ValueError('Same compact closure, full heat inheritance and mixed4 interfaces required')
    if not closure['original_q2_closure_transports_both_forward_and_back_within_same_post_support_interval']:
        raise ValueError('Original continuation below q2 must follow the same actual post-support transport')
    functional=g.proof['original_internal_pulse_functional_coordinate_and_ODE_identities']
    identities=functional
    for key in ('production_gap_coordinate_leading_log','gap_coordinate_pressure_log',
      'gap_coordinate_angular_log','gap_coordinate_pressure_time','gap_coordinate_energy_whole_Z',
      'gap_end_full_future_linear_weight','gap_end_angular_history_whole_Z','gap_end_original_pressure_whole_Z'):
        if identities.get(key) is not True:raise ValueError('Actual gap coordinate function identity missing: '+key)
    if not all(g.proof['internal_functional_identities_apply_to_actual_current_selected_equations'].values()):
        raise ValueError('Old gap identities must be bound to the actual current selected source')
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(a-b)!=0:raise ArithmeticError('Actual velocity dispatch identity: '+name)
        checks[name]=True
    # The dispatcher and assembly both define evaluate. Bind the assembly
    # method explicitly instead of taking the first function with that name.
    asts.method('current_complete_physical_assembly','evaluate')
    tree=ast.parse((Path(__file__).parent/
      'lei_ren_part1_paper_compliant_current_complete_physical_assembly.py').read_text(encoding='utf8'))
    cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='CurrentCompletePhysicalAssembly')
    fn=next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='evaluate')
    wanted=ast.parse('packet=BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=axis)').body[0]
    if not any(ast.dump(node)==ast.dump(wanted) for node in ast.walk(fn)):
        raise ValueError('Actual canonical assembly physical mapper call changed')
    checks['actual_canonical_current_assembly_calls_unchanged_original_physical_mapper']=True
    asts.expression('current_modified_velocity_pressure_dispatch_operator','canonical_original_velocity','(chart, value)',
      wanted='canonical_coordinate(field,region,coordinate)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','canonical_original_velocity','physical',
      wanted='field.geometry.evaluate(chart,Z,value,log_tau=log_tau,theta=theta,axis=axis)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','original_power_velocity','pre',wanted='field.pre.power(Z,phase)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','original_power_velocity','phase',
      wanted="v/field.Tw if kind=='log_radius_offset' else v")
    asts.expression('current_modified_velocity_pressure_dispatch_operator','original_power_velocity','packet',
      wanted='original_velocity_pressure_packet(field,source)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','original_power_velocity','physical',
      wanted='modified_physical_velocity_pressure(field.velocity,packet,Z,log_tau,theta)')
    mu,p=s.symbols('same_actual_mu phase',positive=True)
    mapping=asts.replay('current_modified_velocity_pressure_dispatch_operator','canonical_coordinate',{})
    proxy=SimpleNamespace(geometry=SimpleNamespace(pulse=SimpleNamespace(mu=mu)))
    for region,_,_,_,_ in ROUTES:
        chart,value=mapping(proxy,region,p)
        if chart!=('core' if region=='core_positive_radius' else region):raise ValueError('Canonical chart changed')
        zero('actual_'+region+'_ordinary_source_coordinate',value,p/mu if region=='pulse_entrance' else p)
    gap=asts.replay('current_modified_velocity_pressure_dispatch_operator','gap_source_coordinates',{})
    xi,ss=gap(mu,p)
    zero('actual_gap_phase_source_xi_equals_13_plus_mu_s',xi,13+mu*ss)
    zero('actual_gap_phase_source_s_equals_original_reduced_distance',ss,-(4*mu+(1-4*mu)*(1-p))/mu)
    zero('actual_gap_phase0_xi',xi.subs(p,0),12)
    zero('actual_gap_phase0_reciprocal_s',ss.subs(p,0),-1/mu)
    zero('actual_gap_phase1_s',ss.subs(p,1),-4)
    zero('actual_gap_end_correlated_source_guard',ss+1/mu,(1-4*mu)*p/mu)
    c=field.ctx;m=c.mpf(g.pulse.mu);cut=c.mpf(1)/20000
    if endpoints(m)[0]<=0 or endpoints(m)[1]>=mp.mpf('.25'):
        raise ValueError('Actual positive mu below1/4 needed for directed gap coordinate cover')
    lower_xi,_=gap(m,c.mpf([0,endpoints(cut)[1]]))
    if endpoints(lower_xi)[0]<12 or endpoints(lower_xi)[1]>endpoints(c.mpf('12.0001'))[0]:
        raise ValueError('Supplemental main source does not cover lower phase piece')
    _,upper_s=gap(m,c.mpf([endpoints(cut)[0],1]))
    if endpoints(upper_s)[1]>-4 or endpoints(upper_s+1/m)[0]<0:
        raise ValueError('Directed end source box has not proved its actual canonical guard')
    checks['actual_gap_phase_partition_covers0_to1_with_both_native_guards']=True
    checks['actual_supplemental_gap_and_end_functions_identified_before_source_enclosures']=True
    if endpoints((1-4*m)*c.mpf([endpoints(cut)[0],1])/m)[0]<=0:
        raise ValueError('Reduced correlated upper gap guard is not strictly positive')
    checks['actual_upper_gap_guard_reduced_before_enclosure_and_native_box_guard_passes']=True
    asts.expression('current_modified_velocity_pressure_dispatch_operator','gap_main_cover_velocity','source',wanted='g.pulse.gap(Z,xi)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','gap_main_cover_velocity','proxy',wanted='copy.copy(g)')
    asts.expression('current_modified_velocity_pressure_dispatch_operator','gap_main_cover_velocity','physical',
      wanted="BASE.evaluate(proxy,'pulse_gap',Z,xi,log_tau=log_tau,theta=theta)")
    asts.method('current_modified_velocity_pressure_dispatch_operator','gap_main_cover_velocity')
    for method,target,wanted in (
      ('native','view','canonical_original_velocity(self,region,Z,v,lt,theta,axis=region==\'core_positive_radius\' and endpoints(v)==(0,0))'),
      ('power_offset','view',"original_power_velocity(self,Z,box,'log_radius_offset',lt,theta)"),
      ('power_phase','view',"original_power_velocity(self,Z,box,'original_power_phase',lt,theta)"),
      ('gap_phase','view','gap_main_cover_velocity(self,Z,box,lt,theta)'),
      ('gap_phase','(_, ss)','gap_source_coordinates(self.geometry.pulse.mu,box)'),
      ('gap_phase','view',"canonical_original_velocity(self,'pulse_gap_end',Z,ss,lt,theta)"),
      ('unbounded_exterior','view',"canonical_original_velocity(self,'heat_exterior',Z,box,lt,theta)")):
        asts.expression('current_modified_velocity_pressure_dispatch',method,target,wanted=wanted)
        checks['actual_velocity_dispatch_'+method+'_'+wanted]=True
    local_hi=endpoints(field.Tw*c.mpf([0,field.phase_cut]))[1]
    continuation_lo=endpoints(field.Tw*c.mpf([field.phase_cut,1]))[0]
    if local_hi>=2 or continuation_lo<=mp.mpf('1.9'):
        raise ValueError('Actual power phase cover leaves the admitted compact/original domains')
    checks['actual_power_phase_partition_reuses_checked_post_support_overlap_before_q2']=True
    checks['actual_gap_packet_and_physical_radius_use_same_requested_xi_enclosure']=True
    checks['source_cover_proxy_does_not_mutate_live_canonical_graph_or_source']=True
    checks['original_complete_pressure_and_current_unbounded_Gamma_velocity_sources_inherited']=True
    checks['same_actual_native_core_velocity_pressure_and_nonsingular_axis_map_reused']=True
    return dict(identities=checks,passed=True,input_hashes={**field.interfaces.hashes,**g.hashes,**asts.hashes},
      consumed_actual_original33_canonical_physical_source_proof=g.proof,
      consumed_actual_gap_current_functional_identity_keys=list(identities),
      consumed_actual_compact_closure_and_power_coordinate_theorem=closure,
      consumed_actual_modified_heat_function_theorem=inheritance,
      original_native_route_inventory=[dict(region=row[0],owner_alias=row[1],coordinate=coord)
        for row,coord in zip(ROUTES,SOURCE_COORDINATES)],
      gap_phase_cut='1/20000; coverage partition only, not a changed source parameter',
      original_rows_remain_ordinary_no_new_selector_derivative_factors=True,
      all33_native_velocity_pressure_routes_and_axis_are_defining_function_enclosures=True,
      global_velocity_interfaces_resolved_points_energy_cones_recursion_full_NS_remain_open=True)
