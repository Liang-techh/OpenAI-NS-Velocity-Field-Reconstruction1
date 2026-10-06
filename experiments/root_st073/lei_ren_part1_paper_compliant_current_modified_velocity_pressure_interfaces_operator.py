"""Actual original/modified germs and complete physical mixed4 traces."""
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces import SEAMS
from lei_ren_part1_paper_compliant_current_modified_tensor_interfaces_operator import REPAIR_EDGES
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure_operator import (
    modified_velocity_pressure_packet,modified_physical_velocity_pressure,UT,UZ,UR,P,zero_powers)
from lei_ren_part1_paper_compliant_current_modified_velocity_pressure import endpoints,PREFIX
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import exact_modified_radius

ORIGINAL_LABELS={UT:'Utheta_over_Pstar',UZ:'Uz',UR:'Ur_over_current_sqrt_R_over_2',P:'P_over_Pstar2'}
ORIGINAL_POWERS={UT:1,UZ:0,UR:0,P:2}

def original_velocity_pressure_source(field,region,Z,coordinate):
    """The actual native pre function, preserving its original unit sectors."""
    c=field.ctx;pre_owner=field.velocity.heat.dispatch.pre;q=c.mpf(coordinate)
    if region=='O2_axial':pre=pre_owner.axial(Z,q)
    elif region=='O2_buffer':pre=pre_owner.axial(Z,buffer_offset=q)
    elif region=='quiet_O3_power':pre=pre_owner.power(Z,q/field.velocity.heat.dispatch.Tw)
    else:raise ValueError('One of three actual original boundary providers required')
    return dict(region=region,coordinate=q,Z=c.mpf(Z),current_original_pre_source=pre,
      original_native_velocity_pressure_source=pre,
      original_pre_source_used_with_same_incoming_moments_and_P0=True)

def original_velocity_pressure_packet(field,view):
    """Original packet rows are already ordinary mixed derivatives, not jets."""
    c=field.ctx;pre=view['original_native_velocity_pressure_source']
    raw=pre['physical_velocity_pressure_y_Z_mixed4'];grids={}
    expected={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    for label,name in ORIGINAL_LABELS.items():
        if set(raw[name])!=expected:raise ValueError('Original full mixed4 velocity/pressure rows required')
        powers=list(zero_powers());powers[1]=mp.mpf(ORIGINAL_POWERS[label])
        grids[label]={(j,n):[(tuple(powers),raw[name]['y%d_Z%d'%(j,n)])]
          for j in range(5) for n in range(5-j)}
    if view['region']=='quiet_O3_power':
        logR,radius=exact_modified_radius(field.velocity,view['region'],view['coordinate'],view)
    else:
        logR,radius=field.velocity.geometry.radius(view['region'],view['coordinate'],pre,field.velocity.geometry.pre)
    logs=(c.mpf(0),field.velocity.geometry.logP,c.mpf(0),c.mpf(0),c.mpf(0),logR,c.ln(2))
    return dict(grids=grids,source_log_bases=logs,amplitudes={UT:{},UZ:{},UR:{5:mp.mpf('.5'),6:mp.mpf('-.5')},P:{}},
      signed_Pstar_sector_inventory={label:[power] for label,power in ORIGINAL_POWERS.items()},
      exact_logR=logR,exact_original_geometry_radius_source=radius,
      original_rows_are_already_ordinary_no_second_Taylor_factorial=True,
      source_orders=dict(ordinary_logR=4,total_mixed=4),original_P0_and_incoming_histories_retained=True)

def boundary_source(field,name,Z,side):
    """Resolve both actual named germs; source charts are never added together."""
    if name not in SEAMS or side not in ('left','right'):raise ValueError('One of twelve named two-sided interfaces required')
    original=False
    if name=='buffer_inlet':
        region,q=('O2_axial',1) if side=='left' else ('O2_buffer',0);original=side=='left'
    elif name=='buffer_transition':region,q=('O2_buffer',11) if side=='left' else ('O3_slope_mu',0)
    elif name=='transition_power':region,q=('O3_slope_mu',1) if side=='left' else ('quiet_O3_power',0)
    elif name=='buffer_modulation_start':region,q='O2_buffer',9;original=side=='left'
    elif name in REPAIR_EDGES or name=='transition_modulation_end':
        return 'modified_named_germ',field.interfaces.endpoint_source(name,Z,side)
    else:region,q='quiet_O3_power',2;original=side=='right'
    if original:return 'original_native_germ',original_velocity_pressure_source(field,region,Z,q)
    return 'modified_native_germ',field.velocity.source.history(region,Z,q)

def physical_boundary_germ(field,name,Z,side,log_tau,theta):
    kind,source=boundary_source(field,name,Z,side)
    packet=(original_velocity_pressure_packet(field,source) if kind=='original_native_germ'
      else modified_velocity_pressure_packet(field.velocity,source))
    physical=modified_physical_velocity_pressure(field.velocity,packet,source['Z'],log_tau,theta)
    public=dict(packet,grids={label:{'y%d_Z%d'%index:terms for index,terms in grid.items()}
      for label,grid in packet['grids'].items()},
      amplitudes={label:{physical['source_log_base_names'][index]:power for index,power in parts.items()}
        for label,parts in packet['amplitudes'].items()})
    return dict(physical,actual_boundary_source_kind=kind,actual_boundary_source=source,
      actual_boundary_mixed4_packet=public,named_source_edge=name,named_source_side=side)

def physical_contribution_groups(view):
    rows={}
    for index,components in view['physical_spatial_cartesian_mixed4'].items():
        for component,parts in components.items():
            for label,row in parts.items():rows['spatial/'+index+'/'+component+'/'+label]=row
    for component,parts in view['first_fixed_x_physical_time_derivative'].items():
        for label,row in parts.items():rows['time/'+component+'/'+label]=row
    if len(rows)!=216:raise ValueError('Complete spatial4/time1 physical trace groups required')
    return rows

def coupled_mixed4_packet_source_joins(field,asts):
    """Replay actual recovery/publication/packet programs at each source germ.

    Common germ inputs come from the consumed actual profile, moment and
    partial-FTC joins. This compares defining functions, never interval
    outputs. The original/own branches keep their distinct signed units.
    """
    import ast
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows as full_product
    from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
    from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import SHAPES
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda v:v[0] if isinstance(v,(tuple,list)) else s.Rational(str(v)) if isinstance(v,(str,int,float)) else v,ln=s.log)
    jet=lambda v,order=5:Projection.function(c,v,order)
    delta,Ps,amp,base,rr,lp=s.symbols('delta Pstar original_amp same_Ua same_logR same_logPstar',real=True)
    C=1/(1+Z**2);z=jet(Z);zero=jet(0)
    # In particular the modulation end t=1/2 lies inside the variable
    # original O3 sigmoid. Its higher log-amplitude derivatives do not
    # vanish. Keep every ordinary original derivative independent.
    logjets=[jet(s.Symbol('same_original_logU_y'+str(j+1),real=True)) for j in range(4)]
    p0=jet(s.Function('same_actual_analytic_P0')(Z))
    seeds={key:jet(s.Function('same_original_'+key)(Z)) for key in SHAPES}
    native_env=dict(math=math,IntervalTaylor=Projection,square=lambda v:v*v,
      derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1))
    for stem,name in (('pre_pulse_mixed_C4','rate_rows'),('pre_pulse_mixed_C4','product_rows'),
      ('long_reshape_mixed_C4','exponential_derivatives')):asts.replay(stem,name,native_env)
    native=asts.replay('pre_pulse_mixed_C4','physical_mixed',native_env)(c,Z,delta,jet(amp*C),logjets,[zero]*5,seeds,p0,1/Ps**2)
    pre=dict(native,Utheta_over_Pstar_axial5_coefficients=jet(amp*C).coefficients,
      log_Utheta_ordinary_y_derivatives=logjets,Uz_ordinary_y_derivative_axial5=[zero]*5)
    raw_env=dict(native_env,copy_jet=lambda ctx,v:v,endpoints=lambda v:(v,v),shifted_rows=shifted_rows)
    old=asts.replay('current_pre_pulse_stress_operator','raw_pre_velocity_rows',raw_env)(c,pre)
    increment=asts.replay('current_O3_modulated_histories_operator','increment_rows',
      dict(SHAPES=SHAPES,square=lambda v:v*v,product_rows=full_product,
        derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1)))
    publications={}
    for stem in ('current_O3_modulated_histories','current_O3_repaired_histories'):
        node=asts.method(stem,'history')
        publications[stem]={key:next(item.value for item in ast.walk(node) if isinstance(item,ast.keyword) and item.arg==key)
          for key in ('modified_cylindrical_velocity_source_log_sectors','modified_absolute_pressure_over_Pstar2_ordinary_logR_rows')}
    own_adapter=asts.replay('current_modified_velocity_pressure_operator','modified_velocity_pressure_packet',
      dict(mp=mp,math=math,LABELS={'theta':UT,'axial':UZ,'radial':UR},UT=UT,UZ=UZ,UR=UR,P=P,
        zero_powers=zero_powers,exact_modified_radius=lambda *args:(rr,'same proved original radius')))
    original_adapter=asts.replay('current_modified_velocity_pressure_interfaces_operator','original_velocity_pressure_packet',
      dict(mp=mp,ORIGINAL_LABELS=ORIGINAL_LABELS,ORIGINAL_POWERS=ORIGINAL_POWERS,UT=UT,UZ=UZ,UR=UR,P=P,
        zero_powers=zero_powers,exact_modified_radius=lambda *args:(rr,'same proved original radius')))
    geometry=SimpleNamespace(logP=lp,pre='same actual pre',radius=lambda *args:(rr,'same proved original radius'))
    vpfield=SimpleNamespace(ctx=c,geometry=geometry)
    originalfield=SimpleNamespace(ctx=c,velocity=vpfield)
    def original_packet(region,q):
        return original_adapter(originalfield,dict(region=region,coordinate=q,original_native_velocity_pressure_source=pre))
    def own_packet(stem,region,q,increments,du,V):
        Unew=[a+b for a,b in zip(old['theta'],du)]
        delta_rows,Q=increment(c,z,delta,base,increments,Unew,old['theta'],V,du)
        original_pressure=[p0+pre['actual_normalized_primitive_y_derivative_axial5']['p'][0]]+pre['actual_normalized_primitive_y_derivative_axial5']['p'][1:]
        pressure=[a+b for a,b in zip(original_pressure,delta_rows['p'])]
        env=dict(c=c,original=old,Unew=Unew,Vnew=V,Q=Q,shifted_rows=shifted_rows,
          radial_increment=shifted_rows(Q,c.mpf('.5'),4),pressure=pressure,original_pressure=original_pressure,delta_rows=delta_rows)
        published={key:asts.evaluate(value,env) for key,value in publications[stem].items()}
        return own_adapter(vpfield,dict(published,region=region,coordinate=q))
    common_increment={key:s.Symbol('same_cumulative_'+key,real=True) for key in SHAPES}
    empty_increment=dict.fromkeys(SHAPES,s.Integer(0))
    du=[jet(C*s.Symbol('same_swirl_increment_y'+str(j),real=True)) for j in range(5)]
    V=[jet(C*s.Symbol('same_Vhat_y'+str(j),real=True)) for j in range(5)]
    mod='current_O3_modulated_histories';repair='current_O3_repaired_histories'
    cases={
      'buffer_inlet':(('original','O2_axial',1), (mod,'O2_buffer',0),'zero'),
      'buffer_transition':((mod,'O2_buffer',11),(mod,'O3_slope_mu',0),'shared_modulation'),
      'transition_power':((mod,'O3_slope_mu',1),(repair,'quiet_O3_power',0),'retained_history_flat_profile'),
      'buffer_modulation_start':(('original','O2_buffer',9),(mod,'O2_buffer',9),'zero'),
      'transition_modulation_end':((mod,'O3_slope_mu',s.Rational(1,2)),(mod,'O3_slope_mu',s.Rational(1,2)),'retained_history_flat_profile'),
      'repair_exit':((repair,'quiet_O3_power',2),('original','quiet_O3_power',2),'zero')}
    for name,(a,b,_,_) in REPAIR_EDGES.items():
        spec=(repair,'quiet_O3_power',s.Rational(a,b));cases[name]=(spec,spec,'retained_history_flat_profile')
    checks={};arguments={}
    for name in SEAMS:
        left,right,kind=cases[name]
        inc=empty_increment if kind=='zero' else common_increment
        dU=du if kind=='shared_modulation' else [zero]*5
        v=V if kind=='shared_modulation' else [zero]*5
        packets=[]
        for stem,region,q in (left,right):
            packets.append(original_packet(region,q) if stem=='original' else own_packet(stem,region,q,inc,dU,v))
        for label in (UT,UZ,UR,P):
            for j in range(5):
                for n in range(5-j):
                    values=[sum(Ps**s.Rational(str(powers[1]))*value for powers,value in packet['grids'][label][j,n]) for packet in packets]
                    if s.cancel(values[0]-values[1])!=0:raise ArithmeticError('Actual two-germ mixed4 packet differs: '+name+' '+label+' '+str((j,n)))
                    checks['actual_'+name+'_'+label+'_packet_source_y%d_Z%d'%(j,n)]=True
        arguments[name]=dict(source_programs=[left[0],right[0]],actual_regions=[left[1],right[1]],
          coordinates=[str(left[2]),str(right[2])],common_input_kind=kind,
          nonzero_cumulative_moment_and_pressure_kept=kind!='zero',
          actual_pressure_recurrence_and_radial_recovery_replayed=True)
    return dict(identities=checks,actual_source_germ_input_cases=arguments,
      both_actual_original_and_own_packet_programs_replayed=True,
      ordinary_profile4_and_whole_Z_M_Cp_seed_joins_required_before_projection=True,
      signed_sectors_retained_before_total_function_comparison=True,
      original_boundary_zero_Pstar1_follows_checked_flat_or_implicit_source_identity=True,passed=True)

def exact_modified_velocity_pressure_interfaces_theorem(field):
    """Compose existing actual function/FTC joins with full derivative consumers."""
    from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import CompliantPrePulseMixedC4
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    velocity=field.velocity;interfaces=field.interfaces;dispatch=velocity.heat.dispatch
    velocity.heat.assert_graph()
    if not velocity.acceptance_loaded or not interfaces.acceptance_loaded or interfaces is not velocity.heat.interfaces:
        raise ValueError('Checked same actual velocity and named source interfaces required')
    if field.controls is not velocity.source.repair.controls or interfaces.repaired is not velocity.source:
        raise ValueError('Physical interface germs must use one exact implicit control vector')
    if dispatch.pre.axial.__func__ is not CompliantPrePulseMixedC4.axial or dispatch.pre.power.__func__ is not CompliantPrePulseMixedC4.power:
        raise ValueError('Actual original boundary provider callable changed')
    joins=interfaces.theorem;vp=velocity.theorem
    owners=velocity.heat.registry.owners
    o2=owners['o2'].proof['current_actual_O2_source_and_four_endpoint_theorem']
    o3=owners['o3'].proof['current_actual_transition_source_and_power_endpoint']
    if any(not owners[name].acceptance_loaded or owners[name].physical is not velocity.geometry
      or owners[name].physical.pre is not dispatch.pre for name in ('o2','o3')):
        raise ValueError('Same checked original O2/O3 physical geometry and pre source required')
    if not o2['passed'] or not o3['passed'] or not all(o2['live_original_callable_bindings'].values()) or not all(o3['live_original_callable_bindings'].values()):
        raise ValueError('Checked original O2/O3 mixed4 source-function and geometry theorems required')
    for key in ('axial_buffer_same_actual_radius','buffer_transition_same_actual_radius'):
        if o2['exact_radius_identities'].get(key) is not True:raise ValueError('Original boundary radius equality missing: '+key)
    if not o3['actual_transition_power_same_Rw'] or not o2['ordinary_turnoff_cutoff_theorem']['passed']:
        raise ValueError('Actual ordinary turnoff jets and transition/power radius required')
    if not joins['passed'] or not vp['passed'] or not all(joins['identities'].values()) or not all(vp['identities'].values()):
        raise ValueError('Actual source joins and full velocity/pressure mixed4 theorem required')
    if tuple(joins['source_logR_axial_orders'])!=(4,5):raise ValueError('Insufficient source regularity for physical mixed4')
    profile=velocity.source.histories.candidate.theorem;history=velocity.source.histories.theorem
    if profile['identities'].get('actual_shared_phase_and_all_derivatives_at_old_seam') is not True:
        raise ValueError('Actual finite-N shared ordinary profile jets required')
    for key in ('actual_buffer_Ua_f_C_source_function','actual_transition_Ua_f_C_source_function','actual_continuous_transport_FTC_m',
      'actual_continuous_transport_FTC_p','actual_radial_recovery_same_cumulative_m_not_density','actual_radial_increment_Pstar_units'):
        if history['identities'].get(key) is not True:raise ValueError('Actual shared modulation/cumulative source identity missing: '+key)
    keys={
      'buffer_inlet':('consumed_actual_original_source_join_O2_turnoff_buffer','actual_zero_increment_keeps_original_p_row0'),
      'buffer_transition':('consumed_actual_original_source_join_Rd_buffer_slope_mu','actual_buffer_transition_shared_t'),
      'transition_power':('consumed_actual_original_source_join_Rw_slope_mu_power','actual_transition_power_shared_t','actual_zero_support_keeps_parent_partial_Cp'),
      'buffer_modulation_start':('actual_modulation_start_flat_cutoff_row4','actual_buffer_start_shared_t'),
      'transition_modulation_end':('actual_modulation_end_flat_cutoff_row4','modulation_end_uses_same_nonzero_cumulative_scalars_not_zero_history'),
      'repair_exit':('q2_same_implicit_root_physical_exit_uses_consumed_transport_and_tensor_theorems',)}
    for name in REPAIR_EDGES:keys[name]=('named_'+name+'_partial_values_and_FTC_jets_common_in_both_germs',
      'named_'+name+'_all_local_bump_jets_zero_by_checked_flat_limits')
    if set(keys)!=set(SEAMS):raise ValueError('Twelve source-function proof cases incomplete')
    for name,required in keys.items():
        if any(joins['identities'].get(key) is not True for key in required):raise ValueError('Actual source equality missing: '+name)
    asts=SourceAST();checks={}
    checks['consumed_checked_original_O2_axial_buffer_and_buffer_transition_radii']=True
    checks['consumed_checked_original_O3_transition_power_radius_and_mixed4_functions']=True
    checks['same_actual_original_geometry_pre_and_all_original_callable_bindings']=True
    geometry_proof=dispatch.local.theorem
    for key in ('actual_quiet_log_radius_matches_original_q_over_same_Tw_phase',
      'actual_O2_buffer_original_radius_retained','actual_O3_slope_mu_original_radius_retained',
      'actual_geometry_logPstar_is_same_original_pre_scale_value'):
        if geometry_proof['identities'].get(key) is not True:raise ValueError('Actual own radius/scale binding missing: '+key)
        checks['consumed_'+key]=True
    h=velocity.source.histories
    if h.pre is not dispatch.pre or h.candidate.pre is not dispatch.pre or h.ctx is not velocity.ctx:
        raise ValueError('Original and own amplitude must use the same actual pre/source context')
    asts.expression('current_O3_modulated_histories','__init__','self.Ua',
      wanted="c.mpf(endpoints(self.candidate.direction.bounds['actual_Ua']))")
    asts.expression('current_O3_modulated_histories','history','original',wanted='raw_pre_velocity_rows(c,pre)')
    asts.expression('current_O3_modulated_histories','history','(delta_rows, Q)',
      wanted='increment_rows(c,z,self.delta,self.Ua,scalars,Unew,original[\'theta\'],Vnew,du)')
    asts.expression('current_O3_repaired_histories','history','original',
      wanted="raw_pre_velocity_rows(c,parent['current_original_pre_source'])")
    asts.expression('current_O3_repaired_histories','history','Ahat',
      wanted='(1+square(z)).reciprocal()*(self.histories.Ua*self.f2)')
    asts.expression('current_O3_repaired_histories','history','(delta_rows, Q)',
      wanted="increment_rows(c,z,self.histories.delta,self.histories.Ua,scalar,Unew,original['theta'],Vnew,du)")
    if vp['actual_absolute_pressure_mixed4_source_theorem']['identities'].get('physical_absolute_pressure_has_unchanged_Pstar_squared_units') is not True:
        raise ValueError('Same actual absolute-pressure Pstar squared unit required')
    checks['actual_common_base_is_same_histories_Ua_from_checked_original_inlet']=True
    checks['actual_native_amplitude_rows_are_same_original_pre_Ua_f_C_functions']=True
    checks['actual_repair_Ahat_and_recovery_use_same_Ua_not_new_normalization']=True
    checks['actual_original_and_own_Pstar_and_absolute_pressure_units_same_geometry_source']=True
    def zero(name,a,b):
        if s.cancel(a-b)!=0:raise ArithmeticError('Actual physical interface identity: '+name)
        checks[name]=True
    for target,wanted in (('pre','pre_owner.axial(Z,q)'),('pre','pre_owner.axial(Z,buffer_offset=q)'),
      ('pre','pre_owner.power(Z,q/field.velocity.heat.dispatch.Tw)')):
        asts.expression('current_modified_velocity_pressure_interfaces_operator','original_velocity_pressure_source',target,wanted=wanted)
    asts.expression('current_modified_velocity_pressure_interfaces_operator','original_velocity_pressure_packet','raw',
      wanted="pre['physical_velocity_pressure_y_Z_mixed4']")
    asts.expression('current_modified_velocity_pressure_interfaces_operator','physical_boundary_germ','(kind, source)',
      wanted='boundary_source(field,name,Z,side)')
    asts.expression('current_modified_velocity_pressure_interfaces_operator','physical_boundary_germ','packet',
      wanted="original_velocity_pressure_packet(field,source) if kind=='original_native_germ' else modified_velocity_pressure_packet(field.velocity,source)")
    asts.expression('current_modified_velocity_pressure_interfaces_operator','physical_boundary_germ','physical',
      wanted="modified_physical_velocity_pressure(field.velocity,packet,source['Z'],log_tau,theta)")
    asts.method('current_modified_velocity_pressure_interfaces_operator','boundary_source')
    asts.method('current_modified_tensor_interfaces','endpoint_source')
    # The native germs used in the coupled replay have V=0: the buffer
    # branch explicitly uses five zero cutoff rows, the phase1 turnoff
    # has the same checked flat cutoff, and slope/power pass zero rows.
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='[c.mpf(0)]*5')
    asts.expression('pre_pulse_mixed_C4','axial','V',wanted='[z*(4*b) for b in B]')
    import ast
    for method in ('slope_mu','power'):
        packet_call=next(node for node in ast.walk(asts.method('pre_pulse_mixed_C4',method))
          if isinstance(node,ast.Call) and ast.unparse(node.func)=='self.packet')
        if ast.dump(packet_call.args[6])!=ast.dump(ast.parse('[zero]*5',mode='eval').body):
            raise ValueError('Original native zero axial rows changed: '+method)
    checks['actual_original_native_zero_V_and_all_ordinary_rows_bound_before_coupled_replay']=True
    # Bind the flat transition endpoint t=1 explicitly, in addition to
    # the two modulation support edges. The finite-profile V program
    # uses these actual cutoff rows, not a fresh zero axial placeholder.
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    cutoff=asts.replay('current_O3_finite_frequency_profiles','cutoff_rows',
      dict(sigma_jets=lambda ctx,v:[s.Integer(0 if v<=0 else 1)]+[s.Integer(0)]*4,
        product_rows=product_rows,math=math))
    for name,t in (('buffer_inlet',s.Integer(-11)),('modulation_start',s.Integer(-2)),
      ('modulation_end',s.Rational(1,2)),('transition_power',s.Integer(1))):
        rows=cutoff(None,t)
        for j,row in enumerate(rows):zero('actual_'+name+'_cutoff_flat_row'+str(j),row,0)
    asts.expression('current_O3_finite_frequency_profiles','profile','chi',wanted='cutoff_rows(c,t)')
    asts.expression('current_O3_finite_frequency_profiles','profile','Vnew',
      wanted="[r*(-c.sqrt(mu)/(2*c.pi*n)) for r in product_rows(old['theta'],product_rows(chi,cos2))]")
    asts.expression('current_O3_modulated_histories','history','Vnew',wanted='[row*0 for row in Unew]')
    checks['actual_flat_cutoff_rows_force_original_and_own_zero_axial_germs_before_packet_comparison']=True
    checks['actual_buffer_transition_keeps_shared_nonzero_axial_profile_function']=True
    # Original packet is already ordinary. Replay every row and signed
    # unit through the actual adapter; no artificial radial/axial sectors.
    rr,lp=s.symbols('same_logR same_logPstar',real=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),ln=s.log)
    raw={name:{'y%d_Z%d'%(j,n):s.Symbol('original_'+str(i)+'_y%d_Z%d'%(j,n))
      for j in range(5) for n in range(5-j)} for i,name in enumerate(ORIGINAL_LABELS.values())}
    view=dict(region='quiet_O3_power',coordinate=s.Integer(2),original_native_velocity_pressure_source=dict(physical_velocity_pressure_y_Z_mixed4=raw))
    fn=asts.replay('current_modified_velocity_pressure_interfaces_operator','original_velocity_pressure_packet',
      dict(mp=mp,ORIGINAL_LABELS=ORIGINAL_LABELS,ORIGINAL_POWERS=ORIGINAL_POWERS,UT=UT,UZ=UZ,UR=UR,P=P,
        zero_powers=zero_powers,exact_modified_radius=lambda *args:(rr,'same source radius')))
    packet=fn(SimpleNamespace(ctx=c,velocity=SimpleNamespace(geometry=SimpleNamespace(logP=lp))),view)
    for label,name in ORIGINAL_LABELS.items():
        for j in range(5):
            for n in range(5-j):
                powers,value=packet['grids'][label][j,n][0]
                zero('actual_original_boundary_'+label+'_y%d_Z%d'%(j,n),value,raw[name]['y%d_Z%d'%(j,n)])
                if powers[1]!=ORIGINAL_POWERS[label]:raise ValueError('Original physical source unit changed')
    checks['actual_original_radial_axial_Pstar0_no_fake_Pstar1_sector']=True
    checks['actual_original_rows_not_multiplied_by_Taylor_factorial_twice']=True
    # Both complete radial derivatives are obtained from the actual
    # linear Q recovery of M/V, including the required M_Z fifth row.
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    jet=lambda expr,order=5:Projection.function(c,expr,order)
    y=s.Symbol('same_logR_coordinate',real=True);delta,Ps=s.symbols('same_delta same_Pstar',real=True)
    M0=s.Function('same_M0')(y,Z);M1=s.Function('same_M1')(y,Z)
    V0=s.Function('same_V0')(y,Z);V1=s.Function('same_V1')(y,Z)
    rows=lambda expr:[jet(s.diff(expr,y,j)) for j in range(5)]
    native_history=asts.method('pre_pulse_mixed_C4','physical_mixed')
    own_history=asts.method('current_O3_modulated_histories_operator','increment_rows')
    import ast
    native_append=next(node for node in ast.walk(native_history) if isinstance(node,ast.Call)
      and ast.unparse(node.func)=="rows['m'].append")
    own_append=next(node for node in ast.walk(own_history) if isinstance(node,ast.Call)
      and ast.unparse(node.func)=="d['m'].append")
    native_M=[jet(M0)];own_M=[jet(M1)];complete_M=[jet(M0+Ps*M1)]
    native_V=rows(V0);own_V=rows(V1);complete_V=rows(V0+Ps*V1)
    for j in range(4):
        native_M.append(asts.evaluate(native_append.args[0],dict(V=native_V,j=j,rows=dict(m=native_M))))
        own_M.append(asts.evaluate(own_append.args[0],dict(Vnew=own_V,j=j,d=dict(m=own_M))))
        complete_M.append(asts.evaluate(native_append.args[0],dict(V=complete_V,j=j,rows=dict(m=complete_M))))
    for j in range(5):
        for n in range(6-j):
            zero('actual_common_M_FTC_required_radial_row_y%d_Z%d'%(j,n),
              (native_M[j]+Ps*own_M[j]-complete_M[j])[n]*math.factorial(n),0)
    original_Q=asts.expression('pre_pulse_mixed_C4','physical_mixed','Q')
    own_Q=asts.expression('current_O3_modulated_histories_operator','increment_rows','Q')
    env=dict(z=jet(Z),delta=delta,square=lambda v:v*v,derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1))
    env['L']=asts.evaluate(asts.expression('current_O3_modulated_histories_operator','increment_rows','L'),env)
    Q0=asts.evaluate(original_Q,dict(env,V=rows(V0),rows=dict(m=rows(M0))))
    Q1=asts.evaluate(own_Q,dict(env,Vnew=rows(V1),d=dict(m=rows(M1))))
    total=asts.evaluate(original_Q,dict(env,V=rows(V0+Ps*V1),rows=dict(m=rows(M0+Ps*M1))))
    for j in range(5):
        for n in range(5-j):
            zero('actual_radial_recovery_complete_mixed4_y%d_Z%d'%(j,n),
              (Q0[j]+Ps*Q1[j]-total[j])[n]*math.factorial(n),0)
    checks['same_actual_M_axial5_retains_radial_axial4_dependencies']=True
    pressure=vp['actual_absolute_pressure_mixed4_source_theorem']
    if not pressure['passed'] or not all(pressure['identities'].values()):raise ValueError('Actual absolute P0+own Cp mixed4 source proof required')
    # Already admitted actual function germs identify U,V,M and Cp before
    # differentiation. Their ordinary rows and shared P0 are consumed by
    # the same actual operators. Record every field/order dependency.
    coupled=coupled_mixed4_packet_source_joins(field,asts)
    checks.update(coupled['identities'])
    checks['same_common_source_radius_and_units_precede_Cartesian_and_time_bounds']=True
    checks['actual_common_mixed4_functions_consumed_by_unchanged_spatial4_time1_operators']=True
    asts.expression('current_modified_velocity_pressure_interfaces','interface','left',
      wanted="physical_boundary_germ(self,name,Z,'left',lt,theta)")
    asts.expression('current_modified_velocity_pressure_interfaces','interface','right',
      wanted="physical_boundary_germ(self,name,Z,'right',lt,theta)")
    return dict(identities=checks,passed=True,input_hashes={**joins['input_hashes'],**vp['input_hashes'],**o2['input_hashes'],**o3['input_hashes'],**asts.hashes},
      consumed_original_O2_exact_radius_and_full_mixed4_source_theorem=o2,
      consumed_original_O3_exact_radius_and_full_mixed4_source_theorem=o3,
      consumed_actual_twelve_germ_source_identity_keys=keys,
      consumed_checked_modified_source_interface_identities=len(joins['identities']),
      consumed_checked_velocity_pressure_source_identities=len(vp['identities']),
      actual_original_boundary_sector_inventory=ORIGINAL_POWERS,
      actual_coupled_velocity_pressure_packet_source_joins=coupled,
      coupled_input_function_bindings=dict(original_amp_C='same actual original pre ordinary theta rows; actual Ua*f*C source theorem',
        base='same histories.Ua from checked direction actual original inlet',Pstar='same original pre/geometry logPstar',
        original_zero_V='actual native flat turnoff/buffer/slope/power source programs',
        own_zero_V='actual cutoff at t=-11,-2,1/2,1 or checked named bump flat rows; quiet branch explicit zero',
        shared_V='same actual finite-N profile at buffer11/O3 offset0',
        radius='checked original O2/O3 geometry identities plus actual own exact_radius source theorem'),
      actual_named_partial_FTC_source_theorem_consumed=True,
      exact_original_five_history_and_own_pressure_function_joins_consumed=True,
      mixed4_trace_derivation='actual ordinary profile4 traces + whole-Z M/Cp seeds + actual M and pressure FTC recurrences; radial recovery consumes M axial5',
      complete_velocity_pressure_mixed4_and_fixed_x_time1_trace_scope=True,
      full33_velocity_dispatch_total_energy_cones_recursion_global_NS_remain_open=True)
