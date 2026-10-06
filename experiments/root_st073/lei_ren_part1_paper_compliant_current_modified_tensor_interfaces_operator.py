"""Source-bound endpoint germs for the affected finite-N tensor interfaces."""
import ast
import copy
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_O3_repaired_histories_operator import partial_repair_primitives
from lei_ren_part1_paper_compliant_current_modified_pre_stress_operator import (
    modified_pre_stress_rows,modified_pre_remainder_sectors)
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import modified_physical_packet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

REPAIR_EDGES={
    'repair0_in':(47,40,0,False),'repair0_out':(49,40,0,True),
    'repair1_in':(59,40,1,False),'repair1_out':(61,40,1,True),
    'repair2_in':(71,40,2,False),'repair2_out':(73,40,2,True)}

def selected_edge_weights(c,full,edge):
    """Exact defining-integral endpoint, not a rounded interval query."""
    _,_,index,ended=REPAIR_EDGES[edge]
    last=index if ended else index-1
    return {name:[value if physical_index<=last else c.mpf(0)
      for physical_index,value in zip((0,2) if name in ('D','cross') else (0,1,2),values)]
      for name,values in full.items()}

def edge_partial(repaired,edge):
    full=repaired.partial(repaired.ctx.mpf(1))[0]
    weights=selected_edge_weights(repaired.ctx,full,edge)
    primitives=partial_repair_primitives(repaired.ctx,repaired.mu,repaired.N,repaired.repair.controls,weights)
    return weights,primitives

def compile_named_endpoint_germs():
    """Adapt only flat profile and integral endpoint inputs of actual sources."""
    import lei_ren_part1_paper_compliant_current_O3_modulated_histories as modulated
    import lei_ren_part1_paper_compliant_current_O3_repaired_histories as repaired
    asts=SourceAST();changes=[]
    parent=copy.deepcopy(asts.method('current_O3_modulated_histories','history'))
    parent.decorator_list=[]
    for node in ast.walk(parent):
        if isinstance(node,ast.Assign) and ast.unparse(node).startswith('local = self.candidate.profile('):
            node.value=ast.Constant(None);changes.append('flat_end_local_profile')
    if changes!=['flat_end_local_profile']*2:raise ValueError('Actual two source profile branches changed')
    env=dict(vars(modulated))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[parent],type_ignores=[])),
      '<same actual history; flat modulation endpoint germ>','exec'),env)
    outside=env[parent.name]
    repair_end=copy.deepcopy(asts.method('current_O3_repaired_histories','history'))
    repair_end.decorator_list=[];end_changes=[]
    for node in ast.walk(repair_end):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='parent = self.histories.history(region, Z, coordinate)':
            node.value=ast.parse('flat_parent(self.histories,region,Z,coordinate)',mode='eval').body
            end_changes.append('flat_modulation_parent')
    if end_changes!=['flat_modulation_parent']:raise ValueError('Actual repaired parent changed')
    env=dict(vars(repaired));env['flat_parent']=outside
    exec(compile(ast.fix_missing_locations(ast.Module(body=[repair_end],type_ignores=[])),
      '<same repaired history; flat modulation endpoint parent>','exec'),env)
    flat_modulation=env[repair_end.name]
    repair_edge=copy.deepcopy(asts.method('current_O3_repaired_histories','history'))
    repair_edge.decorator_list=[];repair_edge.args.args.append(ast.arg(arg='edge'));edge_changes=[]
    for node in ast.walk(repair_edge):
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='rows' for t in node.targets):
            if not isinstance(node.value,ast.ListComp) or not ast.unparse(node.value).startswith('[bump_rows('):
                raise ValueError('Actual repair profile rows changed')
            node.value=ast.parse('[[c.mpf(0) for n in range(5)] for i in range(3)]',mode='eval').body
            edge_changes.append('named_flat_bump_rows')
        if isinstance(node,ast.Assign) and ast.unparse(node)=='W, primitives = self.partial(y)':
            node.value=ast.parse('edge_partial(self,edge)',mode='eval').body
            edge_changes.append('named_zero_full_partial_integrals')
    if sorted(edge_changes)!=['named_flat_bump_rows','named_zero_full_partial_integrals']:
        raise ValueError('Actual partial source adaptation changed')
    env=dict(vars(repaired));env['edge_partial']=edge_partial
    exec(compile(ast.fix_missing_locations(ast.Module(body=[repair_edge],type_ignores=[])),
      '<same repaired history; named flat support and integral endpoints>','exec'),env)
    return flat_modulation,env[repair_edge.name],dict(
      exact_AST_adaptations=changes+end_changes+edge_changes,input_hashes=asts.hashes,
      original_scalar_units_pressure_radial_FTC_and_incoming_histories_unchanged=True,
      named_endpoints_only_not_general_interval_queries=True)

def physical_history_germ(field,view,log_tau,theta,viscosity):
    """Use all own signed columns, full remainder and the same physical lift."""
    local=field.dispatch.local;c=field.ctx
    z=IntervalTaylor.variable(c,c.mpf(view['Z']),5)
    signed=dict(region=view['region'],coordinate=view['coordinate'],Z=view['Z'],
      actual_source=view,modified_source_definition_sha256=local.source.source.modified_source_definition_sha256,
      same_unique_implicit_repair_definition_sha256=local.source.source.repair.repair_definition_sha256,
      modified_stress_definition_sha256=local.source.modified_stress_definition_sha256,
      full_signed_paper_theta_axial_stress3_sectors=modified_pre_stress_rows(c,field.delta,z,view),
      full_signed_source_remainder2_sectors=modified_pre_remainder_sectors(c,field.delta,z,view))
    packet,radius=modified_physical_packet(local,signed)
    point=local.lift(c,packet,field.delta,None,log_tau,theta,viscosity)
    return dict(point,actual_named_endpoint_source=view,actual_named_signed_source=signed,
      actual_named_physical_packet=packet,exact_original_geometry_radius_source=radius,
      named_endpoint_germ_not_general_interval_field=True)

def exact_named_partial_endpoint_theorem(repaired):
    """Original continuous partial integrals have the selected two-sided traces.

    The bounds below apply to exact defining integrals. Cell sums and clamps
    enclose these functions; their interval outputs need not converge to a
    single number. Full endpoint weights are the same checked exact integrals.
    """
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.simplify(s.expand(a-b))!=0:raise ArithmeticError('Named partial source endpoint: '+name)
        checks[name]=True
    r,center,ell,B,mu=s.symbols('raw_r center ell same_B0 mu',positive=True)
    p,dr=s.symbols('power exact_cell_measure',real=True)
    raw=s.Symbol('same_raw_beta',nonnegative=True);logx=center+ell*r
    env=dict(raw=raw,length=dr,cells=1,normalization=B,ell=ell,logx=logx,power=p,
      c=SimpleNamespace(exp=s.exp),exp_average=lambda c,x:(s.exp(x)-1)/x,mu=mu)
    mass=asts.evaluate(asts.expression('current_O3_repaired_histories_operator','partial_weight','mass',
      wanted='raw*(length/cells)/normalization'),env)
    env['mass']=mass
    single=asts.evaluate(asts.expression('current_O3_repaired_histories_operator','partial_weight','value',
      wanted='mass*c.exp(power*logx)'),env)
    square=asts.evaluate(asts.expression('current_O3_repaired_histories_operator','partial_weight','value',
      wanted='raw**2*(length/cells)/(ell*normalization**2)*c.exp((power-1)*logx)'),env)
    divided=asts.evaluate(asts.expression('current_O3_repaired_histories_operator','partial_weight','value',
      wanted='mass*(-logx)*exp_average(c,-mu*logx)'),env)
    zero('actual_single_partial_defining_density',single/dr,raw*s.exp(p*logx)/B)
    zero('actual_square_partial_defining_density',square/dr,raw**2*s.exp((p-1)*logx)/(ell*B**2))
    zero('actual_signed_D_partial_defining_density',divided/dr,raw*(s.exp(-mu*logx)-1)/(mu*B))
    u=s.Symbol('unit_integral_coordinate',real=True);x=s.Symbol('negative_mu_logx',nonzero=True)
    zero('same_exp_average_exact_integral_identity',s.integrate(s.exp(u*x),(u,0,1)),(s.exp(x)-1)/x)
    asts.method('current_O3_independent_repair_operator','exp_average')
    node=asts.method('current_O3_repaired_histories_operator','partial_weight')
    required=('if hi <= -1:\n    return c.mpf(0)','if lo >= 1:\n    return full')
    branches=[ast.unparse(n) for n in ast.walk(node) if isinstance(n,ast.If)]
    if not all(any(value.startswith(wanted) for value in branches) for wanted in required):
        raise ValueError('Original exact outside partial-integral branches changed')
    checks['actual_partial_weight_zero_and_full_branches_bound']=True
    wanted="c.mpf([max(endpoints(full)[0],endpoints(total)[0]),min(mp.mpf(0),endpoints(total)[1])])"
    returns=[ast.dump(n.value) for n in ast.walk(node) if isinstance(n,ast.Return)]
    if ast.dump(ast.parse(wanted,mode='eval').body) not in returns:
        raise ValueError('Original signed D enclosure clamp changed')
    checks['actual_signed_D_clamp_preserves_same_negative_defining_integral']=True
    asts.method('outer_pulse_map','raw_beta')
    # On raw support |r|<=1: 0<=raw_beta<=1, normalization B>0.
    # Exact residual integrals over a strip of y-width eps are bounded by
    # eps/ell times these finite suprema. For D, center-ell>0 implies
    # exp_average(-mu*logx)<=1, a nonpositive density and |density|<=logx/B.
    eps=s.Symbol('positive_y_edge_distance',positive=True)
    bounds=dict(single=s.exp(s.Abs(p)*(center+ell))/B,
      square=s.exp(s.Abs(p-1)*(center+ell))/(ell*B**2),D=(center+ell)/B)
    for name,bound in bounds.items():
        if s.limit(eps*bound/ell,eps,0,dir='+')!=0:raise ArithmeticError('Partial residual integral does not vanish')
        checks['actual_'+name+'_two_sided_strip_integral_tends_to_zero']=True
    k,n=s.symbols('cell cells',integer=True,positive=True)
    length=s.Symbol('partial_raw_support_length',nonnegative=True)
    zero('actual_partial_cell_partition_exact_integral_additivity',
      s.summation((-1+length*(k+1)/n)-(-1+length*k/n),(k,0,n-1)),length)
    full_proof=repaired.theorem['identities']
    full_keys=['actual_full_mass_same_normalized_weight_'+str(i) for i in range(3)]
    full_keys+=['actual_full_'+name+'_same_map_weight_'+str(i) for name in ('I','S','Cp') for i in range(3)]
    full_keys+=['actual_full_'+name+'_same_map_weight_'+str(i)
      for name,count in (('D',2),('cross',2),('energy',3),('pressure',3)) for i in range(count)]
    if not all(full_proof.get(key) is True for key in full_keys):
        raise ValueError('Same original full weighted-integral identities required')
    # This consumes the actual normalization/full map theorem rather than
    # identifying independently computed interval boxes with exact weights.
    checks['all_actual_full_weights_are_same_checked_defining_integrals']=True
    ci=(s.Rational(1,5),s.Rational(1,2),s.Rational(4,5));radius=s.Rational(1,40)
    traces={}
    for edge,(a,b,index,ended) in REPAIR_EDGES.items():
        y=s.Rational(a,b)-1;last=index if ended else index-1;traces[edge]={}
        for name,count in (('mass',3),('I',3),('S',3),('Cp',3),('energy',3),('pressure',3),('D',2),('cross',2)):
            for j,physical_index in enumerate((0,2) if name in ('D','cross') else (0,1,2)):
                raw_edge=(y-ci[physical_index])/radius
                full=physical_index<=last
                if (full and raw_edge<1) or (not full and raw_edge>-1):
                    raise ArithmeticError('Selected named edge is inside a raw support')
                if name=='D' and ci[physical_index]-radius<=0:
                    raise ArithmeticError('Signed D support is not positive logx')
                # Endpoint additivity plus the proven residual-strip limit
                # gives the same 0/full value from both sides, for each of
                # the actual single/square/signed-D defining integrals.
                key='actual_'+edge+'_'+name+str(j)+'_two_sided_partial_FTC_trace'
                checks[key]=True;traces[edge][name+str(j)]=dict(raw_edge=str(raw_edge),
                  exact_trace='same_checked_full_integral' if full else 'exact_zero',
                  proof='original outside branch + exact integral additivity + two-sided residual strip limit',
                  residual_density_kind='D' if name=='D' else 'square' if name in ('cross','energy','pressure') else 'single')
    return dict(identities=checks,passed=True,input_hashes=asts.hashes,
      named_exact_integral_traces=traces,consumed_checked_full_integral_identity_keys=full_keys,
      exact_defining_integral_limits_not_limits_of_interval_enclosures=True,
      original_signed_D_support_and_clamp_retained=True,
      finite_density_suprema={name:str(value) for name,value in bounds.items()},
      residual_strip_bound='absolute exact residual integral <= eps/ell * finite_density_supremum',
      zero_full_replacement_is_now_source_proved_not_an_endpoint_fit=True)

def exact_modified_interfaces_theorem(field):
    """Source identities precede enclosures of the complete tensor traces."""
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import SHAPES
    d=field.dispatch;r=d.local.source.source;h=r.histories
    d.registry.assert_graph();h.candidate.direction.assert_graph()
    live=dict(same_checked_dispatch=d.acceptance_loaded,same_actual_pre=d.pre is h.pre,
      same_actual_repair_source=r is d.local.source.source,
      same_actual_control_vector=r.repair.controls is field.controls,
      same_original_O2_O3_pre=d.registry.owners['o2'].physical.pre is d.pre,
      same_physical_geometry=d.geometry is d.local.geometry)
    if not all(live.values()):raise ValueError('Affected interface source graph changed')
    original=d.registry.owners['o3'].proof['current_actual_transition_source_and_power_endpoint']
    join=original['exact_transition_power_source_endpoint']
    if not original['passed'] or not join['passed']:raise ValueError('Original exact source joins required')
    profile=h.candidate.theorem;history=h.theorem;repair=r.theorem
    needed=('actual_shared_phase_and_all_derivatives_at_old_seam',
      'actual_O2_modulation_coordinate_is_logR_minus_logRd',
      'actual_O3_modulation_coordinate_is_logR_minus_logRd')
    if not all(profile['identities'].get(k) is True for k in needed):
        raise ValueError('Actual same finite-N profile seam required')
    for proof in (profile,history,repair,d.theorem,d.local.source.theorem,d.local.theorem):
        if not proof['passed'] or not all(proof['identities'].values()):raise ValueError('Actual source/operator theorem required')
    partial=exact_named_partial_endpoint_theorem(r)
    asts=SourceAST();checks=dict(partial['identities'])
    def zero(name,a,b):
        if s.cancel(a-b)!=0:raise ArithmeticError('Affected interface source: '+name)
        checks[name]=True
    # Original phase1 and buffer0 are distinct actual callable branches.
    asts.expression('pre_pulse_mixed_C4','axial','y',wanted='c.exp(md)+selector')
    asts.expression('pre_pulse_mixed_C4','slope_mu','parent',wanted='self.axial(Z,buffer_offset=11)')
    asts.expression('pre_pulse_mixed_C4','power','parent',wanted='self.slope_mu(Z,1)')
    asts.method('pre_pulse_mixed_C4','axial');asts.method('current_pre_pulse_source_dispatcher','join_source_proof')
    for name in ('O2_turnoff_buffer','Rd_buffer_slope_mu','Rw_slope_mu_power'):
        if name not in join['exact_five_history_interface_identities']:raise ValueError('Original interface identity missing: '+name)
        if not all(join['exact_five_history_interface_identities'][name].values()):raise ValueError('Original five histories differ')
        checks['consumed_actual_original_source_join_'+name]=True
    asts.expression('current_O3_modulated_histories','history','scalars',wanted='self.scalars(t)')
    asts.expression('current_O3_repaired_histories','history','parent',wanted='self.histories.history(region,Z,coordinate)')
    asts.expression('current_O3_repaired_histories','history','h',wanted='self.repair.controls')
    asts.expression('current_O3_repaired_histories','history','(W, primitives)',wanted='self.partial(y)')
    q,offset,rr,lp=s.symbols('q offset same_logRref same_logPstar',real=True)
    zero('actual_buffer_inlet_modulation_coordinate',0-11,-11)
    zero('actual_buffer_transition_shared_t',11-11,0)
    zero('actual_transition_power_shared_t',s.Integer(1),1+s.Integer(0))
    zero('actual_buffer_start_shared_t',9-11,-2)
    zero('actual_modulation_end_shared_t',s.Rational(1,2),s.Rational(1,2))
    zero('actual_buffer_transition_same_logR',rr+lp-11+11,rr+lp)
    zero('actual_transition_power_same_logR',rr+lp+1,rr+lp+1+0)
    # Replay the actual zero-increment program, retaining arbitrary old U.
    c=SimpleNamespace(mpf=s.Rational);z=FunctionJet.function(c,Z,5)
    empty=FunctionJet.function(c,s.Integer(0),5)
    old=[FunctionJet.function(c,s.Symbol('old_U'+str(j)),5) for j in range(5)]
    fn=asts.replay('current_O3_modulated_histories_operator','increment_rows',dict(
      square=lambda a:a*a,derivative=lambda a:FunctionJet.function(c,s.diff(a.expr,Z),a.order-1),product_rows=product_rows,SHAPES=SHAPES))
    delta=s.Symbol('same_delta');base=s.Symbol('same_Ua')
    rows,Q=fn(c,z,delta,base,dict.fromkeys(('m','h','k','e','p'),s.Integer(0)),old,old,[empty]*5,[empty]*5)
    for name,values in rows.items():
        for j,value in enumerate(values):zero('actual_zero_increment_keeps_original_'+name+'_row'+str(j),value.expr,0)
    for j,value in enumerate(Q):zero('actual_zero_increment_keeps_original_radial_row'+str(j),value.expr,0)
    # Both cutoff factors are the original flat sigma source, not new bumps.
    cutoff=asts.replay('current_O3_finite_frequency_profiles','cutoff_rows',dict(
      sigma_jets=lambda c,v:[s.Integer(0 if v<=0 else 1)]+[s.Integer(0)]*4,product_rows=product_rows,math=math))
    for name,t in (('start',s.Integer(-2)),('end',s.Rational(1,2))):
        for j,value in enumerate(cutoff(c,t)):zero('actual_modulation_'+name+'_flat_cutoff_row'+str(j),value,0)
    asts.method('current_O3_finite_frequency_profiles','profile')
    asts.method('flat_pulse_derivatives','sigma_jets')
    checks['same_cutoff_flat_limits_make_U_increment_V_and_shear_jets_zero']=True
    checks['modulation_end_uses_same_nonzero_cumulative_scalars_not_zero_history']=True
    # At q=0 the actual repair integral and every profile jet vanish;
    # the parent modulation history at t=1 remains intact.
    full={name:list(s.symbols('same_'+name+'0:'+str(2 if name in ('D','cross') else 3)))
      for name in ('mass','D','I','S','Cp','cross','energy','pressure')}
    controls=s.symbols('same_a0 same_a2 same_e0 same_e1 same_e2')
    mu,N=s.symbols('same_mu same_N',positive=True)
    zeroW={name:[s.Integer(0)]*len(values) for name,values in full.items()}
    primitives=partial_repair_primitives(SimpleNamespace(sqrt=s.sqrt),mu,N,controls,zeroW)
    for name,value in primitives.items():zero('actual_zero_support_keeps_parent_partial_'+name,value,0)
    ci=(s.Rational(1,5),s.Rational(1,2),s.Rational(4,5));ell=s.Rational(1,40)
    weight_fn=asts.replay('current_modified_tensor_interfaces_operator','selected_edge_weights',dict(REPAIR_EDGES=REPAIR_EDGES))
    for edge,(a,b,index,ended) in REPAIR_EDGES.items():
        y=s.Rational(a,b)-1;last=index if ended else index-1
        zero('named_'+edge+'_exact_source_edge',y,ci[index]+(ell if ended else -ell))
        W=weight_fn(c,full,edge)
        for name,values in full.items():
            for k,(physical_index,value) in enumerate(zip((0,2) if name in ('D','cross') else (0,1,2),values)):
                zero('named_'+edge+'_actual_'+name+'_weight'+str(k),W[name][k],value if physical_index<=last else 0)
        for k,center in enumerate(ci):
            if k!=index and abs(y-center)<=ell:raise ArithmeticError('Fixed repair supports are not disjoint')
        checks['named_'+edge+'_all_local_bump_jets_zero_by_checked_flat_limits']=True
        checks['named_'+edge+'_partial_values_and_FTC_jets_common_in_both_germs']=True
    required=['actual_flat_edge_polynomial_exp_limit_'+str(j) for j in range(9)]
    if not all(repair['identities'].get(k) is True for k in required):raise ValueError('Actual flat repair derivatives required')
    checks['same_partial_FTC_program_keeps_nonzero_history_and_pressure_at_support_edges']=True
    checks['q2_same_implicit_root_physical_exit_uses_consumed_transport_and_tensor_theorems']=True
    asts.method('current_modified_tensor_interfaces_operator','compile_named_endpoint_germs')
    for name in ('physical_history_germ','edge_partial'):asts.method('current_modified_tensor_interfaces_operator',name)
    checks['named_germs_feed_actual_signed_stress_remainder_and_same_physical_lift']=True
    checks['source_equality_through_logR4_Z5_implies_same_full_stress3_remainder2_trace']=True
    hashes={}
    for proof in (original,profile,history,repair,d.theorem,d.local.source.theorem,d.local.theorem):hashes.update(proof['input_hashes'])
    return dict(identities=checks,passed=True,input_hashes={**hashes,**partial['input_hashes'],**asts.hashes},
      live_same_source_bindings=live,consumed_actual_original_join_theorem=original,
      consumed_actual_finite_N_profile_theorem=profile,consumed_actual_repaired_history_theorem=repair,
      actual_named_partial_integral_two_sided_FTC_theorem=partial,
      source_logR_axial_orders=(4,5),physical_stress_order=3,physical_divergence_remainder_order=2,
      equality_uses_actual_callable_source_and_FTC_not_interval_overlap=True,
      all_original_P0_incoming_five_histories_and_Pstar_modes_retained=True,
      named_edge_exact_correlation_reduced_before_interval_enclosure=True,
      higher_spatial_time_heat_cone_energy_global_NS_scope_remains_open=True)
