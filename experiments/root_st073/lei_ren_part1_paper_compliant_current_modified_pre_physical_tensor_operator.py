"""Actual-source packet and original full physical tensor/remainder pullback."""
import ast
import copy
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_modified_pre_stress import PREFIX,sha
from lei_ren_part1_paper_compliant_current_modified_pre_stress import endpoints

def compiled_modified_pre_physical_lift():
    """Change only the full remainder input; physical operators stay original."""
    import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as original
    asts=SourceAST();fn=copy.deepcopy(asts.method('pulse_end_physical_C2','lift_physical_packet'))
    fn.decorator_list=[];changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='source_errors = pulse_remainder_sectors(c, delta, None, z, velocity)':
            node.value=ast.parse("packet['actual_modified_source_remainder_sectors']",mode='eval').body
            changes.append('own_signed_source_remainder_input')
    if changes!=['own_signed_source_remainder_input']:
        raise ValueError('Unreviewed modified full physical lift adaptation')
    env=dict(vars(original))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
      '<original full physical lift; actual modified remainder sectors>','exec'),env)
    name=PREFIX+'current_modified_pre_physical_tensor_operator.py';asts.hashes[name]=sha(name)
    return env[fn.name],dict(exact_AST_adaptations=changes,input_hashes=asts.hashes,
      original_physical_operators_factor_logs_completed_radius_and_Cartesian_programs_unchanged=True,
      all_actual_modified_remainder_Pstar_sectors_retained=True)

def exact_modified_radius(field,region,q,source):
    c=field.ctx;geometry=field.geometry;pre=source['current_original_pre_source']
    if region=='quiet_O3_power':
        phase=q/c.mpf(geometry.params.Tw)
        legacy,radius=geometry.radius('O3_power',phase,pre,geometry.pre)
        # This is the exact Tw*(q/Tw)=q source correlation, before enclosure.
        logR=geometry.logRref+geometry.logP+1+q
        return logR,dict(original_radius_source=radius,original_O3_power_phase_enclosure=phase,
          original_uncorrelated_dispatch_logR_enclosure=legacy,
          exact_radius_function='logRref+logPstar+1+q; original phase=q/same_actual_Tw',
          same_source_Tw_phase_correlation_reduced_before_enclosure=True)
    logR,radius=geometry.radius(region,q,pre,geometry.pre)
    return logR,dict(original_radius_source=radius,same_original_source_radius_unchanged=True)

def modified_physical_packet(field,signed):
    c=field.ctx;region=signed['region'];q=signed['coordinate'];source=signed['actual_source']
    logR,radius=exact_modified_radius(field,region,q,source);sectors={}
    for label,parts in signed['full_signed_paper_theta_axial_stress3_sectors'].items():
        sectors[label]={}
        for name,part in parts.items():
            rp,bp,dp,hp=part['mode']
            if dp!=0 or hp!=0:raise ValueError('Foreign modified pre amplitude modes')
            grid={'s%d_Z%d'%(j,n):row[n]*math.factorial(n)
              for j,row in enumerate(part['full_derivative_rows']) for n in range(4-j)}
            sectors[label][name]=dict(mode=part['mode'],exact_source_log_parts=dict(source_logR=c.mpf(str(rp))*logR,
              logPstar=bp*field.geometry.logP,selected_log_end_scale=c.mpf(0),
              signed_original_memory_log=c.mpf(0),normalization=-c.ln(2)/2),
              full_stress_mixed3_coefficient_enclosures=grid)
    packet=dict(Z=signed['Z'],s=q,exact_logR=logR,
      exact_pulse_reference_logB_parts=dict(logPstar=field.geometry.logP),
      exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors,
      actual_modified_source_remainder_sectors=signed['full_signed_source_remainder2_sectors'],
      parent_modified_source_definition_sha256=signed['modified_source_definition_sha256'],
      parent_modified_stress_definition_sha256=signed['modified_stress_definition_sha256'])
    return packet,radius

def exact_modified_pre_physical_theorem(field):
    """Bind same geometry, source FTC, exact radius and full tensor operator."""
    from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly as BASE
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators,ZSYM,DSYM,BSYM
    from lei_ren_part1_paper_compliant_current_modified_pre_stress import CurrentModifiedPreStress
    if type(field.source) is not CurrentModifiedPreStress or not field.source.acceptance_loaded:
        raise ValueError('Checked actual modified signed source required')
    field.source.source.histories.candidate.direction.assert_graph()
    registry=field.source.source.histories.candidate.direction.registry
    live=dict(same_checked_original_geometry=field.geometry is registry.owners['o3'].physical,
      same_original_pre_source=field.geometry.pre is field.source.source.histories.pre,
      same_O2_pre_source=registry.owners['o2'].physical.pre is field.geometry.pre,
      same_actual_Tw_parameters=field.geometry.params is field.geometry.pre.params,
      same_actual_logPstar_value=endpoints(field.geometry.logP)==endpoints(field.ctx.mpf(endpoints(field.geometry.pre.params.logPstar))),
      original_radius_callable=field.geometry.radius.__func__ is BASE.radius,
      same_context=field.ctx is field.geometry.ctx)
    if not all(live.values()):raise ValueError('Modified source geometry/logPstar lineage changed')
    operator=registry.owners['o3'].proof['consumed_checked_original_full_physical_operator']
    if not all(operator['identities'].values()) or not operator['full_nonzero_meridional_velocity_retained']:
        raise ValueError('Original arbitrary-source full physical operator is required')
    asts=SourceAST();checks={}
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Modified physical source identity: '+name)
        checks[name]=True
    asts.expression('global_physical_assembly','__init__','self.params',wanted='self.pre.params')
    asts.expression('global_physical_assembly','__init__','self.logP',wanted='c.mpf(endpoints(self.params.logPstar))')
    checks['actual_geometry_Tw_is_same_original_pre_parameter_object']=True
    checks['actual_geometry_logPstar_is_same_original_pre_scale_value']=True
    asts.expression('current_modified_pre_physical_tensor','__init__','self.source',
      wanted='source if source is not None else CurrentModifiedPreStress()')
    asts.expression('current_modified_pre_physical_tensor','tensor','signed',wanted='self.source.stress(region,Z,coordinate)')
    asts.expression('current_modified_pre_physical_tensor','tensor','(packet, radius)',wanted='modified_physical_packet(self,signed)')
    asts.expression('current_modified_pre_physical_tensor','tensor','point',wanted='self.lift(c,packet,self.delta,None,lt,theta,nu)')
    checks['actual_same_checked_signed_columns_and_remainder_feed_full_physical_program']=True
    rr,lp,Tw,q=s.symbols('same_logRref same_logPstar same_Tw q',real=True)
    ctx=SimpleNamespace(mpf=s.sympify)
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=()))
    geometry=SimpleNamespace(ctx=ctx,logRref=rr,logP=lp,params=SimpleNamespace(Tw=Tw))
    geometry.pre=object();geometry.radius=lambda chart,value,packet,provider:radius(geometry,chart,value,packet,provider)
    new_radius=asts.replay('current_modified_pre_physical_tensor_operator','exact_modified_radius',{})
    geomfield=SimpleNamespace(ctx=ctx,geometry=geometry)
    actual,_=new_radius(geomfield,'quiet_O3_power',q,dict(current_original_pre_source={}))
    zero('actual_quiet_log_radius_matches_original_q_over_same_Tw_phase',actual,geometry.radius('O3_power',q/Tw,{},None)[0])
    for region,pre in (('O2_buffer',dict(actual_y=lp-11+q,exact_radius_source='same actual y')),
                       ('O3_slope_mu',{})):
        actual,_=new_radius(geomfield,region,q,dict(current_original_pre_source=pre))
        zero('actual_'+region+'_original_radius_retained',actual,rr+lp+q-(11 if region=='O2_buffer' else 0))
    # Normalized delta m is a genuine cumulative moment, with m_y=Vhat-m.
    # Prove physical divergence from the actual source Q and physical H01.
    y,z,delta=s.symbols('logR Z delta',real=True);L=1-delta*z*z
    m=s.Function('same_actual_delta_m')(y,z);V=s.Function('same_actual_Vhat')(y,z)
    node=asts.method('current_O3_modulated_histories_operator','increment_rows')
    expected=ast.dump(ast.parse("d['m'].append(Vnew[j]-d['m'][j])",mode='eval').body)
    if not any(isinstance(n,ast.Call) and ast.dump(n)==expected for n in ast.walk(node)):
        raise ValueError('Actual own M FTC derivative row changed')
    mr=[m,V-m]+[s.Integer(0)]*3;vr=[V,s.diff(V,y)]+[s.Integer(0)]*3
    Q=asts.evaluate(asts.expression('current_O3_modulated_histories_operator','increment_rows','Q'),
      dict(z=z,Vnew=vr,d=dict(m=mr),delta=delta,L=L,square=lambda x:x*x,derivative=lambda x:s.diff(x,z)))
    H01=sum(coefficient.subs({ZSYM:z,DSYM:delta,BSYM:-1-delta})*s.diff(V,y,k,z,n)
      for (k,n),coefficient in physical_operators()[0,1].items())
    zero('actual_own_moment_FTC_Q_makes_modified_physical_velocity_divergence_zero',Q[1]+Q[0]+H01,0)
    checks['actual_original_and_delta_divergence_sum_zero_by_checked_original_operator']=True
    # Taylor coefficient conversion of the actual packet is source-bound.
    # Replay the callable assembler with arbitrary independent sector jets.
    class Row:
        def __init__(self,j):self.j=j
        def __getitem__(self,n):return s.Symbol('same_source_y'+str(self.j)+'_Z'+str(n))/s.factorial(n)
    packetctx=SimpleNamespace(mpf=lambda x:s.Rational(str(x)),ln=s.log)
    abstract=dict(region='quiet_O3_power',coordinate=q,Z=z,actual_source={},
      modified_source_definition_sha256='same source',modified_stress_definition_sha256='same columns',
      full_signed_paper_theta_axial_stress3_sectors=dict(theta=dict(source=dict(mode=(.5,1,0,0),full_derivative_rows=[Row(j) for j in range(4)])),
        axial=dict(source=dict(mode=(-.5,2,0,0),full_derivative_rows=[Row(j) for j in range(4)]))),
      full_signed_source_remainder2_sectors={'same':'same actual signed remainder object'})
    packfn=asts.replay('current_modified_pre_physical_tensor_operator','modified_physical_packet',
      dict(math=math,exact_modified_radius=lambda *args:(rr,'same radius')))
    packed,_=packfn(SimpleNamespace(ctx=packetctx,geometry=SimpleNamespace(logP=lp)),abstract)
    for label,part in packed['full_meridional_stress_log_sectors'].items():
        grid=part['source']['full_stress_mixed3_coefficient_enclosures']
        for j in range(4):
            for n in range(4-j):zero('actual_'+label+'_ordinary_mixed_source_grid_'+str(j)+'_'+str(n),grid['s%d_Z%d'%(j,n)],s.Symbol('same_source_y'+str(j)+'_Z'+str(n)))
        mode=abstract['full_signed_paper_theta_axial_stress3_sectors'][label]['source']['mode']
        parts=part['source']['exact_source_log_parts']
        zero('actual_'+label+'_packet_radius_factor',parts['source_logR'],s.Rational(str(mode[0]))*rr)
        zero('actual_'+label+'_packet_Pstar_factor',parts['logPstar'],mode[1]*lp)
        zero('actual_'+label+'_packet_half_normalization',parts['normalization'],-s.log(2)/2)
    if packed['actual_modified_source_remainder_sectors'] is not abstract['full_signed_source_remainder2_sectors']:
        raise ValueError('Own source remainder replaced in physical packet')
    checks['actual_remainder_object_and_modes_pass_unchanged_to_original_factor_logs']=True
    return dict(identities=checks,passed=True,input_hashes={**operator['input_hashes'],**asts.hashes},
      live_original_geometry_and_checked_source_bindings=live,
      consumed_original_arbitrary_source_full_physical_operator=operator,
      completed_diagonal='Ttheta_theta=r*d_z(Trz); exact log completed-radius=logR/2+log(2)/2',
      own_radial_moment_source_and_physical_incompressibility_bound=True,
      only_remainder_input_changes_in_full_original_lift=True,
      common_cone_N_modified_interfaces_energy_recursion_and_full_NS_remain_open=True)
