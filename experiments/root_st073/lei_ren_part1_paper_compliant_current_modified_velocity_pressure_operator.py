"""Own four-label mixed4 source packet and unchanged physical pullbacks.

Ordinary radius derivatives already include the radial half-power. Units
are fixed at the query basepoint and must not be differentiated again.
"""
import ast
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_modified_pre_physical_tensor_operator import exact_modified_radius
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    cartesian_source_row,time_source_row,log_row,zero_powers,COMPONENTS,BASES)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P,physical_operators
from lei_ren_part1_paper_compliant_cartesian_field import INDICES
from lei_ren_part1_paper_compliant_current_modified_heat_inheritance import sha,PREFIX,endpoints

LABELS={'theta':UT,'axial':UZ,'radial':UR}
_original_velocity_proof=None

def original_velocity_proof():
    """Reuse an original theorem only while every consumed source is unchanged."""
    from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_velocity_remainder_theorem
    global _original_velocity_proof
    if _original_velocity_proof is None or any(sha(name)!=digest
      for name,digest in _original_velocity_proof['input_hashes'].items()):
        _original_velocity_proof=raw_pre_velocity_remainder_theorem()
    return _original_velocity_proof

def modified_velocity_pressure_packet(field,view):
    c=field.ctx;grids={};inventory={};geometry=field.geometry
    for name,label in LABELS.items():
        parts=view['modified_cylindrical_velocity_source_log_sectors'][name]
        grids[label]={(j,n):[] for j in range(5) for n in range(5-j)}
        inventory[label]=[]
        for part in parts:
            power=part['Pstar_power'];rows=part['ordinary_logR_rows']
            if power not in (0,1) or len(rows)!=5:raise ValueError('Full actual own velocity source sectors required')
            if name!='radial' and power!=1:raise ValueError('Theta and modified axial source require Pstar1')
            powers=list(zero_powers());powers[1]=mp.mpf(power)
            inventory[label].append(power)
            for j in range(5):
                if rows[j].order<4-j:raise ValueError('Insufficient actual axial mixed4 source row')
                for n in range(5-j):
                    grids[label][j,n].append((tuple(powers),rows[j][n]*math.factorial(n)))
    if inventory!={UT:[1],UZ:[1],UR:[0,1]}:raise ValueError('Own signed radial Pstar0/1 layout changed')
    pressure=view['modified_absolute_pressure_over_Pstar2_ordinary_logR_rows']
    if len(pressure)!=5 or any(row.order<4-j for j,row in enumerate(pressure)):
        raise ValueError('Own absolute P0+Cp mixed4 pressure source required')
    powers=list(zero_powers());powers[1]=mp.mpf(2)
    grids[P]={(j,n):[(tuple(powers),row[n]*math.factorial(n))]
      for j,row in enumerate(pressure) for n in range(5-j)}
    inventory[P]=[2]
    logR,radius=exact_modified_radius(field,view['region'],view['coordinate'],view)
    logs=(c.mpf(0),geometry.logP,c.mpf(0),c.mpf(0),c.mpf(0),logR,c.ln(2))
    amplitudes={UT:{},UZ:{},UR:{5:mp.mpf('.5'),6:mp.mpf('-.5')},P:{}}
    return dict(grids=grids,source_log_bases=logs,amplitudes=amplitudes,
      signed_Pstar_sector_inventory=inventory,exact_logR=logR,
      exact_original_geometry_radius_source=radius,
      source_orders=dict(ordinary_logR=4,axial_Z=5,total_mixed=4),
      actual_pressure_is_original_P0_plus_own_Cp=True,
      fixed_basepoint_units_not_redifferentiated=True,
      radial_sqrt_R_over_2_rows_already_include_half_power_derivatives=True)

def modified_physical_velocity_pressure(field,packet,Z,log_tau,theta):
    c=field.ctx;delta=field.delta;lt=c.mpf(log_tau);Z=c.mpf(Z)
    cosine=sine=c.mpf([-1,1]) if theta is None else None
    if theta is not None:cosine=c.cos(c.mpf(theta));sine=c.sin(c.mpf(theta))
    grids=packet['grids'];logs=packet['source_log_bases'];amplitudes=packet['amplitudes']
    spatial={}
    for i,j,b in INDICES:
        key='x%d_y%d_z%d'%(i,j,b);spatial[key]={}
        for component in COMPONENTS:
            parts=cartesian_source_row(c,grids,component,i,j,b,Z,delta,cosine,sine,amplitudes)
            spatial[key][component]={label:log_row(c,rows,logs,gamma,lt/2)
              for label,(rows,gamma) in parts.items()}
    time={};basis={'ux':{UR:cosine,UT:-sine},'uy':{UR:sine,UT:cosine},'uz':{UZ:c.mpf(1)},'p':{P:c.mpf(1)}}
    for component,labels in basis.items():
        time[component]={}
        for label,angular in labels.items():
            rows,gamma=time_source_row(c,grids,label,Z,delta,amplitudes)
            time[component][label]=log_row(c,[(powers,value*angular) for powers,value in rows],logs,gamma,lt/2)
    coordinates=dict(exact_log_r='log(lambda)+log(2R)/2',exact_z='Z*lambda^(1-delta)',
      lambda_relation='lambda^2*(1-Z^2)=tau',exact_logR=packet['exact_logR'])
    lo,hi=endpoints(Z)
    if -1<lo and hi<1:
        loglambda=(lt-c.ln(1-Z**2))/2
        coordinates.update(log_lambda_enclosure=loglambda,
          log_r_enclosure=loglambda+(packet['exact_logR']+c.ln(2))/2)
    return dict(physical_spatial_cartesian_mixed4=spatial,
      first_fixed_x_physical_time_derivative=time,physical_coordinates=coordinates,
      requested_log_tau=lt,cos_theta=cosine,sin_theta=sine,
      source_log_base_names=BASES,source_log_bases=logs,
      actual_signed_source_factors_combined_before_final_bound=True,
      positive_source_exponentials_not_materialized=True,
      output_kind='signed factored bounds on actual modified source functions',
      resolved_point_coefficient_values_available=False)

def modified_absolute_pressure_source_theorem(field):
    """Replay P0+old Cp+own Cp and their actual FTC mixed4 recurrence."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    from lei_ren_part1_paper_compliant_current_O3_modulated_histories_operator import SHAPES
    analytic=field.heat.theorem['consumed_original_analytic_preheat_native_pressure_function_theorem']
    if not analytic['passed'] or not analytic['original_P0_not_redefined_or_pressure_patched'] or not analytic['holomorphic_integral_and_axial5_differentiation_justified']:
        raise ValueError('Same actual analytic P0 defining function required')
    if field.source.histories.pre.datum is not field.heat.datum:
        raise ValueError('Own pressure source lost original analytic P0 object')
    required=((field.source.histories.theorem,'actual_continuous_transport_FTC_p'),
      (field.source.histories.theorem,'actual_p_increment_shared_axis_pressure_units'),
      (field.source.theorem,'actual_repair_pressure_FTC'))
    if any(proof['identities'].get(key) is not True for proof,key in required):
        raise ValueError('Actual own Cp continuous/modulation/repair FTC required')
    asts=SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(a-b)!=0:raise ArithmeticError('Actual absolute pressure mixed4 identity: '+name)
        checks[name]=True
    for stem,method,target,wanted in (
      ('current_O3_modulated_histories','history','datum',"IntervalTaylor(c,[c.mpf(endpoints(x)) for x in pre['original_P0_axial5_coefficients']])"),
      ('current_O3_modulated_histories','history','original_P',"[old['p'][0]+datum]+old['p'][1:]"),
      ('current_O3_modulated_histories','history','pressure',"[a+b for a,b in zip(original_P,delta_rows['p'])]"),
      ('current_O3_repaired_histories','history','original_pressure',"parent['original_absolute_pressure_over_Pstar2_ordinary_logR_rows']"),
      ('current_O3_repaired_histories','history','(delta_rows, Q)',"increment_rows(c,z,self.histories.delta,self.histories.Ua,scalar,Unew,original['theta'],Vnew,du)"),
      ('current_O3_modulated_histories_operator','increment_rows','ds',"product_rows(Uold,swirl_increment)[j]*2+product_rows(swirl_increment,swirl_increment)[j]")):
        asts.expression(stem,method,target,wanted=wanted)
        checks['actual_'+stem+'_'+target+'_pressure_source_bound']=True
    repair=asts.method('current_O3_repaired_histories','history')
    expected=ast.parse("[a+b for a,b in zip(original_pressure,delta_rows['p'])]",mode='eval').body
    if not any(isinstance(node,ast.keyword) and node.arg=='modified_absolute_pressure_over_Pstar2_ordinary_logR_rows'
      and ast.dump(node.value)==ast.dump(expected) for node in ast.walk(repair)):
        raise ValueError('Actual repaired absolute pressure publication changed')
    original=asts.method('pre_pulse_mixed_C4','physical_mixed')
    physical=asts.expression('pre_pulse_mixed_C4','physical_mixed','physical')
    value=next(item.value for item in physical.keywords if item.arg=='P_over_Pstar2')
    if ast.dump(value)!=ast.dump(ast.parse("[p0+rows['p'][0]]+rows['p'][1:]",mode='eval').body):
        raise ValueError('Original absolute P0+Cp publication changed')
    append=next(node for node in ast.walk(original) if isinstance(node,ast.Call)
      and ast.unparse(node.func)=="rows['p'].append")
    if ast.dump(append.args[0])!=ast.dump(ast.parse('product_rows(U,U,j)/2',mode='eval').body):
        raise ValueError('Actual original pressure FTC recurrence changed')
    increment=asts.method('current_O3_modulated_histories_operator','increment_rows')
    append_delta=next(node for node in ast.walk(increment) if isinstance(node,ast.Call)
      and ast.unparse(node.func)=="d['p'].append")
    if ast.dump(append_delta.args[0])!=ast.dump(ast.parse('ds/2',mode='eval').body):
        raise ValueError('Actual own pressure FTC recurrence changed')
    checks['actual_original_and_own_pressure_FTC_recurrence_publications_bound']=True
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)))
    jet=lambda v,order=5:Projection.function(c,v,order)
    y=s.Symbol('same_logR',real=True);base=s.Symbol('same_Ua',real=True)
    u=s.Function('same_Uold')(y,Z);du=s.Function('same_swirl_increment')(y,Z)
    old=s.Function('same_original_Cp')(y,Z);own=s.Function('same_own_Cp')(y,Z)
    p0=s.Function('same_analytic_P0')(Z);C=1/(1+Z**2)
    oldrows=[jet(s.diff(u,y,j)) for j in range(5)]
    durows=[jet(s.diff(du,y,j)) for j in range(5)]
    newrows=[a+b for a,b in zip(oldrows,durows)]
    increments={key:s.Symbol('same_'+key) for key in SHAPES};increments['p']=own/(C**2*base**2)
    env=dict(SHAPES=SHAPES,square=lambda v:v*v,product_rows=product_rows,
      derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1))
    fn=asts.replay('current_O3_modulated_histories_operator','increment_rows',env)
    delta_rows,_=fn(c,jet(Z),s.Symbol('delta'),base,increments,newrows,oldrows,
      [jet(s.Function('same_Vhat')(y,Z))]*5,durows)
    original_product=asts.replay('pre_pulse_mixed_C4','product_rows',dict(math=math))
    old_pressure=[jet(old)]+[asts.evaluate(append.args[0],dict(U=oldrows,j=j,product_rows=original_product)) for j in range(4)]
    original_published=asts.evaluate(value,dict(p0=jet(p0),rows=dict(p=old_pressure)))
    published=asts.evaluate(expected,dict(original_pressure=original_published,delta_rows=delta_rows))
    for j in range(5):
        target=p0+old+own if j==0 else s.diff((u+du)**2/2,y,j-1)
        for n in range(5-j):
            zero('actual_absolute_P0_plus_own_Cp_y%d_Z%d'%(j,n),
              published[j][n]*math.factorial(n),s.diff(target,Z,n))
    asts.method('current_pressure_terminal_closure','original_pressure_function_identification')
    checks['actual_P0_holomorphic_derivatives_and_own_partial_FTC_consumed']=True
    checks['P0_only_in_zero_logR_row_with_all_actual_axial_derivatives']=True
    checks['physical_absolute_pressure_has_unchanged_Pstar_squared_units']=True
    return dict(identities=checks,passed=True,input_hashes=asts.hashes,
      consumed_original_analytic_P0_function=True,
      consumed_own_pressure_FTC_keys=[key for _,key in required],
      ordinary_pressure_rows='j=0: P0+old_Cp+own_Cp; j>=1: Dy^(j-1)(Unew^2/2)',
      axial_coefficient_conversion='ordinary Dy^j DZ^n = n!*Z_Taylor_coefficient; j+n<=4',
      same_actual_modulation_and_repair_cumulative_pressure_no_new_datum=True)

def exact_modified_velocity_pressure_theorem(field):
    """Bind actual source programs, ordinary projections and original maps."""
    field.heat.assert_graph()
    if field.source is not field.heat.interfaces.repaired or field.geometry is not field.heat.dispatch.local.geometry:
        raise ValueError('Foreign repaired velocity or physical geometry')
    if field.source.histories.pre is not field.geometry.pre or field.ctx is not field.geometry.ctx:
        raise ValueError('Own velocity and original radius coordinate source differ')
    local=field.heat.dispatch.local.theorem
    if not local['passed'] or not all(local['identities'].values()):
        raise ValueError('Checked original radius/coordinate/unit source theorem required')
    if not field.source.acceptance_loaded or not field.heat.acceptance_loaded:
        raise ValueError('Checked actual repaired history and analytic pressure lineage required')
    raw=original_velocity_proof()
    if not raw['passed'] or not all(raw['identities'].values()):raise ValueError('Original variable mixed4 velocity recovery required')
    for proof in (field.source.theorem,field.source.histories.theorem):
        if not proof['passed'] or not all(proof['identities'].values()):raise ValueError('Own profile and history FTC theorem required')
    pressure=modified_absolute_pressure_source_theorem(field)
    asts=SourceAST();checks=dict(pressure['identities'])
    def zero(name,a,b):
        if s.cancel(a-b)!=0:raise ArithmeticError('Modified physical source identity: '+name)
        checks[name]=True
    for target,wanted in (
      ('parts',"view['modified_cylindrical_velocity_source_log_sectors'][name]"),
      ('pressure',"view['modified_absolute_pressure_over_Pstar2_ordinary_logR_rows']"),
      ('(logR, radius)',"exact_modified_radius(field,view['region'],view['coordinate'],view)")):
        asts.expression('current_modified_velocity_pressure_operator','modified_velocity_pressure_packet',target,wanted=wanted)
    class Row:
        def __init__(self,label,j):self.label=label;self.j=j;self.order=5
        def __getitem__(self,n):return s.Symbol(self.label+'_y%d_Z%d'%(self.j,n))/s.factorial(n)
    part=lambda label,power:dict(Pstar_power=power,ordinary_logR_rows=[Row(label,j) for j in range(5)])
    view=dict(region='quiet_O3_power',coordinate=s.Symbol('q'),
      modified_cylindrical_velocity_source_log_sectors={
        'theta':[part('Ut',1)],'axial':[part('Uz',1)],'radial':[part('Ur0',0),part('Ur1',1)]},
      modified_absolute_pressure_over_Pstar2_ordinary_logR_rows=[Row('Pabs',j) for j in range(5)])
    rr,lp=s.symbols('same_logR same_logPstar',real=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),ln=s.log)
    fn=asts.replay('current_modified_velocity_pressure_operator','modified_velocity_pressure_packet',
      dict(mp=mp,math=math,LABELS=LABELS,UT=UT,UZ=UZ,UR=UR,P=P,zero_powers=zero_powers,
        exact_modified_radius=lambda *args:(rr,'same actual radius')))
    actual=fn(SimpleNamespace(ctx=c,geometry=SimpleNamespace(logP=lp)),view)
    for label,seeds in ((UT,('Ut',)),(UZ,('Uz',)),(UR,('Ur0','Ur1')),(P,('Pabs',))):
        for j in range(5):
            for n in range(5-j):
                terms=actual['grids'][label][j,n]
                if len(terms)!=len(seeds):raise ValueError('Signed source sectors were merged before projection')
                for (powers,value),seed in zip(terms,seeds):
                    zero('actual_ordinary_mixed4_'+seed+'_y%d_Z%d'%(j,n),value,s.Symbol(seed+'_y%d_Z%d'%(j,n)))
    zero('actual_Pstar_log_scale_is_original_source',actual['source_log_bases'][1],lp)
    zero('actual_radius_log_scale_is_original_source',actual['source_log_bases'][5],rr)
    zero('actual_radial_fixed_basepoint_half_unit',s.Rational(str(actual['amplitudes'][UR][5]))*rr+
      s.Rational(str(actual['amplitudes'][UR][6]))*s.log(2),(rr-s.log(2))/2)
    for label,parts in actual['signed_Pstar_sector_inventory'].items():
        for i,power in enumerate(parts):
            zero('actual_'+label+'_Pstar_power_'+str(i),s.Rational(str(actual['grids'][label][0,0][i][0][1])),power)
    # Both original consumers take the ordinary mixed rows without any
    # further phase conversion or moving normalization differentiation.
    for target,wanted in (
      ('parts','cartesian_source_row(c,grids,component,i,j,b,Z,delta,cosine,sine,amplitudes)'),
      ('(rows, gamma)','time_source_row(c,grids,label,Z,delta,amplitudes)')):
        asts.expression('current_modified_velocity_pressure_operator','modified_physical_velocity_pressure',target,wanted=wanted)
    for stem,name,live in (('global_physical_assembly','cartesian_source_row',cartesian_source_row),
      ('global_physical_assembly','time_source_row',time_source_row),('global_physical_assembly','log_row',log_row),
      ('pulse_physical_bounds','physical_operators',physical_operators)):
        asts.method(stem,name)
        if live.__module__!=PREFIX+stem:raise ValueError('Original physical callable replaced: '+name)
        checks['unchanged_actual_physical_callable_'+name]=True
    z,d,beta,g0,gz,gy=s.symbols('Z delta beta same_G same_GZ same_Gy',real=True)
    ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)))
    tf=asts.replay('global_physical_assembly','time_source_row',dict(UR=UR,UT=UT,UZ=UZ,P=P,shift=lambda p,a:p))
    grid={label:{(0,0):[((),g0)],(0,1):[((),gz)],(1,0):[((),gy)]} for label in (UR,UT,UZ,P)}
    for label,bb in ((UR,-1),(UT,-1-d),(UZ,-1-d),(P,-2-2*d)):
        rows,gamma=tf(ctx,grid,label,z,d,dict.fromkeys(grid,{}))
        zero('actual_fixed_x_time_operator_'+label,sum(value for _,value in rows),
          (-bb*g0+(1-d)*z*gz+2*gy)/(2*(1-d*z*z)))
        zero('actual_fixed_x_time_lambda_exponent_'+label,gamma,bb-2)
    lam=s.Symbol('same_lambda',positive=True);L=1-d*z*z
    lamdot=-1/(2*lam*L);zdot=(1-d)*z/(2*lam**2*L);ydot=1/(lam**2*L)
    zero('fixed_x_time_differentiates_actual_tau_relation',
      2*lam*(1-z*z)*lamdot-2*lam**2*z*zdot,-1)
    zero('fixed_x_time_preserves_actual_physical_axial_coordinate',
      zdot+(1-d)*z*lamdot/lam,0)
    zero('fixed_x_time_preserves_actual_physical_radial_coordinate',lamdot/lam+ydot/2,0)
    zero('fixed_x_time_operator_from_actual_implicit_coordinates',
      beta*g0*lamdot/lam+gz*zdot+gy*ydot,
      (-beta*g0+(1-d)*z*gz+2*gy)/(2*lam**2*L))
    dependencies={index for row in physical_operators().values() for index in row}
    if any(j+n>4 for j,n in dependencies):raise ValueError('Original spatial4 consumes unavailable source order')
    checks['all_actual_original_spatial4_operator_dependencies_in_own_packet']=True
    for target,wanted in (('view','self.source.history(region,Z,coordinate)'),
      ('packet','modified_velocity_pressure_packet(self,view)'),
      ('physical','modified_physical_velocity_pressure(self,packet,view[\'Z\'],lt,theta)')):
        asts.expression('current_modified_velocity_pressure','physical',target,wanted=wanted)
    checks['runtime_own_repaired_history_packet_and_original_consumers_bound']=True
    return dict(identities=checks,passed=True,input_hashes={**raw['input_hashes'],**pressure['input_hashes'],**asts.hashes},
      consumed_original_variable_velocity_mixed4_theorem=raw,
      actual_absolute_pressure_mixed4_source_theorem=pressure,
      exact_source_dependency_indices=sorted(dependencies),
      same_original_geometry_and_checked_repaired_history=True,
      checked_original_radius_and_unit_source_theorem_consumed=True,
      only_local_modified_velocity_pressure_packet_admitted=True,
      quantified_affected_interfaces_total_energy_cones_recursion_global_NS_remain_open=True)
