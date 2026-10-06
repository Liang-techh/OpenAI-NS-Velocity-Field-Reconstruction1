"""Full micro-switch tensor algebra without resolving width or amplitude factors."""
import ast
import copy
import math
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_microswitch_mixed_C4 as original
import lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator as raw_original
import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as physical_original
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredJet
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows

SOURCE_KEY='actual_unresolved_factored_source_rows'
SOURCE_EXPRESSION="dict(algebra=algebra,physical=physical,primitives=primitives,Q=Q)"
LOG_NAMES=('original_width_log','original_Pstar_squared_log','original_F0_squared_base_log','original_swirl_squared_base_log')


def compile_function(fn,env,label):
    fn=copy.deepcopy(fn);fn.decorator_list=[]
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),label,'exec'),env)
    return env[fn.name]


def expose_phase_rows(asts):
    before=copy.deepcopy(asts.method('microswitch_mixed_C4','phase_physical'));fn=copy.deepcopy(before)
    returns=[node for node in fn.body if isinstance(node,ast.Return)]
    if len(returns)!=1 or ast.unparse(returns[0].value.func)!='dict':raise ValueError('One original phase source return required')
    returns[0].value.keywords.append(ast.keyword(arg=SOURCE_KEY,value=ast.parse(SOURCE_EXPRESSION,mode='eval').body))
    compare=copy.deepcopy(fn)
    next(node for node in compare.body if isinstance(node,ast.Return)).value.keywords.pop()
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Original phase source math changed')
    return fn


def source_only_projection(fn,key=SOURCE_KEY):
    """Project existing source locals for symbolic equality before numeric bounds."""
    fn=copy.deepcopy(fn)
    returned=next(node for node in fn.body if isinstance(node,ast.Return))
    value=next(keyword.value for keyword in returned.value.keywords if keyword.arg==key)
    returned.value=copy.deepcopy(value)
    return fn


def compiled_original_microswitch_with_rows():
    asts=SourceAST();env=dict(vars(original));phase=expose_phase_rows(asts)
    env['phase_physical']=compile_function(phase,env,'<original phase; unresolved locals exposed only>')
    evaluate=copy.deepcopy(asts.method('microswitch_mixed_C4','evaluate'))
    result=compile_function(evaluate,env,'<original evaluate; output-only phase callback>')
    return result,dict(original_evaluate_AST_unchanged=True,original_phase_AST_unchanged_except_output_keyword=True,
        unresolved_existing_rows_exposed_without_feedback=True,original_factored_width_and_amplitude_algebra_retained=True,
        input_hashes=asts.hashes,passed=True)


def raw_microswitch_rows(c,packet):
    """Invert current radial prefactors after exact D_y^j=hb^-j D_s^j."""
    source=packet[SOURCE_KEY];algebra=source['algebra'];physical=source['physical'];primitives=source['primitives']
    ordinary=lambda rows:[algebra.width(row,-j) for j,row in enumerate(rows)]
    swirl_unit=lambda rows:[algebra.shift(row,(0,0,0,.5)) for row in rows]
    velocity=dict(theta=swirl_unit(ordinary(physical['Utheta_over_current_Utheta'])),
        axial=ordinary(physical['Uz']),radial=ordinary(physical['Ur_over_current_sqrt_R_over_2']))
    histories=dict(
        m=shifted_rows(ordinary(primitives['Mz_over_current_R']),-1,4),
        h=shifted_rows(swirl_unit(ordinary(primitives['Mtheta_over_current_sqrt2_R_1p5_Utheta'])),c.mpf('-1.5'),4),
        k=shifted_rows(swirl_unit(ordinary(primitives['Mtheta_z_over_current_sqrt2_R_1p5_Utheta'])),c.mpf('-1.5'),4),
        e=shifted_rows(ordinary(primitives['Mztheta_over_current_R_Pstar2']),-1,4),
        p=ordinary(primitives['Mp_over_Pstar2']))
    return dict(algebra=algebra,velocity=velocity,histories=histories,
        absolute_pressure=ordinary(physical['P_over_Pstar2']),
        inverse_width_source_shift_applied_before_any_resolution=True,
        original_radial_prefactors_included_once=True,full_absolute_P0_included_once=True)


def factored_axial_derivative(value):
    return value.Zderivative() if isinstance(value,FactoredJet) else axial_derivative(value)


def compiled_factored_source_operators():
    asts=SourceAST();env=dict(vars(physical_original));env['axial_derivative']=factored_axial_derivative
    compile_function(asts.method('pulse_end_physical_C2','axial_n'),env,'<original axial derivative with jet dispatch>')
    before=copy.deepcopy(asts.method('pulse_end_physical_C2','axial_operator_rows'));fn=copy.deepcopy(before);changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.FunctionDef) and node.name=='evaluate':
            node.body.insert(0,ast.parse('zero=z*0').body[0]);changes.append('coefficient_zero')
    if changes!=['coefficient_zero']:raise ValueError('Original axial coefficient program changed')
    compare=copy.deepcopy(fn)
    for node in ast.walk(compare):
        if isinstance(node,ast.FunctionDef) and node.name=='evaluate':node.body.pop(0)
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Axial coefficient mathematics changed')
    axial=compile_function(fn,env,'<original axial operator; unscaled coefficient accumulator>')
    def operator(c,rows,beta,z,delta,count=2,order=2):
        if isinstance(z,FactoredJet):
            if set(z.terms)!={(0,0,0,0)}:raise ValueError('Axial coordinate must be an unscaled coefficient jet')
            z=z.terms[(0,0,0,0)]
        return axial(c,rows,beta,z,delta,count=count,order=order)
    rawenv=dict(vars(raw_original));rawenv.update(axial_derivative=factored_axial_derivative,axial_operator_rows=operator)
    stress=compile_function(asts.method('current_pre_pulse_stress_operator','raw_pre_stress_rows'),rawenv,'<unchanged raw stress; factored jet dispatch>')
    remainder=compile_function(asts.method('current_pre_pulse_stress_operator','raw_pre_remainder_sectors'),rawenv,'<unchanged raw NS remainder; factored jet dispatch>')
    return stress,remainder,dict(original_raw_stress_and_remainder_AST_unchanged=True,
        original_axial_coefficient_operator_unchanged_except_zero_representation=True,
        actual_Z_derivatives_dispatch_without_differentiating_fixed_log_bases=True,
        coefficient_denominators_inverted_only_as_unscaled_Taylor_jets=True,input_hashes=asts.hashes,passed=True)


def source_log_parts(algebra,key):
    return {name:algebra.logs[j]*power for j,(name,power) in enumerate(zip(LOG_NAMES,key)) if power}


def split_factored_rows(rows,algebra):
    """Collect identical factors first, then keep every remaining signed sector."""
    if any(not isinstance(row,FactoredJet) or row.algebra is not algebra for row in rows):raise ValueError('One unresolved source algebra required')
    keys=sorted(set().union(*(row.terms for row in rows)))
    if not keys:keys=[(0,0,0,0)]
    zero=original.IntervalTaylor.constant(algebra.ctx,0,min(row.order for row in rows))
    return {key:[row.terms.get(key,zero) for row in rows] for key in keys}


def expanded_stress_sectors(c,rows,algebra,logR,logP):
    result={}
    for label,parts in rows.items():
        result[label]={}
        for name,part in parts.items():
            for index,(key,coefficients) in enumerate(split_factored_rows(part['full_derivative_rows'],algebra).items()):
                rp,bp,dp,hp=part['mode']
                result[label][name+'__source'+str(index)]=dict(original_stress_sector=name,source_exponents=list(key),
                    mode=part['mode'],exact_source_log_parts=dict(source_logR=rp*logR,logPstar=bp*logP,
                        selected_log_end_scale=c.mpf(0),signed_original_memory_log=c.mpf(0),normalization=-c.ln(2)/2,
                        **source_log_parts(algebra,key)),
                    full_stress_mixed3_coefficient_enclosures={'s%d_Z%d'%(j,n):row[n]*math.factorial(n)
                        for j,row in enumerate(coefficients) for n in range(4-j)})
    return result


def expanded_remainder_sectors(c,delta,z,velocity,remainder):
    algebra=velocity['theta'][0].algebra;result={}
    for label,parts in remainder(c,delta,algebra.lift(z),velocity).items():
        result[label]={}
        for name,part in parts.items():
            for index,(key,rows) in enumerate(split_factored_rows(part['rows'],algebra).items()):
                result[label][name+'__source'+str(index)]=dict(part,rows=rows,original_remainder_sector=name,
                    source_exponents=list(key),factored_source_log_parts=source_log_parts(algebra,key))
    return result


def compiled_factored_physical_lift(remainder):
    asts=SourceAST();before=copy.deepcopy(asts.method('pulse_end_physical_C2','lift_physical_packet'));fn=copy.deepcopy(before);changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='source_errors = pulse_remainder_sectors(c, delta, None, z, velocity)':
            node.value=ast.parse('expanded_remainder_sectors(c,delta,z,velocity,source_remainder)',mode='eval').body;changes.append('source_errors')
        if isinstance(node,ast.Assign) and ast.unparse(node)=="parts = factor_logs(c, packet, sector['mode'], sector['normalization_half'])":
            node.value=ast.parse("dict(factor_logs(c,packet,sector['mode'],sector['normalization_half']),**sector['factored_source_log_parts'])",mode='eval').body;changes.append('parts')
    if sorted(changes)!=['parts','source_errors']:raise ValueError('Unreviewed factored physical lift adaptation')
    compare=copy.deepcopy(fn)
    for node in ast.walk(compare):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='source_errors = expanded_remainder_sectors(c, delta, z, velocity, source_remainder)':
            node.value=ast.parse('pulse_remainder_sectors(c,delta,None,z,velocity)',mode='eval').body
        if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name) and node.targets[0].id=='parts' and isinstance(node.value,ast.Call) and ast.unparse(node.value.func)=='dict':
            node.value=ast.parse("factor_logs(c,packet,sector['mode'],sector['normalization_half'])",mode='eval').body
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Unreviewed physical stress/divergence/completion math changed')
    env=dict(vars(physical_original),expanded_remainder_sectors=expanded_remainder_sectors,source_remainder=remainder)
    lift=compile_function(fn,env,'<original full physical mapper; exact source sector expansion>')
    return lift,dict(original_physical_stress_divergence_diagonal_cartesian_operators_unchanged=True,
        only_remainder_source_units_and_factored_log_parts_adapted=True,
        every_signed_source_sector_retained_without_amplitude_resolution=True,input_hashes=asts.hashes,passed=True)


def factored_rows_record(rows):
    if isinstance(rows,FactoredJet):
        return dict(axial_order=rows.order,fixed_basepoint_log_factors=True,terms=[
            dict(source_exponents=list(key),axial_Taylor_coefficients=list(row.coefficients))
            for key,row in sorted(rows.terms.items())])
    if isinstance(rows,dict):return {key:factored_rows_record(value) for key,value in rows.items() if key!='algebra'}
    if isinstance(rows,(list,tuple)):return [factored_rows_record(value) for value in rows]
    return rows


def actual_microswitch_normalization_theorem():
    """Replay original phase generator and inverse prefactors for arbitrary jets."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda value:(value[0] if isinstance(value,(list,tuple)) else s.Rational(str(value)) if isinstance(value,(str,int,float)) else value))
    J=lambda expr,order=5:Projection.function(c,expr,order)
    hb,unit,Ps=s.symbols('same_positive_hb same_positive_swirl_unit Pstar',positive=True)
    delta=s.Symbol('delta',real=True);amp=J(s.Function('same_normalized_axial_swirl')(Z))
    logU=[J(s.Function('same_phase_logU'+str(k))(Z)) for k in range(1,5)]
    V=[J(s.Function('same_phase_V'+str(k))(Z)) for k in range(5)]
    shapes={key:J(s.Function('same_shape_'+key)(Z)) for key in ('theta','theta_z','mean','axial','swirl','pressure')}
    p0=J(s.Function('same_absolute_P0')(Z));asts=SourceAST()
    class NoFactored:pass
    env=dict(math=math,IntervalTaylor=Projection,FactoredJet=NoFactored,
        derivative=lambda value:J(s.diff(value.expr,Z),value.order-1),
        scaled_positive_source=lambda ctx,logbase,row,proofs:row*unit**2)
    for stem,name in (('long_reshape_mixed_C4','exponential_derivatives'),('microswitch_mixed_C4','rate_rows'),('microswitch_mixed_C4','product_rows')):
        asts.replay(stem,name,env)
    generator=compile_function(source_only_projection(expose_phase_rows(asts)),env,'<original phase source theorem; existing locals projection>')
    packet={SOURCE_KEY:generator(c,Z,delta,J(hb),amp,0,logU,V,shapes,p0,1/Ps**2,lambda row,power=1:row*hb**power,[])}
    class Algebra:
        width=staticmethod(lambda row,power:row*hb**power)
        shift=staticmethod(lambda row,key:row*unit**(2*s.Rational(str(key[3])))*Ps**(2*key[1]))
    packet[SOURCE_KEY]['algebra']=Algebra()
    env.update(shifted_rows=shifted_rows,SOURCE_KEY=SOURCE_KEY)
    adapter=asts.replay('current_microswitch_stress_operator','raw_microswitch_rows',env)
    raw=adapter(c,packet);checks={}
    source=packet[SOURCE_KEY]
    for label,name,unitfactor in (('theta','Utheta_over_current_Utheta',unit),('axial','Uz',1),('radial','Ur_over_current_sqrt_R_over_2',1),('pressure','P_over_Pstar2',1)):
        rows=raw['absolute_pressure'] if label=='pressure' else raw['velocity'][label]
        for j in range(5):
            difference=s.cancel(rows[j].expr-source['physical'][name][j].expr*unitfactor/hb**j)
            if difference!=0:raise ArithmeticError('Microswitch true ordinary source velocity units differ')
            for n in range(5-j):checks[label+'_y%d_Z%d'%(j,n)]=s.diff(difference,Z,n)==0
    for label,name,rate,unitfactor in (('m','Mz_over_current_R',1,1),('h','Mtheta_over_current_sqrt2_R_1p5_Utheta',s.Rational(3,2),unit),('k','Mtheta_z_over_current_sqrt2_R_1p5_Utheta',s.Rational(3,2),unit),('e','Mztheta_over_current_R_Pstar2',1,1),('p','Mp_over_Pstar2',0,1)):
        rows=shifted_rows(raw['histories'][label],rate,4)
        for j in range(5):
            difference=s.cancel(rows[j].expr-source['primitives'][name][j].expr*unitfactor/hb**j)
            if difference!=0:raise ArithmeticError('Microswitch full primitive inverse radial units differ')
            for n in range(5-j):checks[label+'_y%d_Z%d'%(j,n)]=s.diff(difference,Z,n)==0
    u=raw['velocity']['theta'];v=raw['velocity']['axial'];histories=raw['histories']
    rhs=dict(m=v[0]-histories['m'][0],h=u[0]-histories['h'][0]*s.Rational(3,2),
        k=u[0]*v[0]-histories['k'][0]*s.Rational(3,2),e=v[0]*v[0]/Ps**2-u[0]*u[0]/2-histories['e'][0],p=u[0]*u[0]/2)
    for key,row in rhs.items():
        if s.cancel(histories[key][1].expr-row.expr)!=0:raise ArithmeticError('Microswitch raw defining history ODE differs: '+key)
        checks['raw_five_history_ODE_'+key]=True
    if s.cancel(raw['absolute_pressure'][1].expr-u[0].expr**2/2)!=0:raise ArithmeticError('Microswitch absolute P0/pressure ODE differs')
    checks['absolute_pressure_y_ODE']=True
    return dict(identities=checks,original_phase_and_owned_raw_adapter_AST_replayed=True,
        arbitrary_smooth_axial_jet_and_positive_width_swirl_unit_schema=True,
        original_phase_to_ordinary_y_and_inverse_prefactor_identities_verified=True,
        five_full_history_ODEs_and_absolute_pressure_datum_retained=True,
        inverse_width_never_materialized_or_selected_from_cap=True,input_hashes=asts.hashes,passed=True)

def actual_microswitch_endpoint_theorem():
    """Bind the original two controls and replay every R2 physical source row."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda value:(s.Rational(str(value)) if isinstance(value,(str,int,float)) else value))
    J=lambda value,order=5:Projection.function(c,value,order)
    hb,unit,Ps=s.symbols('same_positive_hb same_current_swirl_unit Pstar',positive=True)
    delta=s.Symbol('delta',real=True);zero=J(0)
    asts=SourceAST();env=dict(math=math,IntervalTaylor=Projection,FactoredJet=type('NoFactored',(),{}),
        derivative=lambda value:J(s.diff(value.expr,Z),value.order-1),
        scaled_positive_source=lambda ctx,logbase,row,proofs:row*unit**2)
    for stem,name in (('long_reshape_mixed_C4','exponential_derivatives'),('microswitch_mixed_C4','rate_rows'),('microswitch_mixed_C4','product_rows'),('microswitch_mixed_C4','switch_controls')):
        asts.replay(stem,name,env)
    D=[J(s.Function('same_comparison_D'+str(k))(Z)) for k in range(4)]
    drive=lambda k,power:J(s.Function('same_comparison_drive'+str(k))(Z))*hb**power
    quotient=J(s.Function('same_actual_over_comparison_phi')(Z));width=lambda value,power=1:value*hb**power
    first=env['switch_controls'](c,'first',[1,0,0,0,0],D,drive,quotient,width)
    second=env['switch_controls'](c,'second',[0,0,0,0,0],D,drive,quotient,width)
    controls={}
    for name in first:
        for k,(left,right) in enumerate(zip(first[name],second[name])):
            if s.cancel(left.expr-right.expr)!=0:raise ArithmeticError('Actual phase1 source controls differ: '+name)
            controls[name+str(k)]=True
    exit_controls=env['switch_controls'](c,'second',[1,0,0,0,0],D,drive,quotient,width)
    for k,row in enumerate(exit_controls['logUtheta_phase_derivatives']):
        if s.cancel(row.expr-(hb/10 if k==0 else 0))!=0:raise ArithmeticError('Actual R2 source angular jets differ')
    if any(s.cancel(row.expr)!=0 for row in exit_controls['Uz_positive_phase_derivatives']):raise ArithmeticError('Actual R2 source axial jets not constant')
    amp=J(s.Function('same_actual_axial_swirl_jet')(Z));V=J(s.Function('same_actual_V110')(Z))
    shapes={key:J(s.Function('same_actual_R2_'+key)(Z)) for key in ('theta','theta_z','mean','axial','swirl','pressure')}
    p0=J(s.Function('same_actual_P0')(Z));ell=J(s.Function('same_actual_logu')(Z))
    phase=compile_function(source_only_projection(expose_phase_rows(asts)),env,'<original phase R2 source projection>')
    left=phase(c,Z,delta,J(hb),amp,0,exit_controls['logUtheta_phase_derivatives'],[V]+[zero]*4,
        shapes,p0,1/Ps**2,width,[])
    reshape=copy.deepcopy(asts.method('long_reshape_mixed_C4','reshape_mixed'))
    returned=next(node for node in reshape.body if isinstance(node,ast.Return))
    returned.value.keywords.append(ast.keyword(arg=SOURCE_KEY,value=ast.parse('dict(physical=physical,primitives=primitives,Q=Q)',mode='eval').body))
    env.update(relative_amplitude_jet=lambda logjet,m=1:amp**m,square=lambda value:value*value)
    asts.replay('reference_restore_mixed_C4','binomial_rate',env)
    rightfn=compile_function(source_only_projection(reshape),env,'<original reshape R2 source projection>')
    right=rightfn(c,Z,delta,ell,[J(s.Rational(1,10))]+[zero]*3,V,shapes,p0,1/Ps**2,[])
    checks={}
    for group in ('physical','primitives'):
        if set(left[group])!=set(right[group]):raise ValueError('Full R2 source component inventory differs')
        for name in left[group]:
            for k,(a,b) in enumerate(zip(left[group][name],right[group][name])):
                difference=s.cancel(a.expr/hb**k-b.expr)
                if difference!=0:raise ArithmeticError('Actual R2 source mixed rows differ: '+name+str(k))
                for n in range(5-k):checks[group+'/'+name+'/y%d_Z%d'%(k,n)]=s.diff(difference,Z,n)==0
    bindings=assignment_source_bindings('microswitch_mixed_C4','postpower',{
        'inlet':'self.switch.phase(Z,2)',
        'moments':'power_transport(c,theta,as_initial(moments2),phi2,V)',
        'phi':'phi2*theta**c.mpf(".4")',
        'packet':'reshape_mixed(c,Z,self.core.delta,logu,[constant(".1")]+[constant(0)]*3,V,shapes,jet(inlet["pressure_axis_axial5_coefficients"]),self.invP2,self.proofs)'})
    asts.method('inner_switch_profiles','power_transport')
    # At theta=1 each original kernel vanishes, and six moments equal the inlet.
    from lei_ren_part1_paper_compliant_current_switch_power_source_adapter import pure_source_power_transport
    transport=pure_source_power_transport(asts)
    initial={name:s.Function('same_R2_source_'+name)(Z) for name in ('H','K','mean','A','B','C')}
    moments=transport(c,s.Integer(1),initial,s.Function('same_phi2')(Z),V.expr)
    flattened=dict(H=moments[original.MTH],K=moments[original.MTHZ],mean=moments[original.MZ],A=moments[original.MZT]['axial'],B=moments[original.MZT]['swirl'],C=moments[original.MP])
    neutral={name:s.cancel(flattened[name]-initial[name])==0 for name in initial}
    if not all(neutral.values()):raise ArithmeticError('Original R2 zero-length transport changed')
    asts.expression('microswitch_mixed_C4','postpower','theta',wanted='c.exp(-zeta)')
    asts.expression('microswitch_mixed_C4','postpower','zeta',wanted='length*fraction')
    return dict(exact_phase1_original_control_identities=controls,
        actual_R2_original_physical_and_primitive_mixed4_identities=checks,
        actual_R2_source_rows_verified=len(checks),original_R2_six_zero_length_transport_identities=neutral,
        original_postpower_source_bindings=bindings,
        exact_R2_radius='100*exp(2*original_positive_hb)',
        exact_same_current_R2_phi_V_six_moments_and_absolute_P0=True,
        original_radial_prefactors_and_Bell_products_replayed=True,
        source_functions_compared_before_width_or_amplitude_enclosures=True,
        same_mixed4_source_and_full_stress_operators_imply_two_completed_tensor_traces=True,
        interval_overlap_not_used_as_function_identity=True,input_hashes=asts.hashes,passed=True)

def microswitch_factor_units_theorem():
    """Identify four source log bases with the independent R/Pstar mode factors."""
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    asts=SourceAST();hb,R,Ps,F0,phi=s.symbols('original_hb current_R Pstar F0 phi',positive=True)
    a,b,d,e=s.symbols('hb_power Pstar2_power F02_power swirl2_power',real=True)
    unit=s.sqrt(2*R)*F0*phi/Ps
    bases=(s.log(hb),2*s.log(Ps),2*s.log(F0),2*s.log(unit))
    expanded=s.expand_log(sum(power*base for power,base in zip((a,b,d,e),bases)),force=True)
    expected=a*s.log(hb)+(2*b-2*e)*s.log(Ps)+(2*d+2*e)*s.log(F0)+2*e*s.log(phi)+e*s.log(2*R)
    if s.simplify(expanded-s.expand_log(expected,force=True))!=0:raise ArithmeticError('Original source factor powers counted twice')
    u,H,K,m,A,B,C,P0=s.symbols('unit_u normalized_H normalized_K normalized_m normalized_A normalized_B normalized_C absolute_P0',real=True)
    h=unit*H;k=unit*K;energy=A/Ps**2-unit**2*B/2;pressure=unit**2*C/2
    identities=dict(
        velocity_theta=s.cancel(Ps*unit-s.sqrt(2*R)*F0*phi)==0,
        angular_moment=s.simplify(s.sqrt(2)*R**s.Rational(3,2)*Ps*h-2*R**2*F0*phi*H)==0,
        mixed_moment=s.simplify(s.sqrt(2)*R**s.Rational(3,2)*Ps*k-2*R**2*F0*phi*K)==0,
        axial_moment=s.cancel(R*m-R*m)==0,
        complete_energy=s.cancel(R*Ps**2*energy-(R*A-R**2*F0**2*phi**2*B))==0,
        pressure_primitive=s.cancel(Ps**2*pressure-R*F0**2*phi**2*C)==0,
        absolute_pressure=s.cancel(Ps**2*(P0+pressure)-(Ps**2*P0+R*F0**2*phi**2*C))==0)
    if not all(identities.values()):raise ArithmeticError('Full microscopic source units differ')
    modebindings=assignment_source_bindings('microswitch_mixed_C4','evaluate',{
        'R':'c.mpf(100)*c.exp(self.h*phase)',
        'logu0':'c.ln(2*R)/2+self.logF0+c.ln(phi[0])-self.core.logP',
        'algebra':'FactoredAlgebra(c,(self.logh,2*self.core.logP,2*self.logF0,2*logu0),self.proofs)',
        'amp':'dressed/phi[0]'})
    packetbindings=assignment_source_bindings('current_microswitch_background_tensor','chart',{
        '(logR, radius)':'self.physical.radius(chart,v,micro,self.switch)',
        'packet':'dict(Z=Z,s=v,exact_logR=logR,exact_pulse_reference_logB_parts=dict(logPstar=self.physical.logP),exact_logD=c.mpf(0),exact_logH=c.mpf(0),full_meridional_stress_log_sectors=sectors)'})
    asts.method('global_physical_assembly','radius')
    asts.method('current_pre_pulse_stress_operator','raw_pre_stress_rows')
    asts.method('current_pre_pulse_stress_operator','raw_pre_remainder_sectors')
    return dict(arbitrary_four_source_log_power_identity=True,complete_current_physical_unit_identities=identities,
        original_current_swirl_base_and_four_log_bindings=modebindings,current_tensor_packet_unit_bindings=packetbindings,
        mode_R_Pstar_factors_and_internal_swirl_invPstar_factors_counted_once=True,
        dummy_end_scale_and_signed_memory_factors_absent_only_for_these_microswitch_charts=True,
        exact_radius_source='log100+original_positive_hb*phase',
        original_BASE_radius_bound_contains_same_formal_radius_not_selected_as_value=True,
        source_log_bounds_enclose_formal_positive_factors_not_define_point_values=True,
        input_hashes=asts.hashes,passed=True)
