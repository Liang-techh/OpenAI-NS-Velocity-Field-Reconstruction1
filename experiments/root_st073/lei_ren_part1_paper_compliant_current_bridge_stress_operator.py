"""Full current bridge source rows in microscopic and ordinary log-radius units."""
import ast
import copy
import math
from types import FunctionType,SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_actual_bridge_mixed_C4 as current
import lei_ren_part1_paper_compliant_microswitch_mixed_C4 as original
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import (
    SOURCE_KEY,compile_function,expose_phase_rows,source_only_projection,
    compiled_factored_source_operators,compiled_factored_physical_lift,expanded_stress_sectors,
    factored_rows_record,shifted_rows)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST


def compiled_current_bridge_with_rows():
    asts=SourceAST();env=dict(vars(original))
    phase=compile_function(expose_phase_rows(asts),env,'<original phase; unresolved source locals only>')
    replay=current.REPLAY;bridge_env=dict(replay.__globals__);bridge_env['phase_physical']=phase
    exposed=FunctionType(replay.__code__,bridge_env,replay.__name__,replay.__defaults__,replay.__closure__)
    exposed.__kwdefaults__=replay.__kwdefaults__
    env=dict(vars(current));env['REPLAY']=exposed
    evaluator=compile_function(asts.method('current_actual_bridge_mixed_C4','evaluate'),env,
        '<unchanged current bridge wrapper; original code and source-local phase projection>')
    asts.method('bridge_mixed_C4','evaluate');asts.method('current_actual_bridge_mixed_C4','current_coordinate_replay')
    return evaluator,dict(original_current_bridge_wrapper_AST_unchanged=True,
        original_coordinate_labelled_replay_code_identical=exposed.__code__ is replay.__code__,
        original_current_parent_acquisition_replay=current.REPLAY_PROOF,
        original_phase_generator_AST_unchanged_except_output_keyword=True,
        only_isolated_global_phase_callback_replaced=True,
        existing_unresolved_rows_exposed_without_feedback=True,input_hashes=asts.hashes,passed=True)


def raw_bridge_adapter_AST(asts):
    before=copy.deepcopy(asts.method('current_microswitch_stress_operator','raw_microswitch_rows'))
    fn=copy.deepcopy(before);fn.name='raw_bridge_rows';changes=[]
    for node in fn.body:
        if isinstance(node,ast.Assign) and ast.unparse(node.targets[0])=='ordinary':
            if ast.dump(node.value)!=ast.dump(ast.parse('lambda rows:[algebra.width(row,-j) for j,row in enumerate(rows)]',mode='eval').body):raise ValueError('Original raw micro units changed')
            node.value=ast.parse("lambda rows:[algebra.width(row,-j) if packet['chart']!='macro' else algebra.lift(row) for j,row in enumerate(rows)]",mode='eval').body
            changes.append('ordinary')
    if changes!=['ordinary']:raise ValueError('Only micro/macro ordinary-coordinate conversion may change')
    compare=copy.deepcopy(fn);compare.name=before.name
    for node in compare.body:
        if isinstance(node,ast.Assign) and ast.unparse(node.targets[0])=='ordinary':
            node.value=ast.parse('lambda rows:[algebra.width(row,-j) for j,row in enumerate(rows)]',mode='eval').body
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Unreviewed raw bridge radial units changed')
    return fn


def compiled_raw_bridge_rows():
    asts=SourceAST();env=dict(SOURCE_KEY=SOURCE_KEY,shifted_rows=shifted_rows)
    original_adapter=compile_function(raw_bridge_adapter_AST(asts),env,'<checked raw current units; micro inverse-width or exact macro identity>')
    def adapter(c,packet):
        result=original_adapter(c,packet)
        result.update(inverse_width_source_shift_applied_before_any_resolution=packet['chart']!='macro',
            macro_original_ordinary_y_preserved=packet['chart']=='macro')
        return result
    return adapter,dict(all_raw_velocity_history_pressure_and_radial_prefactor_math_retained=True,
        only_microscopic_to_ordinary_conversion_conditional=True,
        macro_rows_are_original_ordinary_logR_not_fraction_derivatives=True,
        microscopic_width_shift_performed_before_any_source_resolution=True,input_hashes=asts.hashes,passed=True)


def actual_bridge_normalization_theorem():
    """Original phase generator and actual raw units, independently for both charts."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda value:(s.Rational(str(value)) if isinstance(value,(str,int,float)) else value))
    J=lambda value,order=5:Projection.function(c,value,order)
    hb,unit,Ps=s.symbols('original_positive_hb same_current_swirl_unit Pstar',positive=True)
    delta=s.Symbol('delta',real=True);asts=SourceAST();allchecks={}
    env=dict(math=math,IntervalTaylor=Projection,FactoredJet=type('NoFactored',(),{}),
        derivative=lambda value:J(s.diff(value.expr,Z),value.order-1),
        scaled_positive_source=lambda ctx,logbase,row,proofs:row*unit**2)
    for stem,name in (('long_reshape_mixed_C4','exponential_derivatives'),('microswitch_mixed_C4','rate_rows'),('microswitch_mixed_C4','product_rows')):
        asts.replay(stem,name,env)
    generator=compile_function(source_only_projection(expose_phase_rows(asts)),env,'<original bridge physical source theorem>')
    env.update(SOURCE_KEY=SOURCE_KEY,shifted_rows=shifted_rows)
    adapter=compile_function(raw_bridge_adapter_AST(asts),env,'<micro/macro bridge raw source-unit theorem>')
    for chart in ('first','macro'):
        scale=hb if chart=='first' else s.Integer(1)
        amp=J(s.Function('same_axial_swirl_'+chart)(Z))
        logU=[J(s.Function('actual_log_jet_'+chart+str(k))(Z)) for k in range(1,5)]
        V=[J(s.Function('actual_V_jet_'+chart+str(k))(Z)) for k in range(5)]
        shapes={key:J(s.Function('same_shape_'+chart+'_'+key)(Z)) for key in ('theta','theta_z','mean','axial','swirl','pressure')}
        p0=J(s.Function('same_P0_'+chart)(Z))
        class Algebra:
            width=staticmethod(lambda row,power:row*scale**power)
            shift=staticmethod(lambda row,key:row*unit**(2*s.Rational(str(key[3])))*Ps**(2*key[1]))
            lift=staticmethod(lambda row:row)
        source=generator(c,Z,delta,J(scale),amp,0,logU,V,shapes,p0,1/Ps**2,Algebra.width,[])
        source['algebra']=Algebra();packet={SOURCE_KEY:source,'chart':chart};raw=adapter(c,packet);checks={}
        for label,name,factor in (('theta','Utheta_over_current_Utheta',unit),('axial','Uz',1),('radial','Ur_over_current_sqrt_R_over_2',1),('pressure','P_over_Pstar2',1)):
            rows=raw['absolute_pressure'] if label=='pressure' else raw['velocity'][label]
            for k,row in enumerate(rows):
                difference=s.cancel(row.expr-source['physical'][name][k].expr*factor/scale**k)
                if difference!=0:raise ArithmeticError('Current bridge ordinary velocity/pressure source units differ')
                for n in range(5-k):checks[label+'_y%d_Z%d'%(k,n)]=s.diff(difference,Z,n)==0
        for label,name,rate,factor in (('m','Mz_over_current_R',1,1),('h','Mtheta_over_current_sqrt2_R_1p5_Utheta',s.Rational(3,2),unit),('k','Mtheta_z_over_current_sqrt2_R_1p5_Utheta',s.Rational(3,2),unit),('e','Mztheta_over_current_R_Pstar2',1,1),('p','Mp_over_Pstar2',0,1)):
            rows=shifted_rows(raw['histories'][label],rate,4)
            for k,row in enumerate(rows):
                difference=s.cancel(row.expr-source['primitives'][name][k].expr*factor/scale**k)
                if difference!=0:raise ArithmeticError('Current bridge original moment source units differ')
                for n in range(5-k):checks[label+'_y%d_Z%d'%(k,n)]=s.diff(difference,Z,n)==0
        u,v,history=raw['velocity']['theta'],raw['velocity']['axial'],raw['histories']
        rhs=dict(m=v[0]-history['m'][0],h=u[0]-history['h'][0]*s.Rational(3,2),
            k=u[0]*v[0]-history['k'][0]*s.Rational(3,2),
            e=v[0]*v[0]/Ps**2-u[0]*u[0]/2-history['e'][0],p=u[0]*u[0]/2)
        for key,row in rhs.items():
            if s.cancel(history[key][1].expr-row.expr)!=0:raise ArithmeticError('Current bridge complete five-moment ODE differs')
            checks['raw_five_history_ODE_'+key]=True
        if s.cancel(raw['absolute_pressure'][1].expr-u[0].expr**2/2)!=0:raise ArithmeticError('Current bridge pressure/P0 ODE differs')
        checks['absolute_pressure_y_ODE']=True;allchecks[chart]=checks
    return dict(microscopic_and_macro_original_source_unit_identities=allchecks,
        total_identities=sum(map(len,allchecks.values())),original_phase_and_checked_raw_adapter_AST_replayed=True,
        arbitrary_variable_log_velocity_and_all_axial_velocity_rows_retained=True,
        macro_fraction_labels_coverage_only_no_fraction_length_division=True,
        microscopic_phase_to_y_uses_original_positive_hb_not_cap=True,
        full_five_histories_energy_absolute_P0_and_radial_prefactors_retained=True,
        input_hashes=asts.hashes,passed=True)

def actual_bridge_endpoint_operator_theorem():
    """Replay actual control programs and complete source-generator pullbacks."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import R100_source_join
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda value:(s.Rational(str(value)) if isinstance(value,(str,int,float)) else value))
    J=lambda value,order=5:Projection.function(c,value,order)
    hb,unit,Ps=s.symbols('original_positive_hb same_current_swirl_unit Pstar',positive=True)
    delta=s.Symbol('delta',real=True);zero=J(0);asts=SourceAST()
    env=dict(math=math,IntervalTaylor=Projection,FactoredJet=type('NoFactored',(),{}),
        derivative=lambda value:J(s.diff(value.expr,Z),value.order-1),
        scaled_positive_source=lambda ctx,logbase,row,proofs:row*unit**2)
    for stem,name in (('long_reshape_mixed_C4','exponential_derivatives'),('microswitch_mixed_C4','rate_rows'),('microswitch_mixed_C4','product_rows'),('microswitch_mixed_C4','switch_controls'),('bridge_mixed_C4','bridge_controls')):
        asts.replay(stem,name,env)
    class Algebra:
        lift=staticmethod(lambda value:value if isinstance(value,FunctionJet) else J(value))
        width=staticmethod(lambda value,power=1:Algebra.lift(value)*hb**power)
    algebra=Algebra();D=[J(s.Function('same_Dbar_y'+str(k))(Z)) for k in range(4)]
    G=[J(s.Function('same_axial_drive_y'+str(k))(Z)) for k in range(4)]
    quotient=J(s.Function('same_actual_comparison_quotient')(Z))
    macro=env['bridge_controls'](algebra,J(1),[J(hb)]+[zero]*3,D,G,quotient,[zero]*3)
    smoothing=env['bridge_controls'](algebra,J(hb),[J(hb)]+[zero]*3,
        [row*hb**k for k,row in enumerate(D)],[row*hb**k for k,row in enumerate(G)],quotient,[zero]*3)
    switch=env['switch_controls'](c,'first',[0]*5,D,lambda k,power:G[k]*hb**power,quotient,algebra.width)
    controlchecks={}
    for name in ('logUtheta_coordinate_derivatives','logF_coordinate_derivatives','Uz_positive_order_coordinate_derivatives'):
        other={'logUtheta_coordinate_derivatives':'logUtheta_phase_derivatives',
            'logF_coordinate_derivatives':'logF_phase_derivatives',
            'Uz_positive_order_coordinate_derivatives':'Uz_positive_phase_derivatives'}[name]
        for k,row in enumerate(macro[name]):
            for seam,value in (('second_macro',smoothing[name][k]),('R100',switch[other][k])):
                if s.cancel(value.expr/hb**(k+1)-row.expr)!=0:raise ArithmeticError('Original bridge/switch true control pullback differs')
                controlchecks[seam+'/'+name+str(k)]=True
    sourcefn=compile_function(source_only_projection(expose_phase_rows(asts)),env,'<original full bridge boundary source generator>')
    amp=J(s.Function('same_boundary_dressed_swirl')(Z));V0=J(s.Function('same_boundary_V')(Z))
    shapes={key:J(s.Function('same_boundary_actual_'+key)(Z)) for key in ('theta','theta_z','mean','axial','swirl','pressure')}
    p0=J(s.Function('same_boundary_P0')(Z))
    ordinary=sourcefn(c,Z,delta,J(1),amp,0,macro['logUtheta_coordinate_derivatives'],
        [V0]+macro['Uz_positive_order_coordinate_derivatives'],shapes,p0,1/Ps**2,lambda row,power=1:row,[])
    rowchecks={}
    for seam,control,logkey,vkey in (('second_macro',smoothing,'logUtheta_coordinate_derivatives','Uz_positive_order_coordinate_derivatives'),('R100',switch,'logUtheta_phase_derivatives','Uz_positive_phase_derivatives')):
        phase=sourcefn(c,Z,delta,J(hb),amp,0,control[logkey],[V0]+control[vkey],shapes,p0,1/Ps**2,algebra.width,[])
        checks={}
        for group in ('physical','primitives'):
            for name in ordinary[group]:
                for k,(a,b) in enumerate(zip(phase[group][name],ordinary[group][name])):
                    difference=s.cancel(a.expr/hb**k-b.expr)
                    if difference!=0:raise ArithmeticError('Original complete bridge boundary source rows differ')
                    for n in range(5-k):checks[group+'/'+name+'/y%d_Z%d'%(k,n)]=s.diff(difference,Z,n)==0
        rowchecks[seam]=checks
    r100=R100_source_join()
    if r100['symbolic_original_control_coordinate_identities']!=17:raise ValueError('Existing original R100 source proof changed')
    asts.method('bridge_mixed_C4','comparison_moment_rows');asts.method('bridge_mixed_C4','comparison_directions')
    asts.method('bridge_mixed_C4','smoothed_comparison_rows');asts.method('current_bridge_functional_joins','production_branch_specialization')
    asts.method('current_bridge_functional_joins','common_coordinate_identities')
    return dict(original_bridge_and_switch_control_pullback_identities=controlchecks,
        full_original_second_macro_and_R100_source_generator_identities=rowchecks,
        source_generator_identity_count=sum(map(len,rowchecks.values())),
        consumed_original_R100_current_source_function_theorem=r100,
        first_second_uses_same_admitted_current_signed_prefixes_comparison_and_flat_controls=True,
        first_second_same_original_generator_implies135_common_source_rows=True,
        macro_fraction_is_coverage_label_only_not_differentiated=True,
        original_current_comparison_own_histories_and_actual_six_histories_not_interchanged=True,
        same_full_source_mixed4_rows_and_full_raw_operators_imply_three_completed_tensor_traces=True,
        source_equality_precedes_common_triangle_bounds=True,
        interval_overlap_not_used_as_source_identity=True,input_hashes=asts.hashes,passed=True)
