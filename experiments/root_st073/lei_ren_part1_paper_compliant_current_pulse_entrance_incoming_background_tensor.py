"""Current entrance and whole O3 incoming power tensors on one checked graph."""
import ast
import copy
import gzip
import hashlib
import importlib
import json
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_pulse_main_exit_background_tensor as main
from lei_ren_part1_paper_compliant_current_pulse_main_exit_background_tensor import (
    CurrentPulseMainExitBackgroundTensor,HERE,PREFIX,OPEN,SourceAST,sha,pack,encode,
    endpoints,source_precision,accepted,_verify_hashes,canonical_tensor_groups)
from lei_ren_part1_paper_compliant_actual_Rp_source_join import (
    source_recipe_bindings,native_source_graph,native_inlet_projection_proof,
    mixed_source_join_proof,actual_parameter_formula_bindings)
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import source_histories
from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import original_history_source_bindings
from lei_ren_part1_paper_compliant_actual_Rp_source_join import terminal_shape_proof
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,selected_node,keywords_binding
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_velocity_rows,lift_physical_packet
from lei_ren_part1_paper_compliant_pulse_main_exit_physical_C2 import main_exit_to_physical_packet
from lei_ren_part1_paper_compliant_current_pulse_gap_background_tensor import copy_jet,IntervalTaylor
import math
from lei_ren_part1_paper_compliant_global_physical_assembly import PULSE,CompliantGlobalPhysicalAssembly as BASE
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import (
    forward_entrance_energy_rows,gp_energy)

NAME=PREFIX+'current_pulse_entrance_incoming_background_tensor.json.gz'
RECEIPT=PREFIX+'current_pulse_entrance_incoming_background_tensor_check.json'
GATES=('current_actual_entrance_incoming_full_tensors_available',
    'current_actual_entrance_incoming_full_meridional_decomposition_available',
    'current_actual_incoming_entrance_completed_tensor_join_certified',
    'current_actual_entrance_main_completed_tensor_join_certified')
SEAMS=('incoming_entrance','entrance_main')
VIEWS={'whole_entrance':('pulse_entrance',(-1,1),(0,'.02'),('-3','-1'),None,'1'),
    'whole_incoming':('O3_power',(-1,1),(0,1),('-3','-1'),None,'1'),
    'inlet':('pulse_entrance',(-1,1),0,'-1',None,'1'),
    'entrance_main':('pulse_entrance',(-1,1),'.02','-1',None,'1'),
    'incoming_left':('O3_power',(-1,1),0,'-1',None,'1'),
    'fresh_entrance':('pulse_entrance','.537','.01337','-2.6','.41','.8'),
    'fresh_incoming':('O3_power','.731','.537','-2.6','.41','.8'),
    'early_y':('early_y','.537',(0,1),'-2.6','.41','.8')}
LEGACY_ALLOWLIST = {('pulse_entrance_similarity_C4', 'entrance_source_proof'): ('8afd49826cc7835fa0fca307f071033f715ffa347064df37ea2918e9587cee23', '9fc990c72a76e5d45c10d164525b98e8be4c1050ef9ea773f402e7c9311f61e7', '9347a3865d6b75a8b741776bb5e285fb750161389952c81febc5ae975aa9e4c3', '1886b6b0848c2b860eb04743d12102e23a4f80a08607f4d5bca89c1fd1513ba5', 'a6fe9614e038624cb6e0df88a32cf243fe40a8ae519cd2b8283a4216051a802b'), ('pulse_entrance_physical_C2', 'entrance_physical_binding'): ('76081c01f31c56eb6739f3a508fec358f7522e88bdf3c265b074297894f6594d', 'c019e51021b3d0368e1fde4066a2d3ac5ddc847f42620cee04cff4b71dc72444', '13fe00e11825d01d96a193b7f01f6a2698919d92f555707386b9469ef25a3694', 'f470beddaba8d422410d2746864b356cccf661a03fbf9fc81c5f335af685e041', 'd684b801bd8754c02585c65ace6fcaacc6215c0050a8f6a7d7a790c09cf3faff', '9f433d148b2d0a51980e365cf4f809c281fd419ced4355f039340c279bd358a7')}


def generic_entrance_proof(stem,name,original):
    module=importlib.import_module(PREFIX+stem);asts=SourceAST()
    fn=copy.deepcopy(asts.method(stem,name));fn.decorator_list=[]
    legacy=({'records','gates','pressure_gates'} if name=='entrance_source_proof'
        else {'records','gate','original','inlet'})
    kept=[];removed=[];digests=[]
    for node in fn.body:
        names={v.id for v in ast.walk(node) if isinstance(v,ast.Name)}
        if not isinstance(node,ast.Return) and names.intersection(legacy):
            removed.append(ast.unparse(node).splitlines()[0])
            digests.append(hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest())
        else:kept.append(node)
    if tuple(digests)!=LEGACY_ALLOWLIST.get((stem,name)):raise ValueError('Unreviewed entrance legacy receipt removal')
    fn.body=kept;env=dict(vars(module));env['original']=original
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original generic entrance theorem,current source admission separate>','exec'),env)
    result=env[name]({})
    if not all(result['identities'].values()):raise ValueError('Original generic entrance identity failed')
    result.update(omitted_legacy_receipt_checks=removed,omitted_legacy_statement_AST_sha256=digests,
        exact_reviewed_legacy_statement_allowlist_enforced=True,input_hashes={**result['input_hashes'],**asts.hashes})
    return result


def production_power_function_theorem(field):
    """Replay both production recurrences on shared exact kernel functions.

    Equality is inductive over the actual O2/O3 algorithms, not from a hash
    or overlap. Quadrature partitions enclose these defining kernels.
    """
    asts=SourceAST();identities={};recipes=source_recipe_bindings()
    graph=native_source_graph(field.physical.pre,field.pulse)
    parameter=actual_parameter_formula_bindings()
    live_program=actual_live_program_bindings(field)
    def zero(name,a,b):
        if a!=b and s.cancel(s.expand_power_exp(a-b))!=0:raise ArithmeticError('Current production source differs: '+name)
        identities[name]=True
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)) if isinstance(v,(str,int,float)) else v,exp=s.exp)
    z,mu,t,invP=s.symbols('Z mu t invPstar_squared',real=True)
    qi=1/(1+z*z);u1=s.Function('shared_inlet_u')(z)
    names=('m','h','k','e','p');old={k:s.Function('shared_inlet_'+k)(z) for k in names}
    keys=dict(m='Mz_over_R',h='Mtheta_over_sqrt2_R_3half_Pstar',k='Mtheta_z_over_sqrt2_R_3half_Pstar',
        e='Mztheta_over_R_Pstar_squared',p='Mp_over_Pstar_squared')
    get=lambda key:old[next(k for k,v in keys.items() if v==key)]
    J,I0,I1,I2,KB,KB2,KT,KE,KP=s.symbols('J I0 I1 I2 KB KB2 KT KE KP',real=True)
    stages=(('slope','outer_initial'),('axial','outer_initial'),('slope_mu','outer_buffer'),('power','outer_buffer'))
    schemas={}
    for method,native_stem in stages:
        d=s.exp(-t);root=s.exp(-t/2);d3=s.exp(-3*t/2)
        K=dict(B_mass=KB,B_squared_mass=KB2,J=J,theta=KT,energy=KE,pressure=KP)
        env=dict(c=c,z=z,qi=qi,t=t,y=t,yy=t,mu=mu,J=J,integrals=[I0,I1,I2],u1=u1,
            self=SimpleNamespace(invP2=invP),get=get,K=K,kernels=K,mass=[I0,I1,I2],decay=d,root_decay=root,decay3=d3,d3=d3)
        current=dict(z=z,qi=qi,t=t,y=t,J=J,V=4*z,u1=u1,d=d,root=root,d3=d3,mu=mu,
            **{'self.invP2':invP,'mass[0]':I0,'mass[1]':I1,'mass[2]':I2,
                "K['B_mass']":KB,"K['B_squared_mass']":KB2,"K['theta']":KT,"K['energy']":KE,"K['pressure']":KP,
                **{"old['"+k+"']":v for k,v in old.items()}})
        if method=='slope':
            env['factor']=asts.evaluate(asts.expression(native_stem,method,'factor'),env)
            current['factor']=env['factor'];current['h']=asts.evaluate(asts.expression('pre_pulse_mixed_C4',method,'h'),env|{'y':t})
        elif method=='slope_mu':
            env['f']=asts.evaluate(asts.expression(native_stem,method,'f'),env)
            current['factor']=env['f'];env['decay']=d
        elif method=='power':
            env['f']=asts.evaluate(asts.expression(native_stem,method,'f'),env|{'slope':-s.Rational(1,2)-mu})
            env['theta_kernel']=(env['f']-d3)/(1-mu);current['theta']=env['theta_kernel']
            current.update({'decay_integral(c, 2 * mu, t)':(1-s.exp(-2*mu*t))/(2*mu),
                'decay_integral(c, 1 + 2 * mu, t)':(1-s.exp(-(1+2*mu)*t))/(1+2*mu)})
            env['decay_integral']=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate
        native_u=asts.evaluate(asts.expression(native_stem,method,'u'),env)
        native_values={k:asts.evaluate(asts.expression(native_stem,method,k),env|{'h':asts.evaluate(asts.expression(native_stem,method,'h'),env),'V':4*z}) for k in names}
        current_u=(qi*current['factor'] if method=='slope' else u1*root if method=='axial'
            else u1*current['factor'] if method=='slope_mu' else u1*env['f'])
        values=source_histories(method,current)
        schemas[method]=(dict(u=current_u,**values),dict(u=native_u,**native_values))
        zero(method+'_same_u',current_u,native_u)
        for key in names:zero(method+'_same_'+key,values[key],native_values[key])
    # Compose the actual parent chain independently on both sides. The
    # first stage has the defining reference seed, not arbitrary old data.
    Md,Tw,phase=s.symbols('same_Md same_Tw actual_phase',positive=True)
    composed=[{},{}];chain={};ratios={}
    slope_J,transition_J=s.symbols('actual_slope_J_at1 actual_transition_J_at1',real=True)
    for method,_ in stages:
        time=dict(slope=s.Integer(1),axial=s.exp(Md)+10,slope_mu=s.Integer(1),power=Tw*phase)[method]
        for side in range(2):
            previous=composed[side]
            substitutions={t:time,J:slope_J if method=='slope' else transition_J}
            if previous:substitutions.update({u1:previous['u'],**{old[key]:previous[key] for key in names}})
            composed[side]={key:s.factor(expr.subs(substitutions,simultaneous=True)) for key,expr in schemas[method][side].items()}
        for key in ('u',)+names:zero('composed_'+method+'_'+key,composed[0][key],composed[1][key])
        ratio=s.cancel(composed[0]['h']/composed[0]['u'])
        zero('composed_'+method+'_X_is_exact_Z0_ratio',ratio,ratio.subs(z,0))
        ratios[method]=str(ratio)
        chain[method]={key:str(value) for key,value in composed[0].items()}
    terminal={key:s.factor(value.subs(phase,1)) for key,value in composed[0].items()}
    coefficients=dict(U=s.cancel(terminal['u']/qi),M=s.cancel(terminal['m']/z),
        H=s.cancel(terminal['h']/qi),K=s.cancel(terminal['k']/(z*qi)),Pin=s.cancel(terminal['p']/qi**2))
    epoly=s.Poly(s.cancel(terminal['e']/qi**2),z)
    coefficients.update(EZ=epoly.coeff_monomial(z**6),EQ=epoly.coeff_monomial(1))
    if any(z in value.free_symbols for value in coefficients.values()):raise ArithmeticError('Composed terminal coefficient depends on Z')
    canonical=dict(u=coefficients['U']*qi,m=coefficients['M']*z,h=coefficients['H']*qi,
        k=coefficients['K']*z*qi,e=coefficients['EZ']*z*z+coefficients['EQ']*qi**2,p=coefficients['Pin']*qi**2)
    for key in canonical:zero('actual_phase1_live_native_inlet_'+key,terminal[key],canonical[key])
    zero('actual_phase1_Xp',terminal['h']/terminal['u'],coefficients['H']/coefficients['U'])
    # Actual live constants are selected by these original source calls.
    # The composed native packet supplies their defining coefficient
    # functions; quadrature cells only enclose them.
    native_terminal={key:value.subs(phase,1) for key,value in composed[1].items()}
    p0packet={label:[native_terminal[key].subs(z,0),s.diff(native_terminal[key],z).subs(z,0)] for key,label in dict(u='Utheta_over_Pstar',**keys).items()}
    constenv=dict(c=c,box=lambda value:value,p0=p0packet,td=s.exp(Md)+10,
        initial=SimpleNamespace(invP2=invP,params=SimpleNamespace(Tw=Tw,logPstar=s.Symbol('same_logPstar'))),
        kernels=dict(B_squared_mass=KB2))
    source_constants={}
    for key,target in dict(U='U',M='M',K='K',EQ='EQ').items():
        source_constants[key]=asts.evaluate(asts.expression('axial_high_jets','_incoming_constants',target),constenv)
        zero('live_high_constant_'+key,source_constants[key],coefficients[key])
    ez=asts.evaluate(asts.expression('axial_high_jets','_incoming_constants','EZ'),constenv)
    aug=[node for node in ast.walk(asts.method('axial_high_jets','_incoming_constants')) if isinstance(node,ast.AugAssign) and ast.unparse(node.target)=='EZ' and isinstance(node.op,ast.Mult)]
    if len(aug)!=1:raise ValueError('Original direct positive EZ multiplier changed')
    ez*=asts.evaluate(aug[0].value,constenv);zero('live_high_constant_EZ',ez,coefficients['EZ'])
    for key,label in dict(H='Mtheta_over_sqrt2_R_3half_Pstar',Pin='Mp_over_Pstar_squared').items():
        env0=dict(p0=p0packet,box=lambda value:value)
        value=asts.evaluate(asts.expression('pulse_radial_C4','__init__','self.inlet_H' if key=='H' else 'self.inlet_P'),env0)
        zero('live_pulse_ctor_'+key,value,coefficients[key])
    # Replay native data with independent coefficient symbols; substitute
    # the composed coefficient FUNCTIONS into its exact C5 projection.
    U0,M0,H0,K0,EZ0,EQ0,Pin0=s.symbols('U M H K EZ EQ Pin',real=True)
    Ps=s.Symbol('Pstar',positive=True)
    constants=dict(U=U0,M=M0,K=K0,E_Z=EZ0,E_Q=EQ0,C1=M0/(Ps*U0),C2=K0/(Ps*U0**2),C0=EQ0/U0**2,C_E=EZ0/U0**2)
    jet=lambda expr:FunctionJet.function(c,expr,5)
    from lei_ren_part1_paper_candidate_pressure_function import q_jets
    env=dict(IntervalTaylor=FunctionJet,endpoints=lambda value:(0,0),mp=SimpleNamespace(workdps=lambda precision:__import__('contextlib').nullcontext()),q_jets=q_jets,copy_jet=lambda ctx,value:value)
    oldincoming=asts.replay('axial_high_jets','incoming',env)(SimpleNamespace(ctx=c,constants=constants),z)
    append=asts.replay('fifth_axial_jets','append_fifth',env)
    env.update(c=c,Z=z,incoming=oldincoming,k=constants,append_fifth=append)
    moments=asts.evaluate(selected_node('fifth_axial_jets','select','moments'),env)
    energy=asts.evaluate(selected_node('fifth_axial_jets','select','energy'),env)
    env.update(moments=moments,energy=energy)
    publication=[kw.value for node in ast.walk(asts.method('fifth_axial_jets','select')) if isinstance(node,ast.Call) and ast.unparse(node.func)=='result.update' for kw in node.keywords if kw.arg=='incoming']
    if len(publication)!=1:raise ValueError('Original fifth incoming publication changed')
    selected=dict(incoming=asts.evaluate(publication[0],env))
    owner=SimpleNamespace(ctx=c,data_cache={},fifth=SimpleNamespace(select=lambda Z:selected),high=SimpleNamespace(constants=constants),inlet_H=H0,inlet_P=Pin0)
    _,inlet,uu,ee,mm=asts.replay('pulse_radial_C4','data',env)(owner,z)
    substitute={U0:coefficients['U'],M0:coefficients['M'],H0:coefficients['H'],K0:coefficients['K'],EZ0:coefficients['EZ'],EQ0:coefficients['EQ'],Pin0:coefficients['Pin']}
    abstract=dict(u=U0*qi,m=M0*z,h=H0*qi,k=K0*z*qi,e=EZ0*z*z+EQ0*qi**2,p=Pin0*qi**2)
    labels=dict(u='Utheta_over_Pstar',**keys)
    for key,label in labels.items():
        for n,value in enumerate(inlet[label]):
            zero('composed_native_data_'+key+'_Z'+str(n),value,s.diff(abstract[key],z,n)/math.factorial(n))
    # Since the proved substitutions are Z-independent, these derivative
    # equalities remain exact after inserting the actual source functions.
    for key,value,wanted in (('M',mm[0].expr,abstract['m']/(Ps*abstract['u'])),
            ('N',mm[1].expr,abstract['k']/(Ps*abstract['u']**2)),('E',ee.expr,abstract['e']/abstract['u']**2)):
        zero('composed_actual_main_data_'+key,value,wanted)
    # This symbol denotes the live fourteen-atom analytic source already
    # identified by the original callable, definition and complete read
    # state in live_program. It is not a newly chosen pressure function.
    P0=s.Function('current_full_fourteen_atom_preheat_pressure')(z)
    datum=SimpleNamespace(normalized_jets=lambda domain,order:dict(normalized_pressure_coefficients=[s.diff(P0,z,n)/math.factorial(n) for n in range(order+1)]))
    datum_env=dict(c=c,Z=z,endpoints=lambda value:(value,value),self=SimpleNamespace(datum=datum))
    pre_datum=asts.evaluate(asts.expression('pre_pulse_mixed_C4','packet','raw'),datum_env)
    native_env=datum_env|dict(self=SimpleNamespace(selection=SimpleNamespace(future=SimpleNamespace(angular=SimpleNamespace(initial=SimpleNamespace(datum=datum))))))
    native_datum=asts.evaluate(asts.expression('pulse_radial_C4','pressure_moment','rows'),native_env)
    for n,(a,b) in enumerate(zip(pre_datum,native_datum)):zero('actual_live_pre_native_pressure_datum_Z'+str(n),a,b)
    pressure_context=SimpleNamespace(mpf=lambda value:value[0] if isinstance(value,(tuple,list)) else c.mpf(value))
    pressure_env=dict(IntervalTaylor=FunctionJet,c=pressure_context,inlet=inlet,u=uu,kernel=s.Integer(0),rows=native_datum,endpoints=lambda value:(value,value))
    native_p=asts.evaluate(asts.expression('pulse_radial_C4','pressure_moment','p'),pressure_env)
    native_p0=asts.evaluate(asts.expression('pulse_radial_C4','pressure_moment','p0'),pressure_env)
    zero('actual_live_phase1_native_absolute_pressure',native_p.expr+native_p0.expr,abstract['p']+P0)
    asts.expression('current_pulse_main_exit_background_tensor','_data','energy',wanted="copy_jet(c,selected['incoming']['energy_Taylor'])")
    asts.expression('current_pulse_main_exit_background_tensor','_data','incoming',wanted="[copy_jet(c,value) for value in selected['incoming']['moment_Taylor']]")
    asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power','X',wanted="c.mpf(endpoints(scalar['actual_normalized_primitive_y_derivative_axial5']['h'][0][0]))/c.mpf(endpoints(scalar['Utheta_over_Pstar_axial5_coefficients'][0]))")
    asts.method('pre_pulse_mixed_C4','slope_masses')
    projection=native_inlet_projection_proof();mixed=mixed_source_join_proof()
    return dict(identities=identities,source_recipe_bindings=recipes,current_source_graph=graph,parameter_source_bindings=parameter,
        actual_live_program_and_pressure_read_state=live_program,
        actual_parent_call_and_history_chain=original_history_source_bindings(),
        actual_full_canonical_Z_shapes=terminal_shape_proof(),
        independently_composed_actual_parent_chain_source_functions=chain,
        actual_phase1_coefficient_functions={key:str(value) for key,value in coefficients.items()},
        actual_live_high_and_pulse_ctor_source_constant_expressions_replayed=True,
        actual_live_incoming_C4_C5_publication_and_native_data_AST_replayed=True,
        actual_pre_and_native_pressure_datum_calls_replayed_on_same_defining_callable=True,
        checked_main_closed_Rv_pressure_bridge_consumed=field.main_tensor.proof['current_reduced_absolute_Rv_pressure_and_native_Pin_retained'],
        actual_X_source_ratio_by_stage=ratios,
        distinct_slope_and_mu_transition_integral_sources_preserved=True,
        actual_X_Z0_assignment_bound_to_exact_whole_Z_q_cancellation=True,
        actual_phase1_functions_substituted_into_live_native_C5_data=True,
        source_kernel_partitions_do_not_choose_point_values=True,
        same_reference_seed_and_same_exact_scalar_kernels_inductively_identify_both_production_endpoints=True,
        actual_current_O3_power_phase1_equals_live_native_buffer_power_for_whole_Z=True,
        native_C5_incoming_projection=projection,native_physical_mixed_join=mixed,
        input_hashes={**asts.hashes,**projection['input_hashes'],**mixed['input_hashes'],**parameter.get('input_hashes',{}),**live_program['input_hashes']},passed=True)


def actual_live_program_bindings(field):
    """Tie the formal replay to the actual checked production callables.

    For the pressure enclosure algorithm, every self field read is
    compared, not merely its hash. Identical original program, arguments
    and defining input state establish function congruence.
    """
    pre=field.physical.pre;native=field.pulse;asts=SourceAST();checks={}
    owners=(('pre_pulse_mixed_C4','CompliantPrePulseMixedC4',pre,('slope','inlet','axial','slope_mu','power','packet')),
        ('outer_initial','SharedOuterInitial',native.pulse.initial,('slope','axial')),
        ('outer_buffer','SharedOuterBuffer',native.pulse.buffer,('slope_mu','power')),
        ('axial_high_jets','CompliantAxialHighJets',native.fifth.fourth,('_incoming_constants','incoming')),
        ('fifth_axial_jets','CompliantFifthAxialJets',native.fifth,('select',)),
        ('pulse_radial_C4','CompliantPulseRadialC4',native,('data','pressure_moment')))
    for stem,cls,owner,methods in owners:
        original=getattr(importlib.import_module(PREFIX+stem),cls)
        for method in methods:
            if getattr(owner,method).__func__ is not getattr(original,method):raise ValueError('Actual production callable differs: '+stem+'.'+method)
            asts.method(stem,method);checks[stem+'.'+method]=True
    for key in ('Md','mu','Tw','yd','logPstar'):
        if encode(pack(getattr(pre.params,key)))!=encode(pack(getattr(native.pulse.buffer.params,key))):raise ValueError('Actual production parameter input differs: '+key)
        checks['actual_pre_native_parameter_'+key]=True
    if encode(pack(pre.invP2))!=encode(pack(native.pulse.initial.invP2)):raise ValueError('Actual reference inverse Pstar seed differs')
    checks['actual_pre_native_inverse_Pstar_squared_seed']=True
    ctor=asts.method('axial_high_jets','__init__')
    calls=[node for node in ast.walk(ctor) if isinstance(node,ast.Call) and ast.unparse(node.func)=='self._incoming_constants']
    high_class=getattr(importlib.import_module(PREFIX+'axial_high_jets'),'CompliantAxialHighJets')
    if len(calls)!=1 or calls[0].args or calls[0].keywords or type(native.fifth.fourth) is not high_class or native.high.constants is not native.fifth.fourth.constants:
        raise ValueError('Live constructor-created incoming constants witness differs')
    checks['live_constants_object_is_original_constructor_incoming_constants_output']=True
    checks.update({'native_packet_'+key:value for key,value in keywords_binding('outer_initial','packet','result',dict(
        Utheta_over_Pstar='list(u.coefficients)',Mz_over_R='list(m.coefficients)',
        Mtheta_over_sqrt2_R_3half_Pstar='list(h.coefficients)',Mtheta_z_over_sqrt2_R_3half_Pstar='list(k.coefficients)',
        Mztheta_over_R_Pstar_squared='list(e.coefficients)',Mp_over_Pstar_squared='list(p.coefficients)',
        P0_over_Pstar_squared='list(p0.coefficients)',P_over_Pstar_squared='list(P.coefficients)')).items()})
    checks.update({'actual_pre_packet_'+key:value for key,value in keywords_binding('pre_pulse_mixed_C4','packet','result',dict(
        Utheta_over_Pstar_axial5_coefficients='list(u.coefficients)',original_P0_axial5_coefficients='list(p0.coefficients)')).items()})
    physical=asts.method('pre_pulse_mixed_C4','physical_mixed')
    output=next(node.value for node in ast.walk(physical) if isinstance(node,ast.Return))
    fields={kw.arg:ast.unparse(kw.value) for kw in output.keywords}
    if fields['actual_normalized_primitive_y_derivative_axial5']!='rows':raise ValueError('Actual raw history publication differs')
    asts.expression('pre_pulse_mixed_C4','physical_mixed','rows',wanted='{name:[value] for name,value in history.items()}')
    checks['actual_pre_raw_packet_publishes_defining_five_histories']=True
    datums=(pre.datum,native.pulse.initial.datum,native.selection.future.angular.initial.datum)
    original=getattr(importlib.import_module('lei_ren_part1_paper_logarithmic_pressure_datum'),'LogarithmicPressureDatum')
    actual_datum_class=getattr(importlib.import_module(PREFIX+'pressure_source'),'CompliantPressureDatum')
    fn=ast.parse((HERE/'lei_ren_part1_paper_logarithmic_pressure_datum.py').read_text(encoding='utf8'))
    fn=next(node for node in ast.walk(fn) if isinstance(node,ast.FunctionDef) and node.name=='normalized_jets')
    reads={node.attr for node in ast.walk(fn) if isinstance(node,ast.Attribute) and isinstance(node.value,ast.Name) and node.value.id=='self'}
    wanted={'ctx','stages','flatten_complex_upper','rho','m2','m0','parameters','source_sha','datum_sha','input_hashes'}
    if reads!=wanted:raise ValueError('Unreviewed pressure algorithm read state')
    for datum in datums:
        if type(datum) is not actual_datum_class or not isinstance(datum,original) or datum.normalized_jets.__func__ is not original.normalized_jets:raise ValueError('Original actual compliant analytic pressure callable required')
        if datum.definition!=datums[0].definition or len(datum.stages)!=14:raise ValueError('Actual fourteen-atom pressure definition differs')
        for key in wanted-{'ctx','parameters'}:
            if encode(pack(getattr(datum,key)))!=encode(pack(getattr(datums[0],key))):raise ValueError('Actual pressure program input differs: '+key)
        if type(datum.ctx) is not type(datums[0].ctx) or datum.ctx.dps!=datums[0].ctx.dps:raise ValueError('Actual pressure arithmetic differs')
        if encode(pack(datum.parameters.logPstar))!=encode(pack(datums[0].parameters.logPstar)):raise ValueError('Actual pressure log-scale source differs')
    checks['original_pressure_algorithm_identical_complete_read_state_and_definition']=True
    name='lei_ren_part1_paper_logarithmic_pressure_datum.py';asts.hashes[name]=sha(name)
    asts.method('pressure_source','__init__')
    # Consume the actual output packets as well as their defining method
    # ASTs. Their enclosures need not numerically coincide: the source
    # algorithm congruence above identifies the functions before boxes.
    live_power=native.pulse.buffer.power(0,1,cells=128)
    live_data=native.data(pre.ctx.mpf((-1,1)))
    live_pressure=datums[0].normalized_jets((-1,1),5)
    return_fields=('Utheta_over_Pstar','Mz_over_R','Mtheta_over_sqrt2_R_3half_Pstar',
        'Mtheta_z_over_sqrt2_R_3half_Pstar','Mztheta_over_R_Pstar_squared','Mp_over_Pstar_squared','P0_over_Pstar_squared')
    if any(len(live_power[key])!=2 for key in return_fields):raise ValueError('Live native power source packet lost')
    if any(value.order!=5 for value in live_data[4]+[live_data[3],live_data[2]]):raise ValueError('Live native C5 endpoint output lost')
    checks['actual_live_buffer_power_data_incoming_energy_and_pressure_outputs_consumed']=True
    return dict(identities=checks,pressure_algorithm_self_read_fields=sorted(reads),
        live_incoming_constants_witness='Original high constructor calls _incoming_constants once; live high/fourth share that exact output dictionary; its original assignments are replayed against composed source functions.',
        actual_live_native_power_Z0_fields={key:live_power[key] for key in return_fields},
        actual_live_high_source_constants={key:native.high.constants[key] for key in ('U','M','K','E_Q','E_Z','C1','C2','C0','C_E')},
        actual_live_native_C5_inlet=live_data[1],actual_live_native_C5_incoming_moments=live_data[4],actual_live_native_C5_incoming_energy=live_data[3],
        actual_live_pressure_datum_axial5_output=live_pressure['normalized_pressure_coefficients'],
        same_exact_analytic_fourteen_atom_source_before_different_enclosure_layouts=True,
        deterministic_program_congruence_not_interval_overlap=True,input_hashes=asts.hashes,passed=True)


def actual_power_velocity_source_theorem():
    """Bind the lift's full velocity to the actual pre physical source."""
    z=s.Symbol('Z',real=True);U,M,Ps=s.symbols('actual_scalar_U actual_scalar_M Pstar',positive=True)
    mu,delta=s.symbols('mu delta',real=True);C=1/(1+z*z)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)) if isinstance(v,(str,int,float)) else v)
    jet=lambda expr,order=5:FunctionJet.function(c,expr,order)
    env=dict(IntervalTaylor=FunctionJet,square=lambda v:v*v,derivative=lambda v:jet(s.diff(v.expr,z),v.order-1),
        axial_derivative=lambda v:s.diff(v,z),math=math)
    asts=SourceAST()
    for stem,name in (('pre_pulse_mixed_C4','rate_rows'),('pre_pulse_mixed_C4','product_rows'),
            ('long_reshape_mixed_C4','exponential_derivatives'),('collar_stress_C3','shifted_rows')):asts.replay(stem,name,env)
    prefn=asts.replay('pre_pulse_mixed_C4','physical_mixed',env)
    fn=asts.replay('pulse_end_physical_C2','pulse_velocity_rows',env)
    zero=jet(0);hist=dict(m=jet(M*z),h=jet(s.Symbol('H')*C),k=jet(s.Symbol('K')*z*C),
        e=jet(s.Symbol('EZ')*z*z+s.Symbol('EQ')*C*C),p=jet(s.Symbol('Pin')*C*C))
    raw=prefn(c,z,delta,jet(U*C),[zero-(s.Rational(1,2)+mu)]+[zero]*3,[zero]*5,hist,jet(s.Function('same_P0')(z)),1/Ps**2)
    m1=M*z/(Ps*U*C);mm=[m1*(-(s.Rational(1,2)-mu))**j for j in range(5)]
    velocity=fn(c,delta,mu,z,C,[s.Integer(0)]*5,mm);checks={}
    mapping=dict(theta=('Utheta_over_Pstar',U),axial=('Uz',Ps*U),radial=('Ur_over_current_sqrt_R_over_2',Ps*U))
    for key,(label,factor) in mapping.items():
        for j in range(5):
            for n in range(5-j):
                a=raw['physical_velocity_pressure_y_Z_mixed4'][label]['y'+str(j)+'_Z'+str(n)]
                b=s.diff(velocity[key][j]*factor,z,n)
                if s.cancel(a-b)!=0:raise ArithmeticError('Actual pre/source lift velocity differs: '+key+str((j,n)))
                checks[key+'_y'+str(j)+'_Z'+str(n)]=True
    return dict(identities=checks,actual_pre_physical_mixed_and_lift_velocity_AST_replayed=True,
        same_actual_velocity_function_before_BASE_Cartesian_source_map=True,
        scope_velocity_binding_only_pressure_and_NS_operator_admitted_separately=True,
        scope='Same source function with distinct directed sector enclosures; no resolved point field or global cone claimed.',
        input_hashes=asts.hashes,passed=True)


def compiled_current_entrance_incoming_exporter():
    """Extend current main arithmetic to entrance; change only its guard."""
    asts=SourceAST();original=asts.method('current_pulse_main_exit_background_tensor','chart')
    fn=copy.deepcopy(original);fn.decorator_list=[];changes=[]
    replacements={
        "chart not in ('pulse_main', 'pulse_exit')":"chart!='pulse_entrance'",
        "limits = (c.mpf('.02'), c.mpf(10)) if chart == 'pulse_main' else (c.mpf(10), c.mpf(11))":"limits=(c.mpf(0),c.mpf('.02'))"}
    for node in ast.walk(fn):
        if isinstance(node,ast.If) and ast.unparse(node.test) in replacements:
            old=ast.unparse(node.test);node.test=ast.parse(replacements[old],mode='eval').body;changes.append(old)
        elif isinstance(node,ast.Assign) and ast.unparse(node) in replacements:
            old=ast.unparse(node);replacement=ast.parse(replacements[old]).body[0]
            node.value=replacement.value;changes.append(old)
    if set(changes)!=set(replacements):raise ValueError('Current entrance/incoming source extension anchors differ: '+str(changes))
    fn.name='current_extended_source';env=dict(vars(main))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<current source extension>','exec'),env)
    return env[fn.name],dict(exact_reviewed_source_extension_statements=changes,
        all_other_current_source_stress_velocity_logs_and_physical_operators_unchanged=True,
        whole_entrance_original_forward_integrals_and_pressure_unchanged=True,input_hashes=asts.hashes)


def normalized_actual_power_theorem():
    """The pre raw equations convert exactly to the pure-power shape rows."""
    y,z,mu=s.symbols('y Z mu',real=True);C=1/(1+z*z);bp=s.Rational(1,2)+mu
    u,m,h,k,e,p=[s.Function(name)(y,z) for name in ('u','m','h','k','e','Pabsolute')]
    invP=s.Symbol('same_invPstar',positive=True)
    rules={s.diff(u,y):-bp*u,s.diff(m,y):-m,s.diff(h,y):u-s.Rational(3,2)*h,
        s.diff(k,y):-s.Rational(3,2)*k,s.diff(e,y):-e-u*u/2,s.diff(p,y):u*u/2}
    values=dict(M=m*invP/u,N=k*invP/(u*u),E=e/(u*u),X=h/u,P=p*C*C/(u*u))
    rhs=dict(M=-(s.Rational(1,2)-mu)*values['M'],N=-(s.Rational(1,2)-2*mu)*values['N'],
        E=2*mu*values['E']-s.Rational(1,2),X=1-(1-mu)*values['X'],
        P=(1+2*mu)*values['P']+C*C/2)
    checks={}
    for key,value in values.items():
        if s.cancel(s.diff(value,y).xreplace(rules)-rhs[key])!=0:raise ArithmeticError('Actual power shape normalization failed: '+key)
        checks[key+'_source_normalized_ODE']=True
    asts=SourceAST();asts.method('pre_pulse_mixed_C4','physical_mixed');asts.method('global_physical_assembly','normalized_sources')
    for target,wanted in (
        ('m1',"copy_jet(c,raw['m'][0])*invP/u"),('m2',"copy_jet(c,raw['k'][0])*invP/(u*u)"),
        ('energy',"copy_jet(c,raw['e'][0])/(u*u)"),
        ('pressure',"(copy_jet(c,raw['p'][0])+datum)*C*C/(u*u)")):
        asts.expression('current_pulse_entrance_incoming_background_tensor','actual_power',target,wanted=wanted)
    Pstar,R,U=s.symbols('Pstar R current_scalar_U',positive=True)
    raw_m,raw_h,raw_k,raw_e,raw_p=s.symbols('raw_m raw_h raw_k raw_e raw_absolute_p',real=True)
    uu=U*C
    native=dict(theta=Pstar*uu,Mz=R*raw_m,Mtheta=s.sqrt(2)*R**s.Rational(3,2)*Pstar*raw_h,
        Mtheta_z=s.sqrt(2)*R**s.Rational(3,2)*Pstar*raw_k,Mztheta=R*Pstar**2*raw_e,P=Pstar**2*raw_p)
    transformed=dict(theta=Pstar*U*C,Mz=R*Pstar*U*C*(raw_m/(Pstar*uu)),
        Mtheta=s.sqrt(2)*R**s.Rational(3,2)*Pstar*U*C*(raw_h/uu),
        Mtheta_z=s.sqrt(2)*R**s.Rational(3,2)*(Pstar*U*C)**2*(raw_k/(Pstar*uu**2)),
        Mztheta=R*(Pstar*U*C)**2*raw_e/uu**2,P=(Pstar*U)**2*raw_p*C*C/uu**2)
    for key in native:
        if s.cancel(native[key]-transformed[key])!=0:raise ArithmeticError('Actual BASE power unit differs: '+key)
        checks['actual_BASE_and_normalized_full_'+key+'_source_units_identical']=True
    rr,lp,Tw=s.symbols('logRref logPstar Tw',real=True)
    env=dict(self=SimpleNamespace(logRref=rr,logP=lp,params=SimpleNamespace(Tw=Tw)))
    logRp=asts.evaluate(asts.expression('global_physical_assembly','__init__','self.logRp',wanted='self.logRref+self.logP+1+self.params.Tw'),env)
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRref=rr,logP=lp,logRp=logRp,params=SimpleNamespace(Tw=Tw,mu=mu))
    if s.cancel(radius(owner,'O3_power',s.Integer(1),{},None)[0]-radius(owner,'pulse_entrance',s.Integer(0),{},None)[0])!=0:
        raise ArithmeticError('Actual incoming/entrance production radius differs')
    checks['actual_BASE_O3_phase1_and_entrance_y0_same_logRp']=True
    return dict(identities=checks,actual_raw_pre_equations_consumed=True,
        raw_energy_rate_minus1_becomes_normalized_rate_2mu_after_dividing_by_u_squared=True,
        all_positive_actual_U_normalizations_are_source_functions_not_caps=True,input_hashes=asts.hashes,passed=True)


def read_producer():return json.loads(gzip.decompress((HERE/NAME).read_bytes()))


class CurrentPulseEntranceIncomingBackgroundTensor:
    @source_precision
    def __init__(self,main_tensor=None,require_checked=True):
        self.main_tensor=main_tensor if main_tensor is not None else CurrentPulseMainExitBackgroundTensor()
        if not self.main_tensor.acceptance_loaded:raise ValueError('Checked current main/exit tensor required')
        self.physical=self.main_tensor.physical;self.pulse=self.main_tensor.pulse;self.ctx=self.main_tensor.ctx
        self.family=self.main_tensor.family;self.source=self.main_tensor.source;self.datum_sha=self.main_tensor.datum_sha
        self.assert_graph();self.exporter,self.extension=compiled_current_entrance_incoming_exporter()
        if not all(self.physical.operator_bindings.values()) or self.physical.radius.__func__ is not BASE.radius:raise ValueError('Original current BASE radius/normalization/operator required')
        current=self.main_tensor.proof['original_generic_full_main_exit_physical_and_exit_gap_theorem']
        source=generic_entrance_proof('pulse_entrance_similarity_C4','entrance_source_proof',current)
        physical=generic_entrance_proof('pulse_entrance_physical_C2','entrance_physical_binding',current)
        production=production_power_function_theorem(self);normalization=normalized_actual_power_theorem()
        velocity=actual_power_velocity_source_theorem()
        self.proof=dict(original_generic_entrance_source_theorem=source,original_generic_entrance_physical_theorem=physical,
            actual_both_production_O2_O3_function_equality=production,current_main_source_theorem=self.main_tensor.proof,
            current_source_extension=self.extension,current_actual_raw_power_normalization_theorem=normalization,
            actual_pre_and_full_tensor_source_velocity_theorem=velocity,current_selected_graph=self.physical.assert_graph(),
            actual_entrance_main_same_source_function=True,actual_incoming_entrance_full_tensor_source_function_identified=True,
            whole_O3_power_tensor_constructed_directly_from_actual_pre_five_moments_and_pressure=True,
            current_absolute_pressure_and_full_nonzero_incoming_energy_retained=True,
            source_function_equality_precedes_common_bounds=True,passed=True)
        self.hashes=dict(self.main_tensor.hashes)
        for value in (source,physical,production,normalization,velocity,self.extension):self.hashes.update(value['input_hashes'])
        for stem in ('pulse_entrance_similarity_C4','pulse_entrance_physical_C2','current_pulse_entrance_incoming_background_tensor'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):raise ValueError('Current entrance/incoming receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        self.main_tensor.assert_graph()
        if not (self.main_tensor.acceptance_loaded and self.physical is self.main_tensor.physical and self.pulse is self.main_tensor.pulse and self.ctx is self.physical.ctx):raise ValueError('Current entrance/incoming must share checked main/exit graph')

    def actual_power(self,Z,phase,log_tau,theta,viscosity):
        c=self.ctx;Z=c.mpf(Z);mu=c.mpf(self.pulse.mu);delta=c.mpf(self.pulse.delta)
        lt=c.mpf(log_tau);nu=c.mpf(viscosity)
        if not all(mp.isfinite(x) for x in endpoints(Z)+endpoints(lt)+endpoints(nu)) or endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(nu)[0]<=0:raise ValueError('Finite source Z/time and positive viscosity required')
        if theta is not None and not all(mp.isfinite(x) for x in endpoints(c.mpf(theta))):raise ValueError('Finite theta required')
        pre=self.physical.pre.power(Z,phase);scalar=self.physical.pre.power(0,phase)
        raw=pre['actual_normalized_primitive_y_derivative_axial5']
        u=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in pre['Utheta_over_Pstar_axial5_coefficients']])
        if endpoints(u[0])[0]<=0:raise ValueError('Actual pre U positivity required; no cap division')
        z=IntervalTaylor.variable(c,Z,5);C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal();zero=C*0;zeros=[zero]*5
        invP=c.sqrt(c.mpf(self.physical.pre.invP2))
        m1=copy_jet(c,raw['m'][0])*invP/u
        m2=copy_jet(c,raw['k'][0])*invP/(u*u)
        energy=copy_jet(c,raw['e'][0])/(u*u)
        datum=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in pre['original_P0_axial5_coefficients']])
        pressure=(copy_jet(c,raw['p'][0])+datum)*C*C/(u*u)
        if endpoints(c.mpf(endpoints(scalar['Utheta_over_Pstar_axial5_coefficients'][0])))[0]<=0:raise ValueError('Actual scalar U positivity required')
        X=c.mpf(endpoints(scalar['actual_normalized_primitive_y_derivative_axial5']['h'][0][0]))/c.mpf(endpoints(scalar['Utheta_over_Pstar_axial5_coefficients'][0]))
        m=[m1*(-(c.mpf('.5')-mu))**j for j in range(5)];n=[m2*(-(c.mpf('.5')-2*mu))**j for j in range(5)]
        e=[energy];P=[pressure]
        for j in range(4):e.append(e[j]*(2*mu)-(c.mpf('.5') if j==0 else 0));P.append(P[j]*(1+2*mu)+(C*C/2 if j==0 else zero))
        rows=pulse_coefficients(delta,mu,z,C,X,zeros,m,n,e,zeros,P)
        logR=self.physical.radius('O3_power',phase,pre,self.physical.pre)[0]
        logB=dict(logPstar=self.physical.logP,actual_log_current_U=c.mpf(scalar['log_Utheta_over_Pstar_base_source']))
        sectors={}
        for label,parts in rows.items():
            sectors[label]={}
            for name,part in parts.items():
                rp,bp,dp,hp=part['mode']
                grid={'s%d_Z%d'%(j,n):jet[n]*math.factorial(n) for j,jet in enumerate(part['full_derivative_rows']) for n in range(4-j)}
                sectors[label][name]=dict(mode=part['mode'],extra_source='one',
                    exact_source_log_parts=dict(source_logR=rp*logR,**{key:bp*value for key,value in logB.items()},
                        signed_original_memory_log=c.mpf(0),exact_source_history_log=c.mpf(0),normalization=-c.ln(2)/2),
                    full_stress_mixed3_coefficient_enclosures=grid)
        packet=dict(Z=Z,xi=-mu*c.mpf(self.physical.pre.params.Tw)*(1-phase),full_meridional_stress_log_sectors=sectors,
            exact_source_logs=dict(R=logR,B=logB,H=c.mpf(0),D0=c.mpf(0),extra=dict(one=c.mpf(0),incoming1=c.mpf(0))))
        mapped=main_exit_to_physical_packet(packet,mu);velocity=pulse_velocity_rows(c,delta,mu,z,C,zeros,m)
        point=lift_physical_packet(c,mapped,delta,velocity,lt,theta,nu)
        return dict(point,current_actual_source_stress_packet=mapped,current_actual_normalized_power_rows=dict(m=m,n=n,e=e,P=P,X=X),
            actual_upstream_original_pre_power_source=pre,actual_upstream_scalar_U_X_defining_source=scalar,
            actual_upstream_physical_spatial4_time1_packet=self.physical.evaluate('O3_power',Z,phase,log_tau,theta),
            current_source_three_component_velocity_rows=velocity,actual_full_stress_not_local_difference=True,
            actual_full_upstream_power_not_only_unit_left_neighborhood=True,actual_pre_pressure_function_used_in_tensor=True,
            source_bounds_not_resolved_physical_point_values=True,**dict.fromkeys(OPEN,False))

    @source_precision
    def chart(self,chart,Z,coordinate,log_tau='-1',theta=None,viscosity='1'):
        self.assert_graph();c=self.ctx;v=c.mpf(coordinate);mu=c.mpf(self.pulse.mu)
        if chart=='early_y':
            if endpoints(v)[0]<0 or endpoints(v)[1]>1:raise ValueError('Early original y[0,1] required')
            xi=mu*v;internal='pulse_entrance'
        elif chart=='O3_power':
            if endpoints(v)[0]<0 or endpoints(v)[1]>1:raise ValueError('Whole original O3 power phase[0,1] required')
            if not all(mp.isfinite(x) for x in endpoints(v)):raise ValueError('Finite power phase required')
            result=self.actual_power(Z,v,log_tau,theta,viscosity)
            result.update(chart=chart,coverage_coordinate=v,all_nonzero_incoming_histories_and_radial_remainder_retained=True)
            return result
        elif chart=='pulse_entrance':xi=v;internal=chart
        else:raise ValueError('Current entrance or actual O3 power chart required')
        result=self.exporter(self.main_tensor,internal,Z,xi,log_tau,theta,viscosity)
        data=self.main_tensor._data(c.mpf(Z));Bh=result['source_full_main_exit_history_rows']['Bh']
        if internal=='pulse_entrance':
            partial=gp_energy(c,xi,None,self.main_tensor.cells)
            forward=forward_entrance_energy_rows(c,mu,xi,data['ap'],data['incoming_energy'],partial,Bh)
            result.update(original_partial_entrance_energy_bound=partial,current_forward_anchored_full_energy_rows=forward,
                selected_forward_backward_energy_same_source_not_two_added_terms=True)
        result.update(chart=chart,coverage_coordinate=v,source_xi=xi,
            original_ordinary_y_coordinate='mu*early_y or xi=mu*y for entrance;xi=-mu*Tw*(1-phase) for O3 power',
            all_nonzero_incoming_histories_and_radial_remainder_retained=True,
            actual_full_upstream_power_not_only_unit_left_neighborhood=internal=='O3_power',**dict.fromkeys(OPEN,False))
        return result

    @source_precision
    def interface(self,name,Z=(-1,1),log_tau=('-3','-1'),theta=None,viscosity='1'):
        if name=='entrance_main':left=self.chart('pulse_entrance',Z,'.02',log_tau,theta,viscosity);right=self.main_tensor.chart('pulse_main',Z,'.02',log_tau,theta,viscosity)
        elif name=='incoming_entrance':left=self.chart('O3_power',Z,1,log_tau,theta,viscosity);right=self.chart('pulse_entrance',Z,0,log_tau,theta,viscosity)
        else:raise ValueError('Current incoming/entrance or entrance/main seam required')
        a=canonical_tensor_groups(left);b=canonical_tensor_groups(right);common={}
        if set(a)!=set(b):raise ValueError('Actual current entrance/incoming tensor layouts differ')
        for key in a:
            values=[endpoints(row['log_absolute_upper'])[1] for row in a[key]+b[key] if not row['exact_zero']]
            common[key]=dict(exact_zero=not values,log_absolute_upper=self.ctx.mpf(max(values))+self.ctx.ln(len(values)) if values else None)
        return dict(seam=name,common_actual_tensor_rows=common,current_common_tensor_contribution_count=len(common),
            current_actual_both_production_and_full_tensor_function_join_theorem=self.proof,
            source_function_equality_precedes_common_triangle_bounds=True,interval_overlap_not_used_as_function_identity=True,**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            current_actual_entrance_incoming_source_and_tensor_join_theorem=self.proof,
            actual_current_tensor_regions_available=['O3_power','pulse_entrance']+self.main_tensor.manifest()['actual_current_tensor_regions_available'],
            actual_current_completed_tensor_adjacent_interface_count=15,actual_current_completed_tensor_internal_interface_count=4,
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            entrance_domain='xi[0,.02],y[0,.02/mu]; early y[0,1]',incoming_domain='whole O3 power phase[0,1],y=-Tw*(1-phase)',
            producer_storage='Deterministic gzip of complete unpruned UTF-8 JSON',
            remaining_pre_core_axis_angular_internal_global_and_temporal_scope_open=True,
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPulseEntranceIncomingBackgroundTensor(require_checked=False)
    result=field.manifest();result['current_actual_entrance_incoming_tensor_views']={}
    for name,args in VIEWS.items():
        result['current_actual_entrance_incoming_tensor_views'][name]=field.chart(*args)
        print('Current actual entrance/incoming tensor: '+name,flush=True)
    result['current_actual_two_entrance_incoming_tensor_interfaces']={name:field.interface(name) for name in SEAMS}
    payload=(json.dumps(encode(pack(result)),separators=(',',':'))+'\n').encode()
    (HERE/NAME).write_bytes(gzip.compress(payload,compresslevel=6,mtime=0));return result


if __name__=='__main__':run()
