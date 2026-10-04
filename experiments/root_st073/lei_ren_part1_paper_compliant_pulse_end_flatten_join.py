"""Actual original pulse-end / flatten functional moment and stress join.

This receipt closes one similarity interface. Full meridional physical
transfer, pulse admissibility, global bounds and temporal recursion stay open.
"""
import ast
import copy
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets, sigma_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import backward_bump_weights
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(pulse_s=0,flatten_t=0,Z=[-1,1],ordinary_derivative='d_s=d_t=d_logR')
FALSE_FLAGS=('pulse_end_physical_decomposition_constructed','pulse_end_cone_certified',
    'whole_pulse_stress_constructed','global_admissible_stress_lift_constructed',
    'physical_energy_integral_certified','temporal_recursion')


class ExactSourceHalves(ast.NodeTransformer):
    """Preserve exactly representable source halves in symbolic arithmetic.

    Python 0.5 and 1/2 are exact binary halves. SymPy Float arithmetic can
    round large expanded coefficients, so their algebraic replay uses 1/2.
    No parameter, enclosure, selected coefficient or field is rationalized.
    """
    def visit_Constant(self,node):
        if isinstance(node.value,float):
            if abs(node.value)!=.5:raise ValueError('Unexpected symbolic source float')
            return ast.copy_location(ast.Call(func=ast.Name(id='source_Rational',ctx=ast.Load()),
                args=[ast.Constant(value=1 if node.value>0 else -1),ast.Constant(value=2)],keywords=[]),node)
        return node
    def visit_BinOp(self,node):
        node=self.generic_visit(node)
        if isinstance(node.op,ast.Div) and isinstance(node.left,ast.Constant) and isinstance(node.right,ast.Constant):
            if node.left.value==1 and node.right.value==2:
                return ast.copy_location(ast.Call(func=ast.Name(id='source_Rational',ctx=ast.Load()),
                    args=[ast.Constant(value=1),ast.Constant(value=2)],keywords=[]),node)
        return node


def current_sources():
    records={};hashes={};family=None
    for stem in ('pulse_end_stress_C3','flatten_stress_C3','flatten_cone','power_inlet_C4','flatten_mixed_C4','axial_high_jets'):
        name=PREFIX+stem+'_check.json';raw=(HERE/name).read_bytes();check=json.loads(raw)
        if not check['all_passed']:raise ValueError('Current interface prerequisite failed: '+stem)
        pair=(check['actual_five_defect_family_sha256'],check['implicit_source_sha256'])
        if family is None:family=pair
        elif pair!=family:raise ValueError('Interface source family differs: '+stem)
        for source,digest in check['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Interface source changed: '+source)
        hashes.update(check['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes()
        records[stem]=json.loads(raw);records[stem+'_check']=check;hashes[name]=hashlib.sha256(raw).hexdigest()
    return records,hashes,family



def local_inlet_and_datum_source_binding(records):
    """Replay the two local jet producers; do not infer equality from boxes."""
    asts=SourceAST();proofs={}
    z=s.symbols('Z',real=True);U,Pin,Hp=s.symbols('same_U same_Pin same_Hp',real=True)
    datum=s.Function('same_analytic_pressure_datum')(z)
    c=SimpleNamespace(mpf=s.sympify)
    coefficients=[s.diff(datum,z,n)/s.factorial(n) for n in range(6)]
    provider=SimpleNamespace(normalized_jets=lambda Z,order:dict(normalized_pressure_coefficients=coefficients))
    constants={name:s.symbols('same_'+name,real=True) for name in ('M','K','E_Z','E_Q','C1','C2','C_E','C0')}
    constants['U']=U
    gate=records['power_inlet_C4_check']['exact_functional_production_and_join_identities']
    for name in ('canonical_Xp_constant','canonical_Mp_constant',
        'actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes',
        'single_sample_only_encloses_proved_Z_independent_Hp_Pin_constants',
        'source_identity_is_distinct_from_numeric_overlap'):
        if not gate[name]:raise ValueError('Actual shared inlet constants not bound: '+name)
        proofs['consumed_'+name]=True
    if not records['axial_high_jets_check']['actual_incoming_C1_source_agreement_checked']:
        raise ValueError('Original live and serialized inlet constants are not source-bound')
    proofs['consumed_actual_incoming_C1_source_agreement_checked']=True
    # Both U constants are the same fifth provider fourth.constants object.
    # Its producer serializes that defining object in selected.incoming.
    asts.expression('pulse_radial_C4','__init__','self.high',
        wanted='SimpleNamespace(ctx=c,base=self.fifth.fourth.base,constants=self.fifth.fourth.constants,select=self.fifth.select,energy=SimpleNamespace(future=self.fifth.future))')
    asts.expression('fifth_axial_jets','select','k',wanted='self.fourth.constants')
    asts.expression('power_inlet_C4','__init__','self.constants',
        wanted="{k:read_interval(c,v) for k,v in fifth['whole_Z']['selected']['incoming']['Z_independent_constant_definitions'].items() if isinstance(v,dict) and 'lower' in v}")
    fifth=asts.method('fifth_axial_jets','select')
    serialized=[n.value for n in ast.walk(fifth) if isinstance(n,ast.keyword) and n.arg=='Z_independent_constant_definitions']
    if len(serialized)!=1 or ast.dump(serialized[0])!=ast.dump(ast.Name(id='k',ctx=ast.Load())):
        raise ValueError('Fifth provider no longer serializes its actual fourth constants')
    proofs['actual_fifth_serializes_same_live_constants_consumed_by_power']=True
    # Bind the common Pin and U symbols to BOTH actual producer routes.
    # O2/O3 canonical shape gates above bind the slope-mu inlet functions.
    mu,Tw,Uw,Pw=s.symbols('same_mu same_Tw same_slope_inlet_U same_slope_inlet_Pin',real=True)
    bc=SimpleNamespace(mpf=s.Rational,exp=s.exp)
    bq=1+z*z;env=dict(c=bc,mu=mu,t=Tw,
        get=lambda key:{'Utheta_over_Pstar':Uw/bq,'Mp_over_Pstar_squared':Pw/bq**2}[key],
        decay_integral=lambda ctx,rate,length:s.Function('same_original_decay_integral')(rate,length))
    for name in ('slope','f','u1','u','p'):
        env[name]=asts.evaluate(asts.expression('outer_buffer','power',name),env)
    sourcePin0=env['p'].subs(z,0);sourcePinHalf=env['p'].subs(z,s.Rational(1,2))
    pulsePin=asts.evaluate(asts.expression('pulse_radial_C4','__init__','self.inlet_P',
        wanted="box(p0['Mp_over_Pstar_squared'][0])"),dict(box=lambda value:value,p0=dict(Mp_over_Pstar_squared=[sourcePin0])))
    sample=dict(Mp_over_Pstar_squared=[sourcePinHalf])
    constructor_sample=asts.evaluate(asts.expression('power_inlet_C4','__init__','sample',
        wanted="buffer['samples'][-1]"),dict(buffer=dict(samples=[sample])))
    original_q=asts.evaluate(asts.expression('power_inlet_C4','__init__','q',wanted="c.mpf('1.25')"),dict(c=bc))
    powerPin=asts.evaluate(asts.expression('power_inlet_C4','__init__','self.inlet_P',
        wanted="read_interval(c,sample['Mp_over_Pstar_squared'][0])*q*q"),
        dict(read_interval=lambda ctx,value:value,c=bc,sample=constructor_sample,q=original_q))
    if s.cancel(pulsePin-powerPin)!=0:raise ArithmeticError('Actual live/serialized Pin definitions differ')
    proofs['actual_two_Pins_equal_from_buffer_Z0_and_serialized_Zhalf_producers']=True
    asts.expression('pulse_radial_C4','__init__','p0',wanted='self.pulse.buffer.power(0,1)')
    asts.expression('axial_high_jets','_incoming_constants','p0',wanted='buffer.power(0,1,cells=128)')
    liveU=asts.evaluate(asts.expression('axial_high_jets','_incoming_constants','U',
        wanted="box(p0['Utheta_over_Pstar'][0])"),dict(box=lambda value:value,p0=dict(Utheta_over_Pstar=[env['u'].subs(z,0)])))
    if s.cancel(liveU-env['u']*(1+z*z))!=0:raise ArithmeticError('Actual common U definition differs')
    proofs['actual_live_U_equals_same_original_buffer_q_times_swirl']=True
    asts.expression('outer_buffer','report','samples',
        wanted="[self.power('.5',t) for t in ('0','.5','1')]",augmented=True)
    metadata=[]
    for stem in ('outer_buffer','five_moment_repair','power_inlet_C4'):
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();data=json.loads(raw)
        asts.hashes[name]=hashlib.sha256(raw).hexdigest()
        metadata.append((data['implicit_source_sha256'],data['datum_enclosure_sha256']))
        if stem=='outer_buffer':
            last=data['samples'][-1]
            if not last['pulse_inlet'] or mp.mpf(last['Z']['lower'])!=mp.mpf('.5') or mp.mpf(last['Z']['upper'])!=mp.mpf('.5'):
                raise ValueError('Actual serialized terminal buffer sample changed')
    if len(set(metadata))!=1:raise ValueError('Original analytic source/datum hashes differ')
    asts.expression('five_moment_repair','__init__','self.datum',wanted="LogarithmicPressureDatum('40',precision=160)")
    asts.expression('power_inlet_C4','__init__','self.datum',wanted="CompliantPressureDatum('40',precision=160)")
    asts.expression('outer_initial','__init__','self.datum',wanted='self.repair.datum')
    proofs['actual_live_and_serialized_analytic_datum_source_and_enclosure_hashes_equal']=True
    path=HERE/'lei_ren_part1_paper_candidate_pressure_function.py'
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    tree=ast.parse(path.read_text(encoding='utf8'));qfn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='q_jets')
    ns={};exec(compile(ast.Module(body=[qfn],type_ignores=[]),'<actual reciprocal-square coefficients>','exec'),ns)
    pulse=SimpleNamespace(inlet_P=Pin,inlet_H=Hp,high=SimpleNamespace(constants=constants))
    pe=dict(c=c,Z=z,IntervalTaylor=IntervalTaylor,q_jets=ns['q_jets'],self=pulse,k=constants)
    for name in ('z','q','r','z2','u','inlet'):
        pe[name]=asts.evaluate(asts.expression('pulse_radial_C4','data',name),pe)
    power=SimpleNamespace(constants=constants,inlet_P=Pin,Xp=Hp/U,datum=provider)
    fe=dict(c=c,Z=z,IntervalTaylor=IntervalTaylor,endpoints=lambda value:value,self=power)
    for name in ('z','q','invq','k','u','m','e','raw','p0'):
        fe[name]=asts.evaluate(asts.expression('power_inlet_C4','incoming',name),fe)
    fn=asts.method('power_inlet_C4','incoming')
    returns=[n.value for n in ast.walk(fn) if isinstance(n,ast.Return)]
    if len(returns)!=1:raise ValueError('Actual power incoming return changed')
    result=asts.evaluate(returns[0],fe)
    # Replay the actual pulse datum getter and conversion. The source family
    # and canonical pressure theorem bind this analytic datum, not a fitted jet.
    pulse.selection=SimpleNamespace(future=SimpleNamespace(angular=SimpleNamespace(initial=SimpleNamespace(datum=provider))))
    pe['self']=pulse;pe['endpoints']=lambda value:value
    pe['rows']=asts.evaluate(asts.expression('pulse_radial_C4','pressure_moment','rows'),pe)
    pulseP0=asts.evaluate(asts.expression('pulse_radial_C4','pressure_moment','p0'),pe)
    def zero(name,value):
        if s.cancel(value)!=0:raise ArithmeticError('Actual local inlet source differs: '+name)
        proofs[name]=True
    for n in range(6):
        zero('actual_two_U_over_q_local_axial_coefficient'+str(n),pe['u'][n]-result['u'][n])
        zero('actual_two_Pin_over_q_squared_local_axial_coefficient'+str(n),
            pe['inlet']['Mp_over_Pstar_squared'][n]-result['Mp'][n])
        zero('actual_two_same_analytic_datum_local_axial_coefficient'+str(n),pulseP0[n]-result['P0'][n])
        zero('actual_pulse_q_jets_are_derivatives_of_Pin_over_q_squared'+str(n),
            pe['inlet']['Mp_over_Pstar_squared'][n]-s.diff(Pin/(1+z*z)**2,z,n)/s.factorial(n))
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        local_axial_orders=list(range(6)),interval_overlap_used_as_source_proof=False,
        same_canonical_whole_Z_constants_consumed=True,
        original_Pin_defining_source=str(pulsePin),original_U_defining_source=str(liveU),
        shared_analytic_datum_source_and_enclosure_hashes=list(metadata[0]))

def source_endpoint_binding(records):
    asts=SourceAST();proofs={}
    def syntax(stem,method,target,wanted=None,augmented=False):
        return asts.expression(stem,method,target,wanted=wanted,augmented=augmented)
    def zero(name,value):
        if s.cancel(s.expand_power_exp(value))!=0:raise ArithmeticError('Endpoint source differs: '+name)
        proofs[name]=True
    groups=(
        (records['pulse_end_stress_C3']['reference_amplitude_bridge']['identities'],
         ('actual_pulse_B0_equals_C0_times_same_physical_flatten_reference',
          'actual_absolute_pressure_scale_equals_same_physical_Ev0_squared',
          'consumed_actual_flatten_radius_and_right_power_interface_AST_verified')),
        (records['flatten_stress_C3']['source_bridge']['identities'],
         ('actual_flatten_future_energy_half_of_complete_C5_source_AST_verified',
          'actual_complete_future_contains_selected_outer_angular_repairs_and_entire_heat',
          'same_absolute_pressure_identified_by_original_FTC_and_power_datum',
          'actual_flatten_energy_AST_equals_same_power_inlet_Z5')),
        (records['flatten_cone']['source_bridge']['identities'],
         ('actual_original_pulse_Xp_ratio_and_exact_terminal_memory_AST_verified',
          'actual_original_Xv_is_Z_independent_function_not_chosen_box_value',
          'actual_original_sigma_flat_endpoints_and100_unit_scaling_verified',
          'exact_positive_S_source_consumed')))
    for group,names in groups:
        for name in names:
            if not group[name]:raise ValueError('Actual interface source binding missing: '+name)
            proofs['consumed_'+name]=True
    mu,z=s.symbols('mu Z',real=True);r=1-mu;c=SimpleNamespace(exp=s.exp,mpf=s.Rational)
    fu=s.Function('same_complete_future_half')(z)
    nativeE=syntax('axial_pulse_field','end','e',
        'future*c.exp(2*self.mu*s)+baseline-end_energy*(c.exp(2*self.mu*s)*self.E2cap)')
    pulseE=asts.evaluate(nativeE,dict(future=fu,c=c,self=SimpleNamespace(mu=mu,E2cap=s.Symbol('positive_D_squared')),
        s=s.Integer(0),baseline=s.Integer(0),end_energy=s.Integer(0)))
    flatE=asts.evaluate(syntax('flatten_mixed_C4','flatten','energy',
        '(ev-Eint/2)*c.exp(2*self.mu*t)/(F*F)'),
        dict(ev=fu,Eint=s.Integer(0),F=s.Integer(1),c=c,self=SimpleNamespace(mu=mu),t=s.Integer(0)))
    zero('actual_native_terminal_future_energy_functions_equal',pulseE-flatE)
    syntax('axial_pulse_field','end','future',
        "self.selection.future.future(Z)['complete_future_energy_Taylor']/2")
    syntax('flatten_mixed_C4','future_energy','packet',"self.fifth['whole_Z']['energy']")
    syntax('pulse_radial_C4','__init__','self.high',
        'SimpleNamespace(ctx=c,base=self.fifth.fourth.base,constants=self.fifth.fourth.constants,select=self.fifth.select,energy=SimpleNamespace(future=self.fifth.future))')
    syntax('pulse_high_jets','future','r','dict(self.high.energy.future(endpoints(self.high.ctx.mpf(Z))))')
    Xp,H=s.symbols('original_Xp exact_positive_terminal_H',real=True)
    Xv=1/r+(Xp-1/r)*H
    node=syntax('axial_pulse_field','end','X',
        '1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)')
    pulseX=asts.evaluate(node,dict(self=SimpleNamespace(rate=r,mu=mu,Xp=Xp,factor=lambda value:H),s=s.Integer(0)))
    flatX=asts.evaluate(syntax('flatten_mixed_C4','flatten','X',
        '(Xint+self.Xv*c.exp(-self.rate*t))/F'),
        dict(Xint=s.Integer(0),F=s.Integer(1),self=SimpleNamespace(Xv=Xv,rate=r),c=c,t=s.Integer(0)))
    zero('actual_native_terminal_signed_angular_history_functions_equal',pulseX-flatX)
    syntax('flatten_mixed_C4','__init__','self.Xv',
        "read_interval(c,pulse['whole_Z_terminal']['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0])")
    # Actual empty support routines: the result is independent of every
    # selected coefficient and normalization, hence uniform in Z.
    iv=mp.iv;iv.dps=100;zeros=[]
    for center in (-3,-1):
        local=iv.mpf(-center)
        b=beta_jets(iv,local/iv.mpf('.15'))
        w=backward_bump_weights(iv,s.Symbol('mu'),s.Symbol('positive_normalization'),local,1)
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in b.coefficients+tuple(w)):
            raise ArithmeticError('Original terminal beta support is not empty')
        zeros.append(dict(center=center,beta_ordinary_orders=list(range(5)),backward_weights_exact_zero=True))
    sig=sigma_jets(iv,0)
    if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in sig.coefficients):
        raise ArithmeticError('Original flatten sigma inlet is not flat')
    asts.method('flat_pulse_derivatives','beta_jets');asts.method('flat_pulse_derivatives','sigma_jets')
    asts.method('axial_pulse_field','backward_bump_weights')
    proofs['actual_empty_future_beta_supports_uniform_in_axial_selected_data']=True
    proofs['actual_sigma_inlet_derivatives0through4_are_exact_zero']=True
    # Replay both ORIGINAL forward pressure histories. The current negative
    # future representation is identified with this datum by the admitted FTC.
    U,Pin,datum=s.symbols('original_U original_Pin analytic_P0',real=True);q=1+z*z
    inletU=U/q;prate=1+2*mu;decay=s.exp(-13*prate/mu)
    kernel_node=syntax('pulse_radial_C4','pressure_moment','kernel')
    if not isinstance(kernel_node,ast.IfExp):raise ValueError('Actual pressure kernel branch changed')
    kernel=asts.evaluate(kernel_node.orelse,dict(decay=decay,self=SimpleNamespace(prate=prate)))
    node=syntax('pulse_radial_C4','pressure_moment','p',
        "IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)")
    pulseMp=asts.evaluate(node,dict(IntervalTaylor=lambda ctx,value:value,c=c,
        inlet=dict(Mp_over_Pstar_squared=Pin/q**2),u=inletU,kernel=kernel))
    flatMp=asts.evaluate(syntax('flatten_mixed_C4','flatten','Mp_v',
        "data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))"),
        dict(data=dict(Mp=Pin/q**2,u=inletU),self=SimpleNamespace(pressure_decay=s.exp(-13/mu-26),prate=prate)))
    for n in range(6):
        zero('actual_forward_Mp_and_absolute_P_endpoint_Z'+str(n),s.diff((pulseMp+datum)-(flatMp+datum),z,n))
    # Bind the current pulse stress companion's canonical pressure getter,
    # rather than merely assuming a common pressure symbol in the row proof.
    a,Cscale=s.symbols('a exact_positive_C0',real=True);pheat=1+2*a
    canonicalP=s.Function('same_full_flatten_pressure_numerator')(z)
    node=syntax('pulse_end_stress_C3','end','P0',
        '-(IntervalTaylor.constant(c,1,5)/(2*self.heat.prate)+rawP)/(scale*scale)')
    p0=asts.evaluate(node,dict(IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.sympify(value)),
        c=c,self=SimpleNamespace(heat=SimpleNamespace(prate=pheat)),rawP=canonicalP-1/(2*pheat),scale=Cscale))
    zero('actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',p0+canonicalP/Cscale**2)
    syntax('pulse_end_stress_C3','end','rawP0',"packet['flatten_future_defect_zeroth_Taylor']['pressure_defect_rows']")
    syntax('pulse_end_stress_C3','end','scale',
        "2*c.mpf(endpoints(self.flatten.terminal(z0)['Kright']))*c.exp(-100*(c.mpf(endpoints(self.heat.a))-mu))")
    # Exact future decomposition plus the already-bound right power datum
    # identifies the canonical flatten remaining energy with native e(0).
    a,KR=s.symbols('a original_Kright',real=True);f=q/2;d=a-mu
    F0=s.Function('same_original_flatten_energy_integral')(z)
    Rest=s.Function('same_full_post_flatten_energy')(z);total=F0+Rest
    native_right=(total/2-F0/2)*s.exp(200*mu)/f**2
    terminalE=2*KR**2*native_right
    node=syntax('flatten_stress_C3','flatten_defect_rows','Ed',
        "terminal['energy_defect_rows'][0]*e+one*(c.expm1(heat.delta*offset)/heat.delta)")
    env=dict(terminal=dict(energy_defect_rows=[terminalE-1/(2*a)]),e=s.exp(-200*a),
        one=s.Integer(1),c=SimpleNamespace(expm1=lambda x:s.exp(x)-1),heat=SimpleNamespace(delta=2*a),offset=s.Integer(-100))
    Ed=asts.evaluate(node,env)
    extra=syntax('flatten_stress_C3','flatten_defect_rows','Ed',
        "kernels['energy']*(terminal['Kright']**2*ek)",True)
    Ed+=asts.evaluate(extra,dict(kernels=dict(energy=F0/f**2),terminal=dict(Kright=KR),ek=s.exp(-200*d)))
    C0=2*KR*s.exp(-100*d)
    zero('actual_canonical_flatten_complete_E_at_inlet_is_C0_squared_times_complete_future_over_q_squared',
        Ed+1/(2*a)-C0**2*total/q**2)
    # Original production radius, replayed without source-log materialization.
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'));env={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    env[target.id]=ast.literal_eval(node.value)
    asts.replay('global_physical_assembly','radius',env)
    origin=s.symbols('logRp',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    zero('actual_pulse_end_flatten_production_log_radius_equal',
        env['radius'](assembly,'pulse_end',0,{},None)[0]-env['radius'](assembly,'flatten',0,{},None)[0])
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        empty_supports=zeros,signed_original_memory_retained=True,
        pressure_is_same_absolute_forward_and_remaining_function=True,
        endpoint_equality_is_functional_not_interval_overlap=True)


def functional_join_identities():
    asts=SourceAST();proofs={}
    a,mu,z,q0,Xp=s.symbols('a mu Z q0 original_Xp',real=True)
    R,B,C0,D,H=s.symbols('R0 B0 C0 D H',positive=True)
    delta=2*a;k=1-a;r=1-mu;bh=s.Rational(1,2)+a;bp=s.Rational(1,2)+mu
    p=1+2*a;pp=1+2*mu;d=a-mu;C=1/(1+z*z)
    fu=s.Function('same_complete_future_half')(z);P0=s.Function('same_signed_pressure_over_B0_squared')(z)
    Xv=1/r+(Xp-1/r)*H
    c=SimpleNamespace(mpf=s.Rational,exp=s.exp,expm1=lambda v:s.exp(v)-1,ln=s.log)
    stub=lambda ctx,coefficients:s.sympify(coefficients[0])
    stub.constant=lambda ctx,value,order:s.sympify(value)
    stub.variable=lambda ctx,value,order:s.sympify(value)
    env=dict(math=math,mp=SimpleNamespace(mpf=s.Rational),IntervalTaylor=stub,
        product_rows=product_rows,shifted_rows=shifted_rows,axial_derivative=lambda v:s.diff(v,z),
        log_taylor=s.log,c=c)
    for stem,name in (('flatten_stress_C3','flatten_f'),('flatten_stress_C3','flatten_shape'),
        ('collar_stress_C3','collar_stress_rows')):
        asts.replay(stem,name,env)
    fn=ExactSourceHalves().visit(copy.deepcopy(asts.method('pulse_end_stress_C3','pulse_coefficients')))
    env['source_Rational']=s.Rational
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual pulse rows with exact binary halves>','exec'),env)
    asts.replay('heat_pressure_C4','pressure_y_rows',env,remove_pressure_context=True)
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,k=k,delta=delta,prate=p,S=s.exp(q0)/R)
    shape=env['flatten_shape'](heat,dict(Kright=C0*s.exp(100*d)/2),s.Integer(-100),
        dict(Z=z,F_Taylor=s.Integer(1),sigma_y_derivatives=[s.Integer(0)]*5,
            actual_original_flatten_right_endpoint=False))
    K=[s.cancel(s.expand_power_exp(v)) for v in shape['K_rows']]
    shape['K_rows']=K;shape['K_defect_rows']=[K[0]-1]+K[1:]
    square=product_rows(K,K);shape['K_squared_defect_rows']=[square[0]-1]+square[1:]
    fullA=K[0]*Xv;fullE=2*C0**2*C**2*fu;fullP=-C0**2*P0
    A=[fullA-1/k];E=[fullE-1/delta];Pr=[fullP-1/(2*p)]
    rows=dict(A=A,E=E,Pr=Pr,Km=shape['K_defect_rows'],Qm=shape['K_squared_defect_rows'],heat=heat)
    for target in ('A','E','Pr'):
        fn=asts.method('flatten_stress_C3','flatten_defect_rows')
        calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)==target+'.append']
        if len(calls)!=1:raise ValueError('Actual flatten derivative recurrence changed: '+target)
        asts.bindings['flatten_stress_C3.flatten_defect_rows.'+target+'.append']=True
        for j in range(4):
            ns=dict(rows,j=j);exec(compile(ast.fix_missing_locations(ast.Module(body=[ast.Expr(value=calls[0])],type_ignores=[])),'<actual endpoint recurrence>','exec'),ns)
            rows[target][-1]=s.cancel(rows[target][-1])
    defects=dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
        K_defect_rows=shape['K_defect_rows'],K_squared_defect_rows=shape['K_squared_defect_rows'])
    flatStress=env['collar_stress_rows'](heat,shape,defects,z,q0)
    e=[fu];P=[P0];X=[Xv]
    for j in range(4):
        e.append(2*mu*e[j]-(s.Rational(1,2) if j==0 else 0))
        P.append(pp*P[j]+(C*C/2 if j==0 else 0))
        X.append((1 if j==0 else 0)-r*X[j])
    # These are the exact source recurrences in pulse_end_stress_C3.end.
    fn=asts.method('pulse_end_stress_C3','end')
    for target in ('e0','P'):
        calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)==target+'.append']
        if len(calls)!=1:raise ValueError('Actual pulse endpoint derivative recurrence changed: '+target)
        asts.bindings['pulse_end_stress_C3.end.'+target+'.append']=True
        actual=[fu if target=='e0' else P0]
        for j in range(4):
            ns=dict(e0=actual,P=actual,mu=mu,prate=pp,C=C,j=j)
            ns['source_Rational']=s.Rational
            node=ExactSourceHalves().visit(copy.deepcopy(calls[0]))
            exec(compile(ast.fix_missing_locations(ast.Module(body=[ast.Expr(value=node)],type_ignores=[])),'<actual pulse endpoint recurrence>','exec'),ns)
        for j,value in enumerate(actual):
            if s.cancel(value-(e if target=='e0' else P)[j])!=0:
                raise ArithmeticError('Actual endpoint source recurrence differs')
    zeros=[s.Integer(0)]*5
    pulse=env['pulse_coefficients'](delta,mu,z,C,Xp,zeros,zeros,zeros,e,zeros,P)
    def zero(name,value):
        reduced=s.cancel(s.expand_power_exp(value))
        if reduced.is_zero is not True:raise ArithmeticError('Functional interface differs: '+name+' '+str(reduced))
        proofs[name]=True
    def mixed(name,left,right,total):
        for j in range(total+1):
            diff=s.cancel(s.expand_power_exp(left[j]-right[j]))
            zero(name+'_y'+str(j),diff)
            for n in range(1,total-j+1):
                zero(name+'_y'+str(j)+'_Z'+str(n),s.diff(diff,z,n))
    # Equality is in physical similarity units, including the different
    # source reference rates on the two sides of the interface.
    mixed('actual_Utheta_mixed4',
        [B/C0*v for v in shifted_rows(K,-bh,4)],
        [B*C*(-bp)**j for j in range(5)],4)
    fullArows=[fullA]+A[1:];fullErows=[fullE]+E[1:]
    mixed('actual_Mtheta_mixed4',
        [s.sqrt(2)*R**s.Rational(3,2)*B/C0*v for v in shifted_rows(fullArows,k,4)],
        [s.sqrt(2)*R**s.Rational(3,2)*B*C*v for v in shifted_rows(X,r,4)],4)
    mixed('actual_Mztheta_mixed4',
        [R*B**2/(2*C0**2)*v for v in shifted_rows(fullErows,-delta,4)],
        [R*B**2*C**2*v for v in shifted_rows(e,-2*mu,4)],4)
    flatP=env['pressure_y_rows'](K,fullP,p,B**2/C0**2*s.exp(p*q0),q0)
    pulseP=[B**2*v for v in shifted_rows(P,-pp,4)]
    mixed('actual_absolute_pressure_mixed4',flatP,pulseP,4)
    # Both source Mp are the same absolute pressure minus the same analytic
    # P0(Z), which is independent of logR. The native endpoint datum is bound
    # separately by source_endpoint_binding, including all axial derivatives.
    mixed('actual_Mp_mixed4',flatP,pulseP,4)
    for label in ('theta','axial'):
        actual=[]
        for j in range(4):
            value=0
            for part in pulse[label].values():
                rp,bpwr,dp,hp=map(s.Rational,part['mode'])
                value+=R**rp*B**bpwr*D**dp*H**hp/s.sqrt(2)*part['full_derivative_rows'][j]
            actual.append(value)
        factor=s.sqrt(R/2)*(B/C0 if label=='theta' else B**2/C0**2)
        mixed('actual_'+label+'_stress_mixed3',[factor*v for v in flatStress[label]],actual,3)
    # All meridional beta jets, remaining linear moments, and selected
    # backward quadratic loss vanish at s0; their source ODEs then give zero
    # derivatives. The exact empty-support binding above admits these zeros.
    for label in ('Mz','Mtheta_z','Ur','Uz'):
        mixed('actual_'+label+'_mixed4',zeros,zeros,4)
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        stress_total_mixed_order=3,pressure_and_moment_total_mixed_order=4,
        arbitrary_shared_axial_future_and_pressure_functions=True,
        signed_original_memory_H_not_discarded=True,
        different_reference_rates_transferred_before_comparison=True,
        exact_binary_source_halves_preserved_in_symbolic_replay=True,
        endpoint_equality_is_functional_not_interval_overlap=True)


def report():
    records,hashes,family=current_sources()
    inlet=local_inlet_and_datum_source_binding(records)
    source=source_endpoint_binding(records);join=functional_join_identities()
    hashes.update(inlet['input_hashes'])
    hashes.update(source['input_hashes']);hashes.update(join['input_hashes'])
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],domain=DOMAIN,
        local_inlet_and_datum_source_binding=inlet,source_endpoint_binding=source,functional_join=join,input_hashes=hashes,
        pulse_end_flatten_full_moment_stress_pressure_functional_join_verified=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,
        **{flag:False for flag in FALSE_FLAGS})


def run():
    result=report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual pulse-end / flatten five-moment, stress3 and pressure4 functional join generated; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
