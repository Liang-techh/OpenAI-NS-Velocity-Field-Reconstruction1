"""Original sigmoid entry: completed physical stress and nonzero remainder.

The original entry chart, SAME full future and analytic pressure datum are
mapped through the unchanged general-K physical transfer. Only regional
decomposition/physical joins are admitted; entry cone and global flatness wait.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import CompliantSteepEntryStressC3,SourceAST
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    CompliantCollarPhysicalC2,physical_collar_identities,physical_source_row,scale_row,
    collar_velocity_operator_coefficients,stress_source_log_parts)
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import flatten_mixed
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket,ZSYM,DSYM
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def steep_entry_divergence_grids(heat,defects,shape,q,Z):
    """Full positive-factor ordinary q rows of Div(T), through mixed2.

    Homogeneous A/E modes disappear analytically. Pressure remains the same
    full current moment, with unit baselines cancelled before enclosure.
    First/second full-factor axial rows require only the current K^2 jets.
    """
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5)
    L=1-z*z*heat.delta; d=1-z*z; K=shape['K_rows']; p=heat.prate
    viscous=[K[j+2]-K[j+1]*p+K[j]*(heat.a*(1+heat.a)) for j in range(3)]
    shear=shifted_rows(viscous,-1,2)
    theta_raw=[-K[j+1]/L+shear[j]*(2*heat.S*c.exp(-q)) for j in range(3)]
    theta=shifted_rows(theta_raw,-heat.a-c.mpf('.5'),2)
    Pd=defects['pressure_defect_rows']; Qm=shape['K_squared_defect_rows']
    axial=[(d*axial_derivative(Pd)-z*Pd*(2*p)+z*Qm[0])/L,
           z*Qm[1]/L,z*(Qm[2]-Qm[1]*p)/L]
    return {label:{'y'+str(j)+'_Z'+str(n):rows[j][n]*math.factorial(n)
                   for j in range(3) for n in range(3-j)} for label,rows in (('theta',theta),('axial',axial))}


class EntryHeatReference:
    def __init__(self,heat,shape):self.actual_heat=heat; self.actual_entry_shape=shape
    def __getattr__(self,key):return getattr(self.actual_heat,key)
    def shape(self,Z,q):return self.actual_entry_shape


class EntryDispatch:
    def __init__(self,native,stress):self.native=native; self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.steep if chart=='steep_entry' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='steep_entry':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.steep_in(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
                    physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
                    original_steep_entry_absolute_pressure_companion_used=True)


class AtEntryPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['steep_entry_similarity_stress_mixed3_factored'],
                    Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtEntryPhaseAssembly:
    def __init__(self,native,phase):self.native=native; self.phase=phase
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic entry map received an unexpected chart')
        return self.native.evaluate('steep_entry',Z,self.phase,**kwargs)


def steep_entry_physical_identities():
    general=physical_collar_identities(); proofs={}
    a,q,z,S=s.symbols('a q Z S',real=True); delta=2*a; k=1-a; p=1+delta; L=1-delta*z*z; d=1-z*z
    K=s.Function('same_original_entry_K')(q); A=s.Function('same_full_A')(q,z)
    E=s.Function('same_full_E')(q,z); P=s.Function('same_full_P')(q,z); b=(1-delta)/2
    Ct=(k*A-b*z*s.diff(A,z)-K)/L+2*S*s.exp(-q)*(s.diff(K,q)-(1+a)*K)
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*P+d*s.diff(P,z))/L
    rules={s.diff(A,q):K-k*A,s.diff(A,q,z):-k*s.diff(A,z),
        s.diff(E,q):delta*E-K*K,s.diff(E,q,z):delta*s.diff(E,z),
        s.diff(P,q):p*P-K*K/2,s.diff(P,q,z):p*s.diff(P,z)}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Entry physical reduction failed: '+name)
        proofs[name]=True
    Dt=-s.diff(K,q)/L+2*S*s.exp(-q)*(s.diff(K,q,2)-p*s.diff(K,q)+a*(1+a)*K)
    Dz=(d*s.diff(P,z)-2*z*(p*P-K*K/2))/L
    zero('actual_full_theta_homogeneous_angular_mode_cancelled',(s.diff(Ct,q)+k*Ct).subs(rules,simultaneous=True)-Dt)
    zero('actual_full_axial_homogeneous_energy_mode_cancelled',(s.diff(Cz,q)-delta*Cz).subs(rules,simultaneous=True)-Dz)
    Pd=P-1/(2*p); Qm=K*K-1
    zero('actual_pressure_and_squared_K_unit_baselines_cancelled',
         (d*s.diff(Pd,z)-2*p*z*Pd+z*Qm)/L-Dz)
    factor_theta=s.exp(-(a+s.Rational(1,2))*q); factor_axial=s.exp(-p*q)
    expected_axial=[Dz,z*s.diff(K*K,q)/L,z*(s.diff(K*K,q,2)-p*s.diff(K*K,q))/L]
    axial=[Dz]
    for j in range(2):axial.append(s.diff(axial[-1],q).subs(rules,simultaneous=True)-p*axial[-1])
    for j in range(3):zero('actual_pressure_full_factor_q'+str(j),axial[j]-expected_axial[j])
    fn=next(n for n in ast.parse(Path(__file__).read_text(encoding='utf8')).body
            if isinstance(n,ast.FunctionDef) and n.name=='steep_entry_divergence_grids')
    fn.body=fn.body[:-1]+[ast.Return(value=ast.parse('dict(theta=theta,axial=axial)',mode='eval').body)]
    env=dict(math=math,IntervalTaylor=SimpleNamespace(variable=lambda ctx,value,order:value),
             axial_derivative=lambda value:s.diff(value,z),shifted_rows=shifted_rows)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual entry divergence>','exec'),env)
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,a=a,S=S,delta=delta,prate=p)
    shape=dict(K_rows=[s.diff(K,q,j) for j in range(5)],K_squared_defect_rows=[s.diff(K*K-1,q,j) for j in range(5)])
    rows=env['steep_entry_divergence_grids'](heat,dict(pressure_defect_rows=Pd),shape,q,z)
    for label,factor,target in (('theta',factor_theta,Dt),('axial',factor_axial,Dz)):
        for j in range(3):
            expected=s.diff(factor*target,q,j)/factor if label=='theta' else expected_axial[j]
            for n in range(3-j):zero('actual_divergence_AST_'+label+'_q'+str(j)+'_Z'+str(n),s.diff(rows[label][j]-expected,z,n))
    F=1-2*(1-delta)*z*z-delta*delta*z**4
    reduced=2*F/L**3*s.diff(K,q)-4*z*z/L**2*s.diff(K,q,2)
    actual=sum(coef.subs({ZSYM:z,DSYM:delta})*s.diff(K,q,index[0])
        for index,coef in collar_velocity_operator_coefficients(0,2).items() if index[1]==0)
    zero('actual_entry_axial_viscosity_from_unchanged_general_K_operator',-actual-reduced)
    return dict(general_K_physical_transfer=general,reduction_identities=proofs,
        original_K_is_independent_of_source_Z=True,
        actual_regional_physical_steep_entry_stress_remainder_identity_verified=True,
        exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
        source_homogeneous_A_E_and_pressure_baselines_cancelled_before_enclosure=True,
        entry_remainder_generally_nonzero=True,globally_flat_remainder_not_inferred=True)


def entry_native_field_join_binding(stress):
    asts=SourceAST(); proofs={}
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    mu,a,Ts,z=s.symbols('mu a Ts Z',real=True); thetaR=s.Symbol('same_actual_thetaR',positive=True)
    obj=SimpleNamespace(thetaR=thetaR,rate=1-mu,bp=s.Rational(1,2)+mu,k=1-a,Ts=Ts)
    env=dict(c=c,self=obj,one=s.Integer(1),math=math)
    obj.thetaS=asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.thetaS'),env)
    theta_left=asts.evaluate(asts.expression('steep_waiting_C4','steep_in','theta'),dict(env,t=s.Integer(1),J=s.Rational(1,2)))
    theta_right=asts.evaluate(asts.expression('steep_waiting_C4','steep_power','theta'),dict(env,t=s.Integer(0)))
    if not stress.power_source.exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Actual entry flat sigmoid endpoint jets required')
    asts.expression('steep_waiting_C4','steep_in','sig',wanted='sigma_jets(c,t)')
    rates_left=asts.evaluate(asts.expression('steep_waiting_C4','steep_in','rates'),dict(env,sig=[s.Integer(1)]+[s.Integer(0)]*4))
    rates_right=asts.evaluate(asts.expression('steep_waiting_C4','steep_power','rates'),env)
    if s.simplify(theta_left-theta_right)!=0 or rates_left!=rates_right:raise ArithmeticError('Native entry/power theta/rates differ')
    asts.expression('power_angular_C4','_packet','mixed',
        wanted='flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    # The current entry bridge binds the original packet -> outer packet
    # wrapper and the full common X/energy/pressure source endpoint histories.
    if not stress.bridge['actual_entry_source_AST_bindings']['steep_waiting_C4.packet_to_actual_outer_packet']:
        raise ValueError('Native entry packet wrapper required')
    fn=asts.method('flatten_mixed_C4','flatten_mixed'); fn.decorator_list=[]; body=[]
    for node in fn.body:
        if isinstance(node,ast.If):continue
        body.append(node)
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='rows' for t in node.targets):
            body.append(ast.Return(value=ast.Name(id='rows',ctx=ast.Load()))); break
    fn.body=body
    env=dict(math=math,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value)),
        **{label:flatten_mixed.__globals__[label] for label in ('UR','UT','UZ','P')})
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual entry native mixed rows>','exec'),env)
    common=[s.Function('same_native_'+name)(z) for name in ('X','energy','absolute_pressure')]
    left=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_left,theta_left,*common,s.Symbol('same_exact_Ev2'))
    right=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_right,theta_right,*common,s.Symbol('same_exact_Ev2'))
    for label in ('UR','UT','UZ'):
        for j in range(5):
            for n in range(5-j):
                if s.simplify(s.diff(left[env[label]][j]-right[env[label]][j],z,n))!=0:raise ArithmeticError('Native entry velocity mixed4 join failed')
                proofs[label+'_q'+str(j)+'_Z'+str(n)]=True
    for stem,name in (('steep_entry_stress_C3','steep_in'),('steep_power_stress_C3','steep_power')):
        asts.expression(stem,name,'mixed[P]',wanted="{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
        asts.expression(stem,name,'fields[P]',wanted='[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
        keywords={n.arg:n.value for n in ast.walk(asts.method(stem,name)) if isinstance(n,ast.keyword)}
        for key,wanted in (('physical_mixed_derivatives_total_order_le4','mixed'),('physical_velocity_and_pressure_y_derivative_Taylor','fields')):
            if ast.dump(keywords[key])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Native entry pressure packet route changed')
    join=stress.bridge['actual_entry_power_formula_join']
    if not join['actual_entry_power_K4_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Same actual K4/stress3/pressure4 right join required')
    bracket=asts.expression('collar_physical_C2','collar','ebracket',wanted='-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)')
    for i in range(3):
        for j in range(3-i):
            if any(sum(index)>4 for index in collar_velocity_operator_coefficients(i,j+2)):
                raise ValueError('Entry remainder exceeds admitted K4 join')
            proofs['actual_same_remainder_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(identities=proofs,source_bindings=asts.bindings,
        actual_native_velocity_mixed4_AST_join_verified=True,
        actual_full_pressure_mixed4_packet_routing_join_verified=True,
        actual_remainder_mixed2_same_AST_operator_and_K4_join_verified=True,
        actual_original_flat_sigma_right_endpoint_consumed=True,
        source_function_equality_not_interval_overlap=True,input_hashes=asts.hashes)


def steep_entry_adapter_binding(stress):
    tree=ast.parse(Path(__file__).read_text(encoding='utf8')); bindings={}
    requirements={('EntryHeatReference','shape'):'self.actual_entry_shape',
        ('EntryDispatch','provider'):"self.stress.steep if chart=='steep_entry' else self.native.provider(chart)",
        ('EntryDispatch','evaluate'):'self.stress.steep_in(Z,coordinate)',
        ('AtEntryPhaseStress','collar'):"self.packet['steep_entry_similarity_stress_mixed3_factored']",
        ('AtEntryPhaseAssembly','evaluate'):"self.native.evaluate('steep_entry',Z,self.phase,**kwargs)",
        ('CompliantSteepEntryPhysicalC2','steep_in'):'-self.stress.steep.wait-self.stress.steep.Ts-2+phase'}
    for (cls,name),expr in requirements.items():
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==name)
        if not any(ast.dump(n)==ast.dump(ast.parse(expr,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Entry physical adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    asts=SourceAST(); fn=asts.method('global_physical_assembly','radius'); env={}
    path=HERE/(PREFIX+'global_physical_assembly.py')
    for node in ast.parse(path.read_text(encoding='utf8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):env[target.id]=ast.literal_eval(node.value)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual entry physical radius>','exec'),env)
    mu,a,Ts,wait,Lrel,origin=s.symbols('mu a Ts wait Lrel logRp',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    assembly=SimpleNamespace(ctx=c,logRp=origin,params=SimpleNamespace(mu=mu))
    outer=SimpleNamespace(Lrel=Lrel); native=SimpleNamespace(outer=outer,Ts=Ts,wait=wait)
    radii={name:env['radius'](assembly,chart,value,{},provider)[0] for name,chart,value,provider in
        (('left','outer_angular',s.Integer(0),outer),('entry0','steep_entry',s.Integer(0),native),
         ('entry1','steep_entry',s.Integer(1),native),('right','steep_power',s.Integer(0),native))}
    if s.simplify(radii['entry0']-radii['left'])!=0 or s.simplify(radii['entry1']-radii['right'])!=0:
        raise ArithmeticError('Actual native entry physical radius connection failed')
    if not stress.bridge['native_angular_entry_left_field_pressure_mixed4_join_consumed']:
        raise ValueError('Native angular-entry source field/pressure C4 join required')
    nativejoin=entry_native_field_join_binding(stress); nativejoin['input_hashes'].update(asts.hashes)
    return dict(source_bindings=bindings,original_entry_phase_is_ordinary_q_no_Ts_multiplier=True,
        actual_native_angular_entry_left_radius_field_pressure_join_verified=True,
        actual_entry_power_right_radius_join_verified=True,
        same_source_K4_stress3_pressure4_and_common_physical_factors_consumed=True,
        physical_entry_power_stress_diagonal_divergence_remainder_join_verified=True,
        actual_native_velocity_pressure_and_remainder_source_join=nativejoin,
        general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar',
        physical_collar_source_never_evaluated_at_negative_offset=True,
        upstream_angular_stress_companion_constructed=False)


class CompliantSteepEntryPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantSteepEntryStressC3(); self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=EntryDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx; self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Entry physical family mismatch')
        self.proof=steep_entry_physical_identities(); self.adapter_proof=steep_entry_adapter_binding(self.stress)
        self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        self.hashes.update(self.adapter_proof['actual_native_velocity_pressure_and_remainder_source_join']['input_hashes'])
        for stem,gate in (('steep_entry_stress_C3_check','actual_original_steep_entry_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
                          ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
                          ('steep_power_physical_C2_check','actual_regional_physical_steep_power_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Entry physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Entry physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entry physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def steep_in(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); phase=c.mpf(phase); self.stress.steep._phase(phase)
        q=-self.stress.steep.wait-self.stress.steep.Ts-2+phase; packet=self.stress.steep_in(z,phase)
        heat=EntryHeatReference(self.stress.heat,packet['steep_entry_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtEntryPhaseStress(packet),assembly=AtEntryPhaseAssembly(self.assembly,phase))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        grids=steep_entry_divergence_grids(heat,packet['steep_entry_future_defect_zeroth_Taylor'],packet['steep_entry_shape'],q,z)
        logtau=c.mpf(log_tau); nu=c.mpf(viscosity); beta=-2-heat.delta
        parts=point['actual_source_log_factors']; logR=point['source_logR_enclosure']; divmixed={}
        for label,factor in (('theta','Qtheta'),('axial','Qz')):
            divmixed[label]={}
            for i in range(3):
                for j in range(3-i):
                    coefficient=physical_bracket(c,grids[label],i,j,z,heat.delta,beta-1)
                    divmixed[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,coefficient,parts[factor],
                        beta-1-i+j*(heat.delta-1),logtau,nu,i+j,c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        cs=point['cos_theta']; sn=point['sin_theta']
        divcart=dict(x=[scale_row(c,divmixed['theta']['r0_z0'],-sn)],y=[scale_row(c,divmixed['theta']['r0_z0'],cs)],z=[divmixed['axial']['r0_z0']])
        point['physical_cylindrical_stress_divergence_mixed2']=divmixed
        point['physical_completed_stress_divergence_cartesian']=divcart
        point['physical_momentum_residual_decomposition_cartesian']={label:[scale_row(c,row,-1) for row in divcart[label]]
            +point['physical_remainder_cartesian'][label] for label in ('x','y','z')}
        for key in ('actual_regional_physical_collar_stress_remainder_identity_verified','Gamma_endpoint_physical_remainder_exact_zero','collar_cone_certified'):
            point.pop(key,None)
        point.update(chart='steep_entry',steep_entry_phase=phase,steep_entry_offset_from_Rtail=q,
                     actual_regional_physical_steep_entry_stress_remainder_identity_verified=True,
                     steep_entry_homogeneous_A_E_and_pressure_baselines_cancelled_before_enclosure=True,
                     steep_entry_power_stress_mixed3_join_verified=True,steep_entry_power_pressure_mixed4_join_verified=True,
                     steep_entry_power_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                     native_angular_entry_left_radius_field_pressure_join_verified=True,upstream_angular_stress_companion_constructed=False,
                     steep_entry_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
                     steep_entry_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        samples=[summary(self.steep_in(z,t,lt,angle,nu)) for z,t,lt,angle,nu in
                 (('0','0','-1','0','1'),('.5','.5','-10','.7','.01'),('-.5','1','-100','1','.7'))]
        whole=summary(self.steep_in([-1,1],[0,1],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Original steep_entry phase[0,1], Z[-1,1], tau>0, r>0; physical stress mixed3, diagonal/divergence/remainder mixed2',
                    samples=samples,whole_steep_entry=whole,steep_entry_physical_identities=self.proof,steep_entry_adapter_binding=self.adapter_proof,
                    actual_regional_physical_steep_entry_stress_remainder_identity_verified=True,
                    steep_entry_power_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                    native_angular_entry_left_radius_field_pressure_join_verified=True,upstream_angular_stress_companion_constructed=False,
                    steep_entry_regional_remainder_exact_zero=False,steep_entry_cone_certified=False,
                    independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepEntryPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual original steep-entry completed physical stress and nonzero axial-viscosity remainder generated',flush=True)
    return result


if __name__=='__main__':run()
