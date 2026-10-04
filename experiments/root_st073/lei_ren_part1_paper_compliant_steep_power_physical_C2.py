"""Original Ts-long power: physical stress, divergence and nonzero remainder.

The original power chart and complete absolute-pressure datum are mapped by
the unchanged general-K physical transfer. Exact exponential modes cancel
homogeneous angular/energy/pressure divergence before interval enclosure.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_power_stress_C3 import CompliantSteepPowerStressC3
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative
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


def steep_power_divergence_grids(heat,terminal,shape,q,Z):
    """Full positive-factor y rows of Div(T), before physical chain rule."""
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5)
    L=1-z*z*heat.delta; d=1-z*z; K=shape['K_rows'][0]; KQ=terminal['KQ']
    ell=shape['remaining']; k=heat.k; p=heat.prate
    theta=[K*k/L*(-c.mpf('1.5'))**j+K*(4*heat.S*c.exp(-q))*(-c.mpf('2.5'))**j for j in range(3)]
    PdQ=terminal['pressure_defect_rows'][0]
    NpQ=-z-z*PdQ*(2*p)+z*(p*KQ**2/3)+d*axial_derivative(PdQ)
    pure=z*K*K*(2*k)/(3*L)
    axial=[pure+NpQ/L*c.exp(-p*ell),-z*K*K*(2*k)/L,z*K*K*(6*k)/L]
    return {label:{'y'+str(j)+'_Z'+str(n):rows[j][n]*math.factorial(n)
                   for j in range(3) for n in range(3-j)} for label,rows in (('theta',theta),('axial',axial))}


class PowerHeatReference:
    def __init__(self,heat,shape):self.actual_heat=heat; self.actual_power_shape=shape
    def __getattr__(self,key):return getattr(self.actual_heat,key)
    def shape(self,Z,q):return self.actual_power_shape


class PowerDispatch:
    def __init__(self,native,stress):self.native=native; self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.steep if chart=='steep_power' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='steep_power':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.steep_power(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
                    physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
                    original_steep_power_absolute_pressure_companion_used=True)


class AtPowerPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['steep_power_similarity_stress_mixed3_factored'],
                    Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtPowerPhaseAssembly:
    def __init__(self,native,phase):self.native=native; self.phase=phase
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic power map received an unexpected chart')
        return self.native.evaluate('steep_power',Z,self.phase,**kwargs)


def steep_power_physical_identities():
    general=physical_collar_identities(); proofs={}
    a,q,z,qQ,S=s.symbols('a q Z qQ S',real=True)
    delta=2*a; k=1-a; p=1+delta; bh=s.Rational(1,2)+a; L=1-delta*z*z; d=1-z*z
    KQ=s.symbols('actual_KQ',positive=True); ell=qQ-q; K=KQ*s.exp(k*ell)
    A=s.Function('full_A')(q,z); E=s.Function('full_E')(q,z); P=s.Function('full_P')(q,z); b=(1-delta)/2
    Ct=(k*A-b*z*s.diff(A,z)-K)/L+2*S*s.exp(-q)*(s.diff(K,q)-(1+a)*K)
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*P+d*s.diff(P,z))/L
    rules={s.diff(A,q):K-k*A,s.diff(A,q,z):-k*s.diff(A,z),
           s.diff(E,q):delta*E-K*K,s.diff(E,q,z):delta*s.diff(E,z),
           s.diff(P,q):p*P-K*K/2,s.diff(P,q,z):p*s.diff(P,z)}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Power physical reduction failed: '+name)
        proofs[name]=True
    Dt=k*K/L+4*S*s.exp(-q)*K
    Dz=(d*s.diff(P,z)-2*z*(p*P-K*K/2))/L
    zero('actual_full_theta_homogeneous_mode_cancelled',(s.diff(Ct,q)+k*Ct).subs(rules,simultaneous=True)-Dt)
    zero('actual_full_axial_energy_mode_cancelled',(s.diff(Cz,q)-delta*Cz).subs(rules,simultaneous=True)-Dz)
    PdQ=s.Function('same_full_pressure_defect_Q')(z)
    closedP=K*K/6+s.exp(-p*ell)*(1/(2*p)+PdQ-KQ*KQ/6)
    NpQ=-z-2*p*z*PdQ+p*z*KQ*KQ/3+d*s.diff(PdQ,z)
    closedDz=2*k*z*K*K/(3*L)+s.exp(-p*ell)*NpQ/L
    zero('actual_full_pressure_datum_mode_and_power_part',Dz.subs({P:closedP,s.diff(P,z):s.diff(closedP,z)},simultaneous=True)-closedDz)
    factor_theta=s.exp(-bh*q); factor_axial=s.exp(-p*q)
    for j in range(3):
        expected=k*K/L*(-s.Rational(3,2))**j+4*S*s.exp(-q)*K*(-s.Rational(5,2))**j
        zero('actual_theta_divergence_full_factor_y'+str(j),s.diff(factor_theta*Dt,q,j)/factor_theta-expected)
        expected=[closedDz,-2*k*z*K*K/L,6*k*z*K*K/L][j]
        zero('actual_axial_pressure_homogeneous_full_factor_y'+str(j),s.diff(factor_axial*closedDz,q,j)/factor_axial-expected)
    F=1-2*(1-delta)*z*z-delta*delta*z**4
    reduced=2*F/L**3*s.diff(K,q)-4*z*z/L**2*s.diff(K,q,2)
    actual=sum(coefficient.subs({ZSYM:z,DSYM:delta})*s.diff(K,q,index[0])
               for index,coefficient in collar_velocity_operator_coefficients(0,2).items() if index[1]==0)
    zero('actual_nonzero_power_axial_viscosity_from_general_K_operator',-actual-reduced)
    # Recover the ACTUAL divergence producer AST on arbitrary endpoint data.
    fn=next(n for n in ast.parse(Path(__file__).read_text(encoding='utf8')).body
            if isinstance(n,ast.FunctionDef) and n.name=='steep_power_divergence_grids')
    fn.body=fn.body[:-1]+[ast.Return(value=ast.parse("dict(theta=theta,axial=axial)",mode='eval').body)]
    env=dict(math=math,IntervalTaylor=SimpleNamespace(variable=lambda ctx,value,order:value),
             axial_derivative=lambda value:s.diff(value,z))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual power divergence rows>','exec'),env)
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,k=k,delta=delta,prate=p,S=S)
    terminal=dict(KQ=KQ,pressure_defect_rows=[PdQ]); shape=dict(K_rows=[K],remaining=ell)
    rows=env['steep_power_divergence_grids'](heat,terminal,shape,q,z)
    for label,factor,source in (('theta',factor_theta,Dt),('axial',factor_axial,closedDz)):
        for j in range(3):
            for n in range(3-j):
                zero('actual_divergence_AST_'+label+'_y'+str(j)+'_Z'+str(n),
                     s.diff(rows[label][j]-s.diff(factor*source,q,j)/factor,z,n))
    zero('actual_power_source_q_at_exit_endpoint',(-s.Symbol('same_wait')-1-s.Symbol('same_Ts')*(1-s.Integer(1)))-(-s.Symbol('same_wait')-1+s.Integer(0)))
    return dict(general_K_physical_transfer=general,reduction_identities=proofs,
                original_K_is_independent_of_source_Z=True,
                actual_regional_physical_steep_power_stress_remainder_identity_verified=True,
                exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
                source_homogeneous_A_E_and_pressure_divergence_cancelled_before_enclosure=True,
                power_remainder_generally_nonzero=True,globally_flat_remainder_not_inferred=True)


def steep_power_adapter_binding(stress):
    tree=ast.parse(Path(__file__).read_text(encoding='utf8')); bindings={}
    requirements={
        ('PowerHeatReference','shape'):'self.actual_power_shape',
        ('PowerDispatch','provider'):"self.stress.steep if chart=='steep_power' else self.native.provider(chart)",
        ('PowerDispatch','evaluate'):'self.stress.steep_power(Z,coordinate)',
        ('AtPowerPhaseStress','collar'):"self.packet['steep_power_similarity_stress_mixed3_factored']",
        ('AtPowerPhaseAssembly','evaluate'):"self.native.evaluate('steep_power',Z,self.phase,**kwargs)",
        ('CompliantSteepPowerPhysicalC2','steep_power'):'-self.stress.steep.wait-1-self.stress.steep.Ts*(1-phase)'}
    for (cls,name),expression in requirements.items():
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==name)
        if not any(ast.dump(n)==ast.dump(ast.parse(expression,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Power adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    join=stress.bridge['actual_power_exit_formula_join']
    for flag in ('actual_power_exit_K_mixed4_join_verified','actual_power_exit_stress_mixed3_AST_join_verified',
                 'actual_power_exit_pressure_mixed4_AST_join_verified'):
        if not join[flag]:raise ValueError('Power physical requires actual source join: '+flag)
    # Every remainder mixed2 operator uses only K derivatives through four;
    # its positive physical factors and implicit-map coefficients are common.
    if any(sum(index)>4 for i in range(3) for j in range(3-i)
           for index in collar_velocity_operator_coefficients(i,j+2)):
        raise ValueError('Physical remainder operator exceeds the accepted source K join')
    path=HERE/(PREFIX+'steep_exit_physical_C2.py')
    exit_tree=ast.parse(path.read_text(encoding='utf8'))
    exit_cls=next(n for n in exit_tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantSteepExitPhysicalC2')
    exit_fn=next(n for n in exit_cls.body if isinstance(n,ast.FunctionDef) and n.name=='steep_out')
    if not any(ast.dump(n)==ast.dump(ast.parse('-self.stress.steep.wait-1+phase',mode='eval').body) for n in ast.walk(exit_fn)):
        raise ValueError('Actual steep-exit physical q source changed')
    # AST-bound common factor function is applied at the same endpoint q.
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    mu,a,wait,Ts,logone,logu,logp,logrp,Lrel=s.symbols('mu a wait Ts logone logu logP logRp Lrel',real=True)
    assembly=SimpleNamespace(ctx=c,logP=logp,logRp=logrp,
                             dispatch=SimpleNamespace(provider=lambda chart:SimpleNamespace(logEv2_parts=dict(inlet_log=2*logu))))
    heat=SimpleNamespace(a=a,mu=mu,bh=s.Rational(1,2)+a,
                         steep=SimpleNamespace(logone=logone,Ts=Ts,wait=wait,outer=SimpleNamespace(Lrel=Lrel)))
    parts_left=stress_source_log_parts(assembly,heat,-wait-1-Ts*(1-s.Integer(1)))
    parts_right=stress_source_log_parts(assembly,heat,-wait-1+s.Integer(0))
    if parts_left!=parts_right:raise ArithmeticError('Power/exit actual physical source factors differ at join')
    native_join=power_native_field_join_binding(stress)
    return dict(source_bindings=bindings,original_power_phase_not_uncorrelated_offset_division=True,
                same_source_similarity_stress_mixed3_pressure_mixed4_and_K4_join_consumed=True,
                same_actual_source_log_factors_and_original_q_at_right_join=True,
                mixed2_remainder_operator_uses_only_accepted_K4_source_join=True,
                physical_power_exit_stress_diagonal_divergence_remainder_join_verified=True,
                actual_native_velocity_pressure_and_remainder_source_join=native_join,
                actual_shape_pressure_stress_and_source_phase_consumed=True,
                general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar',
                physical_collar_source_never_evaluated_at_negative_offset=True)


def power_native_field_join_binding(stress):
    """Replay original field rows and the common physical remainder operator."""
    hashes={}; bindings={}; proofs={}
    def method(stem,name):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        return next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(n,ast.FunctionDef) and n.name==name)
    def expression(stem,name,target,wanted=None):
        fn=method(stem,name)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        if wanted is not None:values=[n for n in values if ast.dump(n)==ast.dump(ast.parse(wanted,mode='eval').body)]
        if len(values)!=1:raise ValueError('Actual power field join expression not unique: '+stem+'.'+name+'.'+target)
        bindings[stem+'.'+name+'.'+target]=True
        return values[0]
    def evaluate(node,env):return eval(compile(ast.Expression(node),'<actual power field join>','eval'),env,env)
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    mu,k,Ts,z=s.symbols('mu k Ts Z',real=True); thetaS=s.symbols('actual_thetaS',positive=True)
    source=SimpleNamespace(thetaS=thetaS,Ts=Ts,rate=1-mu,k=k)
    source.thetaQ=evaluate(expression('steep_waiting_C4','__init__','self.thetaQ'),dict(c=c,self=source))
    theta_left=evaluate(expression('steep_waiting_C4','steep_power','theta'),dict(c=c,self=source,one=s.Integer(1),t=Ts))
    theta_right=evaluate(expression('steep_waiting_C4','steep_out','theta'),dict(c=c,self=source,one=s.Integer(1),t=s.Integer(0),J=s.Integer(0)))
    rates_left=evaluate(expression('steep_waiting_C4','steep_power','rates'),dict(c=c,self=source))
    expression('steep_waiting_C4','steep_out','sig','sigma_jets(c,t)')
    if not stress.exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Actual original sigma flat endpoint jets required')
    proofs['actual_original_sigma_zero_jets_at_exit_phase0_consumed']=True
    rates_right=evaluate(expression('steep_waiting_C4','steep_out','rates'),dict(c=c,self=source,sig=[s.Integer(0)]*5,math=math))
    if s.simplify(theta_left-theta_right)!=0 or rates_left!=rates_right:raise ArithmeticError('Original power/exit theta/rate source join failed')
    # Bind both original packet wrappers to the same actual derivative helper.
    expression('steep_waiting_C4','packet','result',
               'self.outer._packet(Z,theta,X,energy,pressure,[one*r for r in rates],coordinate,dict(stage=stage,entire_steep_waiting_high_mixed_derivatives_available=True,same_actual_angular_terminal_histories_retained=True,exact_Gamma_and_both_epsilon_atoms_retained=True,exact_relative_velocity_log_parts=logs,no_forward_subtraction_of_unrelated_long_future_energy=True))')
    expression('power_angular_C4','_packet','mixed',
               'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    fn=method('flatten_mixed_C4','flatten_mixed'); fn.decorator_list=[]
    body=[]
    for node in fn.body:
        if isinstance(node,ast.If):continue
        body.append(node)
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='rows' for t in node.targets):
            body.append(ast.Return(value=ast.Name(id='rows',ctx=ast.Load()))); break
    fn.body=body
    env=dict(math=math,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value)),
             **{label:flatten_mixed.__globals__[label] for label in ('UR','UT','UZ','P')})
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual original mixed field rows>','exec'),env)
    common=[s.Function('same_native_'+label)(z) for label in ('X','energy','absolute_pressure')]
    left=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_left,theta_left,*common,s.Symbol('same_exact_Ev2'))
    right=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_right,theta_right,*common,s.Symbol('same_exact_Ev2'))
    for label in ('UR','UT','UZ'):
        key=env[label]
        for j in range(5):
            for n in range(5-j):
                if s.simplify(s.diff(left[key][j]-right[key][j],z,n))!=0:raise ArithmeticError('Original velocity mixed4 source join failed')
                proofs[label+'_y'+str(j)+'_Z'+str(n)]=True
    # Both pressure replacements are consumed explicitly; all mixed4 entries
    # come from the pressure rows already functionally joined by the source.
    for stem,name in (('steep_power_stress_C3','steep_power'),('steep_exit_stress_C3','steep_out')):
        expression(stem,name,'mixed[P]',"{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
        expression(stem,name,'fields[P]','[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
        fn=method(stem,name)
        keywords={n.arg:n.value for n in ast.walk(fn) if isinstance(n,ast.keyword) and n.arg in
                  ('physical_mixed_derivatives_total_order_le4','physical_velocity_and_pressure_y_derivative_Taylor')}
        for key,wanted in (('physical_mixed_derivatives_total_order_le4','mixed'),('physical_velocity_and_pressure_y_derivative_Taylor','fields')):
            if key not in keywords or ast.dump(keywords[key])!=ast.dump(ast.parse(wanted,mode='eval').body):
                raise ValueError('Original absolute-pressure packet routing changed')
        bindings[stem+'.actual_absolute_pressure_packet_routing']=True
    if not stress.bridge['actual_absolute_pressure_datum_and_exact_Ev0_units_consumed']:
        raise ValueError('Original pressure datum and amplitude source required')
    # Explicitly replay the ACTUAL completed mapper's remainder bracket.
    bracket_node=expression('collar_physical_C2','collar','ebracket','-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)')
    fn=method('collar_physical_C2','collar_velocity_bracket'); fn.decorator_list=[]
    env=dict(math=math,collar_velocity_operator_coefficients=collar_velocity_operator_coefficients,
             interval_expression=lambda ctx,value,z,delta,beta:value.subs({ZSYM:z,DSYM:delta}))
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual source remainder operator>','exec'),env)
    delta=s.symbols('same_delta',real=True)
    # The consumed actual K4 source join identifies these rows functionally;
    # no common interval box or cap defines the field.
    Krows=[[s.Symbol('same_K_y'+str(j))]+[s.Integer(0)]*5 for j in range(5)]
    provider=SimpleNamespace(heat=SimpleNamespace(delta=delta))
    for i in range(3):
        for j in range(3-i):
            e=evaluate(bracket_node,dict(c=c,K=Krows,i=i,j=j,z=z,self=provider,
                                       collar_velocity_bracket=env['collar_velocity_bracket']))
            expected=-sum(value.subs({ZSYM:z,DSYM:delta})*Krows[index[0]][index[1]]*math.factorial(index[1])
                          for index,value in collar_velocity_operator_coefficients(i,j+2).items())
            if s.simplify(e-expected)!=0:raise ArithmeticError('Actual remainder mixed2 source operator changed')
            proofs['actual_remainder_source_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(identities=proofs,source_bindings=bindings,
                actual_native_velocity_mixed4_AST_join_verified=True,
                actual_full_pressure_mixed4_packet_routing_join_verified=True,
                actual_remainder_mixed2_AST_operator_and_same_K4_join_verified=True,
                source_function_equality_not_interval_overlap=True,input_hashes=hashes)


class CompliantSteepPowerPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantSteepPowerStressC3(); self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=PowerDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx; self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Power physical family mismatch')
        self.proof=steep_power_physical_identities(); self.adapter_proof=steep_power_adapter_binding(self.stress)
        self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        self.hashes.update(self.adapter_proof['actual_native_velocity_pressure_and_remainder_source_join']['input_hashes'])
        for stem,gate in (('steep_power_stress_C3_check','actual_original_steep_power_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
                          ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
                          ('steep_exit_physical_C2_check','actual_regional_physical_steep_exit_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Power physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Power physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Power physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def steep_power(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); phase=c.mpf(phase); self.stress.steep._phase(phase)
        q=-self.stress.steep.wait-1-self.stress.steep.Ts*(1-phase); packet=self.stress.steep_power(z,phase)
        heat=PowerHeatReference(self.stress.heat,packet['steep_power_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtPowerPhaseStress(packet),assembly=AtPowerPhaseAssembly(self.assembly,phase))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        grids=steep_power_divergence_grids(heat,self.stress.terminal(z),packet['steep_power_shape'],q,z)
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
        point.update(chart='steep_power',steep_power_phase=phase,steep_power_offset_from_Rtail=q,
                     actual_regional_physical_steep_power_stress_remainder_identity_verified=True,
                     steep_power_homogeneous_A_E_and_pressure_divergence_cancelled_before_enclosure=True,
                     steep_power_exit_stress_mixed3_join_verified=True,steep_power_exit_pressure_mixed4_join_verified=True,
                     steep_power_exit_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                     steep_power_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
                     steep_power_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        samples=[summary(self.steep_power(z,t,lt,angle,nu)) for z,t,lt,angle,nu in
                 (('0','0','-1','0','1'),('.5','.5','-10','.7','.01'),('-.5','1','-100','1','.7'))]
        whole=summary(self.steep_power([-1,1],[0,1],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Original steep_power phase[0,1], Z[-1,1], tau>0, r>0; physical stress mixed3, diagonal/divergence/remainder mixed2',
                    samples=samples,whole_steep_power=whole,steep_power_physical_identities=self.proof,steep_power_adapter_binding=self.adapter_proof,
                    actual_regional_physical_steep_power_stress_remainder_identity_verified=True,
                    steep_power_exit_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                    steep_power_regional_remainder_exact_zero=False,steep_power_cone_certified=False,
                    independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepPowerPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual original steep-power completed physical stress and nonzero axial-viscosity remainder generated',flush=True)
    return result


if __name__=='__main__':run()
