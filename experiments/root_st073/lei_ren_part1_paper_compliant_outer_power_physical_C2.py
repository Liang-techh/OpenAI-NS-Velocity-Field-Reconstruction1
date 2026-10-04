"""Original preceding-power completed physical stress and right angular join.

SAME complete moments/absolute datum, original Lrel-4 and physical units.
The unchanged general-K operators retain full moment axial dependence even
when K_Z=0. Only this regional physical layer is admitted; its cone waits.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import CompliantOuterPowerStressC3,SourceAST
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_angular_physical_C2 import AngularHeatReference,angular_divergence_grids
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    CompliantCollarPhysicalC2,physical_source_row,scale_row,
    collar_velocity_operator_coefficients,stress_source_log_parts)
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import flatten_mixed
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def outer_power_physical_identities():
    """Constant-exponential specialization; full axial moment jets survive."""
    a,mu,q,z,S,KR=s.symbols('a mu q Z S KR',real=True)
    delta=2*a;k=1-a;p=1+delta;d=a-mu;b=(1-delta)/2
    L=1-delta*z*z;D=1-z*z;K=KR*s.exp(d*q);Q=K*K-1
    A=s.Function('same_full_A')(q,z);E=s.Function('same_full_E')(q,z);P=s.Function('same_full_P')(q,z)
    Ct=(k*A-b*z*s.diff(A,z)-K)/L+2*S*s.exp(-q)*(s.diff(K,q)-(1+a)*K)
    Cz=(delta*z*E-D*s.diff(E,z)/2-2*p*z*P+D*s.diff(P,z))/L
    rules={s.diff(A,q):K-k*A,s.diff(A,q,z):-k*s.diff(A,z),
        s.diff(E,q):delta*E-K*K,s.diff(E,q,z):delta*s.diff(E,z),
        s.diff(P,q):p*P-K*K/2,s.diff(P,q,z):p*s.diff(P,z)}
    Dt=-d*K/L+2*S*s.exp(-q)*mu*(1+mu)*K
    Dz=(D*s.diff(P,z)-2*p*z*P+z*K*K)/L
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Outer power physical identity failed: '+name)
        proofs[name]=True
    zero('constant_exponential_viscous_coefficient',d*d-p*d+a*(1+a)-mu*(1+mu))
    zero('full_axial_angular_history_cancels_in_divergence',(s.diff(Ct,q)+k*Ct).subs(rules,simultaneous=True)-Dt)
    zero('full_axial_energy_history_cancels_in_divergence',(s.diff(Cz,q)-delta*Cz).subs(rules,simultaneous=True)-Dz)
    Pd=P-1/(2*p)
    zero('same_pressure_and_square_unit_baselines_cancelled',(D*s.diff(Pd,z)-2*p*z*Pd+z*Q)/L-Dz)
    zero('actual_axial_full_factor_q1',(s.diff(Dz,q)-p*Dz).subs(rules,simultaneous=True)-2*d*z*K*K/L)
    zero('actual_axial_full_factor_q2',s.diff(2*d*z*K*K/L,q)-p*2*d*z*K*K/L-(4*d*d-2*d*p)*z*K*K/L)
    return dict(identities=proofs,actual_K_Z_exact_zero=True,full_A_E_P_axial_dependence_retained=True,
        same_general_K_divergence_and_axial_viscosity_operators_used=True,
        exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
        remainder_generally_nonzero_despite_K_Z_zero=True,globally_flat_remainder_not_inferred=True)


class OuterPowerDispatch:
    def __init__(self,native,stress):self.native=native;self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.outer if chart=='outer_power' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='outer_power':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.power(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
            physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
            original_outer_power_absolute_pressure_companion_used=True)


class AtOuterPowerPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['outer_power_similarity_stress_mixed3_factored'],
            Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtOuterPowerPhaseAssembly:
    def __init__(self,native,phase):self.native=native;self.phase=phase
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic outer-power map received an unexpected chart')
        return self.native.evaluate('outer_power',Z,self.phase,**kwargs)


def outer_power_native_field_join_binding(stress):
    asts=SourceAST();proofs={}
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    mu,z,length=s.symbols('mu Z original_Lrel',real=True);bp=s.Rational(1,2)+mu
    obj=SimpleNamespace(bp=bp,rate=1-mu,Lrel=length)
    power_fn=asts.method('power_angular_C4','power')
    packet_call=next(n for n in ast.walk(power_fn) if isinstance(n,ast.Call) and ast.unparse(n.func)=='self._packet')
    wanted='[one*0 for _ in range(4)]'
    if ast.dump(packet_call.args[5])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Native original power rate rows changed')
    rates_left=asts.evaluate(packet_call.args[5],dict(one=s.Integer(1)))
    theta_left=asts.evaluate(asts.expression('power_angular_C4','power','theta'),dict(c=c,self=obj,one=s.Integer(1),y=length-4))
    theta_right=asts.evaluate(asts.expression('power_angular_C4','angular','theta'),dict(c=c,self=obj,F=[s.Integer(1)],s=s.Integer(-4)))
    if s.simplify(theta_left-theta_right)!=0:raise ArithmeticError('Actual original power/angular theta differs')
    asts.expression('power_angular_C4','angular','rates',wanted='quotient_log_rates(F)')
    fn=asts.method('power_angular_C4','quotient_log_rates');fn.body=[n for n in fn.body if not isinstance(n,ast.If)]
    env=dict(math=math)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual angular endpoint rates>','exec'),env)
    rates_right=env['quotient_log_rates']([s.Integer(1)]+[s.Integer(0)]*4)
    if rates_left!=rates_right:raise ArithmeticError('Native power/angular rates differ')
    if not all(stress.angular.bridge['actual_original_endpoint_beta_jets_exact_zero'].values()):
        raise ValueError('Original angular endpoint beta jets required')
    asts.expression('power_angular_C4','_packet','mixed',
        wanted='flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    fn=asts.method('flatten_mixed_C4','flatten_mixed');fn.decorator_list=[];body=[]
    for node in fn.body:
        if isinstance(node,ast.If):continue
        body.append(node)
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='rows' for t in node.targets):
            body.append(ast.Return(value=ast.Name(id='rows',ctx=ast.Load())));break
    fn.body=body
    env=dict(math=math,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value)),
        **{label:flatten_mixed.__globals__[label] for label in ('UR','UT','UZ','P')})
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual native power mixed rows>','exec'),env)
    common=[s.Function('same_native_'+name)(z) for name in ('X','energy','absolute_pressure')]
    left=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_left,theta_left,*common,s.Symbol('same_exact_Ev2'))
    right=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_right,theta_right,*common,s.Symbol('same_exact_Ev2'))
    for label in ('UR','UT','UZ'):
        for j in range(5):
            for n in range(5-j):
                if s.simplify(s.diff(left[env[label]][j]-right[env[label]][j],z,n))!=0:
                    raise ArithmeticError('Actual native power/angular mixed4 velocity join failed')
                proofs[label+'_q'+str(j)+'_Z'+str(n)]=True
    for stem,name in (('outer_power_stress_C3','power'),('angular_stress_C3','angular')):
        asts.expression(stem,name,'mixed[P]',wanted="{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
        asts.expression(stem,name,'fields[P]',wanted='[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
        keywords={n.arg:n.value for n in ast.walk(asts.method(stem,name)) if isinstance(n,ast.keyword)}
        for key,wanted in (('physical_mixed_derivatives_total_order_le4','mixed'),('physical_velocity_and_pressure_y_derivative_Taylor','fields')):
            if ast.dump(keywords[key])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Actual complete pressure routes changed')
    join=stress.bridge['actual_outer_power_angular_formula_join']
    if not join['actual_outer_power_angular_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Actual same-source K/moment/stress/pressure join required')
    asts.expression('collar_physical_C2','collar','ebracket',wanted='-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)')
    for i in range(3):
        for j in range(3-i):
            if any(sum(index)>4 for index in collar_velocity_operator_coefficients(i,j+2)):
                raise ValueError('Power physical remainder exceeds same-source K4 join')
            proofs['same_remainder_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(identities=proofs,source_bindings=asts.bindings,input_hashes=asts.hashes,
        actual_native_velocity_mixed4_AST_join_verified=True,
        actual_full_pressure_mixed4_and_remainder_mixed2_source_join_verified=True,
        source_function_equality_not_interval_overlap=True)


def outer_power_adapter_binding(stress):
    asts=SourceAST();bindings={}
    required={
        ('OuterPowerDispatch','provider'):"self.stress.outer if chart=='outer_power' else self.native.provider(chart)",
        ('OuterPowerDispatch','evaluate'):'self.stress.power(Z,coordinate)',
        ('AtOuterPowerPhaseStress','collar'):"self.packet['outer_power_similarity_stress_mixed3_factored']",
        ('AtOuterPowerPhaseAssembly','evaluate'):"self.native.evaluate('outer_power',Z,self.phase,**kwargs)",
        ('CompliantOuterPowerPhysicalC2','power'):'-self.stress.angular.entry.steep.wait-self.stress.angular.entry.steep.Ts-6+(self.stress.outer.Lrel-4)*(phase-1)',
        ('CompliantOuterPowerPhysicalC2','power_shape'):"AngularHeatReference(self.stress.heat,packet['outer_power_shape'])",
        ('CompliantOuterPowerPhysicalC2','power_packet'):'self.stress.power(z,phase)',
        ('CompliantOuterPowerPhysicalC2','power_mapper'):'CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)',
        ('CompliantOuterPowerPhysicalC2','power_divergence'):"angular_divergence_grids(heat,packet['outer_power_future_defect_zeroth_Taylor'],packet['outer_power_shape'],q,z)"}
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    for (cls,name),expr in required.items():
        lookup='power' if name.startswith('power_') else name
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==lookup)
        if not any(ast.dump(n)==ast.dump(ast.parse(expr,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Outer power physical adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    # Bind the reused source-shape hook to the exact packet passed above.
    fn=asts.method('angular_physical_C2','shape')
    if len(fn.body)!=1 or not isinstance(fn.body[0],ast.Return) or ast.unparse(fn.body[0].value)!='self.actual_angular_shape':
        raise ValueError('Accepted general-K source-shape hook changed')
    # Physical divergence coefficients and units use the same actual source
    # expressions as the accepted angular mapper. No new transfer operator.
    own=asts.method('outer_power_physical_C2','power');other=asts.method('angular_physical_C2','angular')
    for target in ('coefficient',"divmixed[label]['r'+str(i)+'_z'+str(j)]",'divcart',
        "point['physical_momentum_residual_decomposition_cartesian']"):
        canonical=ast.unparse(ast.parse(target,mode='eval').body)
        lv=[n.value for n in ast.walk(own) if isinstance(n,ast.Assign) and any(ast.unparse(t)==canonical for t in n.targets)]
        rv=[n.value for n in ast.walk(other) if isinstance(n,ast.Assign) and any(ast.unparse(t)==canonical for t in n.targets)]
        if len(lv)!=1 or len(rv)!=1 or ast.dump(lv[0])!=ast.dump(rv[0]):raise ValueError('Accepted physical operator/source units changed: '+target)
        bindings['unchanged_angular_physical_transfer.'+target]=True
    mu,a,wait,Ts,length,origin,logP,logone,logu=s.symbols('mu a wait Ts Lrel logRp logP logone logu',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    obj=SimpleNamespace(ctx=c,logRp=origin,params=SimpleNamespace(mu=mu))
    env={}
    for node in ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):env[target.id]=ast.literal_eval(node.value)
    fn=asts.method('global_physical_assembly','radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual power/angular physical radius>','exec'),env)
    outer=SimpleNamespace(Lrel=length)
    if s.simplify(env['radius'](obj,'outer_power',s.Integer(1),{},outer)[0]-env['radius'](obj,'outer_angular',s.Integer(-4),{},outer)[0])!=0:
        raise ArithmeticError('Actual physical power/angular radius differs')
    left_q=asts.expression('outer_power_physical_C2','power','q')
    right_q=asts.expression('angular_physical_C2','angular','q')
    root=SimpleNamespace(stress=SimpleNamespace(outer=outer,angular=SimpleNamespace(entry=SimpleNamespace(steep=SimpleNamespace(wait=wait,Ts=Ts)))))
    qr=asts.evaluate(left_q,dict(self=root,phase=s.Integer(1)))
    ar=asts.evaluate(right_q,dict(self=SimpleNamespace(stress=root.stress.angular),offset=s.Integer(-4)))
    if s.simplify(qr-ar)!=0:raise ArithmeticError('Actual physical q sources differ at power/angular join')
    assembly=SimpleNamespace(ctx=c,logP=logP,logRp=origin,dispatch=SimpleNamespace(provider=lambda chart:SimpleNamespace(logEv2_parts={'inlet_log':2*logu})))
    heat=SimpleNamespace(a=a,mu=mu,bh=s.Rational(1,2)+a,steep=SimpleNamespace(logone=logone,Ts=Ts,wait=wait,outer=outer))
    if stress_source_log_parts(assembly,heat,qr)!=stress_source_log_parts(assembly,heat,ar):
        raise ArithmeticError('Actual physical power/angular source units differ')
    native=outer_power_native_field_join_binding(stress);native['input_hashes'].update(asts.hashes)
    return dict(source_bindings=bindings,actual_original_phase_to_q_factor_retained=True,
        actual_native_power_angular_radius_field_pressure_join_verified=True,
        same_source_log_factors_at_original_right_q=True,
        outer_power_angular_completed_physical_stress_diagonal_divergence_remainder_join_verified=True,
        actual_native_velocity_pressure_and_remainder_source_join=native,
        general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar and angular_divergence_grids',
        physical_collar_source_never_evaluated_at_negative_offset=True,left_flatten_physical_stress_join_verified=False)


class CompliantOuterPowerPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantOuterPowerStressC3();self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=OuterPowerDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx;self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Outer power physical family mismatch')
        self.proof=outer_power_physical_identities();self.adapter_proof=outer_power_adapter_binding(self.stress)
        self.hashes=dict(self.stress.hashes);self.hashes.update(self.assembly.hashes)
        self.hashes.update(self.adapter_proof['actual_native_velocity_pressure_and_remainder_source_join']['input_hashes'])
        for stem,gate in (
            ('outer_power_stress_C3_check','actual_original_outer_power_similarity_stress_recovered'),
            ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
            ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
            ('angular_physical_C2_check','actual_regional_physical_angular_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Outer power physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Outer power physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Outer power physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def power(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx;z=c.mpf(Z);phase=c.mpf(phase)
        q=-self.stress.angular.entry.steep.wait-self.stress.angular.entry.steep.Ts-6+(self.stress.outer.Lrel-4)*(phase-1)
        packet=self.stress.power(z,phase);heat=AngularHeatReference(self.stress.heat,packet['outer_power_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtOuterPowerPhaseStress(packet),assembly=AtOuterPowerPhaseAssembly(self.assembly,phase))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        grids=angular_divergence_grids(heat,packet['outer_power_future_defect_zeroth_Taylor'],packet['outer_power_shape'],q,z)
        logtau=c.mpf(log_tau);nu=c.mpf(viscosity);beta=-2-heat.delta
        parts=point['actual_source_log_factors'];logR=point['source_logR_enclosure'];divmixed={}
        for label,factor in (('theta','Qtheta'),('axial','Qz')):
            divmixed[label]={}
            for i in range(3):
                for j in range(3-i):
                    coefficient=physical_bracket(c,grids[label],i,j,z,heat.delta,beta-1)
                    divmixed[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,coefficient,parts[factor],
                        beta-1-i+j*(heat.delta-1),logtau,nu,i+j,c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        cs=point['cos_theta'];sn=point['sin_theta']
        divcart=dict(x=[scale_row(c,divmixed['theta']['r0_z0'],-sn)],y=[scale_row(c,divmixed['theta']['r0_z0'],cs)],z=[divmixed['axial']['r0_z0']])
        point['physical_cylindrical_stress_divergence_mixed2']=divmixed
        point['physical_completed_stress_divergence_cartesian']=divcart
        point['physical_momentum_residual_decomposition_cartesian']={label:[scale_row(c,row,-1) for row in divcart[label]]
            +point['physical_remainder_cartesian'][label] for label in ('x','y','z')}
        for key in ('actual_regional_physical_collar_stress_remainder_identity_verified','Gamma_endpoint_physical_remainder_exact_zero','collar_cone_certified'):
            point.pop(key,None)
        point.update(chart='outer_power',outer_power_phase=phase,outer_power_offset_from_Rtail=q,
            actual_regional_physical_outer_power_stress_remainder_identity_verified=True,
            full_axial_moment_history_and_pressure_baselines_cancelled_before_enclosure=True,
            actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,
            outer_power_angular_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            native_power_angular_radius_field_pressure_join_verified=True,
            original_phase_to_ordinary_q_factor_retained=True,left_flatten_physical_stress_join_verified=False,
            outer_power_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
            outer_power_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    @source_precision
    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        samples=[summary(self.power(z,v,lt,angle,nu)) for z,v,lt,angle,nu in
            (('0','0','-1','0','1'),('.5','.25','-10','.7','.01'),('-.5','.75','-100','1','.7'),('0','1','-1','0','1'))]
        whole=summary(self.power([-1,1],[0,1],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            scope='Original outer_power phase[0,1],Z[-1,1]; physical stress mixed3, diagonal/divergence/remainder mixed2',
            samples=samples,whole_outer_power=whole,outer_power_physical_identities=self.proof,outer_power_adapter_binding=self.adapter_proof,
            actual_regional_physical_outer_power_stress_remainder_identity_verified=True,
            outer_power_angular_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,
            original_phase_to_ordinary_q_factor_retained=True,left_flatten_physical_stress_join_verified=False,
            outer_power_regional_remainder_exact_zero=False,outer_power_cone_certified=False,
            independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


@source_precision
def run():
    result=CompliantOuterPowerPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original preceding-power completed physical stress, axial-viscosity remainder and angular right join generated',flush=True)
    return result


if __name__=='__main__':run()
