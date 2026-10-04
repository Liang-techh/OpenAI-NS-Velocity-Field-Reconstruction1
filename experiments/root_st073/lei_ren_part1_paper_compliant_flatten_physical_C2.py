"""Original100-unit flatten completed physical stress and right power join.

Same full moments, original variable K and source radius/units are retained.
General-K physical/divergence/remainder operators are unchanged.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_compliant_flatten_stress_C3 import CompliantFlattenStressC3,SourceAST
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_angular_physical_C2 import (
    AngularHeatReference,angular_divergence_grids,angular_physical_identities)
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


def flatten_physical_identities():
    general=angular_physical_identities();proofs={}
    a,mu,t,z,KR=s.symbols('a mu t Z Kright',real=True)
    delta=2*a;p=1+delta;d=a-mu;f=(1+z*z)/2;rho=s.log(f);sig=s.Function('original_sigma')(t)
    K=KR*s.exp(d*(t-100)+(sig-1)*rho);alpha=d+rho*s.diff(sig,t)
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Flatten physical source identity failed: '+name)
        proofs[name]=True
    zero('actual_variable_flatten_K_log_q_rate',s.diff(K,t)-K*alpha)
    zero('actual_variable_flatten_K_log_Z_rate',s.diff(K,z)-K*(sig-1)*2*z/(1+z*z))
    viscous=mu*(1+mu)+rho*(s.diff(sig,t,2)-(1+2*mu)*s.diff(sig,t))+rho*rho*s.diff(sig,t)**2
    zero('actual_variable_flatten_viscous_source_coefficient',s.diff(K,t,2)-p*s.diff(K,t)+a*(1+a)*K-K*viscous)
    zero('actual_flatten_right_K_constant_in_Z',K.subs(sig,1)-KR*s.exp(d*(t-100)))
    zero('actual_flatten_variable_shear_source',delta-2*alpha-(2*mu-2*rho*s.diff(sig,t)))
    return dict(general_K_source_operator=general,source_specialization_identities=proofs,
        actual_K_Z_and_K_ZZ_retained=True,full_A_E_P_axial_dependence_retained=True,
        exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
        generally_nonzero_axial_viscosity_remainder_retained=True,globally_flat_remainder_not_inferred=True)


class FlattenDispatch:
    def __init__(self,native,stress):self.native=native;self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.native if chart=='flatten' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='flatten':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.flatten(Z,coordinate)
        return dict(chart=chart,source_packet=packet,physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
            original_flatten_absolute_pressure_companion_used=True)


class AtFlattenTimeStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['flatten_similarity_stress_mixed3_factored'],
            Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtFlattenTimeAssembly:
    def __init__(self,native,t):self.native=native;self.t=t
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic flatten map received an unexpected chart')
        return self.native.evaluate('flatten',Z,self.t,**kwargs)


def flatten_native_field_join_binding(stress):
    asts=SourceAST();proofs={};c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    mu,z=s.symbols('mu Z',real=True);bp=s.Rational(1,2)+mu
    obj=SimpleNamespace(bp=bp,mu=mu,rate=1-mu)
    theta_left=asts.evaluate(asts.expression('flatten_mixed_C4','flatten','theta',
        wanted='IntervalTaylor.constant(c,c.exp(-100*self.bp)/2,5)'),
        dict(IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:value),c=c,self=obj))
    theta_right=asts.evaluate(asts.expression('power_angular_C4','power','theta'),
        dict(c=c,self=obj,one=s.Integer(1),y=s.Integer(0)))
    if s.simplify(theta_left-theta_right)!=0:raise ArithmeticError('Actual original flatten/power theta differs')
    # Original sigma is exactly1 with zero positive-order endpoint jets.
    endpoint=stress.power.angular.entry.power_source.exit_source.bridge
    if not endpoint['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Original flatten sigma endpoint source required')
    asts.expression('flatten_mixed_C4','flatten','sig',wanted='sigma_jets(c,t/100)')
    sj=asts.evaluate(asts.expression('flatten_mixed_C4','flatten','sj',
        wanted='[sig[k]*math.factorial(k)/100**k for k in range(5)]'),
        dict(sig=[s.Integer(1)]+[s.Integer(0)]*4,math=math))
    asts.expression('flatten_mixed_C4','flatten','mixed',
        wanted='flatten_mixed(c,self.mu,rho,sj,theta,X,energy,pressure,self.Ev2)')
    power_fn=asts.method('power_angular_C4','power')
    call=next(n for n in ast.walk(power_fn) if isinstance(n,ast.Call) and ast.unparse(n.func)=='self._packet')
    wanted='[one*0 for _ in range(4)]'
    if ast.dump(call.args[5])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Original power native rate rows changed')
    rates=asts.evaluate(call.args[5],dict(one=s.Integer(1)))
    asts.expression('power_angular_C4','_packet','mixed',
        wanted='flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    fn=asts.method('flatten_mixed_C4','flatten_mixed');fn.decorator_list=[];body=[]
    for node in fn.body:
        if isinstance(node,ast.If):continue
        body.append(node)
        if isinstance(node,ast.Assign) and any(ast.unparse(target)=='rows' for target in node.targets):
            body.append(ast.Return(value=ast.Name(id='rows',ctx=ast.Load())));break
    fn.body=body
    env=dict(math=math,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value)),
        **{label:flatten_mixed.__globals__[label] for label in ('UR','UT','UZ','P')})
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual native flatten mixed rows>','exec'),env)
    common=[s.Function('same_native_'+name)(z) for name in ('X','energy','absolute_pressure')]
    left=env['flatten_mixed'](c,mu,s.log((1+z*z)/2),sj,theta_left,*common,s.Symbol('same_exact_Ev2'))
    right=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates,theta_right,*common,s.Symbol('same_exact_Ev2'))
    for label in ('UR','UT','UZ'):
        for j in range(5):
            for n in range(5-j):
                if s.simplify(s.diff(left[env[label]][j]-right[env[label]][j],z,n))!=0:
                    raise ArithmeticError('Actual native flatten/power mixed4 velocity join failed')
                proofs[label+'_q'+str(j)+'_Z'+str(n)]=True
    for stem,name in (('flatten_stress_C3','flatten'),('outer_power_stress_C3','power')):
        asts.expression(stem,name,'mixed[P]',wanted="{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
        asts.expression(stem,name,'fields[P]',wanted='[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
        keywords={n.arg:n.value for n in ast.walk(asts.method(stem,name)) if isinstance(n,ast.keyword)}
        for key,wanted in (('physical_mixed_derivatives_total_order_le4','mixed'),('physical_velocity_and_pressure_y_derivative_Taylor','fields')):
            if ast.dump(keywords[key])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Actual full pressure routes changed')
    if not stress.bridge['actual_flatten_power_formula_join']['actual_flatten_power_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Same full source flatten/power moment join missing')
    asts.expression('collar_physical_C2','collar','ebracket',wanted='-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)')
    for i in range(3):
        for j in range(3-i):
            if any(sum(index)>4 for index in collar_velocity_operator_coefficients(i,j+2)):
                raise ValueError('Flatten physical remainder exceeds same-source K4')
            proofs['same_remainder_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(identities=proofs,source_bindings=asts.bindings,input_hashes=asts.hashes,
        actual_native_velocity_mixed4_AST_join_verified=True,
        actual_full_pressure_mixed4_and_remainder_mixed2_source_join_verified=True,
        source_function_equality_not_interval_overlap=True)


def flatten_adapter_binding(stress):
    asts=SourceAST();bindings={}
    required={
        ('FlattenDispatch','provider'):"self.stress.native if chart=='flatten' else self.native.provider(chart)",
        ('FlattenDispatch','evaluate'):'self.stress.flatten(Z,coordinate)',
        ('AtFlattenTimeStress','collar'):"self.packet['flatten_similarity_stress_mixed3_factored']",
        ('AtFlattenTimeAssembly','evaluate'):"self.native.evaluate('flatten',Z,self.t,**kwargs)",
        ('CompliantFlattenPhysicalC2','flatten'):'-self.stress.heat.steep.wait-self.stress.heat.steep.Ts-2-self.stress.power.outer.Lrel+(t-100)',
        ('CompliantFlattenPhysicalC2','flatten_shape'):"AngularHeatReference(self.stress.heat,packet['flatten_shape'])",
        ('CompliantFlattenPhysicalC2','flatten_packet'):'self.stress.flatten(z,t)',
        ('CompliantFlattenPhysicalC2','flatten_mapper'):'CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)',
        ('CompliantFlattenPhysicalC2','flatten_divergence'):"angular_divergence_grids(heat,packet['flatten_future_defect_zeroth_Taylor'],packet['flatten_shape'],q,z)"}
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    for (cls,name),expr in required.items():
        lookup='flatten' if name.startswith('flatten_') else name
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==lookup)
        if not any(ast.dump(n)==ast.dump(ast.parse(expr,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Flatten physical adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    fn=asts.method('angular_physical_C2','shape')
    if len(fn.body)!=1 or not isinstance(fn.body[0],ast.Return) or ast.unparse(fn.body[0].value)!='self.actual_angular_shape':
        raise ValueError('Accepted general-K source-shape hook changed')
    own=asts.method('flatten_physical_C2','flatten');other=asts.method('angular_physical_C2','angular')
    for target in ('coefficient',"divmixed[label]['r'+str(i)+'_z'+str(j)]",'divcart',
        "point['physical_momentum_residual_decomposition_cartesian']"):
        canonical=ast.unparse(ast.parse(target,mode='eval').body)
        lv=[n.value for n in ast.walk(own) if isinstance(n,ast.Assign) and any(ast.unparse(t)==canonical for t in n.targets)]
        rv=[n.value for n in ast.walk(other) if isinstance(n,ast.Assign) and any(ast.unparse(t)==canonical for t in n.targets)]
        if len(lv)!=1 or len(rv)!=1 or ast.dump(lv[0])!=ast.dump(rv[0]):
            raise ValueError('Accepted general-K physical units/operator changed: '+target)
        bindings['unchanged_general_K_transfer.'+target]=True
    mu,a,wait,Ts,length,origin,logP,logone,logu=s.symbols('mu a wait Ts Lrel logRp logP logone logu',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    obj=SimpleNamespace(ctx=c,logRp=origin,params=SimpleNamespace(mu=mu));env={}
    for node in ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):env[target.id]=ast.literal_eval(node.value)
    fn=asts.method('global_physical_assembly','radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual flatten/power physical radius>','exec'),env)
    outer=SimpleNamespace(Lrel=length)
    if s.simplify(env['radius'](obj,'flatten',100,{},None)[0]-env['radius'](obj,'outer_power',0,{},outer)[0])!=0:
        raise ArithmeticError('Actual flatten/power physical radius differs')
    leftq=asts.expression('flatten_physical_C2','flatten','q');rightq=asts.expression('outer_power_physical_C2','power','q')
    steep=SimpleNamespace(wait=wait,Ts=Ts)
    root=SimpleNamespace(stress=SimpleNamespace(heat=SimpleNamespace(steep=steep),power=SimpleNamespace(outer=outer)))
    qr=asts.evaluate(leftq,dict(self=root,t=s.Integer(100)))
    ar=asts.evaluate(rightq,dict(self=SimpleNamespace(stress=SimpleNamespace(angular=SimpleNamespace(entry=SimpleNamespace(steep=steep)),outer=outer)),phase=s.Integer(0)))
    if s.simplify(qr-ar)!=0:raise ArithmeticError('Actual flatten/power physical q differs')
    assembly=SimpleNamespace(ctx=c,logP=logP,logRp=origin,dispatch=SimpleNamespace(provider=lambda chart:SimpleNamespace(logEv2_parts={'inlet_log':2*logu})))
    heat=SimpleNamespace(a=a,mu=mu,bh=s.Rational(1,2)+a,steep=SimpleNamespace(logone=logone,Ts=Ts,wait=wait,outer=outer))
    if stress_source_log_parts(assembly,heat,qr)!=stress_source_log_parts(assembly,heat,ar):
        raise ArithmeticError('Actual physical flatten/power source logs differ')
    native=flatten_native_field_join_binding(stress);native['input_hashes'].update(asts.hashes)
    return dict(source_bindings=bindings,ordinary_t_equals_q_derivatives=True,
        actual_native_flatten_power_radius_field_pressure_join_verified=True,same_source_log_factors_at_original_right_q=True,
        flatten_power_completed_physical_stress_diagonal_divergence_remainder_join_verified=True,
        actual_native_velocity_pressure_and_remainder_source_join=native,
        general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar and angular_divergence_grids',
        physical_collar_source_never_evaluated_at_negative_offset=True,left_pulse_physical_stress_join_verified=False)


class CompliantFlattenPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantFlattenStressC3();self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=FlattenDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx;self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Flatten physical family mismatch')
        self.proof=flatten_physical_identities();self.adapter_proof=flatten_adapter_binding(self.stress)
        self.hashes=dict(self.stress.hashes);self.hashes.update(self.assembly.hashes)
        self.hashes.update(self.adapter_proof['actual_native_velocity_pressure_and_remainder_source_join']['input_hashes'])
        for stem,gate in (
            ('flatten_stress_C3_check','actual_original_flatten_similarity_stress_recovered'),
            ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
            ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
            ('angular_physical_C2_check','actual_regional_physical_angular_stress_remainder_identity_verified'),
            ('outer_power_physical_C2_check','actual_regional_physical_outer_power_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Flatten physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Flatten physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def flatten(self,Z,t,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx;z=c.mpf(Z);t=c.mpf(t)
        q=-self.stress.heat.steep.wait-self.stress.heat.steep.Ts-2-self.stress.power.outer.Lrel+(t-100)
        packet=self.stress.flatten(z,t);heat=AngularHeatReference(self.stress.heat,packet['flatten_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtFlattenTimeStress(packet),assembly=AtFlattenTimeAssembly(self.assembly,t))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        grids=angular_divergence_grids(heat,packet['flatten_future_defect_zeroth_Taylor'],packet['flatten_shape'],q,z)
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
        point.update(chart='flatten',flatten_time=t,flatten_offset_from_Rtail=q,
            actual_regional_physical_flatten_stress_remainder_identity_verified=True,
            full_axial_moment_history_and_pressure_baselines_cancelled_before_enclosure=True,
            actual_K_Z_and_K_ZZ_and_full_moment_axial_dependence_retained=True,
            flatten_power_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            native_flatten_power_radius_field_pressure_join_verified=True,
            ordinary_t_equals_q_derivatives=True,left_pulse_physical_stress_join_verified=False,
            flatten_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
            flatten_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    @source_precision
    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        samples=[summary(self.flatten(z,t,lt,angle,nu)) for z,t,lt,angle,nu in
            (('0','0','-1','0','1'),('.5','50','-10','.7','.01'),('0','100','-100','1','.7'))]
        whole=summary(self.flatten([-1,1],[0,100],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            scope='Original100-unit flatten t[0,100],Z[-1,1]; physical stress mixed3, diagonal/divergence/remainder mixed2',
            samples=samples,whole_flatten=whole,flatten_physical_identities=self.proof,flatten_adapter_binding=self.adapter_proof,
            actual_regional_physical_flatten_stress_remainder_identity_verified=True,
            flatten_power_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            actual_K_Z_and_K_ZZ_and_full_moment_axial_dependence_retained=True,ordinary_t_equals_q_derivatives=True,
            left_pulse_physical_stress_join_verified=False,flatten_regional_remainder_exact_zero=False,flatten_cone_certified=False,
            independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


@source_precision
def run():
    result=CompliantFlattenPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original100-unit flatten completed physical stress, retained axial-viscosity remainder and power right join generated',flush=True)
    return result


if __name__=='__main__':run()
