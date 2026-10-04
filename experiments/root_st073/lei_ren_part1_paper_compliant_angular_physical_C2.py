"""Original angular repair: general-K physical stress and axial viscosity.

Selected axial coefficient functions, original beta pulses, full moments and
absolute pressure are retained. Only regional physical transfer/right joins
are admitted here; the angular cone and preceding-power physical join wait.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_compliant_angular_stress_C3 import CompliantAngularStressC3,SourceAST
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    CompliantCollarPhysicalC2,physical_collar_identities,physical_source_row,scale_row,
    collar_velocity_operator_coefficients)
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import flatten_mixed
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def angular_divergence_grids(heat,defects,shape,q,Z):
    """Full-factor q rows of Div(T), mixed2; every K_Z/Q_Z term remains."""
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5)
    L=1-z*z*heat.delta; D=1-z*z; K=shape['K_rows']; p=heat.prate
    b=(1-heat.delta)/2
    viscous=[K[j+2]-K[j+1]*p+K[j]*(heat.a*(1+heat.a)) for j in range(3)]
    shear=shifted_rows(viscous,-1,2)
    theta_raw=[-(K[j+1]+z*axial_derivative(K[j])*b)/L
               +shear[j]*(2*heat.S*c.exp(-q)) for j in range(3)]
    theta=shifted_rows(theta_raw,-heat.a-c.mpf('.5'),2)
    Pd=defects['pressure_defect_rows']; Q=shape['K_squared_defect_rows']
    axial=[(D*axial_derivative(Pd)-z*Pd*(2*p)+z*Q[0])/L,
           (z*Q[1]-D*axial_derivative(Q[0])/2)/L,
           (z*(Q[2]-Q[1]*p)-D*(axial_derivative(Q[1])-axial_derivative(Q[0])*p)/2)/L]
    return {label:{'y'+str(j)+'_Z'+str(n):rows[j][n]*math.factorial(n)
                   for j in range(3) for n in range(3-j)} for label,rows in (('theta',theta),('axial',axial))}


def angular_physical_identities():
    general=physical_collar_identities(); proofs={}
    a,q,z,S=s.symbols('a q Z S',real=True); delta=2*a; k=1-a; p=1+delta
    b=(1-delta)/2; L=1-delta*z*z; D=1-z*z
    K=s.Function('same_original_angular_K')(q,z)
    A=s.Function('same_full_A')(q,z); E=s.Function('same_full_E')(q,z); P=s.Function('same_full_P')(q,z)
    Ct=(k*A-b*z*s.diff(A,z)-K)/L+2*S*s.exp(-q)*(s.diff(K,q)-(1+a)*K)
    Cz=(delta*z*E-D*s.diff(E,z)/2-2*p*z*P+D*s.diff(P,z))/L
    rules={s.diff(A,q):K-k*A,s.diff(A,q,z):s.diff(K,z)-k*s.diff(A,z),
           s.diff(E,q):delta*E-K*K,s.diff(E,q,z):delta*s.diff(E,z)-s.diff(K*K,z),
           s.diff(P,q):p*P-K*K/2,s.diff(P,q,z):p*s.diff(P,z)-s.diff(K*K,z)/2}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Angular physical identity failed: '+name)
        proofs[name]=True
    Dt=-(s.diff(K,q)+b*z*s.diff(K,z))/L+2*S*s.exp(-q)*(s.diff(K,q,2)-p*s.diff(K,q)+a*(1+a)*K)
    Dz=(D*s.diff(P,z)-2*p*z*P+z*K*K)/L
    zero('actual_axial_dependent_angular_homogeneous_mode_cancelled',(s.diff(Ct,q)+k*Ct).subs(rules,simultaneous=True)-Dt)
    zero('actual_axial_dependent_energy_homogeneous_mode_cancelled',(s.diff(Cz,q)-delta*Cz).subs(rules,simultaneous=True)-Dz)
    Pd=P-1/(2*p); Q=K*K-1
    zero('actual_pressure_and_squared_K_unit_baselines_cancelled',(D*s.diff(Pd,z)-2*p*z*Pd+z*Q)/L-Dz)
    expected_axial=[Dz,(z*s.diff(Q,q)-D*s.diff(Q,z)/2)/L,
        (z*(s.diff(Q,q,2)-p*s.diff(Q,q))-D*(s.diff(Q,q,z)-p*s.diff(Q,z))/2)/L]
    axial=[Dz]
    for j in range(2):axial.append(s.diff(axial[-1],q).subs(rules,simultaneous=True)-p*axial[-1])
    for j in range(3):zero('actual_axial_full_factor_q'+str(j),axial[j]-expected_axial[j])
    fn=next(n for n in ast.parse(Path(__file__).read_text(encoding='utf8')).body
            if isinstance(n,ast.FunctionDef) and n.name=='angular_divergence_grids')
    fn.body=fn.body[:-1]+[ast.Return(value=ast.parse('dict(theta=theta,axial=axial)',mode='eval').body)]
    env=dict(math=math,IntervalTaylor=SimpleNamespace(variable=lambda ctx,value,order:value),
             axial_derivative=lambda value:s.diff(value,z),shifted_rows=shifted_rows)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual angular divergence>','exec'),env)
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,a=a,S=S,delta=delta,prate=p)
    shape=dict(K_rows=[s.diff(K,q,j) for j in range(5)],K_squared_defect_rows=[s.diff(Q,q,j) for j in range(5)])
    rows=env['angular_divergence_grids'](heat,dict(pressure_defect_rows=Pd),shape,q,z)
    for label,factor,target in (('theta',s.exp(-(a+s.Rational(1,2))*q),Dt),('axial',s.exp(-p*q),Dz)):
        for j in range(3):
            expected=s.diff(factor*target,q,j)/factor if label=='theta' else expected_axial[j]
            for n in range(3-j):zero('actual_divergence_AST_'+label+'_q'+str(j)+'_Z'+str(n),s.diff(rows[label][j]-expected,z,n))
    for i in range(3):
        for j in range(3-i):
            if any(sum(index)>4 for index in collar_velocity_operator_coefficients(i,j+2)):
                raise ValueError('Angular axial viscosity exceeds admitted K mixed4')
            proofs['same_general_axial_viscosity_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(general_K_physical_transfer=general,reduction_identities=proofs,
        actual_K_Z_and_K_ZZ_retained=True,entry_Z_independent_reduction_used=False,
        actual_regional_physical_angular_stress_remainder_identity_verified=True,
        exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
        source_homogeneous_A_E_and_pressure_baselines_cancelled_before_enclosure=True,
        angular_remainder_generally_nonzero=True,globally_flat_remainder_not_inferred=True)


class AngularHeatReference:
    def __init__(self,heat,shape):self.actual_heat=heat; self.actual_angular_shape=shape
    def __getattr__(self,key):return getattr(self.actual_heat,key)
    def shape(self,Z,q):return self.actual_angular_shape


class AngularDispatch:
    def __init__(self,native,stress):self.native=native; self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.outer if chart=='outer_angular' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='outer_angular':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.angular(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
            physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
            original_angular_absolute_pressure_companion_used=True)


class AtAngularPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['angular_similarity_stress_mixed3_factored'],
                    Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtAngularPhaseAssembly:
    def __init__(self,native,offset):self.native=native; self.offset=offset
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic angular map received an unexpected chart')
        return self.native.evaluate('outer_angular',Z,self.offset,**kwargs)


def angular_native_field_join_binding(stress):
    asts=SourceAST(); proofs={}; c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    mu,z,length=s.symbols('mu Z Lrel',real=True); bp=s.Rational(1,2)+mu
    obj=SimpleNamespace(bp=bp,rate=1-mu,Lrel=length,outer=SimpleNamespace(Lrel=length))
    obj.thetaR=asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.thetaR'),dict(c=c,self=obj))
    theta_left=asts.evaluate(asts.expression('power_angular_C4','angular','theta'),dict(c=c,self=obj,F=[s.Integer(1)],s=s.Integer(0)))
    theta_right=asts.evaluate(asts.expression('steep_waiting_C4','steep_in','theta'),dict(c=c,self=obj,one=s.Integer(1),t=s.Integer(0),J=s.Integer(0)))
    if s.simplify(theta_left-theta_right)!=0:raise ArithmeticError('Actual angular-entry theta normalization differs')
    asts.expression('power_angular_C4','angular','rates',wanted='quotient_log_rates(F)')
    asts.expression('steep_waiting_C4','steep_in','rates',wanted='[-self.rate*sig[j]*math.factorial(j) for j in range(4)]')
    if not all(stress.bridge['actual_original_endpoint_beta_jets_exact_zero'].values()):
        raise ValueError('Original angular endpoint beta jets required')
    if not stress.entry.power_source.exit_source.bridge['original_sigma_between0and1_reflection_and_flat_endpoints_bound']:
        raise ValueError('Original entry sigma endpoint jets required')
    # Replay both actual rate implementations at their original flat endpoints.
    fn=asts.method('power_angular_C4','quotient_log_rates')
    fn.body=[node for node in fn.body if not isinstance(node,ast.If)]
    env=dict(math=math)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual angular endpoint rates>','exec'),env)
    rates_left=env['quotient_log_rates']([s.Integer(1)]+[s.Integer(0)]*4)
    rates_right=asts.evaluate(asts.expression('steep_waiting_C4','steep_in','rates'),dict(self=obj,sig=[s.Integer(0)]*4,math=math))
    if rates_left!=rates_right:raise ArithmeticError('Actual angular-entry endpoint rate functions differ')
    asts.expression('power_angular_C4','_packet','mixed',
        wanted='flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    fn=asts.method('flatten_mixed_C4','flatten_mixed'); fn.decorator_list=[]; body=[]
    for node in fn.body:
        if isinstance(node,ast.If):continue
        body.append(node)
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='rows' for t in node.targets):
            body.append(ast.Return(value=ast.Name(id='rows',ctx=ast.Load()))); break
    fn.body=body
    env=dict(math=math,IntervalTaylor=SimpleNamespace(constant=lambda ctx,value,order:s.Rational(value)),
        **{label:flatten_mixed.__globals__[label] for label in ('UR','UT','UZ','P')})
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<actual angular native mixed rows>','exec'),env)
    common=[s.Function('same_native_'+name)(z) for name in ('X','energy','absolute_pressure')]
    left=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_left,theta_left,*common,s.Symbol('same_exact_Ev2'))
    right=env['flatten_mixed'](c,mu,s.Integer(1),[s.Integer(0)]+rates_right,theta_right,*common,s.Symbol('same_exact_Ev2'))
    for label in ('UR','UT','UZ'):
        for j in range(5):
            for n in range(5-j):
                if s.simplify(s.diff(left[env[label]][j]-right[env[label]][j],z,n))!=0:raise ArithmeticError('Native angular velocity mixed4 join failed')
                proofs[label+'_q'+str(j)+'_Z'+str(n)]=True
    for stem,name in (('angular_stress_C3','angular'),('steep_entry_stress_C3','steep_in')):
        asts.expression(stem,name,'mixed[P]',wanted="{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
        asts.expression(stem,name,'fields[P]',wanted='[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
        keywords={n.arg:n.value for n in ast.walk(asts.method(stem,name)) if isinstance(n,ast.keyword)}
        for key,wanted in (('physical_mixed_derivatives_total_order_le4','mixed'),('physical_velocity_and_pressure_y_derivative_Taylor','fields')):
            if ast.dump(keywords[key])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Actual angular pressure packet route changed')
    if not stress.bridge['actual_angular_entry_formula_join']['actual_angular_entry_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified']:
        raise ValueError('Same full angular-entry stress/pressure join required')
    asts.expression('collar_physical_C2','collar','ebracket',wanted='-collar_velocity_bracket(c,K,i,j+2,z,self.heat.delta)')
    for i in range(3):
        for j in range(3-i):proofs['same_remainder_operator_r'+str(i)+'_z'+str(j)]=True
    return dict(identities=proofs,source_bindings=asts.bindings,input_hashes=asts.hashes,
        actual_native_velocity_mixed4_AST_join_verified=True,
        actual_full_pressure_mixed4_and_remainder_mixed2_source_join_verified=True,
        source_function_equality_not_interval_overlap=True)


def angular_adapter_binding(stress):
    tree=ast.parse(Path(__file__).read_text(encoding='utf8')); bindings={}
    required={('AngularHeatReference','shape'):'self.actual_angular_shape',
        ('AngularDispatch','provider'):"self.stress.outer if chart=='outer_angular' else self.native.provider(chart)",
        ('AngularDispatch','evaluate'):'self.stress.angular(Z,coordinate)',
        ('AtAngularPhaseStress','collar'):"self.packet['angular_similarity_stress_mixed3_factored']",
        ('AtAngularPhaseAssembly','evaluate'):"self.native.evaluate('outer_angular',Z,self.offset,**kwargs)",
        ('CompliantAngularPhysicalC2','angular'):'-self.stress.entry.steep.wait-self.stress.entry.steep.Ts-2+offset'}
    for (cls,name),expr in required.items():
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==name)
        if not any(ast.dump(n)==ast.dump(ast.parse(expr,mode='eval').body) for n in ast.walk(fn)):
            raise ValueError('Angular physical adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    asts=SourceAST(); fn=asts.method('global_physical_assembly','radius'); env={}
    for node in ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):env[target.id]=ast.literal_eval(node.value)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual angular physical radius>','exec'),env)
    mu,Ts,wait,length,origin=s.symbols('mu Ts wait Lrel logRp',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    assembly=SimpleNamespace(ctx=c,logRp=origin,params=SimpleNamespace(mu=mu))
    outer=SimpleNamespace(Lrel=length); native=SimpleNamespace(outer=outer,Ts=Ts,wait=wait)
    left=env['radius'](assembly,'outer_angular',s.Integer(0),{},outer)[0]
    right=env['radius'](assembly,'steep_entry',s.Integer(0),{},native)[0]
    if s.simplify(left-right)!=0:raise ArithmeticError('Actual angular-entry radius differs')
    nativejoin=angular_native_field_join_binding(stress); nativejoin['input_hashes'].update(asts.hashes)
    return dict(source_bindings=bindings,ordinary_q_equals_original_angular_s=True,
        actual_native_angular_entry_radius_field_pressure_join_verified=True,
        actual_native_velocity_pressure_and_remainder_source_join=nativejoin,
        angular_entry_physical_stress_diagonal_divergence_remainder_join_verified=True,
        general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar',
        physical_collar_source_never_evaluated_at_negative_offset=True,
        preceding_power_physical_stress_join_verified=False)


class CompliantAngularPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantAngularStressC3(); self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=AngularDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx; self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Angular physical family mismatch')
        self.proof=angular_physical_identities(); self.adapter_proof=angular_adapter_binding(self.stress)
        self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        self.hashes.update(self.adapter_proof['actual_native_velocity_pressure_and_remainder_source_join']['input_hashes'])
        for stem,gate in (('angular_stress_C3_check','actual_original_angular_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
                          ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
                          ('steep_entry_physical_C2_check','actual_regional_physical_steep_entry_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Angular physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Angular physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Angular physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def angular(self,Z,offset,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); offset=c.mpf(offset)
        q=-self.stress.entry.steep.wait-self.stress.entry.steep.Ts-2+offset
        packet=self.stress.angular(z,offset); heat=AngularHeatReference(self.stress.heat,packet['angular_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtAngularPhaseStress(packet),assembly=AtAngularPhaseAssembly(self.assembly,offset))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        grids=angular_divergence_grids(heat,packet['angular_future_defect_zeroth_Taylor'],packet['angular_shape'],q,z)
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
        point.update(chart='outer_angular',angular_offset=offset,angular_offset_from_Rtail=q,
            actual_regional_physical_angular_stress_remainder_identity_verified=True,
            angular_homogeneous_A_E_and_pressure_baselines_cancelled_before_enclosure=True,
            actual_K_Z_and_K_ZZ_retained=True,entry_Z_independent_reduction_used=False,
            angular_entry_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            native_angular_entry_radius_field_pressure_join_verified=True,
            preceding_power_physical_stress_join_verified=False,
            angular_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
            angular_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    @source_precision
    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            return result
        samples=[summary(self.angular(z,v,lt,angle,nu)) for z,v,lt,angle,nu in
                 (('0','-4','-1','0','1'),('.5','-3','-10','.7','.01'),('-.5','-1','-100','1','.7'),('0','0','-1','0','1'))]
        whole=summary(self.angular([-1,1],[-4,0],[-1000,-1],None))
        crossings=[summary(self.angular([-1,1],[center+side*.15-.01,center+side*.15+.01],[-1000,-1],None))
                   for center in (-3,-1) for side in (-1,1)]
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            scope='Original angular s[-4,0], Z[-1,1]; physical stress mixed3, diagonal/divergence/remainder mixed2',
            samples=samples,whole_angular=whole,support_crossings=crossings,
            angular_physical_identities=self.proof,angular_adapter_binding=self.adapter_proof,
            actual_regional_physical_angular_stress_remainder_identity_verified=True,
            angular_entry_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
            actual_K_Z_and_K_ZZ_retained=True,entry_Z_independent_reduction_used=False,
            preceding_power_physical_stress_join_verified=False,angular_regional_remainder_exact_zero=False,
            angular_cone_certified=False,independently_bounded_global_flat_remainder=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,
            full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


@source_precision
def run():
    result=CompliantAngularPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original angular general-K completed physical stress and nonzero axial-viscosity remainder generated',flush=True)
    return result


if __name__=='__main__':run()
