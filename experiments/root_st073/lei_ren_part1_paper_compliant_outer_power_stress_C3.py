"""Original preangular power: SAME full moments and analytic pressure.

The whole original Lrel-4 is retained. Native forward angular history is
kept; quadratic full future moments are transported from actual angular
s=-4 with stable factored exponentials. Physical/cone companions wait.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_angular_stress_C3 import CompliantAngularStressC3,SourceAST,quotient_rows
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_stress_rows,axial_derivative
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_heat_pressure_C4 import pressure_y_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(phase=[0,1],Z=[-1,1],offset='q=-wait-Ts-2-Lrel+(Lrel-4)*phase',
            ordinary_derivative='d_q=(Lrel-4)^(-1)*d_phase')


def outer_power_shape(heat,terminal,offset):
    c=heat.ctx; rate=heat.a-heat.mu; one=IntervalTaylor.constant(c,1,5)
    factor=terminal['Kright']*c.exp(rate*offset)
    K=[one*(factor*rate**j) for j in range(5)]
    square=product_rows(K,K)
    return dict(K_rows=K,K_defect_rows=[K[0]-one]+K[1:],
        K_squared_defect_rows=[square[0]-one]+square[1:],
        log_K_rate_rows=[one*rate]+[one*0 for _ in range(3)],
        actual_K_is_independent_of_source_Z=True)


def outer_power_X_rows(heat,X0):
    one=IntervalTaylor.constant(heat.ctx,1,5); r=1-heat.mu; X=[X0]
    for j in range(1,5):X.append((one if j==1 else one*0)-X[j-1]*r)
    return X


def outer_power_defect_rows(heat,terminal,offset,shape,X0):
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5); d=heat.a-heat.mu
    e=c.exp(heat.delta*offset); ep=c.exp(heat.prate*offset); ek=c.exp(2*d*offset)
    kernelE=ek*(-c.expm1(2*heat.mu*offset)/(2*heat.mu))
    kernelP=ek*(-c.expm1((1+2*heat.mu)*offset)/(2*(1+2*heat.mu)))
    Ad=shape['K_rows'][0]*X0-one/heat.k
    Ed=terminal['energy_defect_rows'][0]*e+one*(c.expm1(heat.delta*offset)/heat.delta)
    Ed+=one*(terminal['Kright']**2*kernelE)
    Pd=terminal['pressure_defect_rows'][0]*ep+one*(c.expm1(heat.prate*offset)/(2*heat.prate))
    Pd+=one*(terminal['Kright']**2*kernelP)
    A=[Ad]; E=[Ed]; Pr=[Pd]; Km=shape['K_defect_rows']; Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
        K_defect_rows=Km,K_squared_defect_rows=Qm)


def outer_power_stress_rows(heat,shape,defects,Z,q):
    return collar_stress_rows(heat,shape,defects,Z,q)


def outer_power_pressure_rows(heat,defects,shape,q):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    fullP=[defects['pressure_defect_rows'][0]+one/(2*heat.prate)]+defects['pressure_defect_rows'][1:]
    return pressure_y_rows(shape['K_rows'],fullP[0],heat.prate,heat.pressure_scale,q)


def outer_power_transport_identities():
    a,mu,y,z,KR=s.symbols('a mu y Z Kright',real=True); delta=2*a; k=1-a; p=1+delta; d=a-mu; r=1-mu
    ER=s.Function('same_full_right_E')(z); PR=s.Function('same_full_right_P')(z); XR=s.Function('same_native_right_X')(z)
    K=KR*s.exp(d*y)
    E=s.exp(delta*y)*ER+KR**2*s.exp(2*d*y)*(1-s.exp(2*mu*y))/(2*mu)
    PP=s.exp(p*y)*PR+KR**2*s.exp(2*d*y)*(1-s.exp((1+2*mu)*y))/(2*(1+2*mu))
    X=1/r+(XR-1/r)*s.exp(-r*y)
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Outer power identity failed: '+name)
        proofs[name]=True
    zero('native_forward_X_ODE',s.diff(X,y)-(1-r*X))
    zero('same_cumulative_A_FTC',s.diff(K*X,y)-(K-k*K*X))
    zero('same_full_energy_FTC',s.diff(E,y)-(delta*E-K*K))
    zero('same_full_pressure_FTC',s.diff(PP,y)-(p*PP-K*K/2))
    zero('right_energy_datum',E.subs(y,0)-ER)
    zero('right_pressure_datum',PP.subs(y,0)-PR)
    zero('stable_energy_integral_factoring',s.exp(delta*y)*(s.exp(-2*mu*y)-1)/(2*mu)
         -s.exp(2*d*y)*(1-s.exp(2*mu*y))/(2*mu))
    zero('stable_pressure_integral_factoring',s.exp(p*y)*(s.exp(-(1+2*mu)*y)-1)/(1+2*mu)
         -s.exp(2*d*y)*(1-s.exp((1+2*mu)*y))/(1+2*mu))
    zero('actual_constant_kappa_minus2',delta-2*d-2*mu)
    C,qright=s.symbols('same_pressure_scale qright',real=True)
    absolute=-C*s.exp(-p*(qright+y))*PP
    zero('same_absolute_pressure_derivative',s.diff(absolute,y)-C*s.exp(-p*(qright+y))*K*K/2)
    length,phase=s.symbols('original_Lrel_minus4 phase',real=True)
    zero('original_phase_to_ordinary_q_factor',s.diff(K.subs(y,length*(phase-1)),phase)-length*s.diff(K,y).subs(y,length*(phase-1)))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
        complete_quadratic_future_restored_from_actual_angular_left=True,
        native_forward_angular_history_retained=True,
        actual_source_K_Z_exact_zero=True,full_A_E_P_axial_derivatives_retained=True,
        stable_factored_remaining_integrals=True,source_caps_used_as_fields=False)


def outer_power_angular_join_binding():
    """Actual ordinary-q row implementations, arbitrary full terminal data."""
    a,mu,z,C,wait,Ts=s.symbols('a mu Z C wait Ts',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda v:s.exp(v)-1,mpf=lambda v:s.Rational(str(v)))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,S=s.Symbol('same_positive_inverse_Rtail'),delta=2*a,k=1-a,
        prate=1+2*a,pressure_scale=C,steep=SimpleNamespace(wait=wait,Ts=Ts))
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,c=c,IntervalTaylor=stub,axial_derivative=lambda v:s.diff(v,z),
        decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate)
    asts=SourceAST()
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('collar_stress_C3','collar_stress_rows'),('power_angular_C4','quotient_log_rates'),
        ('angular_stress_C3','angular_shape'),('angular_stress_C3','angular_X_rows'),
        ('angular_stress_C3','angular_defect_rows'),('angular_stress_C3','angular_stress_rows'),
        ('angular_stress_C3','angular_pressure_rows'),('outer_power_stress_C3','outer_power_shape'),
        ('outer_power_stress_C3','outer_power_X_rows'),('outer_power_stress_C3','outer_power_defect_rows'),
        ('outer_power_stress_C3','outer_power_stress_rows'),('outer_power_stress_C3','outer_power_pressure_rows')):
        if name=='quotient_log_rates':
            fn=asts.method(stem,name); fn.decorator_list=[];fn.body=[n for n in fn.body if not isinstance(n,ast.If)]
            exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual quotient algebra>','exec'),env)
        else:asts.replay(stem,name,env)
    asts.replay('heat_pressure_C4','pressure_y_rows',env,True)
    origin=dict(KR=s.Symbol('same_original_KR',positive=True),
        energy_defect_rows=[s.Function('same_entry_Edef')(z)],
        pressure_defect_rows=[s.Function('same_entry_Pdef')(z)])
    packet=dict(swirl_factor_one_plus_h_Taylor=s.Integer(1),actual_angular_bump_y_derivatives=[s.Integer(0)]*5)
    ag=env['angular_shape'](heat,origin,s.Integer(-4),packet)
    XR=s.Function('same_original_power_angular_X')(z)
    ax=env['angular_X_rows'](heat,ag,XR)
    kernels=dict(energy=s.Function('same_original_full_beta_E')(z),pressure=s.Function('same_original_full_beta_P')(z))
    ad=env['angular_defect_rows'](heat,origin,s.Integer(-4),ag,XR,kernels)
    terminal=dict(Kright=ag['K_rows'][0],energy_defect_rows=[ad['energy_defect_rows'][0]],
        pressure_defect_rows=[ad['pressure_defect_rows'][0]])
    ps=env['outer_power_shape'](heat,terminal,s.Integer(0))
    px=env['outer_power_X_rows'](heat,XR)
    pd=env['outer_power_defect_rows'](heat,terminal,s.Integer(0),ps,XR)
    q=-wait-Ts-6
    left=env['outer_power_stress_rows'](heat,ps,pd,z,q)
    right=env['angular_stress_rows'](heat,ag,ad,z,q)
    checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual outer power/angular join failed: '+name)
        checks[name]=True
    for j in range(5):
        zero('actual_K_q'+str(j),ps['K_rows'][j]-ag['K_rows'][j])
        for key in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows'):
            zero(key+'_q'+str(j),pd[key][j]-ad[key][j])
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):zero(label+'_q'+str(j)+'_Z'+str(n),s.diff(left[label][j]-right[label][j],z,n))
    lp=env['outer_power_pressure_rows'](heat,pd,ps,q);rp=env['angular_pressure_rows'](heat,ad,ag,q)
    for j in range(5):
        for n in range(5-j):zero('pressure_q'+str(j)+'_Z'+str(n),s.diff(lp[j]-rp[j],z,n))
    return dict(identities=checks,actual_outer_power_angular_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified=True,
        arbitrary_actual_full_terminal_axial_functions_used=True,source_function_equality_not_interval_overlap=True,input_hashes=asts.hashes)


@source_precision
def outer_power_source_bridge(angular):
    asts=SourceAST();checks={}
    for flag in ('actual_native_theta_B_K_normalization_consumed',
        'actual_same_full_remaining_E_P_identified_by_ODE_and_entry_datum',
        'actual_original_selected_coefficient_functions_retained','actual_full_pressure_both_packet_routes_bound',
        'actual_power_angular_left_source_radius_join_verified'):
        if not angular.bridge['identities'][flag]:raise ValueError('Outer power common angular datum missing: '+flag)
        checks['consumed_'+flag]=True
    if not all(angular.bridge['actual_original_endpoint_beta_jets_exact_zero'].values()):
        raise ValueError('Original angular endpoint beta jets required')
    def syntax(stem,name,target,wanted,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    syntax('outer_power_stress_C3','terminal','point','self.angular.angular(Z,-4)')
    syntax('outer_power_stress_C3','terminal',"data['Kright']","point['angular_shape']['K_rows'][0][0]")
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    value=asts.expression('power_angular_C4','data','value')
    items={kw.arg:kw.value for kw in value.keywords}
    for key,wanted in (('coeff','coeff'),('flatten_exit_X',"f['angular_Taylor']"),
        ('flatten_exit_pressure',"f['pressure']"),('full_angular_change','full_change')):
        if ast.dump(items[key])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Original power data route changed: '+key)
    syntax('power_angular_C4','data','full_change',"(dj*(2*self.weights['E'])+dj*dj*self.weights['F'])*c.exp(-2*self.mu*center)",True)
    syntax('power_angular_C4','power','y','(self.Lrel-4)*phase')
    syntax('power_angular_C4','power','D','4+(self.Lrel-4)*(1-phase)')
    nativeX=syntax('power_angular_C4','power','X',"one/self.rate+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    nativeE=syntax('power_angular_C4','power','energy',"(data['post']+data['full_angular_change'])*(c.exp(-2*self.mu*D)/2)")
    extraE=syntax('power_angular_C4','power','energy','decay_integral(c,2*self.mu,D)/2',True)
    theta=syntax('power_angular_C4','power','theta','one*(c.exp(-100*self.bp-self.bp*y)/2)')
    deltaP=syntax('power_angular_C4','power','deltaP','decay_integral(c,self.prate,y)*(self.flatten.Ev2*c.exp(-100*self.prate)/8)')
    nativeP=syntax('power_angular_C4','power','pressure',"data['flatten_exit_pressure']['P_over_Pstar_squared']+deltaP")
    syntax('power_angular_C4','_packet','mixed','flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    syntax('outer_power_stress_C3','power','mixed[P]',"{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}")
    syntax('outer_power_stress_C3','power','fields[P]','[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
    for target,wanted in (('length','self.outer.Lrel-4'),('y','length*(phase-1)'),
        ('q','-self.angular.entry.steep.wait-self.angular.entry.steep.Ts-6+y')):
        syntax('outer_power_stress_C3','power',target,wanted)
    a,mu,u,z,L,Ev2=s.symbols('a mu ordinary_q_offset Z Lrel same_Ev0_squared',real=True)
    r=1-mu;po=1+2*mu;bp=s.Rational(1,2)+mu
    Xf=s.Function('actual_flatten_Xf')(z);Pf=s.Function('actual_flatten_Pf')(z)
    Post=s.Function('same_full_Post')(z);JE=s.Function('actual_full_beta_E')(z)
    c=SimpleNamespace(exp=s.exp,mpf=lambda val:s.Rational(str(val)))
    obj=SimpleNamespace(mu=mu,rate=r,prate=po,bp=bp,flatten=SimpleNamespace(Ev2=Ev2))
    env=dict(c=c,self=obj,one=s.Integer(1),y=u,D=L-u,
        data=dict(flatten_exit_X=Xf,flatten_exit_pressure={'P_over_Pstar_squared':Pf},post=Post,full_angular_change=JE),
        decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate)
    X=asts.evaluate(nativeX,env);energy=asts.evaluate(nativeE,env)+asts.evaluate(extraE,env)
    th=asts.evaluate(theta,env);env['deltaP']=asts.evaluate(deltaP,env);pressure=asts.evaluate(nativeP,env)
    expressions=dict(actual_native_power_X_ODE_AST_replayed=s.diff(X,u)-(1-r*X),
        actual_native_power_energy_ODE_AST_replayed=s.diff(energy,u)-(2*mu*energy-s.Rational(1,2)),
        actual_native_power_absolute_pressure_ODE_AST_replayed=s.diff(pressure,u)-Ev2*th**2/2,
        actual_native_power_theta_rate_AST_replayed=s.diff(th,u)+bp*th)
    for name,value in expressions.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual outer power source ODE failed: '+name)
        checks[name]=True
    path=HERE/(PREFIX+'global_physical_assembly.py');envr={}
    for node in ast.parse(path.read_text(encoding='utf8')).body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):envr[target.id]=ast.literal_eval(node.value)
    fn=asts.method('global_physical_assembly','radius')
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual outer power radius>','exec'),envr)
    phase,wait,Ts,origin=s.symbols('phase wait Ts logRp',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    out=SimpleNamespace(Lrel=L);steep=SimpleNamespace(outer=out,wait=wait,Ts=Ts)
    qr=envr['radius'](assembly,'outer_power',phase,{},out)[0]-envr['radius'](assembly,'waiting',s.Integer(1),{},steep)[0]
    if s.simplify(qr-(-wait-Ts-2-L+(L-4)*phase))!=0:raise ArithmeticError('Original outer power phase/radius changed')
    right=envr['radius'](assembly,'outer_power',s.Integer(1),{},out)[0]
    left=envr['radius'](assembly,'outer_angular',s.Integer(-4),{},out)[0]
    if s.simplify(right-left)!=0:raise ArithmeticError('Native outer power/angular radius join failed')
    name=PREFIX+'power_angular_C4_check.json';native=json.loads((HERE/name).read_bytes())
    for flag in ('new_provider_exact_power_angular_X','new_provider_exact_power_angular_energy',
        'new_provider_exact_power_angular_pressure','new_provider_exact_power_angular_theta',
        'power_angular_s_minus4_full_future_and_empty_past_supports',
        'same_primitive_ODEs_and_original_flat_beta_identify_mixed_interface_jets'):
        if not native['functional_production_source_identities'][flag]:raise ValueError('Original power/angular source endpoint missing: '+flag)
        checks['consumed_'+flag]=True
    asts.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    join=outer_power_angular_join_binding();asts.hashes.update(join['input_hashes'])
    checks['actual_native_power_complete_history_identified_by_ODE_and_angular_datum']=True
    checks['actual_native_power_theta_equals_same_B_times_K_by_rate_and_right_datum']=True
    checks['actual_original_phase_to_ordinary_q_radius_verified']=True
    checks['actual_outer_power_angular_right_radius_field_pressure_join_verified']=True
    checks['actual_selected_full_angular_future_retained']=True
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,
        actual_outer_power_angular_formula_join=join,input_hashes=asts.hashes,
        whole_original_Lrel_minus4_retained=True,original_forward_angular_history_not_replaced=True,
        source_caps_used_as_fields=False,left_flatten_stress_companion_constructed=False)


class CompliantOuterPowerStressC3:
    @source_precision
    def __init__(self):
        self.angular=CompliantAngularStressC3();self.heat=self.angular.heat;self.outer=self.angular.outer
        self.ctx=self.heat.ctx;self.family,self.source=self.angular.family,self.angular.source
        self.hashes=dict(self.angular.hashes);self.cache={}
        name=PREFIX+'angular_stress_C3_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_angular_similarity_stress_recovered']:
            raise ValueError('Accepted angular complete moment datum required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Outer power prerequisite changed: '+source)
        self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=outer_power_transport_identities();self.bridge=outer_power_source_bridge(self.angular)
        self.hashes.update(self.bridge['input_hashes']);self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def terminal(self,Z):
        key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.cache:
            point=self.angular.angular(Z,-4)
            data={name:[point['angular_future_defect_zeroth_Taylor'][name]]
                for name in ('energy_defect_rows','pressure_defect_rows')}
            data['Kright']=point['angular_shape']['K_rows'][0][0]
            self.cache[key]=data
        return self.cache[key]

    @source_precision
    def power(self,Z,phase):
        c=self.ctx;Z=c.mpf(Z);phase=c.mpf(phase)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(phase)[0]<0 or endpoints(phase)[1]>1:
            raise ValueError('Original outer_power domain Z[-1,1], phase[0,1] required')
        length=self.outer.Lrel-4; y=length*(phase-1)
        q=-self.angular.entry.steep.wait-self.angular.entry.steep.Ts-6+y
        point=dict(self.outer.power(Z,phase));terminal=self.terminal(Z)
        shape=outer_power_shape(self.heat,terminal,y);Xrows=outer_power_X_rows(self.heat,point['angular_Taylor'])
        defects=outer_power_defect_rows(self.heat,terminal,y,shape,Xrows[0])
        stress=outer_power_stress_rows(self.heat,shape,defects,Z,q);pressure=outer_power_pressure_rows(self.heat,defects,shape,q)
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']);fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_mixed_bounds']=mixed[P]
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        one=IntervalTaylor.constant(c,1,5);fullE=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        energy=quotient_rows(fullE,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        point.update(physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
            outer_power_phase=phase,outer_power_offset_from_angular_left=y,outer_power_offset_from_Rtail=q,
            outer_power_shape=shape,
            outer_power_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):stress[label][j][n]*math.factorial(n)
                for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
            outer_power_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(stress[label])]
                for label in ('theta','axial','theta_inertial','theta_shear')},
            outer_power_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
            pressure_over_Pstar_squared_Taylor=pressure[0],pressure_y_derivative_axial5_Taylor=pressure,
            original_native_energy_Taylor=point['energy_Taylor'],same_full_energy_Taylor=energy[0],
            outer_power_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
            actual_original_outer_power_similarity_stress_recovered=True,actual_outer_power_absolute_pressure_same_source_mixed4_available=True,
            outer_power_angular_stress_mixed3_and_pressure_mixed4_join_verified=True,
            original_selected_angular_full_future_retained=True,actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,
            original_power_length=length,ordinary_q_derivatives_not_unscaled_phase_derivatives=True,
            outer_power_physical_decomposition_constructed=False,outer_power_cone_certified=False,
            left_flatten_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    @source_precision
    def report(self):
        keys=('Z','outer_power_phase','outer_power_offset_from_angular_left','outer_power_offset_from_Rtail','outer_power_shape',
            'outer_power_similarity_stress_mixed3_factored','outer_power_future_defect_zeroth_Taylor',
            'pressure_y_derivative_axial5_Taylor','outer_power_meridional_moments_and_velocities',
            'original_power_length','ordinary_q_derivatives_not_unscaled_phase_derivatives')
        def packet(z,phase):
            point=self.power(z,phase);return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            scope='Whole original outer_power phase[0,1],Z[-1,1], actual Lrel-4; similarity stress mixed3 and absolute pressure mixed4',
            samples=[packet(z,t) for z,t in (('0','0'),('.5','.25'),('-.5','.75'),('0','1'))],
            whole_outer_power=packet([-1,1],[0,1]),transport_identities=self.proof,source_bridge=self.bridge,input_hashes=self.hashes,
            actual_original_outer_power_similarity_stress_recovered=True,actual_outer_power_absolute_pressure_same_source_mixed4_available=True,
            outer_power_angular_stress_mixed3_and_pressure_mixed4_join_verified=True,
            original_selected_angular_full_future_retained=True,actual_K_Z_exact_zero_and_full_moment_axial_dependence_retained=True,
            outer_power_physical_decomposition_constructed=False,outer_power_cone_certified=False,
            left_flatten_stress_companion_constructed=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False)


@source_precision
def run():
    result=CompliantOuterPowerStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original preangular power full moments, stress and absolute pressure generated; physical/cone/flatten pending',flush=True)
    return result


if __name__=='__main__':run()
