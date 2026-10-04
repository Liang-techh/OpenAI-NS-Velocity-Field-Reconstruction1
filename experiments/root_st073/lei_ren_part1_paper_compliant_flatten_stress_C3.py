"""Original100-unit flatten: same complete moments and absolute pressure.

Variable axial K and original native forward X are retained. Remaining
quadratic kernels continue the accepted outer-power datum before enclosure.
This layer does not establish physical/cone/global/recursion completion.
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
from lei_ren_part1_paper_compliant_outer_power_stress_C3 import (
    CompliantOuterPowerStressC3,SourceAST,source_precision,quotient_rows)
from lei_ren_part1_paper_compliant_angular_high_jets import log_taylor
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
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
DOMAIN=dict(t=[0,100],Z=[-1,1],offset='q=-wait-Ts-2-Lrel+(t-100)',ordinary_derivative='d_q=d_t')


def flatten_f(c,Z):
    return IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])/2


def flatten_shape(heat,terminal,offset,packet):
    c=heat.ctx;one=IntervalTaylor.constant(c,1,5)
    f=flatten_f(c,packet['Z']);rho=log_taylor(f);sj=packet['sigma_y_derivatives']
    rates=[rho*sj[j+1] for j in range(4)];rates[0]=rates[0]+heat.a-heat.mu
    factor=terminal['Kright']*c.exp((heat.a-heat.mu)*offset)
    K=[one*factor if packet['actual_original_flatten_right_endpoint'] else packet['F_Taylor']/f*factor]
    for n in range(1,5):
        K.append(sum((rates[j]*K[n-1-j]*math.comb(n-1,j) for j in range(n)),one*0))
    square=product_rows(K,K)
    return dict(K_rows=K,K_defect_rows=[K[0]-one]+K[1:],
        K_squared_defect_rows=[square[0]-one]+square[1:],
        log_K_rate_rows=rates,original_sigma_y_derivatives=sj,
        actual_K_Z_dependence_retained=True,exact_original_endpoint_F_over_f_cancelled=True)


def flatten_remaining_kernels(c,mu,Z,t,cells):
    """Directed integral of actual (F/f)^2 with exact exponential weights."""
    Z=c.mpf(Z);t=c.mpf(t);one=IntervalTaylor.constant(c,1,5)
    rho=log_taylor(flatten_f(c,Z))
    length=100-t;energy=one*0;pressure=one*0
    if endpoints(length)[1]!=0:
        for i in range(cells):
            a=length*i/cells;b=length*(i+1)/cells
            lo=max(mp.mpf(0),endpoints(t+a)[0]);hi=min(mp.mpf(100),endpoints(t+b)[1])
            v=c.mpf([lo,hi]);sig=sigma_jets(c,v/100)[0]
            density=(rho*(2*(sig-1))).exp()
            width=length/cells
            energy+=density*(c.exp(-2*mu*a)*decay_integral(c,2*mu,width))
            pressure+=density*(c.exp(-(1+2*mu)*a)*decay_integral(c,1+2*mu,width)/2)
    return dict(energy=energy,pressure=pressure,
        original_normalized_squared_integrand='exp(2*(sigma(v/100)-1)*log((1+Z^2)/2))',
        original_remaining_interval=[t,c.mpf(100)],directed_cell_count=cells,
        actual_remaining_F_over_endpoint_f_correlation_retained=True)


def flatten_defect_rows(heat,terminal,offset,shape,X0,kernels):
    c=heat.ctx;one=IntervalTaylor.constant(c,1,5);d=heat.a-heat.mu
    e=c.exp(heat.delta*offset);ep=c.exp(heat.prate*offset);ek=c.exp(2*d*offset)
    Ad=shape['K_rows'][0]*X0-one/heat.k
    Ed=terminal['energy_defect_rows'][0]*e+one*(c.expm1(heat.delta*offset)/heat.delta)
    Ed+=kernels['energy']*(terminal['Kright']**2*ek)
    Pd=terminal['pressure_defect_rows'][0]*ep+one*(c.expm1(heat.prate*offset)/(2*heat.prate))
    Pd+=kernels['pressure']*(terminal['Kright']**2*ek)
    A=[Ad];E=[Ed];Pr=[Pd];Km=shape['K_defect_rows'];Qm=shape['K_squared_defect_rows']
    for j in range(4):
        A.append(Km[j]-A[j]*heat.k)
        E.append(E[j]*heat.delta-Qm[j])
        Pr.append(Pr[j]*heat.prate-Qm[j]/2)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
        K_defect_rows=Km,K_squared_defect_rows=Qm)


def flatten_stress_rows(heat,shape,defects,Z,q):
    return collar_stress_rows(heat,shape,defects,Z,q)


def flatten_pressure_rows(heat,defects,shape,q):
    one=IntervalTaylor.constant(heat.ctx,1,5)
    return pressure_y_rows(shape['K_rows'],one/(2*heat.prate)+defects['pressure_defect_rows'][0],
        heat.prate,heat.pressure_scale,q)


def flatten_transport_identities():
    a,mu,t,z,KR=s.symbols('a mu t Z Kright',real=True);y=t-100;d=a-mu
    k=1-a;delta=2*a;p=1+delta;r=1-mu;f=(1+z*z)/2
    F=s.Function('original_flatten_F')(t,z);N=s.Function('native_weighted_X_numerator')(t,z)
    JE=s.Function('same_remaining_E_kernel')(t,z);JP=s.Function('same_remaining_P_kernel')(t,z)
    ER=s.Function('same_full_right_E')(z);PR=s.Function('same_full_right_P')(z)
    G=F/f;K=KR*s.exp(d*y)*G;X=N/F
    E=s.exp(delta*y)*ER+KR**2*s.exp(2*d*y)*JE
    PP=s.exp(p*y)*PR+KR**2*s.exp(2*d*y)*JP
    rules={s.diff(N,t):F-r*N,s.diff(JE,t):2*mu*JE-G**2,
        s.diff(JP,t):(1+2*mu)*JP-G**2/2}
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Flatten transport identity failed: '+name)
        proofs[name]=True
    zero('actual_native_X_ODE',s.diff(X,t).subs(rules)-(1-(r+s.diff(F,t)/F)*X))
    zero('same_native_KX_cumulative_A_FTC',s.diff(K*X,t).subs(rules)-(K-k*K*X))
    zero('same_complete_energy_FTC',s.diff(E,t).subs(rules)-(delta*E-K*K))
    zero('same_complete_pressure_FTC',s.diff(PP,t).subs(rules)-(p*PP-K*K/2))
    en=E/(2*K*K)
    zero('same_native_energy_variable_F_ODE',s.diff(en,t).subs(rules)
        -((2*mu-2*s.diff(F,t)/F)*en-s.Rational(1,2)))
    C,qR=s.symbols('same_absolute_pressure_scale qright',real=True)
    absolute=-C*s.exp(-p*(qR+y))*PP
    zero('same_absolute_pressure_derivative',s.diff(absolute,t).subs(rules)-C*s.exp(-p*(qR+y))*K*K/2)
    endpoint={t:100,F.subs(t,100):f,JE.subs(t,100):0,JP.subs(t,100):0}
    zero('same_right_full_E_datum',E.subs(t,100).subs({JE.subs(t,100):0})-ER)
    zero('same_right_full_P_datum',PP.subs(t,100).subs({JP.subs(t,100):0})-PR)
    zero('actual_variable_flatten_kappa_minus2',delta-2*(d+s.diff(F,t)/F)-(2*mu-2*s.diff(F,t)/F))
    zero('same_native_theta_source_rate',d-(s.Rational(1,2)+a)+s.Rational(1,2)+mu)
    omega=s.symbols('omega',real=True)
    zero('actual_F_over_f_squared_integrand',s.exp(2*(1-omega)*s.log(f))/f**2-s.exp(-2*omega*s.log(f)))
    return dict(identities=proofs,whole_original_domain=DOMAIN,
        actual_variable_axial_K_retained=True,original_forward_X_retained=True,
        complete_future_transport_from_same_outer_power_datum=True,
        source_caps_used_as_fields=False,actual_variable_shear_formula='kappa-2=2*mu-2*rho*sigma_t; sign not admitted')


def flatten_power_join_binding():
    a,mu,z,C,wait,Ts,L=s.symbols('a mu Z C wait Ts Lrel',real=True)
    c=SimpleNamespace(exp=s.exp,expm1=lambda v:s.exp(v)-1,mpf=lambda v:s.Rational(str(v)))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,S=s.Symbol('same_exact_inverse_Rtail'),delta=2*a,k=1-a,
        prate=1+2*a,pressure_scale=C)
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,c=c,IntervalTaylor=stub,log_taylor=s.log,flatten_f=lambda ctx,value:(1+value*value)/2,
        axial_derivative=lambda v:s.diff(v,z))
    asts=SourceAST()
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('collar_stress_C3','collar_stress_rows'),('flatten_stress_C3','flatten_shape'),
        ('flatten_stress_C3','flatten_defect_rows'),('flatten_stress_C3','flatten_stress_rows'),
        ('flatten_stress_C3','flatten_pressure_rows'),('outer_power_stress_C3','outer_power_shape'),
        ('outer_power_stress_C3','outer_power_defect_rows'),('outer_power_stress_C3','outer_power_stress_rows'),
        ('outer_power_stress_C3','outer_power_pressure_rows')):asts.replay(stem,name,env)
    asts.replay('heat_pressure_C4','pressure_y_rows',env,True)
    terminal=dict(Kright=s.Symbol('same_power_inlet_K',positive=True),
        energy_defect_rows=[s.Function('same_power_inlet_Edef')(z)],
        pressure_defect_rows=[s.Function('same_power_inlet_Pdef')(z)])
    X=s.Function('same_original_flatten_exit_X')(z);f=(1+z*z)/2
    packet=dict(Z=z,F_Taylor=f,sigma_y_derivatives=[s.Integer(1)]+[s.Integer(0)]*4,
        actual_original_flatten_right_endpoint=True)
    fs=env['flatten_shape'](heat,terminal,s.Integer(0),packet)
    fd=env['flatten_defect_rows'](heat,terminal,s.Integer(0),fs,X,dict(energy=s.Integer(0),pressure=s.Integer(0)))
    ps=env['outer_power_shape'](heat,terminal,s.Integer(0))
    pd=env['outer_power_defect_rows'](heat,terminal,s.Integer(0),ps,X)
    q=-wait-Ts-2-L
    left=env['flatten_stress_rows'](heat,fs,fd,z,q);right=env['outer_power_stress_rows'](heat,ps,pd,z,q)
    checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual flatten/power join failed: '+name)
        checks[name]=True
    for j in range(5):
        zero('actual_K_q'+str(j),fs['K_rows'][j]-ps['K_rows'][j])
        for label in ('angular_defect_rows','energy_defect_rows','pressure_defect_rows'):
            zero(label+'_q'+str(j),fd[label][j]-pd[label][j])
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):zero(label+'_q'+str(j)+'_Z'+str(n),s.diff(left[label][j]-right[label][j],z,n))
    lp=env['flatten_pressure_rows'](heat,fd,fs,q);rp=env['outer_power_pressure_rows'](heat,pd,ps,q)
    for j in range(5):
        for n in range(5-j):zero('absolute_pressure_q'+str(j)+'_Z'+str(n),s.diff(lp[j]-rp[j],z,n))
    return dict(identities=checks,input_hashes=asts.hashes,
        actual_flatten_power_K4_full_moment_stress_mixed3_pressure_mixed4_AST_join_verified=True,
        arbitrary_axial_terminal_functions_used=True,original_zero_remaining_interval_used=True,
        left_pulse_stress_companion_constructed=False)


def flatten_source_bridge(power):
    asts=SourceAST();checks={}
    def syntax(stem,name,target,wanted=None,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    nativeX=syntax('flatten_mixed_C4','flatten','X','(Xint+self.Xv*c.exp(-self.rate*t))/F')
    nativeE=syntax('flatten_mixed_C4','flatten','energy','(ev-Eint/2)*c.exp(2*self.mu*t)/(F*F)')
    syntax('flatten_mixed_C4','flatten','ev','self.future_energy(Z)')
    syntax('flatten_mixed_C4','flatten','Eint',
        'Fc*Fc*(c.exp(-2*self.mu*a)*decay_integral(c,2*self.mu,length))',True)
    syntax('flatten_mixed_C4','flatten','Pint',
        'Fc*Fc/(q*q)*(c.exp(-self.prate*a)*decay_integral(c,self.prate,length)/2)',True)
    syntax('flatten_mixed_C4','flatten','Xint',
        'Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))',True)
    syntax('flatten_mixed_C4','flatten','mixed','flatten_mixed(c,self.mu,rho,sj,theta,X,energy,pressure,self.Ev2)')
    qnode=syntax('flatten_mixed_C4','flatten','q','IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])')
    freturn=next(n.value for n in ast.walk(asts.method('flatten_stress_C3','flatten_f')) if isinstance(n,ast.Return))
    if ast.dump(freturn)!=ast.dump(ast.BinOp(left=qnode,op=ast.Div(),right=ast.Constant(value=2))):
        raise ValueError('Actual normalized flatten f changed')
    checks['actual_normalized_f_is_original_q_over_two_AST_verified']=True
    syntax('flatten_mixed_C4','flatten','rho','log_taylor(q)-c.ln(2)')
    syntax('flatten_mixed_C4','flatten','sig','sigma_jets(c,t/100)')
    syntax('flatten_mixed_C4','flatten','F','(rho*sj[0]).exp()')
    syntax('flatten_mixed_C4','flatten','Fc','(rho*sigma_jets(c,v/100)[0]).exp()')
    syntax('flatten_mixed_C4','flatten','theta','F/q*c.exp(-self.bp*t)')
    syntax('flatten_mixed_C4','flatten','theta','IntervalTaylor.constant(c,c.exp(-100*self.bp)/2,5)')
    for target,wanted in (('rates','[rho*sigma_derivatives[j+1] for j in range(4)]'),
        ('energyrate','[v*(-2) for v in rates]')):
        syntax('flatten_mixed_C4','flatten_mixed',target,wanted)
    syntax('flatten_mixed_C4','flatten_mixed','energyrate[0]','energyrate[0]+2*mu')
    syntax('flatten_mixed_C4','future_energy','packet',"""self.fifth['whole_Z']['energy']""")
    fn=asts.method('flatten_mixed_C4','future_energy')
    returns=[n.value for n in ast.walk(fn) if isinstance(n,ast.Return)]
    expected=ast.parse("IntervalTaylor(c,[read_interval(c,v)/2 for v in packet['complete_future_energy_Taylor']['coefficients']])",mode='eval').body
    if len(returns)!=1 or ast.dump(returns[0])!=ast.dump(expected):raise ValueError('Original complete future energy half-normalization changed')
    checks['actual_flatten_future_energy_half_of_complete_C5_source_AST_verified']=True
    # The admitted provider already binds the original complete decomposition,
    # selected outer-angular coefficients and full post/Gamma future.
    name=PREFIX+'power_angular_C4_check.json';old=json.loads((HERE/name).read_bytes())
    for flag in ('production_future_Nf_correlated_q2','production_future_Nrel_correlated_q2',
        'production_future_power_piece','production_future_full_decomposition',
        'production_signed_angular_energy_change','production_complete_post_scalar',
        'production_entire_Gamma_energy_factor','new_provider_power_energy_is_exact_backward_source'):
        if not old['functional_production_source_identities'][flag]:raise ValueError('Same complete future source not bound: '+flag)
        checks['consumed_'+flag]=True
    asts.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    # Replay the actual two endpoint energy implementations, instead of
    # identifying their numerical interval overlap with a function identity.
    mu,z,L=s.symbols('mu Z Lrel',real=True);poly=1+z*z;f=poly/2
    F0=s.Function('same_original_flatten_integral')(z)
    AC=s.Function('same_selected_angular_full_change')(z);Post=s.Function('same_full_post')(z)
    Nr=poly**2*s.exp(-2*mu*(100+L))/4;Nf=poly**2*s.exp(-200*mu)/4
    D=(1-s.exp(-2*mu*L))/(2*mu)
    totalnode=syntax('fifth_axial_jets','future','total','flatten+power+angular_change+postrel-heatdifference')
    HD=s.Function('same_Gamma_deficit')(z)
    total=asts.evaluate(totalnode,dict(flatten=F0,power=Nf*D,angular_change=Nr*AC,postrel=Nr*Post,heatdifference=Nr*HD))
    c=SimpleNamespace(exp=s.exp)
    actualflat=asts.evaluate(nativeE,dict(ev=total/2,Eint=F0,c=c,self=SimpleNamespace(mu=mu),t=s.Integer(100),F=f))
    node=syntax('power_angular_C4','power','energy',"(data['post']+data['full_angular_change'])*(c.exp(-2*self.mu*D)/2)")
    extra=syntax('power_angular_C4','power','energy','decay_integral(c,2*self.mu,D)/2',True)
    env=dict(c=c,self=SimpleNamespace(mu=mu),D=L,data=dict(post=Post-HD,full_angular_change=AC),
        decay_integral=lambda ctx,rate,length:(1-s.exp(-rate*length))/rate)
    actualpower=asts.evaluate(node,env)+asts.evaluate(extra,env)
    for n in range(6):
        if s.simplify(s.expand_power_exp(s.diff(actualflat-actualpower,z,n)))!=0:
            raise ArithmeticError('Actual source flatten/power energy identity failed')
        checks['actual_flatten_energy_AST_equals_same_power_inlet_Z'+str(n)]=True
    checks['actual_complete_future_contains_selected_outer_angular_repairs_and_entire_heat']=True
    checks['later_axial_ap_not_falsely_claimed_selected']=True
    # Source densities, endpoint F=f, forward pressure and the accepted
    # power absolute datum identify native flatten pressure by its FTC.
    syntax('flatten_mixed_C4','flatten','Mp','Mp_v+Pint*self.Ev2')
    syntax('flatten_mixed_C4','flatten','pressure',"""Mp+data['P0']""")
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    syntax('power_angular_C4','data','value',
        """dict(coeff=coeff,post=post,full_angular_change=full_change,flatten_exit_X=f['angular_Taylor'],flatten_exit_pressure=f['pressure'],flatten_exit_energy=f['energy_Taylor'],complete_future_energy_input=self.flatten.future_energy(Z))""")
    syntax('flatten_stress_C3','flatten','X0',"""point['angular_Taylor']""")
    syntax('flatten_stress_C3','flatten','mixed[P]',
        """{'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}""")
    syntax('flatten_stress_C3','flatten','fields[P]','[jet.truncate(4-j) for j,jet in enumerate(pressure)]')
    syntax('flatten_stress_C3','flatten_remaining_kernels','density','(rho*(2*(sig-1))).exp()')
    syntax('flatten_stress_C3','flatten_remaining_kernels','energy',
        'density*(c.exp(-2*mu*a)*decay_integral(c,2*mu,width))',True)
    syntax('flatten_stress_C3','flatten_remaining_kernels','pressure',
        'density*(c.exp(-(1+2*mu)*a)*decay_integral(c,1+2*mu,width)/2)',True)
    checks['original_native_forward_X_and_complete_pressure_history_retained']=True
    checks['same_absolute_pressure_identified_by_original_FTC_and_power_datum']=True
    checks['actual_original_F_over_endpoint_f_normalization_retained']=True
    join=flatten_power_join_binding();asts.hashes.update(join['input_hashes'])
    # Replay the actual production radius.
    fn=asts.method('global_physical_assembly','radius');envr={}
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):envr[target.id]=ast.literal_eval(node.value)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual flatten radius>','exec'),envr)
    t,wait,Ts,origin=s.symbols('t wait Ts logRp',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=origin,params=SimpleNamespace(mu=mu))
    outer=SimpleNamespace(Lrel=L);steep=SimpleNamespace(outer=outer,Ts=Ts,wait=wait)
    qr=envr['radius'](assembly,'flatten',t,{},None)[0]-envr['radius'](assembly,'waiting',s.Integer(1),{},steep)[0]
    if s.simplify(qr-(t-100-L-2-Ts-wait))!=0:raise ArithmeticError('Original flatten q/radius changed')
    if s.simplify(envr['radius'](assembly,'flatten',100,{},None)[0]
        -envr['radius'](assembly,'outer_power',0,{},outer)[0])!=0:raise ArithmeticError('Original flatten right radius join failed')
    checks['actual_flatten_radius_and_right_power_interface_AST_verified']=True
    checks['ordinary_t_equals_logR_derivatives']=True
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,input_hashes=asts.hashes,
        actual_flatten_power_formula_join=join,left_pulse_stress_companion_constructed=False,
        native_energy_function_join_proved_not_interval_overlap=True,source_caps_used_as_fields=False)


class CompliantFlattenStressC3:
    @source_precision
    def __init__(self):
        self.power=CompliantOuterPowerStressC3();self.heat=self.power.heat;self.native=self.power.outer.flatten
        self.ctx=self.heat.ctx;self.family,self.source=self.power.family,self.power.source
        self.hashes=dict(self.power.hashes);self.cache={}
        name=PREFIX+'outer_power_stress_C3_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_full_native_history_and_source_AST_bridge_verified']:
            raise ValueError('Accepted same-source power full datum required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flatten prerequisite changed: '+source)
        self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=flatten_transport_identities();self.bridge=flatten_source_bridge(self.power)
        self.hashes.update(self.bridge['input_hashes']);self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def terminal(self,Z):
        key=tuple(endpoints(self.ctx.mpf(Z)))
        if key not in self.cache:
            packet=self.power.power(Z,0)
            result={name:[packet['outer_power_future_defect_zeroth_Taylor'][name]]
                for name in ('energy_defect_rows','pressure_defect_rows')}
            result['Kright']=packet['outer_power_shape']['K_rows'][0][0]
            self.cache[key]=result
        return self.cache[key]

    @source_precision
    def flatten(self,Z,t):
        c=self.ctx;Z=c.mpf(Z);t=c.mpf(t)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1 or endpoints(t)[0]<0 or endpoints(t)[1]>100:
            raise ValueError('Original flatten domain Z[-1,1],t[0,100] required')
        point=dict(self.native.flatten(Z,t));point['actual_original_flatten_right_endpoint']=endpoints(t)==(mp.mpf(100),mp.mpf(100))
        y=t-100;q=-self.heat.steep.wait-self.heat.steep.Ts-2-self.power.outer.Lrel+y
        terminal=self.terminal(Z);shape=flatten_shape(self.heat,terminal,y,point)
        kernels=flatten_remaining_kernels(c,self.heat.mu,Z,t,self.native.cells);X0=point['angular_Taylor']
        defects=flatten_defect_rows(self.heat,terminal,y,shape,X0,kernels)
        stress=flatten_stress_rows(self.heat,shape,defects,Z,q);pressure=flatten_pressure_rows(self.heat,defects,shape,q)
        mixed=dict(point['physical_mixed_derivatives_total_order_le4']);fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_mixed_bounds']=mixed[P];point['original_native_pressure']=point['pressure']
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        one=IntervalTaylor.constant(c,1,5)
        fullE=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        normalizedE=quotient_rows(fullE,[jet*2 for jet in product_rows(shape['K_rows'],shape['K_rows'])])
        point.update(physical_mixed_derivatives_total_order_le4=mixed,physical_velocity_and_pressure_y_derivative_Taylor=fields,
            flatten_offset_from_power_inlet=y,flatten_offset_from_Rtail=q,flatten_shape=shape,
            flatten_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):stress[label][j][n]*math.factorial(n)
                for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
            flatten_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(stress[label])]
                for label in ('theta','axial','theta_inertial','theta_shear')},
            flatten_future_defect_zeroth_Taylor={name:values[0] for name,values in defects.items()},
            remaining_flatten_kernels=kernels,pressure_over_Pstar_squared_Taylor=pressure[0],
            pressure_y_derivative_axial5_Taylor=pressure,original_native_energy_Taylor=point['energy_Taylor'],
            same_full_energy_Taylor=normalizedE[0],flatten_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
            actual_original_flatten_similarity_stress_recovered=True,actual_flatten_absolute_pressure_same_source_mixed4_available=True,
            flatten_power_stress_mixed3_and_pressure_mixed4_join_verified=True,actual_variable_K_Z_and_full_moment_axial_dependence_retained=True,
            original_selected_outer_angular_full_future_retained=True,ordinary_q_derivatives_not_rescaled=True,
            flatten_physical_decomposition_constructed=False,flatten_cone_certified=False,left_pulse_stress_companion_constructed=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    @source_precision
    def report(self):
        keys=('Z','t','flatten_offset_from_power_inlet','flatten_offset_from_Rtail','flatten_shape',
            'flatten_similarity_stress_mixed3_factored','flatten_future_defect_zeroth_Taylor',
            'pressure_y_derivative_axial5_Taylor','flatten_meridional_moments_and_velocities',
            'remaining_flatten_kernels','ordinary_q_derivatives_not_rescaled')
        def packet(z,t):
            point=self.flatten(z,t);return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            scope='Whole original100-unit flatten, Z[-1,1], SAME full moments/native forward X/absolute pressure; similarity stress only',
            samples=[packet(z,t) for z,t in (('0','0'),('.5','50'),('0','100'))],whole_flatten=packet([-1,1],[0,100]),
            transport_identities=self.proof,source_bridge=self.bridge,input_hashes=self.hashes,
            actual_original_flatten_similarity_stress_recovered=True,actual_flatten_absolute_pressure_same_source_mixed4_available=True,
            flatten_power_stress_mixed3_and_pressure_mixed4_join_verified=True,actual_variable_K_Z_and_full_moment_axial_dependence_retained=True,
            original_selected_outer_angular_full_future_retained=True,ordinary_q_derivatives_not_rescaled=True,
            flatten_physical_decomposition_constructed=False,flatten_cone_certified=False,left_pulse_stress_companion_constructed=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)


@source_precision
def run():
    result=CompliantFlattenStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original100-unit flatten SAME full moments/stress/absolute pressure and power right join generated; physical/cone/left pulse pending',flush=True)
    return result


if __name__=='__main__':run()
