"""Original selected pulse-end full meridional stress in exact log sectors.

Both narrow supports and all five source moments remain. Unmaterializable
positive B,D,H factors are kept as source logs, not replaced by cap values.
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
from lei_ren_part1_paper_compliant_flatten_stress_C3 import (
    CompliantFlattenStressC3,SourceAST,source_precision,axial_derivative)
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_axial_pulse_field import raw_beta,backward_bump_weights
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(s=[-4,0],Z=[-1,1],ordinary_derivative='d_s=d_logR',
            centers=[-3,-1],width='.15')
PAPER_SHA256='8396b998dcf737cd6a7c12ff1019c6b2fd308a7e6f0907c0c6647a9b2ec64e8d'


def pulse_coefficients(delta,mu,z,C,Xp,Bh,m,n,e0,J,P):
    """Shape rows; modes=(R power,B power,D power,H power)/sqrt(2).

    B is the Z-independent pure-swirl reference, D=exp(log_end_scale),
    H=exp(-13*(1-mu)/mu-(1-mu)*s), C=1/(1+Z^2).
    """
    b=(1-delta)/2;k=1-delta/2;r=1-mu;bp=mu+mp.mpf('.5')
    d=1-z*z;L=1-z*z*delta;zero=C*0
    Cm=[C*v for v in m];Cn=[C*C*v for v in n]
    CE=[C*C*v for v in e0];CJ=[C*C*v for v in J]
    transport=[z*v*(1-delta)+d*axial_derivative(v) for v in Cm]
    mixed=[(z*v*(2*delta-1)-d*axial_derivative(v)) for v in Cn]
    nonlinear=product_rows(Bh,transport)
    equilibrium=(C*((mu-delta/2)/r)-z*axial_derivative(C)*(b/r))/L
    memory=(C*k-z*axial_derivative(C)*b)*(Xp-1/r)/L
    theta={
        'equilibrium':dict(mode=(.5,1,0,0),shape=[equilibrium,zero,zero,zero]),
        'signed_original_memory':dict(mode=(.5,1,0,1),shape=[memory,zero,zero,zero]),
        'meridional_transport':dict(mode=(.5,2,1,0),shape=[(C*transport[j]+mixed[j])/L for j in range(4)]),
        'radial_shear':dict(mode=(-.5,1,0,0),shape=[C*(-2*(1+mu)),zero,zero,zero])}
    axial={
        'full_energy_and_pressure':dict(mode=(.5,2,0,0),shape=[
            (z*CE[j]*(2*delta)-d*axial_derivative(CE[j])+z*P[j]*(2*(1+delta))-d*axial_derivative(P[j]))/L
            for j in range(4)]),
        'axial_transport':dict(mode=(.5,1,1,0),shape=[-C*Bh[j]/L for j in range(4)]),
        'nonlinear_meridional_transport':dict(mode=(.5,2,2,0),shape=[C*nonlinear[j]/L for j in range(4)]),
        'linear_axial_moment':dict(mode=(.5,1,1,0),shape=[
            (Cm[j]-z*axial_derivative(Cm[j]))*(1-delta)/(2*L) for j in range(4)]),
        'selected_backward_energy_loss':dict(mode=(.5,2,2,0),shape=[
            (z*CJ[j]*(-2*delta)+d*axial_derivative(CJ[j]))/L for j in range(4)]),
        'axial_radial_shear':dict(mode=(-.5,1,1,0),shape=[2*C*(Bh[j+1]-Bh[j]*bp) for j in range(4)])}
    for rows in (theta,axial):
        for part in rows.values():
            rp,bpwr,dp,hp=part['mode']
            part['source_logR_rate']=rp-bp*bpwr-r*hp
            part['full_derivative_rows']=shifted_rows(part['shape'],part['source_logR_rate'])
    return dict(theta=theta,axial=axial)


def original_formula_identities():
    """Bind the actual new coefficients to original full (3.16)-(3.18)."""
    asts=SourceAST();z,delta,mu,Xp=s.symbols('Z delta mu Xp',real=True)
    R,B,D,H=s.symbols('R B D H',positive=True)
    C=1/(1+z*z);r=1-mu;d=1-z*z;L=1-delta*z*z
    bh=s.Function('actual_Bhat')(z);m=s.Function('actual_m1hat')(z)
    n=s.Function('actual_m2hat')(z);ev=s.Function('actual_unperturbed_energy')(z)
    J=s.Function('actual_backward_energy_loss')(z);P=s.Function('same_signed_pressure_coefficient')(z)
    X=1/r+(Xp-1/r)*H;Ut=B*C;Uz=B*C*D*bh
    moments=dict(theta=s.sqrt(2)*R**s.Rational(3,2)*Ut*X,
        z=R*Ut*D*m,theta_z=s.sqrt(2)*R**s.Rational(3,2)*Ut**2*D*n,
        z_theta=R*Ut**2*(ev-D**2*J),p=s.Symbol('unused_raw_Mp'))
    env=dict(a=Ut,b=Uz,ay=-(s.Rational(1,2)+mu)*Ut,
        by=B*C*D*(s.Symbol('Bhat_s')-(s.Rational(1,2)+mu)*bh),
        dt=delta,z=z,R=R,root=s.sqrt(2*R),L=L,d=d,m=moments,
        mz={key:s.diff(v,z) for key,v in moments.items()},p=B**2*P,pz=B**2*s.diff(P,z))
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Full pulse stress identity failed: '+name)
        proofs[name]=True
    path=HERE/'lei_ren_part1_paper_mp_stress.py'
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    paperfn=next(v for v in ast.walk(ast.parse(path.read_text(encoding='utf8'))) if isinstance(v,ast.FunctionDef) and v.name=='evaluate_mp_stress')
    def paper_expression(target):
        values=[v.value for v in ast.walk(paperfn) if isinstance(v,ast.Assign) and any(ast.unparse(t)==target for t in v.targets)]
        if len(values)!=1:raise ValueError('Original full paper stress assignment changed: '+target)
        return values[0]
    transport=paper_expression('transport')
    env['transport']=asts.evaluate(transport,env)
    actual={}
    for target in ('Itheta','Iz','Stheta','Sz'):
        node=paper_expression(target)
        # Optional exact shear overrides are not used in this identity.
        if isinstance(node,ast.IfExp):node=node.body
        actual[target]=asts.evaluate(node,env)
    fn=asts.method('pulse_end_stress_C3','pulse_coefficients');ns=dict(math=math,mp=SimpleNamespace(mpf=s.Rational),
        axial_derivative=lambda v:s.diff(v,z),product_rows=product_rows,shifted_rows=shifted_rows)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual full pulse coefficients>','exec'),ns)
    zeros=[s.Integer(0)]*4
    rows=ns['pulse_coefficients'](delta,mu,z,C,Xp,[bh,s.Symbol('Bhat_s')]+zeros[:3],
        [m]+zeros,[n]+zeros,[ev]+zeros,[J]+zeros,[P]+zeros)
    def value(part):
        rp,bp,dp,hp=map(s.Rational,part['mode'])
        return R**rp*B**bp*D**dp*H**hp/s.sqrt(2)*part['shape'][0]
    zero('actual_full_theta_3_16_3_18_all_meridional_terms',
        sum(value(v) for v in rows['theta'].values())-actual['Itheta']-actual['Stheta'])
    zero('actual_full_axial_3_17_3_18_all_meridional_terms',
        sum(value(v) for v in rows['axial'].values())-actual['Iz']-actual['Sz'])
    for label,parts in rows.items():
        for name,part in parts.items():
            rp,bp,dp,hp=map(s.Rational,part['mode'])
            zero(label+'_'+name+'_actual_positive_factor_logR_rate',
                part['source_logR_rate']-(rp-(s.Rational(1,2)+mu)*bp-r*hp))
    return dict(identities=proofs,input_hashes=asts.hashes,
        actual_full_paper_moment_formula_AST_replayed=True,all_meridional_cross_terms_retained=True,
        original_equations=['(3.16)','(3.17)','(3.18)'],paper_pdf_sha256=PAPER_SHA256,
        source_caps_used_as_defining_field_values=False)


def capture_native_end(field,Z,offset,cells):
    """Replay the unchanged native callable, exposing two existing locals."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('axial_pulse_field','end'))
    returns=[node for node in ast.walk(fn) if isinstance(node,ast.Return)]
    if len(returns)!=1 or not isinstance(returns[0].value,ast.Call):raise ValueError('Native end packet changed')
    for key,expression in (
        ('formal_backward_energy_loss_Taylor','end_energy*c.exp(2*self.mu*s)'),
        ('formal_unperturbed_energy_Taylor','future*c.exp(2*self.mu*s)+baseline')):
        returns[0].value.keywords.append(ast.keyword(arg=key,value=ast.parse(expression,mode='eval').body))
    ast.fix_missing_locations(fn)
    env=dict(mp=mp,endpoints=endpoints,raw_beta=raw_beta,backward_bump_weights=backward_bump_weights)
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<unchanged native end with local capture>','exec'),env)
    return env['end'](field,Z,offset,cells),asts.hashes


def pulse_source_bridge(flatten):
    asts=SourceAST();checks={}
    def syntax(stem,name,target,wanted=None,augmented=False):
        return asts.expression(stem,name,target,wanted=wanted,augmented=augmented)
    syntax('axial_pulse_field','end','C',"selected['selected_scaled_end_coefficient_Taylor']")
    syntax('axial_pulse_field','end','Bhat','Cj*beta',True)
    syntax('axial_pulse_field','end','end_energy','Cj*Cj*(c.exp(-2*self.mu*center)*w[2])',True)
    syntax('axial_pulse_field','end','future',"self.selection.future.future(Z)['complete_future_energy_Taylor']/2")
    syntax('axial_pulse_field','end','X','1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)')
    syntax('pulse_mixed_C4','_high_packet','beta',"self.flat.beta(coord['offset_from_Rv']-center)")
    syntax('pulse_mixed_C4','_high_packet','Brows[k]','Cj*(beta[k]*math.factorial(k))',True)
    checks['actual_selected_C5_coefficients_and_original_beta_C4_source_retained']=True
    checks['actual_backward_signed_linear_and_quadratic_energy_source_retained']=True
    checks['actual_original_signed_scalar_X_memory_retained']=True
    for flag in ('actual_flatten_future_energy_half_of_complete_C5_source_AST_verified',
        'same_absolute_pressure_identified_by_original_FTC_and_power_datum',
        'actual_complete_future_contains_selected_outer_angular_repairs_and_entire_heat'):
        if not flatten.bridge['identities'][flag]:raise ValueError('Same absolute/future source unavailable: '+flag)
        checks['consumed_'+flag]=True
    # Actual native pressure endpoint equality, with exponential caps replaced
    # ONLY for symbolic replay by the exact defining exponential.
    mu,z,U,Pin,datum=s.symbols('mu Z U Pin analytic_P0',real=True);q=1+z*z
    c=SimpleNamespace(exp=s.exp)
    node=syntax('flatten_mixed_C4','flatten','Mp_v',
        "data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))")
    expected=Pin/q**2+datum+(U/q)**2*(1-s.exp(-13/mu-26))/(2*(1+2*mu))
    flat=asts.evaluate(node,dict(data=dict(Mp=Pin/q**2,u=U/q),
        self=SimpleNamespace(pressure_decay=s.exp(-13/mu-26),prate=1+2*mu)))+datum
    if s.simplify(flat-expected)!=0:raise ArithmeticError('Native flatten endpoint pressure changed')
    syntax('pulse_radial_C4','pressure_moment','p',"IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)")
    syntax('pulse_radial_C4','pressure_moment','log_decay','sum(log_parts.values(),c.mpf(0))')
    syntax('pulse_radial_C4','pressure_moment','kernel')
    checks['actual_native_pulse_flatten_absolute_pressure_endpoint_formula_bound']=True
    checks['canonical_same_pressure_transported_inward_by_original_FTC_not_fitted']=True
    return dict(identities=checks,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        source_caps_used_as_defining_field_values=False,whole_original_end_interval_retained=True,
        interval_overlap_used_as_functional_join_proof=False)

def pulse_reference_amplitude_bridge(flatten):
    """Original source recipes, not overlap of numerical scale intervals."""
    asts=SourceAST();checks={}
    a,mu,z,L,Ts,wait,eps=s.symbols('a mu Z Lrel Ts wait epsilon',real=True)
    r=1-mu;k=1-a;bh=s.Rational(1,2)+a;bp=s.Rational(1,2)+mu
    c=SimpleNamespace(mpf=s.Rational,exp=s.exp,expm1=lambda x:s.exp(x)-1,ln=s.log)
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,k=k,delta=2*a,bh=bh)
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,variable=lambda ctx,value,order:value)
    env=dict(math=math,IntervalTaylor=stub,product_rows=product_rows,
        quotient_log_rates=lambda rows:[s.Integer(0)]*4,
        sigma_jets=lambda ctx,t:[s.Integer(0)]*5,
        log_taylor=s.log,flatten_f=lambda ctx,Z:(1+Z*Z)/2)
    for stem,name in (('steep_entry_stress_C3','entry_shape'),('angular_stress_C3','angular_shape'),
        ('outer_power_stress_C3','outer_power_shape'),('flatten_stress_C3','flatten_shape')):
        asts.replay(stem,name,env)
    for flag in ('actual_KS_is_original_power_phase0_not_exit_KQ',
        'actual_native_theta_equals_same_reference_B_times_K_verified'):
        if not flatten.power.angular.entry.bridge['identities'][flag]:raise ValueError('Actual entry amplitude binding missing: '+flag)
        checks['consumed_'+flag]=True
    KS=(1-eps)*s.exp(k*(Ts+s.Rational(1,2)))
    entry=env['entry_shape'](heat,dict(KS=KS),s.Integer(0),dict(f=s.Rational(1,2)))
    angular=env['angular_shape'](heat,dict(KR=entry['K_rows'][0]),s.Integer(-4),
        dict(swirl_factor_one_plus_h_Taylor=s.Integer(1),actual_angular_bump_y_derivatives=[0]*5))
    power=env['outer_power_shape'](heat,dict(Kright=angular['K_rows'][0]),-(L-4))
    shape=env['flatten_shape'](heat,dict(Kright=power['K_rows'][0]),s.Integer(-100),
        dict(Z=z,sigma_y_derivatives=[0]*5,F_Taylor=s.Integer(1),actual_original_flatten_right_endpoint=False))
    C0=2*power['K_rows'][0]*s.exp(-100*(a-mu))
    if s.simplify(shape['K_rows'][0]-C0/(1+z*z))!=0:raise ArithmeticError('Original flatten C0/q source changed')
    checks['actual_entry_angular_power_flatten_K0_is_scalar_C0_over_q_AST_replayed']=True
    obj=SimpleNamespace(ctx=c,bp=bp,bh=bh,rate=r,k=k,mu=mu,Ts=Ts,wait=wait,
        logone=s.log(1-eps),outer=SimpleNamespace(Lrel=L))
    for name in ('thetaR','thetaS','thetaQ','thetaT'):
        setattr(obj,name,asts.evaluate(asts.expression('steep_waiting_C4','__init__','self.'+name),dict(self=obj,c=c)))
    h=SimpleNamespace(ctx=c,bh=bh,steep=obj)
    base=asts.evaluate(asts.expression('collar_Gamma_C4','__init__','self.theta_base'),dict(self=h,c=c))
    q0=-wait-Ts-102-L
    if s.simplify(s.expand_power_exp(C0*base*s.exp(-bh*q0)-1))!=0:
        raise ArithmeticError('Original pulse/flatten reference amplitude differs')
    checks['actual_original_C0_times_theta_base_exp_minus_bh_q0_equals1']=True
    if s.simplify(s.expand_power_exp((C0*base*s.exp(-bh*q0))**2-1))!=0:
        raise ArithmeticError('Original absolute pressure scale differs')
    checks['actual_formal_absolute_pressure_scale_ratio_is_same_C0_squared']=True
    flag='actual_flatten_radius_and_right_power_interface_AST_verified'
    if not flatten.bridge['identities'][flag]:raise ValueError('Actual source radius binding missing: '+flag)
    checks['consumed_'+flag]=True
    name=PREFIX+'collar_stress_C3_check.json';receipt=json.loads((HERE/name).read_bytes())
    if not receipt['all_passed'] or (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(flatten.family,flatten.source):
        raise ValueError('Actual collar pressure unit source differs')
    for source,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Actual pressure unit source changed: '+source)
    asts.hashes.update(receipt['input_hashes']);asts.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    for flag in ('actual_Ev0_squared_over_Pstar_squared_is_source_exactlog_not_cap',
        'actual_absolute_pressure_units_use_exact_Ev0_squared','consumed_production_C4_pressure_scale'):
        if not receipt['actual_collar_stress_source_bridge']['identities'][flag]:raise ValueError('Actual pressure unit bridge missing: '+flag)
        checks['consumed_'+flag]=True
    ps,u=s.symbols('Pstar actual_inlet_U',positive=True)
    native=SimpleNamespace(mu=mu,inlet=SimpleNamespace(constants=dict(U=u)))
    native.U=asts.evaluate(asts.expression('flatten_mixed_C4','__init__','self.U'),dict(self=native))
    native.logEv2_parts=asts.evaluate(asts.expression('flatten_mixed_C4','__init__','self.logEv2_parts'),dict(self=native,c=c))
    exactlog=asts.evaluate(asts.expression('flatten_mixed_C4','__init__','exactlog'),dict(self=native,c=c))
    ev=asts.evaluate(asts.expression('collar_stress_C3','collar_stress_source_bridge','ev'),dict(s=s,ps=ps,u=native.U,mu=mu))
    physical=dict(
        actual_Ev0_AST_equals_Pstar_times_inlet_U_times_original_exponent=ev-ps*native.U*s.exp(-13/(2*mu)-13),
        actual_Ev0_squared_equals_production_exactlog_not_cap=ev**2/ps**2-s.exp(exactlog),
        actual_pulse_B0_equals_C0_times_same_physical_flatten_reference=ev-C0*ev*base*s.exp(-bh*q0),
        actual_absolute_pressure_scale_equals_same_physical_Ev0_squared=ev**2-C0**2*(ev*base*s.exp(-bh*q0))**2)
    for name,value in physical.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual pulse physical source units differ: '+name)
        checks[name]=True
    return dict(identities=checks,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        original_positive_parameter_range_required=True,
        pulse_flatten_reference_amplitude_AST_bridge_verified=True,
        source_caps_used_as_defining_field_values=False)


class CompliantPulseEndStressC3:
    @source_precision
    def __init__(self,cells=64):
        self.flatten=CompliantFlattenStressC3();self.native=CompliantPulseMixedC4()
        self.assembly=CompliantGlobalPhysicalAssembly();self.ctx=self.native.ctx;self.cells=cells
        self.heat=self.flatten.heat;self.family,self.source=self.flatten.family,self.flatten.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):
            raise ValueError('Pulse end physical source family differs')
        self.hashes=dict(self.flatten.hashes);self.hashes.update(self.native.hashes);self.hashes.update(self.assembly.hashes)
        for stem,flag in (('pulse_mixed_C4_check','pulse_all_mixed_derivatives_total_order_le4_available'),
            ('flatten_stress_C3_check','actual_original_flatten_similarity_stress_recovered')):
            name=PREFIX+stem+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[flag]:raise ValueError('Pulse stress prerequisite missing: '+flag)
            for source,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Pulse stress source changed: '+source)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=original_formula_identities();self.bridge=pulse_source_bridge(self.flatten)
        self.unit_bridge=pulse_reference_amplitude_bridge(self.flatten)
        self.hashes.update(self.unit_bridge["input_hashes"])
        self.hashes.update(self.proof['input_hashes']);self.hashes.update(self.bridge['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.cache={}

    @source_precision
    def end(self,Z,offset):
        c=self.ctx;z0=c.mpf(Z);v=c.mpf(offset)
        if endpoints(z0)[0]<-1 or endpoints(z0)[1]>1 or endpoints(v)[0]<-4 or endpoints(v)[1]>0:
            raise ValueError('Original end domain Z[-1,1],s[-4,0] required')
        native,hashes=capture_native_end(self.native,z0,v,self.cells);self.hashes.update(hashes)
        selected=self.native.data(z0)[0];z=IntervalTaylor.variable(c,z0,5)
        q=IntervalTaylor(c,[1+z0**2,2*z0,1,0,0,0]);C=q.reciprocal()
        Bh=[C*0 for _ in range(5)]
        for Cj,center in zip(selected['selected_scaled_end_coefficient_Taylor'],(-3,-1)):
            beta=self.native.flat.beta(v-center)
            for j in range(5):Bh[j]+=Cj*(beta[j]*math.factorial(j))
        m,n=([native['formal_scaled_Mz_mixed_moments'][i]] for i in (0,1))
        e0=[native['formal_unperturbed_energy_Taylor']];J=[native['formal_backward_energy_loss_Taylor']]
        mu=self.native.mu;prate=1+2*mu
        key=tuple(endpoints(z0))
        if key not in self.cache:
            packet=self.flatten.flatten(z0,0)
            scale=2*c.mpf(endpoints(self.flatten.terminal(z0)['Kright']))*c.exp(-100*(c.mpf(endpoints(self.heat.a))-mu))
            # K0(Z)=C0/(1+Z^2); scale is the SAME scalar reference ratio.
            # The exact AST unit bridge binds this scalar before enclosing
            # the same absolute pressure; no scale is fitted.
            rawP0=packet['flatten_future_defect_zeroth_Taylor']['pressure_defect_rows']
            rawP=IntervalTaylor(c,[c.mpf(endpoints(v)) for v in rawP0.coefficients])
            P0=-(IntervalTaylor.constant(c,1,5)/(2*self.heat.prate)+rawP)/(scale*scale)
            self.cache[key]=(P0,scale)
        P0,scale=self.cache[key]
        P=[P0*c.exp(prate*v)+C*C*(c.expm1(prate*v)/(2*prate))]
        for j in range(4):
            m.append(Bh[j]-m[j]*(c.mpf('.5')-mu));n.append(Bh[j]-n[j]*(c.mpf('.5')-2*mu))
            square=sum((Bh[l]*Bh[j-l]*math.comb(j,l) for l in range(j+1)),C*0)
            e0.append(e0[j]*(2*mu)-(C*0+1/2 if j==0 else C*0))
            J.append(J[j]*(2*mu)-square)
            P.append(P[j]*prate+(C*C/2 if j==0 else C*0))
        parts=pulse_coefficients(self.native.delta,mu,z,C,self.native.Xp,Bh,m,n,e0,J,P)
        logB=dict(logPstar=self.assembly.logP,actual_log_inlet_U=c.ln(self.flatten.native.U),
            inverse_mu=-13/(2*mu),finite=-13-(c.mpf('.5')+mu)*v)
        logR=self.assembly.logRp+13/mu+v
        logD=self.native.logE;logH=-13*(1-mu)/mu-(1-mu)*v
        packet={}
        for label,rows in parts.items():
            packet[label]={}
            for name,part in rows.items():
                rp,bp,dp,hp=part['mode']
                logs=dict(source_logR=rp*logR,**{key:bp*value for key,value in logB.items()},
                    selected_log_end_scale=dp*logD,signed_original_memory_log=hp*logH,normalization=-c.ln(2)/2)
                grid={'s'+str(j)+'_Z'+str(l):jet[l]*math.factorial(l)
                    for j,jet in enumerate(part['full_derivative_rows']) for l in range(4-j)}
                packet[label][name]=dict(mode=part['mode'],exact_source_log_parts=logs,
                    exact_source_log_factor=sum(logs.values(),c.mpf(0)),source_logR_rate=part['source_logR_rate'],
                    full_stress_mixed3_coefficient_enclosures=grid)
        return dict(Z=z0,s=v,domain=DOMAIN,full_meridional_stress_log_sectors=packet,
            source_formal_Bhat_rows=Bh,source_formal_Mz_rows=m,source_formal_Mtheta_z_rows=n,
            source_unperturbed_energy_rows=e0,source_selected_backward_energy_loss_rows=J,
            signed_pressure_relative_pure_swirl_reference_rows=P,
            native_original_absolute_forward_pressure=native['pressure'],
            original_X_definition='1/(1-mu)+(Xp-1/(1-mu))*exp(-13*(1-mu)/mu-(1-mu)*s)',
            exact_pulse_reference_logB_parts=logB,exact_logR=logR,exact_logD=logD,exact_logH=logH,
            required_source_amplitude_ratio_enclosure=scale,original_source_fields_and_caps_retained=True,
            all_full_meridional_stress_terms_recovered=True,source_caps_used_as_defining_field_values=False,
            pulse_flatten_reference_amplitude_AST_bridge_verified=True,pulse_flatten_full_stress_join_verified=False,
            pulse_end_physical_decomposition_constructed=False,pulse_end_cone_certified=False,
            global_admissible_stress_lift_constructed=False,physical_energy_integral_certified=False,temporal_recursion=False)

    @source_precision
    def report(self):
        keys=('Z','s','full_meridional_stress_log_sectors','source_formal_Bhat_rows','source_formal_Mz_rows',
            'source_formal_Mtheta_z_rows','source_unperturbed_energy_rows','source_selected_backward_energy_loss_rows',
            'signed_pressure_relative_pure_swirl_reference_rows','required_source_amplitude_ratio_enclosure')
        packets=[self.end(z,v) for z,v in (('0','-3'),('.5','-1'),('.5','-2'),('0','0'),([-1,1],[-4,0]))]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            source_bound_full_formula_proof=self.proof,source_bridge=self.bridge,reference_amplitude_bridge=self.unit_bridge,
            samples=[{key:p[key] for key in keys} for p in packets[:-1]],
            whole_original_end={key:packets[-1][key] for key in keys},input_hashes=self.hashes,
            all_full_meridional_stress_terms_recovered=True,whole_original_end_domain_retained=True,
            source_caps_used_as_defining_field_values=False,pulse_flatten_reference_amplitude_AST_bridge_verified=True,
            pulse_flatten_full_stress_join_verified=False,pulse_end_physical_decomposition_constructed=False,
            pulse_end_cone_certified=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False)


@source_precision
def run():
    result=CompliantPulseEndStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original full-meridional pulse end stress in exact log sectors generated; functional join/physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
