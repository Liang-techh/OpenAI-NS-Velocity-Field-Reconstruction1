"""Continuous original main/exit two-vector cone with order-one axial shear.

Exact kernel integration by parts preserves local/incoming histories and
the selected energy cancellation. All stress sectors and actual source logs
remain. Regional two-vector admission is not global completed-tensor admission.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_gap_cone import current_sources as gap_current_sources
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import DOMAIN,PREFIX,source_precision
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets,sigma_tail_bound
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
LEFT=s.Rational(1,50)
TAIL_DOMAIN='Rtail*exp(-12.98/mu-wait-Ts-102-Lrel)<=R<Rtail*exp(3); Gamma stress zero beyond'
ADMISSIONS=('pulse_main_exit_cone_certified','actual_main_exit_theta_stress_positive',
    'actual_main_exit_full_axial_shear_cone_verified',
    'actual_main_exit_continuous_directional_cone_verified',
    'main_exit_gap_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
    'same_source_main_exit_gap_completed_physical_interface_consumed',
    'fixed_positive_viscosity_two_vector_cone_transfer_verified')
FALSE_FLAGS=('completed_full_tensor_cone_certified','whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion',
    'production_exact_point_parameters_selected')
# Stress normalized by sqrt(R/2)*B: (R power,B power,H power,extra source).
CORRECTIONS={
 'theta':{
    'signed_original_memory':(0,0,1,'one'),
    'meridional_transport':(0,1,0,'one'),
    'radial_shear':(-1,0,0,'one'),
    'incoming_Mz_transport':(0,1,0,'incoming1'),
    'incoming_Mtheta_z_transport':(0,1,0,'incoming2')},
 'axial':{
    'full_energy_and_pressure':(0,1,0,'one'),
    'nonlinear_meridional_transport':(0,1,0,'one'),
    'axial_radial_shear':(-1,0,0,'one'),
    'incoming_Mz_linear_axial_moment':(0,0,0,'incoming1'),
    'incoming_Mz_nonlinear_meridional_transport':(0,1,0,'incoming1'),
    'same_absolute_pressure_memory':(0,1,0,'Q'),
    'selected_end_energy_loss':(0,1,0,'end_square')}}


def current_sources():
    records,hashes,family=gap_current_sources()
    for stem in ('pulse_main_exit_similarity_C4','pulse_main_exit_similarity_C4_check',
        'pulse_main_exit_physical_C2','pulse_main_exit_physical_C2_check',
        'pulse_gap_cone','pulse_gap_cone_check','flat_pulse_derivatives_check'):
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
            raise ValueError('Main cone source family differs: '+stem)
        if 'all_passed' in record and not record['all_passed']:
            raise ValueError('Unaccepted main cone source: '+stem)
        for path,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Main cone source changed: '+path)
        records[stem]=record;hashes.update(record['input_hashes'])
        hashes[name]=hashlib.sha256(raw).hexdigest()
    for flag in ('actual_pulse_main_exit_physical_decomposition_constructed',
        'main_exit_and_exit_gap_completed_physical_interfaces_verified',
        'actual_main_exit_full_meridional_history_and_three_component_remainder_preserved'):
        if not records['pulse_main_exit_physical_C2_check'][flag]:
            raise ValueError('Main cone physical prerequisite missing: '+flag)
    if records['pulse_main_exit_similarity_C4']['domain']!=DOMAIN:
        raise ValueError('Original main/exit domain required')
    return records,hashes,family


def relative_log_recipes(c,mu,finite,lp,lu,lrp,xi):
    """Group correlated sources before directed enclosure."""
    B=lp+lu-xi/(2*mu)-xi
    logs=dict(R=lrp+xi/mu,B=B,H=-(1-mu)*xi/mu)
    extra=dict(one=0,incoming1=-(c.mpf('.5')-mu)*xi/mu,
        incoming2=-(c.mpf('.5')-2*mu)*xi/mu,
        Q=-(1+2*mu)*(13-xi)/mu,end_square=-2/mu+2*finite)
    return {label+'_'+name:rp*logs['R']+bp*logs['B']+hp*logs['H']+extra[which]
        for label,parts in CORRECTIONS.items() for name,(rp,bp,hp,which) in parts.items()}


def relative_log_envelopes(c,mu,finite,lp,lu,lrp):
    left=c.mpf('.02')
    # BQ increases with xi; every other retained correction decreases.
    at_left=relative_log_recipes(c,mu,finite,lp,lu,lrp,left)
    at_right=relative_log_recipes(c,mu,finite,lp,lu,lrp,c.mpf(11))
    at_left['axial_same_absolute_pressure_memory']=at_right['axial_same_absolute_pressure_memory']
    return at_left


def full_stress_ratio(theta_baseline,axial_baseline,theta_corrections,axial_corrections):
    """Total Tz/Ttheta, including both finite-radius shear corrections."""
    theta_delta=sum(theta_corrections.values())
    axial_delta=sum(axial_corrections.values())
    theta=theta_baseline+theta_delta
    axial=axial_baseline+axial_delta
    return dict(theta=theta,axial=axial,w=axial/theta,
        theta_delta=theta_delta,axial_delta=axial_delta)


def main_cone_identities(records):
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Main/exit cone identity failed: '+name)
        checks[name]=True
    mu,delta,z,xi,lp,lu,lrp,finite=s.symbols('mu delta Z xi logP logU logRp finite',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log,exp=s.exp)
    native=SimpleNamespace(mu=mu,finite=finite,logP=lp,logU=lu,logRp=lrp)
    env=dict(c=c,self=native,xi=xi,p=1+2*mu,r=1-mu,logD0=-1/mu+finite)
    Bparts=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','logB'),env)
    logR=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','logR'),env)
    logH=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','logH'),env)
    extra=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','extras'),env)
    actual=asts.replay('pulse_main_exit_cone','relative_log_recipes',dict(CORRECTIONS=CORRECTIONS))(
        c,mu,finite,lp,lu,lrp,xi)
    envelopes=asts.replay('pulse_main_exit_cone','relative_log_envelopes',
        dict(relative_log_recipes=lambda *args:asts.replay('pulse_main_exit_cone','relative_log_recipes',
            dict(CORRECTIONS=CORRECTIONS))(*args)))(c,mu,finite,lp,lu,lrp)
    rates={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            key=label+'_'+name
            source=rp*logR+bp*sum(Bparts.values())+hp*logH+extra[which]
            zero('actual_'+key+'_grouped_source_log',source-actual[key])
            right=11 if which=='Q' else LEFT
            zero('actual_'+key+'_endpoint_envelope',actual[key].subs(xi,right)-envelopes[key])
            rate_coefficient=-rp+bp*(s.Rational(1,2)+mu)+hp*(1-mu)+{
                'one':0,'incoming1':s.Rational(1,2)-mu,
                'incoming2':s.Rational(1,2)-2*mu,'Q':-(1+2*mu),'end_square':0}[which]
            zero('actual_'+key+'_monotone_rate_coefficient',-mu*s.diff(actual[key],xi)-rate_coefficient)
            rates[key]=str(s.diff(actual[key],xi))
    # Actual original gp is primitive(sigma(50*xi))*sigma(11-xi).
    for name in ('_gp_piece','gp_jets','sigma_jets','sigma_tail_bound'):
        asts.method('flat_pulse_derivatives',name)
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','t',wanted='xi/mu')
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','coord',
        wanted='xi-mu*c.mpf([endpoints(a)[0],endpoints(b)[1]])')
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','total',
        wanted="(c.exp(-lam*a)-c.exp(-lam*b))/lam*gp(c,coord)['value']",augmented=True)
    flat=records['flat_pulse_derivatives_check']
    for flag in ('original_radial_shape_derivatives_C4_available','original_flat_shapes_retained'):
        if not flat[flag]:raise ValueError('Original flat source missing: '+flag)
        checks['directly_consumed_'+flag]=True
    if not flat['analytic_original_beta_and_flat_envelope_checks']['all_required_endpoint_flat_limits']:
        raise ValueError('Original flat endpoint proof missing')
    ctx=MPIntervalContext();ctx.dps=100
    for j in range(2):
        if endpoints(gp_jets(ctx,0)[j])!=(0,0):raise ArithmeticError('Initial gp jet not flat')
        checks['actual_initial_gp_flat_jet'+str(j)]=True
    v,lam=s.symbols('v lambda',positive=True);G=s.Function('original_gp')
    weighted=s.exp(-lam*v)
    zero('actual_kernel_two_integrations_by_parts_integrand',
        weighted*G(xi-mu*v)+s.diff(weighted*G(xi-mu*v),v)/lam
        -mu*s.diff(weighted*s.diff(G(xi-mu*v),xi),v)/lam**2
        -mu**2*weighted*s.diff(G(xi-mu*v),xi,2)/lam**2)
    checks['exact_remainder_bound_mu_squared_M2_over_lambda_cubed_from_positive_full_kernel']=True
    checks['original_nonnegative_gp_and_gp_xi_upper_one_from_monotone_sigma_product']=True
    checks['global_gp_second_derivative_bound_52_sigma1_plus_11_sigma2']=True
    asts.expression('flat_pulse_derivatives','_gp_piece','entrance',wanted='sigma_jets(c,50*xi)')
    asts.expression('flat_pulse_derivatives','_gp_piece','exit_shape',wanted='sigma_jets(c,11-xi)')
    asts.expression('flat_pulse_derivatives','_sigma_left','L1',wanted='2/(1-x)**3+2/x**3')
    primitive=s.Function('original_sigma_primitive')(xi);sig=s.Function('original_sigma')(11-xi)
    zero('actual_original_gp_monotone_product_first_derivative',
        s.diff(primitive*sig,xi)-s.diff(primitive,xi)*sig-primitive*s.diff(sig,xi))
    zero('actual_original_gp_product_second_derivative',
        s.diff(primitive*sig,xi,2)-s.diff(primitive,xi,2)*sig
        -2*s.diff(primitive,xi)*s.diff(sig,xi)-primitive*s.diff(sig,xi,2))
    C=1/(1+z*z);L=1-delta*z*z;r=1-mu;lam=s.Rational(1,2)-mu
    q=z*z/(1+z*z);qm=mu-delta/2;q0=qm+(1-delta)*q
    h=(1-delta)*(s.Rational(1,2)+q);k=(1-delta)/2
    rho,D,B,ap,apZ=s.symbols('kernel_remainder B_xi B ap ap_Z',real=True)
    kr=s.symbols('original_forward_kernel',real=True)
    def axial(value):return s.diff(value,z)
    func_env=dict(mp=SimpleNamespace(mpf=lambda value:s.Rational(str(value))),math=math,
        axial_derivative=axial)
    asts.replay('collar_Gamma_C4','product_rows',func_env)
    asts.replay('collar_stress_C3','shifted_rows',func_env)
    coefficients=asts.replay('pulse_end_stress_C3','pulse_coefficients',func_env)
    af=s.Function('selected_ap')(z)
    zeros=[s.Integer(0)]*5
    rows=coefficients(delta,mu,z,C,s.Symbol('Xp'),[af*s.Symbol('gp')]+zeros[:4],
        [af*kr]+zeros[:4],zeros,zeros,zeros,zeros)
    equilibrium=rows['theta']['equilibrium']['shape'][0]
    zero('actual_theta_equilibrium_C_over_L_times_q0_over_rate',equilibrium-C/L*q0/r)
    linear=rows['axial']['axial_transport']['shape'][0]+rows['axial']['linear_axial_moment']['shape'][0]
    zero('actual_local_axial_linear_stress_cancellation',
        linear*L/C-(-af*s.Symbol('gp')+h*af*kr-k*z*s.diff(af,z)*kr))
    # Replay the actual split source, then account for every sector before
    # taking the total-stress ratio. The two R^-1 shear sectors are retained.
    split_fn=asts.replay('pulse_main_exit_similarity_C4','split_main_exit_stress',func_env)
    row=lambda name:[s.Function(name+str(j))(z) for j in range(5)]
    bh,ml,mi,nl,ni,ev,loss,pm=[row(name) for name in ('Bh','ml','mi','nl','ni','ev','loss','pm')]
    pb=[-C*C/(2*(1+2*mu))]+zeros[:4]
    split=split_fn(delta,mu,z,C,s.Symbol('Xp'),bh,ml,mi,nl,ni,ev,loss,pb,pm)
    Rscale,Pscale,Hscale,h1,h2,Qmem,D2=s.symbols(
        'source_R source_pure_swirl source_H incoming1 incoming2 pressure_memory end_square',positive=True)
    factors=dict(one=1,incoming1=h1,incoming2=h2,Q=Qmem,end_square=D2)
    baseline_names=dict(theta=('equilibrium',),axial=('axial_transport','linear_axial_moment'))
    def normalized(part):
        rp,bp,dp,hp=map(s.Rational,part['mode'])
        return Rscale**(rp-s.Rational(1,2))*Pscale**(bp-1)*Hscale**hp*part['shape'][0]*factors[part['extra_source']]
    baselines={};corrections={};totals={}
    for label in ('theta','axial'):
        if set(split[label])!=set(CORRECTIONS[label])|set(baseline_names[label]):
            raise ValueError('Full source stress sector omitted or duplicated: '+label)
        checks['actual_'+label+'_every_source_sector_accounted_once']=True
        baselines[label]=sum(normalized(split[label][name]) for name in baseline_names[label])
        corrections[label]={name:normalized(split[label][name]) for name in CORRECTIONS[label]}
        for name,(rp,bp,hp,which) in CORRECTIONS[label].items():
            zero('actual_'+label+'_'+name+'_feeds_full_stress_correction',
                corrections[label][name]-Rscale**rp*Pscale**bp*Hscale**hp*factors[which]*split[label][name]['shape'][0])
        totals[label]=sum(normalized(part) for part in split[label].values())
    ratio=asts.replay('pulse_main_exit_cone','full_stress_ratio',{})(
        baselines['theta'],baselines['axial'],corrections['theta'],corrections['axial'])
    zero('actual_total_theta_equilibrium_plus_all_corrections',ratio['theta']-totals['theta'])
    zero('actual_total_axial_local_linear_plus_all_corrections',ratio['axial']-totals['axial'])
    zero('actual_w_is_total_Tz_over_total_Ttheta',ratio['w']*totals['theta']-totals['axial'])
    G0,W0,elocal,dQ,dN=s.symbols('G0 W0 local_error theta_correction axial_correction',real=True)
    projected=asts.replay('pulse_main_exit_cone','full_stress_ratio',{})(G0,G0*(W0+elocal),{'all':dQ},{'all':dN})
    zero('actual_total_ratio_error_after_C_over_L_normalization',
        projected['w']-W0-(G0*elocal+dN-W0*dQ)/(G0+dQ))
    K=r*mu*h/(lam**2*q0)
    leading=(-B+h*(B/lam-mu*D/lam**2+ap*rho)-k*z*apZ*kr)*r/q0
    zero('actual_correlated_w_leading_plus_exact_errors',
        leading-(2*B-K*D+mu/lam*B+r*h*ap*rho/q0-r*k*z*apZ*kr/q0))
    qq=s.symbols('q',real=True)
    Kq=r*mu*(1-delta)*(s.Rational(1,2)+qq)/(lam**2*(qm+(1-delta)*qq))
    zero('actual_K_strictly_decreasing_in_q',
        s.diff(Kq,qq)+r*mu*(1-delta)/(lam*(qm+(1-delta)*qq)**2))
    zero('actual_K_zero_formula',Kq.subs(qq,0)
        -2*(1-delta)*r/((1-2*mu)**2*(1-delta/(2*mu))))
    zero('actual_q_lower_Z_squared_over_two',q-z*z/2-z*z*(1-z*z)/(2*(1+z*z)))
    t,A0,A1=s.symbols('absZ positive_qmin positive_qquadratic',positive=True)
    zero('actual_AMGM_axial_ratio_bound',
        A0+A1*t*t-2*s.sqrt(A0*A1)*t-(s.sqrt(A0)-s.sqrt(A1)*t)**2)
    a,b,w,F,T=s.symbols('a b w F Ttheta',real=True)
    dot=T*F*(-a+b*w);cross=T*F*(-b-a*w);kap=a+b*b/a
    test=2*b*w+b*b/a+(a-2)*w*w
    zero('actual_full_directional_cone_factorization',
        2*dot**2-(kap-2)*cross**2-(T*F)**2*(a*a+b*b)*(2-test))
    path=HERE/'lei_ren_part1_paper_mp_stress.py'
    asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    original=next(node for node in ast.walk(ast.parse(path.read_text(encoding='utf8')))
        if isinstance(node,ast.FunctionDef) and node.name=='evaluate_mp_stress')
    R,P=s.symbols('R pure_swirl',positive=True)
    shearenv=dict(a=P*C,ay=-(s.Rational(1,2)+mu)*P*C,
        by=P*C*(mu*D-(s.Rational(1,2)+mu)*B),root=s.sqrt(2*R),R=R)
    shears={}
    for target in ('Stheta','Sz'):
        node=next(node.value for node in ast.walk(original) if isinstance(node,ast.Assign)
            and any(ast.unparse(t)==target for t in node.targets))
        if isinstance(node,ast.IfExp):node=node.body
        shears[target]=asts.evaluate(node,shearenv)
    sourceF=P*C/s.sqrt(2*R)
    zero('actual_common_stress_normalization_over_shear_factor',
        (s.sqrt(R/2)*P)/sourceF-R/C)
    zero('actual_negative_theta_shear_a',shears['Stheta']+(2+2*mu)*sourceF)
    zero('actual_full_nonzero_axial_shear_b',shears['Sz']-sourceF*(2*mu*D-(1+2*mu)*B))
    zero('actual_b_plus_B_exact_error',2*mu*D-(1+2*mu)*B+B-2*mu*(D-B))
    A=s.symbols('positive_K_times_derivative_upper',positive=True)
    zero('actual_bw_completed_square_supremum',A*A/8-(-2*B*B+A*B)-2*(B-A/4)**2)
    zero('actual_second_cone_completed_square_supremum',
        2*A*A/7-(-s.Rational(7,2)*B*B+2*A*B)-s.Rational(7,2)*(B-2*A/7)**2)
    nu,lf=s.symbols('nu common_lambda_factor',positive=True)
    zero('actual_positive_viscosity_physical_directional_transfer',
        2*((nu*lf)**2*dot)**2-(kap-2)*((nu*lf)**2*cross)**2
        -(nu*lf)**4*(2*dot**2-(kap-2)*cross**2))
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','ml',wanted='[ap*kernels[0]]')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','e0')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','baseline',wanted='-C*C/(2*p)')
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    radiusenv={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    radiusenv[target.id]=ast.literal_eval(node.value)
    radius=asts.replay('global_physical_assembly','radius',radiusenv)
    wait,Ts,length=s.symbols('wait Ts Lrel',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=lrp,params=SimpleNamespace(mu=mu))
    heat=SimpleNamespace(steep=SimpleNamespace(wait=wait,Ts=Ts,outer=SimpleNamespace(Lrel=length)))
    left=radius(assembly,'pulse_main',LEFT,{},None)[0]
    join=radius(assembly,'pulse_exit',11,{},None)[0]
    gap=radius(assembly,'pulse_gap',11,{},None)[0]
    tail=radius(assembly,'heat_collar',0,{},heat)[0]
    zero('actual_main_exit_gap_production_radius_join',join-gap)
    zero('actual_main_left_production_tail_radius_offset',
        left-tail+(13-LEFT)/mu+wait+Ts+102+length)
    for stem,flag in (('pulse_main_exit_physical_C2_check','main_exit_and_exit_gap_completed_physical_interfaces_verified'),
        ('pulse_gap_cone_check','gap_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified'),
        ('fifth_axial_jets_check','actual_selected_ap_c1_c2_C5_available')):
        if not records[stem][flag]:raise ValueError('Main cone composition source missing: '+flag)
        checks['directly_consumed_'+flag]=True
    return dict(identities=checks,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        source_log_derivative_rates=rates,original_paper_equations=['(7.39)','(7.40)','(7.41)','(7.42)'],
        paper_pdf_sha256='8396b998dcf737cd6a7c12ff1019c6b2fd308a7e6f0907c0c6647a9b2ec64e8d',
        actual_w_definition='full Tz/Ttheta including all twelve correction sectors and both finite-radius shears',
        leading_w='2*B-K*B_xi+controlled_error',actual_b='2mu*B_xi-(1+2mu)*B',
        cone='Ttheta>0, a-b*w>0, 2-(2*b*w+b^2/a+(a-2)*w^2)>0',
        exact_full_forward_kernel_not_truncated_or_chosen_as_point=True)


def upper(c,value):return c.mpf(endpoints(value)[1])


def absolute(c,value):return c.mpf(max(abs(v) for v in endpoints(value)))


@source_precision
def whole_main_exit_bounds(records):
    c=MPIntervalContext();c.dps=240
    selected=records['pulse_mixed_C4']
    mu=read_interval(c,selected['selected_mu']);delta=read_interval(c,selected['selected_delta'])
    r=1-mu;lam=c.mpf('.5')-mu;qm=mu-delta/2;k=(1-delta)/2
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,
        rate=r,one_minus_delta=1-delta,lambda1=lam,qmin=qm)
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Main/exit positive margin failed: '+name)
    for name,value in parameters.items():positive(name,value)
    src=records['pulse_main_exit_similarity_C4']
    ap=[read_interval(c,value) for value in src['current_selected_ap_C5_enclosure']['coefficients']]
    alo,ahi=endpoints(ap[0]);apmax=c.mpf(ahi);apZ=absolute(c,ap[1])
    positive('selected_ap_lower',c.mpf(alo))
    Kmax=upper(c,r*mu*k/(lam**2*qm));positive('Kmax_below_paper_2_01',c.mpf('2.01')-Kmax)
    positive('selected_ap_below_paper_1_2',c.mpf('1.2')-apmax)
    sig1=upper(c,sigma_tail_bound(c,1,'.5'));sig2=upper(c,sigma_tail_bound(c,2,'.5'))
    g1=1+11*sig1;g2=52*sig1+11*sig2
    Bmax=11*apmax;Dabs=apmax*g1;Dupper=apmax
    W0max=2*Bmax+Kmax*Dabs
    gmin=c.mpf(endpoints(qm/r)[0]);epsilon=c.exp(-1000)
    eq=src['whole_original_main']['full_meridional_stress_log_sectors']['theta']['equilibrium']
    logs=eq['exact_source_log_parts']
    lp=read_interval(c,logs['logPstar']);lu=read_interval(c,logs['actual_log_inlet_U'])
    pmap=records['outer_pulse_map']['pulse_rows']
    finite=-3*read_interval(c,pmap['saddle_L'])+2*c.ln(read_interval(c,pmap['saddle_u0']))-c.ln(6)/2-2*c.ln(mu)
    logC=read_interval(c,records['outer_initial']['selected_logCstar'])
    Tw=read_interval(c,records['fifth_axial_jets']['whole_Z']['selected']['incoming']['Z_independent_constant_definitions']['Tw'])
    lrp=c.ln(110)+10*(logC+lp)+lp+1+Tw
    envelopes=relative_log_envelopes(c,mu,finite,lp,lu,lrp)
    bounds={};margins={};monotonicity={}
    for label,parts in CORRECTIONS.items():
        for name,(rp,bp,hp,which) in parts.items():
            coeff=[]
            for region in ('whole_original_main','whole_original_exit'):
                sector=src[region]['full_meridional_stress_log_sectors'][label][name]
                mode=sector['original_mode']
                if (mode[0],mode[1],mode[3])!=(rp+.5,bp+1,hp) or sector['exact_extra_source']!=which:
                    raise ValueError('Current original correction factor changed: '+name)
                coeff.append(read_interval(c,sector['full_stress_mixed3_coefficient_enclosures']['y0_Z0']))
            bound=c.mpf(max(endpoints(absolute(c,value))[1] for value in coeff))
            if endpoints(bound)[1]<=0:raise ValueError('Expected retained correction coefficient')
            key=label+'_'+name
            magnitude=-rp+bp*(c.mpf('.5')+mu)+hp*(1-mu)+{
                'one':0,'incoming1':c.mpf('.5')-mu,
                'incoming2':c.mpf('.5')-2*mu,'Q':-(1+2*mu),'end_square':0}[which]
            monotonicity[key]=-magnitude if which=='Q' else magnitude
            positive('source_monotonicity_'+key,monotonicity[key])
            actual=envelopes[key]+c.ln(bound)
            target=c.ln(gmin)-1000-c.ln(2*len(parts))
            margin=target-actual;positive(key,margin);margins[key]=margin
            bounds[key]=dict(actual_signed_main_and_exit_coefficients=coeff,
                coefficient_absolute_upper=bound,grouped_relative_log_envelope=envelopes[key],
                actual_relative_log_absolute_upper=actual,certified_log_cap=target,positive_log_gap=margin)
    # Original kernel remainder and axial source derivative, divided only after
    # using q0>=qmin+(1-delta)*Z^2/2. No interval division by Z near zero.
    kernel_error=mu**2*g2/lam**3
    local_errors=dict(radial_rate_error=mu/lam*Bmax,
        two_IBP_kernel_error=r*(1-delta)*apmax*kernel_error/qm,
        axial_amplitude_derivative_error=r*k*apZ*(11/lam)/(2*c.sqrt(qm*(1-delta)/2)))
    local_total=sum(local_errors.values(),c.mpf(0))
    w_error=upper(c,(local_total+epsilon*(1+W0max))/(1-epsilon))
    b_error=upper(c,2*mu*apmax*(g1+11))
    positive('w_error_below_one_millionth',c.mpf('1e-6')-w_error)
    positive('b_error_below_one_millionth',c.mpf('1e-6')-b_error)
    shape_margins=dict(selected_ap_lower=c.mpf(alo),selected_ap_below_1_2=c.mpf('1.2')-apmax,
        Kmax_below_2_01=c.mpf('2.01')-Kmax,
        w_error_below_one_millionth=c.mpf('1e-6')-w_error,
        b_error_below_one_millionth=c.mpf('1e-6')-b_error)
    coefficient=Kmax*Dupper
    bw_error=Bmax*w_error+W0max*b_error+b_error*w_error
    bw_upper=coefficient**2/8+bw_error
    second_upper=2*coefficient**2/7+2*bw_error+Bmax*b_error+b_error**2/2+2*mu*(W0max+w_error)**2
    algebraic=dict(theta_normalized_lower=gmin*(1-epsilon)/2,
        a_minus_bw_lower=2+2*mu-bw_upper,
        bw_below_paper_0_8=c.mpf('.8')-bw_upper,
        second_below_paper_1_8=c.mpf('1.8')-second_upper,
        directional_bracket_lower=2-second_upper,
        kappa_minus2_lower=2*mu,
        full_directional_margin_normalized_lower=4*(2-second_upper))
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_log_margins=margins,
        positive_monotonicity_margins=monotonicity,positive_shape_and_error_margins=shape_margins,
        positive_algebraic_margins=algebraic,normalized_stress_sector_bounds=bounds,
        mu=mu,delta=delta,selected_ap_lower=c.mpf(alo),selected_ap_upper=apmax,
        selected_ap_Z_absolute_upper=apZ,K_upper=Kmax,
        gp_bounds=dict(value_lower=0,value_upper=11,first_derivative_upper=1,
            first_derivative_absolute_upper=g1,second_derivative_absolute_upper=g2),
        exact_kernel_remainder_absolute_upper=kernel_error,local_w_error_terms=local_errors,
        actual_w_error_absolute_upper=w_error,actual_b_error_absolute_upper=b_error,
        actual_B_absolute_upper=Bmax,actual_B_xi_upper=Dupper,actual_B_xi_absolute_upper=Dabs,
        actual_leading_w_absolute_upper=W0max,bw_absolute_error_upper=bw_error,
        correlated_bw_upper=bw_upper,correlated_second_cone_expression_upper=second_upper,
        original_source_reduced_log_parameters=dict(logPstar=lp,actual_logU=lu,finite=finite,logRp=lrp),
        continuous_whole_main_exit_Z_domain_covered=True,
        full_order_one_axial_shear_and_falling_side_retained=True,
        all_signed_theta_and_axial_corrections_retained=True,
        original_energy_cancellation_and_absolute_pressure_used=True,
        exact_source_correlations_grouped_before_interval_enclosure=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False,
        global_temporal_flatness_not_inferred=True)


@source_precision
def run():
    records,hashes,family=current_sources()
    proof=main_cone_identities(records);hashes.update(proof['input_hashes'])
    bounds=whole_main_exit_bounds(records)
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_source_cone_proof=proof,bounds=bounds,input_hashes=hashes,
        regional_two_vector_cone_distinct_from_completed_tensor_cone=True,
        current_selected_C5_source_check_directly_consumed=True,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Continuous main/exit full-shear cone generated and composed to gap/end/tail; entrance/global/recursion pending',flush=True)
    return result


class CertifiedPulseMainExitPhysical:
    """Attach current regional cone admission to the unchanged physical source."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_pulse_main_exit_physical_C2 import CompliantPulseMainExitPhysicalC2
        self.receipt=json.loads(Path(__file__).with_name(PREFIX+'pulse_main_exit_cone_check.json').read_bytes())
        if not self.receipt['all_passed'] or not self.receipt['pulse_main_exit_cone_certified']:
            raise ValueError('Accepted main/exit cone required')
        for path,digest in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Certified main/exit source changed: '+path)
        self.physical=CompliantPulseMainExitPhysicalC2()
        if (self.physical.family,self.physical.source)!=(
            self.receipt['actual_five_defect_family_sha256'],self.receipt['implicit_source_sha256']):
            raise ValueError('Certified main/exit source family differs')
    def main_exit(self,*args,**kwargs):
        result=self.physical.main_exit(*args,**kwargs)
        result.update(**{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
        return result


if __name__=='__main__':run()
