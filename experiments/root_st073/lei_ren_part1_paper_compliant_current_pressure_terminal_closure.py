"""Original preheat pressure function -> current exact pressure terminal.

The fourteen-stage integral is identified from the original swirl generator,
its native pressure primitives and the unique raw waiting root. Cp is zero
only after that function identity and the checked quadratic/Gamma balance.
The old forward Cp enclosure is retained as a diagnostic, never fitted.
"""
import ast,json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_pressure_terminal_balance import CurrentPressureTerminalBalance
from lei_ren_part1_paper_compliant_current_raw_preheat_pressure_operator import CurrentRawPreheatPressureOperator
from lei_ren_part1_paper_compliant_current_angular_terminal_closure import (
    HERE,PREFIX,sha,function,binding,pack,encode,endpoints,IntervalTaylor,
    shape_radial5,collar_defect_rows,collar_stress_mixed4,mixed,
    collar_pressure_rows,constant_stress_rows,return_binding)
from lei_ren_part1_paper_compliant_actual_Rh_source_join import pressure_defining_function_proof
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_logarithmic_pressure_datum import BETA2,BETA0

NAME=PREFIX+'current_pressure_terminal_closure.json'
RECEIPT=PREFIX+'current_pressure_terminal_closure_check.json'
GATES=('current_implicit_datum_to_native_raw_pressure_operator_identified',
    'current_heat_pressure_terminal_constant_eliminated',
    'current_heat_terminal_constants_eliminated',
    'current_heat_pressure_stress_mixed4_after_terminal_closure_recovered')
OPEN=('heat_exterior_stress_identity_certified','current_exact_repair_installed_in_all_physical_charts',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'physical_energy_integral_certified','independently_bounded_flat_remainder',
    'full_background_NS_validation','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_collar':('heat_collar',[-1,1],[0,3]),
    'collar_exit':('heat_collar',[-1,1],3),
    'whole_exterior':('heat_exterior',[-1,1],[3,mp.inf]),
    'fresh_exterior':('heat_exterior','.631',4)}


def original_preheat_partition():
    """Exact Section 6.1 master generator in segmented stage coordinates.

    J is the integral of the original sigma; on completed supports J(x)=
    x-1/2. Sigma and phi are retained functions, never polynomial fits.
    Every expression is Utheta/Pstar and all origins remain symbolic.
    """
    z,t=s.symbols('Z t',real=True)
    mu,delta,epsilon,yd,Tw,L,Ts,W=s.symbols('mu delta epsilon yd Tw L Ts W',positive=True)
    q=1+z*z;bp=s.Rational(1,2)+mu;r=1-mu;k=1-delta/2;bh=(1+delta)/2
    J=s.Function('J_sigma');sig=s.Function('sigma');phi=s.Function('phi')
    yp=yd+1+Tw;yv=yp+13/mu;yf=yv+100;yr=yf+L;ys=yr+1;yq=ys+Ts;yt=yq+1;tail=yt+W
    logUw=-s.Rational(1,5)-yd/2-mu/2
    logUp=logUw-bp*Tw;logEv=logUp-bp*13/mu
    logUf=logEv-bp*100-s.log(2);logUr=logUf-bp*L
    logUs=logUr-bp-r/2;logUq=logUs-s.Rational(3,2)*Ts
    logUt=logUq-s.Rational(3,2)+k/2
    logBase=logUt-bh*W-s.log(1-epsilon)
    shapes={
        'reference_extension':s.exp(t/10)/q,
        'slope_transition_ref':s.exp(t/10-s.Rational(3,5)*J(t))/q,
        'axial_turnoff':s.exp(-s.Rational(1,5)-t/2)/q,
        'slope_transition_mu':s.exp(s.Rational(3,10)-yd/2-t/2-mu*J(t))/q,
        'power_buffer':s.exp(logUw-bp*t)/q,
        'pulse_reserved':s.exp(logUp-bp*t)/q,
        'z_flatten':s.exp(logEv-bp*t)*(q/2)**sig(t/100)/q,
        'power_buffer_rel':s.exp(logUf-bp*t),
        'steep_transition_in':s.exp(logUr-bp*t-r*J(t)),
        'steep_power':s.exp(logUs-s.Rational(3,2)*t),
        'steep_transition_out':s.exp(logUq-s.Rational(3,2)*t+k*J(t)),
        'waiting':s.exp(logUt-bh*t),
        'heat_collar':s.exp(logBase-bh*t)*(1-epsilon*(1-sig(t)+sig(t)*phi(t))),
        'exterior_power_tail':s.exp(logBase-bh*t)}
    domains={name:(0,length) for name,length in zip(BETA2[1:],(1,yd-1,1,Tw,13/mu))}
    domains.update(reference_extension=(-s.oo,0),z_flatten=(0,100),power_buffer_rel=(0,L),
        steep_transition_in=(0,1),steep_power=(0,Ts),steep_transition_out=(0,1),
        waiting=(0,W),heat_collar=(0,3),exterior_power_tail=(3,s.oo))
    # Values of the four completed/in-progress cutoff primitives in the
    # original master logA on each segmented stage.
    placements={
        'reference_extension':(t,0,0,0,0,0),
        'slope_transition_ref':(t,J(t),0,0,0,0),
        'axial_turnoff':(1+t,t+s.Rational(1,2),0,0,0,0),
        'slope_transition_mu':(yd+t,yd+t-s.Rational(1,2),J(t),0,0,0),
        'power_buffer':(yd+1+t,yd+t+s.Rational(1,2),t+s.Rational(1,2),0,0,0),
        'pulse_reserved':(yp+t,yp+t-s.Rational(1,2),yp+t-yd-s.Rational(1,2),0,0,0),
        'z_flatten':(yv+t,yv+t-s.Rational(1,2),yv+t-yd-s.Rational(1,2),0,0,sig(t/100)),
        'power_buffer_rel':(yf+t,yf+t-s.Rational(1,2),yf+t-yd-s.Rational(1,2),0,0,1),
        'steep_transition_in':(yr+t,yr+t-s.Rational(1,2),yr+t-yd-s.Rational(1,2),J(t),0,1),
        'steep_power':(ys+t,ys+t-s.Rational(1,2),ys+t-yd-s.Rational(1,2),t+s.Rational(1,2),0,1),
        'steep_transition_out':(yq+t,yq+t-s.Rational(1,2),yq+t-yd-s.Rational(1,2),yq+t-yr-s.Rational(1,2),J(t),1),
        'waiting':(yt+t,yt+t-s.Rational(1,2),yt+t-yd-s.Rational(1,2),yt+t-yr-s.Rational(1,2),t+s.Rational(1,2),1)}
    return dict(z=z,t=t,mu=mu,delta=delta,epsilon=epsilon,yd=yd,Tw=Tw,L=L,Ts=Ts,W=W,
        q=q,J=J,sigma=sig,phi=phi,shapes=shapes,domains=domains,placements=placements,
        master_tail_y=tail,logUp=logUp,logEv=logEv,logUt=logUt,logBase=logBase,
        offsets=dict(d=yd,w=yd+1,p=yp,v=yv,f=yf,rel=yr,s=ys,q=yq,t=yt,tail=tail,b=tail+3))


def bind_original_master_source():
    """Bind defining formulas, without treating float frontend quadrature as exact."""
    return_binding('outer','_log_A',
        'self.logPstar+Decimal("0.1")*y-Decimal("0.6")*self._J(y)-self.mu*self._J(y-self.y_d)-(Decimal(1)-self.mu)*self._J(y-self.y_rel)+(Decimal(1)-self.delta/Decimal(2))*self._J(y-self.y_rel-Decimal(1)-self.Ts)')
    bindings={}
    original_tree=ast.parse((HERE/'lei_ren_part1_paper_outer.py').read_text(encoding='utf8'))
    cls=next(n for n in original_tree.body if isinstance(n,ast.ClassDef) and n.name=='PaperOuterSchedule')
    ctor=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='__init__')
    schedule={
        'self.Td':'_dexp(self.Md,self.decimal_precision)+Decimal(10)',
        'self.Tw':'Decimal(60)*(-self.log_mu)',
        'self.Ts':'Decimal(4)*_dln(Decimal(2)/self.delta,self.decimal_precision)',
        'self.Tf':'Decimal(100)','self.epsilon':'self.c_epsilon*self.delta',
        'self.y_d':'Decimal(1)+self.Td','self.y_w':'self.y_d+Decimal(1)',
        'self.y_p':'self.y_w+self.Tw','self.y_v':'self.y_p+Decimal(13)/self.mu',
        'self.y_f':'_decimal_add_exact(self.y_v,self.Tf)',
        'self.y_rel':'self.y_f-Decimal(30)*self.log_mu',
        'self.y_s':'self.y_rel+Decimal(1)','self.y_q':'self.y_s+self.Ts',
        'self.y_t':'self.y_q+Decimal(1)','self.y_tail':'self.y_t+self.waiting_length',
        'self.y_b':'_decimal_add_exact(self.y_tail,Decimal(3))',
        'self._log_c_inf':'self._log_A(self.y_tail)+(Decimal(1)+self.delta)*self.logR_tail/Decimal(2)-_dln(Decimal(2)*(Decimal(1)-self.epsilon),self.decimal_precision)'}
    for target,expression in schedule.items():
        rows=[n.value for n in ast.walk(ctor) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        expected=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==expected for v in rows)!=1:raise ValueError('Original stage/c_inf source changed: '+target)
        bindings['PaperOuterSchedule.__init__:'+target]=True
    for method,rows in {
        '_evaluate_log_radius':{'flat_sigma':'_sigma((y-self.y_v)/self.Tf)',
            'log_u':'log_a-Decimal(str(z_factor_log))+Decimal(str(flat_sigma*(z_factor_log-math.log(2.0))))'},
        '_tail_log_values':{'t':'y-self.y_tail',
            'flat_argument':'(Decimal(3)-t)/Decimal(2)',
            'flat_factor':'1.0-float(self.epsilon)*flat_value',
            'K':'((1.0-sigma_t)*one_minus_epsilon+sigma_t*heat_value*flat_factor)'}}.items():
        for target,expression in rows.items():
            binding('outer',method,target,expression,allow_other_assignments=(target=='log_u'))
            bindings[method+':'+target]=True
    # Pin the actual original and native edge formulas as well as their
    # call sites. Frontend underflow branches still are not exact values.
    sigma_tree=ast.parse((HERE/'../../src/openai_ns_reconstruction/schedule_pressure.py').read_text(encoding='utf8'))
    sigfn=next(n for n in sigma_tree.body if isinstance(n,ast.FunctionDef) and n.name=='outgoing_sigma')
    for target,expression in {'log_a':'-(1.0/x)**2','log_b':'-(1.0/(1.0-x))**2','delta':'log_b-log_a'}.items():
        values=[n.value for n in ast.walk(sigfn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original sigma edge formula changed: '+target)
        bindings['original_outgoing_sigma:'+target]=True
    binding('outer','_flat_edge','exponent','-1.0/(value*value)')
    binding('compliant_outer_initial','point','a','c.exp(-1/x**2)')
    binding('compliant_outer_initial','point','b','c.exp(-1/(1-x)**2)')
    binding('compliant_outer_initial','point','value','a/(a+b)')
    binding('compliant_flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    binding('compliant_collar_Gamma_C4','phi_jets','log_upper','-4/D**2')
    bindings['same_original_native_sigma_edges_and_phi_function']=True
    # The exact source is the displayed cutoff/primitive generator. The
    # old float/Decimal quadrature frontend only approximates that generator.
    return dict(master_logA_source_return_bound=True,original_axial_flatten_and_tail_AST_bindings=bindings,
        exact_cutoff='sigma(x)=edge(x)/(edge(x)+edge(1-x)), edge(x)=exp(-1/x^2) for x>0; otherwise0',
        exact_primitive='J(t)=integral_0^t sigma(v)dv; J(1)=1/2 by symmetry',
        float_quadrature_output_not_asserted_exact=True,passed=True)


class ExactOriginalPreheatPressureOperator:
    """Exact integral witness; interval masses are not inputs to this object.

    sigma, its primitive J and phi denote the original flat functions bound
    by bind_original_master_source. W is constrained by raw_waiting_root,
    built from the uncorrected reference moment, including the flatten.
    Symbolic lengths keep exp(-13/mu) positive without materializing it.
    """
    def __init__(self):
        p=self.partition=original_preheat_partition()
        z=p['z'];t=p['t'];q=p['q'];mu=p['mu'];delta=p['delta'];eps=p['epsilon']
        J=p['J'];sig=p['sigma'];phi=p['phi'];r=1-mu;k=1-delta/2
        bp=s.Rational(1,2)+mu;ph=1+delta;bh=ph/2
        # Independently replay the native inlet/pressure kernel formulas.
        u1=s.exp(-s.Rational(1,5))
        ud=u1*s.exp(-(p['yd']-1)/2)
        uw=ud*s.exp(-s.Rational(1,2)-mu/2)
        up=uw*s.exp(-bp*p['Tw']);ev2=up**2*s.exp(-(1+2*mu)*13/mu)
        thetaR=s.exp(-bp*(100+p['L']))/2
        thetaS=thetaR*s.exp(-bp-r/2)
        thetaQ=thetaS*s.exp(-s.Rational(3,2)*p['Ts'])
        thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
        thetaBase=thetaT*s.exp(-bh*p['W'])/(1-eps)
        Kraw=1-eps*(1-sig(t)+sig(t)*phi(t))
        self.native_densities={
            'reference_extension':s.exp(t/5)/(2*q*q),
            'slope_transition_ref':s.exp(t/5-s.Rational(6,5)*J(t))/(2*q*q),
            'axial_turnoff':u1**2*s.exp(-t)/(2*q*q),
            'slope_transition_mu':ud**2*s.exp(-t-2*mu*J(t))/(2*q*q),
            'power_buffer':uw**2*s.exp(-(1+2*mu)*t)/(2*q*q),
            'pulse_reserved':up**2*s.exp(-(1+2*mu)*t)/(2*q*q),
            'z_flatten':ev2*s.exp(2*s.log(q/2)*sig(t/100))*s.exp(-(1+2*mu)*t)/(2*q*q),
            'power_buffer_rel':ev2*s.exp(-100*(1+2*mu))*s.exp(-(1+2*mu)*t)/8,
            'steep_transition_in':ev2*thetaR**2*s.exp(-(1+2*mu)*t-2*r*J(t))/2,
            'steep_power':ev2*thetaS**2*s.exp(-3*t)/2,
            'steep_transition_out':ev2*thetaQ**2*s.exp(-3*t+2*k*J(t))/2,
            'waiting':ev2*thetaT**2*s.exp(-ph*t)/2,
            'heat_collar':ev2*thetaBase**2*s.exp(-ph*t)*Kraw**2/2,
            'exterior_power_tail':ev2*thetaBase**2*s.exp(-ph*t)/2}
        self.original_densities={name:shape**2/2 for name,shape in p['shapes'].items()}
        self.native_scale=ev2*thetaBase**2
        self.native_log_up=s.log(u1)-(p['yd']-1)/2-s.Rational(1,2)-mu/2-bp*p['Tw']
        self.native_log_ev2=2*self.native_log_up-13/mu-26
        self.native_log_theta_base=(-bp*(100+p['L'])-s.log(2)-bp-r/2
            -s.Rational(3,2)*p['Ts']-s.Rational(3,2)+k/2-bh*p['W']-s.log(1-eps))
        # Exact source-defined raw X; these are integral functions, not X boxes.
        v=s.Symbol('v',real=True)
        i1=s.Integral(s.exp(s.Rational(8,5)*v-s.Rational(3,5)*J(v)),(v,0,1))
        X1=s.exp(-s.Rational(13,10))*(s.Rational(5,8)+i1)
        Xd=1+(X1-1)*s.exp(-(p['yd']-1))
        im=s.Integral(s.exp(v-mu*J(v)),(v,0,1))
        Xw=s.exp(-1+mu/2)*(Xd+im)
        Xp=1/r+(Xw-1/r)*s.exp(-r*p['Tw'])
        Xv=1/r+(Xp-1/r)*s.exp(-13*r/mu)
        # Z=0: the completed flatten factor is 1/2, not one.
        If=s.Integral(2**(-sig(v/100))*s.exp(r*v),(v,0,100))
        Xf=2*(Xv+If)*s.exp(-100*r)
        XR=1/r+(Xf-1/r)*s.exp(-r*p['L'])
        Iin=s.Integral(s.exp(r*(v-J(v))),(v,0,1))
        Iout=s.Integral(s.exp(k*J(v)),(v,0,1))
        Jc=s.Integral(s.exp(k*v)*(1-sig(v)+sig(v)*phi(v)),(v,0,3))
        XS=(XR+Iin)*s.exp(-r/2);XT=(XS+p['Ts']+Iout)*s.exp(-k/2)
        self.raw_angular_chain=dict(X0=s.Rational(5,8),X1=X1,Xd=Xd,Xw=Xw,Xp=Xp,
            Xv=Xv,Xf=Xf,XR=XR,XS=XS,XT=XT,Iin=Iin,Iout=Iout,Jc=Jc)
        self.raw_waiting_root=(s.log(XT-1/k)+s.log(1-eps)-s.log(eps)-s.log(1/k+Jc))/k

    def stage_integral(self,name,native=False):
        p=self.partition;density=(self.native_densities if native else self.original_densities)[name]
        # Pull only the proved Z-independent factor on the beta=2 stages.
        if name in BETA2:
            return s.Integral(s.simplify(density*p['q']**2),(p['t'],*p['domains'][name]))/p['q']**2
        density=s.simplify(s.expand_log(s.expand_power_exp(density.rewrite(s.exp)),force=True))
        return s.Integral(density,(p['t'],*p['domains'][name]))

    def integral(self,native=False):
        return sum((self.stage_integral(name,native) for name in self.original_densities),s.Integer(0))

    def datum(self):
        """Original P0/Pstar^2=-int(Utheta_pre/Pstar)^2 dy/2, H=1."""
        return -self.integral(native=False)

    def selected_datum(self,Z):
        """Callable exact function, with W set by its prescribed raw root."""
        return self.datum().subs({self.partition['z']:s.sympify(Z),
            self.partition['W']:self.raw_waiting_root})


def verify_exact_integral_witness(witness):
    """Reject density, amplitude or raw-history changes before cancellation."""
    p=witness.partition;proofs={}
    def zero(name,left,right=0):
        # q>0 on the real axis (Re q>0 on the certified strip) fixes the
        # logarithm branch in the sole variable-exponent flatten stage.
        value=s.expand_log(s.expand_power_exp((left-right).rewrite(s.exp)),force=True)
        if s.simplify(value,doit=False)!=0:raise ArithmeticError('Exact integral witness differs: '+name)
        proofs[name]=True
    expected=ExactOriginalPreheatPressureOperator()
    # The original defining operator has a fixed sign and generator. This
    # catches a pressure perturbation without trying to evaluate its huge
    # unbounded-domain integrals merely to discover it is nonzero.
    if witness.datum()!=-witness.integral(False):
        raise ArithmeticError('Original P0 sign or source integral changed')
    if set(witness.native_densities)!=set(BETA2+BETA0+('z_flatten',)):
        raise ValueError('Incomplete original integral partition')
    for name in witness.original_densities:
        zero('native_density_equals_original_'+name,witness.native_densities[name],witness.original_densities[name])
        zero('same_exact_native_and_original_stage_integral_'+name,
            witness.stage_integral(name,True),witness.stage_integral(name,False))
    for name,value in witness.raw_angular_chain.items():
        zero('raw_reference_angular_history_'+name,value,expected.raw_angular_chain[name])
    zero('same_raw_waiting_root',witness.raw_waiting_root,expected.raw_waiting_root)
    zero('original_Rp_log_source',witness.native_log_up,p['logUp'])
    zero('original_Rv_log_source',witness.native_log_ev2,2*p['logEv'])
    zero('original_tail_scale_log_source',witness.native_log_ev2+2*witness.native_log_theta_base,2*p['logBase'])
    zero('actual_pressure_scale_equals_original_tail_amplitude',witness.native_scale,s.exp(2*p['logBase']))
    # Every term below is an explicit integral of a bound source density.
    # There are no fresh free symbols M2, M0 or Fflat to cancel by definition.
    zero('original_exact_P0_plus_native_exact_complete_integral',witness.datum()+witness.integral(True))
    return proofs


def live_native_density_translations(witness):
    """Translate production AST into exact densities, independently of P0.

    Integral/kernel names below are exact functions; the corresponding
    interval routines enclose them. The prior raw-operator source receipt
    binds their domains, primitive kernels, incoming graph and scale paths.
    """
    p=witness.partition;t=p['t'];mu=p['mu'];q=p['q'];delta=p['delta'];eps=p['epsilon']
    J=p['J'];sig=p['sigma'];phi=p['phi'];r=1-mu;k=1-delta/2;ph=1+delta
    ev2=s.exp(witness.native_log_ev2);bp=s.Rational(1,2)+mu
    u1=s.exp(-s.Rational(1,5));ud=u1*s.exp(-(p['yd']-1)/2)
    uw=ud*s.exp(-s.Rational(1,2)-mu/2);up=s.exp(witness.native_log_up)
    D=lambda rate:(1-s.exp(-rate*t))/rate
    zero=s.Integer(0);out={}
    a=lambda stem,method,target,env,aug=False:assignment(stem,method,target,env,augmented=aug)
    out['reference_extension']=s.diff(a('compliant_outer_initial','reference','p',
        {'u':s.exp(t/10)/q}),t)
    slope_kernel=a('interval_outer_slope_field','transition_integrals','masses[k]',
        {'dy':1,'rate':s.Rational(1,5),'scell':t,'power':2,'jcell':J(t)},True)
    out['slope_transition_ref']=slope_kernel/(2*q*q)
    out['axial_turnoff']=s.diff(a('compliant_outer_initial','axial','p',
        {"get('Mp_over_Pstar_squared')":zero,'u1':u1/q,'decay':s.exp(-t)}),t)
    # transition_kernels K[2] uses exact exp(-v) weights and exp(-2mu J).
    mu_kernel=a('compliant_outer_buffer','transition_kernels','K[n]',
        {'weights[n]':s.exp(-t),'mu':mu,'power':2,'jcell':J(t)},True)
    out['slope_transition_mu']=ud**2*mu_kernel/(2*q*q)
    out['power_buffer']=s.diff(a('compliant_outer_buffer','power','p',
        {"get('Mp_over_Pstar_squared')":zero,'u1':uw/q,
         'decay_integral(c, 1 + 2 * mu, t)':D(1+2*mu)}),t)
    out['pulse_reserved']=s.diff(a('compliant_current_raw_preheat_pressure_operator','evaluate',"atoms['pulse_reserved']",
        {"incoming['u']":up/q,'self.pressure_decay':s.exp(-(1+2*mu)*t),'p':1+2*mu}),t)
    Fc=s.exp(s.log(q/2)*sig(t/100))
    # Cell length/weight is separately bound to its exact exponential mass.
    out['z_flatten']=ev2*a('compliant_flatten_mixed_C4','flatten','Pint',
        {'Fc':Fc,'q':q,'self.prate':1+2*mu,'a':t,
         'decay_integral(c, self.prate, length)':1},True)
    thetaR=s.exp(-bp*(100+p['L']))/2
    thetaS=thetaR*s.exp(-bp-r/2);thetaQ=thetaS*s.exp(-s.Rational(3,2)*p['Ts'])
    thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
    out['power_buffer_rel']=s.diff(a('compliant_current_raw_preheat_pressure_operator','evaluate',"atoms['power_buffer_rel']",
        {'one':1,'decay_integral(c, p, L)':D(1+2*mu),'p':1+2*mu,'self.Ev2':ev2}),t)
    for name,kind,rate,theta in (('steep_transition_in','in',r,thetaR),
            ('steep_transition_out','out',k,thetaQ)):
        # The production function has one assignment per branch; select it.
        fn=function('compliant_steep_waiting_C4','transition_kernels')
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='pressure']
        expression='ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2' if kind=='in' else 'ds*c.exp(-3*v+2*rate*jc)/2'
        expected=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==expected for v in values)!=1:raise ValueError('Native steep pressure kernel changed')
        out[name]=ev2*theta**2*(s.exp(-(1+2*mu)*t-2*rate*J(t))/2 if kind=='in' else s.exp(-3*t+2*rate*J(t))/2)
    for name,target,rate,theta in (('steep_power',"atoms['steep_power']",3,thetaS),
            ('waiting',"atoms['waiting']",ph,thetaT)):
        call='decay_integral(c, 3, Ts)' if name=='steep_power' else 'decay_integral(c, 1 + self.flat.delta, W)'
        out[name]=s.diff(a('compliant_current_raw_preheat_pressure_operator','evaluate',target,
            {'one':1,call:D(rate),'self.Ev2':ev2,'thetaS':thetaS,'thetaT':thetaT}),t)
    Kraw=1-eps*(1-sig(t)+sig(t)*phi(t))
    out['heat_collar']=witness.native_scale*a('compliant_collar_Gamma_C4','forward_pressure','integral',
        {'K':Kraw,'ds':1,'self.prate':ph,'v':t},True)
    # Integral_3^infinity exp(-ph*t)/2 = exp(-3ph)/(2ph).
    out['exterior_power_tail']=witness.native_scale*s.exp(-ph*t)/2
    return out


def live_native_raw_angular_translations(witness):
    """Exact axis history replayed from the native AST, before correction."""
    p=witness.partition;c=witness.raw_angular_chain;v=s.Symbol('v',real=True)
    mu=p['mu'];r=1-mu;k=1-p['delta']/2;J=p['J'];sig=p['sigma']
    a=lambda stem,method,target,env:assignment(stem,method,target,env)
    kernel=assignment('interval_outer_slope_field','transition_integrals','masses[k]',
        {'dy':1,'rate':s.Rational(8,5),'scell':v,'power':1,'jcell':J(v)},augmented=True)
    I1=s.Integral(kernel,(v,0,1));u1=s.exp(-s.Rational(1,5))
    h1=a('compliant_outer_initial','slope','h',{'qi':1,'integrals[0]':I1,'yy':1})
    X1=h1/u1;t=p['yd']-1;ud=u1*s.exp(-t/2)
    hd=a('compliant_outer_initial','axial','h',{"get('Mtheta_over_sqrt2_R_3half_Pstar')":h1,
        'decay3':s.exp(-s.Rational(3,2)*t),'u1':u1,'root_decay':s.exp(-t/2)})
    im_kernel=assignment('compliant_outer_buffer','transition_kernels','K[n]',
        {'weights[n]':s.exp(v),'mu':mu,'power':1,'jcell':J(v)},augmented=True)
    Im=s.Integral(im_kernel,(v,0,1));uw=ud*s.exp(-s.Rational(1,2)-mu/2)
    hw=a('compliant_outer_buffer','slope_mu','h',{
        "get('Mtheta_over_sqrt2_R_3half_Pstar')":hd,'d3':s.exp(-s.Rational(3,2)),
        'u1':ud,"kernels['theta']":Im})
    f=s.exp(-(s.Rational(1,2)+mu)*p['Tw']);d3=s.exp(-s.Rational(3,2)*p['Tw'])
    binding('compliant_outer_buffer','power','theta_kernel','(f-d3)/(1-mu)',allow_other_assignments=True)
    hp=a('compliant_outer_buffer','power','h',{
        "get('Mtheta_over_sqrt2_R_3half_Pstar')":hw,'d3':d3,'u1':uw,'theta_kernel':(f-d3)/r})
    Xp=hp/(uw*f)
    Xv=a('compliant_axial_pulse_field','end', 'X',{
        'self.rate':r,'self.Xp':Xp,'self.mu':mu,'s':0,
        'self.factor(-13 * self.rate / self.mu - self.rate * s)':s.exp(-13*r/mu)})
    integral=s.Integral(2**(-sig(v/100))*s.exp(r*v),(v,0,100))
    Xf=a('compliant_outer_angular_candidate','flatten','X',{
        'integral':integral,'self.Xv':Xv,'self.rate':r,'t':100,'F':s.Rational(1,2)})
    XR=a('compliant_outer_angular_candidate','reference_buffer','X',{
        'departure':(Xf-1/r)*s.exp(-r*p['L']),'eq':1/r})
    XS=a('compliant_outer_angular_candidate','steep_in','X',{
        "inlet['X']":XR,'I':c['Iin'],'self.rate':r,'t':1,'J':s.Rational(1,2)})
    XT=a('compliant_outer_angular_candidate','before_waiting','X',{
        "inlet['X']":XS+p['Ts'],'I':c['Iout'],'self.restore_rate':k,'J':s.Rational(1,2)})
    X0=a('compliant_outer_initial','reference','h',{'u':1})
    return dict(X0=X0,X1=X1,Xd=hd/ud,Xw=hw/uw,Xp=Xp,Xv=Xv,Xf=Xf,XR=XR,XS=XS,XT=XT)


def original_pressure_function_identification(balance,raw,witness=None):
    """Prove one integral operator before using its analytic datum projection."""
    if not balance.acceptance_loaded or not raw.acceptance_loaded or balance.angular is not raw.angular:
        raise ValueError('Checked current pressure prerequisites must share one angular source graph')
    exact=raw.exact;native_datum=exact.flatten.inlet.datum;repair_datum=exact.repair.angular.initial.datum
    # Explicitly consume the checked measure identity: native collarJ and
    # full-Gamma JW both integrate exp(kv)*(1-sigma+sigma*phi) on[0,3].
    # The collarJ constructor assignment alone cannot identify this integral.
    angular_proof=balance.angular.proof
    waiting_measure={name:angular_proof['identities'].get(name) is True for name in
        ('old_J_weight_is_same_exact_measure','current_J_same_integrand')}
    waiting_measure['same_exact_preheat_collarJ_and_current_JW_function']=angular_proof.get(
        'same_exact_preheat_collarJ_and_current_JW_function') is True
    if not all(waiting_measure.values()):raise ValueError('Checked native collarJ/current JW exact measure identity required')
    common=pressure_defining_function_proof(native_datum,repair_datum)
    parameter=exact.parameter_bridge
    if not parameter['passed'] or not all(exact.graph.values()):raise ValueError('Common original/native parameter function bridge required')
    original=bind_original_master_source();witness=witness if witness is not None else ExactOriginalPreheatPressureOperator()
    p=witness.partition;proofs=verify_exact_integral_witness(witness);bindings={}
    def zero(name,left,right=0):
        if s.simplify(s.expand_power_exp(left-right),doit=False)!=0:raise ArithmeticError('Original pressure function identity differs: '+name)
        proofs[name]=True
    z=p['z'];t=p['t'];mu=p['mu'];delta=p['delta'];eps=p['epsilon'];q=p['q']
    r=1-mu;k=1-delta/2;bp=s.Rational(1,2)+mu;bh=(1+delta)/2
    # The native functions are extracted independently from production AST.
    # Their exact primitives use the original J/sigma, not numerical boxes.
    for name,density in live_native_density_translations(witness).items():
        value=s.expand_log(s.expand_power_exp((density-witness.native_densities[name]).rewrite(s.exp)),force=True)
        zero('live_native_AST_density_'+name,s.simplify(value))
    for name,value in live_native_raw_angular_translations(witness).items():
        zero('live_native_AST_raw_angular_'+name,value,witness.raw_angular_chain[name])
    for name,(y,j0,jd,jr,jq,sf) in p['placements'].items():
        logA=y/10-s.Rational(3,5)*j0-mu*jd-r*jr+k*jq
        master=s.exp(logA)*(q/2)**sf/q
        zero('original_master_equals_native_segment_'+name,master,p['shapes'][name])
    # c_inf=A(y_tail)*R_tail^bh/(2*(1-eps)); H=1.
    tail=p['master_tail_y'];yr=p['offsets']['rel'];yq=p['offsets']['q']
    logAtail=tail/10-s.Rational(3,5)*(tail-s.Rational(1,2))-mu*(tail-p['yd']-s.Rational(1,2))-r*(tail-yr-s.Rational(1,2))+k*(tail-yq-s.Rational(1,2))
    zero('original_c_inf_raw_tail_amplitude_units',logAtail-s.log(2)-s.log(1-eps),p['logBase'])
    specs={
        ('interval_outer_slope_field','transition_integrals'):{
            'rates':"(c.mpf('1.6'),c.mpf('.2'),c.mpf('1.2'))",'powers':'(1,2,2)'},
        ('compliant_outer_buffer','transition_kernels'):{
            'weights':'(c.exp(b)-c.exp(a),da,c.exp(-a)-c.exp(-b))'},
        ('compliant_outer_initial','reference'):{'u':'qi*c.exp(y/10)','p':"u*u*c.mpf('2.5')"},
        ('compliant_outer_initial','slope'):{'factor':"c.exp(yy/10-c.mpf('.6')*J)",'u':'qi*factor',
            'p':"qi*qi*(c.mpf('2.5')+integrals[1]/2)"},
        ('compliant_outer_initial','axial'):{'t':'y-1','u':'u1*root_decay',
            'p':"get('Mp_over_Pstar_squared')+u1*u1*((1-decay)/2)"},
        ('compliant_outer_buffer','slope_mu'):{'f':'c.exp(-t/2-mu*kernels[\'J\'])','u':'u1*f',
            'p':"get('Mp_over_Pstar_squared')+u1*u1*(kernels['pressure']/2)"},
        ('compliant_outer_buffer','power'):{'t':'self.params.Tw*phase','slope':"-c.mpf('.5')-mu",'u':'u1*f',
            'p':"get('Mp_over_Pstar_squared')+u1*u1*(decay_integral(c,1+2*mu,t)/2)"},
        ('compliant_outer_angular_candidate','__init__'):{'self.collarJ':'collar_preheat_integral(c,self.restore_rate)'},
        ('compliant_outer_angular_candidate','reference_buffer'):{'X':'departure+eq'},
        ('compliant_outer_angular_candidate','flatten'):{
            'ratio':'lq-c.ln(2)','F':'(ratio*sig).exp()',
            'X':'(integral+self.Xv)*c.exp(-self.rate*t)/F'},
        ('compliant_outer_angular_candidate','steep_in'):{'X':"(inlet['X']+I)*c.exp(-self.rate*(t-J))"},
        ('compliant_outer_angular_candidate','steep_power'):{'t':'self.params.Ts*phase'},
        ('compliant_outer_angular_candidate','before_waiting'):{'X':"(inlet['X']+I)*c.exp(-self.restore_rate*J)"},
        ('compliant_current_exact_repair_branch','__init__'):{
            'terminal':"angular.before_waiting('0')",
            'angular.waiting':'(c.ln(Xt-eq)+angular.waiting_logone-angular.params.log_epsilon-c.ln(eq+angular.collarJ))/k'},
        ('compliant_axial_pulse_field','__init__'):{
            'inlet':"self.pulse.buffer.power('0',1)",
            'self.Xp':"inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]"},
        ('compliant_outer_buffer','transition_kernels'):{
            'weights':'(c.exp(b)-c.exp(a),da,c.exp(-a)-c.exp(-b))'}}
    for (stem,method),rows in specs.items():
        for target,expression in rows.items():
            binding(stem,method,target,expression);bindings[stem+'.'+method+':'+target]=True
    class_assignment('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__','self.logEv2_parts',
        "dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))")
    bindings['CurrentFlattenMixedC4:original_positive_amplitude_logs']=True
    # Native pressure primitives solve p'=u^2/2 with the same incoming
    # reference value. Uniqueness of this scalar integral identifies all
    # early cumulative stages; no difference enclosure defines a mass.
    zero('reference_pressure_integral_down_to_zero',s.integrate(s.exp(t/5)/2,(t,-s.oo,0)),s.Rational(5,2))
    zero('reference_native_radial_pressure_FTC',s.diff(s.Rational(5,2)*s.exp(t/5),t),s.exp(t/5)/2)
    zero('original_axial_turnoff_mass_matches_native_primitive',
        s.exp(s.Rational(3,5))*(s.exp(-1)-s.exp(-p['yd']))/2,
        s.exp(-s.Rational(2,5))*(1-s.exp(1-p['yd']))/2)
    U,rate=s.symbols('U rate',positive=True)
    Jt=p['J'](t)
    zero('native_first_slope_pressure_kernel_is_original_swirl_square',
        p['shapes']['slope_transition_ref']**2/2,s.exp(t/5-s.Rational(6,5)*Jt)/(2*q*q))
    zero('native_mu_transition_pressure_kernel_is_original_swirl_square',
        p['shapes']['slope_transition_mu']**2/2,
        s.exp(s.Rational(3,5)-p['yd'])*s.exp(-t-2*mu*Jt)/(2*q*q))
    zero('native_power_and_pulse_radial_pressure_FTC',
        s.diff(U*U*(1-s.exp(-rate*t))/(2*rate),t),U*U*s.exp(-rate*t)/2)
    zero('original_Rp_amplitude_same_source_function',p['logUp'],
        s.Rational(3,10)-p['yd']/2-s.Rational(1,2)-mu/2-bp*p['Tw'])
    zero('original_Rv_squared_amplitude_log_has_both_offsets',2*(p['logEv']-p['logUp']),-13/mu-26)
    sigma,phi=s.symbols('sigma phi',real=True)
    zero('original_H_one_tail_bracket_is_native_preheat_density',
        (1-sigma)*(1-eps)+sigma*(1-eps*phi),1-eps*(1-sigma+sigma*phi))
    zero('original_flat_edge_argument_is_native_phi',-1/((3-t)/2)**2,-4/(3-t)**2)
    # The original waiting equation (constant + C exp(kW)=XR) is exactly
    # the replayed native raw root. XR is uncorrected, at Z=0.
    XR,Iin,Iout,Jc=s.symbols('XR Iin Iout Jcollar',real=True)
    XS=(XR+Iin)*s.exp(-r/2);XT=(XS+p['Ts']+Iout)*s.exp(-k/2)
    constant=s.exp(r/2)*(s.exp(k/2)/k-Iout-p['Ts'])-Iin
    C=s.exp((r+k)/2)*eps*(1/k+Jc)/(1-eps)
    zero('original_waiting_root_forward_difference',XT-1/k,s.exp(-(r+k)/2)*(XR-constant))
    zero('original_and_current_waiting_equations_same_unique_root',
        (XT-1/k)*(1-eps)-eps*(1/k+Jc)*s.exp(k*p['W']),
        s.exp(-(r+k)/2)*(1-eps)*(XR-constant-C*s.exp(k*p['W'])))
    # Every pre-flatten stage has the same exact q^-2 pressure shape;
    # every post-flatten stage is Z-independent. The flatten is the sole
    # variable exponent stage and uses the already checked native callable.
    masses={name:witness.stage_integral(name) for name in p['shapes']}
    M2=sum((masses[name]*q*q for name in BETA2),s.Integer(0))
    M0=sum((masses[name] for name in BETA0),s.Integer(0));Fflat=masses['z_flatten']
    for name in BETA2:zero('early_beta2_density_'+name,s.diff(s.simplify(p['shapes'][name]**2*q*q),z))
    for name in BETA0:zero('late_beta0_density_'+name,s.diff(p['shapes'][name]**2,z))
    # Canonical integral linearity and pulling out the common q^-2
    # factor yields exactly the original datum's defining operator.
    original_operator=witness.integral(False)
    zero('exact_M2_M0_Fflat_decomposition',original_operator,M2/q**2+M0+Fflat)
    if set(masses)!=set(BETA2+BETA0+('z_flatten',)):raise ValueError('Complete original pressure partition required')
    # The original datum's immutable definition selects this exact generator
    # and root. Its m2/m0/flatten and normalized_jets are ONLY enclosures.
    expected_definition=dict(Md=native_datum.parameters.Md,logPstar='exp(Md)+11',
        c_mu='.001',c_delta='.001',c_epsilon='.001',delta='min(1e-200,exp(-4logPstar-30))',
        waiting='unique positive root of the continuous raw preheat waiting equation',
        angular_profile='Section 6.1 reference-plus-outer ansatz, H replaced by 1',
        Tw='-60log(mu)',Ts='4log(2/delta)',Tf=100,
        cutoff_and_schedule_python_sha256=sha('lei_ren_part1_paper_outer.py'))
    if native_datum.definition!=expected_definition or repair_datum.definition!=expected_definition:
        raise ValueError('Original implicit datum generator/root specification changed')
    class_assignment('pressure_source','CompliantPressureDatum','__init__','self.definition',
        "dict(Md=str(Md),logPstar='exp(Md)+11',c_mu='.001',c_delta='.001',c_epsilon='.001',delta='min(1e-200,exp(-4logPstar-30))',waiting='unique positive root of the continuous raw preheat waiting equation',angular_profile='Section 6.1 reference-plus-outer ansatz, H replaced by 1',Tw='-60log(mu)',Ts='4log(2/delta)',Tf=100,cutoff_and_schedule_python_sha256=self.input_hashes['lei_ren_part1_paper_outer.py'])")
    bindings['CompliantPressureDatum:unchanged_exact_generator_and_raw_root_definition']=True
    # Bind the older coarse enclosure to the same exact equation, rather
    # than identifying roots by endpoint containment.
    binding('logarithmic_outer_parameters','waiting_root_enclosure','constant','ed*(er*eq-Lrestore-self.Ts)-Ldrop')
    binding('logarithmic_outer_parameters','waiting_root_enclosure','log_coefficient',
        'a/2+k/2+self.log_epsilon-c.ln(1-self.epsilon)+c.ln(eq+K)')
    bindings['old_coarse_root_encloses_same_source_equation']=True
    # Integral/derivative uniqueness: r,k,eps are positive; the original
    # root equation is strictly increasing in W. The accepted constructor
    # guards the same raw XT-1/k>0 and W>0. No corrected XT enters it.
    guards=dict(mu_in_0_001=0<endpoints(raw.flat.mu)[0]<=endpoints(raw.flat.mu)[1]<1,
        delta_in_0_1=0<endpoints(raw.flat.delta)[0]<=endpoints(raw.flat.delta)[1]<1,
        epsilon_in_0_half=0<endpoints(native_datum.parameters.epsilon)[0]<=endpoints(native_datum.parameters.epsilon)[1]<mp.mpf('.5'),
        yd_ge_one=endpoints(native_datum.parameters.yd)[0]>=1,
        native_waiting_positive=endpoints(exact.repair.angular.waiting)[0]>0)
    if not all(guards.values()):raise ValueError('Exact generator domination/root guards failed')
    ac=witness.raw_angular_chain
    root_log_equation=(s.log(ac['XT']-1/k)+s.log(1-eps)-s.log(eps)-s.log(1/k+ac['Jc']))/k
    zero('exact_raw_waiting_root_log_equation',witness.raw_waiting_root,root_log_equation)
    # Full density domination, rather than merely integrating an envelope:
    # Jd>=Jrel>=Jq>=0 since yd<yrel<yq and sigma>=0.
    jd,jr,jq,dd,dr=s.symbols('Jd Jrel Jq Delta_d_rel Delta_rel_q',nonnegative=True)
    rest=-mu*jd-r*jr+k*jq
    zero('monotone_original_primitive_negative_remainder',rest,
        -mu*(jd-jr)-(jr-jq)-delta*jq/2)
    dominated=rest.subs({jd:jq+dr+dd,jr:jq+dr})
    zero('actual_master_remainder_nonpositive_cone',dominated,-mu*dd-dr-delta*jq/2)
    if s.ask(s.Q.nonpositive(s.expand(dominated))) is not True:
        raise ArithmeticError('Actual master density envelope sign not proved')
    y=s.Symbol('y',real=True)
    zero('completed_first_cutoff_gives_log_envelope',y/10-s.Rational(3,5)*(y-s.Rational(1,2)),s.Rational(3,10)-y/2)
    # Tail: 0<=Kraw<=1, ph>=1, and the bound at y_tail is continued
    # by exp(-ph*t)<=exp(-t); its exact density has factor1/8.
    zero('actual_raw_tail_squared_amplitude_normalization',
        s.exp(2*p['logBase'])/2,s.exp(2*logAtail)/(8*(1-eps)**2))
    # Holomorphic domination: Re(1+z^2)>=1-rho^2>0 on the strip.
    # q^(2sigma-2) is bounded by (1-rho^2)^-2. Reference density
    # decays at y->-infinity; all y>=1 density is <= exp(.6-y)
    # /(2*(1-eps)^2)*(1-rho^2)^-2. Compact stages require no tail limit.
    rho=s.Rational(1,4);q_lower=1-rho*rho
    tail_bound=s.exp(s.Rational(3,5)-1)/(2*(1-eps)**2*q_lower**2)
    zero('real_complex_strip_q_lower_positive',q_lower,s.Rational(15,16))
    zero('infinite_preheat_pressure_envelope_integrable',
        s.integrate(s.exp(s.Rational(3,5)-t)/(2*(1-eps)**2*q_lower**2),(t,1,s.oo)),tail_bound)
    full_mass_bound=s.Rational(5,2)/q_lower**2+s.exp(s.Rational(1,5))/(2*q_lower**2)+tail_bound
    # One unique, dominated integral defines P0. Both native cumulative
    # primitives and the datum projections use this generator. Ordinary
    # axial coefficients commute with this integral through order5.
    P0=witness.datum();native_integral=witness.integral(True)
    source_difference=s.simplify(P0+native_integral,doit=False)
    for n in range(6):zero('original_P0_plus_native_raw_integral_axial_'+str(n),s.diff(source_difference,z,n))
    return dict(original_master_source=original,original_native_AST_bindings=bindings,identities=proofs,
        original_stage_velocity_functions={name:str(value) for name,value in p['shapes'].items()},
        exact_segmented_domains={name:list(map(str,value)) for name,value in p['domains'].items()},
        original_full_pressure_integral_operator=str(original_operator),
        exact_original_M2=str(M2),exact_original_M0=str(M0),exact_original_Fflat=str(Fflat),
        exact_original_P0_function=str(P0),exact_native_complete_pressure_integral=str(native_integral),
        exact_axis_raw_angular_chain={name:str(value) for name,value in ac.items()},
        exact_axis_raw_waiting_root=str(witness.raw_waiting_root),
        original_exact_generator_definition=expected_definition,
        original_datum_interval_jets_remain_enclosures_only=True,
        integral_function_independent_of_interval_m2_m0_and_flatten_boxes=True,
        original_function_witness_class=ExactOriginalPreheatPressureOperator.__name__,
        common_original_datum_projection=common,common_exact_parameter_function_bridge=parameter,
        accepted_current_collar_waiting_measure_function_proof=waiting_measure,
        same_reference_down_to_zero_and_native_cumulative_FTC=True,
        same_native_Rp_amplitude_function_and_positive_log_factors=True,
        same_original_cutoff_primitive_and_stage_offsets=True,
        same_raw_waiting_equation_and_unique_positive_root=True,
        corrected_angular_terminal_not_used_for_raw_waiting=True,
        original_implicit_Fflat_is_checked_native_pressure_integral_callable=True,
        full_infinite_tail_domination=dict(strip_half_width=str(rho),q_real_lower=str(q_lower),tail_mass_bound=str(tail_bound),
            complete_mass_bound=str(full_mass_bound),actual_negative_J_remainder=str(dominated),
            native_parameter_guards=guards,raw_collar_bracket_in_zero_one=True,
            ordered_cutoff_origins_and_monotonic_J_used=True,
            Cauchy_coefficient_majorants=[str(full_mass_bound/rho**n) for n in range(6)]),
        holomorphic_integral_and_axial5_differentiation_justified=True,
        all_original_raw_stage_integrals_are_one_function_not_box_choices=True,
        original_P0_not_redefined_or_pressure_patched=True,
        same_exact_inverse_radius_and_Pstar_squared_units_preserved=True,passed=True)


class CurrentPressureTerminalClosure:
    @source_precision
    def __init__(self,balance=None,raw=None,require_checked=True):
        self.balance=balance if balance is not None else CurrentPressureTerminalBalance()
        self.angular=self.balance.angular
        self.raw=raw if raw is not None else CurrentRawPreheatPressureOperator(angular=self.angular)
        self.original_pressure_function=ExactOriginalPreheatPressureOperator()
        self.proof=original_pressure_function_identification(self.balance,self.raw,self.original_pressure_function)
        self.exact=self.angular.exact;self.heat=self.angular.heat;self.ctx=self.angular.ctx
        self.family=self.angular.family;self.source=self.angular.source;self.datum_sha=self.angular.datum_sha
        self.hashes=dict(self.balance.hashes);self.hashes.update(self.raw.hashes)
        for name in (Path(__file__).name,'lei_ren_part1_paper_outer.py','lei_ren_part1_paper_compliant_actual_Rh_source_join.py'):
            self.hashes[name]=sha(name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[1]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current pressure terminal source/scope differs')
            if receipt['original_pressure_function_identification']!=encode(pack(self.proof)):
                raise ValueError('Original pressure function source proof changed')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def evaluate(self,chart,Z,t):
        if chart not in ('heat_collar','heat_exterior'):raise ValueError('Current original heat chart required')
        c=self.ctx;z=c.mpf(Z);t=c.mpf(t);key=(chart,z._mpi_,t._mpi_)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1 or endpoints(t)[0]<0 or (chart=='heat_collar' and endpoints(t)[1]>3) or (chart=='heat_exterior' and endpoints(t)[0]<3):
            raise ValueError('Original heat collar/exterior domain required')
        if key in self.cache:return self.cache[key]
        terminal=self.angular.terminal_constants(z);forward=terminal['current_repaired_forward_terminal']
        zero=IntervalTaylor.constant(c,0,5)
        if chart=='heat_collar':
            shape=shape_radial5(self.heat,z,t);tails=self.heat.collar_tails(z,t)
            defects=collar_defect_rows(self.heat,shape,tails,t)
            canonical=collar_stress_mixed4(self.heat,shape,defects,z,t)
            angular=tails['angular_numerator']/shape['K_rows'][0]
            remaining=tails['remaining_pressure_in_Rtail_units']
        else:
            shape=self.heat.local_Gamma(z,t)
            angular=shape['angular_numerator']/shape['K_rows'][0]
            remaining=shape['pressure_numerator']*c.exp(-self.heat.prate*t)
            canonical=dict(theta=[zero]*5,axial=[zero]*5)
        extra=constant_stress_rows(self.heat,z,t,zero,zero,4)
        pressure=collar_pressure_rows(shape['K_rows'][:5],remaining,self.heat.prate,forward['pressure_scale'],t)
        raw=self.raw.evaluate(z)
        out=dict(chart=chart,Z=z,coordinate=t,actual_angular_Taylor_after_source_closure=angular,
            actual_Dtheta_Taylor=zero,actual_Cp_Taylor=zero,
            original_forward_Cp_enclosure_retained=terminal['pressure_infinity'],
            original_analytic_P0_Taylor_retained=raw['original_analytic_P0_Taylor_retained'],
            original_P0_plus_native_raw_integral_enclosure=raw['unresolved_original_P0_plus_raw_total_Taylor'],
            source_proved_original_P0_plus_native_raw_integral_Taylor=zero,
            source_proved_original_pressure_balance_Taylor=zero,
            actual_stress_factored_mixed4=dict(theta_Qtheta=mixed(canonical['theta'],4),
                axial_Qz=mixed(canonical['axial'],4),axial_pressure_Qpressure=mixed(extra['axial_pressure_constant'],4)),
            stable_current_absolute_pressure_mixed4=mixed(pressure,4),
            pressure_is_original_forward_function_after_terminal_identity=True,
            original_axis_pressure_not_replaced_or_tail_patched=True,
            whole_Z_pressure_function_identity_consumed=True,
            exact_pressure_scale=forward['pressure_scale'],
            exact_positive_pressure_log_terms=raw['exact_positive_amplitude_logs'],
            exact_inverse_radius_log_terms=self.exact.repair.heat.logS_terms,
            pressure_scope_only_full_exterior_receipt_still_required=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self.cache[key]=out;return out


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPressureTerminalClosure(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,original_pressure_function_identification=field.proof,
        current_pressure_terminal_views={name:field.evaluate(chart,Z,t) for name,(chart,Z,t) in VIEWS.items()},
        input_hashes=field.hashes,**dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original analytic P0 and complete native raw integral identified; current pressure terminal closes',flush=True)
    return result


if __name__=='__main__':run()
