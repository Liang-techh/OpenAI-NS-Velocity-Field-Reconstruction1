"""Current O7 source owners after admitted power/angular, original algorithms."""
import ast
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_power_angular_source import (
    CurrentPowerAngularSourceAssembly, CHARTS as PRIOR_CHARTS,
    UNIFORM,SCOPES,OPEN,HERE,PREFIX,sha)
from lei_ren_part1_paper_compliant_steep_waiting_C4 import CompliantSteepWaitingC4 as BASE
from lei_ren_part1_paper_compliant_steep_waiting_C4_check import functional_source_identities
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

METHODS={'steep_entry':'steep_in','steep_power':'steep_power','steep_exit':'steep_out','waiting':'waiting'}
NEW_CHARTS=tuple(METHODS)
CHARTS=PRIOR_CHARTS+NEW_CHARTS
RECEIPT=PREFIX+'current_steep_waiting_source_check.json'
GATE='current_steep_waiting_source_ownership_certified'
PROVED='current_angular_steep_waiting_source_functional_joins_proved'
THEOREM=PREFIX+'steep_waiting_C4_check.json'


class CurrentSteepWaitingC4(BASE):
    @source_precision
    def __init__(self,source,cells=128):
        self.outer=source.outer;self.future=f=self.outer.future
        self.ctx=c=self.outer.ctx;self.mu=self.outer.mu;self.delta=self.outer.delta
        self.rate=1-self.mu;self.k=1-self.delta/2
        self.bp=c.mpf('.5')+self.mu;self.bh=c.mpf('.5')+self.delta/2
        self.family=source.family;self.source=source.source;self.cells=cells
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive directed integration cell count required')
        self.hashes=dict(source.hashes);self.cache={};self.kernel_cache={}
        box=lambda value:c.mpf(endpoints(value))
        self.wait=box(f.angular.waiting);self.logone=box(f.angular.waiting_logone)
        self.Ts=box(f.params.Ts);self.epsilon=box(f.heat.epsilon)
        atoms=f.heat.preheat_atoms['energy']
        self.atoms=dict(baseline=1/box(f.delta),epsilon_atom=-2*self.epsilon*box(atoms['W']),
            epsilon_squared_atom=self.epsilon**2*box(atoms['W_squared']))
        self.preheat_scalar=sum(self.atoms.values(),c.mpf(0))
        self.tail_normalization=c.exp(-2*self.logone);self.S=box(f.repair.strong_S_cap)
        self.infull=self.kernels(1,'in');self.outfull=self.kernels(1,'out')
        self.inenergy=self.kernels(0,'in')['remaining_energy'];self.outenergy=self.kernels(0,'out')['remaining_energy']
        self.thetaR=c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2
        self.thetaS=self.thetaR*c.exp(-self.bp-self.rate/2)
        self.thetaQ=self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)
        self.thetaT=self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)

    @source_precision
    def data(self,Z):
        c=self.ctx;Z=c.mpf(Z);key=tuple(endpoints(Z))
        if key[0]<-1 or key[1]>1:raise ValueError('Z in[-1,1] required')
        if key in self.cache:return self.cache[key]
        terminal=self.outer.angular(Z,0)
        source=self.outer.fifth.angular(Z)
        if not source['admitted_C4_prefix_preserved']:raise ValueError('Same actual C4/C5 angular prefix required')
        heat=copy_jet(c,source['scaled_Gamma_future_defect_Taylor']['energy'])
        preheat_tail=IntervalTaylor.constant(c,self.preheat_scalar,5)-heat*c.mpf([0,endpoints(self.delta*self.S/2)[1]])
        H=preheat_tail*self.tail_normalization
        waiting_future=H*c.exp(-self.delta*self.wait)+decay_integral(c,self.delta,self.wait)
        after_power=waiting_future*c.exp(-1-self.delta/2)+self.outenergy
        if endpoints(H[0])[0]<=0 or endpoints(after_power[0])[0]<=c.mpf('.5'):
            raise ArithmeticError('Original heat tail / uniform power energy floor not proved')
        after_entry=(after_power*c.exp(-2*self.Ts)+decay_integral(c,2,self.Ts))*c.exp(-1-self.mu)
        XR=terminal['angular_Taylor'];PR=terminal['pressure_over_Pstar_squared_Taylor']
        XS=(XR+self.infull['angular'])*c.exp(-self.rate/2)
        XQ=XS+self.Ts;XT=(XQ+self.outfull['angular'])*c.exp(-self.k/2)
        PS=PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)
        PQ=PS+decay_integral(c,3,self.Ts)*(self.outer.flatten.Ev2*self.thetaS**2/2)
        PT=PQ+self.outfull['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)
        item=dict(H=H,waiting_future=waiting_future,after_power=after_power,after_entry=after_entry,
            XR=XR,XS=XS,XQ=XQ,XT=XT,PR=PR,PS=PS,PQ=PQ,PT=PT,
            original_angular_terminal=terminal,
            current_live_Gamma_C5_and_admitted_C4_prefix_used=True)
        self.cache[key]=item;return item

    def packet(self,*args,**kwargs):
        packet=BASE.packet(self,*args,**kwargs)
        packet.update(full_pulse_C4_installed=False,current_steep_waiting_source_inputs_used=True,
            original_infinite_heat_tail_and_epsilon_atoms_retained=True)
        return packet


def exact_kernel_and_parameter_source_bindings():
    """Identify exact integral definitions before differing cell enclosures."""
    sigma={
        'raw_future_sigma':assignment_source_bindings('outer_initial','point',{
            'a':'c.exp(-1/x**2)','b':'c.exp(-1/(1-x)**2)','value':'a/(a+b)'}),
        'transition_sigma':assignment_source_bindings('axial_pulse_field','point',{
            'odds':'-1/y**2+1/(1-y)**2','e':'c.exp(odds)'})}
    tree=ast.parse((HERE/(PREFIX+'axial_pulse_field.py')).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='sigma_enclosure')
    point=next(n for n in fn.body if isinstance(n,ast.FunctionDef) and n.name=='point')
    if ast.dump(point.body[-1].value)!=ast.dump(ast.parse('e/(1+e)',mode='eval').body):
        raise ValueError('Actual transition logistic source differs')
    x=s.symbols('x',positive=True)
    a=s.exp(-1/x**2);b=s.exp(-1/(1-x)**2);e=a/b
    if s.simplify(a/(a+b)-e/(1+e))!=0:raise ValueError('Exact original sigma source identity differs')
    future=assignment_source_bindings('future_swirl_energy','steep_energy_kernels',{
        'sig':'stable_sigma(c,t)[0]','nextJ':'J+(b-a)*sig'})
    transition=assignment_source_bindings('steep_waiting_C4','transition_kernels',{
        'rate':"1-mu if kind=='in' else 1-delta/2"})
    parsed=ast.parse((HERE/(PREFIX+'steep_waiting_C4.py')).read_text(encoding='utf8'))
    fn=next(n for n in parsed.body if isinstance(n,ast.FunctionDef) and n.name=='transition_kernels')
    wanted=ast.dump(ast.parse('J+ds*sigma_enclosure(c,v)',mode='eval').body)
    if sum(isinstance(n,ast.Assign) and any(ast.unparse(v)=='nextJ' for v in n.targets)
            and ast.dump(n.value)==wanted for n in ast.walk(fn))!=2:
        raise ValueError('Both actual forward/backward J source integrals required')
    transition['both_forward_and_backward_J_source_updates']=True
    expected={
        'future_swirl_energy':("(b-a)*c.exp(-2*mu*t-2*(1-mu)*jc)","(b-a)*c.exp(-2*t+2*(1-delta/2)*jc)"),
        'steep_waiting_C4':('ds*c.exp(-2*mu*v-2*rate*jc)','ds*c.exp(-2*v+2*rate*jc)')}
    for stem,expressions in expected.items():
        parsed=ast.parse((HERE/(PREFIX+stem+'.py')).read_text(encoding='utf8'))
        for expression in expressions:
            wanted=ast.dump(ast.parse(expression,mode='eval').body)
            if sum(isinstance(n,ast.AugAssign) and ast.dump(n.value)==wanted for n in ast.walk(parsed))!=1:
                raise ValueError('Actual exact energy kernel integrand changed: '+stem)
    params=class_assignment('future_swirl_energy','CompliantFutureSwirlEnergy','__init__',
        'self.kernels','steep_energy_kernels(c,self.mu,self.delta,cells)')
    # The actual pinned parameter constructor overrides legacy epsilon.
    name='lei_ren_part1_paper_logarithmic_outer_parameters.py'
    parsed=ast.parse((HERE/name).read_text(encoding='utf8'))
    expected_Ts='4*(c.ln(2)-self.log_delta)'
    if sum(isinstance(n,ast.Assign) and any(ast.unparse(v)=='self.Ts' for v in n.targets)
            and ast.dump(n.value)==ast.dump(ast.parse(expected_Ts,mode='eval').body) for n in ast.walk(parsed))!=1:
        raise ValueError('Actual original Ts definition changed')
    epsilon=class_assignment('pressure_source','CompliantOuterParameters','__init__',
        'self.log_epsilon',"c.ln(c.mpf('.001'))+self.log_delta")
    # A single exact integral is shared by the separately constructed boxes.
    # Replay both PRODUCTION integrands before introducing its abstract name.
    v,t,mu,delta,dt=s.symbols('v t mu delta dt',real=True)
    J=s.Function('J')(v);rate=1-mu;k=1-delta/2
    def formal(node,env):
        label=ast.unparse(node)
        if label in env:return env[label]
        if isinstance(node,ast.Constant):return s.Rational(str(node.value))
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -formal(node.operand,env)
        if isinstance(node,ast.BinOp):
            a,b=formal(node.left,env),formal(node.right,env)
            if isinstance(node.op,ast.Add):return a+b
            if isinstance(node.op,ast.Sub):return a-b
            if isinstance(node.op,ast.Mult):return a*b
            if isinstance(node.op,ast.Div):return a/b
            if isinstance(node.op,ast.Pow):return a**b
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='c.exp':return s.exp(formal(node.args[0],env))
        raise ValueError('Unbound actual kernel integrand: '+label)
    transition_tree=ast.parse((HERE/(PREFIX+'steep_waiting_C4.py')).read_text(encoding='utf8'))
    transition_fn=next(n for n in transition_tree.body if isinstance(n,ast.FunctionDef) and n.name=='transition_kernels')
    actual_transition=[n.value for n in ast.walk(transition_fn) if isinstance(n,ast.AugAssign)
        and ast.unparse(n.target)=='remaining']
    if len(actual_transition)!=2:raise ValueError('Both original backward integral definitions required')
    common={'mu':mu,'delta':delta,'ds':dt,'v':v,'jc':J}
    integrands=[formal(node,dict(common,rate=r))/dt for node,r in zip(actual_transition,(rate,k))]
    future_env={'mu':mu,'delta':delta,'t':v,'jc':J,'b - a':dt}
    future_integrands=[assignment('compliant_future_swirl_energy','steep_energy_kernels',target,future_env,augmented=True)/dt
        for target in ('inside','outside')]
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ValueError('Shared exact kernel source differs: '+name)
        proofs[name]=True
    exact_in=s.exp(-2*mu*v-2*rate*J);exact_out=s.exp(-2*v+2*k*J)
    definitions={}
    for label,actual_future,actual_partial,exact in zip(('Ein','Eout'),future_integrands,integrands,(exact_in,exact_out)):
        zero(label+'_actual_future_integrand',actual_future-exact)
        zero(label+'_actual_transition_integrand',actual_partial-exact)
        partial=s.Integral(actual_partial,(v,t,1));full=s.Integral(actual_future,(v,0,1))
        zero(label+'_same_full_source_at_t0',full-partial.subs(t,0))
        zero(label+'_exact_backward_FTC',s.diff(partial,t)+exact.subs(v,t))
        definitions[label]=str(partial)
    sigma_exact=a/(a+b)
    zero('original_sigma_reflection_identity',sigma_exact+sigma_exact.subs(x,1-x)-1)
    # Reflection identifies 2 int_0^1 sigma = int_0^1 1, hence J(1)=1/2.
    zero('original_J1_half_by_reflected_exact_integral',s.integrate(s.Integer(1),(x,0,1))/2-s.Rational(1,2))
    forward={}
    for target,exact_pair in (
            ('angular',(s.exp(rate*(v-J)),s.exp(k*J))),
            ('pressure',(s.exp(-(1+2*mu)*v-2*rate*J)/2,s.exp(-3*v+2*k*J)/2))):
        nodes=[n.value for n in ast.walk(transition_fn) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)==target]
        if len(nodes)!=2:raise ValueError('Both actual forward integral definitions required: '+target)
        for label,node,r,exact in zip(('in','out'),nodes,(rate,k),exact_pair):
            actual=formal(node,dict(common,rate=r))/dt
            zero(label+'_'+target+'_actual_integrand',actual-exact)
            integral=s.Integral(actual,(v,0,t));forward[label+'_'+target]=str(integral)
            zero(label+'_'+target+'_empty_prefix',integral.subs(t,0))
    current_kernel_inputs={target:class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',target,expression)
        for target,expression in {
            'self.inenergy':"self.kernels(0,'in')['remaining_energy']",
            'self.outenergy':"self.kernels(0,'out')['remaining_energy']",
            'self.infull':"self.kernels(1,'in')",'self.outfull':"self.kernels(1,'out')"}.items()}
    return dict(actual_raw_future_and_transition_sigma_bindings=sigma,
        same_exact_original_logistic_sigma_function=True,
        actual_raw_future_and_transition_integral_assignments=dict(future=future,transition=transition),
        actual_raw_future_full_kernel_constructor=params,
        actual_same_energy_integrands_bound=expected,
        exact_original_Ts_source_bound=expected_Ts,actual_compliant_epsilon_source_bound=epsilon,
        shared_exact_kernel_definitions=dict(J='Integral(original_sigma(v),(v,0,t)); J(0)=0; J(1)=1/2',
            backward=definitions,forward=forward),
        actual_shared_kernel_functional_identities=proofs,
        actual_current_full_kernel_inputs=current_kernel_inputs,
        exact_full_kernel_semantic_bindings={
            "future.kernels['steep_in']":'Ein(0)',"future.kernels['steep_out']":'Eout(0)',
            'self.inenergy':'Ein(0)','self.outenergy':'Eout(0)',
            "self.infull['angular']":'Iin(1)',"self.infull['pressure']":'Pin(1)',
            "self.outfull['angular']":'Iout(1)',"self.outfull['pressure']":'Pout(1)'},
        independent_cell_enclosures_do_not_define_different_integrals=True,
        cells_not_required_numerically_equal=True,passed=True,
        input_hashes={name:sha(name),PREFIX+'pressure_source.py':sha(PREFIX+'pressure_source.py')})


def actual_angular_zero_future_bindings():
    bindings=assignment_source_bindings('power_angular_C4','angular',{
        'local':'s-center','ell':"c.mpf('.15')",'beta':'self.flat.beta(local)',
        'zero':"data['post']*0",'h':'[zero for _ in range(5)]','futureE':'zero',
        'F':'[h[0]+1]+h[1:]',
        'X':'(Xbase+pastA*c.exp(-self.rate*s))/F[0]'})
    tree=ast.parse((HERE/(PREFIX+'power_angular_C4.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CompliantPowerAngularC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='angular')
    loop=next(n for n in fn.body if isinstance(n,ast.For))
    if ast.dump(loop.iter)!=ast.dump(ast.parse("zip(data['coeff'],(-3,-1))",mode='eval').body):
        raise ValueError('Actual angular support centers changed')
    branch=next(n for n in ast.walk(loop) if isinstance(n,ast.If)
        and ast.unparse(n.test)=='lo >= endpoints(ell)[1]')
    assignments={ast.unparse(t):ast.unparse(n.value) for n in branch.body if isinstance(n,ast.Assign) for t in n.targets}
    if assignments!={'past':'self.weights','future':'{k: c.mpf(0) for k in self.weights}'}:
        raise ValueError('Actual angular terminal must retain past and empty future')
    for target,expression in {'h[k]':'dj*(beta[k]*math.factorial(k))',
            'futureE':"(dj*(2*future['E'])+dj*dj*future['F'])*c.exp(-2*self.mu*center)"}.items():
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(isinstance(n,ast.AugAssign) and ast.unparse(n.target)==target
                and ast.dump(n.value)==wanted for n in ast.walk(loop))!=1:
            raise ValueError('Actual empty-support source accumulation changed: '+target)
    wanted=ast.dump(ast.parse('pastP*(self.flatten.Ev2*c.exp(-self.prate*(100+self.Lrel))/4)',mode='eval').body)
    if sum(isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='pressure'
            and ast.dump(n.value)==wanted for n in ast.walk(fn))!=1:
        raise ValueError('Actual past absolute pressure history must be retained')
    beta_tree=ast.parse((HERE/(PREFIX+'flat_pulse_derivatives.py')).read_text(encoding='utf8'))
    beta=next(n for n in beta_tree.body if isinstance(n,ast.FunctionDef) and n.name=='beta_jets')
    branch=next(n for n in beta.body if isinstance(n,ast.If))
    if (ast.unparse(branch.test)!='hi <= -1 or lo >= 1'
            or len(branch.body)!=1 or not isinstance(branch.body[0],ast.Return)
            or ast.unparse(branch.body[0].value)!='IntervalTaylor.constant(c, 0, 4)'):
        raise ValueError('Original beta source must vanish off its actual support')
    scaling=assignment_source_bindings('flat_pulse_derivatives','beta',{
        'ell':"c.mpf('.15')",'raw':'beta_jets(c,c.mpf(s)/ell)'})
    if not all(-center>s.Rational(15,100) and (-center)/s.Rational(15,100)>1 for center in (-3,-1)):
        raise ValueError('Angular s0 must lie after both actual supports')
    return dict(actual_terminal_future_assignments=bindings,
        actual_empty_future_and_retained_past_branch=assignments,actual_beta_coordinate_scaling=scaling,
        actual_beta_empty_support_return_bound=True,
        angular_s0_future_energy_exact_zero=True,angular_s0_F0_exact_one=True,
        angular_s0_nonzero_past_angular_and_pressure_histories_not_reset=True,passed=True)


def current_source_bindings(owner):
    p=owner.steep;f=p.future;outer=owner.before.outer
    graph=dict(checked_current_power_angular=owner.before.acceptance_loaded,
        same_current_outer=p.outer is outer,same_actual_C4_C5_prefix_future=f is outer.fifth.fourth.energy.base,
        same_current_native_context=p.ctx is outer.ctx,same_current_mu=p.mu is outer.mu,
        same_current_delta=p.delta is outer.delta,
        original_kernel_callable=p.kernels.__func__ is BASE.kernels,
        original_entry_callable=p.steep_in.__func__ is BASE.steep_in,
        original_power_callable=p.steep_power.__func__ is BASE.steep_power,
        original_exit_callable=p.steep_out.__func__ is BASE.steep_out,
        original_waiting_callable=p.waiting.__func__ is BASE.waiting,
        same_exact_native_prefix_parameter_source=owner.before.bindings['current_native_parameter_source_bridge']['passed'],
        actual_current_angular_Gamma_prefix=owner.before.bindings['current_prefix_and_future_decomposition']['passed'])
    if not all(graph.values()):raise ValueError('Current steep/waiting graph or original algorithms differ')
    ctor={target:class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',target,expression)
        for target,expression in {'self.outer':'source.outer','self.future':'self.outer.future',
            'self.rate':'1-self.mu','self.k':'1-self.delta/2',
            'self.bp':"c.mpf('.5')+self.mu",'self.bh':"c.mpf('.5')+self.delta/2",
            'self.wait':'box(f.angular.waiting)','self.logone':'box(f.angular.waiting_logone)',
            'self.Ts':'box(f.params.Ts)','self.epsilon':'box(f.heat.epsilon)',
            'self.S':'box(f.repair.strong_S_cap)','self.preheat_scalar':'sum(self.atoms.values(),c.mpf(0))',
            'self.tail_normalization':'c.exp(-2*self.logone)',
            'self.atoms':"dict(baseline=1/box(f.delta),epsilon_atom=-2*self.epsilon*box(atoms['W']),epsilon_squared_atom=self.epsilon**2*box(atoms['W_squared']))",
            'self.thetaR':'c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2',
            'self.thetaS':'self.thetaR*c.exp(-self.bp-self.rate/2)',
            'self.thetaQ':"self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)",
            'self.thetaT':"self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)"}.items()}
    live=assignment_source_bindings('current_steep_waiting_source','data',{
        'terminal':'self.outer.angular(Z,0)','source':'self.outer.fifth.angular(Z)',
        'heat':"copy_jet(c,source['scaled_Gamma_future_defect_Taylor']['energy'])",
        'XR':"terminal['angular_Taylor']",'PR':"terminal['pressure_over_Pstar_squared_Taylor']"})
    identities={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ValueError('Current steep source factor differs: '+name)
        identities[name]=True
    d,eps,A,B,S,lone,W,T,mu,hat=s.symbols('d eps A B S lone W T mu hat',real=True)
    v=s.symbols('v',real=True);J=s.Function('J')(v)
    Kin=s.Integral(s.exp(-2*mu*v-2*(1-mu)*J),(v,0,1))
    Kout=s.Integral(s.exp(-2*v+2*(1-d/2)*J),(v,0,1))
    pre=1/d-2*eps*A+eps**2*B;H=(pre-d*S*hat/2)*s.exp(-2*lone)
    zero('current_complete_Gamma_heat_source',assignment('compliant_current_steep_waiting_source','data','preheat_tail',{
        'IntervalTaylor.constant(c, self.preheat_scalar, 5)':pre,'heat':hat,
        'c.mpf([0, endpoints(self.delta * self.S / 2)[1]])':d*S/2})-(pre-d*S*hat/2))
    zero('current_original_epsilon_normalization',assignment('compliant_current_steep_waiting_source','data','H',{
        'preheat_tail':pre-d*S*hat/2,'self.tail_normalization':s.exp(-2*lone)})-H)
    D=lambda a,b:(1-s.exp(-a*b))/a
    qwait=H*s.exp(-d*W)+D(d,W);qpower=qwait*s.exp(-1-d/2)+Kout
    qentry=(qpower*s.exp(-2*T)+D(2,T))*s.exp(-1-mu)
    env={'H':H,'waiting_future':qwait,'after_power':qpower,'self.delta':d,'self.wait':W,
        'self.Ts':T,'self.mu':mu,'self.outenergy':Kout,
        'decay_integral(c, self.delta, self.wait)':D(d,W),'decay_integral(c, 2, self.Ts)':D(2,T)}
    for target,value in (('waiting_future',qwait),('after_power',qpower),('after_entry',qentry)):
        zero('current_'+target+'_same_complete_backward_source',assignment('compliant_current_steep_waiting_source','data',target,env)-value)
    Ns=s.exp(-1-mu);Nq=Ns*s.exp(-2*T);Nt=Nq*s.exp(-1-d/2)
    post=Kin+Ns*D(2,T)+Nq*Kout+Nt*D(d,W)+Nt*s.exp(-d*W)*H
    tail=Nt*s.exp(-d*W)*s.exp(-2*lone)
    actual_future={'self.mu':mu,'self.delta':d,'self.params.Ts':T,'self.angular.waiting':W,
        'self.angular.waiting_logone':lone,'self.Ns':Ns,'self.Nq':Nq,'self.Nt':Nt,
        'self.Ntail':Nt*s.exp(-d*W),'decay_integral(c, c.mpf(2), self.params.Ts)':D(2,T),
        'decay_integral(c, self.delta, self.angular.waiting)':D(d,W)}
    for target,value in (('self.Ns',Ns),('self.Nq',Nq),('self.Nt',Nt),
            ('self.Ntail',Nt*s.exp(-d*W)),('self.tail_multiplier',tail),
            ('self.steep_power',D(2,T)),('self.waiting_energy',D(d,W))):
        zero('actual_raw_future_'+target[5:]+'_source',assignment('compliant_future_swirl_energy','__init__',target,actual_future)-value)
    actual_scalar=assignment('compliant_current_power_angular_source','__init__','self.post_scalar',{
        "box(future.kernels['steep_in'])":Kin,'box(future.Ns)':Ns,'box(future.steep_power)':D(2,T),
        'box(future.Nq)':Nq,"box(future.kernels['steep_out'])":Kout,'box(future.Nt)':Nt,
        'box(future.waiting_energy)':D(d,W),'box(future.tail_multiplier)':tail,'box(future.delta)':d,
        'box(future.heat.epsilon)':eps,"box(atoms['W'])":A,"box(atoms['W_squared'])":B})
    actual_post=assignment('compliant_current_power_angular_source','data','post',{
        'IntervalTaylor.constant(c, self.post_scalar, 5)':actual_scalar,'heat':hat,
        'c.mpf([0, endpoints(self.heatcap)[1]])':tail*d*S/2})
    zero('actual_current_angular_post_is_shared_exact_future',actual_post-post)
    endpoint=actual_angular_zero_future_bindings()
    zero('actual_current_angular_s0_energy_is_current_post_half',assignment('compliant_power_angular_C4','angular','energy',{
        "data['post']":actual_post,'futureE':s.Integer(0),'self.mu':mu,'s':s.Integer(0),
        'F[0]':s.Integer(1),'decay_integral(c, 2 * self.mu, -s)':s.Integer(0)})-actual_post/2)
    zero('current_angular_terminal_to_steep_inlet_complete_future',post-(Kin+qentry))
    zero('actual_steep_entry_t0_energy_is_current_angular_post_half',assignment('compliant_steep_waiting_C4','steep_in','energy',{
        "data['after_entry']":qentry,"kernels['remaining_energy']":Kin,'self.mu':mu,
        'self.rate':1-mu,'t':s.Integer(0),'J':s.Integer(0)})-actual_post/2)
    # The original live data formulas are copied verbatim. Bind every forward
    # angular/absolute-pressure history, not only a serialized terminal value.
    histories={target:expression for target,expression in {
        'XS':"(XR+self.infull['angular'])*c.exp(-self.rate/2)",'XQ':'XS+self.Ts',
        'XT':"(XQ+self.outfull['angular'])*c.exp(-self.k/2)",
        'PS':"PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)",
        'PQ':'PS+decay_integral(c,3,self.Ts)*(self.outer.flatten.Ev2*self.thetaS**2/2)',
        'PT':"PQ+self.outfull['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)"}.items()}
    history_bindings={stem:assignment_source_bindings(stem,'data',histories)
        for stem in ('current_steep_waiting_source','steep_waiting_C4')}
    return dict(current_defining_object_graph=graph,actual_current_constructor_source_bindings=ctor,
        exact_current_kernel_and_parameter_source_bindings=exact_kernel_and_parameter_source_bindings(),
        actual_current_angular_s0_future_support_bindings=endpoint,
        actual_current_live_angular_Gamma_history_bindings=live,
        current_complete_future_factor_identities=identities,actual_current_and_original_history_assignments=history_bindings,
        current_waiting_root_and_heat_atoms_not_read_from_saved_packets=True,
        current_absolute_pressure_and_nonzero_angular_histories_retained=True,
        same_exact_transition_kernel_definitions_with_independent_directed_cells=True,
        source_caps_used_as_defining_field_values=False,interval_overlap_used_as_join_proof=False,passed=True)


class CurrentSteepWaitingSourceAssembly:
    @source_precision
    def __init__(self,require_checked=True):
        self.before=CurrentPowerAngularSourceAssembly()
        self.family=self.before.family;self.source=self.before.source;self.datum_sha=self.before.datum_sha
        self.steep=CurrentSteepWaitingC4(self.before);self.hashes=dict(self.steep.hashes)
        self.registry=dict(self.before.registry)
        theorem=accepted(THEOREM,self.family,self.source,'angular_steep_and_internal_joins_certified')
        self.functional_proof=functional_source_identities()
        if self.functional_proof!=theorem['functional_production_source_identities'] or not all(self.functional_proof.values()):
            raise ValueError('Canonical original steep/waiting source theorem differs')
        for name,digest in theorem['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Current steep source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[THEOREM]=sha(THEOREM);self.bindings=current_source_bindings(self)
        self.hashes.update(self.bindings['exact_current_kernel_and_parameter_source_bindings']['input_hashes'])
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.registry.update({chart:dict(provider=PREFIX+'current_steep_waiting_source.CurrentSteepWaitingC4',
            method=METHODS[chart],coverage_coordinate='t=logR-logRorigin' if chart in ('steep_entry','steep_exit') else 'phase=offset/original_length',
            domain='[0,1]',acceptance_receipt=RECEIPT) for chart in NEW_CHARTS})
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATE)
            if receipt['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current O7 datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def provider(self,chart):return self.steep if chart in NEW_CHARTS else self.before.provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in NEW_CHARTS:return self.before.evaluate(chart,Z,coordinate)
        packet=getattr(self.steep,METHODS[chart])(Z,coordinate)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_provider=self.registry[chart]['provider'],
            acceptance_receipt=RECEIPT,source_packet=packet,
            physical_mixed_grids={'physical_mixed_derivatives_total_order_le4':packet['physical_mixed_derivatives_total_order_le4']},
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
            current_steep_waiting_physical_owner_installed=False,
            output_kind='current steep/waiting derivative source enclosures; no production point selection',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_downstream_chart_owner_count=27,
            ordered_current_chart_registry=self.registry,current_steep_waiting_source_bindings=self.bindings,
            retained_actual_canonical_functional_identities=self.functional_proof,
            retained_current_power_angular_source_certified=True,
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},full_pulse_C4_installed=False,
            current_steep_waiting_physical_owner_installed=False,all_profile_source_charts_callable=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest();packets={}
        for chart in NEW_CHARTS:
            packets[chart]={}
            for name,value in (('inlet',0),('whole_domain',[0,1]),('exit',1)):
                packets[chart][name]=self.evaluate(chart,[-1,1],value)
                print('Current steep/waiting source: '+chart+' '+name,flush=True)
        result.update(whole_current_steep_waiting_source_maps=packets,
            current_angular_terminal=self.before.evaluate('outer_angular',[-1,1],0),input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentSteepWaitingSourceAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current steep/waiting generated: twenty-seven profile source owners, same complete heat future',flush=True)
    return result


if __name__=='__main__':run()
