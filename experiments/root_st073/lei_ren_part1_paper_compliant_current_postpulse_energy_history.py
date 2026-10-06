"""Current selected energy transported through every post-pulse source owner.

This constructs native algorithms on the new complete future/selection
graph. Integral additivity and original integrating factors identify their
backward energy with the selected forward cumulative moment. No physical
assembler, global cone, flat remainder or time recursion is admitted here.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_selected_energy_source import (
    CurrentSelectedEnergySource,HERE,PREFIX,sha,binding,function,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import CurrentFlattenMixedC4
from lei_ren_part1_paper_compliant_current_power_angular_source import (
    CurrentPowerAngularC4,exact_prefix_and_decomposition_bindings)
from lei_ren_part1_paper_compliant_current_steep_waiting_source import (
    CurrentSteepWaitingC4,exact_kernel_and_parameter_source_bindings,actual_angular_zero_future_bindings)
from lei_ren_part1_paper_compliant_current_heat_source import CurrentCollarGammaC4,same_exact_Gamma_future_bindings
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import CompliantFlattenMixedC4
from lei_ren_part1_paper_compliant_power_angular_C4 import CompliantPowerAngularC4
from lei_ren_part1_paper_compliant_steep_waiting_C4 import CompliantSteepWaitingC4
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4
from lei_ren_part1_paper_compliant_outer_angular_repair import intersect
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_postpulse_energy_history.json'
RECEIPT=PREFIX+'current_postpulse_energy_history_check.json'
GATES=('current_complete_future_postpulse_source_chain_owned',
       'current_selected_forward_energy_equals_full_future_history',
       'current_postpulse_energy_and_zero_meridional_histories_certified')
OPEN=('current_exact_repair_installed_in_all_physical_charts','heat_exterior_stress_identity_certified',
      'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
      'full_background_NS_validation','physical_energy_integral_certified',
      'independently_bounded_flat_remainder','full_cartesian_vector_derivatives_certified','temporal_recursion')
METHODS={'flatten':('flatten','flatten'),'outer_power':('outer','power'),
    'outer_angular':('outer','angular'),'steep_in':('steep','steep_in'),
    'steep_power':('steep','steep_power'),'steep_out':('steep','steep_out'),
    'waiting':('steep','waiting'),'heat_collar':('heat','collar'),'heat_exterior':('heat','exterior')}
VIEWS={'flatten_exit':('flatten',[-1,1],100),'power_inlet':('outer_power',[-1,1],0),
    'angular_terminal':('outer_angular',[-1,1],0),'steep_entry':('steep_in',[-1,1],0),
    'steep_power_exit':('steep_power','.613',1),'steep_exit':('steep_out','.613',1),
    'waiting_terminal':('waiting',[-1,1],1),'collar_inlet':('heat_collar',[-1,1],0),
    'collar_exit':('heat_collar',[-1,1],3),'unbounded_exterior':('heat_exterior',[-1,1],[3,mp.inf]),
    'fresh_exterior':('heat_exterior','.613','4.17')}


def native_energy_density_and_terminal_sources(owner):
    """Bind actual velocity recipes, radial Jacobians and terminal moments.

    Energy units are Rv*Utheta(Rv,Z)^2, not an unnormalized kinetic-energy
    certificate. Every ledger row is the production theta^2 times its
    exact log-radius Jacobian. Existing source bridges identify its actual
    sigma/J, signed beta integrals and complete canonical Gamma integral.
    """
    bindings={};identities={}
    def bind(stem,method,target,value):
        binding(stem,method,target,value);bindings[stem+'.'+method+':'+target]=True
    def keyword(stem,method,target,value):
        found=[n.value for n in ast.walk(function(stem,method)) if isinstance(n,ast.keyword) and n.arg==target]
        if len(found)!=1 or ast.dump(found[0])!=ast.dump(ast.parse(value,mode='eval').body):
            raise ValueError('Actual source publication changed: '+stem+'.'+method+':'+target)
        bindings[stem+'.'+method+':'+target]=True
    def zero(label,left,right=0):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Actual source density differs: '+label)
        identities[label]=True
    specs={
        ('compliant_flatten_mixed_C4','flatten'):{'Fc':'(rho*sigma_jets(c,v/100)[0]).exp()',
            'ev':'self.future_energy(Z)'},
        ('compliant_power_angular_C4','power'):{'theta':'one*(c.exp(-100*self.bp-self.bp*y)/2)'},
        ('compliant_power_angular_C4','angular'):{'theta':'F[0]*(c.exp(-100*self.bp-self.bp*self.Lrel-self.bp*s)/2)'},
        ('compliant_steep_waiting_C4','steep_in'):{'theta':'one*self.thetaR*c.exp(-self.bp*t-self.rate*J)'},
        ('compliant_steep_waiting_C4','steep_power'):{'theta':"one*self.thetaS*c.exp(-c.mpf('1.5')*t)"},
        ('compliant_steep_waiting_C4','steep_out'):{'theta':"one*self.thetaQ*c.exp(-c.mpf('1.5')*t+self.k*J)"},
        ('compliant_steep_waiting_C4','waiting'):{'theta':'one*self.thetaT*c.exp(-self.bh*t)'},
        ('compliant_collar_Gamma_C4','packet'):{'theta':'K[0]*(self.theta_base*c.exp(-self.bh*t))'},
        ('compliant_current_steep_waiting_source','CurrentSteepWaitingC4.__init__'):{
            'self.thetaR':'c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2',
            'self.thetaS':'self.thetaR*c.exp(-self.bp-self.rate/2)',
            'self.thetaQ':"self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)",
            'self.thetaT':"self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)"},
        ('compliant_current_heat_source','CurrentCollarGammaC4.__init__'):{
            'self.theta_base':'self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)'},
        ('compliant_current_pulse_flatten_source','terminal_histories'):{'point':'self.pulse.end(Z,0)'},
        ('compliant_axial_pulse_field','end'):{'mhat':'[zero,zero]','m':'[v*self.Ecap for v in mhat]',
            'B':'Bhat*self.Ecap','future':"self.selection.future.future(Z)['complete_future_energy_Taylor']/2",
            'ell':"c.mpf('.15')",'local':'s-center',
            'w':'backward_bump_weights(c,self.mu,normal,local,cells)'},
        ('compliant_axial_pulse_field','backward_bump_weights'):{
            'begin':'c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])',
            'length':'1-begin','out':'[c.mpf(0),c.mpf(0),c.mpf(0)]'},
        ('compliant_flatten_mixed_C4','flatten_mixed'):{
            'rows':'{UZ:[zero for _ in range(5)],UT:theta_rows,UR:[zero for _ in range(5)],P:prows}'},
        ('compliant_power_angular_C4','_packet'):{
            'mixed':'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)'},
        ('compliant_current_power_angular_source','_packet'):{'packet':'BASE._packet(self,*args,**kwargs)'},
        ('compliant_steep_waiting_C4','packet'):{
            'result':"self.outer._packet(Z,theta,X,energy,pressure,[one*r for r in rates],coordinate,dict(stage=stage,entire_steep_waiting_high_mixed_derivatives_available=True,same_actual_angular_terminal_histories_retained=True,exact_Gamma_and_both_epsilon_atoms_retained=True,exact_relative_velocity_log_parts=logs,no_forward_subtraction_of_unrelated_long_future_energy=True))"}}
    for (stem,method),rows in specs.items():
        for target,value in rows.items():
            if method.endswith('.__init__'):
                class_assignment(stem.removeprefix('compliant_'),method.split('.')[0],'__init__',target,value)
                bindings[stem+'.'+method+':'+target]=True
            else:bind(stem,method,target,value)
    theta_values=[n.value for n in ast.walk(function('compliant_flatten_mixed_C4','flatten'))
        if isinstance(n,ast.Assign) and any(ast.unparse(v)=='theta' for v in n.targets)]
    expected=['F/q*c.exp(-self.bp*t)','IntervalTaylor.constant(c,c.exp(-100*self.bp)/2,5)']
    if [ast.dump(v) for v in theta_values]!=[ast.dump(ast.parse(v,mode='eval').body) for v in expected]:
        raise ValueError('Original flatten theta and exact endpoint correlation changed')
    bindings['original_flatten_theta_and_exact_endpoint_correlation']=True
    units="integral_Rv^infinity Utheta_corrected^2 dR /(Rv*Utheta(Rv,Z)^2)"
    keyword('compliant_future_swirl_energy','future','units',repr(units))
    for stem in ('compliant_future_energy_high_jets','compliant_fifth_axial_jets'):
        keyword(stem,'future','units',"old['units']")
    for target,value in {'Mz_over_R_Utheta':'m[0]',
        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared':'m[1]',
        'Uz_over_Utheta':'B','positive_terminal_future_energy_Taylor':'future'}.items():
        keyword('compliant_axial_pulse_field','end',target,value)
    for target in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        keyword('compliant_current_pulse_flatten_source','terminal_histories',target,"point['"+target+"']")
    for target,expression,op in (
        ('mhat[row - 1]',"Cj*(c.exp((c.mpf('.5')-row*self.mu)*(center-s))*w[row-1])",ast.Sub),
        ('Bhat','Cj*beta',ast.Add)):
        found=[n for n in ast.walk(function('compliant_axial_pulse_field','end'))
            if isinstance(n,ast.AugAssign) and ast.unparse(n.target)==target]
        if len(found)!=1 or not isinstance(found[0].op,op) or ast.dump(found[0].value)!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original end meridional moment recurrence changed: '+target)
        bindings['original_end_meridional_recurrence:'+target]=True
    tail_zero=[n for n in ast.walk(function('compliant_axial_pulse_field','backward_bump_weights'))
        if isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse('endpoints(length)==(mp.mpf(0),mp.mpf(0))',mode='eval').body)]
    if len(tail_zero)!=1 or len(tail_zero[0].body)!=1 or ast.unparse(tail_zero[0].body[0])!='return out':
        raise ValueError('Exact empty end support return changed')
    bindings['actual_empty_end_support_returns_initialized_zero_integrals']=True
    radial_returns=[n.value for n in ast.walk(function('compliant_axial_pulse_field','radial')) if isinstance(n,ast.Return)]
    expression='(2*Z*B[0]-(1-self.delta)*Z*m1[0]-d*(m1[1]-2*Z*m1[0]/q))/L'
    if len(radial_returns)!=1 or ast.dump(radial_returns[0])!=ast.dump(ast.parse(expression,mode='eval').body):
        raise ValueError('Original zero meridional radial recovery changed')
    bindings['actual_original_meridional_radial_recovery']=True
    for label,expression in (
        ('m1.append',"Brows[k]-m1[k]*(c.mpf('.5')-mu)"),
        ('m2.append',"Brows[k]-m2[k]*(c.mpf('.5')-2*mu)")):
        calls=[n for n in ast.walk(function('compliant_pulse_mixed_C4','transport_mixed'))
            if isinstance(n,ast.Call) and ast.unparse(n.func)==label]
        if len(calls)!=1 or len(calls[0].args)!=1 or ast.dump(calls[0].args[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original homogeneous moment derivative recurrence changed: '+label)
        bindings['original_homogeneous_moment_ODE:'+label]=True
    # The Gamma packet delegates the pure-swirl differential recurrence to
    # exactly this outer owner, not to an old receipt or disconnected output.
    calls=[n for n in ast.walk(function('compliant_collar_Gamma_C4','packet'))
        if isinstance(n,ast.Call) and ast.unparse(n.func)=='self.outer._packet']
    if len(calls)!=1 or [ast.unparse(v) for v in calls[0].args[:5]]!=['Z','theta','X','energy','pressure']:
        raise ValueError('Current Gamma meridional transport call changed')
    bindings['canonical_Gamma_packet_same_pure_swirl_recurrence']=True
    mu,delta,L,Ts,W,eps,q,t,J,F,K=s.symbols('mu delta L Ts W eps q t J F K',real=True)
    bp=s.Rational(1,2)+mu;rate=1-mu;k=1-delta/2;bh=(1+delta)/2
    thetaR=s.exp(-bp*(100+L))/2;thetaS=thetaR*s.exp(-bp-rate/2)
    thetaQ=thetaS*s.exp(-s.Rational(3,2)*Ts);thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
    thetaBase=thetaT*s.exp(-bh*W)/(1-eps)
    Nf=q*q*s.exp(-200*mu)/4;Nrel=Nf*s.exp(-2*mu*L)
    Ns=s.exp(-1-mu);Nq=Ns*s.exp(-2*Ts);Nt=Nq*s.exp(-1-delta/2)
    # theta is Utheta/Ev0 and Utheta(Rv)=Ev0/q. R/Rv=e^(origin+t).
    ledger={
        'flatten':(F*s.exp(-bp*t)/q,0,s.exp(-2*mu*t)*F*F),
        'base_power':(s.exp(-bp*(100+t))/2,100,Nf*s.exp(-2*mu*t)),
        'corrected_angular':(F*s.exp(-bp*(100+L+t))/2,100+L,Nrel*s.exp(-2*mu*t)*F*F),
        'steep_in':(thetaR*s.exp(-bp*t-rate*J),100+L,Nrel*s.exp(-2*mu*t-2*rate*J)),
        'steep_power':(thetaS*s.exp(-s.Rational(3,2)*t),101+L,Nrel*Ns*s.exp(-2*t)),
        'steep_out':(thetaQ*s.exp(-s.Rational(3,2)*t+k*J),101+L+Ts,Nrel*Nq*s.exp(-2*t+2*k*J)),
        'waiting':(thetaT*s.exp(-bh*t),102+L+Ts,Nrel*Nt*s.exp(-delta*t)),
        'collar_and_full_exterior':(K*thetaBase*s.exp(-bh*t),102+L+Ts+W,
            Nrel*Nt*s.exp(-delta*W)/(1-eps)**2*s.exp(-delta*t)*K*K)}
    Rv,Ev0,r,z=s.symbols('Rv Ev0 radial_variable Z',positive=True)
    actual_Utheta=s.Function('actual_current_Utheta')
    actual_normalized_density=actual_Utheta(r,z)**2/(Rv*actual_Utheta(Rv,z)**2)
    for label,(theta,origin,density) in ledger.items():
        zero('production_'+label+'_squared_velocity_times_radial_Jacobian',q*q*theta*theta*s.exp(origin+t),density)
        radius=Rv*s.exp(origin+t)
        substituted=actual_normalized_density.subs(r,radius).subs({
            actual_Utheta(radius,z):Ev0*theta,actual_Utheta(Rv,z):Ev0/q})
        zero('actual_'+label+'_density_substitution_into_same_cumulative_integral',
            substituted*s.diff(radius,t),density)
    zero('flatten_exact_endpoint_theta',
        (F*s.exp(-bp*t)/q).subs({t:100,F:q/2}),s.exp(-100*bp)/2)
    Rv,Uv=s.symbols('Rv Utheta_Rv',positive=True);R=Rv*s.exp(t)
    zero('actual_log_radius_Jacobian',s.diff(R,t),Rv*s.exp(t))
    zero('actual_flatten_density_in_Rv_Utheta_Rv_squared_units',
        (Uv*F*s.exp(-bp*t))**2*s.diff(R,t)/(Rv*Uv*Uv),s.exp(-2*mu*t)*F*F)
    d1,d2,b1,b2=s.symbols('d1 d2 beta1 beta2',real=True)
    zero('actual_signed_angular_square_with_disjoint_supports',
        s.expand((1+d1*b1+d2*b2)**2-1-2*d1*d2*b1*b2),
        2*d1*b1+d1*d1*b1*b1+2*d2*b2+d2*d2*b2*b2)
    # The original terminal backward supports are strictly before s=0.
    ell=s.Rational(3,20)
    if not all(0-center>ell for center in (-3,-1)):raise ArithmeticError('Terminal beta support not empty')
    identities['actual_selected_end_m1_m2_and_Uz_zero_by_empty_future_supports']=True
    # When Uz=0, dR Mz=Uz and dR Mtheta_z=sqrt(2R)*Utheta*Uz
    # both vanish. Zero values at Rv therefore remain zero on every chart;
    # (3.8)-(3.9) recover Ur=0 from Mz=Mz_Z=Uz=0.
    Utheta,radial=s.symbols('Utheta R',positive=True)
    zero('zero_Mz_FTC_density',s.Integer(0))
    zero('zero_Mtheta_z_FTC_density',s.sqrt(2*radial)*Utheta*0)
    z,dd,ll=s.symbols('Z d L',real=True)
    zero('zero_radial_velocity_from_original_meridional_recovery',
        (2*z*0-(1-delta)*z*0-dd*(0-2*z*0/(1+z*z)))/ll)
    for n in range(4):
        zero('original_Mz_homogeneous_zero_derivative_induction_'+str(n+1),0-0*(s.Rational(1,2)-mu))
        zero('original_Mtheta_z_homogeneous_zero_derivative_induction_'+str(n+1),0-0*(s.Rational(1,2)-2*mu))
    # All derivative rows are derivatives of these identical zero source
    # functions. The production zero rows are linked above only after the
    # actual selected terminal moments have been bound.
    return dict(actual_source_assignments=bindings,identities=identities,
        native_density_ledger={name:str(v[2]) for name,v in ledger.items()},
        actual_complete_integral_units=units,
        actual_selected_Rv_value='complete_future/2 = integral_Rv^infinity Utheta^2 dR /(2*Rv*Utheta(Rv,Z)^2)',
        stage_integrals_identified_by_actual_J_beta_and_full_Gamma_source_bridges=True,
        zero_meridional_terminal_values_and_downstream_FTC_bound=True,passed=True)


def energy_transport_source_proof(owner):
    """Translate native energy recipes after binding their exact integrands."""
    prefix=exact_prefix_and_decomposition_bindings()
    kernels=exact_kernel_and_parameter_source_bindings()
    gamma=same_exact_Gamma_future_bindings()
    supports=actual_angular_zero_future_bindings()
    if not all(v['passed'] for v in (prefix,kernels,gamma,supports)):raise ValueError('Exact source function bridge missing')
    evidence=owner.selected.exact.companion.source_owner.gamma_evidence
    if not evidence['verified']:raise ValueError('Checked full canonical Gamma enclosure/source required')
    actual_sources=native_energy_density_and_terminal_sources(owner)
    # Use the exact functions represented by these recipes. The same-source
    # proof above has already bound each integral/kernel, before any symbols
    # are used for its value. A cap is never their defining source value.
    mu,delta,L,Ts,W,eps,q=s.symbols('mu delta L Ts W eps q',positive=True)
    Ein,Eout,E0,EF,change=s.symbols('exact_Ein exact_Eout exact_heat_E0 exact_flatten_energy exact_signed_angular_change',real=True)
    I=lambda rate,length:(1-s.exp(-rate*length))/rate
    Ns=s.exp(-1-mu);Nq=Ns*s.exp(-2*Ts);Nt=Nq*s.exp(-1-delta/2)
    H=E0/(1-eps)**2;waiting=H*s.exp(-delta*W)+I(delta,W)
    after_power=waiting*s.exp(-1-delta/2)+Eout
    after_entry=(after_power*s.exp(-2*Ts)+I(2,Ts))*Ns
    post=Ein+after_entry
    Nf=q*q*s.exp(-200*mu)/4;Nrel=Nf*s.exp(-2*mu*L)
    total=EF+Nf*I(2*mu,L)+Nrel*(post+change)
    proofs={};bindings={}
    def zero(name,left,right=0):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Current energy transport differs: '+name)
        proofs[name]=True
    def actual(stem,method,target,env,expected):
        value=assignment(stem,method,target,env)
        zero(stem+'.'+method+':'+target,value,expected)
        bindings[stem+'.'+method+':'+target]=True
        return value
    # Fifth full energy is the exact composition just bound by prefix.
    actual('compliant_fifth_axial_jets','future','total',
        {'flatten':EF,'power':Nf*I(2*mu,L),'angular_change':Nrel*change,
         'postrel':Nrel*(post+Nt*s.exp(-delta*W)/(1-eps)**2*(s.Symbol('exact_Gamma_deficit'))),
         'heatdifference':Nrel*Nt*s.exp(-delta*W)/(1-eps)**2*s.Symbol('exact_Gamma_deficit')},total)
    env={'ev':total/2,'Eint':EF,'self.mu':mu,'t':s.Integer(100),'F':q/2}
    flatten_exit=actual('compliant_flatten_mixed_C4','flatten','energy',env,
        (I(2*mu,L)+s.exp(-2*mu*L)*(post+change))/2)
    # The full flatten integrand is the same original source f(t,Z)^2.
    rho,t,v,dt=s.symbols('rho t v dt',real=True);sig=s.Function('sigma')(v/100)
    Fc=s.exp(s.log(q/2)*sig)
    future_density=s.exp(2*sig*s.log(q/2))*s.exp(-2*mu*v)
    zero('original_flatten_same_exact_energy_density',Fc**2*s.exp(-2*mu*v),future_density)
    zero('flatten_closed_cell_weight_is_exact_integral',s.exp(-2*mu*v)*I(2*mu,dt),
        s.Integral(s.exp(-2*mu*t),(t,v,v+dt)).doit())
    D=s.symbols('remaining_power_length',nonnegative=True)
    # Initial power assignment is augmented by its remaining base integral.
    power0=assignment('compliant_power_angular_C4','power','energy',
        {"data['post']":post,"data['full_angular_change']":change,'self.mu':mu,'D':D})
    poweradd=assignment('compliant_power_angular_C4','power','energy',
        {'decay_integral(c, 2 * self.mu, D)':I(2*mu,D)},augmented=True)
    power=power0+poweradd
    zero('flatten_exit_and_power_inlet_exact_function',flatten_exit,power.subs(D,L))
    futureE,F=s.symbols('remaining_signed_angular_energy actual_angular_factor',real=True)
    angular=actual('compliant_power_angular_C4','angular','energy',
        {"data['post']":post,'futureE':futureE,'self.mu':mu,'s':t,
         'decay_integral(c, 2 * self.mu, -s)':I(2*mu,-t),'F[0]':F},
        ((post+futureE)*s.exp(2*mu*t)+I(2*mu,-t))/(2*F*F))
    zero('power_and_angular_inlet_exact_energy',power.subs(D,4),angular.subs({t:-4,F:1,futureE:change}))
    zero('angular_terminal_positive_complete_post_future',angular.subs({t:0,F:1,futureE:0}),post/2)
    env={'preheat_tail':E0,'self.tail_normalization':1/(1-eps)**2,'H':H,
        'self.delta':delta,'self.wait':W,'decay_integral(c, self.delta, self.wait)':I(delta,W),
        'waiting_future':waiting,'self.outenergy':Eout,'after_power':after_power,
        'self.Ts':Ts,'decay_integral(c, 2, self.Ts)':I(2,Ts),'self.mu':mu}
    actual('compliant_current_steep_waiting_source','data','H',env,H)
    actual('compliant_current_steep_waiting_source','data','waiting_future',env,waiting)
    actual('compliant_current_steep_waiting_source','data','after_power',env,after_power)
    actual('compliant_current_steep_waiting_source','data','after_entry',env,after_entry)
    # Native constructor scalar = full entry kernel + transported remainder.
    scalar=assignment('compliant_current_power_angular_source','__init__','self.post_scalar',{
        "box(future.kernels['steep_in'])":Ein,'box(future.Ns)':Ns,'box(future.steep_power)':I(2,Ts),
        'box(future.Nq)':Nq,"box(future.kernels['steep_out'])":Eout,'box(future.Nt)':Nt,
        'box(future.waiting_energy)':I(delta,W),'box(future.tail_multiplier)':Nt*s.exp(-delta*W)/(1-eps)**2,
        'box(future.delta)':delta,'box(future.heat.epsilon)':eps,
        "box(atoms['W'])":s.Symbol('EW'),"box(atoms['W_squared'])":s.Symbol('EW2')})
    deficit=s.Symbol('exact_Gamma_deficit')
    E0exact=1/delta-2*eps*s.Symbol('EW')+eps**2*s.Symbol('EW2')-deficit
    zero('current_post_scalar_minus_same_Gamma_deficit_equals_transport',
        scalar-Nt*s.exp(-delta*W)/(1-eps)**2*deficit,post.subs(E0,E0exact))
    j,remaining=s.symbols('transition_J remaining_exact_transition_energy',real=True)
    ein=actual('compliant_steep_waiting_C4','steep_in','energy',
        {"data['after_entry']":after_entry,"kernels['remaining_energy']":remaining,
         'self.mu':mu,'self.rate':1-mu,'t':t,'J':j},(after_entry+remaining)*s.exp(2*mu*t+2*(1-mu)*j)/2)
    zero('angular_and_steep_in_energy_join',post/2,ein.subs({t:0,j:0,remaining:Ein}))
    epower=actual('compliant_steep_waiting_C4','steep_power','energy',
        {"data['after_power']":after_power,'left':t},(after_power-s.Rational(1,2))*s.exp(-2*t)/2+s.Rational(1,4))
    zero('steep_in_and_steep_power_energy_join',ein.subs({t:1,j:s.Rational(1,2),remaining:0}),epower.subs(t,Ts))
    eout=actual('compliant_steep_waiting_C4','steep_out','energy',
        {"data['waiting_future']":waiting,'self.delta':delta,"kernels['remaining_energy']":remaining,
         't':t,'self.k':1-delta/2,'J':j},(waiting*s.exp(-1-delta/2)+remaining)*s.exp(2*t-2*(1-delta/2)*j)/2)
    zero('steep_power_and_steep_out_energy_join',epower.subs(t,0),eout.subs({t:0,j:0,remaining:Eout}))
    ewait=actual('compliant_steep_waiting_C4','waiting','energy',
        {"data['H']":H,'self.delta':delta,'left':t,'decay_integral(c, self.delta, left)':I(delta,t)},
        (H*s.exp(-delta*t)+I(delta,t))/2)
    zero('steep_out_and_waiting_energy_join',eout.subs({t:1,j:s.Rational(1,2),remaining:0}),ewait.subs(t,W))
    ecollar=actual('compliant_collar_Gamma_C4','collar','energy',
        {"tails['remaining_energy_in_Rtail_units']":E0,'self.delta':delta,'t':t,'K':F},E0*s.exp(delta*t)/(2*F*F))
    zero('waiting_and_collar_energy_join',ewait.subs(t,0),ecollar.subs({t:0,F:1-eps}))
    xi,h=s.symbols('xi current_Gamma_H',real=True);BE=s.symbols('full_Gamma_energy_future',real=True)
    eext=actual('compliant_collar_Gamma_C4','exterior','energy',
        {"local['energy_numerator']":BE,'K':h},BE/(2*h*h))
    zero('collar_exit_and_full_Gamma_exterior_same_energy',
        ecollar.subs({t:3,F:h,E0:s.exp(-3*delta)*BE}),eext)
    # Same remaining integral at any radial point gives the unique forward
    # moment: Me(R)=Me(Rv)-1/2 int_Rv^R Utheta^2 dR =1/2 int_R^infty Utheta^2 dR.
    R,Rv,r,Z=s.symbols('R Rv radial_integration_variable Z',positive=True)
    actual_Utheta=s.Function('actual_current_Utheta')
    density=actual_Utheta(r,Z)**2
    total_integral=s.Integral(density,(r,Rv,s.oo))
    forward=total_integral/2-s.Integral(density,(r,Rv,R))/2
    backward=s.Integral(density,(r,R,s.oo))/2
    zero('actual_forward_and_remaining_energy_same_radial_derivative',s.diff(forward-backward,R))
    zero('actual_forward_and_remaining_energy_same_selected_Rv_value',
        (forward-backward).subs(R,Rv).doit())
    proofs['same_initial_value_and_original_cumulative_FTC_identify_whole_energy_function']=True
    return dict(source_owner_graph=owner.assert_graph(),source_assignments=bindings,identities=proofs,
        actual_velocity_density_Jacobians_and_terminal_moment_source_proof=actual_sources,
        complete_future_decomposition_source_proof=prefix,shared_transition_integral_proof=kernels,
        full_current_Gamma_future_source_proof=gamma,angular_terminal_empty_future_supports=supports,
        checked_canonical_Gamma_source_evidence=evidence,
        exact_heat_energy_source='E0(Z)=Integral_0^infinity exp(-delta*t)*K(t,Z)^2 dt; K=(1-sigma)*(1-epsilon)+sigma*(1-epsilon*phi)*H',
        exact_inverse_radius_source='S=exp(sum(current exact repair logS terms))=1/Rtail; caps enclose S only',
        cumulative_to_native_energy_units='Me(R)=Integral_R^infinity actual_current_Utheta(r,Z)^2 dr/2; native energy=Me(R)/(R*actual_current_Utheta(R,Z)^2)',
        same_current_full_future_half_normalization_on_every_stage=True,
        selected_forward_cumulative_energy_equals_same_remaining_integral_by_FTC=True,
        zero_meridional_histories_propagate_from_selected_terminal_by_FTC=True,
        positive_physical_energy_not_zeroed=True,passed=True)


class CurrentPostpulseEnergyHistory:
    @source_precision
    def __init__(self,selected=None,require_checked=True):
        self.selected=selected if selected is not None else CurrentSelectedEnergySource()
        if not self.selected.acceptance_loaded:raise ValueError('Checked current complete future/selection required')
        self.family=self.selected.family;self.source=self.selected.source;self.datum_sha=self.selected.datum_sha
        self.hashes=dict(self.selected.hashes);self.ctx=self.selected.pulse.ctx
        proxy=SimpleNamespace(pulse=self.selected.pulse,family=self.family,source=self.source,hashes=self.hashes)
        self.flatten=proxy.flatten=CurrentFlattenMixedC4(proxy.pulse,self.family,self.source,self.hashes)
        self.outer=proxy.outer=CurrentPowerAngularC4(proxy)
        self.steep=proxy.steep=CurrentSteepWaitingC4(proxy)
        self.heat=CurrentCollarGammaC4(proxy)
        c=self.ctx;exact=self.selected.exact
        self.heat.exact_logS_terms=dict(exact.repair.heat.logS_terms)
        self.heat.exact_logRtail_terms=dict(exact.repair.heat.logradius_terms)
        source_S=self.selected.pulse.factor(sum((c.mpf(endpoints(v)) for v in self.heat.exact_logS_terms.values()),c.mpf(0)))
        self.heat.S=intersect(c,source_S,c.mpf([0,endpoints(self.heat.Scap)[1]]))
        self.assert_graph();self.proof=energy_transport_source_proof(self)
        for stem in ('current_pulse_flatten_source','current_power_angular_source','current_steep_waiting_source',
            'current_heat_source','flatten_mixed_C4','power_angular_C4','steep_waiting_C4','collar_Gamma_C4'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current postpulse energy source scope/datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def assert_graph(self):
        sel=self.selected;old=sel.exact.companion.heat
        graph=dict(checked_complete_current_selected_owner=sel.acceptance_loaded,
            current_flatten_terminal_and_future_owner=self.flatten.pulse is sel.pulse,
            current_power_flatten_and_pulse_owners=self.outer.flatten is self.flatten and self.outer.pulse is sel.pulse,
            current_complete_fifth_angular_and_selected_owner=self.outer.fifth is sel.fifth,
            current_complete_future_owner=self.outer.future is self.steep.future is self.heat.future is sel.future,
            current_steep_and_heat_outer_owner=self.steep.outer is self.heat.outer is self.outer,
            current_heat_steep_owner=self.heat.steep is self.steep,
            current_heat_repair_and_exact_Gamma_owner=self.heat.repair is sel.exact.repair and self.heat.exact_heat is sel.future.heat,
            common_native_Taylor_context=self.flatten.ctx is self.outer.ctx is self.steep.ctx is self.heat.ctx is self.ctx,
            same_native_incoming_functions=self.flatten.inlet.constants is sel.pulse.high.constants and self.flatten.inlet.datum is sel.future.angular.initial.datum,
            same_exact_inverse_radius_source=self.heat.exact_logS_terms==sel.exact.repair.heat.logS_terms,
            original_flatten_algorithm=self.flatten.flatten.__func__ is CompliantFlattenMixedC4.flatten,
            original_power_angular_algorithms=self.outer.power.__func__ is CompliantPowerAngularC4.power and self.outer.angular.__func__ is CompliantPowerAngularC4.angular,
            original_steep_waiting_algorithms=all(getattr(self.steep,n).__func__ is getattr(CompliantSteepWaitingC4,n) for n in ('steep_in','steep_power','steep_out','waiting')),
            original_collar_and_exterior_algorithms=self.heat.collar.__func__ is CompliantCollarGammaC4.collar and self.heat.exterior.__func__ is CompliantCollarGammaC4.exterior,
            all_postpulse_source_owners_are_new=self.flatten is not old.outer.flatten and self.outer is not old.outer and self.steep is not old.steep and self.heat is not old,
            all_dependent_caches_are_new=all(getattr(selfobj,name) is not getattr(oldobj,name) for selfobj,oldobj,name in (
                (self.outer,old.outer,'cache'),(self.steep,old.steep,'cache'),(self.steep,old.steep,'kernel_cache'),
                (self.heat,old,'cache'),(self.heat,old,'shape_cache'),(self.heat,old,'tail_cache'),(self.heat,old,'gamma_cache'))))
        if not all(graph.values()):raise ValueError('Current postpulse source graph differs: '+str(graph))
        return graph

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        self.assert_graph();key=(chart,tuple(endpoints(self.ctx.mpf(Z))),tuple(endpoints(self.ctx.mpf(coordinate))))
        if key not in self.cache:
            if chart not in METHODS:raise ValueError('Unknown postpulse chart')
            owner,method=METHODS[chart];point=getattr(getattr(self,owner),method)(Z,coordinate)
            self.cache[key]=dict(chart=chart,Z=self.ctx.mpf(Z),coordinate=self.ctx.mpf(coordinate),
                source_packet=point,current_complete_selected_future_owner_consumed=True,
                same_exact_energy_function_across_stage_units=True,source_owner_graph=self.assert_graph(),
                **dict.fromkeys(OPEN,False))
        return self.cache[key]

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_energy_transport_source_proof=self.proof,
            current_postpulse_energy_views={name:self.evaluate(*args) for name,args in VIEWS.items()},
            installation_scope='current postpulse profile source owners; downstream physical assembly still separate',
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPostpulseEnergyHistory(require_checked=False)
    result=field.report();(HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current selected energy -> native flatten/angular/steep/waiting/full Gamma source chain generated',flush=True)
    return result


if __name__=='__main__':run()
