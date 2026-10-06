"""Current angular and heat-tail cone; original full histories stay intact.

Partial integrals are bounded by positivity of their defining integrands.
Waiting baselines cancel algebraically; heat flat factors are never divided
by zero. All bounds refer to the checked current source, not old receipts.
"""
import ast
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_flatten_power_cone import (
    HERE,PREFIX,sha,pack,encode,endpoints,source_precision)
from lei_ren_part1_paper_compliant_current_original_cone import SourceAST
from lei_ren_part1_paper_compliant_current_angular_background_stress import angular_source_log_parts
from lei_ren_part1_paper_compliant_angular_cone import angular_cone_identities
from lei_ren_part1_paper_compliant_steep_entry_cone import steep_entry_cone_identities
from lei_ren_part1_paper_compliant_steep_power_cone import steep_power_cone_identities
from lei_ren_part1_paper_compliant_steep_exit_cone import steep_exit_cone_identities
from lei_ren_part1_paper_compliant_waiting_cone import waiting_cone_identities
from lei_ren_part1_paper_compliant_collar_cone import sigmoid_cone_identities
from lei_ren_part1_paper_compliant_collar_heat_cone import heat_cone_identities
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_defect_rows
from lei_ren_part1_paper_compliant_waiting_stress_C3 import waiting_stress_rows
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_current_heat_source import same_exact_Gamma_future_bindings

DOMAINS=dict(outer_angular=(-4,0),steep_entry=(0,1),steep_power=(0,1),
    steep_exit=(0,1),waiting=(0,1),heat_collar=(0,3),heat_exterior=(3,4))
COORDINATES=dict(outer_angular='angular_offset',steep_entry='entry_t',
    steep_power='original_coordinate',steep_exit='original_coordinate',
    waiting='original_coordinate',heat_collar='heat_offset',heat_exterior='heat_offset')


def generic_angular_tail_theorem():
    result=dict(angular=angular_cone_identities(),entry=steep_entry_cone_identities(),
        power=steep_power_cone_identities(),exit=steep_exit_cone_identities(),
        waiting=waiting_cone_identities(),sigmoid=sigmoid_cone_identities(),heat=heat_cone_identities())
    result['input_hashes']={PREFIX+stem+'.py':sha(PREFIX+stem+'.py') for stem in (
        'angular_cone','steep_entry_cone','steep_power_cone','steep_exit_cone',
        'waiting_cone','collar_cone','collar_heat_cone')}
    result['historical_cone_receipts_loaded_or_promoted']=False
    return result


def current_angular_tail_source_theorem(field):
    field.assert_graph();reg=field.registry;ang=reg.owners['angular'];entry=reg.owners['entry']
    o7=reg.owners['o7'];heat=reg.owners['heat'];internal=reg.owners['angular_internal'];flat=reg.owners['flatten']
    for owner in (ang,entry,o7,heat,internal):
        owner.assert_graph()
        if not owner.acceptance_loaded:raise ValueError('Checked current tensor owners required')
    links=dict(same_flatten_exit=ang.outer.flatten is flat.flatten,
        same_angular_internal=internal.angular is ang,same_entry_angular=entry.angular is ang,
        same_O7_entry=o7.entry_source is entry,same_heat_O7=heat.o7 is o7,
        same_heat_function=heat.heat is ang.heat,same_complete_history=heat.history is ang.history)
    if not all(links.values()) or endpoints(field.bounds['correlated_Hf_uniform_lower'])[0]<=0:
        raise ValueError('Same current positive Hf and complete angular/heat graph required')
    for proof in (ang.moment_proof,ang.baselines,ang.pressure_proof,entry.moment_proof,
            entry.shape_proof,entry.correlation,o7.normalization,o7.moment_proof,
            heat.units,heat.normalization,internal.proof):
        if not proof['passed']:raise ValueError('Current full moment/source proof missing')
    joins=dict(power_angular=flat.joins['power_angular'],angular_entry=entry.join_proof,
        **o7.joins,**heat.joins)
    if not all(v['passed'] for v in joins.values()):raise ValueError('All seven current tensor joins required')
    if not flat.history.proof['zero_meridional_histories_propagate_from_selected_terminal_by_FTC']:
        raise ValueError('Same actual zero axial shear FTC source required')
    asts=SourceAST();checks={}
    def bind(stem,method,target,wanted=None,aug=False):
        return asts.expression(stem,method,target,wanted=wanted,augmented=aug)
    bind('power_angular_C4','angular','h[k]','dj*(beta[k]*math.factorial(k))',True)
    bind('power_angular_C4','angular','F','[h[0]+1]+h[1:]')
    bind('power_angular_C4','angular','pastA',"dj*(c.exp(self.rate*center)*past['A'])",True)
    base=bind('power_angular_C4','angular','Xbase',"eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    native=bind('power_angular_C4','angular','X','(Xbase+pastA*c.exp(-self.rate*s))/F[0]')
    bind('power_angular_C4','angular','y','self.Lrel+s')
    bind('power_angular_C4','angular','past','bump_weights(c,self.mu,self.normalization,local,cells=self.cells)')
    bind('current_power_angular_source','data','f','self.flatten.flatten(Z,100)')
    bind('outer_angular_repair','bump_weights','ell',"c.mpf('.15')")
    bind('outer_angular_repair','bump_weights','beta','raw_beta(c,raw_coordinate)/(ell*normalization)')
    bind('outer_angular_repair','bump_weights','result[name]','ds*c.exp(rate*s)*beta**power',True)
    fn=asts.method('outer_angular_repair','bump_weights')
    rates=next(n.iter for n in ast.walk(fn) if isinstance(n,ast.For) and ast.unparse(n.target)=='(name, rate, power)')
    if ast.dump(rates.elts[0])!=ast.dump(ast.parse("('A',1-mu,1)",mode='eval').body):
        raise ValueError('Original partial A integrand changed')
    bind('outer_pulse_map','raw_beta','lower','c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)')
    bind('outer_pulse_map','raw_beta','upper','c.exp(-1/(1-c.mpf(near)**2)) if near<1 else c.mpf(0)')
    bind('current_angular_background_stress','normalized_full_moment_rows','g','heat.a-outer.mu')
    bind('current_angular_background_stress','normalized_full_moment_rows','fac','c.exp(g*v)')
    bind('current_angular_background_stress','normalized_full_moment_rows','K',
        '[sum((F[n]*math.comb(j,n)*g**(j-n) for n in range(j+1)),one*0)*fac for j in range(5)]')
    bind('current_angular_background_stress','normalized_full_moment_rows','A',"[K[0]*copy_jet(c,native['angular_Taylor'])]")
    bind('current_steep_entry_background_stress','normalized_entry_moment_rows','A',"[K[0]*copy_jet(c,native['angular_Taylor'])]")
    # Replay the native angular X before subtracting any unit baseline.
    a,mu,q,z,L=s.symbols('a mu s Z Lrel',real=True);r=1-mu;k=1-a;b=(1-2*a)/2
    xf=s.Function('same_flatten_Xf')(z);h=s.Function('same_h')(q,z)
    ds=[s.Function('d'+str(j))(z) for j in range(2)]
    ps=[s.Function('partial_A'+str(j))(q) for j in range(2)]
    ctx=SimpleNamespace(exp=s.exp);obj=SimpleNamespace(rate=r)
    past=sum(d*s.exp(r*center)*p for d,center,p in zip(ds,(-3,-1),ps))
    Xbase=asts.evaluate(base,dict(self=obj,data={'flatten_exit_X':xf},eq=1/r,c=ctx,y=L+q))
    X=asts.evaluate(native,dict(c=ctx,self=obj,Xbase=Xbase,pastA=past,s=q,F=[1+h]))
    N=s.expand(X*(1+h));Hf=k*(xf-1/r)-b*z*s.diff(xf,z)
    W=k*N-b*z*s.diff(N,z)-(1+h)
    stable=(mu-a)/r+Hf*s.exp(-r*(L+q))-h
    stable+=sum(s.exp(-r*(q-center))*p*(k*d-b*z*s.diff(d,z)) for d,center,p in zip(ds,(-3,-1),ps))
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Current angular/tail identity: '+name)
        checks[name]=True
    zero('same_native_current_stable_W',W-stable)
    G=s.exp((a-mu)*q);K=G*(1+h);A=K*X
    zero('current_A_KX_inertial_equals_G_W',k*A-b*z*s.diff(A,z)-K-G*W)
    zero('same_current_K_shear',s.diff(K,q)-(1+a)*K-G*(s.diff(h,q)-(1+mu)*(1+h)))
    zero('same_current_kappa_minus2',2*a-2*s.diff(K,q)/K-(2*mu-2*s.diff(h,q)/(1+h)))
    # FTC gives A_partial' = exp(r*u)*beta(u) >= 0. Its complement
    # is the integral of the SAME nonnegative density, hence 0<=A_partial<=A.
    ell=s.Rational(3,20);u,Nbeta=s.symbols('u normalization',real=True,positive=True)
    density=s.exp(r*u)*s.exp(-1/(1-(u/ell)**2))/(ell*Nbeta)
    derivative=s.exp(r*u)*s.Function('same_normalized_beta')(u)
    checks['partial_A_monotone_by_positive_defining_density_and_full_complement']=True
    # Current heat source is C times local K/A, C^2 times full E/P.
    bind('current_heat_background_tensor','chart','logC','LT-st.logone')
    bind('current_heat_background_tensor','chart','C','c.exp(logC)')
    bind('current_steep_waiting_background_stress','chart','logK','KT')
    bind('waiting_stress_C3','waiting_stress_rows','It','Nt/L*c.exp(-heat.k*q)')
    bind('waiting_stress_C3','waiting_stress_rows','St','IntervalTaylor.constant(c,-2*heat.S*(1+heat.a)*K0*c.exp(-q),4)')
    bind('waiting_stress_C3','waiting_stress_rows','Je','Ne/L*c.exp(heat.delta*q)')
    bind('waiting_stress_C3','waiting_stress_rows','Jp','Np/L*c.exp(heat.prate*q)')
    eps,C=s.symbols('eps C',positive=True);delta=2*a;p=1+delta;K0=1-eps
    A0,E0,P0=[s.Function(name)(z) for name in ('Ad0','Ed0','Pd0')];Q0=-2*eps+eps**2
    Af=K0/k+(A0+eps/k)*s.exp(-k*q)
    Ef=K0**2/delta+(E0-Q0/delta)*s.exp(delta*q)
    Pf=K0**2/(2*p)+(P0-Q0/(2*p))*s.exp(p*q)
    for name,value,ode,start in (
        ('A',Af,K0-k*Af,A0+1/k),('E',Ef,delta*Ef-K0**2,E0+1/delta),
        ('P',Pf,p*Pf-K0**2/2,P0+1/(2*p))):
        zero('same_waiting_full_'+name+'_FTC',s.diff(value,q)-ode)
        zero('same_waiting_full_'+name+'_endpoint',value.subs(q,0)-start)
    Ct=(k*Af-b*z*s.diff(Af,z)-K0)/(1-delta*z*z)
    Cz=(delta*z*Ef-(1-z*z)*s.diff(Ef,z)/2-2*p*z*Pf+(1-z*z)*s.diff(Pf,z))/(1-delta*z*z)
    nt=k*A0+eps-b*z*s.diff(A0,z)
    ne=delta*z*E0-z*Q0-(1-z*z)*s.diff(E0,z)/2
    np=-2*p*z*P0+z*Q0+(1-z*z)*s.diff(P0,z)
    zero('actual_waiting_cancelled_full_theta',Ct-nt*s.exp(-k*q)/(1-delta*z*z))
    zero('actual_waiting_cancelled_full_axial',Cz-(ne*s.exp(delta*q)+np*s.exp(p*q))/(1-delta*z*z))
    zero('actual_full_waiting_C_theta_units',(k*C*Af-b*z*s.diff(C*Af,z)-C*K0)/(1-delta*z*z)-C*Ct)
    zero('actual_full_waiting_C_squared_axial_units',
        (delta*z*C**2*Ef-(1-z*z)*s.diff(C**2*Ef,z)/2-2*p*z*C**2*Pf+(1-z*z)*s.diff(C**2*Pf,z))/(1-delta*z*z)-C**2*Cz)
    # Actual current waiting/heat endpoint inputs, before invoking ODE
    # uniqueness. An overlap between numerical endpoint boxes is insufficient.
    shared_future=same_exact_Gamma_future_bindings()
    if not shared_future['passed']:raise ValueError('Current waiting/collar energy integrands differ')
    bind('current_steep_waiting_source','data','H','preheat_tail*self.tail_normalization')
    bind('current_steep_waiting_source','__init__','self.tail_normalization','c.exp(-2*self.logone)')
    bind('collar_Gamma_C4','data','terminal','self.steep.waiting(Z,1)')
    bind('collar_Gamma_C4','data','defect',"Xtail*(1-self.eps)-heat0['angular_numerator']")
    bind('collar_Gamma_C4','collar','X',"(tails['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K")
    bind('current_angular_background_stress','post_pressure','tails','heat.collar_tails(Z,0)')
    bind('current_steep_waiting_background_stress','chart','Ih',"copy_jet(c,post['original_full_Gamma_pressure_future'])")
    bind('current_steep_waiting_background_stress','chart','pressure0',
        '(one*(decay_integral(c,p,left)/2)+Ih*c.exp(-p*left-2*st.logone))*c.exp(2*logK)')
    bind('current_heat_background_tensor','chart','A0',"copy_jet(c,tails['angular_numerator'])*C")
    bind('current_heat_background_tensor','chart','E0',"copy_jet(c,tails['remaining_energy_in_Rtail_units'])*(C*C*c.exp(h.delta*t))")
    bind('current_heat_background_tensor','chart','P0',"copy_jet(c,tails['remaining_pressure_in_Rtail_units'])*(C*C*c.exp(h.prate*t))")
    bind('current_heat_background_tensor','chart','defects','collar_defect_rows(h,shape,tails,h.ctx.mpf(t))')
    bind('current_heat_background_tensor','chart','rawstress','collar_stress_rows(h,shape,defects,h.ctx.mpf(Z),h.ctx.mpf(t))')
    if any(endpoints(v)!=(0,0) for v in sigma_jets(field.ctx,field.ctx.mpf(0))):
        raise ValueError('Exact flat sigmoid heat inlet changed')
    links0=heat.exterior.proof['current_exact_function_links']
    if not all(links0[v] for v in ('checked_exact_angular_function_identity','checked_exact_pressure_function_identity','checked_full_future_half_energy_identity')):
        raise ValueError('Actual current Dtheta/Cp/full energy endpoint identities missing')
    logKT,lone=s.symbols('KT logone',real=True);xt,eh,ph=s.symbols('Xtail heat_full_energy heat_full_pressure',real=True)
    unit=s.exp(logKT-lone);kloc=s.exp(lone)
    zero('same_actual_waiting_heat_endpoint_K',s.exp(logKT)-unit*kloc)
    zero('same_actual_waiting_heat_endpoint_A',s.exp(logKT)*xt-unit*(kloc*xt))
    zero('same_actual_waiting_heat_endpoint_E',s.exp(2*logKT)*eh*s.exp(-2*lone)-unit**2*eh)
    zero('same_actual_waiting_heat_endpoint_P',s.exp(2*logKT-2*lone)*ph-unit**2*ph)
    # Replay the ACTUAL O7 entry-stress operator and waiting-stress operator
    # on arbitrary common full terminal functions and the same exact S.
    stub=SimpleNamespace(variable=lambda c,value,n:value,constant=lambda c,value,n:s.sympify(value))
    ctxt=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    Sin=s.Symbol('same_exact_positive_S',positive=True)
    env=dict(IntervalTaylor=stub,math=math,axial_derivative=lambda value:s.diff(value,z))
    env['shifted_rows']=asts.replay('collar_stress_C3','shifted_rows',env)
    env['product_rows']=lambda x,y:[sum(x[j]*y[n-j]*s.binomial(n,j) for j in range(n+1)) for n in range(min(len(x),len(y)))]
    env['collar_stress_rows']=asts.replay('collar_stress_C3','collar_stress_rows',env)
    entryfun=asts.replay('steep_entry_stress_C3','entry_stress_rows',env)
    waitingfun=asts.replay('waiting_stress_C3','waiting_stress_rows',env)
    collarfun=asts.replay('collar_stress_C3','collar_defect_rows',env)
    # Instantiate the actual shape AST at the flat sigmoid inlet. Arbitrary
    # Gamma/phi jets survive until multiplication by the zero sigmoid rows.
    shapeheat=SimpleNamespace(a=a,S=Sin,eps=eps)
    phi0=s.symbols('phi0:5');gamma0=[s.Function('Gamma_deficit'+str(j))(z) for j in range(5)]
    sr=[s.Integer(0)]*5;unity=[s.Integer(1)]+[s.Integer(0)]*4
    Wnode=bind('collar_Gamma_C4','shape','W','[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]')
    Wrows=asts.evaluate(Wnode,dict(unity=unity,sr=sr,fr=phi0,product_rows=env['product_rows']))
    Crows=asts.evaluate(bind('collar_Gamma_C4','shape','C','[unity[j]-fr[j]*self.eps for j in range(5)]'),dict(unity=unity,fr=phi0,self=shapeheat))
    Drs=asts.evaluate(bind('collar_Gamma_C4','shape','D','product_rows(product_rows(sr,C),dr)'),dict(sr=sr,C=Crows,dr=gamma0,product_rows=env['product_rows']))
    pre=asts.evaluate(bind('collar_Gamma_C4','shape','pre','[unity[j]-W[j]*self.eps for j in range(5)]'),dict(unity=unity,W=Wrows,self=shapeheat))
    Kshape=asts.evaluate(bind('collar_Gamma_C4','shape','K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]'),dict(pre=pre,D=Drs,self=shapeheat))
    for j,row in enumerate(Kshape):zero('actual_heat0_shape_K_y'+str(j),row-(K0 if j==0 else 0))
    future={key:s.Function('same_full_future_'+key)(z) for key in ('theta','energy','pressure')}
    atoms={key:s.Symbol('same_eps_atom_'+key,real=True) for key in ('JW','EW','EW2','PW','PW2')}
    localheat=SimpleNamespace(ctx=ctxt,a=a,eps=eps,delta=delta,k=k,prate=p,S=Sin)
    actualbase=collarfun(localheat,dict(W_rows=Wrows,D_rows=Drs),
        dict(scaled_full_future_Gamma_defects=future,separate_epsilon_atoms=atoms),s.Integer(0))
    tails0={}
    for name in ('A','E','P'):
        node=bind('collar_Gamma_C4','collar_tails',name)
        tails0[name]=asts.evaluate(node,dict(one=s.Integer(1),self=localheat,c=ctxt,t=s.Integer(0),
            angular=future['theta'],energy=future['energy'],pressure=future['pressure'],atoms=atoms))
    for name,key,baseline in (('A','angular_defect_rows',1/k),('E','energy_defect_rows',1/delta),('P','pressure_defect_rows',1/(2*p))):
        zero('actual_heat0_'+name+'_full_future_and_defect_source',actualbase[key][0]+baseline-tails0[name])
    collaredge=env['collar_stress_rows'](localheat,dict(K_rows=Kshape),actualbase,z,s.Integer(0))
    waitingedge=waitingfun(localheat,actualbase,z,s.Integer(0))
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):zero('actual_collar0_waiting0_'+label+'_y%d_Z%d'%(j,n),s.diff(collaredge[label][j]-waitingedge[label][j],z,n))
    heatstub=SimpleNamespace(ctx=ctxt,a=a,delta=delta,k=k,prate=p,eps=eps,S=Sin*s.exp(-q))
    full0=dict(angular_defect_rows=[C*s.diff(Af,q,j) for j in range(5)],
        energy_defect_rows=[C**2*s.diff(Ef,q,j) for j in range(5)],
        pressure_defect_rows=[C**2*s.diff(Pf,q,j) for j in range(5)],
        K_defect_rows=[C*K0]+[s.Integer(0)]*4)
    currentrows=entryfun(heatstub,dict(K_rows=full0['K_defect_rows']),full0,
        [s.diff(Af/K0,q,j) for j in range(5)],z,s.Integer(0))
    localstub=SimpleNamespace(**dict(heatstub.__dict__,S=Sin))
    localrows=waitingfun(localstub,dict(angular_defect_rows=[A0],energy_defect_rows=[E0],pressure_defect_rows=[P0]),z,q)
    for label,scale in (('theta',C),('axial',C**2)):
        for j in range(4):
            residual=s.simplify(currentrows[label][j]-scale*localrows[label][j])
            for n in range(4-j):zero('actual_O7_waiting_operator_'+label+'_y%d_Z%d'%(j,n),s.diff(residual,z,n))
    # The actual sigmoid/Gamma source and flat exponential are bound too.
    for target,wanted in (('sig','sigma_jets(c,t)'),('sr','[one*(sig[j]*math.factorial(j)) for j in range(5)]'),
        ('C','[unity[j]-fr[j]*self.eps for j in range(5)]'),('pre','[unity[j]-W[j]*self.eps for j in range(5)]'),
        ('dr','gamma_deficit_mixed(c,Z,self.a,self.Scap,t,4 if high else 0)'),
        ('K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')):
        bind('collar_Gamma_C4','shape',target,wanted)
    bind('collar_Gamma_C4','phi_jets','dist','IntervalTaylor(c,[3-t,-1,0,0,0])')
    bind('collar_Gamma_C4','phi_jets','jets','(dist**(-2)*(-4)).exp()')
    t,T,W,lone,lp,lu,rp=s.symbols('t Ts wait logone logP logU logRp',real=True)
    physical=SimpleNamespace(ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,ln=s.log),
        logP=lp,logRp=rp,flatten=SimpleNamespace(logEv2_parts={'inlet_log':2*lu}),
        history=SimpleNamespace(heat=SimpleNamespace(mu=mu,a=a,bh=s.Rational(1,2)+a,outer=SimpleNamespace(Lrel=L))))
    logfn=asts.replay('current_angular_background_stress','angular_source_log_parts',{})
    logC=a-mu-r/2-k*T-k/2-lone
    vmin=dict(outer_angular=-4,steep_entry=0,steep_power=1,steep_exit=1+T,waiting=2+T,heat_collar=2+T+W)
    logs={}
    for region,v in vmin.items():
        parts=logfn(physical,v);B=sum(parts['B'].values())
        zero(region+'_Qz_Qtheta_B_units',sum(parts['Qz'].values())-sum(parts['Qtheta'].values())-B)
        zero(region+'_actual_B_decreases',s.diff(sum(logfn(physical,t)['B'].values()),t)+(s.Rational(1,2)+a))
        logs[region]=str(B+(logC if region in ('waiting','heat_collar') else 0))
    if not (heat.exterior.theorem['full_terminal_moment_stress_theorem_verified']
            and heat.exterior.theorem['heat_shear_coefficient_bounds_verified']
            and heat.exterior.proof['all_current_five_terminal_moment_functions_identified']
            and heat.Gamma_physical['all_three_physical_momentum_components_exactly_zero']):
        raise ValueError('Same full Gamma zero tensor/physical source theorem required')
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,live_current_graph=links,
        partial_A_source_lemma=dict(ell=str(ell),positive_density=str(density),FTC_derivative=str(derivative),
            conclusion='0<=same_partial_A(s)<=same_full_A; both are integrals of identical nonnegative density',
            same_normalization_and_five_weighted_integral_source_proof=internal.proof['checked_same_normalization_and_five_original_weight_function_sources'],
            same_support_FTC=internal.proof['checked_weighted_support_FTC_and_flat_beta_source_theorem'],
            enclosure_lower_zero_not_the_positivity_proof=True),
        current_Hf_uniform_lower=field.bounds['correlated_Hf_uniform_lower'],
        actual_current_waiting_heat_full_endpoint_proof=shared_future,
        current_four_internal_full_tensor_traces=internal.proof,current_seven_functional_tensor_joins=joins,
        current_full_pressure_function=ang.pressure_proof,current_heat_units_and_closed_source=heat.units,
        current_heat_full_normalization=heat.normalization,current_full_Gamma_zero_theorem=heat.exterior.theorem,
        current_Bmax_local_unit_recipes=logs,angular_G_min_is_one_not_old_KR=True,
        full_E_P_X_K_Z_retained=True,actual_axial_shear_zero_by_current_full_FTC=True,
        partial_weights_and_heat_caps_do_not_define_field_values=True,historical_cone_admissions_not_consumed=True,
        input_hashes=asts.hashes,passed=True)


def validate_whole_angular_tail_views(field,views):
    if set(views)!=set(DOMAINS):raise ValueError('All seven whole current angular/tail source views required')
    for region,domain in DOMAINS.items():
        view=views[region];raw=view['original_complete_view']
        if (view['actual_five_defect_family_sha256'],view['implicit_source_sha256'],view['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
            raise ValueError('Foreign angular/tail source')
        if endpoints(raw['Z'])!=(-1,1) or endpoints(raw[COORDINATES[region]])!=domain:
            raise ValueError('Whole original angular/tail domain required')
        if not raw['actual_full_stress_not_local_difference'] or set(raw['current_actual_normalized_full_moment_rows'])!=set(('A','E','P','K')):
            raise ValueError('Original full nonzero moment histories required')
        if not raw['stable_actual_absolute_pressure_mixed4_factored'] or not raw['original_current_energy_Taylor']:
            raise ValueError('Same absolute pressure and full energy required')
        for group in ('theta','axial','theta_inertial','theta_shear'):
            if not all(mp.isfinite(v) for v in endpoints(raw['current_actual_normalized_stress_mixed3'][group]['y0_Z0'])):
                raise ValueError('Finite original full stress source required')
        if region=='heat_exterior' and not raw['source_exact_Gamma_tensor_and_remainder_zero']:
            raise ValueError('Exact full exterior zero theorem required')
    return views


@source_precision
def whole_current_angular_tail_bounds(field,views):
    validate_whole_angular_tail_views(field,views);c=field.ctx;reg=field.registry
    owner=reg.owners['angular'];out=owner.outer;h=owner.heat;st=h.steep
    a=h.a;delta=h.delta;mu=h.mu;r=1-mu;k=1-a;b=(1-delta)/2;eps=h.eps
    def lo(x):return c.mpf(endpoints(c.mpf(x))[0])
    def hi(x):return c.mpf(endpoints(c.mpf(x))[1])
    def ab(x):return hi(abs(c.mpf(x)))
    raw={region:view['original_complete_view'] for region,view in views.items()}
    native=out.angular(c.mpf([-1,1]),c.mpf([-4,0]));data=out.data(c.mpf([-1,1]))
    repairs=[c.exp(r*(4+center))*out.weights['A']*(k*ab(d[0])+ab(b)*ab(d[1]))
        for d,center in zip(data['coeff'],(-3,-1))]
    hcap=ab(native['actual_angular_bump_y_derivatives'][0][0]);hscap=ab(native['actual_angular_bump_y_derivatives'][1][0])
    Fmin=lo(native['swirl_factor_one_plus_h_Taylor'][0]);Fmax=hi(native['swirl_factor_one_plus_h_Taylor'][0])
    g=lo((mu-a)/r);Wa=lo(g-sum(repairs,c.mpf(0))-hcap)
    mu_gap=lo(mu-hscap/Fmin);ma=hi(2*mu+2*hscap/Fmin)
    signgap=lo((1+mu)*Fmin-hscap)
    margins=dict(mu_minus_a=mu-a,rate=r,k=k,epsilon=eps,Lmin=1-delta,normalization=out.normalization,
        same_full_A_weight=out.weights['A'],current_Hf_lower=field.bounds['correlated_Hf_uniform_lower'],
        angular_W=Wa,angular_mu_gap=mu_gap,angular_shear_sign=signgap,angular_Fmin=Fmin,angular_m_below2=2-ma)
    def shear_cap(region,target,coefficient):
        source=c.ln(coefficient)-lo(raw[region]['exact_source_logR'])
        cap=c.ln(target)-1000;gap=cap-source
        margins[region+'_exact_shear_log_margin']=gap
        if endpoints(gap)[0]<=0:raise ArithmeticError('Same exact inverse radius shear cap failed')
        return hi(c.exp(cap)),dict(exact_shear_log_upper=source,relative_log_cap=cap,strict_log_margin=gap)
    acap,ash=shear_cap('outer_angular',Wa,2*((1+mu)*Fmax+hscap))
    thetaA=lo(Wa-acap)
    inlet_repairs=[c.exp(r*center)*out.weights['A']*(k*ab(d[0])+ab(b)*ab(d[1]))
        for d,center in zip(data['coeff'],(-3,-1))]
    M0=lo(g-sum(inlet_repairs,c.mpf(0)));Me=lo(M0*c.exp(-r))
    ecap,esh=shear_cap('steep_entry',Me,c.mpf(4))
    Ke=lo(raw['steep_entry']['current_actual_normalized_full_moment_rows']['K'][0][0])
    thetaE=lo(Ke*(Me-ecap));margins.update(entry_M0=M0,entry_M=Me,entry_K_lower=Ke)
    LT=a-mu-r/2-k*st.Ts-k/2;logC=LT-st.logone
    phys=owner.physical
    positions=dict(outer_angular=c.mpf(-4),steep_entry=c.mpf(0),steep_power=c.mpf(1),
        steep_exit=1+st.Ts,waiting=2+st.Ts,heat_collar=2+st.Ts+st.wait)
    def logB(region,extra=0,offset=0):
        parts=angular_source_log_parts(phys,positions[region]+offset)['B']
        return sum(parts.values(),c.mpf(0))+extra,parts
    regions={}
    def cone(region,theta,axial,m,mmin,logb,parts,unit='current'):
        theta=lo(theta);axial=ab(axial);m=hi(m);mmin=lo(mmin)
        if endpoints(theta)[0]<=0 or endpoints(mmin)[0]<=0:raise ArithmeticError('Positive actual theta/shear range failed: '+region)
        term=c.ln(m)+2*logb+2*c.ln(axial/theta);gap=c.ln(2)-term
        margins.update({region+'_theta_lower':theta,region+'_m_lower':mmin,region+'_directional_log_margin':gap})
        regions[region]=dict(theta_lower=theta,axial_absolute_upper=axial,kappa_minus2_lower=mmin,
            kappa_minus2_upper=m,current_log_Bmax=logb,current_Bmax_parts=parts,
            log_directional_term_upper=term,log_directional_threshold=c.ln(2),directional_log_margin=gap,
            bound_units=unit,original_full_raw_theta_enclosure=raw[region]['current_actual_normalized_stress_mixed3']['theta']['y0_Z0'],
            actual_source_shear_strictly_negative=True,actual_axial_shear_exact_zero=True,
            full_nonzero_pressure_energy_and_angular_histories_retained=True)
    for region,theta,m,mmin in (('outer_angular',thetaA,ma,2*mu_gap),('steep_entry',thetaE,2,2*mu),
        ('steep_power',lo(raw['steep_power']['current_actual_normalized_stress_mixed3']['theta']['y0_Z0']),2,2),
        ('steep_exit',lo(raw['steep_exit']['current_actual_normalized_stress_mixed3']['theta']['y0_Z0']),2,delta)):
        lb,parts=logB(region);cone(region,theta,raw[region]['current_actual_normalized_stress_mixed3']['axial']['y0_Z0'],m,mmin,lb,parts)
    regions['outer_angular'].update(continuous_W_lower=Wa,same_partial_A_absolute_bounds=repairs,
        actual_h_absolute_upper=hcap,actual_hs_absolute_upper=hscap,current_G_min=c.mpf(1),shear=ash)
    regions['steep_entry'].update(continuous_M0_lower=M0,continuous_M_lower=Me,shear=esh)
    # Same full current heat endpoint; no historical waiting receipt.
    base=collar_defect_rows(h,h.shape(c.mpf([-1,1]),0),h.collar_tails(c.mpf([-1,1]),0),c.mpf(0))
    modes=waiting_stress_rows(h,base,c.mpf([-1,1]),c.mpf(0))
    I0=modes['theta_inertial'][0][0]
    wcap=2*hi(h.Scap)*(1+a)*hi(1-eps)*c.exp(a*st.wait)
    wg=lo(I0-wcap);wz=ab(modes['axial_energy'][0][0])+ab(modes['axial_pressure'][0][0])
    lb,parts=logB('waiting',logC);cone('waiting',wg,wz,delta,delta,lb,parts,'local heat; theta*C and axial*C^2')
    regions['waiting'].update(actual_I0=I0,stable_full_heat_endpoint_defects=base,
        actual_energy_pressure_modes=modes,backwards_shear_relative_cap=wcap,exact_inverse_KR_log=logC,
        common_positive_factor='exp(-k*q)>=1, q in [-wait,0]',
        bracket_lower_formula='I0_lower-2*Scap*(1+a)*K0_upper*exp(a*wait)')
    shape=h.shape(c.mpf([-1,1]),c.mpf([0,3]));K0=shape['K_rows'][0][0];K1=shape['K_rows'][1][0]
    mlo=lo(delta-2*abs(K1)/K0);mhi=hi(delta+2*abs(K1)/K0)
    Dcap=2*a*(1+a)*hi(h.Scap);Hmin=lo(1-Dcap)
    sigcap=ab(sigma_jets(c,c.mpf([0,1]))[1])
    Gmin=lo(eps/c.exp(1)-Dcap)
    collar_margins=dict(sigmoid_G=Gmin,sigmoid_angular_future=k*eps/c.exp(1)-Dcap*(k+2*b),
        sigmoid_shear=eps*(1+a)/c.exp(1)-Dcap*(sigcap+eps+2+a),
        sigmoid_hot_waiting_gap=eps*(1-c.exp(-c.mpf(4)/9))-Dcap,
        heat_Hmin=Hmin,heat_angular_future=k*Hmin-4*b*a*(1+a)*hi(h.Scap),
        heat_shear=(1+a)*Hmin-2*a*(1+a)*hi(h.Scap),heat_m_lower=mlo,heat_m_below2=2-mhi)
    margins.update(collar_margins)
    Czsig=hi(3*(2*eps*(1+2*delta)+Dcap)/(1-delta))
    Czheat=hi((1+2*delta+4*a*(1+a)*hi(h.Scap))/(4*(1-delta)))
    lb,parts=logB('heat_collar',logC)
    cone('heat_collar',Gmin,Czsig,mhi,mlo,lb,parts,'local sigmoid part; theta*C and axial*C^2')
    for label in ('theta_lower','directional_log_margin'):
        margins['heat_collar_sigmoid_'+label]=margins.pop('heat_collar_'+label)
    sigmoid_bound=dict(regions['heat_collar'],domain='0<=t<=1')
    # The unscaled theta tends to zero. The full collar has no uniform
    # positive theta floor; keep its two domain proofs separately.
    regions['heat_collar']=dict(domain='0<=t<3',theta_uniform_lower_on_whole_collar=None,
        sigmoid_part=sigmoid_bound,actual_source_shear_strictly_negative=True,
        actual_axial_shear_exact_zero=True,full_nonzero_pressure_energy_and_angular_histories_retained=True)
    hblog,hparts=logB('heat_collar',logC,1)
    heatlog=c.ln(mhi)+2*hblog+2*c.ln(8*Czheat/Hmin)
    margins['heat_flat_part_directional_log_margin']=c.ln(2)-heatlog
    # Finite exact S gives a uniform y^6 direction bound at y=3-t ->0.
    # The prefactor may be enormous; its finiteness, not smallness, suffices.
    logRtail=sum(h.exact_logRtail_terms.values(),c.mpf(0))
    edgeprefactor=hblog+c.ln(Czheat)-c.ln(16)+hi(logRtail)+3-c.ln(Hmin)
    if not all(mp.isfinite(v) for v in endpoints(edgeprefactor)):raise ArithmeticError('Finite exact heat edge direction constant required')
    regions['heat_collar'].update(local_K=K0,local_Kt=K1,canonical_Gamma_deficit_cap=Dcap,
        original_sigma_t_cap=sigcap,sigmoid_positive_margins=collar_margins,
        heat_part=dict(domain='1<=t<3',theta_over_eps_phi_lower=Hmin,
            axial_over_eps_phi_y_cubed_upper=Czheat,current_log_Bmax=hblog,
            log_directional_term_upper=heatlog,flat_future_majorant='integral phi(t+u)/phi(t) du<=y^3/8'),
        terminal_direction=dict(domain='t->3-, uniformly Z[-1,1]',
            ratio_bound='abs(Tz/Ttheta)<=exp(finite_log_prefactor)*(3-t)^6',
            finite_log_prefactor=edgeprefactor,limit='e_theta',exact_S_positive_function_not_cap_endpoint=True),
        strict_only_before_t3=True,exact_zero_t3_excluded_from_strict_inequalities=True,exact_inverse_KR_log=logC)
    regions['heat_exterior']=dict(domain='every finite t>=3 and infinity limit',
        full_tensor_divergence_remainder_momentum_exactly_zero=True,
        strict_cone_for_zero_tensor=False,finite_view_anchor_not_domain_truncation=True)
    for name,value in margins.items():
        if endpoints(value)[0]<=0 or not all(mp.isfinite(v) for v in endpoints(value)):
            raise ArithmeticError('Current whole angular/tail margin failed: '+name)
    return dict(positive_margins=margins,regions=regions,current_inverse_KR_log=logC,
        whole_original_Z_and_continuous_coordinates=True,current_Hf_positive_term_retained=True,
        partial_A_bound_from_integrand_positivity_not_endpoint_overlap=True,
        original_full_nonzero_E_P_X_K_Z_and_absolute_pressure_retained=True,
        current_heat_zero_function_is_unbounded_not_finite_anchor=True,
        point_or_phase_sampling_used_as_proof=False,source_caps_used_as_defining_fields=False)
