"""Independent original steep kernels, functional O7 joins and mixed jets."""
import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_steep_waiting_C4 import transition_kernels
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_power_angular_C4_check import independent_nonconstant_rate_fixture
from lei_ren_part1_paper_compliant_pulse_physical_bounds import UZ,UT,UR,P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
NAME=PREFIX+'compliant_steep_waiting_C4.json'


def independent_original_kernel_fixture():
    """Point integration of the original sigmoid, separate from directed cells.

    An elementary certified fixture J is not substituted for the actual
    sigmoid. High-precision quadrature is a numerical diagnostic only;
    production enclosures rely on interval cells and original monotonicity.
    """
    with mp.workdps(60):
        c=MPIntervalContext(); c.dps=90; mu=mp.mpf('.04'); delta=mp.mpf('.003')
        def sig(t):
            if t<=0:return mp.mpf(0)
            if t>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/t**2-1/(1-t)**2))
        # Solve J'=sigma once with high-precision polynomial quadrature
        # evaluation. Source J(1)=1/2 is used only at the exact endpoint.
        cache={}
        def J(t):
            key=mp.nstr(t,65)
            if key not in cache:cache[key]=mp.quad(sig,[0,t],method='gauss-legendre') if t else mp.mpf(0)
            return cache[key]
        checks=0
        for kind in ('in','out'):
            rate=1-mu if kind=='in' else 1-delta/2
            angular=lambda v:mp.exp(rate*(v-J(v))) if kind=='in' else mp.exp(rate*J(v))
            energy=lambda v:mp.exp(-2*mu*v-2*rate*J(v)) if kind=='in' else mp.exp(-2*v+2*rate*J(v))
            pressure=lambda v:mp.exp(-(1+2*mu)*v-2*rate*J(v))/2 if kind=='in' else mp.exp(-3*v+2*rate*J(v))/2
            for t in (mp.mpf(0),mp.mpf('.37'),mp.mpf(1)):
                out=transition_kernels(c,c.mpf(t),c.mpf(mu),c.mpf(delta),kind,128)
                values=dict(J=J(t),angular=mp.quad(angular,[0,t],method='gauss-legendre') if t else 0,
                    pressure=mp.quad(pressure,[0,t],method='gauss-legendre') if t else 0,
                    remaining_energy=mp.quad(energy,[t,1],method='gauss-legendre') if t<1 else 0)
                for key,value in values.items():
                    lo,hi=endpoints(out[key])
                    # J(1) is an exact source identity, not quadrature's
                    # finite rounding. Zero/empty source integrals are exact.
                    if key=='J' and t==1:value=mp.mpf('.5')
                    if not lo<=value<=hi:raise ArithmeticError('Original steep kernel diagnostic failed: '+str((kind,t,key)))
                    checks+=1
        return dict(original_sigma_direct_quadrature_diagnostics=checks,
            numerical_diagnostic_not_interval_integral_proof=True,passed=True)


def functional_source_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Steep source identity failed: '+name)
        proofs[name]=True
    mu,d,t,T,W,J=s.symbols('mu delta t Ts W J',real=True)
    bp=s.Rational(1,2)+mu; r=1-mu; k=1-d/2; bh=s.Rational(1,2)+d/2
    XR,PR,thetaR,Ev2,H,Kentry,Kout,Ientry,Iout,Pentry,Pout=s.symbols('XR PR thetaR Ev2 H Kentry Kout Ientry Iout Pentry Pout',real=True)
    D=lambda a,l:(1-s.exp(-a*l))/a
    Qwait=D(d,W)+s.exp(-d*W)*H
    Qafter=Kout+s.exp(-1-d/2)*Qwait
    Qentry=s.exp(-1-mu)*(D(2,T)+s.exp(-2*T)*Qafter)
    XS=(XR+Ientry)*s.exp(-r/2); XQ=XS+T; XT=(XQ+Iout)*s.exp(-k/2)
    thetaS=thetaR*s.exp(-bp-r/2); thetaQ=thetaS*s.exp(-s.Rational(3,2)*T)
    thetaT=thetaQ*s.exp(-s.Rational(3,2)+k/2)
    PS=PR+Ev2*thetaR**2*Pentry; PQ=PS+Ev2*thetaS**2*D(3,T)/2
    PT=PQ+Ev2*thetaQ**2*Pout
    values={}
    # Bind actual NEW provider expressions to the original source formulas.
    # Treat angular/pressure initial data and the entire Gamma tail as
    # arbitrary functions, then differentiate interfaces in Z below.
    base={'c':s.Integer(0),'self.mu':mu,'self.delta':d,'self.rate':r,'self.k':k,'self.bp':bp,'self.bh':bh,
        'self.Ts':T,'self.wait':W,'t':t,'J':J,'one':s.Integer(1),'self.outer.flatten.Ev2':Ev2,
        'self.thetaR':thetaR,'self.thetaS':thetaS,'self.thetaQ':thetaQ,'self.thetaT':thetaT,
        "data['XR']":XR,"data['XS']":XS,"data['XQ']":XQ,"data['XT']":XT,
        "data['PR']":PR,"data['PS']":PS,"data['PQ']":PQ,"data['PT']":PT,
        "data['H']":H,"data['after_entry']":Qentry,"data['after_power']":Qafter,
        "data['waiting_future']":Qwait}
    L=s.symbols('Lrel',real=True)
    zero('new_exact_steep_power_length',assignment('compliant_steep_waiting_C4','__init__','self.Ts',base)-4*(s.log(2)-s.log(d)))
    zero('new_actual_epsilon_relation',assignment('compliant_steep_waiting_C4','__init__','self.epsilon',base)-d/1000)
    ibase=dict(base,**{'self.outer.Lrel':L})
    zero('new_angular_terminal_velocity_scale',assignment('compliant_steep_waiting_C4','__init__','self.thetaR',ibase)-s.exp(-bp*(100+L))/2)
    for target,expected in (('self.thetaS',thetaS),('self.thetaQ',thetaQ),('self.thetaT',thetaT)):
        zero('new_'+target[5:]+'_velocity_history',assignment('compliant_steep_waiting_C4','__init__',target,ibase)-expected)
    method='steep_in'; env=dict(base,**{"kernels['angular']":s.symbols('It'),"kernels['remaining_energy']":s.symbols('Kt'),"kernels['pressure']":s.symbols('Pt')})
    for target,expected in (
        ('theta',thetaR*s.exp(-bp*t-r*J)),('X',(XR+env["kernels['angular']"])*s.exp(-r*(t-J))),
        ('energy',(Qentry+env["kernels['remaining_energy']"])*s.exp(2*mu*t+2*r*J)/2),
        ('pressure',PR+env["kernels['pressure']"]*Ev2*thetaR**2)):
        actual=assignment('compliant_steep_waiting_C4',method,target,env); zero('new_entry_'+target,actual-expected); values['entry_'+target]=actual
    left=s.symbols('left',real=True); env=dict(base,left=left,**{'decay_integral(c, 3, t)':D(3,t)})
    for target,expected in (('theta',thetaS*s.exp(-s.Rational(3,2)*t)),('X',XS+t),
        ('energy',(D(2,left)+s.exp(-2*left)*Qafter)/2),('pressure',PS+Ev2*thetaS**2*D(3,t)/2)):
        actual=assignment('compliant_steep_waiting_C4','steep_power',target,env); zero('new_power_'+target,actual-expected); values['power_'+target]=actual
    env=dict(base,**{"kernels['angular']":s.symbols('Ot'),"kernels['remaining_energy']":s.symbols('Lt'),"kernels['pressure']":s.symbols('Qt')})
    for target,expected in (('theta',thetaQ*s.exp(-s.Rational(3,2)*t+k*J)),
        ('X',(XQ+env["kernels['angular']"])*s.exp(-k*J)),
        ('energy',(Qwait*s.exp(-1-d/2)+env["kernels['remaining_energy']"])*s.exp(2*t-2*k*J)/2),
        ('pressure',PQ+env["kernels['pressure']"]*Ev2*thetaQ**2)):
        actual=assignment('compliant_steep_waiting_C4','steep_out',target,env); zero('new_exit_'+target,actual-expected); values['exit_'+target]=actual
    env=dict(base,left=left,**{'decay_integral(c, self.delta, left)':D(d,left),'decay_integral(c, 1 + self.delta, t)':D(1+d,t)})
    for target,expected in (('theta',thetaT*s.exp(-bh*t)),('X',(XT-1/k)*s.exp(-k*t)+1/k),
        ('energy',(D(d,left)+s.exp(-d*left)*H)/2),('pressure',PT+Ev2*thetaT**2*D(1+d,t)/2)):
        actual=assignment('compliant_steep_waiting_C4','waiting',target,env); zero('new_waiting_'+target,actual-expected); values['waiting_'+target]=actual
    # Source normalized future decomposition retains every epsilon atom
    # and the full Gamma factor, with no finite-tail replacement.
    eps,A,B,S,hat,lone=s.symbols('eps A B S hat lone',real=True)
    heat=(1/d-2*eps*A+eps**2*B-d*S*hat/2)*s.exp(-2*lone)
    rawheat=heat/s.exp(-2*lone)
    heatcap_label='c.mpf([0, endpoints(self.delta * self.S / 2)[1]])'
    henv={'IntervalTaylor.constant(c, self.preheat_scalar, 5)':1/d-2*eps*A+eps**2*B,
        'heat':hat,heatcap_label:d*S/2}
    zero('new_exact_preheat_Gamma_source',assignment('compliant_steep_waiting_C4','data','preheat_tail',henv)-rawheat)
    zero('new_preheat_epsilon_normalization',assignment('compliant_steep_waiting_C4','data','H',
        {'preheat_tail':rawheat,'self.tail_normalization':s.exp(-2*lone)})-heat)
    zero('new_exact_waiting_epsilon_normalization',assignment('compliant_steep_waiting_C4','__init__','self.tail_normalization',
        {'self.logone':lone})-s.exp(-2*lone))
    Ns=s.exp(-1-mu); Nq=Ns*s.exp(-2*T); Nt=Nq*s.exp(-1-d/2)
    post_original=Kentry+Ns*D(2,T)+Nq*Kout+Nt*D(d,W)+Nt*s.exp(-d*W)*heat
    zero('same_original_postangular_future_decomposition',post_original-(Kentry+Qentry).subs(H,heat))
    denv=dict(base,**{'H':H,'decay_integral(c, self.delta, self.wait)':D(d,W),
        'waiting_future':Qwait,'self.outenergy':Kout,'after_power':Qafter,'decay_integral(c, 2, self.Ts)':D(2,T)})
    for target,expected in (('waiting_future',Qwait),('after_power',Qafter),('after_entry',Qentry)):
        zero('new_data_'+target,assignment('compliant_steep_waiting_C4','data',target,denv)-expected)
    historyenv=dict(base,XR=XR,PR=PR,XS=XS,XQ=XQ,XT=XT,PS=PS,PQ=PQ,
        **{"self.infull['angular']":Ientry,"self.outfull['angular']":Iout,
           "self.infull['pressure']":Pentry,"self.outfull['pressure']":Pout,
           'decay_integral(c, 3, self.Ts)':D(3,T)})
    for target,expected in (('XS',XS),('XQ',XQ),('XT',XT),('PS',PS),('PQ',PQ),('PT',PT)):
        zero('new_data_'+target+'_history',assignment('compliant_steep_waiting_C4','data',target,historyenv)-expected)
    # Independently normalize ORIGINAL corrected outer energy expressions.
    Nr=s.symbols('Nrel',positive=True); LF=s.symbols('LR',real=True)
    y0=s.symbols('y0',real=True); tr=s.exp(LF); relation=Nr-s.exp(y0)*tr**2
    zero('Nrel_is_Rrel_UthetaR_squared_in_Rv_Ev0_units',relation.subs({Nr:s.exp(y0+2*LF)}))
    oldenv={'self.Nrel':Nr,'self.Ns':Ns,'self.Nq':Nq,'self.Nt':Nt,
        "data['heat0']['energy']":H/s.exp(-2*lone),'self.Ntail':Nt*s.exp(-d*W),
        'self.tailmult':s.exp(-2*lone),"tails['energy']":s.symbols('Kt'),
        "self.fullout['energy']":Kout,'decay_integral(c, 2, self.params.Ts)':D(2,T),
        'decay_integral(c, self.delta, self.angular.waiting)':D(d,W),
        'self.mu':mu,'self.delta':d,'t':t,'left':left,
        'decay_integral(c, 2, left)':D(2,left),'decay_integral(c, self.delta, left)':D(d,left)}
    oldentry=assignment('compliant_corrected_outer_field','steep_in','E',oldenv)/Nr*s.exp(2*mu*t+2*r*J)/2
    zero('original_entry_backward_energy_units',oldentry-values['entry_energy'])
    oldpower=assignment('compliant_corrected_outer_field','steep_power','E',oldenv)/Nr*s.exp(1+mu+2*t)/2
    zero('original_power_backward_energy_units',oldpower.subs(left,T-t)-values['power_energy'].subs(left,T-t))
    oldenv["tails['energy']"]=s.symbols('Lt')
    oldexit=assignment('compliant_corrected_outer_field','steep_out','E',oldenv)/Nr*s.exp(1+mu+2*T+2*t-2*k*J)/2
    zero('original_exit_backward_energy_units',oldexit-values['exit_energy'])
    oldwait=assignment('compliant_corrected_outer_field','waiting','E',oldenv)/Nr*s.exp(2+mu+d/2+2*T+d*t)/2
    zero('original_waiting_backward_energy_units',oldwait.subs(left,W-t)-values['waiting_energy'].subs(left,W-t))
    oldbase=dict(base,**{"data['XR']":XR,"data['Xs']":XS,"data['Xq']":XQ,"data['Xt']":XT,
        'I':s.symbols('It'),'self.LR':s.log(thetaR),'self.Ls':s.log(thetaR)-bp-r/2,
        'self.Lq':s.log(thetaR)-bp-r/2-s.Rational(3,2)*T,
        'self.Lt':s.log(thetaR)-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)+k/2})
    for method,stage,log_expected in (
        ('steep_in','entry',s.log(thetaR)-bp*t-r*J),
        ('steep_power','power',s.log(thetaR)-bp-r/2-s.Rational(3,2)*t),
        ('steep_out','exit',s.log(thetaR)-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)*t+k*J),
        ('waiting','waiting',s.log(thetaR)-bp-r/2-s.Rational(3,2)*T-s.Rational(3,2)+k/2-bh*t)):
        oenv=dict(oldbase)
        if stage=='exit':oenv['I']=s.symbols('Ot')
        zero('original_'+stage+'_X_history',assignment('compliant_corrected_outer_field',method,'X',oenv)-values[stage+'_X'])
        zero('original_'+stage+'_theta_history',assignment('compliant_corrected_outer_field',method,'L',oenv)-log_expected)
    # Bind rate lists as production syntax, then prove their coefficients
    # against the original slope equations. The factorials convert sigma
    # Taylor coefficients to ordinary y derivatives, with no missing chain.
    tree=ast.parse((HERE/(PREFIX+'compliant_steep_waiting_C4.py')).read_text(encoding='utf8'))
    expected_rates={
        'steep_in':'[-self.rate*sig[j]*math.factorial(j) for j in range(4)]',
        'steep_power':'[-self.rate,c.mpf(0),c.mpf(0),c.mpf(0)]',
        'steep_out':'[-self.rate+self.k*sig[0]]+[self.k*sig[j]*math.factorial(j) for j in range(1,4)]',
        'waiting':'[self.mu-self.delta/2,c.mpf(0),c.mpf(0),c.mpf(0)]'}
    for method,expected in expected_rates.items():
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        found=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)=='rates' for v in n.targets)]
        if len(found)!=1 or ast.dump(found[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Actual O7 ordinary-rate source changed: '+method)
        proofs['production_'+method+'_ordinary_rate_list_bound']=True
    sig0=s.symbols('sigma0',real=True)
    for name,a0,theta_rate,angular_rate,energy_rate in (
        ('entry',-r*sig0,-bp-r*sig0,r*(1-sig0),2*mu+2*r*sig0),
        ('power',-r,-s.Rational(3,2),0,2),
        ('exit',-r+k*sig0,-s.Rational(3,2)+k*sig0,k*sig0,2-2*k*sig0),
        ('waiting',mu-d/2,-bh,k,d)):
        zero(name+'_production_theta_rate',a0-bp-theta_rate)
        zero(name+'_production_angular_rate',r+a0-angular_rate)
        zero(name+'_production_energy_rate',2*mu-2*a0-energy_rate)
    # Each original unit kernel is defined by its exact integral. The
    # following substitutions are the fundamental theorem of calculus,
    # not finite differences or an interval overlap criterion.
    Jt=s.Function('J')(t); sig=s.Function('sigma')(t)
    ft=lambda expr:s.diff(expr,t).subs(s.diff(Jt,t),sig)
    for stage,theta,X,e,pressure,a0,subs in (
        ('entry',values['entry_theta'],values['entry_X'],values['entry_energy'],values['entry_pressure'],-r*sig,
            {J:Jt,s.symbols('It'):s.Function('I')(t),s.symbols('Kt'):s.Function('K')(t),s.symbols('Pt'):s.Function('Pint')(t)}),
        ('exit',values['exit_theta'],values['exit_X'],values['exit_energy'],values['exit_pressure'],-r+k*sig,
            {J:Jt,s.symbols('Ot'):s.Function('I')(t),s.symbols('Lt'):s.Function('K')(t),s.symbols('Qt'):s.Function('Pint')(t)})):
        th,xx,ee,pp=[v.subs(subs) for v in (theta,X,e,pressure)]
        if stage=='entry':
            deriv={s.diff(s.Function('I')(t),t):s.exp(r*(t-Jt)),s.diff(s.Function('K')(t),t):-s.exp(-2*mu*t-2*r*Jt),
                s.diff(s.Function('Pint')(t),t):s.exp(-(1+2*mu)*t-2*r*Jt)/2}
        else:
            deriv={s.diff(s.Function('I')(t),t):s.exp(k*Jt),s.diff(s.Function('K')(t),t):-s.exp(-2*t+2*k*Jt),
                s.diff(s.Function('Pint')(t),t):s.exp(-3*t+2*k*Jt)/2}
        zero(stage+'_theta_source_ODE',ft(th)-(a0-bp)*th)
        zero(stage+'_angular_source_ODE',ft(xx).subs(deriv)-1+(r+a0)*xx)
        zero(stage+'_energy_source_ODE',ft(ee).subs(deriv)+s.Rational(1,2)-(2*mu-2*a0)*ee)
        zero(stage+'_pressure_source_ODE',ft(pp).subs(deriv)-Ev2*th**2/2)
    for stage,length,a0 in (('power',T,-r),('waiting',W,mu-d/2)):
        th,xx,ee,pp=[values[stage+'_'+label].subs(left,length-t) for label in ('theta','X','energy','pressure')]
        zero(stage+'_theta_source_ODE',s.diff(th,t)-(a0-bp)*th)
        zero(stage+'_angular_source_ODE',s.diff(xx,t)-1+(r+a0)*xx)
        zero(stage+'_energy_source_ODE',s.diff(ee,t)+s.Rational(1,2)-(2*mu-2*a0)*ee)
        zero(stage+'_pressure_source_ODE',s.diff(pp,t)-Ev2*th**2/2)
    # Exact endpoint function equalities and identical ordinary-y rates
    # determine all mixed spatial jets at the four local interfaces.
    ep_in0={t:0,J:0,s.symbols('It'):0,s.symbols('Kt'):Kentry,s.symbols('Pt'):0}
    ep_in1={t:1,J:s.Rational(1,2),s.symbols('It'):Ientry,s.symbols('Kt'):0,s.symbols('Pt'):Pentry}
    ep_out0={t:0,J:0,s.symbols('Ot'):0,s.symbols('Lt'):Kout,s.symbols('Qt'):0}
    ep_out1={t:1,J:s.Rational(1,2),s.symbols('Ot'):Iout,s.symbols('Lt'):0,s.symbols('Qt'):Pout}
    labels=('theta','X','energy','pressure')
    Z=s.symbols('Z',real=True)
    history={v:s.Function(str(v))(Z) for v in (XR,PR,H)}
    angular=[thetaR,XR,post_original.subs(heat,H)/2,PR]
    # post_original contains H only after replacing the complete source
    # heat expression. Use its already-proven identical source factor.
    angular[2]=(Kentry+Qentry)/2
    interfaces=(('angular_entry',angular,[values['entry_'+a].subs(ep_in0) for a in labels],0),
        ('entry_power',[values['entry_'+a].subs(ep_in1) for a in labels],
            [values['power_'+a].subs({t:0,left:T}) for a in labels],-r),
        ('power_exit',[values['power_'+a].subs({t:T,left:0}) for a in labels],
            [values['exit_'+a].subs(ep_out0) for a in labels],-r),
        ('exit_waiting',[values['exit_'+a].subs(ep_out1) for a in labels],
            [values['waiting_'+a].subs({t:0,left:W}) for a in labels],mu-d/2))
    for name,a,b,a0 in interfaces:
        for label,lhs,rhs in zip(labels,a,b):
            zero(name+'_exact_'+label,lhs-rhs)
            for n in range(6):zero(name+'_'+label+'_axial'+str(n),s.diff((lhs-rhs).subs(history),Z,n))
        leftrows=list(a); rightrows=list(b)
        for order in range(5):
            for n in range(5-order):
                for label,lhs,rhs in zip(labels,leftrows,rightrows):
                    zero(name+'_'+label+'_y'+str(order)+'_Z'+str(n),s.diff((lhs-rhs).subs(history),Z,n))
            # sigma is flat at both endpoints, so the boundary logarithmic
            # rates are constant through the derivative orders required.
            leftrows=[(a0-bp)*leftrows[0],(1 if order==0 else 0)-(r+a0)*leftrows[1],
                (-s.Rational(1,2) if order==0 else 0)+(2*mu-2*a0)*leftrows[2],
                Ev2*a[0]**2/2*(2*(a0-bp))**order]
            rightrows=[(a0-bp)*rightrows[0],(1 if order==0 else 0)-(r+a0)*rightrows[1],
                (-s.Rational(1,2) if order==0 else 0)+(2*mu-2*a0)*rightrows[2],
                Ev2*b[0]**2/2*(2*(a0-bp))**order]
    zero('waiting_exit_full_infinite_Gamma_future',values['waiting_energy'].subs({t:W,left:0})-H/2)
    return proofs


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Steep source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    checks=zeros=0
    points=[p for group in r['samples'].values() for p in group]+[p for group in r['whole_Z'].values() for p in group.values()]
    for p in points:
        if endpoints(read(p['energy_Taylor']['coefficients'][0]))[0]<=0:raise ArithmeticError('Whole-domain future energy positivity lost')
        for flag in ('same_actual_angular_terminal_histories_retained','exact_Gamma_and_both_epsilon_atoms_retained',
            'no_forward_subtraction_of_unrelated_long_future_energy'):
            if not p[flag]:raise ValueError('Original histories/positive correlation missing')
        for label,grid in p['physical_mixed_derivatives_total_order_le4'].items():
            if len(grid)!=15:raise ValueError('Mixed-order coverage incomplete')
            for value in grid.values():
                ends=endpoints(read(value))
                if any(not mp.isfinite(v) for v in ends):raise ArithmeticError('Nonfinite O7 derivative')
                checks+=1
                if label in (UZ,UR):
                    if ends!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Actual inherited axial/radial zero lost')
                    zeros+=1
        for flag in ('full_outer_C4_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
            if p[flag]:raise ValueError('O7 scope promoted')
    for kind in ('entry','exit'):
        for end in ('inlet','terminal'):
            if endpoints(read(r['whole_Z'][kind][end]['original_transition_kernels']['J']))!=(mp.mpf(0 if end=='inlet' else '.5'),)*2:
                raise ArithmeticError('Original exact transition J endpoint lost')
    if endpoints(read(r['whole_Z']['power']['entire']['energy_Taylor']['coefficients'][0]))[0]<mp.mpf('.25'):
        raise ArithmeticError('Uniform steep-power energy floor lost')
    overlaps=prefix_checks=0
    def overlap(a,b,label):
        nonlocal overlaps
        a,b=endpoints(read(a)),endpoints(read(b))
        if max(a[0],b[0])>min(a[1],b[1]):raise ArithmeticError('Independent same-source diagnostic disjoint: '+label)
        overlaps+=1
    prior=json.loads((HERE/(PREFIX+'compliant_power_angular_C4.json')).read_bytes())
    interfaces=[('angular_entry',prior['whole_Z_angular_terminal'],r['whole_Z']['entry']['inlet'])]
    interfaces +=[(a+'_'+b,r['whole_Z'][a]['terminal'],r['whole_Z'][b]['inlet']) for a,b in (('entry','power'),('power','exit'),('exit','waiting'))]
    for name,left,right in interfaces:
        for label,grid in left['physical_mixed_derivatives_total_order_le4'].items():
            for index,value in grid.items():overlap(value,right['physical_mixed_derivatives_total_order_le4'][label][index],name+' '+label+' '+index)
        for key in ('angular_y_derivative_Taylor','energy_y_derivative_Taylor'):
            for a,b in zip(left[key],right[key]):
                for av,bv in zip(a['coefficients'],b['coefficients']):overlap(av,bv,name+' '+key)
    old=json.loads((HERE/(PREFIX+'compliant_corrected_outer_field.json')).read_bytes())
    stages={'corrected O.7 steep entry':'entry','corrected O.7 steep power':'power','corrected O.7 steep exit':'exit','corrected O.7 waiting':'waiting'}
    for p in old['samples']:
        if p['stage'] not in stages:continue
        name=stages[p['stage']]; coord=p['coordinate']; key='phase' if name in ('power','waiting') else 'offset'
        current=next(q for q in r['samples'][name] if endpoints(read(q['Z']))==endpoints(read(p['Z']))
            and endpoints(read(q['coordinate'][key]))==endpoints(read(coord[key])))
        for key,prior_key in (('theta_over_Ev0_Taylor','Utheta_over_Ev0_enclosure'),
            ('angular_Taylor','Mtheta_over_sqrt2_R_3half_Utheta'),('energy_Taylor','Mztheta_over_R_Utheta_squared'),
            ('pressure_over_Pstar_squared_Taylor','P_over_Pstar_squared')):
            for a,b in zip(current[key]['coefficients'],p[prior_key]['coefficients']):overlap(a,b,name+' earlier C1 '+key); prefix_checks+=1
    source_proofs=functional_source_identities()
    fixture=independent_nonconstant_rate_fixture()
    kernel_fixture=independent_original_kernel_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        all_passed=True,functional_production_source_identities=source_proofs,
        independent_nonconstant_rate_fixture=fixture,independent_original_kernel_fixture=kernel_fixture,
        finite_actual_mixed_bounds_checked=checks,inherited_exact_zero_bounds_checked=zeros,
        independent_interface_and_prefix_overlap_diagnostics=overlaps,earlier_actual_C1_prefix_diagnostics=prefix_checks,
        entire_steep_waiting_high_mixed_derivatives_available=True,uniform_complete_future_energy_strictly_positive=True,
        angular_steep_and_internal_joins_certified=True,waiting_exit_Gamma_future_source_retained=True,
        join_scope='Leading spatial/profile mixed<=4 and primitive axial5, angular/O7 and three internal O7 joins; collar external join excluded',
        full_pulse_C4_installed=True,full_outer_C4_certified=False,
        collar_and_Gamma_high_mixed_derivatives_available=False,physical_energy_integral_certified=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Actual O7 steep/waiting: original kernels, exact source ODEs, mixed interface identities and positive future energy PASS',flush=True)
    return out


if __name__=='__main__':run()
