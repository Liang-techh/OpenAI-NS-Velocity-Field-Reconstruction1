"""Exact source joins, angular equations and direct radial-energy checks."""
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, roots_genlaguerre
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions_check as Gamma_check
from lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions_check import evaluate as bump_integral_value


def angular_fixture():
    g=current.source.FunctionTransportGraph();z=g.symbol('Z')
    c=g.constant;Z=sy.Symbol('Z')
    expression=lambda f:current.source.expression(g,sy.sympify(f),{'Z':z})
    pair=lambda f:current.source.C1Function(expression(f),expression(sy.diff(f,Z)))
    # Known nonconstant solution of the ORIGINAL two equations. Build
    # its RHSs independently and test the small-root and implicit derivative.
    p,q,ks=sy.Rational(1,7),sy.Rational(1,9),sy.Rational(1,10)
    x=Z/100+sy.Rational(1,200);y=Z**2/150-sy.Rational(1,250)
    b1=p*x+y;b2=x+q*y+ks*(x*x+q*y*y)
    selected,eq=current.small_quadratic_functions(g,pair(b1),pair(b2),
        c('1/7'),c('1/9'),c(1),c('1/10'))
    evaluate=Gamma_check.FixtureEvaluator(g,.3)
    for point in (-.8,0.,.7):
        env={z.node:point}
        xx,yy=evaluate(selected['x'].value,env),evaluate(selected['y'].value,env)
        assert abs(xx-float(x.subs(Z,point)))<1e-14
        assert abs(yy-float(y.subs(Z,point)))<1e-14
        matrix=np.array([[float(p),1.],[1+2*float(ks)*xx,float(q)*(1+2*float(ks)*yy)]])
        rhs=np.array([float(sy.diff(b1,Z).subs(Z,point)),float(sy.diff(b2,Z).subs(Z,point))])
        want=np.linalg.solve(matrix,rhs)
        assert abs(evaluate(selected['x'].Z,env)-want[0])<1e-14
        assert abs(evaluate(selected['y'].Z,env)-want[1])<1e-14
    P,Q,K,S,X,B1,B2=sy.symbols('p q k scale x b1 b2')
    Y=B1-P*X
    expanded=X+Q*Y+K*S*(X*X+Q*Y*Y)-B2
    expected=K*S*(1+Q*P*P)*X*X+(1-P*Q-2*K*S*Q*P*B1)*X+Q*B1+K*S*Q*B1*B1-B2
    assert sy.expand(expanded-expected)==0
    # Independent actual normalized beta integrals, with finite mu.
    norm=quad(lambda t:np.exp(-1/(1-t*t)) if abs(t)<1 else 0.,-1,1,epsabs=1e-13)[0]
    W=current.angular_weights(g,c('7/100'),c(str(norm)))
    beta=lambda s:np.exp(-1/(1-(s/.15)**2))/(.15*norm) if abs(s)<.15 else 0.
    for name,rate,power in (('A',.93,1),('B',-1.14,1),('D',-1.14,2),('E',-.14,1),('F',-.14,2)):
        want=quad(lambda s:np.exp(rate*s)*beta(s)**power,-.15,.15,epsabs=1e-12)[0]
        assert abs(float(bump_integral_value(g,W[name]))-want)<2e-11,name
    return dict(passed=True,mathematical_fixture_only=True,
        original_two_angular_equations_solved=True,nonconstant_C1_matches_independent_implicit_inverse=True,
        exact_eliminated_polynomial_identity=True,all_five_finite_mu_bump_integrals_checked=True)


def energy_fixture():
    g=current.source.FunctionTransportGraph();c=g.constant;z=g.symbol('Z')
    mu,delta,eps,wait,S=.02,.3,.02,1.3,.12
    r,k,bp,bh=1-mu,1-delta/2,.5+mu,.5+delta/2
    Lrel,Ts=-30*np.log(mu),4*np.log(2/delta)
    schedule=dict(mu=c('1/50'),delta=c('3/10'),epsilon=c('1/50'),r=c('49/50'),k=c('17/20'),
        waiting=c('13/10'),S=c('3/25'),Lrel=g.mul(c(-30),g.unary('log',c('1/50'))),
        Ts=g.mul(c(4),g.unary('log',c('20/3'))),log_one_minus_epsilon=g.unary('log1p',c('-1/50')))
    t=g.symbol('fixture_preheat_t');sig=current.Gamma.sigma(g,t);ph=current.Gamma.phi(g,t)
    W=g.add(g.sub(g.one,sig),g.mul(sig,ph));weight=g.unary('exp',g.mul(c('-3/10'),t))
    schedule['preheat_atoms']={'energy':{key:current.controls.integral(g,g.mul(weight,body),
        g.nodes[t.node]['name'],g.zero,c(3)) for key,body in (('W',W),('W_squared',g.mul(W,W)))}}
    heat,_=current.Gamma.Gamma_future_functions(g,schedule['delta'],schedule['epsilon'],schedule['S'])
    Z=sy.Symbol('Z');expr=lambda f:current.source.expression(g,sy.sympify(f),{'Z':z})
    pair=lambda f:current.source.C1Function(expr(f),expr(sy.diff(f,Z)))
    ds=[-sy.Rational(1,1000)+Z/5000,sy.Rational(3,2000)+Z**2/5000]
    angular=dict(d1=pair(ds[0]),d2=pair(ds[1]))
    norm=quad(lambda t:np.exp(-1/(1-t*t)) if abs(t)<1 else 0.,-1,1,epsabs=1e-13)[0]
    weights=current.angular_weights(g,schedule['mu'],c(str(norm)))
    functions,pieces=current.complete_energy_functions(g,schedule,heat,angular,weights)
    evaluate=Gamma_check.FixtureEvaluator(g,delta)
    # These fixture bindings are independently integrated fixed weights;
    # they are not production substitutions for the actual integral nodes.
    weight_values={ref.node:float(bump_integral_value(g,ref)) for ref in weights.values()}
    vv,ww=roots_genlaguerre(128,delta/2)
    H=lambda Z,tt:np.dot(ww,np.exp(-delta/2*np.log1p(2*(1-Z*Z)*S*np.exp(-tt)*vv)))/gamma(1+delta/2)
    sigma=lambda tt:float(Gamma_check.flat_sigma(tt))
    phi=lambda tt:float(Gamma_check.flat_phi(tt))
    J=lambda tt:quad(sigma,0,tt,epsabs=2e-12)[0]
    beta=lambda s:np.exp(-1/(1-(s/.15)**2))/(.15*norm) if abs(s)<.15 else 0.
    qint=lambda fn,lo,hi:quad(fn,lo,hi,epsabs=2e-12,epsrel=2e-12,limit=160)[0]
    def direct(Z):
        Q=1+Z*Z;yrel=100+Lrel
        d1,d2=[float(d.subs(sy.Symbol('Z'),Z)) for d in ds]
        def U_power(y):
            bump=d1*beta(y-yrel+3)+d2*beta(y-yrel+1)
            return Q/2*np.exp(-bp*y)*(1+bump)
        flatten=qint(lambda y:np.exp(y)*(np.exp(-bp*y)*(Q/2)**sigma(y/100))**2,0,100)
        breaks=[100,yrel-3.15,yrel-2.85,yrel-1.15,yrel-.85,yrel]
        power=sum(qint(lambda y:np.exp(y)*U_power(y)**2,lo,hi) for lo,hi in zip(breaks,breaks[1:]))
        base=qint(lambda y:np.exp(y)*(Q/2*np.exp(-bp*y))**2,100,yrel)
        Erel=Q/2*np.exp(-bp*yrel)
        inside=qint(lambda tt:np.exp(yrel+tt)*(Erel*np.exp(-bp*tt-r*J(tt)))**2,0,1)
        ys=yrel+1;Es=Erel*np.exp(-bp-r/2)
        steep=qint(lambda tt:np.exp(ys+tt)*(Es*np.exp(-1.5*tt))**2,0,Ts)
        yq=ys+Ts;Eq=Es*np.exp(-1.5*Ts)
        outside=qint(lambda tt:np.exp(yq+tt)*(Eq*np.exp(-1.5*tt+k*J(tt)))**2,0,1)
        yt=yq+1;Et=Eq*np.exp(-1.5+k/2)
        waiting=qint(lambda tt:np.exp(yt+tt)*(Et*np.exp(-bh*tt))**2,0,wait)
        ytail=yt+wait;Etail=Et*np.exp(-bh*wait)
        factor=np.exp(ytail)*(Etail/(1-eps))**2
        def brackets(tt):
            sig,ph=sigma(tt),phi(tt)
            pre=1-eps*(1-sig+sig*ph)
            actual=(1-sig)*(1-eps)+sig*(1-eps*ph)*H(Z,tt)
            return pre,actual
        tail=factor*(qint(lambda tt:np.exp(-delta*tt)*brackets(tt)[1]**2,0,3)
            +qint(lambda tt:np.exp(-delta*tt)*H(Z,tt)**2,3,np.inf))
        pre_tail=factor*(qint(lambda tt:np.exp(-delta*tt)*brackets(tt)[0]**2,0,3)+np.exp(-3*delta)/delta)
        heat=factor*(qint(lambda tt:np.exp(-delta*tt)*(brackets(tt)[0]**2-brackets(tt)[1]**2),0,3)
            +qint(lambda tt:np.exp(-delta*tt)*(1-H(Z,tt)**2),3,np.inf))
        total=flatten+power+inside+steep+outside+waiting+tail
        return dict(complete=total,flatten=flatten,power=base,angular_change=power-base,
            steep_waiting_preheat_tail=inside+steep+outside+waiting+pre_tail,Gamma_heat_difference=heat)
    errors={}
    for point in (-.4,.65):
        env={**weight_values,z.node:point};want=direct(point)
        for name,pair_ref in pieces.items():
            if name not in want:continue
            actual=evaluate(pair_ref.value,env)
            error=abs(actual-want[name])/max(abs(want[name]),1e-40)
            errors[name]=max(errors.get(name,0.),error)
            assert error<2e-7,(name,error)
        value=evaluate(functions['complete'].value,env)
        assert abs(value-want['complete'])<2e-9
        h=1e-4
        derivative=(direct(point+h)['complete']-direct(point-h)['complete'])/(2*h)
        assert abs(evaluate(functions['complete'].Z,env)-derivative)<2e-7
        assert abs(evaluate(functions['weighted'].value,env)-mu*np.exp(-26)/2*want['complete'])<2e-22
    return dict(passed=True,mathematical_fixture_only=True,original_parameters_not_replaced=True,
        actual_piecewise_swirl_integrated_with_physical_radial_measure=True,
        signed_two_bump_cross_and_square_terms_retained=True,
        independent_whole_future_value_and_first_Z_checked=True,
        exact_nonzero_heat_prefactor_checked_separately=True,
        entire_positive_Gamma_exterior_integrated=True,maximum_relative_piece_errors=errors)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1CompleteFutureFunctions(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert owner.Gamma.acceptance_loaded and owner.Gamma.pulse.acceptance_loaded
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['function_bindings']==owner.snapshot
    assert not any(raw[key] for key in current.OPEN)
    assert owner.source_identity['passed'] and owner.source_identity['actual_native_Xv_inside_archived_candidate_interval']
    assert owner.branch_certificate['passed'] and owner.branch_certificate['exact_root_in_original_unique_small_branch']
    assert len(owner.angular)==4 and set(owner.functions)=={'complete','weighted'}
    for pair in [*owner.angular.values(),*owner.functions.values(),*owner.pulse_functions.values()]:
        assert type(pair) is current.source.C1Function
        owner.g.ids([pair.value,pair.Z])
    fixtures=dict(angular_equations_and_bump_integrals=angular_fixture(),
        independent_full_physical_energy=energy_fixture(),
        actual_source_identity=owner.source_identity,uniform_original_small_branch=owner.branch_certificate,
        fresh_actual_pulse_branch=owner.pulse_branch_certificate)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=fixtures,
        **{current.GATE:True,current.ANGULAR:True,current.PULSE:True},**dict.fromkeys(current.OPEN,False),
        original_pulse_positive_branch_proved=True,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_COMPLETE_FUTURE_AND_SMALL_ANGULAR_C1_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
