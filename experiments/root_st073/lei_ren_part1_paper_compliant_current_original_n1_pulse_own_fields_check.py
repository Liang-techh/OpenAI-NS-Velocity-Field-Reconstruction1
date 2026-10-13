"""Independent complete radial histories, physical units and pulse joins.

Finite fixtures test the actual defining formulas with nonunit Pstar and
nonzero finite mu. They do not replace the original microscopic parameters
or certify uniform original numerical accuracy or full five-moment repair.
"""
import json
import math
from pathlib import Path
import time
import numpy as np
import sympy as sy
from scipy.integrate import quad, solve_ivp
from scipy.special import expit, exprel

import lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields as current


def sigma(x):
    if x <= 0:return 0.
    if x >= 1:return 1.
    return float(expit(-1/x**2+1/(1-x)**2))


def sigma_prime(x):
    if not 0 < x < 1:return 0.
    s=sigma(x)
    return s*(1-s)*(2/x**3+2/(1-x)**3)


def beta(x):
    return math.exp(-1/(1-x*x)) if abs(x)<1 else 0.


def beta_prime(x):
    return -2*x*math.exp(-1/(1-x*x))/(1-x*x)**2 if abs(x)<1 else 0.


def primitive(x):
    if x <= 0:return 0.
    # Exact symmetry sigma(u)+sigma(1-u)=1 on the whole transition.
    if x >= .02:return x-.01
    return quad(lambda v:sigma(50*v),0,x,epsabs=2e-13,epsrel=2e-13)[0]


def gp(x):return primitive(x)*sigma(11-x)
def gp_prime(x):return sigma(50*x)*sigma(11-x)-primitive(x)*sigma_prime(11-x)


class FixtureEvaluator:
    """Evaluate the defining finite DAG, including complete own histories."""
    def __init__(self,g,mu):
        self.g,self.mu=g,mu
        self.variables={r['name']:i for i,r in enumerate(g.nodes) if r['operation']=='bound_variable'}

    def __call__(self,ref,env=None):return float(self.value(ref.node,env or {}))

    def value(self,index,env):
        cache={}
        def at(i):
            if i in env:return env[i]
            if i in cache:return cache[i]
            row=self.g.nodes[i];op=row['operation']
            if op=='exact_rational':out=row['numerator']/row['denominator']
            elif op=='sum':out=sum(at(j) for j in row['arguments'])
            elif op=='product':out=math.prod(at(j) for j in row['arguments'])
            elif op=='negative':out=-at(row['argument'])
            elif op=='positive_quotient':out=at(row['numerator'])/at(row['denominator'])
            elif op=='analytic_unary':
                x=at(row['argument']);name=row['name']
                if name=='original_flat_sigma_prime':out=sigma_prime(x)
                elif name=='compact_raw_beta_derivative':out=beta_prime(x)
                else:out=float(exprel(x)) if name=='exprel' else getattr(math,name)(x)
            elif op=='exact_original_source_flat_sigma':out=sigma(at(row['argument']))
            elif op=='exact_original_source_flat_sigma_derivative':out=sigma_prime(at(row['argument']))
            elif op=='compact_raw_beta':out=beta(at(row['argument']))
            elif op=='compact_raw_beta_derivative':out=beta_prime(at(row['argument']))
            elif op=='function_substitution':
                out=self.value(row['expression'],{**env,row['variable']:at(row['value'])})
            elif op=='definite_integral':
                name=row['variable'];lo,hi=at(row['lower']),at(row['upper'])
                variable=self.variables[name]
                if name=='original_pulse_gp_primitive_v':out=primitive(hi)-primitive(lo)
                else:
                    if name=='original_own_main_history_v':breaks=[.02/self.mu,10/self.mu,11/self.mu]
                    elif name=='original_own_end_history_v':breaks=[-3.15,-3.,-2.85,-1.15,-1.,-.85]
                    elif name=='original_pulse_xi':breaks=[.02,10.]
                    else:breaks=[0.]
                    points=[lo]+[v for v in breaks if lo<v<hi]+[hi]
                    body=lambda v:self.value(row['integrand'],{**env,variable:v})
                    out=sum(quad(body,l,r,epsabs=2e-12,epsrel=2e-12,limit=150)[0]
                        for l,r in zip(points,points[1:]))
            else:raise ValueError('No independent own-field fixture evaluator for '+op)
            cache[i]=out
            return out
        return at(index)


def rows(pair):return [pair.value,pair.Z,pair.ZZ]+([pair.ZZZ] if hasattr(pair,'ZZZ') else [])


def convolve(a,b):
    return np.array([sum(math.comb(n,k)*a[k]*b[n-k] for k in range(n+1)) for n in range(len(b))])


def Q_reference(Z,delta,U,m):
    """Independent Taylor-series quotient; inputs are ordinary ODE rows."""
    mul=lambda a,b:np.convolve(a,b)[:3]
    z=np.array([Z,1.,0.]);square=mul(z,z);one=np.array([1.,0.,0.])
    us=np.array(U[:3])/np.array([1.,1.,2.]);ms=np.array(m[:3])/np.array([1.,1.,2.])
    mz=np.array(m[1:4])/np.array([1.,1.,2.])
    N=2*mul(z,us)-(1-delta)*mul(z,ms)-mul(one-square,mz)
    L=one-delta*square;out=[]
    for n in range(3):out.append((N[n]-sum(L[k]*out[n-k] for k in range(1,n+1)))/L[0])
    return np.array(out)*np.array([1.,1.,2.])


def own_fields_fixture():
    mu,delta,Pstar,logRp=.07,.3,2.3,1.2
    g=current.source.FunctionTransportGraph();c=g.constant;z=g.symbol('Z');Z=sy.Symbol('Z')
    a=current.c3.C3Algebra(g)
    expressions=dict(u=sy.Rational(2,5)/(1+Z**2),m1=sy.Rational(3,100)*Z*(1+Z**2),
        m2=sy.Rational(1,20)*Z*(1+Z**2),P0=sy.Rational(3,10)*(1+Z**2),
        Mp=sy.Rational(1,5)/(1+Z**2)**2)
    ap_expr=sy.Rational(51,50)+sy.Rational(3,200)*Z+Z**2/1000
    make=lambda expr:current.c3.C3Function(*(current.source.expression(g,sy.diff(expr,Z,n),{'Z':z}) for n in range(4)))
    native={key:make(expr) for key,expr in expressions.items()}
    ap=make(ap_expr);W=current.jets.complete.Gamma.pulse.kernel_functions(g,c('7/100'))
    A,D,det=W['first_row'],W['divided_row'],W['positive_determinant_absolute']
    inv=lambda r1,r2:[a.div(a.add(a.scale(r2,A[1]),a.neg(a.scale(r1,D[1]))),a.fixed(det)),
        a.div(a.add(a.scale(r1,D[0]),a.neg(a.scale(r2,A[0]))),a.fixed(det))]
    incoming=[a.scale(native[key],factor) for key,factor in zip(('m1','m2'),W['incoming_factors'])]
    u=inv(a.neg(a.scale(incoming[0],W['mu'])),a.neg(a.add(incoming[1],a.neg(incoming[0]))))
    v=inv(a.fixed(g.neg(g.mul(W['mu'],W['scaled_main_rows'][0]))),
        a.fixed(g.neg(g.sub(W['scaled_main_rows'][1],W['scaled_main_rows'][0]))))
    C=[a.add(p,a.mul(ap,q)) for p,q in zip(u,v)]
    selected=dict(ap=ap,C1=C[0],C2=C[1])
    fields=current.pulse_fields(g,c('7/100'),c('3/10'),c('23/10'),native,selected,W,c('6/5'))
    derivative=current.jets.OrdinaryZDerivative(g)
    for i,node in enumerate(g.nodes):
        if node['operation']=='analytic_unary' and node['name'] in ('original_flat_sigma_prime','compact_raw_beta_derivative'):
            assert node['derivative_source']
            assert derivative.derivative(current.source.FunctionRef(g,i))==g.zero
    evaluate=FixtureEvaluator(g,mu)
    # Fixture-only scalar cache: exact kernel bodies are evaluated first.
    scalar_refs=[W['normalization'],det,W['log_end_scale'],*A,*D,
        *W['incoming_factors'],*W['scaled_main_rows']]
    scalar_bindings={ref.node:evaluate(ref) for ref in scalar_refs}
    norm=quad(beta,-1,1,epsabs=2e-13,epsrel=2e-13)[0];ell=.15
    beta_j=lambda s,center:beta((s-center)/ell)/(ell*norm)
    Yv=13/mu;bp=.5+mu;lams=np.array([.5-mu,.5-2*mu])
    logs=-1/mu-3*(4*mu)**(-2/3)+(2/3)*math.log(4*mu)-.5*math.log(6)-math.log(mu)
    endscale=math.exp(logs-math.log(mu))
    matrix=np.array([[quad(lambda s:math.exp(lam*s)*beta_j(s,center),center-ell,center+ell,
        epsabs=2e-13,epsrel=2e-13)[0] for center in (-3,-1)] for lam in lams])
    main_tail=np.array([sum(quad(lambda x:math.exp(-lam*(13-x)/mu)*gp(x)/mu,l,r,
        epsabs=2e-13,epsrel=2e-13)[0] for l,r in zip([0,.02,10],[.02,10,11])) for lam in lams])
    native_refs={key:sy.lambdify(Z,[sy.diff(expr,Z,n) for n in range(4)],'numpy') for key,expr in expressions.items()}
    ap_ref=sy.lambdify(Z,[sy.diff(ap_expr,Z,n) for n in range(4)],'numpy')
    max_errors=dict(coefficients=0.,U=0.,U_y=0.,m=0.,histories=0.,Q=0.,Q_y=0.,Ur=0.,pressure=0.,joins=0.,terminal_ODE=0.)
    gap_seen=False
    test_points=[('main',0.),('main',.008),('main',.7),('main',10.5),('main',11.),
        ('gap',11.),('gap',12.),('gap',13-4*mu),('end',-4.),('end',-3.),('end',-1.),('end',0.)]
    for zz in (-.7,0.,.65,1.):
        nj={key:np.array(fn(zz),dtype=float) for key,fn in native_refs.items()};aj=np.array(ap_ref(zz),dtype=float)
        mj=np.array([nj['m1'],nj['m2']])
        rhs=-(mj*np.exp(-lams[:,None]*Yv)+main_tail[:,None]*aj)/endscale
        cj=np.linalg.solve(matrix,rhs)
        env={**scalar_bindings,z.node:zz}
        actual_c=np.array([[evaluate(ref,env) for ref in rows(q)] for q in C])
        coefficient_error=float(np.max(np.abs(actual_c-cj)/np.maximum(1,np.abs(cj))))
        max_errors['coefficients']=max(max_errors['coefficients'],coefficient_error)
        assert coefficient_error<2e-10
        def Bjet(y):
            if y<=11/mu:return aj*gp(mu*y)
            s=y-Yv
            return endscale*sum((row*beta_j(s,center) for row,center in zip(cj,(-3,-1))),np.zeros(4))
        def Byjet(y):
            if y<=11/mu:return aj*mu*gp_prime(mu*y)
            s=y-Yv
            return endscale*sum((row*beta_prime((s-center)/ell)/(ell*ell*norm)
                for row,center in zip(cj,(-3,-1))),np.zeros(4))
        boundaries=sorted(set([0,.02/mu,10/mu,11/mu,Yv-4,Yv-3.15,Yv-3,Yv-2.85,
            Yv-1.15,Yv-1,Yv-.85,Yv]))
        pieces=[];state=mj.flatten()
        for lo,hi in zip(boundaries,boundaries[1:]):
            def rhs_ODE(y,flat):return (Bjet(y)[None,:]-lams[:,None]*flat.reshape(2,4)).flatten()
            sol=solve_ivp(rhs_ODE,(lo,hi),state,dense_output=True,rtol=2e-12,atol=2e-13,
                max_step=.025 if lo>=Yv-4 else .25)
            assert sol.success
            state=sol.y[:,-1];pieces.append((lo,hi,sol))
        max_errors['terminal_ODE']=max(max_errors['terminal_ODE'],float(np.max(np.abs(state))))
        assert np.max(np.abs(state))<3e-9
        def history(y):
            if y==0:return mj
            return next(sol.sol(y).reshape(2,4) for lo,hi,sol in pieces if lo-1e-12<=y<=hi+1e-12)
        def checked(pair,expect,environment,scale,key,tol=4e-8):
            actual=np.array([evaluate(ref,environment) for ref in rows(pair)])
            error=float(np.max(np.abs(actual-expect)/scale));max_errors[key]=max(max_errors[key],error)
            assert error<tol,(zz,key,error,actual/scale,expect/scale)
            return actual
        for chart,x in test_points:
            f=fields[chart];y=Yv+x if chart=='end' else x/mu
            environment={**env,f['query'].node:x}
            theta=nj['u']*Pstar*math.exp(-bp*y);scale=theta[0]
            U=convolve(theta,Bjet(y));Uy=convolve(theta,Byjet(y)-bp*Bjet(y))
            means=history(y);m=convolve(theta,means[0]);my=U-m
            for pair,expected in zip(f['normalized_linear_histories'],means):
                checked(pair,expected,environment,1.,'histories')
            checked(f['Utheta'],theta,environment,scale,'U')
            checked(f['U'],U,environment,scale,'U')
            checked(f['U_y'],Uy,environment,scale,'U_y')
            actual_m=checked(f['m'],m,environment,scale,'m')
            checked(f['m_y'],my,environment,scale,'m')
            expected_q=Q_reference(zz,delta,U,m)
            q=checked(f['Q'],expected_q,environment,scale,'Q')
            checked(f['Q_y'],Q_reference(zz,delta,Uy,my),environment,scale,'Q_y')
            radial_scale=scale*math.sqrt(math.exp(logRp+y)/2)
            checked(f['Ur'],expected_q*math.sqrt(math.exp(logRp+y)/2),environment,radial_scale,'Ur')
            mass=-math.expm1(-(1+2*mu)*y)/(1+2*mu)
            pressure=Pstar**2*(nj['P0']+nj['Mp']+.5*mass*convolve(nj['u'],nj['u']))
            checked(f['absolute_swirl_base_pressure'],pressure,environment,Pstar**2,'pressure',2e-10)
            if chart=='gap' and x==12 and zz==.65:
                assert evaluate(f['U'].value,environment)==0
                assert abs(actual_m[0]/scale)>1e-9 and abs(q[0]/scale)>1e-9
                gap_seen=True
            if chart=='end' and x==0:
                for key in ('U','U_y','m','m_y','Q','Q_y','Ur'):
                    assert all(evaluate(ref,environment)==0 for ref in rows(f[key])),key
                assert evaluate(f['absolute_swirl_base_pressure'].value,environment)!=0
        for left,xl,right,xr in [('main',11.,'gap',11.),('gap',13-4*mu,'end',-4.)]:
            lf,rf=fields[left],fields[right];le={**env,lf['query'].node:xl};re={**env,rf['query'].node:xr}
            scale=evaluate(lf['Utheta'].value,le)
            for key in ('Utheta','U','U_y','m','m_y','Q','Q_y'):
                error=max(abs(evaluate(l,le)-evaluate(r,re))/scale for l,r in zip(rows(lf[key]),rows(rf[key])))
                max_errors['joins']=max(max_errors['joins'],error)
                assert error<4e-8,(left,right,key,error)
            for lpair,rpair in zip(lf['normalized_linear_histories'],rf['normalized_linear_histories']):
                error=max(abs(evaluate(l,le)-evaluate(r,re)) for l,r in zip(rows(lpair),rows(rpair)))
                max_errors['joins']=max(max_errors['joins'],error)
                assert error<4e-8,(left,right,'normalized history',error)
    assert gap_seen
    # Radial derivative primitives are checked independently by finite differences.
    for x in (.004,.5,10.3,10.7):
        h=1e-6
        assert abs(gp_prime(x)-(gp(x+h)-gp(x-h))/(2*h))<3e-7
    for x in (-.65,0.,.65):
        h=1e-6
        assert abs(beta_prime(x)-(beta(x+h)-beta(x-h))/(2*h))<2e-9
    return dict(passed=True,mathematical_fixture_only=True,fixture_mu=mu,fixture_delta=delta,
        fixture_Pstar=Pstar,fixture_logRp=logRp,original_parameters_not_replaced=True,
        independent_two_original_moment_matrix_solve=True,
        independent_full_radial_ODE_for_both_mean_histories_Z0_to_Z3=True,
        nonunit_physical_scale_applied_once=True,
        main_gap_end_joins_checked_through_available_Z_orders=True,
        inactive_gap_retains_nonzero_own_mean_and_radial_velocity=True,
        exact_end_empty_support_zero_velocity_mean_Q=True,
        flat_and_beta_derivative_AST_bindings_and_Z_independence_checked=True,
        absolute_pressure_datum_and_swirl_history_retained=True,
        Q_Z0_Z1_Z2_checked_against_independent_ODE_Taylor_quotient=True,
        Q_y_uses_mean_radial_FTC=True,maximum_normalized_errors=max_errors)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1PulseOwnFields(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['own_main_gap_end_fields']==owner.snapshot
    assert owner.jets.acceptance_loaded and not any(raw[key] for key in current.OPEN)
    assert set(owner.fields)=={'main','gap','end'}
    for key,field in owner.fields.items():
        for name in ('Utheta','U','U_y','m','m_y','absolute_swirl_base_pressure'):
            assert isinstance(field[name],current.c3.C3Function)
            owner.g.ids(rows(field[name]))
        for name in ('Q','Q_y','Ur'):
            assert type(field[name]) is current.c2.C2Function
            assert not hasattr(field[name],'ZZZ')
            owner.g.ids(rows(field[name]))
        assert field['pressure_is_not_n1_pressure_forcing']
        assert field['physical_Pstar_applied_exactly_once'] and field['own_history_not_zeroed_on_gap']
    fixtures=dict(own_fields=own_fields_fixture(),source=dict(passed=True,
        selected_actual_C3_ap_and_scaled_end_functions_consumed=True,
        original_native_and_parent_graph_prefix_preserved=True,
        exact_nonzero_end_scale_retained=True,missing_fourth_mean_derivative_not_zero_padded=True,
        all_absolute_cumulative_histories_retained=True))
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=fixtures,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        uniform_original_field_numerical_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_ACTUAL_OWN_PULSE_FIELDS_CHECKED',flush=True)
    print(json.dumps(fixtures['own_fields']['maximum_normalized_errors']),flush=True)
    return receipt


if __name__=='__main__':run()
