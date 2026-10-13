"""Independent collar ODE/history and physical pressure checks."""
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions as current


def symbolic(g,ref,cache=None):
    cache={} if cache is None else cache
    if ref.node in cache:return cache[ref.node]
    row=g.nodes[ref.node];op=row['operation']
    child=lambda i:symbolic(g,current.source.FunctionRef(g,i),cache)
    if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
    elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
    elif op=='sum':value=sy.Add(*(child(i) for i in row['arguments']))
    elif op=='product':value=sy.Mul(*(child(i) for i in row['arguments']))
    elif op=='negative':value=-child(row['argument'])
    elif op=='positive_quotient':value=sy.cancel(child(row['numerator'])/child(row['denominator']))
    elif op=='analytic_unary':value=getattr(sy,row['name'])(child(row['argument']))
    elif op=='function_substitution':value=child(row['expression']).subs(child(row['variable']),child(row['value']))
    elif op=='definite_integral':
        body=sy.expand(child(row['integrand']));variable=sy.Symbol(row['variable'],real=True)
        value=sy.integrate(body,(variable,child(row['lower']),child(row['upper'])),conds='none')
    else:raise AssertionError('A source/cap node is not a fixture value: '+op)
    cache[ref.node]=value
    return value


def fixture(exponential=False):
    g=current.source.FunctionTransportGraph();a=current.c3.C3Algebra(g)
    z,r,s=sy.symbols('Z fixture_core_rho original_lower_collar_s',real=True)
    variables={str(v):g.symbol(str(v)) for v in (z,r,s)}
    bind=lambda p:current.source.expression(g,sy.sympify(p),variables)
    jet=lambda p:current.c3.C3Function(*(bind(sy.diff(p,z,k)) for k in range(4)))
    hb=sy.Rational(1,5);eps=sy.Rational(1,9);delta=sy.Rational(1,5);sc=sy.Rational(3,5)
    b=(3+z)/7;alpha=sy.Rational(1,7)+z/11
    Phi=sy.exp(alpha*r) if exponential else 2+z*z/10
    Psi=0 if exponential else b*r*r
    mean=0 if exponential else b*r*r/3
    rows={key:[jet(sy.diff(expr,r,i)) for i in range(3)]
        for key,expr in (('Phi',Phi),('Psi',Psi),('mean_Psi',mean))}
    params={key:bind(value) for key,value in dict(hb=hb,epsilon=eps,delta=delta,
        logLambda=sy.log(1/eps),sc=sc,sc_quarter=sc/4,sc_half=sc/2).items()}
    functions=dict(U0=jet(4*z+sy.Rational(1,10)),F0=jet(1+z*z/20),P0=jet(sy.Rational(3,2)-z*z/4))
    sigma=(lambda x:x) if exponential else (lambda x:g.constant('1/3'))
    sigma_prime=(lambda x:g.one) if exponential else (lambda x:g.zero)
    fields=current.collar_fields(g,variables[str(r)],rows,functions,params,sigma,sigma_prime)
    return g,fields,params,dict(z=z,r=r,s=s,hb=hb,epsilon=eps,delta=delta,sc=sc,b=b,alpha=alpha,
        Phi=Phi,functions=functions)


def ratio_fixture():
    g,f,p,v=fixture(exponential=True);z,r,s,hb=v['z'],v['r'],v['s'],v['hb'];alpha=v['alpha']
    rs=4*sy.exp(hb*s)
    logratio=4*(hb-1)*alpha*((s-1/hb)*sy.exp(hb*s)+1/hb)
    # Derive rho rows from the independently integrated exponential ODE.
    phase=sy.log(r/4)/hb
    exact=sy.exp(alpha*r+4*(hb-1)*alpha*((phase-1/hb)*r/4+1/hb))
    cache={}
    for k,ref in enumerate(current.c3.rows(f['logratio'])):
        assert sy.simplify(symbolic(g,ref,cache)-sy.diff(logratio,z,k))==0
    for i,jet in enumerate(f['Phi_rho_rows']):
        for k,ref in enumerate(current.c3.rows(jet)):
            want=sy.diff(exact,r,i,z,k).subs(r,rs)
            assert sy.simplify(symbolic(g,ref,cache)-want)==0,(i,k)
    for k,ref in enumerate(current.c3.rows(f['ratio'])):
        assert sy.simplify(symbolic(g,ref,cache).subs(s,0)-(1 if k==0 else 0))==0
    return dict(passed=True,mathematical_fixture_only=True,nonconstant_phase_sigma_fixture=True,
        independent_integrated_logarithmic_ODE=True,nonzero_chi_rho_and_C3_rows_checked=True,
        rho4_ratio_inlet_exact_one=True)


def history_and_pressure_fixture():
    g,f,p,v=fixture();z,r,s=v['z'],v['r'],v['s'];hb,eps,de,b=v['hb'],v['epsilon'],v['delta'],v['b']
    chi=1+(hb-1)/3;rs=4*sy.exp(hb*s)
    psi=b*(16+chi*(r*r-16))
    mean=b*(chi*r*r/3+16*(1-chi)-128*(1-chi)/(3*r))
    U0=4*z+sy.Rational(1,10);U=U0+eps*psi;m=U0+eps*mean
    d=1-z*z;L=1-de*z*z
    Q=sy.cancel((2*z*U-(1-de)*z*m-d*sy.diff(m,z))/L)
    cache={}
    for key,expr in (('Psi_rho_rows',psi),('mean_Psi_rho_rows',mean),('U_rho_rows',U),('m_rho_rows',m)):
        for i,jet in enumerate(f[key]):
            for k,ref in enumerate(current.c3.rows(jet)):
                assert sy.cancel(symbolic(g,ref,cache)-sy.diff(expr,r,i,z,k).subs(r,rs))==0,(key,i,k)
    for i,jet in enumerate(f['Q_rho_rows']):
        for k,ref in enumerate((jet.value,jet.Z,jet.ZZ)):
            assert sy.cancel(symbolic(g,ref,cache)-sy.diff(Q,r,i,z,k).subs(r,rs))==0,('Q',i,k)
    assert sy.cancel(r*sy.diff(mean,r)+mean-psi)==0
    assert sy.cancel(mean.subs(r,4)-16*b/3)==0
    qy=r*sy.diff(Q,r);qyy=r*sy.diff(qy,r)
    for key,expr in (('Q_y',qy),('Q_yy',qyy)):
        for k,ref in enumerate((f[key].value,f[key].Z,f[key].ZZ)):
            assert sy.cancel(symbolic(g,ref,cache)-sy.diff(expr,z,k).subs(r,rs))==0,(key,k)
    base=(sy.Rational(3,2)-z*z/4)+eps*(1+z*z/20)**2*(2+z*z/10)**2*r
    for k,ref in enumerate(current.c3.rows(f['absolute_swirl_base_pressure'])):
        assert sy.cancel(symbolic(g,ref,cache)-sy.diff(base,z,k).subs(r,rs))==0,('base pressure',k)
    pressure,bindings=current.direct_pressure_functions(g,f,p)
    # Independent physical-R differential operator, including both viscous
    # derivatives. It does not consume the producer's reduced S or measures.
    R=sy.Symbol('physical_fixture_R',positive=True)
    qr=Q.subs(r,R/eps);ur=U.subs(r,R/eps)
    N=sy.cancel(-((qr+(1-de)*z*sy.diff(qr,z)/2+R*sy.diff(qr,R))/L
        +qr*(qr/2+R*sy.diff(qr,R))
        +ur*(-2*z*qr+d*sy.diff(qr,z)-2*z*R*sy.diff(qr,R))/L)/2
        +2*sy.diff(qr,R)+R*sy.diff(qr,R,2))
    lower=4*eps*sy.exp(hb*v['sc']/4);upper=4*eps*sy.exp(hb*v['sc']/2)
    independent={key:sign*sy.integrate(sy.expand(sy.cancel(R**power*N)),(R,lower,upper),conds='none')
        for key,power,sign in (('m13',0,1),('m15',1,-1))}
    for key,pair in pressure.items():
        for k,ref in enumerate((pair.value,pair.Z)):
            got=symbolic(g,ref,cache);want=sy.diff(independent[key],z,k)
            # Every exponential endpoint is an integer power of exp(3/100).
            # Treating those as one exact generator avoids an expensive
            # generic transcendental simplifier while preserving equality.
            E=sy.Symbol('exact_exp_three_over_one_hundred',positive=True)
            normal=lambda x:sy.expand_log(x,force=True).replace(
                lambda q:q.func==sy.exp and q.args[0].is_Rational,
                lambda q:E**(q.args[0]*100/3))
            assert sy.cancel(normal(got)-normal(want))==0,(key,k)
    assert sy.cancel(symbolic(g,bindings['N1p'].value,cache)-N.subs(R,eps*rs))==0
    return dict(passed=True,mathematical_fixture_only=True,epsilon=str(eps),hb=str(hb),delta=str(de),
        independent_own_mean_and_rho4_history=True,ordinary_C3_and_Q_C2_checked=True,
        core_plus_collar_swirl_pressure_history_checked=True,
        both_signed_pressure_C1_equal_direct_physical_R_integrals=True,
        both_viscous_radial_derivatives_retained=True,original_m15_minus_applied_once=True,
        nonunit_physical_measures_checked=True)


def source_fixture(owner):
    g=owner.g;f=owner.fields;p=owner.parameters;c=g.constant
    assert owner.geometry_certificate['sc_exact_mpf_tuple']==owner.sc_bits
    assert p['sc']==current.exact_binary_scalar(g,owner.sc_bits)
    assert p['sc_quarter']==g.mul(c('1/4'),p['sc']) and p['sc_half']==g.mul(c('1/2'),p['sc'])
    assert owner.geometry_certificate['original_modification_support']['increments_exact_zero']
    assert g.nodes[p['hb'].node]['name']=='hbB'
    assert f['radius']==g.mul(p['epsilon'],f['rho'])
    def factors(ref):
        if ref.node in (p['epsilon'].node,p['hb'].node,f['rho'].node):return [ref.node]
        n=g.nodes[ref.node]
        assert n['operation']=='product'
        return [j for i in n['arguments'] for j in factors(current.source.FunctionRef(g,i))]
    assert sorted(factors(f['physical_measure']))==sorted([p['epsilon'].node,p['hb'].node,f['rho'].node])
    assert sorted(factors(f['weighted_physical_measure']))==sorted([p['epsilon'].node,p['epsilon'].node,p['hb'].node,f['rho'].node,f['rho'].node])
    for quantity,jets in owner.core_radial_rows.items():
        for i,jet in enumerate(jets):
            for k,ref in enumerate(current.c3.rows(jet)):
                n=g.nodes[ref.node]
                assert n['operation']=='current_original_unique_analytic_core_component'
                assert n['quantity']==quantity and n['radial_order']==i and n['ordinary_slow_Z_derivative_order']==k
                assert n['unique_analytic_Banach_solution'] and n['source_family']==owner.family
    for pair in owner.functions.values():
        for k,ref in enumerate((pair.value,pair.Z)):
            n=g.nodes[ref.node]
            assert n['operation']=='definite_integral' and n['ordinary_slow_Z_derivative_order']==k
            assert n['lower']==p['sc_quarter'].node and n['upper']==p['sc_half'].node
            assert n['signed_integrand_not_a_cap'] and n['all_finite_viscous_derivatives_retained']
    assert f['all_histories_start_at_true_core_rho4']
    return dict(passed=True,checked_true_core_nodes_only=True,exact_original_binary_sc_retained=True,
        original_nonzero_hb_and_positive_physical_measures_retained=True,
        exact_Z_independent_endpoints_and_original_zero_increment_support=True,
        genuine_signed_value_and_first_Z_integral_bodies=True,
        no_enclosure_center_or_comparison_selected_as_field=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1LowerCollarFunctions(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['lower_collar_bindings']==owner.snapshot
    checks=dict(actual_source=source_fixture(owner))
    print('ACTUAL_COLLAR_SOURCE_AND_MEASURES_CHECKED',flush=True)
    checks['ratio']=ratio_fixture()
    print('INDEPENDENT_COLLAR_RATIO_C3_CHECKED',flush=True)
    checks['history_and_physical_pressure']=history_and_pressure_fixture()
    print('INDEPENDENT_COLLAR_HISTORY_AND_PRESSURE_C1_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        numerical_original_collar_point_accuracy_claimed=False,
        complete_inner_cutoff_five_defect_adapter_installed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_TRUE_LOWER_COLLAR_C1_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
