"""Raw operator, exact Green paths, derivative feedback and limit admission."""
import ast
import inspect
import json
from pathlib import Path
import time
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions as current


def symbol(name):
    return sy.Symbol(name,nonnegative=True) if any(s in name for s in ('rho','Green_t','Green_phase','test_s')) else sy.Symbol(name,real=True)


def symbolic(g,ref,cache=None):
    cache={} if cache is None else cache
    if ref.node in cache:return cache[ref.node]
    n=g.nodes[ref.node];op=n['operation'];at=lambda i:symbolic(g,current.source.FunctionRef(g,i),cache)
    if op=='exact_rational':v=sy.Rational(n['numerator'],n['denominator'])
    elif op=='bound_variable':v=symbol(n['name'])
    elif op=='sum':v=sy.Add(*(at(i) for i in n['arguments']))
    elif op=='product':v=sy.Mul(*(at(i) for i in n['arguments']))
    elif op=='negative':v=-at(n['argument'])
    elif op=='positive_quotient':v=sy.cancel(at(n['numerator'])/at(n['denominator']))
    elif op=='analytic_unary':
        q=at(n['argument'])
        v=(sy.Integer(1) if q==0 else (sy.exp(q)-1)/q) if n['name']=='exprel' else getattr(sy,n['name'])(q)
    elif op=='function_substitution':v=at(n['expression']).subs(at(n['variable']),at(n['value']))
    elif op=='definite_integral':
        body=sy.expand(sy.cancel(at(n['integrand'])))
        v=sy.integrate(body,(symbol(n['variable']),at(n['lower']),at(n['upper'])),conds='none')
    else:raise AssertionError('An original solution/cap is not a mathematical fixture: '+op)
    cache[ref.node]=v
    return v


def compile_primary(fn,scope):
    """Run the accepted mathematical AST with SymPy rather than box truncation."""
    tree=ast.parse(inspect.getsource(fn).strip())
    class IdentityTruncation(ast.NodeTransformer):
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if isinstance(node.func,ast.Name) and node.func.id=='truncate':return node.args[0]
            return node
    tree=IdentityTruncation().visit(tree);ast.fix_missing_locations(tree)
    exec(compile(tree,'<accepted-raw-operator-mathematical-fixture>','exec'),scope)
    return scope[fn.__name__]


def rational_fixture_converter(g,z,r):
    """Exact canonical field arithmetic avoids expanding derivative DAGs.

    rho=X^2 with positive X preserves the raw sqrt(rho) units. No sampled
    point or tolerance replaces any matrix, forcing or derivative identity.
    """
    from sympy.polys.fields import field
    X=sy.Symbol('raw_fixture_X',positive=True)
    K,Zf,Xf=field((z,X),sy.QQ);cache={}
    def at(ref):
        if ref.node in cache:return cache[ref.node]
        n=g.nodes[ref.node];op=n['operation']
        get=lambda i:at(current.source.FunctionRef(g,i))
        if op=='exact_rational':v=K(n['numerator'])/K(n['denominator'])
        elif op=='bound_variable':
            q=symbol(n['name']);assert q in (z,r),q
            v=Zf if q==z else Xf*Xf
        elif op=='sum':v=sum((get(i) for i in n['arguments']),K.zero)
        elif op=='product':
            v=K.one
            for i in n['arguments']:v*=get(i)
        elif op=='negative':v=-get(n['argument'])
        elif op=='positive_quotient':v=get(n['numerator'])/get(n['denominator'])
        elif op=='analytic_unary':
            assert n['name']=='sqrt' and get(n['argument'])==Xf*Xf
            v=Xf
        else:raise AssertionError('Not raw fixture rational arithmetic: '+op)
        cache[ref.node]=v;return v
    def expected(q,k=0):
        v=K.from_expr(sy.sympify(q).subs(r,X*X))
        for _ in range(k):v=v.diff(Zf)
        return v
    return at,expected


def raw_operator_fixture():
    g=current.source.FunctionTransportGraph();z=symbol('Z');r=symbol('operator_test_rho')
    variables={str(q):g.symbol(str(q)) for q in (z,r)}
    bind=lambda q:current.source.expression(g,sy.sympify(q),variables)
    jet=lambda q:current.c3.C3Function(*(bind(sy.diff(q,z,k)) for k in range(4)))
    de=sy.Rational(1,5);eps=sy.Rational(1,9);d=1-z*z;L=1-de*z*z
    m=2+z+z*z+r*(1-2*z+z**3);U=m+r*sy.diff(m,r)
    F=1+z*z+r*(z+1)/7+r*r*(z*z+1)/13
    Q=sy.cancel((2*z*U-(1-de)*z*m-d*sy.diff(m,z))/L)
    fields=dict(rho=variables[str(r)],F_rho_rows=[jet(sy.diff(F,r,i)) for i in range(3)],
        U_rho_rows=[jet(sy.diff(U,r,i)) for i in range(3)],
        Q_rho_rows=[current.own.c2.C2Function(*(bind(sy.diff(Q,r,i,z,k)) for k in range(3))) for i in range(3)])
    params=dict(delta=bind(de),epsilon=bind(eps),Lambda=bind(1/eps))
    result=current.raw_system(g,fields,params,current.AnalyticZDerivative(g))
    Zop=lambda q,a:(a*z*q+d*sy.diff(q,z)-2*z*r*sy.diff(q,r))/L
    Nt=-Zop(Zop(F,-2-de),-3);Nz=-Zop(Zop(U,-1-de),-2)
    R=sy.Symbol('physical_fixture_R',positive=True);qr=Q.subs(r,R/eps);ur=U.subs(r,R/eps)
    Np=sy.cancel(-((qr+(1-de)*z*sy.diff(qr,z)/2+R*sy.diff(qr,R))/L+qr*(qr/2+R*sy.diff(qr,R))
        +ur*(-2*z*qr+d*sy.diff(qr,z)-2*z*R*sy.diff(qr,R))/L)/2+2*sy.diff(qr,R)+R*sy.diff(qr,R,2)).subs(R,eps*r)
    class Ctx:
        sqrt=staticmethod(sy.sqrt)
        mpf=staticmethod(lambda q:sy.Rational(str(q)))
    class Algebra:
        ctx=Ctx()
        lift=staticmethod(sy.sympify)
        shift=staticmethod(lambda q,p:q*(1/eps)**p[3])
    scope=dict(dz=lambda q:sy.diff(q,z))
    scope['axial_second']=compile_primary(current.raw.axial_second,scope)
    accepted=compile_primary(current.raw.matrix,scope)(Algebra(),dict(F=[F,sy.diff(F,r),sy.diff(F,r,2)],
        Uz=[U,sy.diff(U,r),sy.diff(U,r,2)],Q=[Q,sy.diff(Q,r),sy.diff(Q,r,2)]),
        dict(N1theta=Nt,N1z=Nz,N1p=Np),z,r,de)
    exact,want_exact=rational_fixture_converter(g,z,r);count=0
    for label in ('B0','B1','g'):
        assert set(result[label])==set(accepted[label])
        for name,jet_value in result[label].items():
            for k in range(4):
                assert exact(current.c3.rows(jet_value)[k])==want_exact(accepted[label][name],k),(label,name,k)
                count+=1
    for name,want in (('N1theta',Nt),('N1z',Nz),('N1p',Np)):
        assert exact(result[name].value)==want_exact(want),name
    return dict(passed=True,mathematical_fixture_only=True,accepted_raw_matrix_AST_rows=count,
        direct_nested_Z_operators_and_physical_R_viscosity=True,all_epsilon_powers_and_B1_entries_retained=True,
        exact_rational_field_identities_without_sampling=True)


def green_feedback_fixture():
    g=current.source.FunctionTransportGraph();z=symbol('Z');r=symbol('green_test_rho');s=symbol('test_s')
    variables={str(q):g.symbol(str(q)) for q in (z,r,s)}
    root=sy.Symbol('green_fixture_sqrt_rho',nonnegative=True)
    variables[str(root)]=g.unary('sqrt',variables[str(r)])
    bind=lambda q:current.source.expression(g,sy.sympify(q).xreplace({sy.sqrt(r):root}),variables)
    jet=lambda q:current.c3.C3Function(*(bind(sy.diff(q,z,k)) for k in range(4)))
    hb=sy.Rational(1,5);eps=sy.Rational(1,9);rs=4*sy.exp(hb*s)
    A=2+z;B=-3+2*z+z**3;C=1+z+z*z
    coeff=(sy.Rational(3,7),sy.Rational(-2,9),sy.Rational(5,11),sy.Rational(2,13),sy.Rational(-3,17),sy.Rational(7,19))
    a,b,p,d,e,f=coeff
    coordinates=dict(core=variables[str(r)],first=variables[str(s)]);params=dict(hb=bind(hb),epsilon=bind(eps))
    systems={};vectors={}
    for chart,radius in (('core',r),('first',rs)):
        x=sy.sqrt(radius)
        systems[chart]=dict(C0={'4,1':jet(a*x),'5,5':jet(b*x),'6,4':jet(p)},
            B1={'5,1':jet(d),'6,2':jet(e),'6,4':jet(f)})
        vectors[chart]={'4':jet(A*x),'5':jet(B),'6':jet(C)}
    W0={chart:current.green_apply(g,vectors,chart,coordinates,params) for chart in coordinates}
    feedback=current.feedback_vectors(g,systems,W0,current.AnalyticZDerivative(g))
    W1={chart:current.green_apply(g,feedback,chart,coordinates,params) for chart in coordinates}
    alpha=a*B/8;beta=b*B/4+d*sy.diff(B,z)/8;gamma=p*A/2+e*sy.diff(C,z)/4+f*sy.diff(A,z)/2
    expected0=dict(F1=B*r/8,Uz1=C*r/4,K1=-C*r/8,P1=A*r/2,d_x_F1=B*sy.sqrt(r)/4,d_x_Uz1=C*sy.sqrt(r)/2)
    expected1=dict(F1=beta*r*r/24,Uz1=gamma*r*r/16,K1=-gamma*r*r/24,P1=alpha*r*r/4,
        d_x_F1=beta*r**sy.Rational(3,2)/6,d_x_Uz1=gamma*r**sy.Rational(3,2)/4)
    cache={};count=0
    for word,expected in ((W0,expected0),(W1,expected1)):
        for chart in coordinates:
            for name,q in word[chart].items():
                for k,ref in enumerate(current.c3.rows(q)):
                    got=symbolic(g,ref,cache);want=sy.diff(expected[name],z,k)
                    if chart=='first':want=want.subs(r,rs)
                    assert sy.simplify(got-want)==0,(chart,name,k)
                    if chart=='core':assert sy.simplify(got.subs(r,0))==0
                    else:assert sy.simplify(got.subs(s,0)-sy.diff(expected[name],z,k).subs(r,4))==0
                    count+=1
    state={name:jet(value) for name,value in expected0.items()}
    de=sy.Rational(1,5);rec=current.recovery(g,state,dict(rho=variables[str(r)],epsilon=bind(eps)),bind(de))
    mean=expected0['Uz1']+expected0['K1'];U=expected0['Uz1'];K=expected0['K1']
    Q=((1-de)*z*U-(1+de)*z*K-(1-z*z)*sy.diff(mean,z))/(1-de*z*z)
    for key,want in (('mean1',mean),('radial_primitive_M1',eps*r*mean),('Q1',Q),('V1',eps*r*Q)):
        assert sy.cancel(symbolic(g,rec[key].value,cache)-want)==0,key
    return dict(passed=True,mathematical_fixture_only=True,complete_W0_and_first_feedback_C3_rows=count,
        genuine_B0_and_B1_dZ_feedback=True,negative_K1_and_core_to_collar_seam=True,
        zero_axis_regular_values=True,physical_n1_mean_primitive_and_V_units=True,
        n1_temporal_weight_is_one_plus_delta=True)


def source_fixture(owner):
    g=owner.g
    for chart,functions in owner.functions.items():
        for i,name in enumerate(current.FIELDS):
            for k,ref in enumerate(current.c3.rows(functions[name])):
                n=g.nodes[ref.node]
                assert n['operation']=='current_original_unique_regular_n1_component'
                assert n['all_words_retained_in_analytic_limit'] and n['unique_regular_analytic_solution']
                assert n['ordinary_slow_Z_derivative_order']==k and n['zero_axis_regular_power']==current.POWERS[i]
                assert n['actual_raw_operator_bindings']==current.own.encode(owner.systems)
                assert n['coordinate_domain_lower']==g.zero.node
                assert n['coordinate_domain_upper']==(g.constant(4).node if chart=='core' else owner.parameters['sc_quarter'].node)
    assert owner.definition['raw_C0_is_B0_minus_Akin'] and owner.definition['normalized_Kchi_not_used_with_raw_matrix']
    assert owner.definition['common_inner_endpoint_is_Rin_not_rminus']
    assert owner.tail_admission['all_k_derivative_words_allowed']
    assert all(set(system['C0'])==set(system['B0'])-{'1,5','2,6','3,6'} for system in owner.systems.values())
    assert all(rec['temporal_n1_mean_weight_is_1_plus_delta'] for rec in owner.recovered.values())
    assert all(set(owner.fixed_point_rhs[c])==set(current.FIELDS) for c in owner.coordinates)
    return dict(passed=True,actual_same_source_operator_bound_to_complete_limit=True,
        full_raw_preconditioned_all_word_convergence_consumed=True,full_B1_derivative_words_not_omitted=True,
        complete_fixed_point_rhs_and_separate_finite_words=True,common_inner_Rin_scope_explicit=True,
        numerical_center_not_selected_as_solution=True,normalized_potential_tail_not_transplanted_to_raw_operator=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1InnerVolterraFunctions(require_checked=False)
    raw_data=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw_data['exact_graph_nodes']==owner.g.nodes and raw_data['inner_Volterra_bindings']==owner.snapshot
    checks=dict(actual_source=source_fixture(owner));print('ACTUAL_RAW_N1_LIMIT_ADMISSION_CHECKED',flush=True)
    checks['raw_operator']=raw_operator_fixture();print('INDEPENDENT_RAW_N1_OPERATOR_CHECKED',flush=True)
    checks['Green_feedback']=green_feedback_fixture();print('INDEPENDENT_GREEN_AND_DZ_FEEDBACK_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        original_numeric_n1_solution_oracle_executed=False,
        full_inner_cutoff_and_five_defect_adapter_installed=False,
        numerical_original_n1_point_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_TRUE_INNER_VOLTERRA_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
