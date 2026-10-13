"""Independent manufactured solution, exact regular axis and physical units."""
import copy
import json
from pathlib import Path
import time
import sympy as sy
from sympy.polys.fields import field

import lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions as current


def exact_fixture(g,z,r):
    X=sy.Symbol('radial_fixture_X',positive=True)
    K,Zf,Xf=field((z,X),sy.QQ);cache={}
    def at(ref):
        if ref.node in cache:return cache[ref.node]
        n=g.nodes[ref.node];op=n['operation'];get=lambda i:at(current.source.FunctionRef(g,i))
        if op=='exact_rational':v=K(n['numerator'])/K(n['denominator'])
        elif op=='bound_variable':
            assert n['name'] in (str(z),str(r)),n['name']
            v=Zf if n['name']==str(z) else Xf*Xf
        elif op=='sum':v=sum((get(i) for i in n['arguments']),K.zero)
        elif op=='product':
            v=K.one
            for i in n['arguments']:v*=get(i)
        elif op=='negative':v=-get(n['argument'])
        elif op=='positive_quotient':v=get(n['numerator'])/get(n['denominator'])
        elif op=='analytic_unary':
            assert n['name']=='sqrt' and get(n['argument'])==Xf*Xf
            v=Xf
        elif op=='function_substitution':
            assert g.nodes[n['variable']]['name']==str(r)
            target=get(n['value']);assert target==K.zero
            v=K.from_expr(get(n['expression']).as_expr().subs(X,0))
        else:raise AssertionError('Not an independent rational fixture: '+op)
        cache[ref.node]=v;return v
    want=lambda q:K.from_expr(sy.sympify(q).subs(r,X*X))
    return at,want


def manufactured_solution_fixture():
    g=current.source.FunctionTransportGraph();z=sy.Symbol('Z',real=True);r=sy.Symbol('radial_test_rho',nonnegative=True)
    X=sy.Symbol('radial_compile_X',positive=True)
    variables={str(z):g.symbol(str(z)),str(X):g.unary('sqrt',g.symbol(str(r)))}
    bind=lambda q:current.source.expression(g,sy.sympify(q).subs(r,X*X),variables)
    jet=lambda q:current.c3.C3Function(*(bind(sy.diff(q,z,k)) for k in range(4)))
    eps=sy.Rational(1,9);de=sy.Rational(1,5);L=1-de*z*z
    leadingF=1+z*z+r*(1+z)/7+r*r/13
    leadingM=2+z+z*z+r*(1-2*z+z**3)
    leadingU=leadingM+r*sy.diff(leadingM,r)
    leadingQ=sy.cancel((2*z*leadingU-(1-de)*z*leadingM-(1-z*z)*sy.diff(leadingM,z))/L)
    leading=dict(F_rho_rows=[jet(sy.diff(leadingF,r,i)) for i in range(3)],
        U_rho_rows=[jet(sy.diff(leadingU,r,i)) for i in range(3)],
        Q_rho_rows=[jet(sy.diff(leadingQ,r,i)) for i in range(3)])
    F=(1+z)*r+(2+z*z)*r*r+(3+z**3)*r**3
    U=(2-z+z**3)*r+(1+z*z)*r*r+(2+z)*r**3
    P=(1+2*z+z**3)*r+(3+z*z)*r*r
    mean=sum(sy.expand(U).coeff(r,n)*r**n/(n+1) for n in range(1,4));K1=mean-U
    expressions=dict(F1=F,Uz1=U,K1=K1,P1=P,d_x_F1=2*sy.sqrt(r)*sy.diff(F,r),d_x_Uz1=2*sy.sqrt(r)*sy.diff(U,r))
    state={name:jet(q) for name,q in expressions.items()}
    values=dict(rho=r,x=sy.sqrt(r),Z=z,delta=de,epsilon=eps,F=leadingF,Fr=sy.diff(leadingF,r),
        FZ=sy.diff(leadingF,z),U=leadingU,Ur=sy.diff(leadingU,r),UZ=sy.diff(leadingU,z),Q=leadingQ)
    matrices={label:{name:q.subs({sy.Symbol(k):v for k,v in values.items()}) for name,q in rows.items()}
        for label,rows in current.inner.raw_matrix_expressions().items()}
    lhs={'5':8*sy.diff(F,r)+4*r*sy.diff(F,r,2),'6':4*sy.diff(U,r)+4*r*sy.diff(U,r,2)}
    forcing={}
    for row in ('5','6'):
        q=lhs[row]
        for label in ('B0','B1'):
            for key,B in matrices[label].items():
                i,j=key.split(',')
                if i!=row:continue
                w=expressions[current.inner.FIELDS[int(j)-1]]
                q-=B*(sy.diff(w,z) if label=='B1' else w)
        forcing[row]=sy.cancel(q)
    Np=sy.diff(P,r)/eps-2*leadingF*F
    system=dict(**{label:{name:jet(q) for name,q in rows.items()} for label,rows in matrices.items()},
        g={**{row:jet(q) for row,q in forcing.items()},'4':jet(2*eps*sy.sqrt(r)*Np)},
        N1p=jet(Np),rho=g.symbol(str(r)),x=g.unary('sqrt',g.symbol(str(r))))
    zd=current.inner.AnalyticZDerivative(g);rd=current.CoreRhoDerivative(g,system['rho'])
    velocity,H=current.equation_velocity_rows(g,system,state,zd)
    rows={**velocity,**current.own_mean_rows(g,state,velocity,system['rho'])}
    axis,axis_binding=current.regular_axis_rows(g,system,leading,dict(epsilon=bind(eps),delta=bind(de)),rd,zd)
    recovered=current.recovered_radial_rows(g,rows,system['rho'],bind(eps),bind(de))
    pressure=current.pressure_radial_function(g,jet(leadingF),state['F1'],system['N1p'],bind(eps))
    exact,want=exact_fixture(g,z,r);count=0
    for name,q in (('F1',F),('Uz1',U),('mean1',mean),('K1',K1)):
        for i in range(3):
            for k,ref in enumerate(current.c3.rows(rows[name][i])):
                expected=sy.diff(q,r,i,z,k)
                assert exact(ref)==want(expected),(name,i,k,'interior')
                assert exact(current.c3.rows(axis[name][i])[k])==want(expected.subs(r,0)),(name,i,k,'axis')
                count+=2
    Q=((1-de)*z*U-(1+de)*z*K1-(1-z*z)*sy.diff(mean,z))/L;V=eps*r*Q
    for name,q in (('Q1',Q),('V1',V)):
        for i in range(3):
            for k,ref in enumerate((recovered[name][i].value,recovered[name][i].Z,recovered[name][i].ZZ)):
                assert exact(ref)==want(sy.diff(q,r,i,z,k)),(name,i,k)
                count+=1
    for k,ref in enumerate(current.c3.rows(pressure)):
        assert exact(ref)==want(sy.diff(P,r,z,k)),('pressure',k)
        assert exact(current.c3.rows(axis['P1_rho'])[k])==want(sy.diff(P,r,z,k).subs(r,0))
        count+=2
    # A radial derivative must preserve ordinary Z order and increase only
    # the radial order of an actual leading analytic component.
    component=g.node('current_original_unique_analytic_core_component',rho=system['rho'].node,
        radial_order=2,ordinary_slow_Z_derivative_order=3,quantity='Psi')
    n=g.nodes[rd.derivative(component).node]
    assert n['radial_order']==3 and n['ordinary_slow_Z_derivative_order']==3
    assert rd.derivative(g.symbol('Z'))==g.zero
    return dict(passed=True,mathematical_fixture_only=True,exact_mixed_derivative_identities=count,
        full_nonzero_raw_B0_B1_manufactured_feedback=True,independent_cubic_regular_solution=True,
        exact_first_second_axis_limits=True,nonunit_epsilon_delta_and_pressure_units=True,
        same_own_mean_K_Q_V_and_n1_weight=True,no_odd_axis_rho_derivative=True)


def source_fixture(owner):
    g=owner.g
    for chart,rows in owner.rows.items():
        assert set(rows)=={'F1','Uz1','mean1','K1'}
        assert all(len(jets)==3 for jets in rows.values())
        assert owner.definition['full_raw_B0_and_B1_dZ_retained']
        if chart=='core':
            for key,jets in rows.items():
                for i in (1,2):
                    for k,ref in enumerate(current.c3.rows(jets[i])):
                        n=g.nodes[ref.node]
                        assert n['operation']=='regular_even_radial_axis_continuation'
                        assert n['axis_expression']==current.c3.rows(owner.axis[key][i])[k].node
                        assert n['interior_expression']==current.c3.rows(owner.interior_rows[chart][key][i])[k].node
                        assert n['ordinary_rho_derivative_order']==i and n['ordinary_slow_Z_derivative_order']==k
        for name in ('Q1','V1'):assert len(owner.recovered[chart][name])==3
    assert not owner.definition['pressure_second_radial_function_installed']
    assert not owner.definition['original_numeric_n1_solution_oracle_executed']
    extra=g.nodes[len(owner.inner.prefix):]
    assert any(n['operation']=='current_original_unique_analytic_core_component' and n['radial_order']>=3 for n in extra)
    assert owner.inner.definition['common_inner_endpoint_is_Rin_not_rminus']
    return dict(passed=True,true_complete_regular_solution_used=True,
        source_bound_core_radial_derivatives_not_error_envelopes=True,
        explicit_regular_axis_functions_and_positive_rho_bodies=True,
        original_common_inner_Rin_scope_preserved=True,pressure_first_derivative_only=True,
        numerical_original_point_oracle_not_claimed=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1InnerRadialFunctions(require_checked=False)
    data=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert data['exact_graph_nodes']==owner.g.nodes and data['inner_radial_bindings']==owner.snapshot
    checks=dict(source=source_fixture(owner));print('TRUE_N1_RADIAL_SOURCE_CHECKED',flush=True)
    checks['manufactured']=manufactured_solution_fixture();print('EXACT_REGULAR_RADIAL_AND_AXIS_FIXTURE_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,**{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        original_numeric_n1_solution_oracle_executed=False,full_inner_cutoff_and_five_defect_adapter_installed=False,
        pressure_second_radial_function_installed=False,numerical_original_n1_point_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_TRUE_INNER_RADIAL_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
