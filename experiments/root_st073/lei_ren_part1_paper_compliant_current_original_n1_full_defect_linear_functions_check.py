"""Full signed partition and universal physical C1 linear cancellation."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator_check as linear_check


def symbolic(g,ref,bindings,cache):
    if ref.node in bindings:return bindings[ref.node]
    if ref.node in cache:return cache[ref.node]
    n=g.nodes[ref.node];op=n['operation'];at=lambda i:symbolic(g,current.source.FunctionRef(g,i),bindings,cache)
    if op=='exact_rational':v=sy.Rational(n['numerator'],n['denominator'])
    elif op=='sum':v=sy.Add(*(at(i) for i in n['arguments']))
    elif op=='product':v=sy.Mul(*(at(i) for i in n['arguments']))
    elif op=='negative':v=-at(n['argument'])
    elif op=='positive_quotient':v=sy.cancel(at(n['numerator'])/at(n['denominator']))
    elif op=='analytic_unary' and n['name'] in ('exp','log','sqrt'):
        v=getattr(sy,n['name'])(sy.expand(at(n['argument'])))
    else:raise AssertionError('Unbound exact linear algebra source '+op)
    cache[ref.node]=v;return v


def signed_partition_fixture():
    g=current.source.FunctionTransportGraph();z=sy.Symbol('Z',real=True);zr=g.symbol('Z')
    bind=lambda q:current.source.expression(g,sy.sympify(q),{'Z':zr})
    jet=lambda q:current.cutoff.c3.C3Function(*(bind(sy.diff(q,z,k)) for k in range(4)))
    inner=dict(zip(current.MOMENTS,(1+z,-1+2*z,3+z*z,z**3,2-z*z)))
    lower=dict(m13=-sy.Rational(1,10)+z,m15=sy.Rational(3,10)+2*z*z)
    preRp=dict(m13=2-z**3,m15=-3+z)
    pulse=dict(m13=-sy.Rational(7,10)+z*z,m15=sy.Rational(9,10)-z**3)
    pair=lambda q:current.source.C1Function(bind(q),bind(sy.diff(q,z)))
    moments=current.compose_full_moments(g,{k:jet(v) for k,v in inner.items()},
        *({k:pair(v) for k,v in fields.items()} for fields in (lower,preRp,pulse)))
    cache={}
    for key,value in moments.items():
        expected=inner[key]+sum((row.get(key,0) for row in (lower,preRp,pulse)),sy.Integer(0))
        assert sy.expand(symbolic(g,value.value,{zr.node:z},cache)-expected)==0,key
        assert sy.expand(symbolic(g,value.Z,{zr.node:z},cache)-sy.diff(expected,z))==0,(key,'Z')
    return dict(passed=True,mathematical_fixture_only=True,all_five_signed_C1_rows=True,
        all_three_exterior_pressure_intervals_counted_once=True,m15_not_negated_twice=True)


def universal_linear_identity(owner):
    """Actual selected graph expressions, arbitrary C1 defects and amplitude.

    The accepted matrix determinant identities identify the same-graph
    positive denominators with -det(A) and det(T). This substitutes exact
    identities, not numerical matrix values or defect enclosures.
    """
    g=owner.g;W=owner.linear.functions;b=owner.linear.background
    mu=sy.Symbol('original_mu',positive=True);Rc2=sy.Symbol('two_Rc',positive=True)
    Av=sy.Symbol('Astar',positive=True);Az=sy.Symbol('Astar_Z',real=True)
    A=sy.Matrix(2,2,lambda i,j:sy.Symbol(f'axial_{i}_{j}',real=True))
    T=sy.Matrix(3,3,lambda i,j:sy.Symbol(f'angular_{i}_{j}',real=True))
    bindings={W['mu'].node:mu,W['log2Rc'].node:sy.log(Rc2),
        b['Astar_ordinary_Z0_Z3'][0]:Av,b['Astar_ordinary_Z0_Z3'][1]:Az,
        W['positive_axial_divided_determinant'].node:-A.det(),
        W['positive_angular_determinant'].node:T.det()}
    for i,row in enumerate(W['axial_rows']):
        for j,ref in enumerate(row):bindings[ref.node]=A[i,j]
    for i,row in enumerate(W['angular_rows']):
        for j,ref in enumerate(row):bindings[ref.node]=T[i,j]
    for key,pair in owner.moments.items():
        for label,ref in (('value',pair.value),('Z',pair.Z)):
            if ref.node not in (g.zero.node,g.one.node):bindings[ref.node]=sy.Symbol(key+'_'+label,real=True)
    cache={};identities=[]
    for key,pair in owner.linear_response_residuals.items():
        for label,ref in (('value',pair.value),('Z',pair.Z)):
            value=symbolic(g,ref,bindings,cache)
            assert sy.cancel(sy.powsimp(value,force=True))==0,(key,label)
            identities.append(key+'_'+label)
    return dict(passed=True,actual_selected_graph_expressions=True,
        arbitrary_defect_and_amplitude_C1_identity_rows=identities,
        accepted_exact_determinant_identities_consumed=True,
        original_mu_Rc_Astar_factors_retained=True,Astar_Z_not_omitted=True,
        algebraic_response_only_not_repaired_own_terminal_fields=True)


def physical_response_fixture():
    with mp.workdps(75):
        g=current.source.FunctionTransportGraph();c=g.constant;z=sy.Symbol('Z',real=True);zr=g.symbol('Z')
        mu=g.symbol('finite_fixture_mu');logRc=g.unary('log',c('37/10'));log2Rc=g.unary('log',c('37/5'))
        bind=lambda q:current.source.expression(g,sy.sympify(q),{'Z':zr})
        A=sy.Rational(23,10)*(1+z*z)
        background=dict(parameters={'mu':mu.node},logRc=logRc.node,log2Rc=log2Rc.node,
            Astar_ordinary_Z0_Z3=[bind(sy.diff(A,z,k)).node for k in range(4)])
        W=current.linear.graph_weights(g,background)
        definitions=dict(zip(current.MOMENTS,(sy.Rational(37,100)+z*z,-sy.Rational(29,100)+2*z,
            sy.Rational(41,100)-z,-sy.Rational(23,100)+3*z*z,sy.Rational(19,100)+z**3)))
        moments={key:current.source.C1Function(bind(q),bind(sy.diff(q,z))) for key,q in definitions.items()}
        coefficients=current.linear.coefficient_action(g,W,background,moments)
        matrix=current.physical_response_matrix(g,W,background)
        residuals=current.response_residuals(g,moments,coefficients,matrix)
        patch=current.repair_velocity_functions(g,W,coefficients,background)
        zz=mp.mpf('.371');muv=mp.mpf('.08');Rc=mp.mpf('3.7');Av=mp.mpf('2.3')*(1+zz*zz)
        interp=linear_check.Interpreter(g,{mu.node:muv,zr.node:zz})
        beta=lambda t:mp.exp(-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)
        ell=mp.log(2)/400;centers=[j*mp.log(2)/40 for j in current.linear.CENTERS]
        integ=lambda fn:mp.quad(fn,[-1,-.75,0,.75,1])
        normal=integ(lambda t:beta(t)*mp.exp(ell*t/2))
        def physical_entry(j,k):
            ci=centers[j]
            def value(t):
                R=Rc*mp.exp(ci+ell*t)
                bump=mp.sqrt(2/Rc)*beta(t)/(ell*mp.exp(ci/2)*normal)
                return bump*(2*R)**k*R*ell
            return integ(value)
        expected=dict(m11=[physical_entry(j,0) for j in current.linear.AXIAL]+[mp.mpf(0)]*3,
            m12=[mp.mpf(0)]*2+[physical_entry(j,mp.mpf('.5')) for j in current.linear.ANGULAR],
            m13=[mp.mpf(0)]*2+[2*Av*physical_entry(j,-mp.mpf('1.5')-muv) for j in current.linear.ANGULAR],
            m14=[Av*physical_entry(j,-muv) for j in current.linear.AXIAL]+[mp.mpf(0)]*3,
            m15=[mp.mpf(0)]*2+[-Av*physical_entry(j,-mp.mpf('.5')-muv) for j in current.linear.ANGULAR])
        worst=mp.mpf(0)
        for key,row in matrix.items():
            for j,pair in enumerate(row):
                error=abs(interp.at(pair.value.node)-expected[key][j]);worst=max(worst,error)
                assert error<mp.mpf('1e-60'),(key,j,error)
        for key,pair in residuals.items():
            for ref in (pair.value,pair.Z):
                error=abs(interp.at(ref.node));worst=max(worst,error)
                assert error<mp.mpf('1e-57'),(key,error)
        for j in range(5):
            point=linear_check.Interpreter(g,{mu.node:muv,zr.node:zz,W['bump_log_radius_variable'].node:centers[j]})
            Ut=point.at(patch['Utheta1'].value.node);F=point.at(patch['F1'].value.node)
            assert abs(Ut-mp.sqrt(2*Rc*mp.exp(centers[j]))*F)<mp.mpf('1e-60')
        return dict(passed=True,mathematical_fixture_only=True,independent_physical_dR_bump_integrals=True,
            nonunit_Rc_Astar_and_finite_mu=True,m13_factor_two_and_m15_minus_preserved=True,
            five_value_and_Z_response_cancellations=True,physical_F_equals_Utheta_over_sqrt_2R=True,
            max_fixture_difference=str(worst))


def source_fixture(owner):
    assert set(owner.moments)==set(current.MOMENTS) and set(owner.coefficients)==set(current.COEFFICIENTS)
    g=owner.g
    for key in ('m11','m12','m14'):assert owner.moments[key]==owner.segments['inner'][key]
    expected=current.compose_full_moments(g,owner.cutoff.endpoint,owner.lower.functions,owner.preRp.functions,owner.pulse.functions)
    assert owner.moments==expected
    assert list(owner.coefficients.values())==owner.linear.coefficient_functions(owner.moments)
    assert all((p.value==g.zero and p.Z==g.zero) for p in owner.segments['beyond_Rv'].values())
    assert owner.reserved_background==current.reserved_background_guard(owner.linear)
    assert len(owner.reserved_background['actual_local_Uz_and_angular_increment_zero_rows'])==24
    receipt=current.parent.linear.read(current.linear.RECEIPT)
    assert receipt['determinant_identities']['axial_exact_difference_identity']
    assert receipt['determinant_identities']['angular_exact_Vandermonde_identity']
    assert receipt['actual_mu_positive_and_retained'] and receipt['original_amplitude_positive_on_all_four_Z_cells']
    assert not owner.definition['coefficient_norms_and_higher_derivative_bounds_transferred']
    assert not owner.definition['repaired_own_M_V_P_installed']
    return dict(passed=True,exact_retained_same_graph_source_owners=True,
        complete_signed_pressure_partition_and_zero_support=True,
        actual_complete_function_handles_selected_not_enclosure_centers=True,
        original_positive_invertibility_and_amplitude_guards_consumed=True,
        actual_leading_local_Uz_zero_and_pure_swirl_support_guard=True,
        nonlocal_histories_and_fixed_N1p_retained=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1FullDefectLinearFunctions(require_checked=False)
    raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['full_defect_linear_bindings']==owner.snapshot
    checks=dict(source=source_fixture(owner));print('TRUE_FULL_N1_DEFECT_SOURCES_CHECKED',flush=True)
    checks['partition']=signed_partition_fixture()
    checks['universal_linear']=universal_linear_identity(owner);print('FIVE_UNIVERSAL_LINEAR_C1_IDENTITIES_CHECKED',flush=True)
    checks['physical_response']=physical_response_fixture();print('PHYSICAL_BUMP_RESPONSE_AND_PATCH_UNITS_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        algebraic_linear_five_moment_response_identities_checked=True,
        coefficient_norms_and_higher_derivative_bounds_transferred=False,
        original_numeric_n1_solution_oracle_executed=False,numerical_original_n1_point_accuracy_claimed=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_FULL_DEFECT_AND_LINEAR_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
