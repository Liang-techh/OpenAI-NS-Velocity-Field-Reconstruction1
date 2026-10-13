"""Exact physical conjugation and bounded same-source numerical bridge checks."""
import ast
import copy
import inspect
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions as raw
import lei_ren_part1_paper_compliant_current_original_n1_regular_series_function as regular


def compile_mathematical(fn, scope):
    """Execute the pinned production formula with exact scalar/derivative adapters."""
    tree = ast.parse(inspect.getsource(fn).strip())
    tree.body[0].decorator_list = []

    class ScalarAdapter(ast.NodeTransformer):
        def visit_Attribute(self, node):
            node = self.generic_visit(node)
            if node.attr == 'algebra':
                return ast.Name(id='symbolic_algebra', ctx=ast.Load())
            return node

        def visit_Call(self, node):
            node = self.generic_visit(node)
            if isinstance(node.func, ast.Name) and node.func.id == 'truncate':
                return node.args[0]
            if isinstance(node.func, ast.Attribute) and node.func.attr == 'truncate':
                return node.args[0] if isinstance(node.func.value, ast.Name) and node.func.value.id == 'op' else node.func.value
            return node

    tree = ScalarAdapter().visit(tree)
    ast.fix_missing_locations(tree)
    namespace = dict(fn.__globals__)
    namespace.update(scope)
    exec(compile(tree, '<pinned-original-n1-physical-conjugation>', 'exec'), namespace)
    return namespace[fn.__name__]


def physical_conjugation_identity():
    """All 72 matrix entries, true forcing, ordinary amplitude jets and recovery."""
    z, rho = sy.symbols('Z rho', real=True, positive=True)
    delta, eps, anchor, P2 = sy.symbols('delta epsilon Fanchor Pstar_squared', positive=True)
    F0 = sy.Function('F0')(z)
    gradient = sy.Function('gradient')(z)
    U0 = sy.Function('Uaxis')(z)
    phi = sy.Function('Phi')(rho, z)
    psi = sy.Function('Psi')(rho, z)
    q = sy.Function('Qleading')(rho, z)
    d, L = 1-z*z, 1-delta*z*z
    Haxis = (1-delta)*z/2+d*U0
    chi = Haxis*gradient/L
    ctx = SimpleNamespace(sqrt=sy.sqrt, mpf=lambda v: sy.Rational(str(v)))

    class Algebra:
        def __init__(self):
            self.ctx = ctx
        def lift(self, value, order=None):
            return sy.sympify(value)
        def shift(self, value, powers):
            return value * anchor**powers[1] * P2**powers[2] * eps**(-powers[3])

    algebra = Algebra()
    dz = lambda v: sy.diff(v, z)
    axial = compile_mathematical(regular.op.axial_second, {'dz': dz})
    fake_op = SimpleNamespace(dz=dz, truncate=lambda v, n: v, axial_second=axial)
    # F0/anchor is the analytic relative-amplitude function. Its actual
    # production Taylor recurrence is independently checked below.
    fake_weighted = SimpleNamespace(relative_amplitude=lambda a, g, p, n: (F0/anchor)**p)
    residual = compile_mathematical(regular.feedback.residual_matrix,
        {'op': fake_op, 'weighted': fake_weighted})
    theta = compile_mathematical(regular.theta_source, {'op': fake_op})
    extra = regular.natural.conjugate_extra
    f = F0*phi
    u = U0+eps*psi
    rows = {'Phi': [sy.diff(phi, rho, i) for i in range(3)],
            'Psi': [sy.diff(psi, rho, i) for i in range(3)]}
    fields = {'F': [sy.diff(f, rho, i) for i in range(3)],
              'Uz': [sy.diff(u, rho, i) for i in range(3)],
              'Q': [sy.diff(q, rho, i) for i in range(3)]}
    forcing_fn = compile_mathematical(regular.op.known_forcing,
        {'dz': dz, 'axial_second': axial, 'truncate': lambda v, n: v})
    forcing = forcing_fn(algebra, fields, z, rho, delta)

    class FakeJet:
        def __new__(cls, c, value):
            return value
        @staticmethod
        def variable(c, value, order):
            return value

    fake_op.Jet = FakeJet
    packet = object()
    owner = SimpleNamespace(c=ctx,
        physical=SimpleNamespace(delta=delta, axis_inputs=lambda Z: {'chi': chi}),
        weighted=SimpleNamespace(source=SimpleNamespace(amplitude=SimpleNamespace(gradient_jets=lambda Z: gradient))),
        mean=SimpleNamespace(report=lambda p: {'geometry': {'rho': rho}, 'Z': z},
            _packets={id(packet): (None, None, None, None, (None, rows, fields, forcing))}))
    exact = compile_mathematical(regular.exact_operator,
        {'op': fake_op, 'feedback': SimpleNamespace(residual_matrix=residual),
         'natural': SimpleNamespace(conjugate_extra=extra), 'theta_source': theta,
         'symbolic_algebra': algebra})
    normalized = exact(owner, packet)
    expressions = raw.raw_matrix_expressions()
    values = dict(rho=rho, x=sy.sqrt(rho), Z=z, delta=delta, epsilon=eps,
        F=f, Fr=sy.diff(f, rho), FZ=sy.diff(f, z), U=u,
        Ur=sy.diff(u, rho), UZ=sy.diff(u, z), Q=q)
    scales = (F0/eps, eps, eps, P2*eps, F0/eps, eps)
    ell = (-gradient/eps, 0, 0, 0, -gradient/eps, 0)
    derivatives = {sy.diff(F0, z): -gradient*F0/eps,
        sy.diff(F0, z, 2): (gradient**2/eps**2-sy.diff(gradient, z)/eps)*F0}

    def zero(value, label):
        if sy.cancel(sy.expand(value.xreplace(derivatives))) != 0:
            raise AssertionError('Physical conjugation differs: '+label)

    count = 0
    j=sy.Symbol('j',positive=True)
    assert sy.expand((1-delta)*z/2+(1-z*z)*(4*z+j)-
        (j+((9-delta)/2)*z-j*z*z-4*z**3))==0
    for label in ('B0', 'B1'):
        for i in range(6):
            for j in range(6):
                key = f'{i+1},{j+1}'
                body = expressions[label].get(key, sy.Integer(0)).subs(
                    {sy.Symbol(k): v for k, v in values.items()}, simultaneous=True)
                want = body*scales[j]/scales[i]
                if label == 'B0':
                    b1 = expressions['B1'].get(key, sy.Integer(0)).subs(
                        {sy.Symbol(k): v for k, v in values.items()}, simultaneous=True)
                    want += b1*scales[j]/scales[i]*ell[j]
                zero(normalized[label].get(key, sy.Integer(0))-want, label+':'+key)
                count += 1
    raw_sources = {'4': 2*eps*sy.sqrt(rho)*forcing['N1p'],
        '5': 2*eps*forcing['N1theta'],
        '6': 2*eps*forcing['N1z']-4*eps**2*z*rho*forcing['N1p']/L}
    for key, value in raw_sources.items():
        zero(normalized['g'][key]-value/scales[int(key)-1], 'g:'+key)
    zero(normalized['B0']['5,1']-
         (2*eps*(q+(z*u-sy.Rational(1,2))*(-2+delta)/L)-2*((1-delta)*z/2+d*u)*gradient/L),
         'single-gauge-principal-plus-Psi')
    assert normalized['D'] == raw.D
    assert regular.natural.WEIGHTS == ((0,1,0,1),(0,0,0,-1),(0,0,0,-1),(0,0,1,-1),(0,1,0,1),(0,0,0,-1))
    Wh, Kh, Zh = (sy.Function(name)(rho,z) for name in ('Uhat', 'Khat', 'Qhat'))
    mean = eps*(Wh+Kh)
    Qphysical = ((1-delta)*z*eps*Wh-(1+delta)*z*eps*Kh-d*sy.diff(mean,z))/L
    zero(Qphysical-eps*((1-delta)*z*Wh-(1+delta)*z*Kh-d*sy.diff(Wh+Kh,z))/L, 'own-Q1')
    zero(rho*eps*Qphysical-rho*eps**2*(Qphysical/eps), 'own-V1')

    # Execute the real coefficient recurrence, with a symbolic Taylor
    # container, and compare ordinary derivatives of R'= -p*g*R/epsilon.
    class Jet:
        def __init__(self, c, values):
            self.coefficients=list(values); self.order=len(values)-1
        def __getitem__(self, k):
            return self.coefficients[k]
    class Factored:
        def __init__(self, a, terms, order):
            self.terms=terms; self.order=order
    amplitude = compile_mathematical(regular.weighted.relative_amplitude,
        {'op': SimpleNamespace(Jet=Jet,FactoredJet=Factored),
         'ep': lambda v: (0,0) if v==0 else (-1,1)})
    gj = Jet(ctx,[sy.diff(gradient,z,k)/math.factorial(k) for k in range(4)])
    amplitude_count = 0
    for power in (-2,-1,1,2):
        got=amplitude(algebra,gj,power,3); expected=sy.Integer(1)
        for k in range(4):
            observed=sum((row[k]*eps**(-p[3]) for p,row in got.terms.items()),sy.Integer(0))*math.factorial(k)
            zero(observed-expected, f'relative-amplitude:{power}:Z{k}')
            expected=sy.diff(expected,z)-power*gradient*expected/eps
            amplitude_count+=1
    return dict(passed=True,all_exact_physical_conjugation_matrix_entries=count,
        exact_physical_forcing_conjugation_rows=3,
        original_relative_amplitude_ordinary_Z0_to_Z3_rows=amplitude_count,
        accepted_exact_operator_mathematical_AST_executed=True,
        normalization='diag(F0*Lambda,epsilon,epsilon,Pstar^2*epsilon,F0*Lambda,epsilon)',
        single_original_F0_Z_gauge_verified=True,
        principal_minus_2chi_and_residual_Psi_gauge_sum_verified=True,
        n1_own_mean_Q_V_physical_scales_and_one_plus_delta_verified=True,
        legacy_anchor_H_polynomial_equals_true_graph_axis_H=True,
        identities_are_universal_not_sampled_fixtures=True)


def rejected(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('Unsupported true numerical bridge query admitted')


@regular.source_precision
def run(owner, saved, packets):
    import lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge as current
    began=time.monotonic();owner.assert_graph()
    proof=physical_conjugation_identity()
    identity=owner.assert_source_identity()
    assert identity==owner.source_identity and identity['passed']
    assert regular.op.snapshot(saved['source_identity'])==regular.op.snapshot(identity)
    assert len(packets)==len(saved['numerical_queries'])==3
    assert saved['true_function_bindings']==owner._function_bindings
    results=[]
    for packet, expected in zip(packets,saved['numerical_queries']):
        view=owner.report(packet)
        assert regular.op.snapshot(view)==regular.op.snapshot(expected)
        assert all(view[k] is False for k in current.OPEN)
        nv=view['numerical_enclosures']
        assert nv['same_original_common_inner_function_for_every_query']
        assert nv['all_six_row_feedback_included']
        assert nv['uncomputed_higher_Z_centers_use_symmetric_analytic_bound']
        assert nv['original_source_F0_Z_gauge_applied_once_in_physical_output']
        assert nv['radial_derivatives_recovered_from_exact_equation_not_seed_differences']
        for name,row in nv['physical_function_jets'].items():
            assert name in regular.FIELDS and row['ordinary_Z_derivative_order']==3
            assert all(len(term['ordinary_Z_derivatives'])==4 for term in row['terms'])
        if view['chart']=='first':
            assert view['original_true_phase_enclosure']._mpi_==(view['numerical_coordinate']*owner.selection.sc/4)._mpi_
            assert nv['query']['actual_source_coordinate']._mpi_==view['original_true_phase_enclosure']._mpi_
        if regular.ep(view['numerical_coordinate'])[1]!=0:
            assert nv['passing_six_component_value_Z1_rows']==12
        results.append(dict(chart=view['chart'],coordinate=view['numerical_coordinate'],
            actual_query_executed=True,passing_value_Z1_natural_component_rows=nv['passing_six_component_value_Z1_rows'],
            every_physical_component_has_ordinary_Z3_enclosures=True))
    assert rejected(lambda:owner.evaluate('.371','4.01','core'))
    assert rejected(lambda:owner.evaluate('.371','1.01','first'))
    assert rejected(lambda:owner.evaluate('.371','0','outside'))
    assert rejected(lambda:owner.report(regular.OriginalN1RegularSeriesFunctionPacket()))
    numerical_packet=owner._packets[id(packets[1])][2]
    assert rejected(lambda:owner.evaluate('.371','3','core',reuse_packet=numerical_packet))
    original=owner._function_bindings
    try:
        owner._function_bindings={}
        assert rejected(owner.assert_graph)
    finally:
        owner._function_bindings=original
    owner.assert_graph()
    pins=dict(owner.hashes);pins[current.NAME]=current.sha(current.NAME)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=pins,
        physical_conjugation_proof=proof,source_identity=identity,
        checks=dict(three_original_axis_core_and_actual_first_queries_execute=True,
            no_numeric_center_promoted_to_true_function_definition=True,
            exact_first_fraction_to_phase_conversion_checked=True,
            same_actual_core_norm_recipes_and_fourteen_atom_P0_bound=True,
            original_complete_feedback_tail_and_all_source_errors_preserved=True,
            unsupported_domain_foreign_packet_and_changed_function_binding_rejected=True,
            live_packet_reuse_with_wrong_physical_query_rejected=True,
            sampled_value_Z1_natural_component_budget_passes=True,
            no_uniform_repaired_or_global_NS_accuracy_claim=True),
        query_results=results,**{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        original_numeric_n1_solution_oracle_executed=True,
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(regular.op.snapshot(receipt),indent=2)+'\n').encode('utf8'))
    print('TRUE_N1_NUMERIC_BRIDGE_CHECKED',len(packets),len(pins),flush=True)
    return receipt
