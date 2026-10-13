"""Direct physical-pressure identity and exact original DAG relocation checks."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as current


def symbolic(g, ref, bindings=None, cache=None):
    bindings, cache = bindings or {}, cache if cache is not None else {}
    if ref.node in bindings:
        return bindings[ref.node]
    if ref.node in cache:
        return cache[ref.node]
    row = g.nodes[ref.node]
    at = lambda i: symbolic(g, current.source.FunctionRef(g, i), bindings, cache)
    op = row['operation']
    if op == 'exact_rational':
        value = s.Rational(row['numerator'], row['denominator'])
    elif op == 'bound_variable':
        value = s.Symbol(row['name'])
    elif op == 'sum':
        value = s.Add(*map(at, row['arguments']))
    elif op == 'product':
        value = s.Mul(*map(at, row['arguments']))
    elif op == 'negative':
        value = -at(row['argument'])
    elif op == 'positive_quotient':
        value = at(row['numerator']) / at(row['denominator'])
    elif op == 'analytic_unary':
        value = getattr(s, row['name'])(at(row['argument']))
    elif op == 'function_substitution':
        value = at(row['expression']).subs(at(row['variable']), at(row['value']))
    elif op == 'definite_integral':
        value = s.integrate(at(row['integrand']),
                            (s.Symbol(row['variable']), at(row['lower']), at(row['upper'])))
    else:
        raise AssertionError('Fixture cannot turn a source/cap node into a value: ' + op)
    cache[ref.node] = value
    return value


def independent_pressure_fixture():
    """Direct N1p integration checks signs, finite endpoints and true C1 rows."""
    g = current.source.FunctionTransportGraph()
    Z, R = s.symbols('Z pressure_fixture_R')
    variables = {str(v): g.symbol(str(v)) for v in (Z, R)}
    bind = lambda p: current.source.expression(g, s.sympify(p), variables)
    c1 = lambda p: current.source.C1Function(bind(p), bind(s.diff(p, Z)))
    jet = lambda p: current.active.c3.C3Function(*(bind(s.diff(p, Z, k)) for k in range(4)))
    delta, physical_unit = s.Rational(1, 5), s.Rational(3, 2)
    lo, hi = s.Rational(37, 10), s.Rational(74, 10)
    # U=m+R*m_R supplies the true own-mean equation. The unit and physical
    # radius are independently nonunit, making omitted scale factors visible.
    m = physical_unit * (2 + Z + Z**2 + (1 - 2*Z + Z**3)*R)
    U = m + R*s.diff(m, R)
    d, L = 1-Z**2, 1-delta*Z**2
    q = s.cancel((2*Z*U-(1-delta)*Z*m-d*s.diff(m, Z))/L)
    qy = R*s.diff(q, R)
    source_qplus = current.q_plus_q_y(g, bind(delta), jet(U), c1(R*s.diff(U, R)))
    cache = {}
    assert s.cancel(symbolic(g, source_qplus.value, cache=cache) - q - qy) == 0
    assert s.cancel(symbolic(g, source_qplus.Z, cache=cache) - s.diff(q+qy, Z)) == 0
    products = dict(I0_q_Z=s.diff(q, Z), I0_U_q_Z=U*s.diff(q, Z),
        I0_U_squared=U**2, I0_U_U_Z=U*s.diff(U, Z), I1_q=q,
        I1_q_Z=s.diff(q, Z), I1_q_squared=q**2,
        I1_U_q_Z=U*s.diff(q, Z), I1_U_squared=U**2, I1_U_U_Z=U*s.diff(U, Z))
    integrate = lambda p, power: s.integrate(s.Poly(s.cancel(p*R**power), R).as_expr(), (R, lo, hi))
    integrals = {name: c1(integrate(body, 0 if name.startswith('I0_') else 1))
                 for name, body in products.items()}
    a = current.C1Algebra(g)
    boundary = []
    for edge in (lo, hi):
        at = lambda p: a.at(c1(p), variables[str(R)], bind(edge))
        boundary.append(current.endpoint_functions(g, at(U), at(q), at(qy), bind(edge)))
    got = current.pressure_functions(g, bind(delta), integrals, *boundary)
    # Independently integrate the original physical forcing; do not reuse
    # the reduced-S formula or producer's pressure moment implementation.
    T = (q + Z*(1-delta)*s.diff(q, Z)/2 + R*s.diff(q, R))/L
    convective = q*(q/2+R*s.diff(q, R))
    axial = U*(-2*Z*q+d*s.diff(q, Z)-2*Z*R*s.diff(q, R))/L
    N1p = -(T+convective+axial)/2 + 2*s.diff(q, R) + R*s.diff(q, R, 2)
    direct = dict(m13=integrate(N1p, 0), m15=-integrate(N1p, 1))
    for name, pair in got.items():
        assert s.cancel(symbolic(g, pair.value, cache=cache)-direct[name]) == 0
        assert s.cancel(symbolic(g, pair.Z, cache=cache)-s.diff(direct[name], Z)) == 0
    # Also exercise the production compact builder: an independent nonunit
    # Rc fixture detects omitted or duplicated Rc and Rc^2*x measures.
    x = s.Symbol('pressure_fixture_x')
    variables[str(x)] = g.symbol(str(x))
    common_m = (m/physical_unit).subs(R, lo*x)
    common_u = (U/physical_unit).subs(R, lo*x)
    encode_jet = lambda p: [v.node for v in current.active.c3.rows(jet(p))]
    original_band = dict(original_band_variable=variables[str(x)].node,
        profiles_at_band_x=dict(radial_y_rows=[dict(V=encode_jet(common_u)),
            dict(V=encode_jet(x*s.diff(common_u, x)))]),
        complete_histories=dict(m=encode_jet(common_m)),
        complete_history_y=dict(m=encode_jet(x*s.diff(common_m, x))))
    fake = SimpleNamespace(graph=g,
        background=dict(logRc=g.unary('log', bind(lo)).node),
        reports={current.linear.quiet.band.NAME: dict(actual_C3_repair_band_functions=original_band)})
    compact, _ = current.compact_functions(fake, bind(physical_unit), bind(delta))
    for name, pair in compact.items():
        assert s.cancel(symbolic(g, pair.value, cache=cache)-direct[name]) == 0
        assert s.cancel(symbolic(g, pair.Z, cache=cache)-s.diff(direct[name], Z)) == 0
    assert all(s.cancel(v) != 0 for v in direct.values())
    return dict(passed=True, independent_mathematical_fixture=True,
        original_source_evaluation=False, delta=str(delta), physical_unit=str(physical_unit),
        physical_R_endpoints=[str(lo), str(hi)], both_finite_viscous_endpoints_retained=True,
        correct_m15_minus_applied_once=True, both_true_ordinary_Z_rows_checked_symbolically=True)


def relocation_fixture():
    original = current.source.FunctionTransportGraph()
    x = original.symbol('x')
    seven = original.constant(7)
    row = original.node('original_lazy_flat_branch', active_body=seven.node,
        Delta=x.node, eta=original.one.node, flat_value=7,
        flat_predicate=original.node('exact_real_comparison', operator='>=',
            left=x.node, right=seven.node).node,
        original_collar_coordinate=None, original_collar_upper=None)
    g = current.source.FunctionTransportGraph()
    for j in range(11):
        g.symbol('unrelated_' + str(j))
    importer = current.GraphImport(g, original)
    imported = importer(row)
    node = g.nodes[imported.node]
    assert node['flat_value'] == 7
    assert node['original_collar_coordinate'] is None
    assert g.nodes[node['active_body']] == dict(operation='exact_rational', numerator=7, denominator=1)
    assert g.nodes[node['Delta']] == dict(operation='bound_variable', name='x')
    assert node['active_body'] != seven.node
    try:
        importer(g.one)
    except TypeError:
        pass
    else:
        raise AssertionError('Foreign graph accepted')
    return dict(passed=True, literal_flat_value_and_optional_null_preserved=True,
                exact_rational_fields_not_relocated=True, foreign_graph_rejected=True)


def actual_graph_relocation(owner):
    mapping = owner.bindings['active']['import_mapping']
    old, new = owner.fields.g, owner.g
    moved = 0
    for i, ref in mapping.items():
        expected = copy.deepcopy(old.nodes[i])
        for key in current.REFERENCE_KEYS[expected['operation']]:
            value = expected.get(key)
            if value is not None:
                expected[key] = ([mapping[v].node for v in value] if isinstance(value, list)
                                 else mapping[value].node)
        assert new.nodes[ref.node] == expected
        moved += i != ref.node
    assert moved > 0 and len(owner.fields.charts) == 17
    # The actual_patch logarithmic Jacobian keeps both variable substitutions.
    candidates = [row for row in new.nodes if row['operation'] == 'substitute_free_query_at_native_coordinate'
                  and row['native_variable'] == 'pressure_moment_native_actual_patch']
    assert candidates and all(row['inner_Duhamel_native_variable_preserved'] for row in candidates)
    lower = [row for row in new.nodes if row['operation'] == 'current_original_physical_source_quantity'
             and row.get('quantity') == 'V_y']
    assert {row['Z_order'] for row in lower} == {0, 1}
    assert all(row['raw_quantity_path'].startswith('actual_background_source.V_y[') for row in lower)
    return dict(passed=True, exact_source_nodes_checked=len(mapping), relocated_nodes=moved,
        all_17_source_charts_retained=True, actual_patch_inner_Duhamel_variable_preserved=True,
        true_lower_Uy_and_UyZ_provider_bound=True)


def run(owner):
    began = time.monotonic()
    owner.assert_graph()
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert raw['source_family'] == owner.family
    assert raw['exact_graph_nodes'] == owner.g.nodes
    assert raw['signed_pressure_C1_functions'] == current.encode(owner.functions)
    assert set(owner.pieces) == {'active', 'compact', 'quiet'}
    assert raw['inner_lower_collar_and_outer_pulse_not_included']
    assert raw['complete_five_defects_not_installed']
    assert all(not raw[key] for key in current.OPEN)
    physical = independent_pressure_fixture()
    relocation = relocation_fixture()
    live = actual_graph_relocation(owner)
    for name, digest in owner.hashes.items():
        if current.sha(name) != digest:
            raise AssertionError('Changed original source ' + name)
    checks = dict(physical_pressure=physical, synthetic_relocation=relocation,
                  actual_graph_relocation=live)
    receipt = dict(all_passed=True, source_family=owner.family, checks=checks,
        radial_scope=raw['radial_scope'], signed_function_graph_only=True,
        actual_full_defect_coefficients_not_selected=True,
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
                     Path(__file__).name: current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_PRE_RP_SIGNED_C1_FUNCTIONS_CHECKED', live['exact_source_nodes_checked'], flush=True)
    return receipt


if __name__ == '__main__':
    independent_pressure_fixture(); relocation_fixture()
    print('INDEPENDENT_PRE_RP_PRESSURE_FUNCTION_FIXTURES_PASSED', flush=True)
