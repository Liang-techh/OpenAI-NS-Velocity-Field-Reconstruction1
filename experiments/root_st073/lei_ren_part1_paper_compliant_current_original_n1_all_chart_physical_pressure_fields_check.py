"""Independent acceptance checks for the all-chart physical pressure adapter."""
import json
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as current
import lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields as physical


class SymbolicGraph:
    """Independent smooth-expression graph for the production Z-row helpers."""
    def __init__(self):
        self.zero = self.ref(s.Integer(0))
        self.one = self.ref(s.Integer(1))

    def ref(self, value):
        return current.source.FunctionRef(self, value)

    def constant(self, value):
        return self.ref(value if isinstance(value, s.Basic) else s.Rational(str(value)))

    def symbol(self, name):
        return self.ref(s.Symbol(name))

    def ids(self, values):
        assert all(value.graph is self for value in values)

    def add(self, *values):
        return self.ref(s.Add(*(value.node for value in values)))

    def neg(self, value):
        return self.ref(-value.node)

    def sub(self, left, right):
        return self.ref(left.node - right.node)

    def mul(self, *values):
        return self.ref(s.Mul(*(value.node for value in values)))

    def quotient(self, numerator, denominator, proof):
        return self.ref(numerator.node / denominator.node)

    def unary(self, name, value):
        return self.ref({'exp': s.exp}[name](value.node))

    def node(self, operation, **attrs):
        raise AssertionError('Fixture did not expect graph operation ' + operation)


def _jet(graph, expression, Z):
    return current.c3.C3Function(*(graph.ref(s.diff(expression, Z, k))
                                   for k in range(4)))


def _assert_rows(value, expression, Z, through=3):
    keys = ('value', 'Z', 'ZZ', 'ZZZ')
    for k, key in enumerate(keys):
        if k > through or not hasattr(value, key):
            break
        assert s.simplify(getattr(value, key).node - s.diff(expression, Z, k)) == 0


def independent_smooth_physical_fixture():
    """Exercise production q/products with an independently chosen physical mean."""
    graph = SymbolicGraph()
    Z, y = s.symbols('Z y')
    delta, scale = s.Rational(1, 5), s.Rational(3, 2)
    assert scale != 1
    m0 = (2 + y) + (1 - 2 * y) * Z + (y + 1) * Z**2 + y * Z**3
    delta_m = (1 + y) * Z - y * Z**2 + (2 - y) * Z**3
    V, delta_V = m0 + s.diff(m0, y), delta_m + s.diff(delta_m, y)
    U_expr, m_expr = scale * (V + delta_V), scale * (m0 + delta_m)
    U, m = _jet(graph, U_expr, Z), _jet(graph, m_expr, Z)
    q_expr = (2 * Z * U_expr - (1 - delta) * Z * m_expr
              - (1 - Z**2) * s.diff(m_expr, Z)) / (1 - delta * Z**2)
    q = current.partial.recover_Q_increment(graph, U, m, graph.ref(delta))
    _assert_rows(q, q_expr, Z, through=2)
    algebra = current.c3.C3Algebra(graph)
    mean_y = algebra.add(U, algebra.neg(m))
    assert s.simplify(U_expr - m_expr - s.diff(m_expr, y)) == 0
    _assert_rows(mean_y, s.diff(m_expr, y), Z)
    expected = {
        'I0_q_Z': s.diff(q_expr, Z),
        'I0_U_q_Z': U_expr * s.diff(q_expr, Z),
        'I0_U_squared': U_expr**2,
        'I0_U_U_Z': U_expr * s.diff(U_expr, Z),
        'I1_q': q_expr,
        'I1_q_Z': s.diff(q_expr, Z),
        'I1_q_squared': q_expr**2,
        'I1_U_q_Z': U_expr * s.diff(q_expr, Z),
        'I1_U_squared': U_expr**2,
        'I1_U_U_Z': U_expr * s.diff(U_expr, Z),
    }
    products = physical.nonradial_integrands(graph, U, q)
    for name, expression in expected.items():
        _assert_rows(products[name], expression, Z, through=1)
    # R = exp(y), so dR = R dy and R dR = R^2 dy.  Both measures are
    # checked on value/Z1 coefficient functions and against their y forms.
    R = s.exp(y)
    count = 7 + 2 * len(expected)
    for coefficient in (1 + y, 2 + 3 * y):
        dR = s.integrate(R * coefficient, (y, 0, 1))
        R_dR = s.integrate(R**2 * coefficient, (y, 0, 1))
        assert s.simplify(dR - s.integrate(s.diff(R, y) * coefficient, (y, 0, 1))) == 0
        assert s.simplify(R_dR - s.integrate(R * s.diff(R, y) * coefficient, (y, 0, 1))) == 0
        assert s.simplify(dR - s.integrate(coefficient, (y, 0, 1))) != 0
        assert s.simplify(R_dR - s.integrate(R * coefficient, (y, 0, 1))) != 0
        count += 4
    return count


def _reject(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('All-chart physical scope guard admitted an invalid query')


def _actual_patch_jacobian_binding_regression(owner):
    """Walk the production signed graph for nonconstant J=1/native_actual_patch."""
    chart = 'actual_patch'
    graph = owner.g
    base = owner.modified.target.phase.windows[chart]
    raw_jacobian = base['Jacobian']
    raw_node = graph.nodes[raw_jacobian]
    assert raw_node['operation'] == 'positive_quotient'
    denominator = raw_node['denominator']
    assert graph.nodes[denominator] == {
        'operation': 'bound_variable', 'name': 'native_actual_patch'}
    query = owner.functions[chart]['query']
    # The first exact source-cell defining integral is production output.
    first_cell = base['actual_current_source_cells'][0]
    key = current._interval_key(first_cell['exact_left'], first_cell['exact_right'])
    integral = owner.functions[chart]['defining_integrals_by_cell']['I0_q_Z'][key]
    definite = graph.nodes[integral.value.node]
    assert definite['operation'] == 'definite_integral'
    outer = graph.nodes[definite['integrand']]
    assert outer['operation'] == 'substitute_free_query_at_native_coordinate'
    assert outer['free_query'] == query.node
    assert outer['native_variable'] == 'pressure_moment_native_actual_patch'
    assert outer['bound_partial_mean_native_variable'] == 'native_actual_patch'

    # Search only the production integrand's dependency cone for the unique
    # native-J substitution.  This rejects a raw native variable left beside
    # the outer wrapper while preserving the inner Duhamel variable metadata.
    seen = set()
    stack = [outer['body']]
    bindings = []
    while stack:
        node_id = stack.pop()
        if node_id in seen:
            continue
        seen.add(node_id)
        node = graph.nodes[node_id]
        if (node.get('operation') == 'substitute_original_native_coordinate' and
                node.get('body') == raw_jacobian):
            bindings.append(node)
        for key in ('argument', 'body', 'integrand', 'query_coordinate'):
            value = node.get(key)
            if isinstance(value, int) and 0 <= value < len(graph.nodes):
                stack.append(value)
        for key in ('arguments',):
            values = node.get(key, ())
            if isinstance(values, list):
                stack.extend(value for value in values
                             if isinstance(value, int) and 0 <= value < len(graph.nodes))
    assert len(bindings) == 1
    binding = bindings[0]
    assert binding['native_variable'] == 'native_actual_patch'
    assert binding['query_coordinate'] == query.node
    return dict(raw_jacobian_node=raw_jacobian,
                raw_native_denominator=denominator,
                free_query_binding=binding,
                outer_pressure_native_variable=outer['native_variable'],
                inner_duhamel_native_variable=outer['bound_partial_mean_native_variable'])


@current.source_precision
def run(owner, packets):
    began = time.monotonic()
    fixture = independent_smooth_physical_fixture()
    owner.assert_graph()
    jacobian_check = _actual_patch_jacobian_binding_regression(owner)
    assert len(owner.charts) == 17
    assert owner.partitions['first_micro'][0] == ('sc_half',)
    assert owner.partitions['actual_patch'][-1] == 'Rh'
    expected_names = ({'I0_' + name for name in current.I0_NAMES} |
                      {'I1_' + name for name in current.I1_NAMES})
    assert set(packets) == set(owner.charts)
    for chart, packet in packets.items():
        view = owner.report(packet)
        partition = owner.partitions[chart]
        assert view['source_family'] == owner.family
        assert view['exact_native_partition'] == list(partition)
        assert len(view['exact_native_query_union']) == len(partition) - 1
        assert view['exact_original_source_cell_union']
        assert view['baseline_rows_are_source_rows_not_range_values']
        assert view['ordinary_baseline_factorial_conversion_once']
        assert view['physical_common_scale_Pstar_applied_once']
        assert view['physical_radius_original_geometry']
        assert view['retained_parent_partial_packet']
        assert view['no_finite_endpoint_q_R_or_support_claim']
        assert all(not view[key] for key in current.OPEN)
        assert len(view['whole_cell_physical_U_Z3_magnitude']) == 4
        assert len(view['whole_cell_physical_m_Z3_magnitude']) == 4
        assert len(view['whole_cell_physical_q_Z2_magnitude']) == 3
        assert set(view['whole_cell_value_Z1_integrand_magnitudes']) == expected_names
        assert set(view['whole_cell_integrated_value_Z1_magnitudes']) == expected_names
        assert all(len(rows) == 2 for rows in
                   view['whole_cell_value_Z1_integrand_magnitudes'].values())
        assert all(len(rows) == 2 for rows in
                   view['whole_cell_integrated_value_Z1_magnitudes'].values())
        assert set(view['signed_defining_integrals']) == expected_names

    # Packet fingerprints and exact-union endpoint guards remain live after a
    # child appends one graph symbol.
    token = packets['first_micro']
    entry = owner._packets[id(token)]
    key = 'physical_common_scale_Pstar_applied_once'
    old = entry[1]['view'][key]
    entry[1]['view'][key] = False
    _reject(lambda: owner.report(token))
    entry[1]['view'][key] = old
    _reject(lambda: owner.range_cell(('0', '.5'), 'first_micro', ('sc_half',), (3, 8)))
    owner.g.symbol('all_chart_physical_pressure_child')
    owner.report(token)

    # The receipt is written before this public checked-owner query.  This
    # verifies that a newly loaded checked owner can serve one nonmicro chart.
    hashes = dict(owner.hashes)
    for name in (current.NAME, Path(__file__).name, Path(physical.__file__).name):
        hashes[name] = current.sha(name)
    result = dict(
        all_passed=True, source_family=owner.family, input_hashes=hashes,
        checks=dict(
            independent_smooth_physical_scale_not_one_Q_Z2_mean_y_and_products=fixture,
            all_ten_integrand_value_and_Z1_rows_verified=True,
            physical_dR_and_R_dR_measures_verified=True,
            all_17_original_chart_full_unions_report=True,
            symbolic_sc_half_and_Rh_partitions_preserved=True,
            source_rows_factorials_scale_and_original_geometry_guards=True,
            nonconstant_actual_patch_jacobian_bound_to_free_query=jacobian_check,
            exact_union_and_packet_fingerprint_guards=True,
            support_endpoint_q_R_repair_n2_and_global_claims_not_installed=True),
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(current.c3.encoded(result), separators=(',', ':')) + '\n',
        encoding='utf-8', newline='\n')
    checked = current.CurrentOriginalN1AllChartPhysicalPressureFields(
        owner.parent, require_checked=True)
    chart = 'actual_patch'
    public = checked.range_chart(('0', '.5'), chart)
    checked_view = checked.report(public)
    assert checked_view['chart'] == chart
    assert checked_view[current.GATE]
    print('PASS_ALL_CHART_PHYSICAL_PRESSURE_FIELDS', fixture, len(packets), flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit('Use run(checked_all_chart_owner, packets)')
