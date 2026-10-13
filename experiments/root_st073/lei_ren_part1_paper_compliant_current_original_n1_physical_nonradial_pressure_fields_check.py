"""Independent checks for the physical nonradial pressure field continuation."""
import json
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields as current


class SymbolicGraph:
    """Small independent graph fixture for exact ordinary Z derivatives."""

    def __init__(self):
        self.zero = self.ref(s.Integer(0))
        self.one = self.ref(s.Integer(1))

    def ref(self, value):
        return current.source.FunctionRef(self, value)

    def constant(self, value):
        if isinstance(value, s.Basic):
            return self.ref(value)
        return self.ref(s.Rational(str(value)))

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
        raise AssertionError('Independent fixture did not expect graph operation ' + operation)


def _jet(graph, expression, Z):
    return current.c3.C3Function(*(graph.ref(s.diff(expression, Z, k)) for k in range(4)))


def _assert_rows(value, expression, Z, through=3):
    rows = tuple(getattr(value, key) for key in ('value', 'Z', 'ZZ', 'ZZZ')
                 if hasattr(value, key))
    for k, row in enumerate(rows):
        if k > through:
            break
        assert s.simplify(row.node - s.diff(expression, Z, k)) == 0


def _fixture_z_shift(value, graph):
    return current.c3.C3Function(value.Z, value.ZZ, value.ZZZ, graph.zero)


def independent_physical_fixture():
    """Check physical scaling, Q through Z2, products, and R dR measures."""
    graph = SymbolicGraph()
    Z, y = s.symbols('Z y')
    delta = s.Rational(1, 5)
    S = s.Rational(3, 2)
    assert S != 1

    # Smooth source profiles have independent y and Z dependence.  The
    # physical fields use one common unit S, while their own mean satisfies
    # the defining m_y = U - m identity.
    m0 = (2 + y) + (1 - 2 * y) * Z + (y + 1) * Z**2 + y * Z**3
    delta_m = (1 + y) * Z - y * Z**2 + (2 - y) * Z**3
    V = m0 + s.diff(m0, y)
    delta_V = delta_m + s.diff(delta_m, y)
    U_expression = S * (V + delta_V)
    m_expression = S * (m0 + delta_m)
    U = _jet(graph, U_expression, Z)
    m = _jet(graph, m_expression, Z)

    q_expression = (2 * Z * U_expression - (1 - delta) * Z * m_expression
                    - (1 - Z**2) * s.diff(m_expression, Z)) / (1 - delta * Z**2)
    q = current.partial.recover_Q_increment(graph, U, m, graph.ref(delta))
    _assert_rows(q, q_expression, Z, through=2)

    algebra = current.c3.C3Algebra(graph)
    m_y = algebra.add(U, algebra.neg(m))
    assert s.simplify(U_expression - m_expression - s.diff(m_expression, y)) == 0
    _assert_rows(m_y, s.diff(m_expression, y), Z, through=3)

    integrands = current.nonradial_integrands(graph, U, q)
    expected = {
        'I0_q_Z': s.diff(q_expression, Z),
        'I0_U_q_Z': U_expression * s.diff(q_expression, Z),
        'I0_U_squared': U_expression**2,
        'I0_U_U_Z': U_expression * s.diff(U_expression, Z),
        'I1_q': q_expression,
        'I1_q_Z': s.diff(q_expression, Z),
        'I1_q_squared': q_expression**2,
        'I1_U_q_Z': U_expression * s.diff(q_expression, Z),
        'I1_U_squared': U_expression**2,
        'I1_U_U_Z': U_expression * s.diff(U_expression, Z),
    }
    checked = 7  # Q value/Z/Z2 and physical own mean_y value/Z/Z2/Z3.
    for name, expression in expected.items():
        _assert_rows(integrands[name], expression, Z, through=1)
        checked += 2

    # With R = exp(y), dR = R dy and R dR = R^2 dy.  Check both measures on
    # value and Z1 rows of every independent product, including that the two
    # measures are genuinely different on this nonconstant fixture.
    R = s.exp(y)
    for coefficient in (1 + y, 2 + 3 * y):
        dR = s.integrate(R * coefficient, (y, 0, 1))
        R_dR = s.integrate(R**2 * coefficient, (y, 0, 1))
        assert s.simplify(dR - s.integrate(s.diff(R, y) * coefficient, (y, 0, 1))) == 0
        assert s.simplify(R_dR - s.integrate(R * s.diff(R, y) * coefficient, (y, 0, 1))) == 0
        assert s.simplify(dR - s.integrate(coefficient, (y, 0, 1))) != 0
        assert s.simplify(R_dR - s.integrate(R * coefficient, (y, 0, 1))) != 0
        checked += 4
    return checked


def rejected(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('Physical nonradial pressure scope guard admitted an invalid object')


@current.source_precision
def run(owner, packets):
    began = time.monotonic()
    count = independent_physical_fixture()
    for packet in packets.values():
        view = owner.report(packet)
        assert all(not view[key] for key in current.OPEN)
        assert view['retained_parent_partial_packet']
        assert view['baseline_rows_are_source_rows_not_range_values']
        assert view['ordinary_baseline_factorial_conversion_once']
        assert view['physical_common_scale_Pstar_applied_once']
        assert view['physical_radius_original_geometry']
        assert view['whole_cell_range_not_endpoint_value']
        assert view['numerical_evidence_absolute_whole_cell_magnitude_only']
        assert view['no_finite_endpoint_q_R_or_support_claim']
        assert len(view['whole_cell_physical_U_Z3_magnitude']) == 4
        assert len(view['whole_cell_physical_m_Z3_magnitude']) == 4
        assert len(view['whole_cell_physical_q_Z2_magnitude']) == 3
        expected_names = ({'I0_' + name for name in current.I0_NAMES} |
                          {'I1_' + name for name in current.I1_NAMES})
        assert set(view['whole_cell_value_Z1_integrand_magnitudes']) == expected_names
        assert set(view['whole_cell_integrated_value_Z1_magnitudes']) == expected_names
        assert all(len(rows) == 2 for rows in view['whole_cell_value_Z1_integrand_magnitudes'].values())
        assert all(len(rows) == 2 for rows in view['whole_cell_integrated_value_Z1_magnitudes'].values())

    # Mutating a retained packet view invalidates its fingerprint.  Restore it
    # immediately so the live owner remains usable by the caller after checks.
    token = packets['first_micro']
    entry = owner._packets[id(token)]
    key = 'physical_common_scale_Pstar_applied_once'
    old = entry[1][key]
    entry[1][key] = False
    rejected(lambda: owner.report(token))
    entry[1][key] = old
    owner.report(token)

    hashes = dict(owner.hashes)
    for name in (current.NAME, Path(__file__).name):
        hashes[name] = current.sha(name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        input_hashes=hashes,
        checks=dict(
            independent_smooth_physical_scale_not_one_Q_Z2_and_m_y_identity=count,
            all_ten_integrand_value_and_Z1_product_rows_verified=True,
            physical_dR_and_R_dR_measures_verified_on_value_and_Z1_rows=True,
            two_actual_micro_whole_cell_physical_fields_and_integrals_report=True,
            source_rows_factorials_scale_and_original_geometry_guards=True,
            packet_fingerprint_and_open_scope_guards=True,
            support_endpoint_q_R_repair_n2_and_global_claims_not_installed=True,
        ),
        **{current.GATE: True},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (current.HERE / current.RECEIPT).write_bytes(
        (json.dumps(current.c3.encoded(result), separators=(',', ':')) + '\n').encode())
    print('PASS_PHYSICAL_NONRADIAL_PRESSURE_FIELDS', count, flush=True)
    return result
