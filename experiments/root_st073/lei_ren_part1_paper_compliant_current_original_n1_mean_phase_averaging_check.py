"""Independent checks for signed finite-N own-mean phase averaging."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_mean_phase_averaging as current


class SymbolicGraph:
    """Small expression graph used only to check the production IBP builder."""
    def __init__(self):
        self.zero = self.ref(sy.Integer(0))
        self.one = self.ref(sy.Integer(1))

    def ref(self, value):
        return current.source.FunctionRef(self, value)

    def constant(self, value):
        if isinstance(value, sy.Basic):
            return self.ref(value)
        return self.ref(sy.Rational(str(value)))

    def symbol(self, name):
        return self.ref(sy.Symbol(name))

    def ids(self, values):
        assert all(value.graph is self for value in values)

    def add(self, *values):
        self.ids(values)
        return self.ref(sy.Add(*(value.node for value in values)))

    def neg(self, value):
        self.ids([value])
        return self.ref(-value.node)

    def sub(self, left, right):
        self.ids([left, right])
        return self.ref(left.node - right.node)

    def mul(self, *values):
        self.ids(values)
        return self.ref(sy.Mul(*(value.node for value in values)))

    def quotient(self, numerator, denominator, proof):
        self.ids([numerator, denominator])
        return self.ref(numerator.node / denominator.node)

    def unary(self, name, value):
        self.ids([value])
        if name != 'exp':
            raise AssertionError('Unexpected fixture unary ' + name)
        return self.ref(sy.powsimp(sy.exp(value.node), force=True))

    def node(self, operation, **attrs):
        if operation == 'substitute_original_native_coordinate':
            body = attrs['body']
            native = sy.Symbol(attrs['native_variable'])
            return self.ref(body.subs(native, attrs['query_coordinate']))
        if operation == 'definite_integral':
            variable = sy.Symbol(attrs['variable'])
            value = sy.integrate(attrs['integrand'],
                                 (variable, attrs['lower'], attrs['upper']))
            return self.ref(sy.simplify(value))
        raise AssertionError('Unexpected production operation ' + operation)


def independent_ibp_fixture():
    """Check all four ordinary Z rows against direct moving-query integration.

    The independently chosen G has a nonzero inlet and third Z derivative;
    H is a nonconstant fixed-phase slow-y envelope.  The Jacobian is also
    nonconstant, so the fixture detects omission or duplication of the native
    physical measure.
    """
    g = SymbolicGraph()
    n, q, z = sy.symbols('native_fixture q Z')
    N = sy.Integer(7)
    offset_expr = sy.log(1 + n)
    jacobian_expr = 1 / (1 + n)
    G_expr = (3 + n + n**2) * (1 + z + z**2 + z**3)
    H_expr = (1 + 2*n) * (1 + 2*z + z**2 + 3*z**3)
    query = g.ref(q)
    offset = g.ref(offset_expr)
    left_offset = g.zero
    jacobian = g.ref(jacobian_expr)
    direct_kernel = sy.exp(-(sy.log(1 + q) - sy.log(1 + n)))
    count = 0
    for order in range(4):
        G_row = g.ref(sy.diff(G_expr, z, order))
        H_row = g.ref(sy.diff(H_expr, z, order))
        got = current._ibp_row_from_handles(
            g, 'fixture', query, offset, left_offset, G_row, H_row,
            jacobian, g.zero,
            g.quotient(g.one, g.constant(N), 'same selected exact positive N'),
            order, rate=1)
        # The source density is defined by the total native derivative of G:
        # B = (G_y_total - G_y_slow) / N.  The production result represents
        # integral K*B/N dy, so this is a direct independent comparison.
        dG_dy = sy.diff(G_expr, n) / jacobian_expr
        B_row = (sy.diff(G_expr, z, order).diff(n) / jacobian_expr -
                 sy.diff(H_expr, z, order)) / N
        direct = sy.integrate(
            direct_kernel * B_row / N * jacobian_expr, (n, 0, q))
        assert sy.simplify(got.node - direct) == 0, order
        assert sy.simplify(sy.diff(G_expr, z, 3)) != 0
        count += 1
    assert sy.simplify(G_expr.subs(n, 0)) != 0
    assert sy.simplify(H_expr.diff(n)) != 0
    assert sy.simplify(dG_dy - (sy.diff(G_expr, n) / jacobian_expr)) == 0
    return dict(ordinary_Z_rows=count, nonzero_inlet=True,
                nonconstant_slow_y_envelope=True,
                nonconstant_original_Jacobian=True,
                third_Z_derivative_verified=True,
                moving_query_endpoint_retained=True,
                direct_measure_Jacobian_once=True, N_squared_factor_verified=True)


def _reject(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('Mean phase averaging scope guard admitted invalid input')


@current.source_precision
def run(owner, packets):
    began = time.monotonic()
    fixture = independent_ibp_fixture()
    owner.assert_graph()
    assert set(packets) == {'first_micro'}
    packet = packets['first_micro']
    view = owner.report(packet)
    assert view['source_family'] == owner.family
    assert len(view['signed_phase_primitive_Z0_to_Z3']) == 4
    assert len(view['signed_slow_y_phase_primitive_Z0_to_Z3']) == 4
    assert len(view['signed_moving_query_IBP_mean_Z0_to_Z3']) == 4
    assert view['original_inlet_memory_retained']
    assert view['endpoint_terms_retained']
    assert view['same_N_squared_IBP_factor']
    assert view['original_Jacobian_applied_once']
    assert view['inner_mean_native_variable'] != view['outer_IBP_native_variable']
    assert not view['genuine_slow_y_source_cap_available']
    assert not view['numerical_phase_averaging_gate']
    assert view['numerical_fallback_direct_old_whole_mean_bounds']
    assert all(not view[key] for key in current.OPEN)

    operations = [node.get('operation') for node in owner.g.nodes]
    assert 'current_original_n1_periodic_phase_primitive' in operations
    assert 'current_original_n1_periodic_phase_primitive_slow_y' in operations
    integrals = [node for node in owner.g.nodes
                 if node.get('operation') == 'definite_integral']
    assert any(node.get('variable') == 'mean_phase_native_first_micro'
               for node in integrals)
    assert any(node.get('bound_inner_mean_native_variable') ==
               'native_first_micro' for node in integrals)
    assert any(node.get('ordinary_slow_Z_derivative_order') == 3
               for node in integrals)

    entry = owner._packets[id(packet)]
    key = 'endpoint_terms_retained'
    old = entry[1]['view'][key]
    entry[1]['view'][key] = False
    _reject(lambda: owner.report(packet))
    entry[1]['view'][key] = old
    _reject(lambda: owner.range_cell(('0', '.5'), 'first_micro',
                                     ('sc_half',), (3, 8)))
    owner.g.symbol('mean_phase_averaging_downstream_child')
    owner.report(packet)

    hashes = dict(owner.hashes)
    for name in (current.NAME, Path(__file__).name):
        hashes[name] = current.sha(name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        input_hashes=hashes,
        checks=dict(
            independent_moving_query_IBP_Z0_to_Z3=fixture,
            signed_original_B_phase_primitive_rows_Z0_to_Z3=True,
            endpoint_and_inlet_memory_terms_retained=True,
            distinct_inner_mean_and_outer_IBP_native_variables=True,
            original_Jacobian_applied_once_in_graph=True,
            same_selected_N_squared_factor_and_global_phase=True,
            no_genuine_mixed_slow_y_cap_fabricated=True,
            direct_old_whole_mean_caps_retained_as_fallback=True,
            complete_pressure_tails_repair_and_global_claims_not_installed=True,
        ),
        **{current.GATE: True},
        **{current.NUMERIC_GATE: False},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(current.c3.encoded(result), separators=(',', ':')) + '\n',
        encoding='utf-8', newline='\n')
    print('PASS_N1_MEAN_PHASE_AVERAGING', fixture, flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit('Use run(owner, packets)')
