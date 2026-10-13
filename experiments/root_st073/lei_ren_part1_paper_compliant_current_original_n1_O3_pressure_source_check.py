"""Independent checks for the O3 pressure-source adapter.

The production companion binds the accepted O3 source graph.  This checker
keeps its algebraic fixture independent: it builds ``m=C(Z)/R``, recovers
``q`` from the physical equation with ``U=0``, evaluates the pressure source
directly, and integrates the resulting ``N1p``.  It never evaluates a radius,
``Pstar``, ``N``, or a tiny signed coefficient.
"""

import json
from fractions import Fraction
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_O3_pressure_source as current


def _canonical(value):
    """Normalize rational/log expressions without using the production formula.

    SymPy leaves integer logarithms such as ``log(35)`` and ``log(12)``
    unexpanded after ``expand_log``.  The direct integral therefore contains
    an algebraically equal prime-log basis to the production primitive, but a
    structural comparison sees different terms.  Factoring only positive
    rational/integer log arguments gives an exact canonical form and does not
    perform numerical sampling or weaken the identity check.
    """
    expression = s.together(s.expand_log(s.sympify(value), force=True))

    def factor_integer_log(node):
        if node.func != s.log:
            return node
        argument = node.args[0]
        if not (argument.is_Rational or argument.is_Integer):
            return node
        numerator = int(s.numer(argument))
        denominator = int(s.denom(argument))
        if numerator <= 0:
            return node
        result = s.Integer(0)
        for prime, exponent in s.factorint(numerator).items():
            result += exponent * s.log(prime)
        for prime, exponent in s.factorint(denominator).items():
            result -= exponent * s.log(prime)
        return result

    expression = expression.replace(
        lambda node: node.func == s.log, factor_integer_log)
    return s.factor(s.together(expression))


def _equal(left, right):
    return _canonical(left - right) == 0


def _direct_zero_axial(C, Z, R, delta, a, b):
    """Directly recover q and integrate the pressure equation for m=C/R."""
    d = 1 - Z**2
    L = 1 - delta * Z**2
    m = C / R
    U = s.Integer(0)
    q = s.cancel((2 * Z * U - (1 - delta) * Z * m - d * s.diff(m, Z)) / L)
    q_Z = s.diff(q, Z)
    q_y = R * s.diff(q, R)
    S = s.cancel(
        (q + (1 - delta) * Z * q_Z / 2 + q_y) / L
        + q * (q / 2 + q_y)
    )
    N1p = s.cancel(-S / 2 + 2 * s.diff(q, R) + R * s.diff(q, R, 2))
    B = s.cancel(R * q)
    m13 = s.integrate(N1p, (R, a, b))
    m15 = -s.integrate(R * N1p, (R, a, b))
    return dict(C=C, m=m, U=U, q=q, q_Z=q_Z, q_y=q_y, S=S,
                N1p=N1p, B=B, B_Z=s.diff(B, Z), m13=m13, m15=m15)


def _check_general_case(C, Z, R, delta, a, b):
    direct = _direct_zero_axial(C, Z, R, delta, a, b)
    production = current.general_zero_axial_primitives(
        direct['B'], direct['B_Z'], Z, delta, a, b)
    assert _equal(production['m13'], direct['m13'])
    assert _equal(production['m15'], direct['m15'])
    assert _equal(s.diff(production['m13'], Z), s.diff(direct['m13'], Z))
    assert _equal(s.diff(production['m15'], Z), s.diff(direct['m15'], Z))
    # The radial cancellations used by the primitive must also hold in the
    # independently recovered field.
    assert _equal(direct['q_y'], R * s.diff(direct['q'], R))
    assert _equal(direct['N1p'],
                  -direct['S'] / 2 + 2 * s.diff(direct['q'], R)
                  + R * s.diff(direct['q'], R, 2))
    return direct


def independent_pressure_fixture():
    """Use unrelated cubic/quartic own means and exact rational endpoints."""
    Z, R = s.symbols('Z R', real=True, positive=True)
    cases = (
        (s.Rational(1, 7),
         2 + 3 * Z - 5 * Z**2 + Z**3 + Z**4 / 11,
         s.Rational(3, 5), s.Rational(7, 4)),
        (s.Rational(2, 11),
         1 - 2 * Z + 4 * Z**2 + 3 * Z**3 - Z**4 / 13,
         s.Rational(2, 3), s.Rational(5, 3)),
    )
    checked = 0
    normalized = []
    for delta, C, a, b in cases:
        direct = _check_general_case(C, Z, R, delta, a, b)
        checked += 8
        normalized.append((delta, C, direct))

    # Bind the same [Rc exp(-2), Rc] geometry as the accepted O3 source.
    Rc = s.symbols('Rc', positive=True)
    delta, C0, _ = normalized[0]
    C = Rc * C0
    a, b = Rc * s.exp(-2), Rc
    direct = _direct_zero_axial(C, Z, R, delta, a, b)
    btilde = s.cancel(direct['B'] / Rc)
    production = current.normalized_pressure_functions(Z, delta, btilde)
    expected13 = s.cancel(direct['m13'] / Rc)
    expected15 = s.cancel(direct['m15'] / Rc**2)
    assert _equal(production['m13'], expected13)
    assert _equal(production['m15'], expected15)
    assert _equal(s.diff(production['m13'], Z), s.diff(expected13, Z))
    assert _equal(s.diff(production['m15'], Z), s.diff(expected15, Z))
    assert s.simplify(s.diff(btilde, Z)) != 0
    # At Z=0 the derivative term vanishes, so these are direct signed checks
    # for a nonconstant btilde without evaluating its tiny source scale.
    assert s.simplify(expected13.subs(Z, 0)) > 0
    assert s.simplify(expected15.subs(Z, 0)) < 0
    checked += 7
    return dict(general_cases=2, normalized_cases=1, checked_rows=checked,
                normalized_btilde=s.sstr(btilde))


def _find(value, key):
    """Find a source fact in the view/identity packet without inventing it."""
    if isinstance(value, dict):
        if key in value:
            return value[key]
        for child in value.values():
            found = _find(child, key)
            if found is not None:
                return found
    elif isinstance(value, (list, tuple)):
        for child in value:
            found = _find(child, key)
            if found is not None:
                return found
    return None


def _assert_source_facts(owner, view):
    sources = (view, getattr(owner, 'identity', None),
               getattr(owner, 'functions', None))

    def fact(key):
        for source in sources:
            found = _find(source, key)
            if found is not None:
                return found
        raise AssertionError('Missing accepted O3 source fact: ' + key)

    incoming = fact('incoming_node_rows')
    if isinstance(incoming, dict):
        incoming = list(incoming.values())
    assert any(value not in (0, '0', None, False) for value in incoming)
    local = fact('local_density_node_rows')
    assert list(local) == [0, 0, 0, 0]
    assert str(fact('geometry_dlogR_width')) == '2'
    for key in ('original_Rc_phase_offset_shared',
                'original_U_Z3_exact_zero',
                'original_selected_N_division_once',
                'physical_Pstar_applied_once'):
        assert fact(key) is True
    return True


def _assert_signed_view(view):
    assert view['m13']['sign'] == 'positive'
    assert view['m15']['sign'] == 'negative'
    assert 'log_coefficient' in view['m13']
    assert 'log_absolute_coefficient' in view['m15']
    assert view['normalized_coefficient_relative_width_passed'] is True
    assert 'normalized_coefficient_relative_width' in view
    assert view['value_error_relative_to_baseline_upper'] == '1e-60'
    assert _find(view, 'true_log_error') is not None
    assert view['source_scaled_Z1_accuracy'] is False
    assert isinstance(view['source_scaled_Z1_error_magnitudes'], dict)
    assert view['finite_N_incoming_mean_retained'] is True
    return True


def _read_exact_interval(c, row, label):
    """Read a saved interval through its exact MPF tuples.

    The decimal fields in these packets are display data.  Rechecking the
    tuple shape here makes the checker use the directed interval that the
    producer actually saved, without constructing a source value or radius.
    """
    assert isinstance(row, dict), label
    for side in ('lower', 'upper'):
        key = side + '_exact_mpf_tuple'
        value = row.get(key)
        assert isinstance(value, list) and len(value) == 4, label + ':' + key
        assert all(isinstance(item, int) for item in value), label + ':' + key
    return current.read_interval(c, row)


def _read_log_cap(c, row, label):
    assert row.get('exact_zero') is False, label + ': unexpected exact zero'
    return _read_exact_interval(c, row['log_absolute_upper'], label + ':log')


def _exact_mpf_fraction(record):
    """Decode one mpmath binary tuple without lowering its precision."""
    sign, mantissa, exponent, _ = record
    value = Fraction(mantissa) * (Fraction(2) ** exponent)
    return -value if sign else value


def _interval_contains(outer, inner):
    outer_lower, outer_upper = outer._mpi_
    inner_lower, inner_upper = inner._mpi_
    return (_exact_mpf_fraction(outer_lower)
            <= _exact_mpf_fraction(inner_lower)
            and _exact_mpf_fraction(inner_upper)
            <= _exact_mpf_fraction(outer_upper))


def _assert_interval_evidence(owner, view):
    """Recompute signed-value and source-cap evidence from exported tuples."""
    c = MPIntervalContext()
    c.dps = 180
    assert view['source_family'] == owner.family
    for key, raw in owner.sources.items():
        assert raw['source_family'] == owner.family, 'source-family mismatch: ' + key

    # The independent baseline is rebuilt from the accepted quiet K record,
    # rather than trusting the baseline packet copied into the O3 view.
    quiet_raw = owner.sources['quiet']
    quiet_actual = quiet_raw['actual_quiet_source']
    logk = _read_exact_interval(c, quiet_actual['K']['log_coefficient'], 'quiet:K')
    baseline = dict(
        m13=2 * logk + c.ln((c.exp(2) - 1) / 4),
        m15=2 * logk - c.ln(2),
    )
    threshold = c.ln(c.mpf('1e-60'))
    baseline_rows = {}
    value_rows = {}
    for name in ('m13', 'm15'):
        base = baseline[name]
        actual = _read_exact_interval(c, view[name]['log_coefficient'], name + ':value')
        base_lo, base_hi = current.ep(base)
        actual_lo, actual_hi = current.ep(actual)
        assert actual_lo <= base_lo and base_hi <= actual_hi, \
            name + ': signed interval does not extend independent baseline'
        exported_base = _read_exact_interval(
            c, view['baseline_log_coefficients'][name], name + ':baseline')
        exported_lo, exported_hi = current.ep(exported_base)
        assert exported_lo <= base_lo and base_hi <= exported_hi, \
            name + ': exported baseline disagrees with quiet K baseline'

        proof = view['value_error_proofs'][name]
        error = _read_exact_interval(
            c, proof['true_log_error']['log_absolute_upper'],
            name + ':true-log-error')
        recomputed_difference = error - c.mpf(base_lo)
        stored_difference = _read_exact_interval(
            c, proof['log_error_relative_to_baseline_lower'],
            name + ':stored-difference')
        assert current.ep(recomputed_difference)[1] < current.ep(threshold)[0], \
            name + ': recomputed true-error/baseline-lower bound is too large'
        assert current.ep(stored_difference)[1] < current.ep(threshold)[0], \
            name + ': retained true-error/baseline-lower bound is too large'
        baseline_rows[name] = dict(
            independent_lower=str(base_lo), independent_upper=str(base_hi),
            exported_lower=str(exported_lo), exported_upper=str(exported_hi),
            signed_lower=str(actual_lo), signed_upper=str(actual_hi),
            exact_tuple_records_checked=True,
        )
        value_rows[name] = dict(
            recomputed_error_difference_upper=str(current.ep(recomputed_difference)[1]),
            retained_error_difference_upper=str(current.ep(stored_difference)[1]),
            threshold_upper=str(current.ep(threshold)[0]),
        )

    delta_rows = view['source_scaled_delta_b_Z2_magnitudes']
    assert len(delta_rows) == 3
    delta_signs = []
    for index, row in enumerate(delta_rows):
        interval = _read_log_cap(c, row, 'delta_b[%d]' % index)
        lo, hi = current.ep(interval)
        if index < 2:
            assert hi < 0, 'delta_b[%d] cap log is not negative' % index
            delta_signs.append('negative')
        else:
            assert lo > 0, 'delta_b[2] cap log is not positive'
            delta_signs.append('positive')

    # Recompute width from the two exported signed log intervals.  The
    # producer's boolean is checked only after this independent calculation.
    coefficient_intervals = [
        _read_exact_interval(c, view[name]['log_coefficient'], name + ':width')
        for name in ('m13', 'm15')]
    # Keep the subtraction in the interval context.  Converting endpoints via
    # quiet.ep would first round them through the process-wide mpmath context,
    # which can lose the last directed bits at this 180 digit scale.
    endpoint_widths = [interval.b - interval.a
                       for interval in coefficient_intervals]
    recomputed_widths = [c.expm1(width) for width in endpoint_widths]
    exported_width = _read_exact_interval(
        c, view['normalized_coefficient_relative_width'], 'relative-width')
    for recomputed_width in recomputed_widths:
        assert _interval_contains(exported_width, recomputed_width), \
            'saved normalized width does not enclose endpoint recomputation'
    recomputed_upper = max(
        _exact_mpf_fraction(width._mpi_[1]) for width in recomputed_widths)
    threshold_upper = _exact_mpf_fraction(c.mpf('.001')._mpi_[1])
    assert recomputed_upper < threshold_upper
    assert view['normalized_coefficient_relative_width_passed'] is True

    exported_lo, exported_hi = current.ep(exported_width)
    recomputed_lo = min(
        current.ep(width)[0] for width in recomputed_widths)
    recomputed_hi = max(
        current.ep(width)[1] for width in recomputed_widths)

    return dict(
        independent_baseline_logs=baseline_rows,
        true_log_error_vs_baseline_lower=value_rows,
        source_scaled_delta_b_log_signs=delta_signs,
        normalized_width_from_exported_endpoints=dict(
            lower=str(recomputed_lo), upper=str(recomputed_hi),
            exported_lower=str(exported_lo), exported_upper=str(exported_hi),
            target='0.001'),
        exact_mpf_tuple_records_checked=True,
    )


def _assert_parent_scope_and_baseline_binding(owner, view):
    """Check the accepted parent scope and the actual C3/frame baseline link."""
    identity = view['source_identity']
    functions = view['functions']
    assert identity['physical_radius_and_N_not_materialized'] is True
    assert view['physical_radius_scale_retained_as_original_function'] is True
    assert functions['actual_signed_functions_not_caps_or_selected_centers'] is True

    inlet = owner.sources['atlas']['samples']['O3_power']
    scope = inlet['whole_chart_scope']
    assert scope['complete_original_partition'] is True
    assert scope['query_is_contiguous_union_of_original_cells'] is True
    assert scope['own_mean_prefix_starts_at_original_inlet'] is True
    assert scope['incoming_history_is_from_saved_whole_Z_transport'] is True
    assert scope['original_phase_and_Jacobian_Z_independent'] is True
    for key in (
            'same_selected_N_division_once_no_integer_materialization',
            'whole_cell_range_not_outgoing_endpoint_value',
            'magnitude_bounds_not_point_values_and_no_cap_differentiation',
            'partial_mean_free_query_is_distinct_from_bound_integration_coordinate'):
        assert inlet[key] is True, 'parent scope fact missing: ' + key

    # The baseline kernel node is the accepted original definite integral in
    # the Rp/native graph.  The target graph is an exact prefix of that graph;
    # the frame metadata also records that the compact repaired exit is used.
    quiet_identity = owner.sources['quiet']['source_identity']
    frame = owner.sources['frame']
    frame_functions = frame['actual_Rp_native_frame_functions']
    target_nodes = owner.sources['targets']['exact_graph_nodes']
    frame_nodes = frame['exact_graph_nodes']
    kernel_node = functions['baseline_kernel_source_node']
    assert identity['baseline_K_from_same_original_leading_power_and_quiet_source'] is True
    assert identity['baseline_source_identity'] == quiet_identity
    assert functions['baseline_kernel_graph_namespace'] == current.quiet.prior.FRAME_NAME
    assert kernel_node == quiet_identity['actual_kernel_node']
    assert frame_nodes[:len(target_nodes)] == target_nodes
    assert 0 <= kernel_node < len(frame_nodes)
    assert frame_nodes[kernel_node]['original_kernel_function'] is True
    assert frame_nodes[kernel_node]['coordinate_and_kernel_Z_independent'] is True
    assert quiet_identity['actual_compact_repaired_quiet_source_consumed'] is True
    assert frame_functions['source_identity_requires_actual_compact_exit_not_an_outer_shortcut'] is True
    assert frame_functions['original_Rc_phase_offset'] == identity['original_Rc_phase_offset_node']
    assert owner.sources['quiet']['actual_quiet_source']['original_log_Rc_source_node'] == \
        frame_functions['original_logRc']
    leading_histories = frame_functions['actual_leading_2Rc_histories_C3']
    assert all(isinstance(rows, list) and rows for rows in leading_histories.values())
    return dict(
        source_family_matches_all_dependencies=True,
        raw_parent_scope_checked=True,
        no_radius_or_N_materialization=True,
        baseline_kernel_matches_quiet_and_shared_C3_frame_prefix=True,
        baseline_kernel_node=kernel_node,
    )


def _cell_key(cell):
    return tuple(str(value) for value in cell)


def _collect_evaluated_cells(owner, default_view):
    """Normalize the in-memory result of evaluate_all without rebuilding it."""
    raw = getattr(owner, 'evaluated_cells', None)
    if raw is None:
        return [(None, default_view)]
    if isinstance(raw, dict):
        items = list(raw.items())
    else:
        items = list(raw)
    collected = []
    for item in items:
        key = None
        candidate = item
        if isinstance(raw, dict):
            key, candidate = item
        elif isinstance(item, (tuple, list)) and len(item) == 2 and isinstance(item[1], dict):
            key, candidate = item
        assert isinstance(candidate, dict), 'evaluate_all cell must be a view dictionary'
        cell = candidate.get('exact_Z_cell', key)
        assert cell is not None, 'evaluate_all view must retain exact_Z_cell'
        collected.append((cell, candidate))
    return collected


def _assert_four_cell_ranges(owner, views):
    """Check the accepted four-Z-cell V/B zero source evidence."""
    expected = (('-1', '-.5'), ('-.5', '0'), ('0', '.5'), ('.5', '1'))
    expected_set = set(expected)
    seen = [_cell_key(cell) for cell, _ in views]
    assert set(seen) == expected_set and len(seen) == len(expected), \
        'evaluate_all must cover the four accepted O3 Z cells exactly once'

    identity_cover = owner.identity.get('exact_Z_cover')
    if identity_cover is not None:
        assert {_cell_key(cell) for cell in identity_cover} == expected_set
        assert len(identity_cover) == len(expected)

    ranges = owner.sources['ranges']
    assert ranges['source_family'] == owner.family
    raw_cells = ranges['actual_four_Z_C3_target_transports']
    assert len(raw_cells) == len(expected)
    raw_seen = set()
    for row in raw_cells:
        cell = _cell_key(row['exact_Z_cell'])
        raw_seen.add(cell)
        assert cell in expected_set
        assert row['range_caps_not_point_function_values'] is True
        assert row['original_independent_P0_in_source_and_not_added_twice'] is True
        transport = next(item for item in row['actual_17_chart_C3_transports']
                         if item['chart'] == 'O3_power')
        assert len(transport['cells']) == 2
        for source_cell, expected_lr in zip(transport['cells'],
                                            (((0, 1), (1, 1)), ((1, 1), (2, 1)))):
            assert (tuple(source_cell['exact_left']), tuple(source_cell['exact_right'])) == expected_lr
            ranges_for_source = source_cell['actual_source_and_primitive_range']
            for name in ('V', 'B'):
                rows = ranges_for_source[name]
                assert len(rows) == 4
                assert all(bound['exact_zero'] for bound in rows), \
                    'O3 four-Z source %s row is not exact zero' % name
    assert raw_seen == expected_set
    # The raw range packet deliberately records candidate gates as false.  The
    # accepted state is the independently hashed checker receipt consumed by
    # the production constructor.
    ranges_receipt = json.loads(
        (current.rg.HERE / current.rg.RECEIPT).read_text(encoding='utf-8'))
    assert ranges_receipt['all_passed'] is True
    assert ranges_receipt[
        'current_original_actual_four_Z_C3_five_target_ranges_installed'] is True
    return dict(
        exact_four_Z_cover=[list(cell) for cell in expected],
        all_four_O3_V_and_B_Z3_rows_exact_zero=True,
        range_caps_are_not_point_values=True,
    )


def _assert_cell_view_scope(cell, view):
    key = _cell_key(cell)
    assert _cell_key(view['exact_Z_cell']) == key
    scope = str(view['interval_scope'])
    assert 'full actual O3_power segment' in scope
    for endpoint in key:
        assert endpoint in scope, 'interval scope lost endpoint ' + endpoint
    identity_cell = view.get('source_identity', {}).get('exact_Z_cell')
    if identity_cell is not None:
        assert _cell_key(identity_cell) == key
    assert not view.get(current.DERIVATIVE_GATE, False)
    assert all(not view.get(name, False) for name in current.OPEN)


@current.source_precision
def run(owner, view):
    began = time.monotonic()
    owner.assert_source()
    fixture = independent_pressure_fixture()
    evaluated = _collect_evaluated_cells(owner, view)
    range_evidence = _assert_four_cell_ranges(owner, evaluated)
    default_key = ('0', '.5')
    assert _cell_key(view['exact_Z_cell']) == default_key
    per_cell_evidence = {}
    for cell, cell_view in evaluated:
        _assert_cell_view_scope(cell, cell_view)
        _assert_source_facts(owner, cell_view)
        _assert_signed_view(cell_view)
        per_cell_evidence['[%s,%s]' % _cell_key(cell)] = _assert_interval_evidence(owner, cell_view)
    interval_evidence = per_cell_evidence['[0,.5]']
    scope_evidence = _assert_parent_scope_and_baseline_binding(owner, view)

    hashes = dict(owner.hashes)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    checks = dict(
        independent_cubic_quartic_direct_pressure_fixture=fixture,
        direct_physical_S_and_N1p_equation_integrated_dR_and_RdR=True,
        general_zero_axial_primitive_value_and_Z1_rows_verified=True,
        actual_O3_Rc_exp_minus2_to_Rc_normalization_value_and_Z1_verified=True,
        nonconstant_btilde_Z_derivative_and_m13_plus_m15_minus_verified=True,
        accepted_incoming_rows_nonzero_and_local_density_rows_exact_zero=True,
        accepted_geometry_width_two_and_Rc_phase_offset_shared=True,
        selected_N_division_once_and_Pstar_applied_once=True,
        signed_value_records_checked_without_materializing_tiny_coefficients=True,
        source_scaled_Z1_gate_remains_open_with_error_magnitudes=True,
        interval_evidence=interval_evidence,
        all_evaluated_cell_interval_evidence=per_cell_evidence,
        four_cell_source_range_evidence=range_evidence,
        parent_scope_and_baseline_binding=scope_evidence,
    )
    result = dict(
        all_passed=True,
        source_family=owner.family,
        input_hashes=hashes,
        checks=checks,
        **{current.GATE: True},
        **{current.DERIVATIVE_GATE: False},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(result, separators=(',', ':')) + '\n',
        encoding='utf-8', newline='\n')
    print('PASS_O3_PRESSURE_SOURCE', fixture['checked_rows'], flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit('Use run(accepted_o3_owner, accepted_o3_view)')
