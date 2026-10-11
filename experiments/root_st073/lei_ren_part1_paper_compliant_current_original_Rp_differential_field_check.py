"""Independent derivative semantics, primitive cancellations and linear bounds."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_differential_field as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError, KeyError):
        return True
    raise AssertionError('Invalid original differential field input admitted')


def expression(graph, node, bindings, aliases):
    node = node.node if hasattr(node, 'node') else node
    if node in aliases:
        return expression(graph, aliases[node], bindings, aliases)
    if node in bindings:
        return s.Symbol('source_' + str(node), real=True)
    data = graph.nodes[node]
    op = data['operation']
    recurse = lambda n: expression(graph, n, bindings, aliases)
    if op == 'exact_rational':
        return s.Rational(data['numerator'], data['denominator'])
    if op == 'sum':
        return s.Add(*(recurse(n) for n in data['arguments']))
    if op == 'product':
        return s.Mul(*(recurse(n) for n in data['arguments']))
    if op == 'negative':
        return -recurse(data['argument'])
    if op == 'positive_quotient':
        return recurse(data['numerator']) / recurse(data['denominator'])
    if op == 'analytic_unary' and data['name'] in ('sin', 'cos', 'exp', 'log'):
        return getattr(s, data['name'])(recurse(data['argument']))
    if op == 'exact_operator_integer_power':
        return recurse(data['base']) ** data['exponent']
    return s.Symbol('graph_' + str(node), real=True)


def semantic_forms():
    x, y, z, t = s.symbols('x y z t')
    variables = (x, y, z)
    U = (x*x*y + y*z*z + t*z, x*y*y + z*x*x + t*x, x*y*z + z**3 + t*y)
    jac = s.Matrix(3, 3, lambda i, j: s.diff(U[i], variables[j]))
    expected = dict(divergence=s.trace(jac),
        vorticity_x=jac[2, 1]-jac[1, 2], vorticity_y=jac[0, 2]-jac[2, 0],
        vorticity_z=jac[1, 0]-jac[0, 1])
    for i, component in enumerate(current.COMPONENTS[:3]):
        expected['Laplacian_' + component] = sum(s.diff(U[i], v, 2) for v in variables)
        for j in range(i, 3):
            expected['strain_' + str(i) + str(j)] = (jac[i, j] + jac[j, i]) / 2
    for name, terms in current.standard_forms().items():
        actual = 0
        for component, index, weight in terms:
            value = U[current.COMPONENTS.index(component)]
            for v, order in zip(variables, index):
                value = s.diff(value, v, order)
            actual += s.Rational(weight.numerator, weight.denominator) * value
        assert s.expand(actual - expected[name]) == 0, name
    return len(expected)


@source_precision
def run(owner, fields, deliveries):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[k] for k in current.GATES + current.OPEN)
    assert raw['source_family'] == owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name] == digest for name, digest in owner.before.hashes.items())
    forms_checked = semantic_forms()
    c = MPIntervalContext(); c.dps = 800
    observations = {}; rows_checked = bounds_checked = cancellations_checked = 0
    for name, field in fields.items():
        record = owner._require(field)
        assert record['delivery'] is deliveries[name]
        report = owner.report(field)
        assert not report['unrestricted_physical_point_API'] and not report['full_certified_physical_accuracy']
        assert report['background_momentum_products_stress_and_flat_remainder_not_installed']
        assert report['original_analytic_divergence_identity']['exact_original_Mz_recovery_divergence_identity']
        indices = {(i, j, k) for i in range(5) for j in range(5-i) for k in range(5-i-j)}
        assert set(record['entries']) == {(component, index) for component in current.COMPONENTS for index in indices | {('t',)}}
        for (component, index), entry in record['entries'].items():
            row = entry['row']; handle = entry['handle']
            assert entry['raw'].component == component and entry['raw'].derivative == index
            assert row['physical_derivative'] == index and row['component'] == component
            assert entry['status'] == 'nonzero' and handle is not None
            issued = owner.before._values[id(handle)][1]
            assert issued['unit'] == current.unit(component, index)
            assert issued['scale'].node == row['exact_reference_log_scale_function']
            assert issued['coefficient_tuple'] == row['common_scale_coefficient_enclosure']._mpi_
            assert issued['sign'] == row['physical_accuracy']['sign']
            assert issued['absolute_width_satisfied'] == row['physical_accuracy']['ordinary_numeric_relative_width_satisfied']
            assert not issued['original_materialized']
            assert issued['original_row'] is row
            rows_checked += 1
        independent_reader = current.correlated.CorrelatedGraphBounds(current.correlated.locator.CancelledGraphBounds(
            owner.graph, c, dict(record['bindings']), record['allowed']), record['aliases'])
        for operator_name, terms in current.standard_forms().items():
            result = report['actual_linear_differential_fields'][operator_name]
            primitives = {}
            # Derive source-correlated sums directly from the typed source rows.
            for component, index, q in terms:
                for term in record['entries'][component, index]['raw'].terms:
                    log = owner.graph.add(*(ref for _, ref in term.log_scale_parts))
                    key = (term.source_label, term.source_row.derivative, term.source_row.source_units,
                        term.source_row.powers, tuple(sorted(independent_reader.polynomial(log).items())))
                    part = primitives.setdefault(key, dict(log=log, exact=s.Integer(0),
                        operator=c.mpf(0), source=c.mpf(current.ends(term.source_row.coefficients[0]))))
                    part['exact'] += s.Rational(q.numerator, q.denominator) * expression(
                        owner.graph, term.operator_function, record['bindings'], record['aliases'])
                    part['operator'] += (c.mpf(q.numerator) / q.denominator) * c.mpf(current.ends(term.operator_coefficient[0]))
            for cancelled in result['exact_cancelled_primitive_operators']:
                e = expression(owner.graph, cancelled['exact_combined_operator'], record['bindings'], record['aliases'])
                assert s.cancel(s.expand(e)) == 0
                cancellations_checked += 1
            if result['exact_zero_enclosure']:
                assert all(s.cancel(s.expand(part['exact'])) == 0 for part in primitives.values())
                continue
            scale_groups = {}
            for key, part in primitives.items():
                if s.cancel(s.expand(part['exact'])) == 0:
                    continue
                scale_key = (key[-1], key[2], key[3])
                group = scale_groups.setdefault(scale_key, dict(log=part['log'], coefficient=c.mpf(0)))
                group['coefficient'] += part['source'] * part['operator']
            independent_sum = c.mpf(0)
            for part in scale_groups.values():
                diff = owner.graph.sub(part['log'], current.box.pulse.radius.FunctionRef(
                    owner.graph, result['exact_reference_log_scale_function']))
                log_ratio = independent_reader.at(diff)
                lo, hi = current.ends(log_ratio)
                if not independent_reader.polynomial(diff):
                    ratio = c.mpf(1)
                elif hi <= -1000:
                    ratio = c.mpf([0, current.ends(c.exp(-1000))[1]])
                else:
                    assert -1000 <= lo <= hi <= 1000, operator_name
                    ratio = c.exp(log_ratio)
                independent_sum += part['coefficient'] * ratio
            assert current.product.contains(result['common_scale_coefficient_enclosure'], independent_sum), operator_name
            assert not result['ordinary_numeric_materialized']
            bounds_checked += 1
        differential = report['actual_linear_differential_fields']
        divergence = differential['divergence']['common_scale_coefficient_enclosure']
        assert current.ends(divergence)[0] <= 0 <= current.ends(divergence)[1]
        assert not differential['divergence']['exact_zero_enclosure']
        assert differential['vorticity_z']['exact_primitive_cancellation_count'] == 2
        assert {v['source_label'] for v in differential['vorticity_z']['exact_cancelled_primitive_operators']} == {current.correlated.physical.mixed.UR}
        observations[name] = dict(derivative_rows=144, spatial_rows=140, fixed_x_time_rows=4,
            counts=report['derivative_counts'], linear_fields=len(differential),
            vorticity_signs={key: differential[key]['signed_log_value']['sign'] for key in
                ('vorticity_x', 'vorticity_y', 'vorticity_z')},
            exact_radial_primitives_cancelled_in_axial_curl=2,
            analytic_divergence_identity_retained=True, numeric_divergence_enclosure_contains_zero=True,
            numeric_divergence_not_promoted_to_exact_zero=True,
            nonlinear_background_stress_flat_remainder_recursion_open=True)
    field = next(iter(fields.values())); record = owner._fields[id(field)][1]
    u = owner.derivative(field, 'ux')
    assert current.ends(owner.before.ratio(u, u)['ordinary_numeric_ratio_enclosure']) == (1, 1)
    zero = owner.linear_combination(field, (('ux', (1, 0, 0), 1), ('ux', (1, 0, 0), -1)))
    assert zero['exact_zero_enclosure'] and current.ends(zero['ordinary_numeric_enclosure']) == (0, 0)
    invalid = dict(copied_field=rejected(lambda: owner.report(copy.copy(field))),
        unavailable_order=rejected(lambda: owner.derivative(field, 'ux', (5, 0, 0))),
        unavailable_component=rejected(lambda: owner.derivative(field, 'u_bad')),
        mismatched_units=rejected(lambda: owner.linear_combination(field,
            (('ux', (1, 0, 0), 1), ('p', (1, 0, 0), 1)))))
    key = ('ux', (1, 0, 0)); saved = record['entries'][key]['handle']
    try:
        record['entries'][key]['handle'] = None
        invalid['removed_nonzero_handle'] = rejected(lambda: owner.report(field))
    finally:
        record['entries'][key]['handle'] = saved
    saved = record['aliases']
    try:
        record['aliases'] = {**saved, owner.graph.constant(314159).node: owner.graph.one.node}
        invalid['changed_aliases'] = rejected(lambda: owner.report(field))
    finally:
        record['aliases'] = saved
    saved = record['entries'][key]['raw']
    try:
        record['entries'][key]['raw'] = record['entries']['uy', (1, 0, 0)]['raw']
        invalid['substituted_source_primitive'] = rejected(lambda: owner.report(field))
    finally:
        record['entries'][key]['raw'] = saved
    saved = record['entries'].pop(('p', ('t',)))
    try:
        invalid['missing_time_row'] = rejected(lambda: owner.report(field))
    finally:
        record['entries']['p', ('t',)] = saved
    saved = record['entries'][key]['raw'].terms[0].source_row.coefficients
    try:
        saved.order = 1
        invalid['Taylor_instead_of_ordinary_derivative'] = rejected(lambda: owner.report(field))
    finally:
        saved.order = 0
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes = dict(owner.hashes); hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family_record,
        actual_differential_field_observations=observations, independently_checked_derivative_rows=rows_checked,
        independently_checked_linear_bounds=bounds_checked,
        independently_checked_cancelled_source_operators=cancellations_checked,
        independent_symbolic_cartesian_semantic_forms=forms_checked,
        independent_800_digit_source_correlated_linear_sums=True,
        linear_operator_definitions=current.forms_definition(),
        original_physical_map_definitions=current.correlated.physical.operator_definitions(),
        original_parameters_pressure_scales_derivative_units_and_absolute_flags_retained=True,
        divergence_interval_overlap_not_used_as_exact_identity=True,
        invalid_inputs_rejected=invalid, input_hashes=hashes, execution_seconds=time.monotonic()-began,
        **dict.fromkeys(current.GATES, True), **dict.fromkeys(current.OPEN, False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt), indent=2)+'\n',
        encoding='utf8', newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_DIFFERENTIAL_FIELD original derivative and correlated linear fields', flush=True)
    return receipt
