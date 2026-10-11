"""Independent source partitions, nonzero values and same-scale arithmetic."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_scale_value_arithmetic as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError, KeyError):
        return True
    raise AssertionError('Invalid source-scaled arithmetic input admitted')


def symbolic(poly):
    return s.Add(*(s.Rational(q.numerator, q.denominator) * s.Mul(*(s.Symbol('n' + str(n)) for n in atoms))
        for atoms, q in poly.items()))


def interval(ctx, value):
    return ctx.mpf(current.ends(value))


@source_precision
def run(owner, deliveries):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES + current.OPEN)
    assert raw['source_family'] == owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name] == digest for name, digest in owner.before.hashes.items())
    c = MPIntervalContext()
    c.dps = 800
    observations = {}
    partition_count = 0
    first_values = None
    for name, delivery in deliveries.items():
        original = owner.before.velocity_pressure(delivery, '1/1000')
        packet = owner.velocity_pressure(delivery, '1/1000')
        values = packet['values']
        if first_values is None:
            first_values = values
        reader = owner.before.reader(delivery)
        exported = {key: owner.export(value) for key, value in values.items()}
        for key, record in exported.items():
            row = original['values'][key]
            assert record['exact_original_log_scale_function'] == row['exact_reference_log_scale_function']
            assert record['sign'] == row['physical_accuracy']['sign'] and record['proved_nonzero']
            assert record['signed_common_scale_coefficient_enclosure']._mpi_ == row['common_scale_coefficient_enclosure']._mpi_
            assert record['original_absolute_width_target_satisfied'] == row['physical_accuracy']['ordinary_numeric_relative_width_satisfied']
            assert not record['ordinary_numeric_materialized'] and record['absolute_scale_accuracy_not_promoted']
            inverse = record['exact_source_inverse_mu']
            assert not inverse['zero_error_variation'] and inverse['directed_uncertainty_retained']
            q = Fraction(inverse['coefficient']['numerator'], inverse['coefficient']['denominator'])
            assert q == (Fraction(21, 4) if key in ('u', 'v') else 0)
            scale = reader.polynomial(record['exact_original_log_scale_function'])
            fixed = reader.polynomial(record['exact_fixed_origin']['source_function'])
            inv = reader.polynomial(inverse['term_function'])
            residual = reader.polynomial(record['exact_remaining_log_scale_function'])
            assert s.expand(symbolic(scale) - symbolic(fixed) - symbolic(inv) - symbolic(residual)) == 0
            # Inspect the actual reciprocal definition, not just node equality.
            actual = owner.graph.nodes[inverse['source_function']]
            assert actual['operation'] == 'positive_quotient'
            one = owner.graph.nodes[actual['numerator']]
            assert one['operation'] == 'exact_rational' and Fraction(one['numerator'], one['denominator']) == 1
            assert actual['denominator'] == owner.before.amplitude.mu.node
            if q:
                independent = (c.mpf(q.numerator) / q.denominator) / interval(c, reader.at(owner.before.amplitude.mu))
                assert current.product.contains(inverse['directed_term_enclosure'], independent)
                lo, hi = current.ends(inverse['directed_term_enclosure'])
                assert hi > lo > 0
            partition_count += 1
        u, v, w, p = (values[key] for key in ('u', 'v', 'w', 'p'))
        uv = owner.ratio(u, v)
        ur = original['values']['u']; vr = original['values']['v']
        assert not s.expand(symbolic(reader.polynomial(ur['exact_reference_log_scale_function'])) -
                            symbolic(reader.polynomial(vr['exact_reference_log_scale_function'])))
        independent_uv = interval(c, ur['common_scale_coefficient_enclosure']) / interval(c, vr['common_scale_coefficient_enclosure'])
        assert current.product.contains(uv['ordinary_numeric_ratio_enclosure'], independent_uv)
        assert uv['common_scale_cancelled_exactly'] and uv['scale_polynomial_terms'] == 0
        assert current.ends(uv['directed_scale_log_difference']) == (0, 0)
        assert current.ends(uv['directed_scale_ratio_enclosure']) == (1, 1)
        assert uv['sign'] == 1 and uv['exact_ratio_is_nonzero']
        assert current.ends(uv['ordinary_numeric_ratio_enclosure'])[0] > 0
        # An 800-digit direct cot(theta) is a useful independent sanity value;
        # it is not used as the original ratio's defining expression.
        cot = c.cos(c.mpf(7) / 10) / c.sin(c.mpf(7) / 10)
        assert current.product.contains(uv['ordinary_numeric_ratio_enclosure'], cot)
        self_ratio = owner.ratio(u, u)
        twice = owner.multiply_rational(u, 2)
        twice_ratio = owner.ratio(twice, u)
        negative = owner.multiply_rational(p, '-3/2')
        negative_ratio = owner.ratio(negative, p)
        assert current.ends(self_ratio['ordinary_numeric_ratio_enclosure']) == (1, 1)
        assert current.ends(twice_ratio['ordinary_numeric_ratio_enclosure']) == (2, 2)
        assert current.ends(negative_ratio['ordinary_numeric_ratio_enclosure']) == (Fraction(-3, 2), Fraction(-3, 2))
        assert owner.compare(u, u)['order'] == 0
        assert owner.compare(twice, u)['order'] == 1 and owner.compare(u, twice)['order'] == -1
        assert owner.compare(negative, p)['order'] == 1
        assert owner.compare(u, w)['order'] == 1
        assert owner.compare(w, u)['order'] == -1
        tail = owner.ratio(w, u)
        assert tail['ratio_method']['kind'] == 'retained_positive_tail'
        assert tail['exact_ratio_is_nonzero'] and tail['sign'] == 1
        assert current.ends(tail['ordinary_numeric_ratio_enclosure'])[0] == 0
        assert current.ends(tail['ordinary_numeric_ratio_enclosure'])[1] > 0
        huge = owner.ratio(u, p)
        assert huge['ordinary_numeric_ratio_enclosure'] is None and huge['sign'] == -1
        assert huge['exact_ratio_is_nonzero'] and huge['ratio_method']['kind'] == 'unresolved_ratio'
        assert not exported['u']['original_absolute_width_target_satisfied']
        assert not exported['v']['original_absolute_width_target_satisfied']
        assert exported['w']['original_absolute_width_target_satisfied'] and exported['p']['original_absolute_width_target_satisfied']
        observations[name] = dict(nonzero_value_count=4, exact_source_partitions=4,
            u_over_v=uv, inverse_mu_coefficient='21/4',
            shared_source_scale_cancels_before_directed_arithmetic=True,
            independent_800_digit_finite_ratio_included=True,
            rational_multiple_correlation_retained=True,
            positive_tail_ratio_not_replaced_by_zero=True,
            huge_ratio_not_materialized=True,
            absolute_width_flags_retained=dict(u=False, v=False, w=True, p=True))
    u, v, w, p = (first_values[key] for key in ('u', 'v', 'w', 'p'))
    other = current.CurrentOriginalRpScaleValueArithmetic(owner.before, require_checked=False)
    invalid = dict(copied_value=rejected(lambda: owner.export(copy.copy(u))),
        unissued_value=rejected(lambda: owner.export(current.SourceScaledPhysicalValue('u'))),
        foreign_owner=rejected(lambda: other.export(u)),
        mismatched_units=rejected(lambda: owner.compare(u, p)),
        zero_multiplier=rejected(lambda: owner.multiply_rational(u, 0)),
        invalid_target=rejected(lambda: owner.velocity_pressure(next(iter(deliveries.values())), '0')))
    record = owner._values[id(u)][1]
    # A freshly reconstructed reader ignores poisoned caller caches.
    reader = owner._fresh_reader(record)
    reader.memo[record['scale'].node] = owner.ctx.mpf(0)
    reader.polys[record['scale'].node] = {}
    assert owner.ratio(u, v)['common_scale_cancelled_exactly']
    node = record['scale'].node
    saved = copy.deepcopy(owner.graph.nodes[node])
    try:
        owner.graph.nodes[node]['operation'] = 'forged_source_scale'
        invalid['changed_scale_DAG'] = rejected(lambda: owner.export(u))
    finally:
        owner.graph.nodes[node].clear(); owner.graph.nodes[node].update(saved)
    saved = record['bindings'][owner.before.delta.node]
    try:
        record['bindings'][owner.before.delta.node] = owner.ctx.mpf(0)
        invalid['changed_source_binding'] = rejected(lambda: owner.ratio(u, v))
    finally:
        record['bindings'][owner.before.delta.node] = saved
    saved = dict(record['aliases'])
    try:
        record['aliases'][owner.graph.constant(987654321).node] = owner.graph.one.node
        invalid['changed_issued_aliases'] = rejected(lambda: owner.ratio(u, v))
    finally:
        record['aliases'].clear(); record['aliases'].update(saved)
    saved = record['allowed']
    try:
        extra = owner.graph.unary('exp', owner.graph.constant(-30)).node
        assert extra not in saved
        record['allowed'] = frozenset(saved) | {extra}
        invalid['changed_exponential_whitelist'] = rejected(lambda: owner.export(u))
    finally:
        record['allowed'] = saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family_record,
        actual_value_operations=observations, independently_checked_source_partitions=partition_count,
        independent_exact_source_polynomial_cancellation=True,
        independent_800_digit_finite_ratio_arithmetic=True,
        actual_nonzero_values_and_rational_operations_checked=True,
        inverse_mu_uncertainty_and_absolute_width_flags_retained=True,
        positive_tails_and_unresolved_huge_ratios_retained=True,
        invalid_inputs_rejected=invalid, input_hashes=hashes,
        execution_seconds=time.monotonic() - began,
        **dict.fromkeys(current.GATES, True), **dict.fromkeys(current.OPEN, False))
    (current.HERE / current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt), indent=2) + '\n',
        encoding='utf8', newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_SCALE_VALUE_ARITHMETIC nonzero source values and exact shared-scale ratios', flush=True)
    return receipt
