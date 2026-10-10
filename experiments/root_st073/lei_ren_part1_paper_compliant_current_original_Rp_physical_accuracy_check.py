"""Independent finite physical ratios, exact affine source identities, scopes."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_physical_accuracy as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected(fn):
    try:
        fn()
    except (ValueError, TypeError, ArithmeticError):
        return True
    raise AssertionError('Invalid physical accuracy input admitted')


def symbolic_polynomial(polynomial):
    return s.Add(*(s.Rational(value.numerator, value.denominator)*
        s.Mul(*(s.Symbol('node_'+str(node)) for node in atoms))
        for atoms, value in polynomial.items()))


def fixture(ctx, coefficient, scale, materialized=False):
    return dict(exact_zero_enclosure=current.ends(coefficient) == (0, 0),
        common_scale_coefficient_enclosure=coefficient, directed_reference_log_scale=scale,
        exact_reference_log_scale_function=-1, ordinary_numeric_materialized=materialized,
        requested_relative_width_satisfied=materialized)


@source_precision
def run(owner, deliveries):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert not owner.acceptance_loaded and not any(raw[key] for key in current.GATES+current.OPEN)
    assert raw['source_family'] == owner.family_record and all(owner.assert_graph().values())
    assert all(raw['input_hashes'][name] == digest for name, digest in owner.before.hashes.items())
    for name in (current.pressure.NAME, current.pressure.RECEIPT, Path(current.__file__).name):
        assert current.sha(name) == raw['input_hashes'][name]
    ctx, independent = owner.ctx, MPIntervalContext()
    independent.dps = 500
    # Direct physical exponential ratios, independent of the logarithmic test.
    finite = []
    for sign in (1, -1):
        coefficient = ctx.mpf([2, '2.01']) if sign > 0 else ctx.mpf(['-2.01', -2])
        scale = ctx.mpf([10, '10.0001'])
        row = fixture(ctx, coefficient, scale, materialized=True)
        result = current.relative_error_budget(ctx, row, Fraction(1, 100))
        a, b = current.ends(coefficient)
        low, high = min(abs(a), abs(b)), max(abs(a), abs(b))
        physical_ratio = independent.mpf(high)*independent.exp(independent.mpf(current.ends(scale)[1]))
        physical_ratio /= independent.mpf(low)*independent.exp(independent.mpf(current.ends(scale)[0]))
        width = physical_ratio-1
        assert current.pressure.amplitude.contains(result['ordinary_relative_width_bound'], width)
        assert result['sign'] == sign and result['factored_relative_width_satisfied']
        assert result['ordinary_numeric_relative_width_satisfied'] and result['ordinary_numeric_delivery_target_satisfied']
        finite.append(dict(sign=sign, independent_physical_ratio_included=True))
    # Actual binary-exact huge log points/boxes need no enormous exponential.
    enormous_point = ctx.mpf(2)**10000
    huge_exact = current.relative_error_budget(ctx, fixture(ctx, ctx.mpf(2), enormous_point), Fraction(1, 100))
    assert huge_exact['ordinary_numeric_relative_width_satisfied']
    assert not huge_exact['ordinary_numeric_materialized'] and not huge_exact['ordinary_numeric_delivery_target_satisfied']
    huge_box = ctx.mpf([2**10000, 2**10001])
    factored = current.relative_error_budget(ctx, fixture(ctx, ctx.mpf([2, '2.01']), huge_box), Fraction(1, 100))
    assert factored['factored_relative_width_satisfied'] and not factored['ordinary_numeric_relative_width_satisfied']
    assert factored['ordinary_relative_width_upper_representation']['physical_exponential_not_materialized']
    crossing = current.relative_error_budget(ctx, fixture(ctx, ctx.mpf([-1, 1]), huge_box), Fraction(1, 100))
    assert crossing['sign'] is None and not crossing['factored_relative_width_satisfied']
    zero = current.relative_error_budget(ctx, fixture(ctx, ctx.mpf(0), huge_box, materialized=True), Fraction(1, 100))
    assert zero['exact_zero_enclosure'] and zero['ordinary_numeric_delivery_target_satisfied']
    unresolved = current.relative_error_budget(ctx, dict(exact_zero_enclosure=False,
        common_scale_coefficient_enclosure=None, ordinary_numeric_materialized=False), Fraction(1, 100))
    assert not unresolved['factored_relative_width_satisfied']
    strict = current.relative_error_budget(ctx, fixture(ctx, ctx.mpf([2, '2.01']), ctx.mpf(0)), Fraction(1, 100000000))
    assert not strict['factored_relative_width_satisfied']
    observations = {}
    totals = dict(rows=0, exact_zero_rows=0, nonzero_factored_target_rows=0,
        ordinary_width_target_rows=0, ordinary_numeric_target_rows=0,
        coefficient_sign_unresolved_rows=0, unresolved_ratio_rows=0)
    proofs = 0
    targets = {name: raw['relative_width_target'] for name in deliveries}
    for name, delivery in deliveries.items():
        observed = owner.evaluate(delivery, targets[name])
        assert current.correlated.report(observed) == raw['actual_physical_accuracy_views'][name]
        original = owner.before.before.evaluate(delivery, targets[name])
        reader = owner.before.amplitude.reader(delivery)
        for section, components in observed['physical_value_rows'].items():
            for component, rows in components.items():
                for label, row in rows.items():
                    old = original['physical_value_rows'][section][component][label]
                    assert row['physical_derivative'] == old['physical_derivative']
                    assert current.correlated.report(row['signed_canonical_scale_groups']) == current.correlated.report(old['signed_canonical_scale_groups'])
                    accuracy = row['physical_accuracy']
                    assert accuracy['exact_reference_scale_retained'] and accuracy['failure_does_not_prove_true_field_error']
                    assert accuracy['coefficient_ledger_keeps_shared_source_dependencies']
                    for term, budget in zip(row['ratio_terms'], accuracy['signed_coefficient_and_ratio_variation_ledger']):
                        assert budget['exact_ratio_function'] == term['exact_ratio_function']
                        if term['directed_ratio_enclosure'] is None:
                            assert budget['variation_bound_unresolved']
                            continue
                        coeff = independent.mpf(term['signed_coefficient'])
                        ratio = independent.mpf(term['directed_ratio_enclosure'])
                        contribution = coeff*ratio
                        assert current.pressure.amplitude.contains(budget['normalized_contribution_enclosure'], contribution)
                        width = independent.mpf(current.ends(contribution)[1])-independent.mpf(current.ends(contribution)[0])
                        assert current.ends(width)[1] <= current.ends(budget['normalized_variation_budget'])[1]
                        if term['method']['kind'] == 'retained_positive_tail':
                            assert current.ends(budget['ratio_enclosure'])[1] > 0 and budget['strictly_positive_tail_upper_retained']
                    if row['exact_zero_enclosure']:
                        assert accuracy['ordinary_numeric_delivery_target_satisfied']
                        continue
                    ledger = accuracy['actual_common_scale_source_ledger']
                    proof = owner.graph.nodes[ledger['exact_affine_logRp_identity']]
                    full = symbolic_polynomial(reader.polynomial(proof['original_full_scale_log']))
                    a = symbolic_polynomial(reader.polynomial(proof['exact_affine_coefficient']))
                    origin = symbolic_polynomial(reader.polynomial(proof['original_logRp']))
                    residual = symbolic_polynomial(reader.polynomial(proof['exact_residual_function']))
                    assert s.expand(a*origin+residual-full) == 0
                    assert s.Symbol('node_'+str(owner.pivot.node)) not in residual.free_symbols
                    expected = []
                    for term in ledger['source_term_width_ledger']:
                        assert current.pressure.amplitude.contains(term['directed_term'], reader.at(term['exact_term_function']))
                        expected.append(symbolic_polynomial(reader.polynomial(term['exact_term_function'])))
                        for leaf in term['original_bound_leaves']:
                            assert leaf['directed_source_bound']._mpi_ == reader.bindings[leaf['node']]._mpi_
                    assert s.expand(s.Add(*expected)-full) == 0
                    assert ledger['shared_dependencies_and_reader_rounding_prevent_additive_error_attribution']
                    proofs += 1
                    common = row.get('common_scale_coefficient_enclosure')
                    if common is not None and accuracy.get('sign'):
                        lo, hi = current.ends(common)
                        m, M = min(abs(lo), abs(hi)), max(abs(lo), abs(hi))
                        independent_relative = (independent.mpf(M)-independent.mpf(m))/independent.mpf(m)
                        assert current.pressure.amplitude.contains(accuracy['coefficient_relative_width_upper'], independent_relative)
                        assert accuracy['absolute_error_representation']['coefficient_width_is_not_absolute_physical_error']
                    assert not accuracy['ordinary_numeric_delivery_target_satisfied']
        for key, value in observed['physical_accuracy_counts'].items():
            totals[key] += value
        packet = owner.velocity_pressure(delivery, targets[name])
        assert set(packet['values']) == {'u', 'v', 'w', 'p'}
        assert not packet['current_P0_numeric_inclusion_pending']
        assert not packet['full_certified_physical_accuracy'] and not packet['unrestricted_physical_point_API']
        assert not any(packet[key] for key in current.OPEN)
        observations[name] = dict(chart=packet['chart'], exact_Z=packet['exact_Z'],
            counts=observed['physical_accuracy_counts'], base_component_accuracy={
                key: {field: value for field, value in row['physical_accuracy'].items()
                      if field != 'actual_common_scale_source_ledger'} for key, row in packet['values'].items()})
    invalid = dict(copied_delivery=rejected(lambda: owner.evaluate(copy.copy(next(iter(deliveries.values()))))),
        caller_error_oracle=rejected(lambda: owner.velocity_pressure(next(iter(deliveries.values())), {'error': [0, 0]})),
        invalid_width_target=rejected(lambda: owner.velocity_pressure(next(iter(deliveries.values())), '0')),
        nonlinear_origin=rejected(lambda: current.affine_origin({(7, 7): Fraction(1)}, {(7,): Fraction(10)}, 7)))
    saved = owner.origin
    try:
        owner.origin = owner.graph.zero
        invalid['substituted_original_radius'] = rejected(owner.assert_graph)
    finally:
        owner.origin = saved
    assert all(invalid.values()) and all(owner.assert_graph().values())
    hashes = dict(owner.hashes)
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(all_passed=True, source_family=owner.family_record,
        actual_row_counts=totals, actual_exact_affine_scale_identities=proofs,
        actual_supported_physical_error_observations=observations,
        independent_500_digit_finite_physical_ratio_checks=finite,
        exact_huge_scale_kept_separate_from_numeric_materialization=True,
        factored_certificate_survives_huge_scale_box_but_not_sign_crossing=True,
        exact_zero_and_unresolved_ratios_preserved=True,
        source_scale_dependency_terms_and_exact_bound_leaves_checked=True,
        actual_signed_coefficient_and_positive_ratio_tail_variation_budgets_checked=True,
        requested_target_is_enclosure_diameter_not_NS_residual_threshold=True,
        failed_width_test_not_claimed_as_true_field_error=True,
        invalid_inputs_rejected=invalid, **dict.fromkeys(current.GATES, True), **dict.fromkeys(current.OPEN, False),
        input_hashes=hashes, execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.correlated.report(receipt), indent=2)+'\n', encoding='utf8', newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_PHYSICAL_ACCURACY exact factors, ordinary error bounds and actual source dependencies', flush=True)
    return receipt
