"""Fresh synthetic checks for rational exponential enclosures."""
from fractions import Fraction


def _raises(fn):
    try:
        fn()
    except Exception:
        return True
    return False


def run_fixture(api):
    I = api.Interval
    backend = api.Backend()
    checks = {}
    zero = backend.unary('exp', I(0, 0))
    checks['exp_zero_exact'] = zero.lo == zero.hi == 1

    # Independently use a 40-term positive Taylor sum as a lower bound and
    # twice its next term as an upper-tail bound for exp(1).
    term = Fraction(1)
    partial = Fraction(1)
    for k in range(1, 41):
        term /= k
        partial += term
    next_term = term / 41
    reference_lo, reference_hi = partial, partial + 2 * next_term
    one = backend.unary('exp', I(1, 1))
    checks['exp_one_contains_independent_bounds'] = one.lo <= reference_lo and one.hi >= reference_hi

    minus_one = backend.unary('exp', I(-1, -1))
    product = minus_one * one
    checks['opposite_exponents_enclose_one'] = product.lo <= 1 <= product.hi
    spanning = backend.unary('exp', I(Fraction(-1, 2), Fraction(1, 2)))
    checks['spanning_interval_contains_one'] = spanning.lo < 1 < spanning.hi
    checks['monotonic_endpoints'] = spanning.lo < zero.lo and spanning.hi > zero.hi
    checks['float_endpoint_rejected'] = _raises(lambda: I(0.0, 1))
    checks['exp_large_argument_budget_rejected'] = _raises(lambda: backend.unary('exp', I(1025, 1025)))
    checks['unsupported_unary_rejected'] = _raises(lambda: backend.unary('log', I(1, 1)))
    checks['interval_not_scalar'] = _raises(lambda: float(one))

    # A new grid and function, distinct from the v1 synthetic quadrature.
    def polynomial(cell):
        return cell * cell + 3 * cell - 2
    coarse = api.range_integrate(polynomial, 1, 3, partitions=5)
    fine = api.range_integrate(polynomial, 1, 3, partitions=20)
    exact = Fraction(50, 3)
    checks['fresh_polynomial_integral_enclosed'] = (
        coarse.interval.lo <= exact <= coarse.interval.hi
        and fine.interval.lo <= exact <= fine.interval.hi
    )
    checks['fresh_grid_refinement_narrows'] = fine.interval.width < coarse.interval.width
    checks['conditional_receipt_retained'] = (
        fine.kind == 'conditional_enclosure' and fine.coverage_count == 20
        and fine.provider_certified is False and fine.actual_value_selected is False
    )
    return {'schema': 'original_q4_fraction_interval_fixture_v2',
            'passed': all(checks.values()), 'checks': checks,
            'limitations': ['Synthetic exact rational cases only.',
                            'No original source or integral certification is claimed.']}
