"""Fresh directed-endpoint checks; no old provider or interval helper."""
from fractions import Fraction


class FakeInterval:
    def __init__(self, lower, upper):
        self._mpi_ = (lower, upper)


def _raises(fn):
    try:
        fn()
    except Exception:
        return True
    return False


def run_fixture(api):
    converter = api.make_converter(version='1.3.0', max_bits=80)
    negative = converter(FakeInterval((1, 7, -2, 3), (1, 5, -2, 3)))
    crossing = converter(FakeInterval((1, 9, -3, 4), (0, 11, -3, 4)))
    zero = converter(FakeInterval((0, 0, 0, 0), (0, 1, -4, 1)))
    checks = {
        'negative_directed_order_correct': negative == (Fraction(-7, 4), Fraction(-5, 4)),
        'cross_zero_exact': crossing == (Fraction(-9, 8), Fraction(11, 8)),
        'zero_and_positive_exact': zero == (Fraction(0), Fraction(1, 16)),
        'fraction_outputs_only': all(isinstance(x, Fraction) for x in negative + crossing + zero),
        'reversed_new_inputs_rejected': _raises(lambda: converter(FakeInterval((0, 11, -3, 4), (1, 9, -3, 4)))),
        'fractional_fields_rejected': _raises(lambda: converter(FakeInterval((0, Fraction(3, 2), 0, 1), (0, 1, 0, 1)))),
        'conversion_scope_not_certified': converter.receipt['runtime_version_verified'] is False
            and converter.receipt['source_identity_verified'] is False
            and converter.receipt['endpoint_rounding_verified'] is False,
    }
    return {'schema': 'original_mpmath_directed_endpoints_fixture_v2',
            'passed': all(checks.values()), 'checks': checks, 'synthetic_only': True,
            'actual_value_selected': False,
            'limitations': ['Fresh synthetic endpoint tuples only; mpmath and source provider were not invoked.']}
