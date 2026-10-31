"""Adapt already integrated inlet moment enclosures without recomputation.

Caller authenticates all raw bytes and their lineage. Binary mpf endpoints are
rounded outwards to a dyadic grid, including exponents too small to expand.
"""
from fractions import Fraction

BITS = 1024
FAMILY_FIELDS = ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256')

def outward_mpf_tuple(row):
    if not isinstance(row, list) or len(row) != 4 or any(type(x) is not int for x in row):
        raise ValueError('finite exact mpf tuple required')
    sign, mantissa, exponent, count = row
    if sign not in (0, 1) or mantissa < 0 or count != mantissa.bit_length():
        raise ValueError('invalid exact mpf tuple')
    if mantissa == 0:
        return Fraction(0), Fraction(0)
    shift = exponent + BITS
    if shift >= 0:
        if count + shift > 8192: raise ValueError('endpoint exceeds conversion limit')
        q, remainder = mantissa << shift, False
    else:
        right = -shift
        if right >= count:
            q, remainder = 0, True
        else:
            q = mantissa >> right
            remainder = (q << right) != mantissa
    if sign:
        lower, upper = -(q + int(remainder)), -q
    else:
        lower, upper = q, q + int(remainder)
    denominator = 1 << BITS
    return Fraction(lower, denominator), Fraction(upper, denominator)

def interval(row):
    if not isinstance(row, dict): raise ValueError('directed endpoint record required')
    lo = outward_mpf_tuple(row['lower_exact_mpf_tuple'])[0]
    hi = outward_mpf_tuple(row['upper_exact_mpf_tuple'])[1]
    if lo > hi: raise ValueError('reversed directed endpoints')
    return lo, hi

def encode_interval(row):
    lo, hi = interval(row)
    return [str(lo), str(hi)]

def saved_inlet(document, check, view='inlet'):
    """Read selected saved radial-zero ordinary Z0/Z1 inlet rows only."""
    if check.get('all_passed') is not True or check.get('current_core_interior_moment_rho4_axial6_and_Euler4_available') is not True:
        raise ValueError('checked inlet moment and tail rows required')
    family = {k: document[k] for k in FAMILY_FIELDS}
    if any(check.get(k) != value for k, value in family.items()): raise ValueError('source family mismatch')
    if view not in ('inlet', 'fresh_inlet'): raise ValueError('selected fixed-radius inlet view required')
    packet = document['current_core_interior_views'][view]
    if interval(packet['rho']) != (Fraction(4), Fraction(4)):
        raise ValueError('rho=4 required')
    if packet['radial_degree'] != 24 or packet['axial_depth'] != 6:
        raise ValueError('unexpected stored truncation')
    if packet['finite_polynomial_is_not_complete_solution'] is not True or packet['source_enclosures_not_resolved_point_values'] is not True:
        raise ValueError('stored enclosure semantics required')
    moments = {}
    for name in ('H', 'M', 'K', 'A', 'B', 'C'):
        moment = {}
        for k in (0, 1):
            key = 'rho0_Z' + str(k)
            finite = packet['finite_integrated_axial_Taylor_coefficient_grids'][name][key]
            error = packet['directed_integrated_axial_Taylor_tail_bounds'][name][key]
            error_lo, error_hi = interval(error)
            if error_lo < 0 or error_hi < 0: raise ValueError('nonnegative saved tail required')
            # k! = 1 for k=0,1: these rows are already ordinary derivatives.
            moment['z' + str(k)] = {'finite': encode_interval(finite), 'error_upper': str(error_hi)}
        moments[name] = moment
    ratios = packet['original_F0_relative_axial6_ordinary_derivatives']
    if len(ratios) < 2: raise ValueError('missing F0 relative Z derivative')
    return {'view': view, 'Z': encode_interval(packet['Z']), 'rho': '4', 'family': family,
            'moments': moments, 'log_derivative': encode_interval(ratios[1]),
            'dyadic_grid_bits': BITS, 'outward_rounding': True,
            'F0_value_materialized': False, 'Q4_capture_safe_substitution_proved': False}
