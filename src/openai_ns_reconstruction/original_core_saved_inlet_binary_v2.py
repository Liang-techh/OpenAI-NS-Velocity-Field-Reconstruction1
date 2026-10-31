"""Extract exact binary endpoints from authenticated, already integrated rows."""

FAMILY_FIELDS = ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256')

def _tuple(row, key):
    value = row[key + '_exact_mpf_tuple']
    if not isinstance(value, list) or len(value) != 4 or any(type(x) is not int for x in value):
        raise ValueError('exact binary endpoint required')
    sign, mantissa, exponent, count = value
    if sign not in (0, 1) or mantissa < 0 or count != mantissa.bit_length():
        raise ValueError('invalid exact binary endpoint')
    return list(value)

def _interval(row):
    if not isinstance(row, dict): raise ValueError('directed interval required')
    return [_tuple(row, 'lower'), _tuple(row, 'upper')]

def saved_inlet(document, check, view):
    if check.get('all_passed') is not True or check.get('current_core_interior_moment_rho4_axial6_and_Euler4_available') is not True:
        raise ValueError('checked saved inlet rows required')
    family = {k: document[k] for k in FAMILY_FIELDS}
    if any(check.get(k) != value for k, value in family.items()): raise ValueError('source family mismatch')
    if view not in ('inlet', 'fresh_inlet'): raise ValueError('invalid saved inlet view')
    packet = document['current_core_interior_views'][view]
    if _interval(packet['rho']) != [[0, 1, 2, 1], [0, 1, 2, 1]]:
        raise ValueError('rho must equal 4')
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
            lo, hi = _interval(error)
            if (lo[0] and lo[1]) or (hi[0] and hi[1]): raise ValueError('nonnegative saved tail required')
            # k! = 1 for k=0,1, so these are ordinary Z derivatives.
            moment['z' + str(k)] = {'finite': _interval(finite), 'error_upper': hi}
        moments[name] = moment
    ratios = packet['original_F0_relative_axial6_ordinary_derivatives']
    if len(ratios) < 2: raise ValueError('missing F0 relative Z1')
    return {'view': view, 'Z': _interval(packet['Z']), 'rho': '4', 'family': family,
            'moments': moments, 'log_derivative': _interval(ratios[1]),
            'F0_value_materialized': False, 'Q4_capture_safe_substitution_proved': False}
