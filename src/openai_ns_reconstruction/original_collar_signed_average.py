"""Enclose the signed collar primitive and its first radial correction.

Inputs are already authenticated receipts for one selected axial point. The
exponential source atoms remain formal. The returned C0/C1 coefficients are
intervals in explicitly factored units.
"""
from __future__ import annotations

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

BASES = ['hb', 'F0_anchor', 'Pstar_squared', 'Lambda', 'sigma_a', 'sigma_2a']


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _read(iv, row):
    lo, hi = (tuple(row[key]) for key in ('lower_exact_mpf_tuple', 'upper_exact_mpf_tuple'))
    _require(len(lo) == len(hi) == 4, 'expected exact endpoints')
    value = iv.make_mpf((lo, hi))
    _require(value.a <= value.b, 'reversed saved interval')
    return value


def _hull(iv, lo, hi):
    return iv.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _record(value):
    mp = MPContext()
    return {
        'lower_exact_mpf_tuple': list(value._mpi_[0]), 'upper_exact_mpf_tuple': list(value._mpi_[1]),
        'display_endpoints_rounded': [mp.nstr(mp.make_mpf(bits), 16) for bits in value._mpi_],
    }


def _model(iv, row, inverse_pstar2):
    """Read the saved directed inlet model -A0/Pstar², before its radial error."""
    target = (1, 0, 1, -2, 0, 1)
    secondary = (1, 0, 0, -2, 0, 1)
    coefficients = {tuple(t['source_exponents']): _read(iv, t['coefficient_interval'])
                    for t in row['signed_approximation']}
    _require(len(row['signed_approximation']) == 2 and set(coefficients) == {target, secondary},
             'unexpected scaling of the inlet model')
    return coefficients[target] + inverse_pstar2 * coefficients[secondary]


def _psi_bound(iv, terms, S, inverse_pstar2, sigma):
    # The saved positive variation majorant is at least 128*S*B1.
    # Indeed 4*Cw>=64*S and its common multiplier is >=2.
    total = iv.mpf(0)
    for term in terms:
        h, f, p, l, a, b = term['source_exponents']
        coefficient = _read(iv, term['coefficient_interval'])
        _require((h, f, l, a) == (0, 0, 0, 0) and p in (0, 1) and b in (0, 1)
                 and coefficient.a >= 0, 'unexpected positive majorant sector')
        total += coefficient * (1 if p else inverse_pstar2) * (sigma if b else 1)
    return (total / (128*S)).b


def evaluate(*, q2_report, ftc_receipt):
    q2 = q2_report['actual_correlated_q2_inlet']
    lower = q2_report['actual_signed_lower_m15']
    ftc = ftc_receipt['actual_FTC']
    packet = q2['shared_q2_inlet']
    lower_sha = ftc_receipt['raw_input_sha256'][next(p for p in ftc_receipt['raw_input_sha256']
                                                    if p.endswith('n1_lower_collar_functions.json.gz'))]
    archive_name = 'lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz'
    _require(ftc_receipt['passed'] is True and q2_report['checks']['passed'] is True and q2['passed'] is True
             and lower['passed'] is True and lower_sha == q2_report['input_hashes'][archive_name]
             and q2['source_family'] == lower['source_family'] == packet['source_family']
             and ftc['signed_primitive']['integral_C0_C1'] == [157, 158]
             and ftc['radial_FTC']['normalized_A_minus_A0_over_hb_C0_C1'] == [145, 146]
             and q2['actual_primitive_binding']['radial_Psi_order'] == 1
             and q2['actual_primitive_binding']['substitute_exact_sc_half'] is True
             and q2['source_log_basis_names'] == BASES, 'incompatible selected source receipts')
    atom = q2['sigma_S_atom']
    large = q2['variation_proof']
    ratio = lower['complex_error_proofs'][1]
    _require(large['bound_per_hb'] == '4*Cw*B1+32*S*rho_max*Ac*sigma(S)*B1+16*Cd*B2'
             and q2['signed_remainder_model'] == '-hb*epsilon^2*A0*sigma(S)'
             and atom['exact_function'] == 'sigma(s)=exp(odds)/(1+exp(odds)), odds=(1-s)^(-2)-s^(-2)'
             and atom['ordinary_Z_derivative_is_zero'] is True and atom['sigma_is_strictly_positive'] is True
             and atom['monotonicity_proof']['flat_left_endpoint'] == 'sigma(0)=0 from pinned source branch'
             and lower['lower_sc_and_hb_Z_independent'] is True
             and packet['Z'] == '371/1000' and packet['q'] == '2'
             and q2['selected_Z_interval'] == lower['selected_Z_interval'] == packet['selected_Z_interval'],
             'selected point, endpoint, or source convention changed')

    iv = MPIntervalContext()
    iv.prec = 1800
    S = _read(iv, q2['selected_phase'])
    _require(S.a == S.b and 0 < S.a < iv.mpf(1)/4096 and list(S._mpi_[0]) == packet['exact_phase_tuple']
             and atom['argument']['lower_exact_mpf_tuple'] == packet['exact_phase_tuple'],
             'exact positive collar endpoint required')
    zero, one = iv.mpf(0), iv.mpf(1)
    log_hb, _, log_pstar2, _, _, _ = (_read(iv, t) for t in q2['source_log_bases'])
    _require(log_hb.b < -1024 and log_pstar2.a > 1024, 'positive factors require log-scale admission')
    tiny = iv.mpf(2)**(-1024)
    hb = inverse_pstar2 = _hull(iv, zero, tiny)
    sigma = _hull(iv, zero, _read(iv, atom['positive_cap']).b)
    _require(sigma.b < S and (hb*S).b < one/8, 'analytic radius and sigma margin required')
    Lmin = _read(iv, large['complex_L_modulus_lower']).a
    gap_interval = _read(iv, large['complex_Z_gap'])
    gap = gap_interval.a
    _require(0 < Lmin <= 1 and gap > 0 and large['whole_error_Z1_uses_Cauchy'] is True
             and large['complex_Z_gap'] == ratio['whole_error_Z1_gap'], 'complex quotient or C1 tube missing')
    Ac = _read(iv, ratio['core_log_derivative_upper']).b
    B2 = _read(iv, large['core_Psi_rho2_bound']).b * inverse_pstar2
    B1 = _psi_bound(iv, q2['variation_bound_per_hb'], S, inverse_pstar2, sigma)
    _require(Ac > 0 and B1 > 0 and B2.b > 0 and (5*hb*S*Ac*sigma).b < one/2,
             'ratio exponential or core derivative bounds missing')

    # K-2T0 on the same complex axial tube. |rho|<=5, |q|<=2,
    # |rho-4|<=5*hb*S and |q-1|<=10*hb*S*Ac*sigma(S).
    K_error = (180*hb*S + 320*hb*S*Ac*sigma + 250*Ac*sigma)*B1 + (160*hb*S + 250)*B2
    pref_bound = 2*(1+16*gap_interval.b)/Lmin   # eta = 16*gap in the accepted tube
    error = (pref_bound*K_error).b
    # Integrate by parts against sigma', retaining sigma(S) as its own atom.
    # W=1-sigma(t)/sigma(S), M=∫W, 0<=1-M/S<=S² by log-sigma monotonic growth.
    mass = _hull(iv, (one-S*S).a, one)
    z_results = {}
    signs = []
    for order in ('0', '1'):
        inlet = _model(iv, q2['ordinary_Z'][order]['actual_signed_A_minus_A0'], inverse_pstar2)
        radius = error if order == '0' else (error/gap).b
        correction = -2*(one-hb)*mass*inlet + _hull(iv, -radius, radius)
        full = inlet + S*correction
        _require(inlet.b < 0 and correction.a > 0 and full.b < 0, 'signed coefficient unresolved')
        signs.append('negative total, positive radial correction')
        z_results[order] = {'inlet_in_total_unit': _record(inlet),
                            'radial_correction_in_scaled_unit': _record(correction),
                            'radial_correction_relative_width_upper': _record(((correction.b-correction.a)/correction.a).b),
                            'radial_K_error_after_axial_factors': _record(radius),
                            'complete_primitive_in_total_unit': _record(full)}
    return {
        'passed': True, 'selected_Z': packet['Z'], 'selected_Z_interval': q2['selected_Z_interval'],
        'ordinary_Z_orders': [0, 1], 'signs': signs,
        'total_unit': 'epsilon^2*hb*Pstar^2*sigma(S)',
        'radial_correction_unit': 'epsilon^2*hb*Pstar^2*sigma(S)*S',
        'radial_correction_definition': 'signed primitive minus its exact same-source inlet model -epsilon^2*hb*A0*sigma(S)',
        'geometry': {'S': _record(S), 'M_over_S_deficit_upper': _record((S*S).b), 'sigma_S_log_atom': atom},
        'bounds': {'B1_over_Pstar_squared_upper': _record(B1), 'B2_over_Pstar_squared': _record(B2),
                   'K_minus_2T0_over_Pstar_squared_upper': _record(K_error.b)},
        'ordinary_Z': z_results, 'ancestor_calculations_replayed': False,
        'complete_axial_function_admission': False, 'forward212_identity': False,
        'legacy_error_terms_removed': False, 'terminal_matching_closed': False,
    }
