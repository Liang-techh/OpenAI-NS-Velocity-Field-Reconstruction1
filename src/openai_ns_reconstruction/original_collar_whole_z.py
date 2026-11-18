"""Cellwise C1 enclosures for the canonical signed collar m15 primitive.

The pressure contribution is formed from the checked positive beta measures.
Everything is in factored physical units, and all numerical rows are directed
intervals. No live historical source or solver is called.
"""
from __future__ import annotations

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(iv, row):
    endpoints = tuple(tuple(row[f'{side}_exact_mpf_tuple']) for side in ('lower', 'upper'))
    require(all(len(bits) == 4 for bits in endpoints), 'invalid saved interval endpoints')
    value = iv.make_mpf(endpoints)
    require(value.a <= value.b, 'reversed saved endpoints')
    return value


def hull(iv, lo, hi):
    return iv.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def record(value):
    return dict(lower_exact_mpf_tuple=list(value._mpi_[0]), upper_exact_mpf_tuple=list(value._mpi_[1]), directed=True)


def sign(value):
    if value.b < 0:
        return 'negative'
    if value.a > 0:
        return 'positive'
    return 'exact_zero' if value.a == value.b == 0 else 'unresolved'


def _dyadic_inv_p2(iv, log_p):
    """Positive enclosure of Pstar^-2 without evaluating its exponential."""
    def floor_ceil(bits):
        s, man, exp, _ = bits
        require(s == 0 and man > 0, 'positive scale logarithm required')
        if exp >= 0:
            return man << exp, man << exp
        denominator = 1 << (-exp)
        whole, rest = divmod(man, denominator)
        return whole, whole + bool(rest)

    low, high = (2*log_p/iv.ln(2))._mpi_
    floor = floor_ceil(low)[0]
    ceiling = floor_ceil(high)[1]
    require(0 < floor <= ceiling <= floor+2, 'unresolved reciprocal log scale')
    return iv.make_mpf(((0, 1, -ceiling, 1), (0, 1, -floor, 1)))


def _cell(z, data):
    d0 = 1-z**2
    ratio = 1+z**2
    delta, j = data['delta'], data['j']
    gamma = 1+delta
    L, L1 = 1-delta*z**2, -2*delta*z
    require(L.a > 0 and d0.a >= 0, 'real core cell is outside its analytic domain')
    m2, m0, F = data['m2'], data['m0'], data['flat_mass']
    # F, beta_moment and the 6*F bound have beta weights 1, beta,
    # beta*(beta+1) respectively in the same positive integral.
    beta_moment = 2*F
    slope = m2*(2*d0/ratio**3+gamma/ratio**2) + gamma*m0 + d0*beta_moment + gamma*F
    slope1 = -2*z*(m2*(2*(1+gamma)/ratio**3+6*d0/ratio**4) + (1+gamma)*beta_moment + 6*d0*F)
    pressure = z*slope/L
    pressure1 = (slope+z*slope1-pressure*L1)/L

    u = 4*z+j
    H = z*(1-delta)/2+d0*u
    H1 = (1-delta)/2-2*z*u+4*d0
    w = 1-2*z*u
    hydro = (w*u*gamma/2+4*H)/(2*L)
    hydro1 = (((-2*u-8*z)*u+4*w)*gamma/2+4*H1-2*hydro*L1)/(2*L)
    hydro += data['N0']
    hydro1 += data['N1']
    inlet = pressure+data['inv_p2']*hydro
    inlet1 = pressure1+data['inv_p2']*hydro1

    # P/Q = Z*g, and the separate radial correction is Z*r in units Q*S.
    # g is regular at Z=0; no evaluation divides by Z.
    g = -32*inlet/L
    g1 = (-32*inlet1-g*L1)/L
    factor = -2*(1-data['hb'])*data['mass']
    r = factor*g+data['radial_error0']
    r1 = factor*g1+data['radial_error1']
    total, total1 = g+data['S']*r, g1+data['S']*r1
    return {
        'inlet': (inlet, inlet1), 'inlet_pressure': (pressure, pressure1),
        'inlet_hydro_and_nonlinear': (hydro, hydro1),
        'radial': (z*r, r+z*r1), 'primitive': (z*total, total+z*total1),
    }


def _prepare(iv, pressure, pressure_check, core, tube, domain, q2, point):
    p = pressure['compliant_source']
    family = {key: pressure_check[key] for key in ('implicit_source_sha256', 'datum_enclosure_sha256')}
    proof = pressure_check['original_pressure_function_identification']
    require(pressure_check['all_passed'] and proof['passed'] and proof['holomorphic_integral_and_axial5_differentiation_justified']
            and proof['original_implicit_Fflat_is_checked_native_pressure_integral_callable']
            and proof['identities']['same_exact_native_and_original_stage_integral_z_flatten']
            and 'sigma(t/100)*log((Z**2 + 1)**2/4)' in proof['exact_original_Fflat']
            and p['stages']['z_flatten']['beta'] == 'variable [0,2]' and p['stage_count_total'] == 14,
            'actual positive pressure beta measure required')
    require(all(p[k] == core[k] == domain['source_family'][k] == q2['source_family'][k] == value
                for k, value in family.items()) and q2['source_family'] == domain['source_family']
            and tube['analytic_core_family_sha256'] == core['analytic_core_family_sha256']
            and core['contraction_proved'] and tube['common_complex_axis_poles_excluded'],
            'different core or pressure family')
    require(domain['joint_core_domain']['Z_real_centers'] == ['-1', '1'] and q2['checks']['passed'] and point['passed']
            and point['frozen_inputs']['experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json']
            == '420eebe81d021418b869ec4df291409fccb85ce49de35602e740ec72c2727b13', 'whole-axis core tube missing')
    selected = q2['actual_correlated_q2_inlet']
    variation = selected['variation_proof']
    eta = read(iv, domain['common_tubes']['outer_Z_radius'])
    gap = read(iv, variation['complex_Z_gap'])
    lower_L = read(iv, variation['complex_L_modulus_lower']).a
    # The saved K bound used global analytic core/ratio norms. The q2 input
    # selected only the inlet model. For every real center the same axial
    # loss is available on the uniform tube.
    require(0 < gap.a and gap.b*4 < eta.a and lower_L > 0
            and variation['whole_error_Z1_uses_Cauchy']
            and domain['joint_core_domain']['Cauchy_rho_radius'] == '.5'
            and domain['actual_inner_geometry']['actual_bridge_not_replaced_by_comparison']
            and read(iv, domain['actual_inner_geometry']['real_rho_upper_bound']).b < 5,
            'insufficient uniform Cauchy neighborhood')
    S = read(iv, point['geometry']['S'])
    accepted_endpoint = tuple(selected['selected_phase']['lower_exact_mpf_tuple'])
    require(S.a == S.b and S._mpi_[0] == accepted_endpoint and 0 < S.a < iv.mpf(1)/4096,
            'radial endpoint not the accepted exact sc/2')
    log_hb, log_p2 = (read(iv, selected['source_log_bases'][k]) for k in (0, 2))
    require(log_hb.b < -1024 and log_p2.a > 1024, 'inadmissible physical scales')
    K = read(iv, point['bounds']['K_minus_2T0_over_Pstar_squared_upper'])
    require(K.a > 0, 'no certified K remainder')
    e0 = (2*K/lower_L).b  # bound the factored prefactor 2/L; no division by Z
    e1 = (e0/gap.a).b
    true_delta = read(iv, q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    require(0 < true_delta.a and true_delta.b < read(iv, tube['delta']).a, 'actual delta outside common parameter gate')
    core_norm = read(iv, core['scaled_map_size_upper'])/(1-read(iv, core['scaled_map_Lipschitz_upper']))
    h = read(iv, tube['Xh_parameter'])
    denominator = 1-iv.mpf(4)/20
    # Apply the X_h coefficient embedding at rho=4, with ordinary orders (1,0),(1,1).
    N0 = (core_norm/(20*denominator**2)).b
    N1 = (2*core_norm/(80*h*denominator**3)).b
    floor = iv.mpf(0)
    j = read(iv, core['required_j'])
    require(j.a > 0 and read(iv, tube['j'])._mpi_ == j._mpi_ and N0 < j.a, 'selected j or core error differs')
    Fflat = read(iv, p['stages']['z_flatten']['mass'])
    m2, m0 = (read(iv, p[f'fixed_beta{k}_mass_normalized']) for k in (2, 0))
    require(Fflat.a == 0 and Fflat.b > 0 and m2.a > 0 and m0.a >= 0, 'invalid pressure measure enclosures')
    return dict(S=S, delta=true_delta, j=j, m2=m2, m0=m0, flat_mass=Fflat,
                inv_p2=_dyadic_inv_p2(iv, read(iv, core['logPstar'])), N0=hull(iv, -N0, N0), N1=hull(iv, -N1, N1),
                hb=hull(iv, floor, iv.mpf(2)**(-1024)), mass=hull(iv, (1-S*S).a, iv.mpf(1)),
                radial_error0=hull(iv, -e0, e0), radial_error1=hull(iv, -e1, e1))


def evaluate(*, pressure, pressure_check, core, tube, domain, q2, point):
    iv = MPIntervalContext()
    iv.prec = 1800
    data = _prepare(iv, pressure, pressure_check, core, tube, domain, q2, point)
    def encode(z, result):
        return dict(axial_cell=record(z), ordinary_C0_C1={key: [record(v) for v in row] for key,row in result.items()},
                    signs={key: [sign(v) for v in row] for key,row in result.items()})

    cells = []
    for i in range(-32,32):
        z = hull(iv, iv.mpf(i)/32, iv.mpf(i+1)/32)
        cells.append(encode(z, _cell(z, data)))
    queries = {}
    mp = MPContext()
    for key, n, d in (('-0.3',-3,10), ('0',0,1), ('0.3',3,10), ('0.371',371,1000)):
        z = iv.mpf(n)/d
        answer = _cell(z, data)
        queries[key] = encode(z, answer)
        queries[key]['rounded_displays'] = {k: [[mp.nstr(mp.make_mpf(bits), 16) for bits in row._mpi_] for row in answer[k]]
                                            for k in ('primitive','radial')}
    require(queries['0']['signs']['primitive'] == ['exact_zero','negative']
            and queries['0']['signs']['radial'] == ['exact_zero','positive'], 'analytic axis cancellation or slope unresolved')
    require(queries['-0.3']['signs']['primitive'] == ['negative','positive']
            and queries['0.3']['signs']['primitive'] == ['negative','negative'], 'legacy sample signs unresolved')
    return dict(passed=True, admitted_domain=['-1','1'], covers_every_real_Z=True, axial_cell_count=len(cells), cells=cells,
                exact_samples=queries, ordinary_Z_orders=[0,1], units={
                    'inlet': 'Psi_core,rho(4,Z)/Pstar^2', 'primitive': 'epsilon^2*hb*Pstar^2*sigma(S)',
                    'radial': 'epsilon^2*hb*Pstar^2*sigma(S)*S'},
                pressure_beta_masses_retained=True, source_exponentials_materialized=False,
                uniform_ratio_and_core_remainder_bound_reused=True, forward212_identity=False,
                error_terms_removed=False, terminal_matching_closed=False, historical_calculations_replayed=False)
