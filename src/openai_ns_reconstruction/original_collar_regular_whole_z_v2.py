"""Directed whole-axis enclosures of the fourth forward m15 source."""
from __future__ import annotations

import hashlib
import json
from math import factorial

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from .original_collar_whole_z import _dyadic_inv_p2, hull, read, require, sign


BASIS = ['hb', 'F0_anchor', 'Pstar_squared', 'Lambda', 'sigma_a', 'sigma_2a']
UNIT = 'epsilon^2*hb*Pstar^2*S, S=sc/2'


def product(a, b):
    """Product jet in ordinary derivative convention."""
    a0, a1, a2 = a
    b0, b1, b2 = b
    return a0*b0, a1*b0+a0*b1, a2*b0+2*a1*b1+a0*b2


def quotient(a, b):
    """Solve a=b*q at orders zero, one and two."""
    b0, b1, b2 = b
    require(b0.a > 0, 'axial quotient denominator is not positive')
    q0 = a[0]/b0
    q1 = (a[1]-q0*b1)/b0
    q2 = (a[2]-2*q1*b1-q0*b2)/b0
    return q0, q1, q2


def inlet_second_jet(z, data):
    """The actual pressure, hydro and fixed-point contributions at rho=4."""
    zero, one = data['zero'], data['one']
    delta, j = data['delta'], data['j']
    d = (1-z**2, -2*z, -2*one)
    x = (1+z**2, 2*z, 2*one)
    L = (1-delta*z**2, -2*delta*z, -2*delta)
    gamma = 1+delta
    x2 = product(x, x)
    x3 = product(x2, x)
    i0, i1, i2, i3 = data['flat_moments']
    # I_k = integral (beta)_k x^(-beta-k) dmu, I_k' = -2 Z I_(k+1).
    flat = (i0, -2*z*i1, 4*z**2*i2-2*i1)
    flat1 = (i1, -2*z*i2, 4*z**2*i3-2*i2)
    rational = quotient(d, x3)
    reciprocal = quotient((one, zero, zero), x2)
    weighted_flat = product(d, flat1)
    slope = tuple(data['m2']*(2*r+gamma*q)+f+gamma*g
                  for r,q,f,g in zip(rational, reciprocal, weighted_flat, flat, strict=True))
    slope = (slope[0]+gamma*data['m0'], *slope[1:])
    z_jet = (z, one, zero)
    pressure = quotient(product(z_jet, slope), L)

    u = (4*z+j, 4*one, zero)
    du = product(d, u)
    H = tuple(v+(1-delta)*w/2 for v,w in zip(du, z_jet, strict=True))
    zu = product(z_jet, u)
    weight = (1-2*zu[0], -2*zu[1], -2*zu[2])
    wu = product(weight, u)
    hydro = quotient(tuple(gamma*v/4+2*w for v,w in zip(wu, H, strict=True)), L)
    return tuple(p+data['inv_p2']*(h+n)
                 for p,h,n in zip(pressure, hydro, data['nonlinear'], strict=True))


def _viscous_model(z, inlet, delta):
    p, pz, pzz = inlet
    L = 1-delta*z**2
    d = 1-z**2
    value = 16*(d*pz-(3+delta)*z*p)/L
    derivative = (16*(d*pzz-2*z*pz-(3+delta)*(p+z*pz))+2*delta*z*value)/L
    return value, derivative


def _viscous_radii(iv, data, packet, point):
    """Apply two nested Cauchy losses to the regular viscous difference."""
    bounds = point['bounds']
    T0 = 16*read(iv, bounds['B1_over_Pstar_squared_upper']).b
    k_error = read(iv, bounds['K_minus_2T0_over_Pstar_squared_upper']).b
    gap, zmax, dmax, Lmin = (data[k] for k in ('gap', 'zmax', 'dmax', 'Lmin'))
    sigma = read(iv, packet['sigma_S_atom']['positive_cap']).b
    hb_upper = iv.mpf(2)**(-1024)
    require(min(T0.a, k_error.a, sigma.a) > 0 and sigma < data['S'].a, 'inadmissible viscous bounds')
    # T_s=hb K. On a smaller complex disk, differentiate T-T0 once
    # before applying Cauchy to the entire difference for its ordinary Z1.
    diff_T = hb_upper*data['S']*(2*T0+k_error)
    inlet_numerator = (3+data['delta'].b)*zmax*T0+dmax*T0/gap
    remainder = ((1+data['delta'].b)*zmax*diff_T+dmax*diff_T/gap
                 +2*zmax*k_error+sigma*inlet_numerator)/Lmin
    return remainder.b, (remainder/gap).b


def _prepare(iv, inputs):
    pressure, core, tube = (inputs[name] for name in ('pressure', 'core', 'tube'))
    domain, q2, point = (inputs[name] for name in ('domain', 'q2', 'point'))
    proof = inputs['pressure_check']
    accepted_measure = proof['original_pressure_function_identification']
    family = q2['source_family']
    p = pressure['compliant_source']
    whole = inputs['whole']
    require(whole['passed'] and whole['covers_every_real_Z'] and point['passed']
            and core['contraction_proved'] and q2['checks']['passed']
            and domain['source_family'] == family == inputs['domain_check']['source_family']
            and tube['analytic_core_family_sha256'] == core['analytic_core_family_sha256']
            and all(p[k] == core[k] == family[k] == proof[k]
                    for k in ('implicit_source_sha256', 'datum_enclosure_sha256')),
            'accepted canonical family mismatch')
    require(proof['all_passed'] and accepted_measure['passed']
            and accepted_measure['holomorphic_integral_and_axial5_differentiation_justified']
            and accepted_measure['identities']['same_exact_native_and_original_stage_integral_z_flatten']
            and accepted_measure['original_implicit_Fflat_is_checked_native_pressure_integral_callable']
            and p['stage_count_total'] == 14 and p['stages']['z_flatten']['beta'] == 'variable [0,2]',
            'the actual positive pressure measure is required')
    packet = q2['actual_correlated_q2_inlet']
    require(packet['source_log_basis_names'] == BASIS, 'physical scale order changed')
    log_hb, _, log_p2, log_lambda, _, _ = (read(iv, r) for r in packet['source_log_bases'])
    S = read(iv, point['geometry']['S'])
    zero, one = iv.mpf(0), iv.mpf(1)
    require(S.a == S.b and 0 < S.a < one/4096 and S._mpi_[0] == tuple(packet['selected_phase']['lower_exact_mpf_tuple'])
            and log_hb.b < -1024 and log_p2.a > 1024, 'wrong collar endpoint or log scales')
    eta = read(iv, domain['common_tubes']['outer_Z_radius'])
    gap = read(iv, packet['variation_proof']['complex_Z_gap']).a
    Lmin = read(iv, domain['denominator_bounds']['L_modulus_lower']).a
    rho_max = read(iv, domain['actual_inner_geometry']['real_rho_upper_bound']).b
    delta = read(iv, q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    require(0 < gap and 4*gap < eta.a and 0 < Lmin and 4 < rho_max < 5
            and 0 < delta.a and delta.b < read(iv, tube['delta']).a
            and domain['joint_core_domain']['Z_real_centers'] == ['-1', '1']
            and domain['actual_inner_geometry']['actual_bridge_not_replaced_by_comparison']
            and packet['variation_proof']['whole_error_Z1_uses_Cauchy'],
            'insufficient whole-axis analytic neighborhood')

    core_norm = read(iv, core['scaled_map_size_upper'])/(1-read(iv, core['scaled_map_Lipschitz_upper']))
    h = read(iv, tube['Xh_parameter'])
    nonlinear = []
    for k in range(3):
        # Cauchy coefficient embedding for ordinary (rho,Z) orders (1,k).
        prefactor = factorial(1+k)/(20*(k+1)**2*h**k*(1-one/5)**(2+k))
        radius = (core_norm*prefactor).b
        nonlinear.append(hull(iv, -radius, radius))
    flat_mass = read(iv, p['stages']['z_flatten']['mass'])
    require(flat_mass.a == 0 and flat_mass.b > 0, 'invalid pressure measure magnitude')
    moments = tuple(hull(iv, zero, (factorial(k+1)*flat_mass).b) for k in range(4))
    zmax = (1+eta).b
    return dict(zero=zero, one=one, S=S, gap=gap, Lmin=Lmin, eta=eta, zmax=zmax,
                dmax=(1+zmax**2).b, delta=delta, rho_max=rho_max,
                j=read(iv, core['required_j']), m2=read(iv, p['fixed_beta2_mass_normalized']),
                m0=read(iv, p['fixed_beta0_mass_normalized']), nonlinear=tuple(nonlinear),
                flat_moments=moments, log_hb=log_hb, log_p2=log_p2, log_lambda=log_lambda,
                inv_p2=_dyadic_inv_p2(iv, read(iv, core['logPstar'])))


def _compare_forward(iv, encode, data, inputs, exact_samples):
    quiet, inner, singular = (inputs[k] for k in ('quiet', 'inner', 'singular'))
    family = inputs['q2']['source_family']
    require(quiet['checks']['passed'] and inner['checks']['passed'] and singular['passed']
            and quiet['source_family'] == inner['source_family'] == family, 'forward receipts differ')
    comparisons = {}
    for key, label in (('-.3', '-0.3'), ('0', '0'), ('.3', '0.3')):
        terminal, = (p for p in quiet['samples'] if
                     (p['Z'], p['chart'], p['native_coordinate']) == (key, 'quiet', '1'))
        ledger = terminal['m15_C0_C1_terms']
        parents = [p for p in inner['samples'] if (p['Z'], p['q']) == (key, '2')
                   and p['m15_C0_C1_terms'] == ledger[:6]]
        require(parents, 'missing archived q2 prefix')
        parent = parents[0]
        raw_digest = hashlib.sha256(json.dumps(ledger, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        require(len(ledger) == 212 and ledger[:6] == parent['m15_C0_C1_terms']
                and raw_digest == singular['samples'][key]['all_forward_terms_digest'], 'wrong forward ledger')
        step = parent['exact_phase_tuple']
        require(data['S']._mpi_[0] == (step[0], step[1], step[2]+1, step[3])
                and parent['exact_sc_tuple'] == [step[0], step[1], step[2]+2, step[3]], 'wrong chart width')
        row = ledger[3]
        require(row['name'] == 'first_lower_regular_pressure_width'
                and len(row['log_scale_terms']) == 4
                and row['log_scale_terms'][2].get('source_log_atom_only') is True
                and row['actual_function_binding']['leading_N1p_after_velocity_support_retained']
                and row['actual_function_binding']['constant_core_pressure_cancels_exactly'],
                'wrong archived regular source')
        norm_logs = (-2*data['log_lambda'], data['log_p2'], data['log_hb'], iv.ln(data['S']/2))
        records = []
        for index, (old, new) in enumerate(zip(row['coefficient_C0_C1'], exact_samples[label], strict=True)):
            box = iv.make_mpf(tuple(tuple(bits) for bits in old['source_enclosure_exact_mpf_tuples']))/2
            model = iv.mpf(old['approximation'])/2
            records.append(dict(ordinary_Z_order=index, saved_finite_coefficient=encode(model),
                saved_directed_coefficient=encode(box), canonical_minus_saved_finite=encode(new-model),
                difference_between_enclosures=encode(new-box),
                directed_intervals_overlap=bool(new.a <= box.b and box.a <= new.b),
                saved_enclosure_contains_canonical=bool(box.a <= new.a and new.b <= box.b)))
        log_diagnostics = []
        for atom, exact in zip(row['log_scale_terms'], norm_logs, strict=True):
            old_log = iv.make_mpf(tuple(tuple(bits) for bits in atom['source_enclosure_exact_mpf_tuples']))
            log_diagnostics.append(encode(exact-old_log))
        comparisons[label] = dict(ledger_count=len(ledger), ledger_sha256=raw_digest, raw_row=row,
            canonical_and_archived_coefficient_unit=UNIT, archived_scale_divisor=2,
            exact_minus_saved_source_log_enclosures=log_diagnostics, ordinary_Z=records)
    return comparisons


def _positive_log_cap(iv, log_value):
    """Round a log upper bound up to a power of two without exponentiating it."""
    binary_log = log_value/iv.ln(2)
    # Powers of hb can put even log_value outside any materializable
    # integer range. Keep a positive, directed dyadic cap in that case.
    minimum_power = -10**20
    if binary_log.b < minimum_power:
        return iv.make_mpf(((0, 1, minimum_power, 1),)*2)
    require(binary_log.b < 10**20, 'normalized source bound did not remain finite')
    negative, mantissa, exponent, _ = binary_log._mpi_[1]
    integer = -mantissa if negative else mantissa
    power = integer << exponent if exponent >= 0 else -((-integer) >> -exponent)
    return iv.make_mpf(((0, 1, power, 1),)*2)


def _inertial_radii(iv, data, inputs):
    p = inputs['q2']['actual_selected_leading_pressure_primitive']
    proof = p['inertial_volume_bound_proof']
    gap = read(iv, p['radius_factors']['complex_S1_S2_derivative_gap']).a
    require(p['passed'] and p['source_family'] == inputs['q2']['source_family']
            and p['source_log_basis_names'] == BASIS[:4]
            and proof['S1_S2_norm_tube'] == 'eta/8 on the real radial interval'
            and proof['positive_polynomial_coefficients_are_error_bounds_not_signed_values']
            and gap == data['gap'], 'inertial polynomial or its Cauchy gap changed')
    residual = iv.mpf(0)
    for order in (1, 2):
        for row in p['inertial_volume_norms'][f'S{order}']:
            hp, fp, pp, lp = row['source_exponents']
            bound = read(iv, row['coefficient_interval']).b
            require(hp >= 0 and fp == lp == 0 and 0 <= pp <= order and bound > 0
                    and row['coefficient_interval']['directed'], 'invalid physical inertial radius')
            log_value = iv.ln(bound) + hp*data['log_hb'] + (pp-1)*data['log_p2'] - order*data['log_lambda']
            residual += _positive_log_cap(iv, log_value)
    e0 = (data['rho_max']**2*residual/2).b
    return e0, (e0/gap).b


def _axis_inertia(z, data):
    """Analytic radial-constant part of the actual inertial integrand."""
    one, zero = data['one'], data['zero']
    delta = data['delta']
    z_jet = (z, one, zero)
    z2 = z**2
    d = (1-z2, -2*z, -2*one)
    L = (1-delta*z2, -2*delta*z, -2*delta)
    u = (4*z+data['j'], 4*one, zero)
    qa = quotient(tuple((1+delta)*v-4*w for v,w in zip(product(z_jet,u), d, strict=True)), L)
    qz = (qa[1], qa[2], zero)
    # Truncation here only discards unused order two in the inertial density.
    time = quotient(tuple(a+(1-delta)*b/2 for a,b in zip(qa,product(z_jet,qz),strict=True)), L)
    B = tuple(-2*a+b for a,b in zip(product(z_jet,qa),product(d,qz),strict=True))
    advective = quotient(product(u,B), L)
    sq = product(qa, qa)
    # rho^2 is real, Z independent and lies between 16 and the accepted cap.
    # Its radial average carries the same positive enclosure.
    rho_avg = data['rho_square_average']
    return tuple(rho_avg*data['inv_p2']*(a+b/2+c)/2
                 for a,b,c in zip(time[:2], sq[:2], advective[:2], strict=True))


def _cutoff_radii(iv, data, inputs):
    domain = inputs['domain']
    q2 = inputs['q2']
    check = inputs['domain_check']
    convergence = q2['actual_first_raw_n1_function']
    ports = inputs['regular']['canonical_function_bridge']['ports']
    require(check['all_passed'] and check['common_tube_nozero_guards_passed']
            and check['accepted_core_matrix_norm_containment'] and convergence['passed']
            and convergence['all_word_binding']['raw_all_word_convergence_passed']
            and convergence['actual_function_handles']['first']['F1'][:2] == ports['raw_F1_source_nodes']
            and convergence['scale_conjugacy']['canonical_state'][0] == 'F1'
            and domain['source_inner_common_holomorphic_bounds_constructed']
            and domain['nonadjacent_B1_support_inherited'],
            'no common actual all-word F1 bound')
    tail = domain['regular_neumann_tail_after_degree64']
    require(tail['retained_iteration_indices'] == [0,64]
            and tail['source_measure_is_integrated_before_simplex_bound']
            and tail['formula'] == 'M0*(2*a*C0)^k/k!*max(1,alpha*ceil(k/2))^ceil(k/2)',
            'undefined raw regular series majorant')
    M0 = read(iv, domain['initial_Gg_common_tube_bound']).b
    C0 = read(iv, domain['regular_B0_infinity_norm_upper']).b
    # The tail receipt uses directed alpha and an independent certified omitted sum.
    alpha = read(iv, tail['alpha']).b
    a = iv.sqrt(data['rho_max']).b
    omitted = read(iv, tail['omitted_sum_bound']).b
    head = iv.mpf(0)
    for k in range(65):
        m = (k+1)//2
        head += M0*(2*a*C0)**k*max(data['one'],(alpha*m).b)**m/factorial(k)
    raw_f1 = (head+omitted).b
    F0 = read(iv, domain['original_axis_amplitude']['complex_modulus_upper']).b
    Phi = read(iv, domain['actual_collar_radial_bounds']['Phi'][0]['0']['coefficient']).b
    gap = read(iv, domain['common_tubes']['solution_Z_radius']).a/2
    require(min(F0.a, Phi.a, raw_f1.a, gap.a) > 0
            and domain['original_axis_amplitude']['no_source_amplitude_materialized'],
            'nonpositive analytic bound on the cutoff factors')
    # |eta_cut|<=1 and its support has length at most 3S/8. Put the
    # leading and all-word factors on the same complex tube before Cauchy.
    log_product = iv.ln(F0)+iv.ln(Phi)+iv.ln(raw_f1)-data['log_p2']
    cap = _positive_log_cap(iv, log_product)
    e0 = (3*data['rho_max']**2*cap/4).b
    return (e0, (e0/gap).b), {'raw_F1_norm_upper': raw_f1, 'Z_Cauchy_radius': gap, 'leading_F0_norm_upper': F0}


def evaluate(**inputs):
    iv, mp = MPIntervalContext(), MPContext()
    iv.prec = 1800
    data = _prepare(iv, inputs)
    data['rho_square_average'] = hull(iv, iv.mpf(16), (16/(1-2*iv.mpf(2)**(-1024)*data['S'])).b)

    def encode(row):
        return dict(lower_exact_mpf_tuple=list(row._mpi_[0]), upper_exact_mpf_tuple=list(row._mpi_[1]),
                    display=[mp.nstr(mp.make_mpf(bits), 16) for bits in row._mpi_], sign=sign(row))

    q2 = inputs['q2']['actual_correlated_q2_inlet']
    regular = inputs['regular']
    require(regular['passed'] and regular['original_measure'] == 'epsilon^2*hb*rho^2 ds'
            and regular['cutoff_support_q'] == '[0,3/4]'
            and regular['canonical_function_bridge']['ports']['row4_regular_plus_cut'] == [109,111],
            'actual regular source receipt changed')
    visc = _viscous_radii(iv, data, q2, inputs['point'])
    inertial = _inertial_radii(iv, data, inputs)
    cutoff, cut_meta = _cutoff_radii(iv, data, inputs)
    errors = dict(regular_viscous=visc, inertial_nonlinear=inertial, cutoff_swirl=cutoff)
    c0c1_error = tuple(hull(iv, -e, e) for e in (sum(t[k] for t in errors.values()).b for k in range(2)))

    def cell(z):
        inlet = inlet_second_jet(z, data)
        V = _viscous_model(z, inlet, data['delta'])
        A = _axis_inertia(z, data)
        value = tuple(v+a+r for v,a,r in zip(V,A,c0c1_error,strict=True))
        return value, dict(axial_cell=encode(z), inlet_ordinary_C0_C1_C2=[encode(v) for v in inlet],
            viscous_inlet_ordinary_C0_C1=[encode(v) for v in V],
            signed_axis_inertial_ordinary_C0_C1=[encode(v) for v in A],
            row_4_ordinary_C0_C1=[encode(v) for v in value])

    cells = [cell(hull(iv,iv.mpf(k)/64,iv.mpf(k+1)/64))[1] for k in range(-64,64)]
    samples = {}
    exact_samples = {}
    for label in ('-0.3','0','0.3','0.371'):
        exact_samples[label], samples[label] = cell(iv.mpf(label))
    previous = regular['selected_point']['ordinary_Z']
    point_agreement = []
    for k, new in enumerate(exact_samples['0.371']):
        old = read(iv, previous[str(k)]['interval'])
        point_agreement.append(dict(ordinary_Z_order=k, difference_between_enclosures=encode(new-old),
                                    overlap=bool(new.a <= old.b and old.a <= new.b)))
    compare = _compare_forward(iv, encode, data, inputs, exact_samples)
    return dict(passed=True, normalized_unit=UNIT, admitted_axial_domain=['-1','1'],
        ordinary_Z_orders=[0,1], cell_count=len(cells), cells=cells, samples=samples,
        directed_radii={key:[encode(v) for v in row] for key,row in errors.items()},
        cutoff_norm_data={key:encode(value) for key,value in cut_meta.items()},
        independent_saved_selected_point_diagnostic=point_agreement, legacy_fourth_row_diagnostics=compare,
        legacy_212_terms_preserved=True, inherited_error_packets_preserved=True,
        legacy_finite_coefficients_are_function_values=False, functional_equality_from_overlap=False,
        whole_six_forward_source_identity=False, independent_terminal_matching_closed=False,
        ancestor_calculations_replayed=False)
