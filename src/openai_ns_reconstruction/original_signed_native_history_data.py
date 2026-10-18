"""New signed post-restoration histories from immutable saved source data.

Only new endpoint memory transport and signed A/b/c are evaluated here.
Absolute caps enter error radii, never signed values. Physical radii and
exponentially small memories remain in logarithmic form until a proved
finite enclosure is available. No historical numerical provider is imported.
"""
from __future__ import annotations

import hashlib
import json
import math

from mpmath import make_mpf, nstr
from mpmath.libmp import mpf_cmp


NATIVE_SHA = '10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a'
NATIVE_BYTES = 42177499


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(c, row):
    require('log_absolute_upper_bound' not in row, 'A cap cannot be a signed value')
    if 'value' in row:
        row = row['value']
    require('log_absolute_upper_bound' not in row, 'A cap cannot be a signed value')
    return c.mpf([make_mpf(tuple(row[k + '_exact_mpf_tuple'])) for k in ('lower', 'upper')])


def record(value):
    if value is None:
        return None
    return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),
                upper_exact_mpf_tuple=list(value._mpi_[1]),
                lower_display=nstr(make_mpf(value._mpi_[0]), 26),
                upper_display=nstr(make_mpf(value._mpi_[1]), 26), directed=True)


def upper(c, value):
    v = make_mpf(value._mpi_[1])
    return c.mpf([v, v])


def abslog(c, value):
    ends = [(0, *t[1:]) for t in value._mpi_]
    t = ends[0] if mpf_cmp(ends[0], ends[1]) >= 0 else ends[1]
    return None if t[1] == 0 else c.ln(c.mpf([make_mpf(t), make_mpf(t)]))


def logsum(c, terms):
    terms = [t for t in terms if t is not None]
    if not terms:
        return None
    best = max(terms, key=lambda t: make_mpf(t._mpi_[1]))
    return upper(c, upper(c, best) + c.ln(len(terms)))


def logscale(value, scale):
    return None if value is None else value+scale


def logproduct(c, left, right, order):
    return [logsum(c, [None if left[j] is None or right[k-j] is None else
                      left[j] + right[k-j] + c.ln(math.comb(k, j))
                      for j in range(k+1)]) for k in range(order+1)]


def caps(c, rows, order=2):
    result = []
    for k, row in enumerate(rows[:order+1]):
        require(row['ordinary_Z_order'] == k, 'Ordinary derivative order changed')
        result.append(None if row['exact_zero'] else upper(c, read(c, row['log_absolute_upper_bound'])))
    require(len(result) == order+1, 'Incomplete ordinary derivative rows')
    return result


def add(*jets):
    return [sum(row) for row in zip(*jets)]


def scale(jet, factor):
    return [v*factor for v in jet]


def mul(left, right):
    return [sum(math.comb(k, j)*left[j]*right[k-j] for j in range(k+1))
            for k in range(min(len(left), len(right)))]


def signed_rows(c, rows):
    require([r['ordinary_Z_order'] for r in rows[:3]] == [0, 1, 2], 'Signed row order changed')
    return [read(c, r) for r in rows[:3]]


def hull(c, left, right):
    lo = left._mpi_[0] if mpf_cmp(left._mpi_[0], right._mpi_[0]) <= 0 else right._mpi_[0]
    hi = left._mpi_[1] if mpf_cmp(left._mpi_[1], right._mpi_[1]) >= 0 else right._mpi_[1]
    return c.mpf([make_mpf(lo), make_mpf(hi)])


def rows(jet):
    return [dict(ordinary_Z_order=k, value=record(v),
                 kind='signed_function_enclosure_not_selected_value') for k, v in enumerate(jet)]


def small_error(c, logs, label):
    threshold = -1024*c.ln(2)
    for value in logs:
        require(value is None or mpf_cmp(value._mpi_[1], threshold._mpi_[0]) < 0,
                label + ' is not below the chosen finite enclosure radius')
    radius = c.mpf(2)**-1024
    return [c.mpf(0) if value is None else c.mpf([-radius.b, radius.b]) for value in logs]


def compute(c, native_raw, gamma2_record):
    require(len(native_raw) == NATIVE_BYTES and hashlib.sha256(native_raw).hexdigest() == NATIVE_SHA,
            'Native source packet changed before decoding')
    data = json.loads(native_raw)
    raw = data['actual_native_directed_controls_v2']['raw_interval_tree']
    shared, long = raw['shared'], raw['long_reference']
    restoration = raw['restoration']['charts']['restoration']
    z, zero = read(c, shared['z']), c.mpf(0)
    Z = [z, c.mpf(1), zero]
    T = read(c, long['slope']['T'])
    margin = read(c, long['slope']['margin'])
    require(margin.b < 0 and long['slope']['ratio_bound'] == '|B|/T<=1/200',
            'Source long-slope premise is unavailable')
    require(data['actual_long_reference_native_controls']['slope_proof']['sigma_prime_upper'] == 8,
            'Source sigma derivative premise changed')
    Lref = read(c, long['charts']['reference']['scale']['reference_length'])
    require(T.a > 0 and Lref.a > 0, 'Positive source lengths required')
    # Rz lies in the restoration radius domain. Using its upper endpoint
    # here enlarges an error bound, without selecting any source value.
    logRz_upper = upper(c, read(c, restoration['scales']['log_R_upper']))
    phi110 = signed_rows(c, data['actual_switch_post_background_I1_numeric']
                           ['charts']['post_switch']['native_fields']['phi'])
    require(phi110[0].a > 0, 'Positive actual R110 phi required')
    logphi110_squared = 2*c.ln(upper(c, phi110[0]))
    relative = caps(c, long['charts']['long_reshape']['phi_relative'])
    require(relative[0]._mpi_ == c.mpf(0)._mpi_, 'The relative phi value bound must be one')
    density = logproduct(c, relative, relative, 2)
    incoming = data['actual_long_reference_native_controls']['charts']['long_reshape']['memory_bound']['incoming']
    nu = c.mpf(18)/25
    endpoint_logs = {}
    endpoint_parts = {}
    for name, rate in (('B', 2), ('C', 1)):
        old = caps(c, incoming[name])
        inlet = [logscale(value, -rate*T) for value in old]
        forcing = [logscale(value, logphi110_squared-nu*T-c.ln(rate-nu)) for value in density]
        endpoint_logs[name] = [logsum(c, [a, b]) for a, b in zip(inlet, forcing)]
        endpoint_parts[name] = dict(incoming=[record(v) for v in inlet],
                                    density=[record(v) for v in forcing],
                                    total=[record(v) for v in endpoint_logs[name]])
    # The identity invphi^2=2 Rz F0^2/Utheta^2 is differentiated with
    # ordinary binomial factors; all factors here are error majorants.
    inverse_utheta = caps(c, restoration['inverse_Utheta'])
    inverse_utheta_squared = logproduct(c, inverse_utheta, inverse_utheta, 2)
    f2 = [upper(c, read(c, v)) for v in shared['F2_logs5'][:3]]
    inverse_phi = [v+c.ln(2)+logRz_upper
                   for v in logproduct(c, f2, inverse_utheta_squared, 2)]
    deviations, components = {}, {}
    for name, rate, canonical in (('B', 2, c.mpf(5)/6), ('C', 1, c.mpf(5))):
        memory = [logscale(v, -rate*Lref) for v in logproduct(c, endpoint_logs[name], inverse_phi, 2)]
        subtract = c.ln(canonical)-(c.mpf(rate)-c.mpf(4)/5)*Lref
        bound = [logsum(c, [memory[0], subtract]), memory[1], memory[2]]
        deviations[name] = small_error(c, bound, name+' normalized Rz deviation')
        components[name] = dict(incoming_Rsh=[record(v) for v in memory],
                                canonical_inlet_subtraction=record(subtract),
                                deviation_log_caps=[record(v) for v in bound],
                                signed_terms_retained=True)
    inlet_source = data['actual_restoration_background_I1_numeric']['charts']['restoration']['native_source']
    U110, E = signed_rows(c, inlet_source['U_R110']), signed_rows(c, inlet_source['E'])
    initialA = caps(c, incoming['A'])
    Usquared_logs = [abslog(c, value) for value in mul(U110, U110)]
    total_length = T+Lref
    Alogs = [logscale(logsum(c, [a, b]), -total_length) for a, b in zip(initialA, Usquared_logs)]
    epsilon_A = small_error(c, Alogs, 'A_Rz incoming deviation')
    gamma1 = read(c, data['actual_restoration_background_I1_numeric']['charts']['restoration']
                         ['arithmetic']['source']['gamma_proof']['endpoint'])
    gamma2 = read(c, gamma2_record)
    baseA = scale(mul(Z, Z), 16)
    memoryA = add(scale(mul(Z, E), 8*gamma1), scale(mul(E, E), gamma2))
    results = {}
    post_fields = data['actual_restoration_post_jet_repair']['corrected_background']['charts']['post_restore']
    operations = data['actual_long_reference_background_I1_numeric']['scalar_graph']['node_operations']
    require(operations['reference_length'] == 'node 81519 = logref - T - 8',
            'Fixed source total-length identity changed')
    require(post_fields['native_source']['formula'] == 'U=4Z; M=4Z+(M_restore-4Z)*exp(-q)',
            'Saved ordinary-Z M rows need the same affine exponential source')
    for chart, p in (('post_restore', c.mpf([0, 1])), ('Rm', c.mpf(1))):
        A = add(baseA, scale(memoryA, c.exp(-p)), scale(epsilon_A, c.exp(-(p+1))))
        b = add([c.mpf(5)/6, zero, zero], deviations['B'])
        cc = add([c.mpf(5), zero, zero], deviations['C'])
        M = signed_rows(c, post_fields['endpoint_fields']['M_end'])
        if chart == 'post_restore':
            # Every ordinary-Z row is affine in exp(-p). Thus the union of
            # the saved endpoint enclosures covers p in [0,1]. The old M
            # provider is neither called nor numerically reproduced.
            start = signed_rows(c, post_fields['endpoint_fields']['M_start'])
            M = [hull(c, a, b) for a, b in zip(start, M)]
        results[chart] = dict(A=rows(A), b=rows(b), c=rows(cc), U=rows(scale(Z, 4)), M=rows(M),
                             coordinate=record(p), canonical_centres_selected_as_actual=False)
    return dict(passed=True, source_report_sha256=NATIVE_SHA, ordinary_Z_orders=[0, 1, 2],
                signed_histories=results, new_long_endpoint_memory_caps=endpoint_parts,
                normalized_Rz_memory_transport=components, inverse_phi_Rz_squared_log_caps=[record(v) for v in inverse_phi],
                A_Rz_deviation_log_caps=[record(v) for v in Alogs],
                A_Rz_identity='U110^2+exp(-(T+Lref))*(A_R110-U110^2)',
                A_R110_unknown_signed_memory_retained=True,
                finite_enclosure_radius=record(c.mpf(2)**-1024), gamma1_saved=record(gamma1),
                gamma2_new=record(gamma2), incoming_A_B_C_retained=True,
                long_density_decay_nu='18/25', long_density_proof='phi_y/phi<=-2/5+8/200=-9/25',
                long_endpoint_formula='exp(-rate*T)*inlet+phi110_upper^2*Dk*exp(-18*T/25)/(rate-18/25)',
                normalized_memory_formula='exp(-rate*Lref)*J_Rsh/phi_Rz^2-canonical*exp(-(rate-4/5)*Lref)',
                A_post_formula='16Z^2+exp(-p)*(8Z E gamma1+E^2 gamma2)+exp(-(p+1))*epsilon_A_Rz',
                derivative_of_error_bound_used=False, caps_promoted_to_values=False,
                old_gamma1_recomputed=False, physical_radius_materialized=False,
                accepted_mathematical_producers_called=False)


def compute_p2(c, native_raw, history_packet):
    """Bind the new signed histories and saved independent pressure once.

    K=p2/(R*Utheta/2) is a finite signed function enclosure. Its derivative
    K_Z differs from the physical p2_Z coefficient K_Z+(log Utheta)_Z*K.
    Both are retained. The enormous positive R scale is never exponentiated.
    """
    from .original_signed_native_p2_kernel import evaluate_scaled_signed_p2_source

    require(len(native_raw) == NATIVE_BYTES and hashlib.sha256(native_raw).hexdigest() == NATIVE_SHA,
            'Native source packet changed before pressure decoding')
    data = json.loads(native_raw)
    raw = data['actual_native_directed_controls_v2']['raw_interval_tree']
    shared = raw['shared']
    z, delta, zero, one = read(c, shared['z']), read(c, shared['delta']), c.mpf(0), c.mpf(1)
    Z, d = [z, one, zero], [1-z*z, -2*z, -2*one]
    q = 1+z*z
    q_squared = [q*q, 4*z*q, 4+12*z*z]
    g = [-4*z/q, -4*(1-z*z)/(q*q), 8*z*(3-z*z)/(q*q*q)]
    logP = read(c, raw['long_reference']['charts']['reference']['scale']['logP'])
    logPstar = read(c, shared['logPstar'])
    pressure_source = data['actual_micro_b_p2_mixed_bounds']['whole_tube_pressure']
    require(pressure_source['physical_scale']['Pstar_squared_not_applied'] is True,
            'Independent pressure normalization changed')
    P0 = [read(c, row) for row in pressure_source['p0_jets'][:3]]
    results = {}
    for chart, histories in history_packet['signed_histories'].items():
        U, M, A, b, cc = [signed_rows(c, histories[name]) for name in ('U', 'M', 'A', 'b', 'c')]
        Mz, Az = M[1:], A[1:]
        hydro = add(scale(U, -1), scale(add(M, scale(mul(Z, Mz), -1)), (1-delta)/2),
                    mul(U, add(scale(mul(Z, M), 1-delta), mul(d, Mz))),
                    scale(mul(Z, A), 2*delta), scale(mul(d, Az), -1))
        require(len(hydro) == 2, 'Only source-supported hydro Z0/Z1 may be supplied')
        p = read(c, histories['coordinate'])
        offset = -7+p
        inverse_u2 = scale(q_squared, c.exp(-2*logP-offset/5))
        # Distinct logP and logPstar records are retained, rather than
        # silently identifying their numerical enclosures by overlap.
        ratio = scale(q_squared, c.exp(2*(logPstar-logP)-offset/5))
        result = evaluate_scaled_signed_p2_source(context=c, Z=Z, delta=delta,
                    hydro=hydro, P0=P0, b=b, c=cc, inverse_Utheta_squared=inverse_u2,
                    pressure_ratio=ratio, g=g, pressure_is_physical=False)
        K = list(result.K)
        loggradient = -2*z/q
        physical_Z1 = K[1]+loggradient*K[0]
        require(K[0].b < 0 or K[0].a > 0, 'Signed native p2 value is unresolved')
        logR = read(c, data['actual_restoration_post_jet_repair']['corrected_background']
                         ['charts']['post_restore']['log_R_upper'])
        # Source offsets differ by p; the saved scale is an upper radius
        # bound for the full chart. It is explicitly scale metadata here.
        log_scale_upper = logR+logP+offset/10-c.ln(q)-c.ln(2)
        results[chart] = dict(K=rows(K), p2_Z1_over_positive_scale=record(physical_Z1),
            physical_positive_scale='R*Utheta/2', scale_source='R=110*exp(logref-7+p); Utheta=exp(logP+(-7+p)/10)/(1+Z^2)',
            log_positive_scale_upper=record(log_scale_upper),
            p2_value_sign='negative' if K[0].b < 0 else 'positive',
            p2_Z_over_p2=record(physical_Z1/K[0]),
            ordinary_Z_derivatives_of_scaled_K=True,
            hydro=rows(hydro), normalized_swirl=rows(list(result.S)[:2]),
            scaled_hydro=rows(list(result.hydro_scaled)[:2]),
            scaled_pressure=rows(list(result.pressure_scaled)[:2]),
            independent_p0=rows(P0), pressure_ratio=rows(ratio),
            inverse_Utheta_squared=rows(inverse_u2), g=rows(g),
            source_coordinate=histories['coordinate'], physical_scale_materialized=False)
    return dict(passed=True, charts=results, ordinary_Z_orders=[0, 1],
                source_report_sha256=NATIVE_SHA, physical_R_applied_once=True,
                normalized_P0_used_with_explicit_Pstar_squared_ratio=True,
                F0_derivative_cancelled_with_same_source_phi=True,
                canonical_history_centres_selected=False, incoming_memory_retained=True,
                midpoint_or_logcap_promoted=False, full_native_p2_signed_scaled_enclosure=True,
                full_NS_residual_validated=False)
