"""Same-original-source Fourier C3 bounds consumed by pre-Rp pressure C1.

The signed functions, selected N and physical scales are unchanged. Only
their directed derivative enclosures are improved. The resulting nonzero
pre-Rp errors enter the accepted Rp-to-Rv derivative models, giving actual
r_minus-to-infinity pressure derivative enclosures, not repair coefficients.
"""
import copy
import gzip
import importlib
import json
import math
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_C3_target_ranges as rg
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as partial
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source as previous
import lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1 as endpoint
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1 as pulse_C1
import lei_ren_part1_paper_compliant_current_original_C3_limit_controls as control_C3
import lei_ren_part1_paper_compliant_current_original_C3_repair_band as repair_C3

HERE, PREFIX, sha, ep = previous.HERE, previous.PREFIX, previous.sha, previous.ep
source_precision = previous.source_precision
NAME = PREFIX + 'current_original_n1_preRp_pressure_C1.json.gz'
RECEIPT = PREFIX + 'current_original_n1_preRp_pressure_C1_check.json'
GATE = 'current_original_n1_preRp_pressure_C1_installed'
UPPER_GATE = 'current_original_n1_rminus_to_infinity_pressure_C1_accuracy_installed'
OPEN = tuple(dict.fromkeys((*endpoint.OPEN, *pulse_C1.OPEN,
    'complete_signed_defect_C1_functions_bound', 'actual_n1_repair_coefficients_computed')))
CELLS = (('-1', '-.5'), ('-.5', '0'), ('0', '.5'), ('.5', '1'))
interval = previous.interval


def logcap(c,value):
    if isinstance(value,rg.LogUpper):
        return rg.LogUpper(c) if value.log is None else rg.LogUpper(c,interval(c,value.log))
    return previous.log_cap(c,value)


def record(value):
    if isinstance(value,rg.JetBound):
        return record(rg.rows(value))
    if isinstance(value,dict):
        return {key:record(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):
        return [record(item) for item in value]
    if hasattr(value,'_mpf_'):
        return dict(lower=mp.nstr(value,90),upper=mp.nstr(value,90),width='0.0',
            lower_exact_mpf_tuple=list(value._mpf_),upper_exact_mpf_tuple=list(value._mpf_))
    return previous.record(value)


def jet(c, records):
    if len(records) != 4:
        raise ValueError('Four ordinary Z derivative rows required')
    return rg.JetBound(*(logcap(c, value) for value in records))


def fourier_theorem():
    return dict(
        chi='2*h*sum(r^n*cos((n+1)*psi), n>=0)',
        h='(1+u^2)^(-1/2)', r='u*h', hyperbolic_coordinate='u=sinh(alpha)',
        exact_covariance='integral chi(u)*chi(v)/(4*pi)=sech(alpha-beta)',
        squared_hyperbolic_derivative_norms=[1, 1, 5, 61],
        L2_u_derivatives_over_2sqrtpi=[1, 1, 4, 17],
        point_bounds=dict(chi='2+4*abs(u)', chi_u=12, chi_uu=80,
            chi_psi='8*(1+abs(u))^3', chi_psiu='40*(1+abs(u))^2',
            chi_psipsi='32*(1+abs(u))^5'),
        unnormalized_inverse='F=psi+integral_0^psi t^2-2*pi*phi*nu=0',
        t='t0+q*chi', nu='1+t0^2+2*q^2',
        inverse_monotonicity='F_psi=1+t^2>=1',
        mixed_partial_division='2*abs(t)/(1+t^2)<=1',
        phase_domain=[0, 1], inverse_angle_domain=['0', '2*pi'],
        genuine_fixed_angle_derivatives_before_composition=True,
        same_original_q_root_and_sigma_C3_source=True,
        no_derivative_of_a_magnitude_cap=True,
        same_signed_A_and_B_function_bounds_only=True)


def sharpen_primitive(c, cell):
    """Bound the original B derivatives using Parseval and implicit jets."""
    bd = rg.Bounds(c)
    proof = cell['proof']
    if not (proof['full_original_phase_box'] == [0, 1]
            and proof['live_source_phase_Z_exact_zero']
            and proof['exact_same_source_P0']):
        raise ValueError('Original full phase and independent P0 binding required')
    old = jet(c, cell['B'])
    if proof['branch'] == 'exact_original_flat':
        if any(v.log is not None for v in rg.rows(old)):
            raise ValueError('Accepted original flat B must be exact zero')
        return dict(selected_B=rg.rows(old), candidate_B=rg.rows(old), selected_A=rg.rows(jet(c,cell['A'])),
                    branch='exact_original_flat', original_function_handles=proof['actual_C3_function_handles'])
    raw = proof['actual_source_ordinary_Z3_rows']
    jets = {key: jet(c, values) for key, values in raw.items()}
    q = jet(c, proof['q_bounds'])
    a, b, E, t0, p2 = [jets[key] for key in ('a', 'b', 'E', 't0', 'p2')]
    u = bd.scaled(bd.product(q, p2), rg.LogUpper(c, -interval(c, proof['actual_dstar_log'])))
    u0, u1, u2, u3 = rg.rows(u)
    q0, q1, q2, q3 = rg.rows(q)
    t00, t01, t02, t03 = rg.rows(t0)
    const, plus, power = bd.constant, bd.sum, bd.power
    sp, s2 = const(2*c.sqrt(c.pi)), const(c.sqrt(2*c.pi))
    W = [bd.one, u1, plus(bd.scale(power(u1, 2), 4), u2),
         plus(bd.scale(power(u1, 3), 17), bd.scale(u1*u2, 12), u3)]
    T = [plus(s2*rg.rows(t0)[j], sp*plus(*(
        bd.scale(rg.rows(q)[k]*W[j-k], math.comb(j, k)) for k in range(j+1))))
        for j in range(4)]
    nu = bd.add(bd.fixed(bd.one), bd.product(t0, t0), bd.scaled(bd.product(q, q), const(2)))
    D1 = plus(bd.scale(T[0]*T[1], 2), const(2*c.pi)*nu.Z)
    D2 = plus(bd.scale(plus(power(T[1], 2), T[0]*T[2]), 2), const(2*c.pi)*nu.ZZ)
    D3 = plus(bd.scale(plus(bd.scale(T[1]*T[2], 3), T[0]*T[3]), 2), const(2*c.pi)*nu.ZZZ)
    peak = plus(const(2), bd.scale(u0, 4))
    TZ = plus(t01, q1*peak, bd.scale(q0*u1, 12))
    TZZ = plus(t02, q2*peak, bd.scale(q1*u1, 24),
               q0*plus(bd.scale(power(u1, 2), 80), bd.scale(u2, 12)))
    Tp = bd.scale(q0*power(plus(bd.one, u0), 3), 8)
    TpZ = plus(bd.scale(q1*power(plus(bd.one, u0), 3), 8),
               bd.scale(q0*u1*power(plus(bd.one, u0), 2), 40))
    Tpp = bd.scale(q0*power(plus(bd.one, u0), 5), 32)
    P1 = D1
    P2 = plus(D2, bd.scale(TZ*P1, 2), Tp*power(P1, 2))
    P3 = plus(D3, bd.scale(plus(bd.scale(power(TZ, 2), 2), TZZ)*P1, 3),
        bd.scale(plus(bd.scale(TZ*Tp, 2), TpZ)*power(P1, 2), 3),
        plus(bd.scale(power(Tp, 2), 2), Tpp)*power(P1, 3),
        bd.scale(plus(TZ, Tp*P1)*P2, 3))
    Tpoint0 = plus(t00, q0*peak)
    C = rg.JetBound(s2*T[0], plus(s2*T[1], Tpoint0*P1),
        plus(s2*T[2], bd.scale(TZ*P1, 2), Tp*power(P1, 2), Tpoint0*P2),
        plus(s2*T[3], bd.scale(TZZ*P1, 3), bd.scale(TpZ*power(P1, 2), 3),
             Tpp*power(P1, 3), bd.scale(TZ*P2, 3), bd.scale(Tp*P1*P2, 3), Tpoint0*P3))
    active_a = rg.JetBound(bd.minimum(a.value, const('2.5')), a.Z, a.ZZ, a.ZZZ)
    M = bd.add(bd.scaled(bd.product(active_a, C), const(1/(2*c.pi))), b)
    candidate = bd.scaled(bd.product(E, M), const('.5'))
    selected = rg.JetBound(*(bd.minimum(x,y) for x,y in zip(rg.rows(old),rg.rows(candidate))))
    chi = rg.JetBound(bd.one, P1*const(1/(2*c.pi)), P2*const(1/(2*c.pi)), P3*const(1/(2*c.pi)))
    candidate_A = bd.scaled(bd.product(active_a, chi), const('.5'))
    selected_A = rg.JetBound(*(bd.minimum(x,y) for x,y in zip(rg.rows(jet(c,cell['A'])),rg.rows(candidate_A))))
    return dict(selected_B=rg.rows(selected), candidate_B=rg.rows(candidate),
        selected_A=rg.rows(selected_A), candidate_A=rg.rows(candidate_A),
        original_function_handles=proof['actual_C3_function_handles'], branch=proof['branch'],
        L2_t_Z3=T, inverse_Z3=[P1, P2, P3], composed_T1_Z3=rg.rows(C),
        original_signed_C0_through_C3_handles_unchanged=True)


def updated_density(c, cell, sharp):
    """Original density algebra with the accepted epsilon-domain constants."""
    bd = rg.Bounds(c)
    original = cell['actual_source_and_primitive_range']
    A, B = rg.JetBound(*sharp['selected_A']), rg.JetBound(*sharp['selected_B'])
    E, V = jet(c, original['E']), jet(c, original['V'])
    old = cell['actual_density_ranges']
    factor = logcap(c, old['full_exp_exprel_factor_upper'])
    epsilon = logcap(c, old['epsilon_upper'])
    F = rg.JetBound(E.value*A.value*factor,
        bd.sum(E.Z*A.value,E.value*A.Z)*factor,
        bd.sum(E.ZZ*A.value,bd.scale(E.Z*A.Z,2),E.value*bd.sum(A.ZZ,bd.power(A.Z,2)*epsilon))*factor,
        bd.sum(E.ZZZ*A.value,bd.scale(E.ZZ*A.Z,3),
            bd.scale(E.Z*bd.sum(A.ZZ,bd.power(A.Z,2)*epsilon),3),
            E.value*bd.sum(A.ZZZ,bd.scale(A.Z*A.ZZ*epsilon,3),bd.power(A.Z,3)*bd.power(epsilon,2)))*factor)
    EF, VF, EB, VB = [bd.product(x,y) for x,y in ((E,F),(V,F),(E,B),(V,B))]
    F2, B2 = bd.product(F,F), bd.product(B,B)
    first = dict(m=B,h=F,k=bd.add(VF,EB),e=bd.add(bd.scaled(VB,bd.constant(2)),EF),p=EF)
    second = dict(m=bd.fixed(bd.zero),h=bd.fixed(bd.zero),k=bd.product(F,B),
        e=bd.add(B2,bd.scaled(F2,bd.constant('.5'))),p=bd.scaled(F2,bd.constant('.5')))
    return {key:bd.add(first[key],bd.scaled(second[key],epsilon)) for key in first}


def updated_targets_and_compact(c, ranges, sources, old_control, old_band):
    """Propagate to the SAME target functions and SAME implicit control branch."""
    bd = rg.Bounds(c)
    result = []
    new_histories = {}
    maximum = bd.zero
    for transport in ranges['actual_four_Z_C3_target_transports']:
        ends = tuple(transport['exact_Z_cell'])
        history = {key:bd.fixed(bd.zero) for key in ('m','h','k','e','p')}
        for chart in transport['actual_17_chart_C3_transports']:
            local = {key:bd.fixed(bd.zero) for key in history}
            for cell in chart['cells']:
                sharp = sharpen_primitive(c,cell['actual_source_and_primitive_range'])
                density = updated_density(c,cell,sharp)
                for key in history:
                    w = cell['original_positive_own_rate_weights'][key]
                    measure = logcap(c,w['original_mass'])*logcap(c,w['suffix_decay'])
                    local[key] = bd.add(local[key],bd.scaled(density[key],measure))
            history = {key:bd.add(bd.scaled(history[key],logcap(c,chart['original_own_rate_memory'][key])),local[key])
                       for key in history}
        new_histories[ends] = history
        A = jet(c,transport['actual_terminal_amplitude_C3'])
        A2 = bd.product(A,A)
        logA = interval(c,transport['actual_positive_terminal_amplitude_log_lower'])
        logmu = interval(c,transport['actual_positive_mu_log'])
        targets = {}
        for row,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
            targets[row] = bd.quotient(history[key],A if degree==1 else A2,degree*logA)
        joint = bd.add(history['k'],bd.product(A,history['m']))
        angular_key = next(key for key in transport['actual_N_scaled_five_target_C3_ranges'] if key not in targets)
        targets[angular_key] = bd.quotient(joint,bd.scaled(A2,rg.LogUpper(c,logmu)),2*logA+logmu)
        cell_max = previous.maximum(c,[value.ZZZ for value in targets.values()])
        maximum = previous.maximum(c,[maximum,cell_max])
        result.append(dict(ends=ends,histories=history,targets=targets,maximum_ZZZ=cell_max,
            same_signed_target_functions=transport['actual_five_target_C3_function_handles'],
            original_mu_amplitude_and_rate_memories_retained=True))
    bounds = old_control['actual_quantitative_C3_bounds']
    rd = lambda key:logcap(c,bounds[key])
    CQ,rho,H2,epsilon,inverse = [rd(k) for k in ('original_quadratic_upper','original_C1_ball_radius_upper',
        'actual_limit_and_iterate_ZZ_upper','actual_inverse_epsilon_upper','actual_J_inverse_upper')]
    H3 = bd.minimum(rd('actual_limit_and_iterate_ZZZ_upper'),inverse*bd.sum(maximum,bd.scale(CQ*rho*H2*epsilon,6)))
    H = rg.JetBound(rho,rho,H2,H3)
    inputs = importlib.import_module(PREFIX+'current_original_n1_preRp_C1_inputs')
    compact = {}
    for old in old_band['actual_four_Z_C3_band_ranges']:
        ends = tuple(old['exact_Z_cell'])
        row = copy.deepcopy(old)
        A = jet(c,old['actual_endpoint_source']['ordinary_terminal_amplitude_C3'])
        for profile,shape in zip(row['actual_selected_N_profile_y0_y1_y2_C3_bounds'],old['original_bump_y0_y1_y2_caps']):
            cap = logcap(c,shape)
            newF = bd.scaled(H,bd.scale(cap*epsilon,3))
            newG = bd.scaled(H,bd.scale(cap*epsilon,2))
            for key,value in (('F',newF),('G',newG),('V',bd.product(A,newG)),('delta_V',bd.product(A,newG))):
                if key in profile:
                    profile[key] = record(rg.JetBound(*(bd.minimum(x,y) for x,y in zip(rg.rows(jet(c,profile[key])),rg.rows(value)))))
        V = jet(c,row['actual_selected_N_profile_y0_y1_y2_C3_bounds'][0]['V'])
        incoming = bd.scaled(new_histories[ends]['m'],epsilon)
        mass = bd.constant(interval(c,old['actual_own_rate_positive_masses']['m']))
        mean = bd.add(incoming,bd.scaled(V,mass))
        oldm = jet(c,row['actual_partial_correction_C3_bounds']['m'])
        row['actual_partial_correction_C3_bounds']['m'] = record(rg.JetBound(*(bd.minimum(x,y) for x,y in zip(rg.rows(oldm),rg.rows(mean)))))
        baseline = sources['baseline']['|'.join(ends)+'|first_micro']
        compact['|'.join(ends)] = record(inputs._compact_c1_cell(c,ends,row,logcap(c,baseline['Pstar_record']),
            c.exp(-4*(c.exp(40)+11)-30),interval(c,sources['quiet_source']['K']['log_coefficient']),
            interval(c,sources['absolute_Rc_source']['log_Rc'])))
    return dict(four_Z_targets=result,global_target_ZZZ=maximum,control_ZZZ=H3,
        old_control_ZZZ=rd('actual_limit_and_iterate_ZZZ_upper'),
        same_original_implicit_J_inverse=inverse,same_selected_N=True,
        derivative_equation='J*h_ZZZ=-d_ZZZ-3*D2Q[h_Z,h_ZZ]/N',
        compact=compact,original_signed_target_and_control_functions_unchanged=True)


def _read(module):
    blob = (HERE/module.NAME).read_bytes()
    return json.loads(gzip.decompress(blob) if module.NAME.endswith('.gz') else blob)


def load_checked(module, hashes, family=None):
    data = _read(module)
    receipt = json.loads((HERE/module.RECEIPT).read_bytes())
    gates = module.GATES if hasattr(module, 'GATES') else (module.GATE,)
    if not receipt.get('all_passed') or not all(receipt.get(key) for key in gates):
        raise ValueError('Checked source required: '+module.NAME)
    if receipt['source_family'] != data['source_family'] or (family is not None and family != data['source_family']):
        raise ValueError('Same original source family required: '+module.NAME)
    pins = dict(receipt['input_hashes'])
    for name in (module.NAME, module.RECEIPT, Path(module.__file__).name):
        pins[name] = sha(name)
    for name, digest in pins.items():
        if name in hashes:
            if hashes[name] != digest:
                raise ValueError('Conflicting immutable source '+name)
        elif sha(name) != digest:
            raise ValueError('Changed immutable source '+name)
        hashes[name] = digest
    return data


def upgraded_chart_data(c, ranges, old_data, source_data, partial_data):
    """Sum new B_ZZZ along the whole rate-one mean history, without resets."""
    bd = rg.Bounds(c)
    old_views = {(tuple(v['ends']), v['chart']): v['view'] for v in old_data['source_data']['actual_chart_views']}
    logs = [interval(c, v['selected_N_log']) for v in partial_data['samples'].values()]
    if not logs or any(value._mpi_ != logs[0]._mpi_ for value in logs):
        raise ValueError('One actual original selected N required')
    logN = logs[0]
    out, native = [], []
    for transport in ranges['actual_four_Z_C3_target_transports']:
        ends = tuple(transport['exact_Z_cell'])
        if ends not in CELLS:
            raise ValueError('Original four-Z domain required')
        charts = transport['actual_17_chart_C3_transports']
        first = charts[0]['actual_N_scaled_incoming_C3']['m'][3]
        if not first.get('exact_zero'):
            raise ValueError('Original initial finite-N mean increment must be zero')
        prefix = bd.zero
        for chart in charts:
            name = chart['chart']
            baseline = source_data['baseline']['|'.join(ends)+'|'+name]
            Pstar = logcap(c, baseline['Pstar_record'])
            original_view = old_views[(ends, name)]
            for cell in chart['cells']:
                if not cell['physical_Jacobian_applied_once']:
                    raise ValueError('Original native cell measure required')
                sharp = sharpen_primitive(c, cell['actual_source_and_primitive_range'])
                # A whole-chart rate-zero width also bounds every native
                # subcell. This bounds all free-query prefixes, including a
                # query inside the current cell, without endpoint kernels.
                mass = logcap(c, original_view['original_positive_dlogR_width'])
                prefix = bd.sum(prefix, sharp['selected_B'][3]*mass)
                native.append(dict(ends=ends, chart=name, exact_left=cell['exact_left'], exact_right=cell['exact_right'],
                    selected_B=sharp['selected_B'], candidate_B=sharp['candidate_B'],
                    original_function_handles=sharp['original_function_handles'],
                    branch=sharp['branch'], original_dlogR_width_upper=mass,
                    actual_whole_chart_rate_zero_width_used_as_subcell_upper=True,
                    full_prefix_B_ZZZ=prefix, rate_one_kernel_upper=1, predecessor_never_reset=True,
                    new_physical_local_B_ZZZ_over_N=sharp['selected_B'][3]*rg.LogUpper(c, Pstar.log-logN)
                    if sharp['selected_B'][3].log is not None else bd.zero))
            U = jet(c, original_view['whole_cell_physical_U_Z3_magnitude'])
            oldm = jet(c, original_view['whole_cell_physical_m_Z3_magnitude'])
            baseline_m3 = logcap(c, baseline['physical_m_Z3_magnitude'][3])
            increment = prefix*rg.LogUpper(c, Pstar.log-logN)
            m3 = bd.minimum(oldm.ZZZ, bd.sum(baseline_m3, increment))
            m = rg.JetBound(oldm.value, oldm.Z, oldm.ZZ, m3)
            delta = c.exp(-4*(c.exp(40)+11)-30)
            fresh_q = physical.physical._physical_q_bounds(c, ends, delta, U, m)
            oldq = original_view['whole_cell_physical_q_Z2_magnitude']
            q = [bd.minimum(v, logcap(c, old)) for v, old in zip(fresh_q, oldq)]
            out.append(dict(ends=ends, chart=name, U=rg.rows(U), m=rg.rows(m), q=q,
                baseline_m_ZZZ=baseline_m3, full_prefix_B_ZZZ=prefix,
                physical_mean_increment_ZZZ=increment, actual_selected_N_log=logN,
                q_recovery_source='physical._physical_q_bounds',
                original_leading_history_retained=True, all_preceding_charts_in_mean=True,
                common_Pstar_applied_once=True, original_B_function_unchanged=True))
    if len(out) != 68 or len(native) != 228:
        raise ValueError('Complete 17-chart/four-Z/57-native-cell cover required')
    return out, native


def active_pressure_C1(c, ends, views, endpoint_view, delta=None):
    """Differentiate the defining reduced-S coefficients and true products."""
    bd = rg.Bounds(c)
    maxima = lambda key, k: previous.maximum(c, [logcap(c, v[key][k]) for v in views])
    U, UZ, UZZ = [maxima('U', k) for k in range(3)]
    q, qZ, qZZ = [maxima('q', k) for k in range(3)]
    product = lambda x, y: (x[0]*y[0], bd.sum(x[1]*y[0], x[0]*y[1]))
    square = product((U, UZ), (U, UZ))
    pairs = dict(q=(q, qZ), qZ=(qZ, qZZ), q2=product((q, qZ), (q, qZ)),
        U2=square, UqZ=product((U, UZ), (qZ, qZZ)), UUZ=product((U, UZ), (UZ, UZZ)))
    z = IntervalTaylor.variable(c, c.mpf(ends), 1)
    delta = c.exp(-4*(c.exp(40)+11)-30) if delta is None else c.mpf(delta)
    if not 0 <= ep(delta)[0] <= ep(delta)[1] < c.mpf('.5'):
        raise ValueError('Original positive L majorant requires 0<=delta<1/2')
    L, d, a = 1-z*z*delta, 1-z*z, 1-delta
    invL = L.reciprocal()
    def term(coefficient, pair, radius_factor=1):
        caps = [bd.constant(c.mpf(ep(abs(coefficient[k]))[1])) for k in range(2)]
        return (bd.scale(caps[0]*pair[0], radius_factor),
                bd.scale(bd.sum(caps[1]*pair[0], caps[0]*pair[1]), radius_factor))
    def sum_terms(items):
        return tuple(bd.sum(*(v[k] for v in items)) for k in range(2))
    half = IntervalTaylor.constant(c, c.mpf('.5'), 1)
    S0 = sum_terms([term(invL, pairs['q'], 2), term(z*a*invL/2, pairs['qZ']),
        term(half, pairs['q2'], 2), term(d*invL, pairs['UqZ']),
        term(-2*z*z*invL*invL, pairs['U2'], 2), term(z*z*(-2*delta)*invL*invL, pairs['U2']),
        term(2*z*d*invL*invL, pairs['UUZ'])])
    S1 = sum_terms([term(invL, pairs['q'], 2), term(-invL, pairs['q']),
        term(z*a*invL/2, pairs['qZ']), term(half, pairs['q2'], 2), term(-half, pairs['q2']),
        term(d*invL, pairs['UqZ']), term(-2*z*z*invL*invL, pairs['U2'], 2),
        term(z*z*(2*a)*invL*invL, pairs['U2']), term(2*z*d*invL*invL, pairs['UUZ'])])
    lower, upper = endpoint_view['lower']['caps'], endpoint_view['upper']['caps']
    visc = dict(m13=bd.sum(*(logcap(c, v['q_plus_q_y_Z']) for v in (lower, upper))),
                m15=bd.sum(*(logcap(c, v['q_y_Z']) for v in (lower, upper))))
    return dict(fields=dict(U=U, U_Z=UZ, U_ZZ=UZZ, q=q, q_Z=qZ, q_ZZ=qZZ),
        differentiated_products=pairs, S0_Z_over_Rc=S0[1], S1_Z_over_Rc_squared=S1[1],
        viscous_Z_over_radius_power=visc,
        m13=bd.sum(bd.scale(S0[1], '.5'), visc['m13']),
        m15=bd.sum(bd.scale(S1[1], '.5'), visc['m15']),
        original_Rc_ge_one_used_only_for_endpoint_bound=True,
        boundary_radius_upper='R<=Rc', integral_dR_upper='Rc', integral_RdR_upper='Rc^2',
        all_coefficient_derivatives_retained=True, genuine_U_ZZ_and_q_ZZ_used=True)


@source_precision
def build_data(c, ranges, old, sources, partial_data, endpoints, old_control, old_band):
    if not old['source_data']['geometry']['original_Rc_ge_one']:
        raise ValueError('Original fixed Rc>=1 geometry required')
    control_logN = interval(c,old_control['actual_quantitative_C3_bounds']['actual_log_N'])
    if any(interval(c,v['selected_N_log'])._mpi_ != control_logN._mpi_ for v in partial_data['samples'].values()):
        raise ValueError('Physical mean and implicit controls must use the same original selected N')
    if old_band['actual_selected_repair_integer'] != old_control['actual_quantitative_C3_bounds']['actual_same_repair_integer']:
        raise ValueError('Compact profiles and implicit branch selected integer differ')
    control_update = updated_targets_and_compact(c,ranges,sources,old_control,old_band)
    effective_sources = copy.deepcopy(sources)
    effective_sources['compact'] = control_update['compact']
    upgraded, native = upgraded_chart_data(c, ranges, old, sources, partial_data)
    cells = []
    for ends in CELLS:
        views = [v for v in upgraded if tuple(v['ends']) == ends]
        endpoint_view = next(v for v in endpoints['samples'].values() if tuple(v['exact_Z_cell']) == ends)
        active = active_pressure_C1(c, ends, views, endpoint_view)
        rows = {}
        for key in ('m13', 'm15'):
            unit = interval(c, sources['scale_ratios'][key])
            compact = logcap(c, effective_sources['compact']['|'.join(ends)]['normalized_first_Z_upper'][key])
            logs = dict(active=active[key].log+unit,
                compact=compact.log+unit if compact.log is not None else None)
            if any(value is not None and ep(value)[1] >= -1000 for value in logs.values()):
                raise ValueError('Original pre-Rp derivative error exceeds budget: '+key+' '+str(ends))
            if not sources['quiet_source']['source_scaled_Z1_exact_zero']:
                raise ValueError('Actual quiet derivative source must be Z independent')
            rows[key] = dict(normalized_log_caps=logs, quiet_Z_exact_zero=True,
                active_nonzero_absolute_upper=c.exp(-1000), compact_nonzero_absolute_upper=c.exp(-1000),
                total_nonzero_absolute_upper=2*c.exp(-1000), counted_once=True,
                original_pressure_units_fixed=True, no_signed_zero_assigned_to_bounded_terms=True)
        cells.append(dict(ends=ends, active=active, fields=rows))
    return dict(source_family=ranges['source_family'], cells=cells, upgraded_chart_views=upgraded,
        native_primitive_bounds=native, fourier=fourier_theorem(),
        original_target_and_control_C3_update=control_update,
        exact_shared_selected_N_checked=True,
        radial_scope='original r_minus to Rp', original_functions_and_selected_N_unchanged=True)


class CurrentOriginalN1PreRpPressureC1:
    @source_precision
    def __init__(self, require_checked=True):
        self.c = MPIntervalContext()
        self.c.dps = 540
        self.hashes = {}
        inputs = importlib.import_module(PREFIX+'current_original_n1_preRp_C1_inputs')
        ranges = load_checked(rg, self.hashes)
        self.family = copy.deepcopy(ranges['source_family'])
        self.parents = dict(
            ranges=ranges, old=load_checked(previous, self.hashes, self.family),
            sources=load_checked(inputs, self.hashes, self.family),
            partial_data=load_checked(partial, self.hashes, self.family),
            endpoints=load_checked(endpoint, self.hashes, self.family),
            old_control=load_checked(control_C3,self.hashes,self.family),
            old_band=load_checked(repair_C3,self.hashes,self.family))
        self.pulse_data = load_checked(pulse_C1, self.hashes, self.family)
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            if not (receipt.get('all_passed') and receipt.get(GATE) and receipt.get(UPPER_GATE)
                    and receipt['source_family'] == self.family and not any(receipt.get(k) for k in OPEN)):
                raise ValueError('Checked original pre-Rp C1 receipt required')
            for name, digest in receipt['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('Changed accepted pre-Rp C1 source '+name)
            raw = _read(importlib.import_module(Path(__file__).stem))
            self.data = raw['source_data']
            self.hashes.update(receipt['input_hashes'])
        else:
            self.data = record(build_data(self.c, **self.parents))
        self.error = {key: 2*self.c.exp(-1000) for key in ('m13', 'm15')}
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.acceptance_loaded = require_checked
        self._snapshot = json.dumps(record((self.data, self.family, self.hashes)), sort_keys=True)

    @source_precision
    def assert_graph(self):
        if self._snapshot != json.dumps(record((self.data, self.family, self.hashes)), sort_keys=True):
            raise ValueError('Original pre-Rp C1 data or identity changed')

    @source_precision
    def upper_C1_grid(self):
        self.assert_graph()
        grid = copy.deepcopy(self.pulse_data['c1_grid'])
        for model in grid:
            for key, field in model['fields'].items():
                field['preRp_nonzero_source_error'] = self.error[key]
                field['total_error'] = interval(self.c, field['total_error'])+self.error[key]
                field['absolute_accuracy_certified'] = ep(field['total_error'])[1] < ep(self.c.mpf('.001'))[0]
            model['radial_scope'] = 'original r_minus to infinity'
        return grid

    @source_precision
    def derivative_enclosure(self, Z):
        """Directed genuine pressure derivative enclosure, in original units."""
        self.assert_graph()
        z = self.c.mpf(Z)
        lo, hi = ep(z)
        if not (mp.isfinite(lo) and mp.isfinite(hi) and -1 <= lo <= hi <= 1):
            raise ValueError('Z must stay within the original axial domain [-1,1]')
        values = {key: [] for key in self.error}
        for model in self.pulse_data['c1_grid']:
            left, right = ep(interval(self.c, model['left']))[0], ep(interval(self.c, model['right']))[1]
            if right < lo or left > hi:
                continue
            for key, field in model['fields'].items():
                e = interval(self.c, field['total_error'])+self.error[key]
                values[key].append(interval(self.c, field['center_slope'])+self.c.mpf([-ep(e)[1], ep(e)[1]]))
        if any(not rows for rows in values.values()):
            raise ValueError('Original derivative grid coverage missing')
        return dict(Z=z, source_family=self.family, radial_scope='original r_minus to infinity',
            normalized_pressure_Z={key: self.c.mpf([min(ep(v)[0] for v in rows), max(ep(v)[1] for v in rows)])
                                   for key, rows in values.items()},
            preRp_nonzero_source_error=self.error, numerical_enclosures_not_coefficient_functions=True,
            **{GATE: self.acceptance_loaded, UPPER_GATE: self.acceptance_loaded}, **dict.fromkeys(OPEN, False))


@source_precision
def run():
    began = time.monotonic()
    owner = CurrentOriginalN1PreRpPressureC1(False)
    maxima = {key: max(ep(interval(owner.c, v['fields'][key]['total_error'])+owner.error[key])[1]
        for v in owner.pulse_data['c1_grid']) for key in owner.error}
    samples = {z: record(owner.derivative_enclosure(z)) for z in ('-1', '-.639', '0', '.211', '.371', '.639', '1')}
    raw = dict(source_family=owner.family, source_data=owner.data,
        pressure_C1_composition=dict(parent=pulse_C1.NAME, parent_sha256=sha(pulse_C1.NAME),
            grid_count=len(owner.pulse_data['c1_grid']), added_preRp_error=owner.error,
            maximum_normalized_derivative_error=maxima, radial_scope='original r_minus to infinity',
            unchanged_signed_center_slopes=True), samples=samples,
        source_bindings=dict(original_C3_signed_functions=rg.target.NAME,
            primitive_bound=Path(__file__).name+':sharpen_primitive',
            complete_physical_q=Path(physical.physical.__file__).name+':_physical_q_bounds',
            pressure_identity=Path(previous.reduced.__file__).name+':reduced_S_integrals',
            original_finite_endpoint_C1=endpoint.NAME),
        input_hashes=owner.hashes, **{GATE: False, UPPER_GATE: False}, **dict.fromkeys(OPEN, False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(record(raw), separators=(',', ':'))+'\n').encode(), mtime=0))
    print('PRE_RP_C1_CONSUMED', len(owner.data['upgraded_chart_views']), len(owner.data['native_primitive_bounds']), flush=True)
    return owner, raw


if __name__ == '__main__':
    run()
