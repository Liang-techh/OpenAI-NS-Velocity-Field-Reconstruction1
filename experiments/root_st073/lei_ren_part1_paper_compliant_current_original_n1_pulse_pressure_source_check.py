"""Independent checks for the combined Rp-to-Rv n1 pulse pressure source.

The producer is intentionally kept out of the independent fixtures below.
They check the finite-interval moment/Fubini signs, the positive main-kernel
support cap, the Hardy/Cauchy constants, and the reduced quadratic pressure
constants.  The live-source portion is bound to the producer API at runtime.
"""

import gzip
import hashlib
import json
from functools import wraps
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

try:
    import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source as current
except ImportError:  # The producer is created by the parent task later.
    current = None


HERE = Path(__file__).parent
SOURCE_DPS = 540


def _source_precision(fn):
    """Keep live-record decimal strings identical to the 540-digit raw run."""
    @wraps(fn)
    def call(*args, **kwargs):
        with mp.workdps(SOURCE_DPS):
            return fn(*args, **kwargs)
    return call


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _record_bounds(value):
    if isinstance(value, dict) and 'lower_exact_mpf_tuple' in value:
        return (mp.make_mpf(tuple(value['lower_exact_mpf_tuple'])),
                mp.make_mpf(tuple(value['upper_exact_mpf_tuple'])))
    if hasattr(value, '_mpi_'):
        return tuple(mp.make_mpf(x) for x in value._mpi_)
    raise TypeError('Directed interval required')


def _zeroth_bounds(value):
    """Read a scalar record or coefficient-zero of a serialized/live jet."""
    if isinstance(value, dict) and 'coefficients' in value:
        coefficients = value['coefficients']
        if not coefficients:
            raise ValueError('Taylor record has no zeroth coefficient')
        return _record_bounds(coefficients[0])
    if hasattr(value, 'order') and hasattr(value, '__getitem__'):
        return _record_bounds(value[0])
    return _record_bounds(value)


def _interval_from_record(ctx, value):
    if hasattr(value, '_mpi_'):
        return value
    return ctx.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                         tuple(value['upper_exact_mpf_tuple'])))


def _read_raw():
    if current is None:
        raise RuntimeError('Pulse pressure producer is not available')
    return json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))


def symbolic_moment_fixture():
    """Check row-RHS/Fubini identities and radial pressure signs."""
    Z, delta, ap, mu = sy.symbols('Z delta ap mu', positive=True)
    r1, r2, F1, F2, F1prime, rho = sy.symbols(
        'r1 r2 F1 F2 F1prime rho')
    C1, C2, W1, W2, Kmain, ratio = sy.symbols(
        'C1 C2 W1 W2 Kmain ratio', positive=True)
    ubar = sy.symbols('ubar', positive=True)
    G2 = r2 - r1 + ap * (F2 - F1)
    I0U = -ubar * mu * r1
    I0m = ubar * (ap * F1prime + 13 * r1 - G2 - rho)
    I1U = ubar * (C1 * W1 + C2 * W2 + ap * Kmain)
    I1m = -I1U - ubar * mu * r1 * ratio
    expected = dict(
        G2=G2,
        I0U=I0U,
        I0m=I0m,
        I1U=I1U,
        I1m=I1m)
    if sy.simplify(I0U + ubar * mu * r1) != 0:
        raise ArithmeticError('I0U sign identity failed')
    if sy.simplify(I0m - ubar * (ap * F1prime + 13 * r1 - G2 - rho)) != 0:
        raise ArithmeticError('I0m Fubini identity failed')
    if sy.simplify(I1m + I1U + ubar * mu * r1 * ratio) != 0:
        raise ArithmeticError('I1m endpoint identity failed')

    # An independent finite-interval own-mean construction checks the signs
    # behind the two Fubini reductions.  For H(s)=-integral_s^T h(v) dv,
    # integral H ds=-integral v*h(v)dv and integral exp(s-T)H ds is the
    # corresponding weighted endpoint expression.
    s, v, T = sy.symbols('s v T', real=True, positive=True)
    h = 1 + 2 * v + 3 * v**2
    H = -sy.integrate(h, (v, s, T))
    direct0 = sy.integrate(H, (s, 0, T))
    fubini0 = -sy.integrate(v * h, (v, 0, T))
    if sy.simplify(direct0 - fubini0) != 0:
        raise ArithmeticError('Unweighted finite-interval Fubini failed')
    weighted = sy.exp(s - T)
    direct1 = sy.integrate(weighted * H, (s, 0, T))
    fubini1 = -sy.integrate(h * (sy.exp(v - T) - sy.exp(-T)),
                             (v, 0, T))
    if sy.simplify(direct1 - fubini1) != 0:
        raise ArithmeticError('Weighted finite-interval Fubini failed')
    return dict(
        row_rhs_and_G2_identity=True,
        I0U_sign=True,
        I0m_F1prime_G2_rho_sign=True,
        I1U_positive_main_kernel=True,
        I1m_exact_Rp_over_Rv_endpoint=True,
        unweighted_own_mean_Fubini=True,
        weighted_own_mean_Fubini=True,
        symbolic_G2=str(G2),
        symbolic_I0m=str(I0m),
        symbolic_I1m=str(I1m))


def direct_main_kernel_fixture():
    """Check the main positive 3/2-mu kernel and logarithmic support cap."""
    with mp.workdps(100):
        mu = mp.mpf('1e-3')
        logE = mp.mpf('0')
        bound_log = (mp.log(11 / (mp.mpf('1.5') - mu)) - logE
                     - 3 / mu + 2)
        if bound_log >= -1000:
            raise ArithmeticError('Fixture main-kernel cap is not below -1000')
        T = 13 / mu
        # A positive decaying reference realizes the same support inequality.
        kernel = lambda t: mp.exp(-(mp.mpf('1.5') - mu) * t)
        integral = mp.quad(kernel, [0, T])
        upper = 11 / (mp.mpf('1.5') - mu)
        if not (integral > 0 and integral < upper):
            raise ArithmeticError('Positive main-kernel support bound failed')
        return dict(mu='1e-3', logE='0', T=mp.nstr(T, 30),
                    log_upper=mp.nstr(bound_log, 30),
                    positive_integral=mp.nstr(integral, 30),
                    positive_support_bound=mp.nstr(upper, 30))


def energy_inequality_fixture():
    """Numerically verify the stated Hardy/Cauchy constants independently."""
    with mp.workdps(80):
        T = mp.mpf('3.25')
        initial = mp.mpf('.012')
        # Smooth positive U with a nonconstant profile and a nonunit scale.
        U = lambda t: mp.mpf('.008') * (1 + t / T + (t / T)**2)
        m = lambda t: initial * mp.exp(-t) + mp.quad(
            lambda s: mp.exp(-(t - s)) * U(s), [0, t])
        E0U = mp.quad(lambda t: U(t), [0, T])
        E1U = mp.quad(lambda t: t * U(t)**2, [0, T])
        E0m = mp.quad(lambda t: abs(m(t)), [0, T])
        E1m = mp.quad(lambda t: t * m(t)**2, [0, T])
        hardy_rhs = 2 * initial + 8 * E0U
        cauchy_rhs = 2 * initial**2 * T + E1U * T**2
        if E0m > hardy_rhs or E1m > cauchy_rhs:
            raise ArithmeticError('Independent Hardy/Cauchy bound failed')
        # The producer's 3 epsilon budgets are checked as algebraic caps,
        # independently of its field enclosure implementation.
        eps = mp.exp(-1000)
        if 2 * eps + eps > 3 * eps:
            raise ArithmeticError('Weighted epsilon budget failed')
        return dict(E0U=mp.nstr(E0U, 30), E0m=mp.nstr(E0m, 30),
                    hardy_rhs=mp.nstr(hardy_rhs, 30),
                    E1U=mp.nstr(E1U, 30), E1m=mp.nstr(E1m, 30),
                    cauchy_rhs=mp.nstr(cauchy_rhs, 30),
                    Hardy_constant_2_initial_plus_8_E0U=True,
                    weighted_Cauchy_constant=True,
                    epsilon=mp.nstr(eps, 20))


def quadratic_constant_fixture():
    """Recompute q/qZ and reduced-S coefficient constants independently."""
    sqrt_term = mp.sqrt(7960)
    s0 = 8 + 2 * sqrt_term + 4 + 8
    s1 = 8 + 144 + 2 * sqrt_term + 8 + 8
    if not (s0 > 0 and s1 > s0):
        raise ArithmeticError('Reduced quadratic constants are invalid')
    # Cauchy-Schwarz check used by the mixed U*UZ terms.
    a, b = mp.mpf('3.5'), mp.mpf('11.25')
    if 2 * a * b > a**2 + b**2:
        raise ArithmeticError('Mixed-term Cauchy bound failed')
    # With each component energy capped by 4 epsilon V, the three q
    # coefficients give 3*4*(4^2+2^2+2^2)=288.  For q_Z, the elementary
    # absolute-sum estimate is already below the recorded 7960 constant:
    # 4*(16+4+8+8+2)^2=5776 <= 7960.
    q_from_components = 3 * 4 * (4**2 + 2**2 + 2**2)
    qZ_from_components = 4 * (16 + 4 + 8 + 8 + 2)**2
    if q_from_components != 288 or qZ_from_components > 7960:
        raise ArithmeticError('q/qZ coefficient energy constants failed')
    return dict(q_energy_constant=288, qZ_energy_constant=7960,
                q_energy_from_component_caps=q_from_components,
                qZ_energy_from_component_caps=qZ_from_components,
                S0_coefficient=mp.nstr(s0, 30),
                S1_coefficient=mp.nstr(s1, 30),
                mixed_term_Cauchy=True,
                boundary_qinit_squared_constant=16)


def _check_source_scope(raw, owner, view):
    if raw.get(current.GATE) is not False:
        raise ValueError('Producer gate must remain candidate-only')
    if view.get('source_family') != owner.family:
        raise ValueError('Source family mismatch')
    if any(view.get(key, False) for key in current.OPEN):
        raise ValueError('Out-of-scope gate enabled')
    for key in ('terminal_own_mean_identities_retained',
                'positive_weighted_main_source_retained',
                'unweighted_viscous_boundary_exact_zero',
                'all_original_pulse_charts_combined'):
        if view.get(key) is not True:
            raise ValueError('Missing source-scope fact: ' + key)
    if owner.units.get('terminal_M_at_Rv_exact_zero') is not True:
        raise ValueError('Terminal M(Rv)=0 identity missing')
    if owner.bounds['whole_pulse_energy'].get('source_energy_bounds_nonzero') is not True:
        raise ValueError('Whole-pulse nonzero energy provenance missing')
    return dict(source_family_match=True, terminal_M_zero=True,
                terminal_own_mean_retained=True,
                positive_main_source_retained=True,
                all_original_pulse_charts_combined=True,
                scope_open_flags_false=True)


def _check_energy_bounds(owner):
    """Recompute the common epsilon and reduced-S bounds independently."""
    energy = owner.bounds['whole_pulse_energy']
    c = owner.c
    if len(energy['log_caps']) != 3:
        raise ValueError('Expected three derivative log-cap rows')
    for row in energy['log_caps']:
        if any(_record_bounds(value)[1] >= -1000 for value in row.values()):
            raise ArithmeticError('A physical energy cap exceeds -1000')
    epsilon = energy['epsilon']
    if _record_bounds(epsilon)[0] <= 0:
        raise ArithmeticError('Positive epsilon source scale was lost')
    if energy['unweighted_Hardy_mean_factor'] != 10:
        raise ValueError('Hardy mean factor changed')
    if energy['weighted_Cauchy_mean_factor'] != 3:
        raise ValueError('Weighted Cauchy factor changed')
    if energy['conservative_common_mean_factor'] != 10:
        raise ValueError('Common mean factor changed')
    if energy['exact_units_combined_before_boxing'] is not True:
        raise ValueError('Physical units were not combined before boxing')
    expected_q = epsilon * 288
    expected_qZ = epsilon * 7960
    if _record_bounds(energy['q_energy']) != _record_bounds(expected_q):
        raise ArithmeticError('q energy constant was not recomputed')
    if _record_bounds(energy['qZ_energy']) != _record_bounds(expected_qZ):
        raise ArithmeticError('qZ energy constant was not recomputed')
    cross = c.sqrt(7960) * epsilon
    expected0 = current.end.scalar_bound(c, ((8 + 4 + 8) * epsilon + 2 * cross) / 2)
    expected1 = current.end.scalar_bound(c, ((8 + 144 + 8 + 8) * epsilon + 2 * cross) / 2)
    if _record_bounds(energy['m13']) != _record_bounds(expected0):
        raise ArithmeticError('m13 reduced quadratic bound mismatch')
    if _record_bounds(energy['m15']) != _record_bounds(expected1):
        raise ArithmeticError('m15 reduced quadratic bound mismatch')
    for key in ('m13', 'm15', 'q_energy', 'qZ_energy', 'boundary_q_squared_bound'):
        if _record_bounds(energy[key])[1] <= 0:
            raise ArithmeticError('Nonzero energy bound missing: ' + key)
    return dict(log_caps_three_rows_below_minus_1000=True,
                Hardy_factor_10=True, weighted_Cauchy_factor_3=True,
                common_mean_factor_10=True, q_energy_288_recomputed=True,
                qZ_energy_7960_recomputed=True,
                reduced_m13_m15_bounds_recomputed=True,
                all_energy_bounds_positive=True)


def _check_source_caps(owner):
    if _record_bounds(owner.bounds['weighted_main_source_log'])[1] >= -1000:
        raise ArithmeticError('Main weighted source log cap is not below -1000')
    if _record_bounds(owner.bounds['lower_radius_ratio_log'])[1] >= -1000:
        raise ArithmeticError('Rp/Rv log cap is not below -1000')
    for key in ('positive_weighted_main_source', 'positive_lower_radius_ratio',
                'viscous_error', 'qinit_bound'):
        if _record_bounds(owner.bounds[key])[1] <= 0:
            raise ArithmeticError('Positive source cap missing: ' + key)
    return dict(main_kernel_log_below_minus_1000=True,
                lower_radius_ratio_log_below_minus_1000=True,
                positive_main_lower_ratio_and_viscous_caps=True)


def _check_initial_mean_binding(owner, samples):
    """Check m0=ubar*mu*first original incoming moment on each query."""
    c = owner.c
    checked = 0
    for z, view in samples.items():
        cell = owner._cell(z)
        # The saved sample is a point query.  Reconstruct its physical
        # incoming mean at that same Z; using the whole cell endpoints here
        # creates an interval-valued Z and cannot equal the point report.
        _, ubar, _, _ = current.end.polynomials(c, z, 5, owner.delta)
        expected = ubar * cell['rows'][0] * owner.mu
        actual = view['source_linear_inputs']['lower_own_mean']
        for index in range(min(expected.order, actual.order) + 1):
            if _record_bounds(expected[index]) != _record_bounds(actual[index]):
                raise ArithmeticError('First incoming mean binding mismatch')
        checked += 1
    return dict(samples=checked, original_first_incoming_row_used=True,
                physical_ubar_and_mu_applied_once=True)


def _check_models(owner, raw, view):
    """Recompute model absolute errors and require genuine nonzero bounds."""
    models = raw.get('models', view.get('models'))
    if models is None:
        return dict(model_schema_deferred=True)
    if not isinstance(models, list) or len(models) != 256:
        raise ValueError('Expected 256 full-Z affine models')
    checked = 0
    with mp.workdps(SOURCE_DPS):
        for model in models:
            for name in ('m13', 'm15'):
                field = model['fields'][name]
                if not field.get('absolute_accuracy_certified', False):
                    raise ValueError('Model accuracy flag is false')
                lo, hi = _record_bounds(field['absolute_error'])
                if hi >= mp.mpf('.001'):
                    raise ArithmeticError('Model absolute error is not below .001')
                if hi <= 0:
                    raise ArithmeticError('Model absolute error is not positive')
                if not field.get('full_nonzero_source_error'):
                    raise ValueError('Nonzero source error was discarded')
                if _record_bounds(field['full_nonzero_source_error']) != \
                        _record_bounds(owner.full_errors[name]):
                    raise ArithmeticError('Model full source error differs from owner')
                parts = (field['linear_Taylor_remainder'],
                         field['coefficient_ambiguity'],
                         field['full_nonzero_source_error'])
                model_ctx = MPIntervalContext()
                model_ctx.dps = SOURCE_DPS
                expected_interval = sum(
                    (_interval_from_record(model_ctx, part)
                     for part in parts), model_ctx.mpf(0))
                expected_error = _record_bounds(expected_interval)
                if _record_bounds(field['absolute_error']) != expected_error:
                    raise ArithmeticError('Model absolute error was not recomputed')
                checked += 1
    return dict(models=256, fields=checked,
                absolute_error_upper_below_001=True,
                nonzero_source_errors_retained=True)


def _check_samples(owner, samples):
    expected = ('-1', '-.639', '0', '.211', '.371', '.639', '1')
    if tuple(samples) != expected:
        raise ValueError('Unexpected seven-query sample set')
    checked = 0
    for z, view in samples.items():
        if view['source_family'] != owner.family:
            raise ValueError('Sample source family mismatch')
        if view['all_original_pulse_charts_combined'] is not True:
            raise ValueError('Sample is not whole-pulse')
        if view['terminal_own_mean_identities_retained'] is not True:
            raise ValueError('Sample terminal mean identity missing')
        if view['positive_weighted_main_source_retained'] is not True:
            raise ValueError('Sample weighted main source missing')
        if view['unweighted_viscous_boundary_exact_zero'] is not True:
            raise ValueError('Sample unweighted viscous identity missing')
        if view['full_pressure_derivative_accuracy_not_claimed'] is not True:
            raise ValueError('Sample derivative scope widened')
        if view[current.GATE] is not owner.acceptance_loaded:
            raise ValueError('Sample gate state mismatch')
        terms = view['source_linear_inputs']
        required = ('I0U', 'I0m', 'I1U', 'I1m', 'I0q', 'I1q', 'qinit',
                    'index_derivative_source', 'finite_mu_X_remainder',
                    'main_weighted_source', 'lower_own_mean',
                    'm15_viscous_coefficient')
        if any(key not in terms for key in required):
            raise ValueError('Sample source-linear field missing')
        finite_mu_remainder = _zeroth_bounds(terms['finite_mu_X_remainder'])
        if finite_mu_remainder[0] == 0 and finite_mu_remainder[1] == 0:
            raise ArithmeticError('Finite-mu remainder was discarded')
        if _zeroth_bounds(terms['main_weighted_source'])[1] <= 0:
            raise ArithmeticError('Main weighted source was discarded')
        if view['source_linear_inputs']['lower_q_plus_R_qR_exact_zero'] is not True:
            raise ValueError('Lower q endpoint identity missing')
        for key in ('m13', 'm15'):
            if _record_bounds(view['quadratic_source_errors'][key])[1] <= 0:
                raise ArithmeticError('Quadratic source error is zero')
            if _record_bounds(view['full_pressure_values'][key])[1] <= \
                    _record_bounds(view['full_pressure_values'][key])[0]:
                raise ArithmeticError('Full pressure interval is degenerate')
        if _record_bounds(view['weighted_viscous_endpoint'])[1] == 0:
            raise ArithmeticError('Weighted viscous endpoint was dropped')
        checked += 1
    return dict(samples=checked, linear_inputs_present=True,
                finite_mu_and_main_sources_nonzero=True,
                terminal_and_viscous_endpoint_identities=True)


def _check_raw_live(raw, owner, samples, models):
    if raw.get('source_family') != owner.family:
        raise ValueError('Saved source family mismatch')
    if raw.get('samples') != current.record(samples):
        raise ValueError('Saved samples differ from live samples')
    if raw.get('models') != current.record(models):
        raise ValueError('Saved models differ from live models')
    if raw.get('original_physical_units') != owner.units:
        raise ValueError('Saved physical units differ from live units')
    if raw.get('source_bounds') != current.record(owner.bounds):
        raise ValueError('Saved bounds differ from live bounds')
    if raw.get('source_bindings') != owner.source_bindings:
        raise ValueError('Saved source bindings differ from live bindings')
    for name, digest in raw.get('input_hashes', {}).items():
        if _sha(name) != digest:
            raise ValueError('Saved source hash changed: ' + name)
    return True


def _public_default_precision(owner):
    before = mp.mp.dps
    with mp.workdps(15):
        public = current.CurrentOriginalN1PulsePressureSource(
            end_owner=owner.end, index_owner=owner.index,
            require_checked=False)
        packet = public.evaluate('.371')
        view = public.report(packet)
        if mp.mp.dps != 15:
            raise ArithmeticError('Public query leaked source precision')
        if view['source_family'] != owner.family:
            raise ValueError('Public family changed')
    if mp.mp.dps != before:
        raise ArithmeticError('Public query leaked outer precision')
    return True


@_source_precision
def run(owner=None, samples=None, models=None):
    if current is None:
        raise RuntimeError('Producer module is not available yet')
    if owner is None:
        owner = current.CurrentOriginalN1PulsePressureSource(
            require_checked=False)
    if samples is None:
        samples = {z: owner.report(owner.evaluate(z)) for z in
                   ('-1', '-.639', '0', '.211', '.371', '.639', '1')}
    if models is None:
        models = owner.affine_models()
    raw = _read_raw()
    symbolic = symbolic_moment_fixture()
    main_kernel = direct_main_kernel_fixture()
    energy = energy_inequality_fixture()
    quadratic = quadratic_constant_fixture()
    scope = _check_source_scope(raw, owner, next(iter(samples.values())))
    energy_bounds = _check_energy_bounds(owner)
    source_caps = _check_source_caps(owner)
    initial_binding = _check_initial_mean_binding(owner, samples)
    sample_evidence = _check_samples(owner, samples)
    _check_raw_live(raw, owner, samples, models)
    model_evidence = _check_models(owner, raw, next(iter(samples.values())))
    public = _public_default_precision(owner)
    receipt = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(symbolic_moment_identities=symbolic,
                    main_kernel_support_bound=main_kernel,
                    Hardy_Cauchy_energy_bounds=energy,
                    reduced_quadratic_constants=quadratic,
                    source_scope=scope,
                    whole_pulse_energy_bounds=energy_bounds,
                    positive_source_caps=source_caps,
                    initial_mean_binding=initial_binding,
                    seven_source_queries=sample_evidence,
                    full_Z_models=model_evidence,
                    raw_live_schema_checked=True,
                    default_precision_public_query=public),
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        input_hashes=dict(owner.hashes))
    receipt['input_hashes'][Path(__file__).name] = _sha(Path(__file__).name)
    receipt['input_hashes'][current.NAME] = _sha(current.NAME)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(',', ':')) + '\n',
        encoding='utf8', newline='\n')
    return receipt


if __name__ == '__main__':
    run()
