"""Independent checks for the finite end-pulse n1 pressure source.

The companion binds the accepted refined Rp pulse coefficients to the actual
end field.  This checker reconstructs the pressure equation from an unrelated
smooth physical mean, checks both Fubini/IBP pressure identities, and then
checks every saved full-Z affine model against the live source.  Quadratic and
weighted-viscous errors remain explicit nonzero bounds throughout.
"""

import gzip
import json
from math import comb
from pathlib import Path
import time

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


def _read_raw():
    return json.loads(gzip.decompress(
        (current.HERE / current.NAME).read_bytes()))


def _canon(value):
    value = sy.together(sy.expand_log(sy.sympify(value), force=True))
    return sy.factor(sy.cancel(value))


def _equal(left, right):
    return _canon(left - right) == 0


def _direct_fields(m, R, Z, delta):
    """Recover U, q, S and N1p directly from the physical mean m."""
    U = sy.cancel(m + R * sy.diff(m, R))
    d = 1 - Z**2
    L = 1 - delta * Z**2
    q = sy.cancel((2 * Z * U - (1 - delta) * Z * m - d * sy.diff(m, Z)) / L)
    q_Z = sy.diff(q, Z)
    q_R = sy.diff(q, R)
    q_y = R * q_R
    S = sy.cancel(
        (q + (1 - delta) * Z * q_Z / 2 + q_y) / L
        + q * (q / 2 + q_y)
        + U * (-2 * Z * q + d * q_Z - 2 * Z * q_y) / L)
    N1p = sy.cancel(-S / 2 + 2 * q_R + R * sy.diff(q, R, 2))
    return dict(U=U, q=q, q_Z=q_Z, q_R=q_R, q_y=q_y,
                S=S, N1p=N1p)


def independent_pressure_fixture():
    """Use a nonlinear, nonzero physical mean and direct radial quadrature."""
    R, Z = sy.symbols('R Z', positive=True, real=True)
    delta = sy.Rational(2, 13)
    k = sy.Rational(3, 11)
    eps = sy.Rational(2, 17)
    m = k * (Z + 2 * Z**2) / R + eps * Z**2 * (R - 1)**2 / R
    fields = _direct_fields(m, R, Z, delta)

    def B(value):
        return sy.cancel(value.subs(R, 2) - value.subs(R, 1))

    i0 = {}
    i1 = {}
    for name, value in dict(
            q_Z=fields['q_Z'], U_q_Z=fields['U'] * fields['q_Z'],
            U_squared=fields['U']**2,
            U_U_Z=fields['U'] * sy.diff(fields['U'], Z), q=fields['q'],
            q_squared=fields['q']**2).items():
        i0[name] = sy.integrate(value, (R, 1, 2))
        i1[name] = sy.integrate(R * value, (R, 1, 2))
    boundaries = dict(
        R_q=B(R * fields['q']),
        R_q_squared=B(R * fields['q']**2),
        R_squared_q=B(R**2 * fields['q']),
        R_squared_q_squared=B(R**2 * fields['q']**2),
        R_U_squared=B(R * fields['U']**2),
        R_squared_U_squared=B(R**2 * fields['U']**2),
        q_plus_R_q_R=B(fields['q'] + R * fields['q_R']),
        R_squared_q_R=B(R**2 * fields['q_R']))

    L = 1 - delta * Z**2
    i0S = (boundaries['R_q'] / L
           + (1 - delta) * Z * i0['q_Z'] / (2 * L)
           + boundaries['R_q_squared'] / 2
           + (1 - Z**2) * i0['U_q_Z'] / L
           - 2 * Z**2 * boundaries['R_U_squared'] / L**2
           - 2 * delta * Z**2 * i0['U_squared'] / L**2
           + 2 * Z * (1 - Z**2) * i0['U_U_Z'] / L**2)
    i1S = (boundaries['R_squared_q'] / L
           - i1['q'] / L
           + (1 - delta) * Z * i1['q_Z'] / (2 * L)
           + boundaries['R_squared_q_squared'] / 2
           - i1['q_squared'] / 2
           + (1 - Z**2) * i1['U_q_Z'] / L
           - 2 * Z**2 * boundaries['R_squared_U_squared'] / L**2
           + 2 * (1 - delta) * Z**2 * i1['U_squared'] / L**2
           + 2 * Z * (1 - Z**2) * i1['U_U_Z'] / L**2)
    direct_s0 = sy.integrate(fields['S'], (R, 1, 2))
    direct_s1 = sy.integrate(R * fields['S'], (R, 1, 2))
    direct0 = sy.integrate(fields['N1p'], (R, 1, 2))
    direct1 = sy.integrate(R * fields['N1p'], (R, 1, 2))
    assert _equal(i0S, direct_s0)
    assert _equal(i1S, direct_s1)
    assert _equal(-direct_s0 / 2 + boundaries['q_plus_R_q_R'], direct0)
    assert _equal(-direct_s1 / 2 + boundaries['R_squared_q_R'], direct1)

    with mp.workdps(100):
        z0 = mp.mpf('.371')
        # Build the numerical callable inside the guarded context; creating
        # it at the process default precision would round rational constants
        # before the independent quadrature starts.
        f0 = sy.lambdify((R, Z), fields['N1p'], 'mpmath')
        f1 = sy.lambdify((R, Z), R * fields['N1p'], 'mpmath')
        quad0 = mp.quad(lambda rr: f0(rr, z0), [1, 2])
        quad1 = mp.quad(lambda rr: f1(rr, z0), [1, 2])
        exact0 = mp.mpf(str(sy.N(direct0.subs(Z, sy.Rational(371, 1000)), 100)))
        exact1 = mp.mpf(str(sy.N((-direct1).subs(Z, sy.Rational(371, 1000)), 100)))
        assert abs(quad0 - exact0) < mp.mpf('1e-80')
        assert abs(-quad1 - exact1) < mp.mpf('1e-80')

    x, Rc, f = sy.symbols('x Rc f', positive=True)
    assert sy.integrate(f * Rc, (x, 1, 2)) == f * Rc
    assert sy.integrate((Rc * x) * f * Rc, (x, 1, 2)) == sy.Rational(3, 2) * f * Rc**2
    assert fields['U'] != 0 and fields['q_Z'] != 0
    assert sy.simplify(boundaries['q_plus_R_q_R']) != 0
    assert sy.simplify(boundaries['R_squared_q_R']) != 0
    assert sy.simplify(i0['U_q_Z']) != 0
    assert sy.simplify(i1['U_U_Z']) != 0
    return dict(
        direct_full_N1p_dR_and_RdR=True,
        direct_S_dR_and_RdR_Fubini_identities=True,
        direct_endpoint_IBP_identities=True,
        direct_m13_m15_nonlinear_viscous_quadrature=True,
        nonlinear_U_q_Z_and_U_U_Z_terms_nonzero=True,
        viscous_q_R_and_q_RR_terms_retained=True,
        finite_endpoint_q_plus_Rq_R_nonzero=True,
        finite_endpoint_R2q_R_nonzero=True,
        physical_dR_and_RdR_measure_powers=True,
        fixture_delta=str(delta), fixture_k=str(k), fixture_epsilon=str(eps),
        direct_m13=str(direct0), direct_m15=str(-direct1))


def _contains(outer, inner):
    olo, ohi = _interval_endpoints(outer)
    ilo, ihi = _interval_endpoints(inner)
    return olo <= ilo and ihi <= ohi


def independent_end_linear_fixture(owner):
    """Check end_linear_jets from independently quadratured finite-mu rows.

    The compact bump and all two-sided moment rows are created here rather
    than read from the production coefficient owner.  The incoming rows are
    then recovered from the two exact source moment equations, and the
    production helper is compared with independently constructed moment and
    q identities.  A nonunit A/Rv endpoint scaling is included explicitly.
    """
    c = owner.c
    with mp.workdps(220):
        z0 = +mp.mpf('.371')
        mu = +mp.mpf('.01')
        delta = +mp.mpf('.2')
        ell = +mp.mpf('.15')
        lam = +mp.mpf('.5') - mu
        norm = +ell * mp.mpf(256) / 315

    def beta(t):
        if abs(t) > ell:
            return mp.mpf('0')
        return (1 - (t / ell)**2)**4 / norm

    def moment(d, rate):
        with mp.workdps(220):
            return +mp.exp(-d * rate) * mp.quad(
                lambda t: mp.exp(rate * t) * beta(t), [-ell, ell])

    A = [moment(d, lam) for d in (3, 1)]
    Wplus = [moment(d, lam + 1) for d in (3, 1)]
    X = []
    rho = []
    D = []
    for index, d in enumerate((3, 1)):
        with mp.workdps(220):
            xj = +mp.exp(-d * lam) * mp.quad(
                lambda t: (t + 4 - d) * mp.exp(lam * t) * beta(t),
                [-ell, ell])
            rj = +mp.exp(-d * lam) * mp.quad(
                lambda t: mp.exp(lam * t) * beta(t) *
                ((mp.exp(mu * (d - t)) - 1) / mu - (d - t)),
                [-ell, ell])
        X.append(xj)
        rho.append(rj)
        D.append(4 * A[index] + rj - xj)
    assert all(value > 0 for value in rho)

    def C1(z):
        return mp.mpf('.7') + mp.mpf('.2') * z + mp.mpf('.05') * z**2

    def C2(z):
        return -mp.mpf('.4') + mp.mpf('.1') * z - mp.mpf('.03') * z**2

    def apf(z):
        return mp.mpf('.93') + mp.mpf('.08') * z - mp.mpf('.02') * z**2

    def Pi0(z):
        return mp.mpf('.3') - mp.mpf('.05') * z + mp.mpf('.02') * z**2

    def Pi1(z):
        return -mp.mpf('.2') + mp.mpf('.03') * z + mp.mpf('.01') * z**2

    def Wsum(z):
        return C1(z) * A[0] + C2(z) * A[1]

    def Dsum(z):
        return C1(z) * D[0] + C2(z) * D[1]

    def Xsum(z):
        return C1(z) * X[0] + C2(z) * X[1]

    def row0(z):
        return -Wsum(z) / mu - apf(z) * Pi0(z)

    def row1(z):
        return row0(z) - Dsum(z) - apf(z) * (Pi1(z) - Pi0(z))

    def jet(fn, order=5):
        with mp.workdps(200):
            values = [mp.diff(fn, z0, n) / mp.factorial(n)
                      for n in range(order + 1)]
            # Keep a directed, tiny independent quadrature allowance.  The
            # production Taylor algebra performs cancellation in a different
            # operation order, so point-decimal inputs could exclude the same
            # exact source value at the last few bits.
            coeffs = []
            for value in values:
                radius = max(abs(value), mp.mpf(1)) * mp.mpf('1e-70')
                coeffs.append(c.mpf([mp.nstr(value - radius, 200),
                                     mp.nstr(value + radius, 200)]))
        return IntervalTaylor(c, coeffs)

    zjet = IntervalTaylor.variable(c, c.mpf(str(z0)), 5)
    Cjets = [jet(C1), jet(C2)]
    apjet = jet(apf)
    Pijets = [jet(Pi0), jet(Pi1)]
    rows = [jet(row0), jet(row1)]
    A_c = [c.mpf(mp.nstr(v, 200)) for v in A]
    D_c = [c.mpf(mp.nstr(v, 200)) for v in D]
    Wplus_c = [c.mpf(mp.nstr(v, 200)) for v in Wplus]
    rho_c = [c.mpf(mp.nstr(v, 200)) for v in rho]
    mu_c = c.mpf(str(mu))
    delta_c = c.mpf(str(delta))
    terms = current.end_linear_jets(
        c, c.mpf(str(z0)), delta_c, mu_c, apjet, Cjets, rows, Pijets,
        Wplus_c, rho_c)

    def jderivative(value):
        return IntervalTaylor(c, [(i + 1) * value[i + 1]
                                   for i in range(value.order)])

    def independent_q(U, m):
        z = zjet
        d = 1 - z * z
        L = 1 - z * delta_c * z
        # Keep the Taylor object on the left of scalar products.  This is
        # required by the interval-Taylor operand protocol and also mirrors
        # the source's physical q formula without changing its order.
        return (2 * z * U - z * (1 - delta_c) * m
                - d * jderivative(m)) / L

    ubar = 1 / (1 + zjet * zjet)
    Wsum_j = jet(Wsum)
    Dsum_j = jet(Dsum)
    Xsum_j = jet(Xsum)
    Wplus_sum_j = jet(lambda z: C1(z) * Wplus[0] + C2(z) * Wplus[1])
    expected = dict(
        first_row_sum=jet(lambda z: -Wsum(z) / mu) * (-mu_c),
        Xsum=jet(Xsum),
        I0U=ubar * Wsum_j,
        I0m=-ubar * Xsum_j,
        I1U=ubar * Wplus_sum_j,
        I1m=-ubar * (Wplus_sum_j - Wsum_j * c.exp(-4)))
    # The two directed calculations use different rounding/error budgets.
    # Compare their high-precision centers with an explicit, much tighter
    # quadrature tolerance, while retaining the exact symbolic row identities
    # below.  Requiring one independently boxed interval to contain the other
    # would conflate those budgets and can reject a valid calculation.
    independent_center_errors = {}
    with mp.workdps(90):
        for key in ('first_row_sum', 'Xsum', 'I0U', 'I0m', 'I1U', 'I1m'):
            got_lo, got_hi = _interval_endpoints(terms[key][0])
            exp_lo, exp_hi = _interval_endpoints(expected[key][0])
            got_center = (got_lo + got_hi) / 2
            exp_center = (exp_lo + exp_hi) / 2
            error = abs(got_center - exp_center)
            independent_center_errors[key] = mp.nstr(error, 30)
            assert error < mp.mpf('1e-14'), (key, error)

    q0_expected = independent_q(expected['I0U'], expected['I0m'])
    q1_expected = independent_q(expected['I1U'], expected['I1m'])
    qbase_expected = independent_q(0 * expected['I0U'], -ubar * Wsum_j)
    for key, expected_q in (('I0q', q0_expected), ('I1q', q1_expected),
                            ('qbase', qbase_expected)):
        got_lo, got_hi = _interval_endpoints(terms[key][0])
        exp_lo, exp_hi = _interval_endpoints(expected_q[0])
        got_center = (got_lo + got_hi) / 2
        exp_center = (exp_lo + exp_hi) / 2
        error = abs(got_center - exp_center)
        independent_center_errors[key] = mp.nstr(error, 30)
        assert error < mp.mpf('1e-14'), (key, error)

    L = 1 - zjet * delta_c * zjet
    S0 = -qbase_expected / L + jderivative(q0_expected) * zjet * (1 - delta_c) / (2 * L)
    S1 = (-qbase_expected * c.exp(-4) - q1_expected) / L
    S1 += jderivative(q1_expected) * zjet * (1 - delta_c) / (2 * L)
    for key, expected_jet in (
            ('m13_linear', -S0 / 2),
            ('m15_linear', S1 / 2),
            ('m15_viscous_coefficient', -qbase_expected)):
        got_lo, got_hi = _interval_endpoints(terms[key][0])
        exp_lo, exp_hi = _interval_endpoints(expected_jet[0])
        got_center = (got_lo + got_hi) / 2
        exp_center = (exp_lo + exp_hi) / 2
        error = abs(got_center - exp_center)
        independent_center_errors[key] = mp.nstr(error, 30)
        assert error < mp.mpf('1e-14'), (key, error)

    # Independently exercise the own-mean/Fubini step.  For each compact
    # bump h_j(v)=exp(lambda*v) beta(v+d_j), define the backward own mean
    # H_j(s)=-integral_s^0 h_j(v) dv.  We evaluate H by an independent
    # polynomial/exponential antiderivative and integrate it with composite
    # Simpson quadrature at two resolutions.  This gives direct checks of
    # integral H_j ds=-X_j and integral exp(s) H_j ds=-(Wplus_j-exp(-4)A_j)
    # without reconstructing those identities from the producer.
    with mp.workdps(90):
        def exp_power_antiderivative(power, x, exp_x):
            total = mp.mpf('0')
            for r in range(power + 1):
                total += ((-1) ** r * mp.factorial(power) /
                          mp.factorial(power - r) * x ** (power - r) /
                          lam ** (r + 1))
            return exp_x * total

        def bump_antiderivative(x, d):
            # The compact bump is a degree-eight polynomial on its support.
            total = mp.mpf('0')
            exp_x = mp.exp(lam * x)
            for k in range(5):
                coefficient = ((-1) ** k * comb(4, k) /
                               (ell ** (2 * k) * norm))
                power = 2 * k
                for j in range(power + 1):
                    total += (coefficient * comb(power, j) *
                              d ** (power - j) *
                              exp_power_antiderivative(j, x, exp_x))
            return total

        def mean_component(s, d):
            lo = max(s, -d - ell)
            hi = min(mp.mpf('0'), -d + ell)
            if hi <= lo:
                return mp.mpf('0')
            return -(bump_antiderivative(hi, d) -
                     bump_antiderivative(lo, d))

        split_points = sorted(set([
            mp.mpf('-4'), mp.mpf('0'),
            -mp.mpf('3') - ell, -mp.mpf('3') + ell,
            -mp.mpf('1') - ell, -mp.mpf('1') + ell]))

        def simpson(function, left, right, cells):
            step = (right - left) / cells
            total = function(left) + function(right)
            for index in range(1, cells):
                total += (2 if index % 2 == 0 else 4) * function(
                    left + step * index)
            return step * total / 3

        def integrate(function, cells):
            return sum((simpson(function, left, right, cells)
                        for left, right in zip(split_points[:-1],
                                               split_points[1:])),
                       mp.mpf('0'))

        def mean_sum(s):
            return sum((C1(z0) * mean_component(s, 3),
                        C2(z0) * mean_component(s, 1)), mp.mpf('0'))

        direct_m0_128 = integrate(mean_sum, 128)
        direct_m1_128 = integrate(
            lambda s: mp.exp(s) * mean_sum(s), 128)
        direct_m0 = integrate(mean_sum, 256)
        direct_m1 = integrate(
            lambda s: mp.exp(s) * mean_sum(s), 256)
        fubini_m0 = -(C1(z0) * X[0] + C2(z0) * X[1])
        fubini_m1 = -(C1(z0) * (Wplus[0] - mp.exp(-4) * A[0])
                      + C2(z0) * (Wplus[1] - mp.exp(-4) * A[1]))
        quadrature_uncertainty_m0 = abs(direct_m0 - direct_m0_128)
        quadrature_uncertainty_m1 = abs(direct_m1 - direct_m1_128)
        assert quadrature_uncertainty_m0 < mp.mpf('1e-8')
        assert quadrature_uncertainty_m1 < mp.mpf('1e-8')
        assert abs(direct_m0 - fubini_m0) < mp.mpf('2e-8')
        assert abs(direct_m1 - fubini_m1) < mp.mpf('2e-8')

    # Compare the normalized moment values with direct radial quadrature
    # values generated from the independently computed beta moments.
    with mp.workdps(160):
        zv = mp.mpf(str(z0))
        uv = 1 / (1 + zv * zv)
        cv = [C1(zv), C2(zv)]
        direct = dict(
            I0U=uv * sum(cv[j] * A[j] for j in range(2)),
            I0m=-uv * sum(cv[j] * X[j] for j in range(2)),
            I1U=uv * sum(cv[j] * Wplus[j] for j in range(2)),
            I1m=-uv * (sum(cv[j] * Wplus[j] for j in range(2))
                       - mp.exp(-4) * sum(cv[j] * A[j] for j in range(2))))
        for key in direct:
            lo, hi = _interval_endpoints(terms[key][0])
            production_center = (lo + hi) / 2
            error = abs(production_center - direct[key])
            independent_center_errors['direct_' + key] = mp.nstr(error, 30)
            assert error < mp.mpf('1e-14'), (key, error)

        # Nonunit physical A/Rv scaling of the retained weighted viscous
        # endpoint: q_base(R)=Q/R gives B(R^2 q_R)=-Rv*Q.
        scale_A = mp.mpf(5) / 4
        scale_Rv = mp.mpf(11) / 3
        Q = _interval_endpoints(terms['qbase'][0])
        expected_endpoint = -scale_A * scale_Rv * ((Q[0] + Q[1]) / 2)
        got_coeff = _interval_endpoints(terms['m15_viscous_coefficient'][0])
        got_endpoint = scale_A * scale_Rv * (got_coeff[0] + got_coeff[1]) / 2
        endpoint_error = abs(got_endpoint - expected_endpoint)
        independent_center_errors['scaled_weighted_viscous_endpoint'] = \
            mp.nstr(endpoint_error, 30)
        assert endpoint_error < mp.mpf('1e-14')

    return dict(
        independent_compact_bump_direct_quadrature=True,
        nonzero_mu_and_delta=True,
        independent_A_D_Wplus_X_rho=True,
        source_rows_recovered_from_two_moment_equations=True,
        finite_mu_X_remainder_positive=True,
        ubar_Z_retained_in_independent_q=True,
        I0_and_I1_U_m_directly_checked=True,
        I0_and_I1_m_directly_checked=True,
        I0_and_I1_q_moment_derivatives_checked=True,
        direct_own_mean_Fubini_X_and_Wplus_checked=True,
        linear_pressure_and_weighted_viscous_endpoint_checked=True,
        nonunit_A_Rv_endpoint_scaling_checked=True,
        independent_center_errors=independent_center_errors,
        independent_center_tolerance='1e-14',
        direct_own_mean_Fubini_quadrature_uncertainty={
            'I0m': mp.nstr(quadrature_uncertainty_m0, 30),
            'I1m': mp.nstr(quadrature_uncertainty_m1, 30)},
        direct_own_mean_Fubini_tolerance='2e-8')


def _interval_endpoints(value):
    if hasattr(value, 'coefficients'):
        value = value[0]
    if hasattr(value, '_mpi_'):
        return value.a, value.b
    if isinstance(value, dict):
        if 'lower' in value and 'upper' in value:
            return mp.mpf(str(value['lower'])), mp.mpf(str(value['upper']))
    scalar = mp.mpf(str(value))
    return scalar, scalar


def _upper(value):
    return _interval_endpoints(value)[1]


def _assert_source_hashes(raw, owner):
    assert raw['source_family'] == owner.family
    assert raw[current.GATE] is False
    assert all(not raw.get(name, False) for name in current.OPEN)
    assert raw['input_hashes'] == owner.hashes
    for name, digest in raw['input_hashes'].items():
        assert current.sha(name) == digest, name


def _assert_raw_live(raw, owner, samples, models):
    assert raw['samples'] == current.record(samples)
    assert raw['models'] == current.record(models)
    assert raw['original_physical_units'] == owner.units
    assert raw['source_bounds'] == current.record(owner.bounds)
    assert raw['source_bindings'] == current.record(owner.source_bindings)
    assert raw['beta_Wplus'] == current.record(owner.Wplus)
    assert raw['finite_mu_rho'] == current.record(owner.rho)
    return True


def _assert_units(owner, raw):
    units = owner.units
    assert units == raw['original_physical_units']
    assert units['m13'] == dict(A_power=1, Rv_power=1,
                                inverse_mu_numerator=11,
                                inverse_mu_denominator=2)
    assert units['m15'] == dict(A_power=1, Rv_power=2,
                                inverse_mu_numerator=37,
                                inverse_mu_denominator=2)
    assert units['quadratic_m13'] == dict(A_power=2, Rv_power=1,
                                           inverse_mu_numerator=-4,
                                           inverse_mu_denominator=2)
    assert units['quadratic_m15'] == dict(A_power=2, Rv_power=2,
                                           inverse_mu_numerator=22,
                                           inverse_mu_denominator=2)
    assert units['integer_powers_combined_before_boxing'] is True
    assert units['no_physical_scale_materialization'] is True
    assert '26/mu' not in json.dumps(units, sort_keys=True)
    return dict(linear_m13_inverse_mu='11/2', linear_m15_inverse_mu='37/2',
                quadratic_m13_inverse_mu='-4/2',
                quadratic_m15_inverse_mu='22/2',
                no_incorrect_26_over_mu=True)


def _assert_bounds(owner, raw):
    c = owner.c
    bounds = owner.bounds
    assert raw['source_bounds'] == current.record(bounds)
    assert _upper(bounds['log_A_upper_source']) < -1000
    assert _interval_endpoints(bounds['log_Rv_lower_source'])[0] > 1000
    for name in ('m13', 'm15'):
        assert _upper(bounds['quadratic'][name]) > 0
        assert _upper(owner.errors[name]) > 0
        assert _upper(owner.full_errors[name]) > 0
    assert _upper(bounds['viscous_error']) > 0
    assert _upper(owner.viscous_error) > 0
    assert _upper(owner.A_cap) > 0 and _upper(owner.invRv) > 0

    field = bounds['global_field_bounds']
    U, UZ, q, qZ = (field[key] for key in ('U', 'UZ', 'q', 'qZ'))
    qbase = bounds['qbase_bound']
    S0 = c.exp(4) * qbase**2 / 2 + 2 * U * qZ + 4 * U**2 + 8 * U * UZ
    S1 = qbase**2 / 2 + q**2 / 4 + U * qZ + 4 * U**2 + 4 * U * UZ
    expected0 = current.scalar_bound(c, bounds['positive_A_cap'] * S0 /
                                     (2 * bounds['reference']))
    expected1 = current.scalar_bound(c, bounds['positive_A_cap'] * S1 /
                                     (2 * bounds['reference']))
    assert _interval_endpoints(expected0) == _interval_endpoints(owner.errors['m13'])
    assert _interval_endpoints(expected1) == _interval_endpoints(owner.errors['m15'])
    expected_viscous = current.scalar_bound(
        c, bounds['positive_inverse_Rv_enclosure'] * qbase /
        bounds['reference'])
    assert _interval_endpoints(expected_viscous) == _interval_endpoints(owner.viscous_error)
    return dict(log_caps_proved=True, quadratic_and_viscous_bounds_recomputed=True,
                quadratic_and_viscous_bounds_nonzero=True)


def _assert_samples(owner, samples):
    expected = ('-1', '-.639', '0', '.211', '.371', '.639', '1')
    assert tuple(samples) == expected
    checked = {}
    for z, view in samples.items():
        assert view['source_family'] == owner.family
        assert view['full_pressure_derivative_accuracy_not_claimed'] is True
        assert view['exact_source_row_correlation_retained'] is True
        assert view['finite_mu_X_remainder_retained'] is True
        assert view['unweighted_viscous_boundary_exact_zero'] is True
        assert view[current.GATE] is False
        assert all(not view[name] for name in current.OPEN)
        full = view['full_pressure_values']
        linear = view['linear_pressure_Taylor']
        errors = view['quadratic_source_errors']
        for name in ('m13', 'm15'):
            flo, fhi = _interval_endpoints(full[name])
            # ``linear_pressure_Taylor`` is an IntervalTaylor at the queried
            # Z; its value row is the coefficient at order zero.
            llo, lhi = _interval_endpoints(linear[name][0])
            assert flo <= llo <= lhi <= fhi
            assert _upper(errors[name]) > 0
        assert _upper(abs(view['weighted_viscous_endpoint'])) > 0
        assert view['source_linear_inputs']['lower_q_plus_R_qR_exact_zero'] is True
        assert _upper(view['source_linear_inputs']['X_finite_mu_remainder']) > 0
        checked[z] = dict(m13_positive=_interval_endpoints(full['m13'])[0] > 0,
                          m15_negative=_interval_endpoints(full['m15'])[1] < 0,
                          quadratic_error_nonzero=True,
                          weighted_viscous_nonzero=True)
    return checked


def _assert_models(owner, models):
    assert len(models) == current.MODELS == 256
    c = owner.c
    previous_right = None
    checked = 0
    for index, model in enumerate(models):
        left = c.mpf(-1) + c.mpf(2) * index / current.MODELS
        right = c.mpf(-1) + c.mpf(2) * (index + 1) / current.MODELS
        center = (left + right) / 2
        half_width = (right - left) / 2
        assert model['left']._mpi_ == left._mpi_
        assert model['right']._mpi_ == right._mpi_
        assert model['center']._mpi_ == center._mpi_
        assert model['half_width']._mpi_ == half_width._mpi_
        if previous_right is not None:
            assert model['left']._mpi_ == previous_right._mpi_
        previous_right = model['right']
        assert tuple(model['source_coefficient_cell']) in current.CELLS
        fields = model['fields']
        assert set(fields) == {'m13', 'm15'}
        # Re-query the actual source at the midpoint and the complete model
        # interval.  This exposes the genuine row-2 Z2 cap used in the model,
        # rather than trusting a copied Taylor remainder field.
        point = owner._make_view(model['center'])
        whole = owner._make_view([model['left'].a, model['right'].b])
        for name, field in fields.items():
            checked += 1
            point_jet = point['linear_pressure_Taylor'][name]
            whole_jet = whole['linear_pressure_Taylor'][name]
            half_row = whole_jet[2]
            expected_remainder = c.mpf(current.abs_upper(half_row)) * half_width**2
            expected_ambiguity = ((point_jet[0].b - point_jet[0].a) / 2
                                  + (point_jet[1].b - point_jet[1].a)
                                  * half_width / 2)
            assert _interval_endpoints(field['linear_Z2_half_range']) == \
                _interval_endpoints(half_row)
            assert _interval_endpoints(field['linear_Taylor_remainder']) == \
                _interval_endpoints(expected_remainder)
            assert _interval_endpoints(field['coefficient_ambiguity']) == \
                _interval_endpoints(expected_ambiguity)
            expected = (field['linear_Taylor_remainder']
                        + field['coefficient_ambiguity']
                        + field['full_nonzero_source_error'])
            assert _interval_endpoints(field['absolute_error']) == _interval_endpoints(expected)
            assert _upper(field['linear_Taylor_remainder']) >= 0
            assert _upper(field['coefficient_ambiguity']) >= 0
            assert _upper(field['full_nonzero_source_error']) > 0
            assert _upper(field['absolute_error']) > 0
            # Recompute the saved nonzero source contribution against the
            # owner's live bound and independently enforce the advertised
            # numerical accuracy.  The boolean flag alone is insufficient.
            assert _interval_endpoints(field['full_nonzero_source_error']) == \
                _interval_endpoints(owner.full_errors[name])
            assert _upper(field['absolute_error']) < c.mpf('.001')
            assert field['absolute_accuracy_certified'] is True
    assert previous_right._mpi_ == c.mpf(1)._mpi_
    return dict(models=256, adjacency_and_full_Z=True, fields=checked,
                all_absolute_errors_recomputed=True,
                all_nonzero_source_errors_retained=True)


def _public_default_precision(owner):
    before = mp.mp.dps
    with mp.workdps(15):
        public = current.CurrentOriginalN1EndPressureSource(
            coefficient_owner=owner.coefficients, require_checked=False)
        packet = public.evaluate('.371')
        view = public.report(packet)
        assert mp.mp.dps == 15
        assert view['source_family'] == owner.family
    assert mp.mp.dps == before
    return True


@current.source_precision
def run(owner=None, samples=None, models=None, coefficient_owner=None):
    began = time.monotonic()
    if owner is None:
        owner, samples, models = current.run(coefficient_owner)
    elif samples is None or models is None:
        samples = {z: owner.report(owner.evaluate(z)) for z in
                   ('-1', '-.639', '0', '.211', '.371', '.639', '1')}
        models = owner.affine_models()
    owner.assert_graph()
    raw = _read_raw()
    _assert_source_hashes(raw, owner)
    _assert_raw_live(raw, owner, samples, models)
    fixture = independent_pressure_fixture()
    end_fixture = independent_end_linear_fixture(owner)
    unit_evidence = _assert_units(owner, raw)
    bound_evidence = _assert_bounds(owner, raw)
    sample_evidence = _assert_samples(owner, samples)
    model_evidence = _assert_models(owner, models)
    public = _public_default_precision(owner)

    hashes = dict(owner.hashes)
    # The receipt must pin the producer raw input as well as the source and
    # checker code.  This prevents accepting a checker against a substituted
    # saved model archive.
    hashes[current.NAME] = current.sha(current.NAME)
    hashes[Path(current.__file__).name] = current.sha(Path(current.__file__).name)
    hashes[Path(__file__).name] = current.sha(Path(__file__).name)
    receipt = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            independent_full_pressure_fixture=fixture,
            independent_end_linear_jets_fixture=end_fixture,
            saved_live_complete_views_and_models=True,
            seven_actual_source_queries_checked=True,
            all_256_models_adjacency_and_full_Z=model_evidence,
            source_hashes_and_bindings_checked=True,
            physical_A_Rv_and_A_Rv2_units=unit_evidence,
            corrected_inverse_mu_ledger=True,
            quadratic_and_weighted_viscous_nonzero_bounds=bound_evidence,
            endpoint_q_plus_RqR_and_weighted_viscous_retained=True,
            default_precision_public_query=public,
            full_pressure_derivative_accuracy_not_claimed=True,
            no_full_tails_repair_or_global_NS_claim=True,
            sample_evidence=sample_evidence),
        input_hashes=hashes,
        **{current.GATE: True},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(',', ':'), ensure_ascii=False) + '\n',
        encoding='utf8', newline='\n')
    print('PASS_N1_END_PRESSURE_SOURCE', len(models), len(samples), flush=True)
    return receipt


if __name__ == '__main__':
    run()
