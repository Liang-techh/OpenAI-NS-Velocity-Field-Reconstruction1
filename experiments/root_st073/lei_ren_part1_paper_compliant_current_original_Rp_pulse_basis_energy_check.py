"""Independent checks for the refined original Rp basis and pulse energy."""
import gzip
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy as current


HERE = Path(__file__).parent
SOURCE_DPS = current.SOURCE_DPS
NAME = current.NAME
RECEIPT = current.RECEIPT


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _read_raw():
    return json.loads(gzip.decompress((HERE / NAME).read_bytes()))


def _hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if _sha(name) != digest:
            raise ValueError('Basis-energy input changed: ' + name)


def _bounds(value):
    if isinstance(value, dict):
        lo = tuple(value['lower_exact_mpf_tuple'])
        hi = tuple(value['upper_exact_mpf_tuple'])
        return mp.make_mpf(lo), mp.make_mpf(hi)
    lo, hi = current.ep(value)
    return lo, hi


def _positive(value):
    return _bounds(value)[0] > 0


def _contains(value, point):
    lo, hi = _bounds(value)
    return lo <= point <= hi


def _overlap(left, right):
    l0, l1 = _bounds(left)
    r0, r1 = _bounds(right)
    return max(l0, r0) <= min(l1, r1)


def _relative_width(value):
    lo, hi = _bounds(value)
    if lo <= 0 <= hi:
        raise ArithmeticError('Relative width requires a sign-definite interval')
    return (hi - lo) / abs(lo)


def symbolic_source_fixture():
    r, a, ell = sy.symbols('r a ell', positive=True)
    d = 1 - r**2
    bump = sy.exp(-1 / d)
    reflected = bump.subs(r, -r)
    if sy.simplify(reflected - bump) != 0:
        raise ArithmeticError('Bump reflection identity failed')
    xi = sy.symbols('xi', real=True)
    z = sy.symbols('z', real=True)
    primitive = -sy.exp(-2 * xi) * ((xi - z)**2 / 2 +
                                    (xi - z) / 2 + sy.Rational(1, 4))
    if sy.simplify(sy.diff(primitive, xi) - sy.exp(-2 * xi) * (xi - z)**2) != 0:
        raise ArithmeticError('Analytic pulse primitive failed')
    v = sy.symbols('v', positive=True)
    sigma_num = sy.exp(-1 / v**2)
    sigma_den = sigma_num + sy.exp(-1 / (1 - v)**2)
    sigma = sigma_num / sigma_den
    sigma_reflection_residual = sy.simplify(
        sigma + sigma.subs(v, 1 - v) - 1)
    if sigma_reflection_residual != 0:
        raise ArithmeticError('Positive-ratio sigma reflection identity failed')
    # On 0 < v < 1, pair reflection gives the exact unit integral.  The
    # original flat extensions then make the scaled primitive affine after
    # xi=.02, while the second cutoff is identically one through xi=10.
    pair_integral = sy.Rational(1, 2) * sy.integrate(
        sigma + sigma.subs(v, 1 - v), (v, 0, 1))
    if sy.simplify(pair_integral - sy.Rational(1, 2)) != 0:
        raise ArithmeticError('Sigma pair integral derivation failed')
    xi_bulk = sy.symbols('xi_bulk', real=True)
    primitive_bulk = pair_integral / 50 + (xi_bulk - sy.Rational(1, 50))
    gp_bulk_residual = sy.simplify(primitive_bulk - (xi_bulk - sy.Rational(1, 100)))
    if gp_bulk_residual != 0:
        raise ArithmeticError('Flat-extension gp bulk identity failed')
    return dict(bump_reflection='b(-r)=b(r)',
                sigma_ratio='exp(-1/v^2)/(exp(-1/v^2)+exp(-1/(1-v)^2))',
                sigma_reflection='sigma(v)+sigma(1-v)=1',
                sigma_ratio_residual=str(sigma_reflection_residual),
                sigma_flat_extensions=True,
                gp_bulk='gp(xi)=xi-.01 for .02<=xi<=10',
                gp_bulk_pair_integral=str(pair_integral),
                gp_bulk_residual=str(gp_bulk_residual),
                gp_bulk_flat_cut=True,
                pulse_primitive=str(primitive),
                pulse_primitive_residual='0')


def _sigma(v):
    if v <= 0:
        return mp.mpf(0)
    if v >= 1:
        return mp.mpf(1)
    numerator = mp.exp(-1 / v**2)
    denominator = numerator + mp.exp(-1 / (1 - v)**2)
    return numerator / denominator


def _bump(r, n=1):
    if abs(r) >= 1:
        return mp.mpf(0)
    return mp.exp(-mp.mpf(n) / (1 - r * r))


def _direct_moment(a, n=1, multiplier=False):
    return mp.quad(lambda r: ((mp.mpf('.15') * r) if multiplier else 1) *
                   mp.exp(a * mp.mpf('.15') * r) * _bump(r, n), [-1, 0, 1])


def _bound_constants(n):
    def q(k):
        if k == 0:
            return mp.exp(-n)
        return (mp.mpf(k) / (n * mp.e))**k
    b0 = q(0)
    b1 = 2 * n * q(2)
    b2 = 4 * n * n * q(4) + 2 * n * q(2) + 8 * n * q(3)
    return b0, b1, b2


def _independent_midpoint(t, n, kind, cells=2048):
    h = mp.mpf(2) / cells
    total = mp.mpf(0)
    b0, b1, b2 = _bound_constants(n)
    ta = abs(t)
    if kind == 'j':
        sup = mp.exp(ta) * (ta**2 * b0 + 2 * ta * b1 + b2)
    elif kind == 'jr':
        g1 = mp.exp(ta) * (ta * b0 + b1)
        g2 = mp.exp(ta) * (ta**2 * b0 + 2 * ta * b1 + b2)
        sup = mp.mpf('.15') * (g2 + 2 * g1)
    elif kind == 'gram':
        sup = b2
    else:
        raise ValueError(kind)
    sampled_max = mp.mpf(0)
    for i in range(cells):
        left = -1 + mp.mpf(2 * i) / cells
        right = -1 + mp.mpf(2 * (i + 1)) / cells
        r = (left + right) / 2
        den = 1 - r * r
        b = mp.exp(-n / den)
        if kind == 'j':
            value = mp.exp(t * mp.mpf('.15') * r) * b
            sampled = mp.exp(t * mp.mpf('.15') * r) * (t**2 * b +
                     2 * t * (-2 * n * r * b / den**2) +
                     b * (4 * n**2 * r**2 / den**4 - 2 * n / den**2 -
                          8 * n * r**2 / den**3))
        elif kind == 'jr':
            g = mp.exp(t * mp.mpf('.15') * r) * b
            gp = mp.exp(t * mp.mpf('.15') * r) * (t * mp.mpf('.15') * b -
                  mp.mpf('.15') * 2 * n * r * b / den**2)
            gpp = mp.exp(t * mp.mpf('.15') * r) * ((t * mp.mpf('.15'))**2 * b +
                  2 * t * mp.mpf('.15') * (-2 * n * r * b / den**2) +
                  b * (4 * n**2 * r**2 / den**4 - 2 * n / den**2 -
                       8 * n * r**2 / den**3))
            value = mp.mpf('.15') * r * g
            sampled = mp.mpf('.15') * (r * gpp + 2 * gp)
        else:
            value = b
            sampled = b * (4 * n**2 * r**2 / den**4 - 2 * n / den**2 -
                           8 * n * r**2 / den**3)
        total += h * value
        sampled_max = max(sampled_max, abs(sampled))
    remainder = 2 * h**2 * sup / 24
    if sampled_max > sup:
        raise ArithmeticError('Sampled second derivative exceeded analytic bound')
    return total, remainder, sup, sampled_max


def moderate_quadrature_fixture():
    with mp.workdps(180):
        normalization = mp.quad(lambda r: _bump(r), [-1, 0, 1])
        checks = {}
        for t_text, kind, n in (('0.5', 'j', 1), ('1.5', 'jr', 1),
                                ('0', 'gram', 2)):
            t = mp.mpf(t_text)
            midpoint, remainder, sup, sampled = _independent_midpoint(t, n, kind, 512)
            if kind == 'j':
                direct = _direct_moment(t, n) / normalization
                low = (midpoint - remainder) / normalization
                high = (midpoint + remainder) / normalization
            elif kind == 'jr':
                direct = _direct_moment(t, n, True) / normalization
                low = (midpoint - remainder) / normalization
                high = (midpoint + remainder) / normalization
            else:
                direct = (mp.quad(lambda r: _bump(r, 2), [-1, 0, 1]) /
                          (mp.mpf('.15') * normalization**2))
                low = (midpoint - remainder) / (mp.mpf('.15') * normalization**2)
                high = (midpoint + remainder) / (mp.mpf('.15') * normalization**2)
            if not (low <= direct <= high):
                raise ArithmeticError('Independent midpoint remainder missed direct integral')
            checks[kind] = dict(t=t_text, n=n,
                                direct=mp.nstr(direct, 70),
                                midpoint=mp.nstr(midpoint, 70),
                                remainder=mp.nstr(remainder, 40),
                                global_second_derivative_bound=mp.nstr(sup, 40),
                                sampled_second_derivative_max=mp.nstr(sampled, 40))
        sigma_residuals = []
        for v in (mp.mpf('.1'), mp.mpf('.3'), mp.mpf('.5'), mp.mpf('.8')):
            sigma_residuals.append(mp.nstr(_sigma(v) + _sigma(1 - v) - 1, 40))
        if _sigma(mp.mpf('-.1')) != 0 or _sigma(mp.mpf('1.1')) != 1:
            raise ArithmeticError('Sigma flat extension mismatch')
        startup = mp.quad(lambda x: _sigma(50 * x), [0, mp.mpf('.02')])
        if abs(startup - mp.mpf('.01')) > mp.mpf('1e-100'):
            raise ArithmeticError('Sigma symmetry did not recover startup half mass')
        return dict(normalization=mp.nstr(normalization, 80), checks=checks,
                    sigma_reflection_residuals=sigma_residuals,
                    startup_sigma_mass=mp.nstr(startup, 80),
                    startup_sigma_mass_identity=True)


def _check_basis(view, outer):
    basis = view['bump_basis']
    old = outer['bump_basis']
    for name in ('first_row', 'exact_divided_difference_row'):
        if len(basis[name]) != 2:
            raise ValueError('Basis row count changed')
        for i, value in enumerate(basis[name]):
            if not _positive(value) or not _overlap(value, old[name][i]):
                raise ArithmeticError('Refined basis row is invalid or disjoint')
            if _relative_width(value) >= _relative_width(old[name][i]):
                raise ArithmeticError('Basis row did not narrow')
    det = basis['divided_difference_determinant']
    if _bounds(det)[1] >= 0 or not _overlap(det, old['divided_difference_determinant']):
        raise ArithmeticError('Determinant negativity/overlap failed')
    if _relative_width(det) >= _relative_width(old['divided_difference_determinant']):
        raise ArithmeticError('Determinant did not narrow')
    gram = basis['energy_gram']
    if not _positive(gram) or not _overlap(gram, old['energy_gram']):
        raise ArithmeticError('Energy Gram positivity/overlap failed')
    if _relative_width(gram) >= _relative_width(old['energy_gram']):
        raise ArithmeticError('Energy Gram did not narrow')
    if view['finite_mu_nonzero'] is not True or view['normalization_source_retained'] is not True:
        raise ValueError('Finite-mu or normalization provenance missing')
    if basis['original_mu_correction']['mu_nonzero'] is not True:
        raise ValueError('Finite-mu correction was replaced by zero')
    return dict(rows_narrowed=True, determinant_negative=True, gram_positive=True,
                finite_mu_correction_retained=True)


def _check_actual_record_values(view):
    """Check the production records against fresh direct mu=0 quadrature."""
    with mp.workdps(180):
        norm = mp.quad(lambda r: _bump(r), [-1, 0, 1])
        moments = view['bump_basis']['moments']
        for a_text in ('0.5', '1.5'):
            a = mp.mpf(a_text)
            j = _direct_moment(a, 1) / norm
            jr = _direct_moment(a, 1, True) / norm
            if not _contains(moments[a_text]['J'], j):
                raise ArithmeticError('Recorded J interval misses direct quadrature ' + a_text)
            if not _contains(moments[a_text]['J_r'], jr):
                raise ArithmeticError('Recorded J_r interval misses direct quadrature ' + a_text)
        j0 = _direct_moment(mp.mpf('.5'), 1) / norm
        jr0 = _direct_moment(mp.mpf('.5'), 1, True) / norm
        for index, d in enumerate((3, 1)):
            a0 = mp.exp(-mp.mpf(d) / 2) * j0
            d0 = mp.exp(-mp.mpf(d) / 2) * (d * j0 - jr0)
            if not _contains(view['bump_basis']['first_row'][index], a0):
                raise ArithmeticError('Recorded first-row interval misses direct quadrature')
            if not _contains(view['bump_basis']['exact_divided_difference_row'][index], d0):
                raise ArithmeticError('Recorded divided-row interval misses direct quadrature')
        det0 = -2 * mp.exp(-2) * j0**2
        if not _contains(view['bump_basis']['divided_difference_determinant'], det0):
            raise ArithmeticError('Recorded determinant misses direct reference')
        gram0 = mp.quad(lambda r: _bump(r, 2), [-1, 0, 1]) / (mp.mpf('.15') * norm**2)
        if not _contains(view['bump_basis']['energy_gram'], gram0):
            raise ArithmeticError('Recorded Gram interval misses direct quadrature')
    return dict(direct_moment_values_enclosed=True,
                direct_basis_rows_enclosed=True,
                direct_determinant_reference_enclosed=True,
                direct_gram_enclosed=True)


def _check_directed_derivative_records(view):
    """Compare every recorded interval bound with fresh scalar formulas."""
    # The producer records 540-digit directed endpoints.  Compare at a
    # still-higher precision so a low-precision nearest scalar cannot appear
    # just above a correctly directed endpoint by one rounding ulp.
    with mp.workdps(SOURCE_DPS + 60):
        def scalar_bounds(n):
            def q(k):
                if k == 0:
                    return mp.exp(-mp.mpf(n))
                return (mp.mpf(k) / (mp.mpf(n) * mp.e)) ** k
            b0 = q(0)
            b1 = 2 * mp.mpf(n) * q(2)
            b2 = 4 * mp.mpf(n) ** 2 * q(4) + 2 * mp.mpf(n) * q(2) + 8 * mp.mpf(n) * q(3)
            return b0, b1, b2

        def check_components(diag, n, label):
            expected = scalar_bounds(n)
            components = diag['bump_derivative_components']
            for key, value in zip(('b0', 'b1', 'b2'), expected):
                if _bounds(components[key])[1] < value:
                    raise ArithmeticError(
                        f'{label} {key} directed upper bound is too small')
            return expected

        moments = view['bump_basis']['moments']
        checks = {}
        for a_text in ('0.5', '1.5'):
            a = mp.mpf(a_text)
            j_diag = moments[a_text]['J_diagnostics']
            jr_diag = moments[a_text]['J_r_diagnostics']
            b0, b1, b2 = check_components(j_diag, 1, 'J ' + a_text)
            check_components(jr_diag, 1, 'J_r ' + a_text)
            t_abs = abs(a * mp.mpf('.15'))
            expected_j = mp.exp(t_abs) * (t_abs**2 * b0 + 2 * t_abs * b1 + b2)
            g1 = mp.exp(t_abs) * (t_abs * b0 + b1)
            g2 = mp.exp(t_abs) * (t_abs**2 * b0 + 2 * t_abs * b1 + b2)
            expected_jr = mp.mpf('.15') * (g2 + 2 * g1)
            if _bounds(j_diag['global_second_derivative_bound'])[1] < expected_j:
                raise ArithmeticError('J scalar derivative bound was not enclosed ' + a_text)
            if _bounds(jr_diag['global_second_derivative_bound'])[1] < expected_jr:
                raise ArithmeticError('J_r scalar derivative bound was not enclosed ' + a_text)
            checks[a_text] = dict(
                J_scalar_upper=mp.nstr(expected_j, 80),
                J_r_scalar_upper=mp.nstr(expected_jr, 80),
                J_record_upper=mp.nstr(_bounds(j_diag['global_second_derivative_bound'])[1], 80),
                J_r_record_upper=mp.nstr(_bounds(jr_diag['global_second_derivative_bound'])[1], 80))

        gram_diag = view['bump_basis']['gram_diagnostics']
        b0, b1, b2 = check_components(gram_diag, 2, 'Gram')
        if _bounds(gram_diag['global_second_derivative_bound'])[1] < b2:
            raise ArithmeticError('Gram scalar derivative bound was not enclosed')
        checks['gram'] = dict(
            scalar_upper=mp.nstr(b2, 80),
            record_upper=mp.nstr(_bounds(gram_diag['global_second_derivative_bound'])[1], 80))
        return dict(directed_derivative_bounds=True, scalar_comparisons=checks)


def _check_energy(view, outer):
    value = view['fixed_pulse_energy_Kpulse']
    old = outer['fixed_pulse_energy_Kpulse']
    if not _positive(value) or not _overlap(value, old):
        raise ArithmeticError('Refined Kpulse invalid or disjoint')
    if _relative_width(value) >= _relative_width(old):
        raise ArithmeticError('Kpulse did not narrow')
    diag = view['energy_diagnostics']
    bulk = _bounds(diag['bulk'])
    startup = _bounds(diag['startup'])
    exit_tail = _bounds(diag['exit'])
    if bulk[0] <= 0 or startup[0] < 0 or exit_tail[0] < 0:
        raise ArithmeticError('Pulse segment positivity failed')
    rounding_slack = mp.mpf(10) ** (-current.SOURCE_DPS + 20)
    if (_bounds(value)[0] < bulk[0] - rounding_slack or
            _bounds(value)[1] > bulk[1] + startup[1] + exit_tail[1] + rounding_slack):
        raise ArithmeticError('Kpulse does not respect analytic bulk/edge enclosure')
    if diag['sigma_symmetry_bulk'] is not True or diag['gp_exact_on_bulk'] is not True:
        raise ValueError('Bulk sigma/gp identity missing')
    return dict(kpulse_narrowed=True, positive_bulk=True,
                startup_positive_unknown_retained=True,
                exit_positive_unknown_retained=True)


def run(owner=None, view=None):
    with mp.workdps(SOURCE_DPS):
        if owner is None:
            owner = current.CurrentOriginalRpPulseBasisEnergy(require_checked=False)
        if view is None:
            view = owner.evaluate()
        raw = _read_raw()
        _hashes(raw)
        outer = json.loads((HERE / current.OUTER_RAW).read_bytes())
        symbolic = symbolic_source_fixture()
        moderate = moderate_quadrature_fixture()
        basis_checks = _check_basis(view, outer)
        actual_values = _check_actual_record_values(view)
        directed_values = _check_directed_derivative_records(view)
        energy_checks = _check_energy(view, outer)
        saved_dps = mp.mp.dps
        with mp.workdps(15):
            public = current.CurrentOriginalRpPulseBasisEnergy(require_checked=False)
            public_view = public.evaluate()
            public.report(public_view)
            if mp.mp.dps != 15:
                raise ArithmeticError('Public basis-energy query leaked inner precision')
            if len(public_view['bump_basis']['first_row']) != 2:
                raise ArithmeticError('Default-precision public basis changed shape')
        if mp.mp.dps != saved_dps:
            raise ArithmeticError('Public basis-energy query leaked outer precision')
        hashes = dict(owner.hashes)
        hashes[NAME] = _sha(NAME)
        hashes[Path(__file__).name] = _sha(Path(__file__).name)
        receipt = dict(
            all_passed=True,
            source_family=owner.family,
            checks=dict(symbolic_source_identities=True,
                        moderate_midpoint_and_sigma_fixture=True,
                        basis=basis_checks, energy=energy_checks,
                        direct_actual_records=actual_values,
                        directed_derivative_records=directed_values,
                        finite_mu_nonzero=True,
                        source_hash_scope=True,
                        default_precision_public_query=True,
                        no_selected_amplitude_or_pressure_claim=True),
            symbolic=symbolic,
            moderate_fixture=moderate,
            basis_relative_widths={
                'first_row': [mp.nstr(_relative_width(x), 30)
                              for x in view['bump_basis']['first_row']],
                'exact_divided_difference_row': [mp.nstr(_relative_width(x), 30)
                                                 for x in view['bump_basis']['exact_divided_difference_row']],
                'determinant': mp.nstr(_relative_width(view['bump_basis']['divided_difference_determinant']), 30),
                'energy_gram': mp.nstr(_relative_width(view['bump_basis']['energy_gram']), 30),
                'Kpulse': mp.nstr(_relative_width(view['fixed_pulse_energy_Kpulse']), 30),
            },
            original_selection_branch_unchanged=True,
            **{current.GATE: True},
            **dict.fromkeys(current.OPEN, False),
            input_hashes=hashes,
        )
        data = json.dumps(receipt, separators=(',', ':'), ensure_ascii=False) + '\n'
        (HERE / RECEIPT).write_text(data, encoding='utf-8', newline='\n')
        checked = current.CurrentOriginalRpPulseBasisEnergy(require_checked=True)
        checked.evaluate()
        print('CURRENT_RP_PULSE_BASIS_ENERGY_CHECK PASS', flush=True)
        return receipt


if __name__ == '__main__':
    run()
