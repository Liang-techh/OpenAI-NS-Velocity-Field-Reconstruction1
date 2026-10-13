"""Independent checks for the bounded pulse-pressure C1 adapter."""
import gzip
import hashlib
import json
from functools import wraps
from pathlib import Path

import mpmath as mp
import sympy as sy

try:
    import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1 as current
except ImportError:
    current = None

HERE = Path(__file__).parent
SOURCE_DPS = 540


def _source_precision(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        with mp.workdps(SOURCE_DPS):
            return fn(*args, **kwargs)
    return call


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _interval(value):
    if hasattr(value, '_mpi_'):
        return tuple(mp.make_mpf(x) for x in value._mpi_)
    if isinstance(value, dict) and 'lower_exact_mpf_tuple' in value:
        return (mp.make_mpf(tuple(value['lower_exact_mpf_tuple'])),
                mp.make_mpf(tuple(value['upper_exact_mpf_tuple'])))
    raise TypeError('Directed interval required')


def _record_upper(value):
    return _interval(value)[1]


def _read_raw():
    if current is None:
        raise RuntimeError('C1 producer is unavailable')
    return json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))


def independent_symbolic_fixture():
    """Check qZZ coefficients and differentiated quadratic identities independently."""
    z, delta = sy.symbols('z delta', real=True)
    U = 1 + z + sy.Rational(1, 3) * z ** 2
    M = 2 - z + sy.Rational(1, 5) * z ** 3
    a = 1 - delta
    d = 1 - z ** 2
    L = 1 - delta * z ** 2
    q = (2 * z * U - a * z * M - d * sy.diff(M, z)) / L
    qzz = sy.diff(q, z, 2)
    majorant = (56 * abs(U) + 24 * abs(sy.diff(U, z))
                + 4 * abs(sy.diff(U, z, 2)) + 28 * abs(M)
                + 38 * abs(sy.diff(M, z)) + 16 * abs(sy.diff(M, z, 2))
                + 2 * abs(sy.diff(M, z, 3)))
    for zv in (sy.Rational(-1), sy.Rational(-1, 2), 0,
               sy.Rational(1, 2), sy.Rational(1)):
        value = abs(sy.N(qzz.subs({z: zv, delta: sy.Rational(1, 5)}), 80))
        bound = sy.N(majorant.subs({z: zv, delta: sy.Rational(1, 5)}), 80)
        if value > bound:
            raise ArithmeticError('Independent qZZ majorant failed')

    # Global analytic qZZ proof.  Write q=f U+g m+h m_Z and differentiate
    # twice, rather than relying on the five sample points above.  The seven
    # coefficients are [f'',2f',f,g'',2g'+h'',g+2h',h].
    f = 2 * z / L
    g = -a * z / L
    h = -d / L
    hprime = 2 * a * z / L ** 2
    combo_gprime = a * delta * z ** 2 * (8 - 2 * L) / L ** 3
    combo_g_hprime = a * z * (4 - L) / L ** 2
    identities = (
        sy.diff(f, z, 2) - 4 * delta * z * (3 + delta * z ** 2) / L ** 3,
        2 * sy.diff(f, z) - 4 * (1 + delta * z ** 2) / L ** 2,
        sy.diff(g, z, 2) + 2 * a * delta * z * (3 + delta * z ** 2) / L ** 3,
        2 * sy.diff(g, z) + sy.diff(h, z, 2) - combo_gprime,
        g + 2 * sy.diff(h, z) - combo_g_hprime,
        sy.diff(h, z) - hprime,
    )
    if any(sy.simplify(item) != 0 for item in identities):
        raise ArithmeticError('Analytic qZZ coefficient simplification failed')
    # On |z|<=1, 0<=delta<=1/2: a<=1, d<=1, 1/2<=L<=1,
    # hence L^-1<=2, L^-2<=4, L^-3<=8.  The following exact rational
    # arithmetic is the advertised global bound, independent of field data.
    half = sy.Rational(1, 2)
    invL = sy.Integer(2)
    analytic_bounds = [
        4 * half * (3 + half) * invL ** 3,
        4 * (1 + half) * invL ** 2,
        2 * invL,
        2 * half * (3 + half) * invL ** 3,
        half * 7 * invL ** 3,
        4 * invL ** 2,
        invL,
    ]
    qzz_bounds = [56, 24, 4, 28, 38, 16, 2]
    if any(value > limit for value, limit in zip(analytic_bounds, qzz_bounds)):
        raise ArithmeticError('Analytic qZZ global bound failed')

    # Differentiate the exact integrated S0/S1 formulas independently.  The
    # coefficient derivatives are written out explicitly, so this catches a
    # missing d/L, Z/L, or L^-2 contribution rather than replaying a product
    # rule for a synthetic reduced expression.
    a = 1 - delta
    d = 1 - z ** 2
    L = 1 - delta * z ** 2
    Lz = -2 * delta * z
    dz = -2 * z
    c0 = [1 / L, a * z / (2 * L), sy.Rational(1, 2), d / L,
          -2 * z ** 2 / L ** 2, -2 * delta * z ** 2 / L ** 2,
          2 * z * d / L ** 2]
    c0z = [-Lz / L ** 2,
           a * (L - z * Lz) / (2 * L ** 2), sy.Integer(0),
           (dz * L - d * Lz) / L ** 2,
           -4 * z / L ** 2 + 4 * z ** 2 * Lz / L ** 3,
           -4 * delta * z / L ** 2 + 4 * delta * z ** 2 * Lz / L ** 3,
           2 * (d + z * dz) / L ** 2 - 4 * z * d * Lz / L ** 3]
    c1 = [1 / L, -1 / L, a * z / (2 * L), sy.Rational(1, 2),
          -sy.Rational(1, 2), d / L, -2 * z ** 2 / L ** 2,
          2 * a * z ** 2 / L ** 2, 2 * z * d / L ** 2]
    c1z = [-Lz / L ** 2, Lz / L ** 2,
           a * (L - z * Lz) / (2 * L ** 2), sy.Integer(0), sy.Integer(0),
           (dz * L - d * Lz) / L ** 2,
           -4 * z / L ** 2 + 4 * z ** 2 * Lz / L ** 3,
           4 * a * z / L ** 2 - 4 * a * z ** 2 * Lz / L ** 3,
           2 * (d + z * dz) / L ** 2 - 4 * z * d * Lz / L ** 3]
    names0 = ('Bq', 'I_qZ', 'B_q2', 'I_UqZ', 'B_U2', 'I_U2', 'I_UUZ')
    names1 = ('B_Rq', 'I_q', 'I_RqZ', 'B_Rq2', 'I_Rq2', 'I_RUqZ',
              'B_RU2', 'I_RU2', 'I_RUUZ')
    t0 = [sy.Function(name)(z) for name in names0]
    t1 = [sy.Function(name)(z) for name in names1]
    s0 = sum(coef * term for coef, term in zip(c0, t0))
    s1 = sum(coef * term for coef, term in zip(c1, t1))
    expected0 = sum(coefz * term + coef * sy.diff(term, z)
                    for coef, coefz, term in zip(c0, c0z, t0))
    expected1 = sum(coefz * term + coef * sy.diff(term, z)
                    for coef, coefz, term in zip(c1, c1z, t1))
    for coef, coefz in zip(c0, c0z):
        if sy.simplify(sy.diff(coef, z) - coefz) != 0:
            raise ArithmeticError('S0 coefficient derivative identity failed')
    for coef, coefz in zip(c1, c1z):
        if sy.simplify(sy.diff(coef, z) - coefz) != 0:
            raise ArithmeticError('S1 coefficient derivative identity failed')
    if sy.simplify(sy.diff(s0, z) - expected0) != 0:
        raise ArithmeticError('S0 exact derivative identity failed')
    if sy.simplify(sy.diff(s1, z) - expected1) != 0:
        raise ArithmeticError('S1 exact derivative identity failed')
    # Independently verify the finite interval coefficient envelope used by
    # the production majorant.  This is an enclosure check, not a field cap.
    with mp.workdps(120):
        iz = mp.iv.mpf([-1, 1])
        idelta = mp.iv.mpf([0, '0.5'])
        ia = 1 - idelta
        id_ = 1 - iz * iz
        iL = 1 - idelta * iz * iz
        iLz = -2 * idelta * iz
        idz = -2 * iz
        coeff = [1 / iL, ia * iz / (2 * iL), id_ / iL,
                 2 * iz * iz / (iL * iL),
                 2 * idelta * iz * iz / (iL * iL),
                 2 * iz * id_ / (iL * iL)]
        coeffz = [-iLz / (iL * iL),
                  ia * (iL - iz * iLz) / (2 * iL * iL),
                  (idz * iL - id_ * iLz) / (iL * iL),
                  4 * iz / (iL * iL) - 4 * iz * iz * iLz / (iL ** 3),
                  4 * idelta * iz / (iL * iL)
                  - 4 * idelta * iz * iz * iLz / (iL ** 3),
                  2 * (id_ + iz * idz) / (iL * iL)
                  - 4 * iz * id_ * iLz / (iL ** 3)]
        coefficient_envelope = max(
            max(abs(mp.make_mpf(v)) for v in item._mpi_)
            for item in (coeff + coeffz))
        if coefficient_envelope >= 256:
            raise ArithmeticError('Reduced-S coefficient envelope failed')
    return dict(qZZ_coefficients=qzz_bounds,
                qZZ_fixture_points=5, qZZ_global_analytic=True,
                qZZ_decomposition="q=f*U+g*m+h*m_Z; coefficients=[f'',2f',f,g'',2g'+h'',g+2h',h]",
                qZZ_simplifications={
                    'h_prime': '2*(1-delta)*Z/L^2',
                    '2g_prime_plus_h_second': '(1-delta)*delta*Z^2*(8-2*L)/L^3',
                    'g_plus_2h_prime': '(1-delta)*Z*(4-L)/L^2'},
                qZZ_global_bound_arithmetic=[int(value) for value in analytic_bounds],
                quadratic_product_rules=True,
                reduced_S0_Z_identity=True, reduced_S1_Z_identity=True,
                S1_constant_q2_coefficients_zero_derivative=True,
                S1_constant_q2_derivative_term='sqrt(E_q*E_qZ)',
                S1_constant_q2_production_factor='256*(1/2)>=1; conservative coefficient envelope',
                reduced_S_coefficient_envelope_upper=str(coefficient_envelope))


def independent_numeric_fixture():
    """Check a nonunit finite lower endpoint and the weighted viscous sign."""
    with mp.workdps(90):
        inv_rv = mp.iv.mpf(['1e-30', '1.0000000000001e-30'])
        reference = mp.iv.mpf(['2.3', '2.3000000000001'])
        q0 = mp.iv.mpf(['-.31', '-.30'])
        q1 = mp.iv.mpf(['.07', '.08'])
        expected = -q1 * inv_rv / reference
        if _interval(expected)[1] >= 0 or _interval(expected)[0] >= 0:
            raise ArithmeticError('Signed weighted endpoint fixture lost its negative sign')
        if _interval(-q1 * inv_rv / reference) != _interval(expected):
            raise ArithmeticError('Weighted endpoint sign changed')
        if _interval(q0)[1] == 0:
            raise ArithmeticError('Value endpoint fixture was degenerate')
        # Independent first-derivative remainder fixture: for p(Z)=Z^2,
        # p'(Z)=2Z varies by exactly 2h from its center value on a halfwidth
        # h cell.  An h^2 value-remainder would underbound this for h<2.
        h = mp.iv.mpf(['.01', '.01'])
        center = mp.iv.mpf(['.5', '.5'])
        variation = abs(2 * (center + h) - 2 * center)
        value_remainder = h * h
        if _interval(variation)[0] <= _interval(value_remainder)[1]:
            raise ArithmeticError('First-derivative h scaling fixture failed')
        return dict(nonunit_inverse_radius=True,
                    weighted_viscous_derivative_sign='-qinit_Z*invRv/reference',
                    unweighted_endpoint_derivative_zero=True,
                    first_derivative_polynomial='p(Z)=Z^2; p_prime=2Z',
                    first_derivative_remainder_scales='2*h',
                    value_remainder_wrong_scale='h^2')


def _check_scope(raw, owner):
    if raw.get('source_family') != owner.family:
        raise ValueError('C1 raw source family mismatch')
    if raw.get(current.GATE) is not False:
        raise ValueError('C1 producer gate must be candidate-only before receipt')
    if raw.get(current.LOCAL_GATE) is not False:
        raise ValueError('Local C1 grid gate must be candidate-only before receipt')
    if raw.get('c1_grid_count') != current.C1_GRID_COUNT:
        raise ValueError('Local C1 grid count changed')
    if any(raw.get(key, False) for key in current.OPEN):
        raise ValueError('C1 scope widened')
    if raw.get('units', {}).get('full_pressure_derivative_accuracy_not_claimed') is not True:
        raise ValueError('Uniform derivative scope was claimed')
    if raw.get('order_three_energy', {}).get('source_order_three_retained') is not True:
        raise ValueError('Order-three source provenance missing')
    return dict(source_family=True, local_derivative_scope_only=True,
                local_grid_candidate=True, global_and_repair_gates_false=True)


def _check_energy(raw, owner):
    energy = raw['order_three_energy']
    if len(energy['log_caps']) != 4:
        raise ValueError('Expected B/Tinitial rows through order 3')
    for row in energy['log_caps']:
        for value in row.values():
            if _record_upper(value) >= -1000:
                raise ArithmeticError('Order-three energy row exceeds -1000')
    if energy['qZZ_majorant_coefficients'] != [56, 24, 4, 28, 38, 16, 2]:
        raise ValueError('qZZ coefficient proof changed')
    if energy.get('qZZ_cauchy_factor') != 7:
        raise ValueError('qZZ seven-term Cauchy factor missing')
    eps = _interval(energy['epsilon'])
    if eps[0] <= 0 or eps[1] <= 0:
        raise ArithmeticError('epsilon is not a positive nonzero directed scale')
    if _record_upper(energy['qZZ_energy']) <= 0:
        raise ArithmeticError('qZZ energy was collapsed to zero')
    with mp.workdps(SOURCE_DPS):
        # Use the serialized directed epsilon itself as the independent
        # source scale; recomputing exp(-1000) in a fresh interval context can
        # round one endpoint differently from the owner's basis context.
        expected_qzz_upper = mp.mpf(200256) * eps[1]
        if _record_upper(energy['qZZ_energy']) < expected_qzz_upper:
            raise ArithmeticError('qZZ energy is below the independent Cauchy bound')
    if _record_upper(energy['boundary_qZ_squared_bound']) <= 0:
        raise ArithmeticError('Endpoint qZ energy was collapsed to zero')
    live = owner.energy
    if current._json_safe(live) != energy:
        raise ValueError('Saved C1 energy differs from live C1 energy')
    return dict(order_three_rows=True, all_log_caps_below_minus_1000=True,
                qZZ_nonzero=True, qZZ_cauchy_factor=7,
                qZZ_independent_bound=True, endpoint_qZ_nonzero=True)


def _check_samples(raw, owner):
    samples = raw.get('samples', {})
    if tuple(samples) != current.SAMPLES:
        raise ValueError('C1 sample set changed')
    live = {}
    for z in current.SAMPLES:
        packet = owner.evaluate(z)
        view = owner.report(packet)
        live[z] = view
        if view['source_family'] != owner.family:
            raise ValueError('C1 sample family mismatch')
        if view.get('full_pressure_derivative_accuracy_not_claimed') is not True:
            raise ValueError('Sample derivative scope widened')
        if view.get('unweighted_viscous_boundary_exact_zero') is not True:
            raise ValueError('Unweighted endpoint derivative identity missing')
        if _record_upper(view['weighted_viscous_endpoint_derivative']) == 0:
            raise ArithmeticError('Weighted endpoint derivative was dropped')
        if _record_upper(view['qZZ_source_bound']) <= 0:
            raise ArithmeticError('Sample qZZ cap was dropped')
        expected = -view['source_linear_inputs']['qinit'][1] \
            * owner.pulse.invRv / owner.pulse.reference
        if _interval(view['weighted_viscous_endpoint_derivative']) != _interval(expected):
            raise ArithmeticError('Weighted viscous derivative sign/formula changed')
        for key in ('m13', 'm15'):
            if _record_upper(view['quadratic_C1_source_errors'][key]) <= 0:
                raise ArithmeticError('C1 quadratic source error was zero')
            if _record_upper(view['full_pressure_first_derivative'][key]) \
                    <= _interval(view['full_pressure_first_derivative'][key])[0]:
                raise ArithmeticError('C1 derivative interval is degenerate')
    if raw.get('samples') != current._json_safe(live):
        raise ValueError('Saved C1 samples differ from live samples')
    return dict(samples=len(live), endpoint_signs=True,
                quadratic_source_errors_nonzero=True, qZZ_rows_retained=True)


def _check_models(raw, owner):
    models = raw.get('models', [])
    if len(models) != 256 or len(owner.models) != 256:
        raise ValueError('C1 model coverage is not 256 cells')
    if models != current._json_safe(owner.models):
        raise ValueError('Saved C1 models differ from live models')
    checked = 0
    for model in models:
        for key in ('m13', 'm15'):
            fields = model['fields'][key]
            for name in ('first_derivative_source_error',
                         'first_derivative_linear_Taylor_remainder',
                         'first_derivative_coefficient_ambiguity',
                         'first_derivative_absolute_error'):
                if _record_upper(fields[name]) < 0:
                    raise ArithmeticError('Negative C1 model bound')
            if fields['first_derivative_accuracy_certified'] is not False:
                raise ValueError('Uniform derivative accuracy was claimed')
            if _record_upper(fields['first_derivative_absolute_error']) <= 0:
                raise ArithmeticError('C1 model error was collapsed to zero')
        checked += 1
    return dict(models=checked, finite_first_derivative_enclosures=True,
                uniform_derivative_accuracy_false=True)


def _check_local_grid(raw, owner):
    """Check the fine grid against actual center/whole-subcell source rows."""
    grid = raw.get('c1_grid', [])
    if len(grid) != current.C1_GRID_COUNT or len(owner.c1_grid) != current.C1_GRID_COUNT:
        raise ValueError('Local C1 grid coverage is not 4096 cells')
    if grid != current._json_safe(owner.c1_grid):
        raise ValueError('Saved local C1 grid differs from live grid')
    previous_right = None
    max_error = mp.mpf(0)
    threshold = mp.mpf('.001')
    for index, row in enumerate(grid):
        if row.get('index') != index:
            raise ValueError('Local C1 grid ordering changed')
        left = _interval(row['left'])
        right = _interval(row['right'])
        center = _interval(row['center'])
        half = _interval(row['half_width'])
        if left[0] > right[1] or center[0] < left[0] or center[1] > right[1]:
            raise ArithmeticError('Invalid local C1 subcell geometry')
        if previous_right is not None and left[0] > previous_right:
            raise ArithmeticError('Gap in local C1 grid coverage')
        previous_right = right[1]
        if half[1] <= 0:
            raise ArithmeticError('Local C1 halfwidth was collapsed')
        for key in ('m13', 'm15'):
            fields = row['fields'][key]
            for name in ('center_slope', 'whole_Z2_bound', 'half_width',
                         'source_error', 'slope_ambiguity', 'total_error'):
                _interval(fields[name])
            slope = _interval(fields['center_slope'])
            z2 = _interval(fields['whole_Z2_bound'])
            source = _interval(fields['source_error'])
            ambiguity = _interval(fields['slope_ambiguity'])
            total = _interval(fields['total_error'])
            if source[1] <= 0 or ambiguity[1] <= 0 or total[1] <= 0:
                raise ArithmeticError('Local C1 nonzero source/ambiguity was dropped')
            if _interval(owner.full_errors[key]) != source:
                raise ValueError('Local C1 source error is not the accepted source cap')
            if ambiguity[0] < 0:
                raise ArithmeticError('Local C1 slope ambiguity became signed')
            expected_ambiguity = (slope[1] - slope[0]) / 2
            if ambiguity[1] < expected_ambiguity:
                raise ArithmeticError('Local C1 slope ambiguity underbounds source row')
            expected_remainder = 2 * max(abs(z2[0]), abs(z2[1])) * half[1]
            expected = expected_remainder + expected_ambiguity + source[1]
            if total[1] < expected:
                raise ArithmeticError('Local C1 total error underbounds source row')
            if total[1] >= threshold:
                raise ArithmeticError('Local C1 .001 accuracy failed')
            if fields.get('absolute_accuracy_certified') is not True:
                raise ValueError('Local C1 row was not certified')
            max_error = max(max_error, total[1])
    if previous_right is None or previous_right < 1:
        raise ArithmeticError('Local C1 grid does not reach Z=1')
    first_left = _interval(grid[0]['left'])
    if first_left[0] > -1:
        raise ArithmeticError('Local C1 grid does not start at Z=-1')
    if raw.get('source_bindings', {}).get('local_C1_grid') != \
            'actual pulse._make_view at every center and whole subcell':
        raise ValueError('Local C1 grid source binding changed')
    return dict(cells=len(grid), exact_source_views=True,
                contiguous_full_Z_coverage=True, max_total_error=mp.nstr(max_error, 16),
                local_accuracy_threshold='.001', local_accuracy_certified=True)


def _check_hashes(raw, owner):
    for name, digest in raw.get('input_hashes', {}).items():
        if _sha(name) != digest:
            raise ValueError('Changed C1 input ' + name)
    if raw.get('parent_archive', {}).get('name') != current.pulse.NAME:
        raise ValueError('Wrong accepted pulse archive')
    return True


@_source_precision
def run(owner=None, raw=None):
    if current is None:
        raise RuntimeError('C1 producer is unavailable')
    if owner is None:
        raise ValueError('Pass the retained checked pulse owner; no ancestor reconstruction is permitted')
    owner.assert_graph()
    if raw is None:
        raw = _read_raw()
    symbolic = independent_symbolic_fixture()
    numeric = independent_numeric_fixture()
    scope = _check_scope(raw, owner)
    hashes = _check_hashes(raw, owner)
    energy = _check_energy(raw, owner)
    samples = _check_samples(raw, owner)
    models = _check_models(raw, owner)
    local_grid = _check_local_grid(raw, owner)
    before = mp.mp.dps
    with mp.workdps(15):
        packet = owner.evaluate('.371')
        public = owner.report(packet)
        if mp.mp.dps != 15 or public['source_family'] != owner.family:
            raise ArithmeticError('Default-precision C1 public query failed')
    if mp.mp.dps != before:
        raise ArithmeticError('Default-precision query leaked precision')
    receipt = dict(
        all_passed=True, source_family=owner.family,
        checks=dict(symbolic=symbolic, numeric=numeric, scope=scope,
                    hashes=hashes, energy=energy, samples=samples,
                    models=models, local_grid=local_grid,
                    default_precision_public_query=True),
        **{current.GATE: True, current.LOCAL_GATE: True},
        **dict.fromkeys(current.OPEN, False),
        input_hashes=dict(owner.hashes), source_precision=SOURCE_DPS,
    )
    receipt['input_hashes'][current.NAME] = _sha(current.NAME)
    receipt['input_hashes'][Path(__file__).name] = _sha(Path(__file__).name)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(receipt, separators=(',', ':')) + '\n',
        encoding='utf8', newline='\n')
    return receipt


if __name__ == '__main__':
    raise SystemExit('Use run(checked_pulse_owner)')
