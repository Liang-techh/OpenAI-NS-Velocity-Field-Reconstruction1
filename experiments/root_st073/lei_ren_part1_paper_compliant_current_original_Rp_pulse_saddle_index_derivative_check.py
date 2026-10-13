"""Independent checks for the original Rp saddle F1 index derivative."""

import gzip
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_index_derivative as current


HERE = Path(__file__).parent
NAME = current.NAME
RECEIPT = current.RECEIPT
SOURCE_DPS = current.SOURCE_DPS


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _read_raw():
    return json.loads(gzip.decompress((HERE / NAME).read_bytes()))


def _record_bounds(value):
    if isinstance(value, dict):
        return (mp.make_mpf(tuple(value['lower_exact_mpf_tuple'])),
                mp.make_mpf(tuple(value['upper_exact_mpf_tuple'])))
    return current.ep(value)


def _positive(value):
    return _record_bounds(value)[0] > 0


def symbolic_derivative_fixture():
    """Derive the index derivative and the positive u remainder independently."""
    i, u, x, w = sy.symbols('i u x w', real=True)
    v = 1 + x * w
    phi = x**2 * (2 * v + 1) / (6 * v**2)
    denominator = 1 + sy.exp(-1 / u**2 + (1 - u)**(-2))
    finite = ((sy.Rational(1099, 100) - u) *
              sy.exp((1 - u)**(-2) - phi + i * (2 + u)) /
              denominator)
    if sy.simplify(sy.diff(finite, i) - (2 + u) * finite) != 0:
        raise ArithmeticError('Index derivative identity failed')
    f1 = finite.subs(i, 1)
    if sy.simplify(sy.diff(finite, i).subs(i, 1) - (2 * f1 + u * f1)) != 0:
        raise ArithmeticError('F1prime decomposition failed')
    phase = sy.factor(phi - x**2 / 2 +
                      x**3 * w * (3 * v + 1) / (6 * v**2))
    if phase != 0:
        raise ArithmeticError('Saddle phase identity failed')
    return dict(index_derivative_residual='0',
                F1prime_decomposition_residual='0',
                phase_identity=str(phase),
                exact_definition='F1prime=2*F1+integral(u*kernel_1)',
                positive_kernel='kernel_1>0 on the accepted original support')


def _direct_kernel(mu, row, x):
    u0 = (4 * mu) ** (mp.mpf(1) / 3)
    w = u0 / mp.sqrt(6)
    v = 1 + x * w
    u = u0 * v
    phi = x**2 * (2 * v + 1) / (6 * v**2)
    return ((mp.mpf('10.99') - u) *
            mp.exp((1 - u)**(-2) - phi + row * (2 + u)) /
            (1 + mp.exp(-1 / u**2 + (1 - u)**(-2))))


def moderate_mu_fixture():
    """Directly quadrature the actual finite-mu core kernel and derivative."""
    with mp.workdps(90):
        mu = mp.mpf('1e-24')
        band = mp.mpf(48)
        pieces = [-band, -24, 0, 24, band]
        f1 = lambda x: _direct_kernel(mu, 1, x)

        def weighted(x):
            u0 = (4 * mu) ** (mp.mpf(1) / 3)
            w = u0 / mp.sqrt(6)
            return (2 + u0 * (1 + x * w)) * f1(x)

        actual = mp.quad(f1, pieces)
        weighted_u = mp.quad(lambda x: (
            (lambda u0, w: u0 * (1 + x * w))(
                (4 * mu) ** (mp.mpf(1) / 3),
                (4 * mu) ** (mp.mpf(1) / 3) / mp.sqrt(6))) * f1(x), pieces)
        derivative = mp.quad(weighted, pieces)
        residual = abs(derivative - (2 * actual + weighted_u))
        if residual >= mp.mpf('1e-50'):
            raise ArithmeticError('Direct finite-mu derivative residual too large')
        if actual <= 0 or weighted_u <= 0:
            raise ArithmeticError('Direct positive kernel integral lost positivity')
        return dict(mu='1e-24', band='48', F1=mp.nstr(actual, 60),
                    weighted_u=mp.nstr(weighted_u, 60),
                    F1prime=mp.nstr(derivative, 60),
                    identity_residual=mp.nstr(residual, 20),
                    direct_positive_integrals=True)


def _check_hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if _sha(name) != digest:
            raise ValueError('Index-derivative input changed: ' + name)


def _check_scope(raw, view, owner):
    if raw.get('candidate_original_Rp_pulse_saddle_index_derivative_constructed') is not True:
        raise ValueError('Candidate index-derivative record missing')
    if raw.get(current.GATE) is not False:
        raise ValueError('Producer gate must remain candidate-only')
    if view['source_family'] != owner.family:
        raise ValueError('Source-family mismatch')
    binding = view['source_binding']
    if binding['selected_branch_unchanged'] is not True:
        raise ValueError('Selection branch provenance missing')
    if binding['common_logpref_preserved'] is not True:
        raise ValueError('Common log pref provenance missing')
    if view['original_source_rows_unchanged'] is not True:
        raise ValueError('Original source rows were changed')
    if view['original_positive_tail_reused'] is not True:
        raise ValueError('Accepted positive tail was not reused')
    if view['no_radius_or_N_materialization'] is not True:
        raise ValueError('Forbidden radius/N materialization claim')
    if any(view.get(key, False) for key in current.OPEN):
        raise ValueError('Out-of-scope gate enabled')
    if view['common_logpref_index_independent'] is not True:
        raise ValueError('Index-independent common prefactor missing')
    return dict(source_family_match=True, selected_branch_unchanged=True,
                accepted_tail_reused=True, no_radius_or_N_materialization=True,
                scope_open_flags_false=True)


def _check_bounds(view, owner):
    c = owner.ctx
    row = owner.saddle.refined_rows[0]
    core = row['core']
    tail = row['tail_bound']
    u_lo = mp.mpf(_record_bounds(owner.bounds['weighted_core'])[0] /
                  _record_bounds(core)[0])
    u_hi = current.ep(owner.u0 * (1 + owner.band * owner.w))[1]
    core_lo, core_hi = current.ep(core)
    tail_hi = current.ep(tail)[1]
    expected_core = c.mpf([
        current.ep(owner.u0 * (1 - owner.band * owner.w))[0] * core_lo,
        u_hi * core_hi])
    expected_tail = c.mpf([0, c.mpf(11) * tail_hi])
    expected_total = c.mpf([
        current.ep(expected_core)[0] + current.ep(expected_tail)[0],
        current.ep(expected_core)[1] + current.ep(expected_tail)[1]])
    if _record_bounds(view['weighted_core']) != current.ep(expected_core):
        raise ArithmeticError('Weighted core bound was not independently reproduced')
    if _record_bounds(view['weighted_tail']) != current.ep(expected_tail):
        raise ArithmeticError('Weighted tail bound was not independently reproduced')
    if _record_bounds(view['weighted_u_kernel']) != current.ep(expected_total):
        raise ArithmeticError('Weighted u-kernel bound was not independently reproduced')
    f1 = owner.f1
    expected_derivative = f1 * 2 + expected_total
    if _record_bounds(view['F1prime']) != current.ep(expected_derivative):
        raise ArithmeticError('F1prime bound does not equal 2F1 plus weighted bound')
    if not (current.ep(owner.u0 * (1 - owner.band * owner.w))[0] > 0):
        raise ArithmeticError('Core u lower bound is not positive')
    if u_hi > mp.mpf('.25'):
        raise ArithmeticError('Core u upper guard failed')
    if abs(owner.band * current.ep(owner.w)[1]) > mp.mpf('.25'):
        raise ArithmeticError('Core x-to-u guard failed')
    if current.ep(expected_core)[0] <= 0:
        raise ArithmeticError('Positive weighted core remainder was lost')
    if current.ep(expected_total)[1] <= 0:
        raise ArithmeticError('Positive weighted remainder was lost')
    if current.ep(expected_tail)[1] <= 0:
        raise ArithmeticError('Accepted positive tail was lost')
    return dict(weighted_core_recomputed=True, weighted_tail_recomputed=True,
                F1prime_decomposition_recomputed=True,
                positive_core_and_tail_retained=True,
                core_u_upper_guard=True, support_u_upper='11')


def run(owner=None, view=None):
    with mp.workdps(SOURCE_DPS):
        if owner is None:
            owner = current.CurrentOriginalRpPulseSaddleIndexDerivative(
                require_checked=False)
        if view is None:
            view = owner.evaluate()
        raw = _read_raw()
        _check_hashes(raw)
        if raw['view'] != current.record(view):
            raise ValueError('Saved and live derivative views differ')
        symbolic = symbolic_derivative_fixture()
        moderate = moderate_mu_fixture()
        scope = _check_scope(raw, view, owner)
        bounds = _check_bounds(view, owner)
        saved_dps = mp.mp.dps
        with mp.workdps(15):
            public = current.CurrentOriginalRpPulseSaddleIndexDerivative(
                saddle_owner=owner.saddle, require_checked=False)
            public_view = public.evaluate()
            public.report(public_view)
            if mp.mp.dps != 15:
                raise ArithmeticError('Public query leaked source precision')
        if mp.mp.dps != saved_dps:
            raise ArithmeticError('Public query leaked outer precision')
        hashes = dict(owner.hashes)
        hashes[NAME] = _sha(NAME)
        hashes[Path(__file__).name] = _sha(Path(__file__).name)
        receipt = dict(
            all_passed=True,
            source_family=owner.family,
            checks=dict(symbolic_derivative_identity=True,
                        moderate_mu_direct_quadrature=True,
                        source_scope=scope,
                        weighted_core_and_tail_bounds=bounds,
                        raw_live_view_checked=True,
                        default_precision_public_query=True,
                        original_positive_source_preserved=True,
                        no_radius_or_N_materialization=True),
            symbolic=symbolic,
            moderate_mu_fixture=moderate,
            view=view,
            **{current.GATE: True},
            **dict.fromkeys(current.OPEN, False),
            input_hashes=hashes)
        (HERE / RECEIPT).write_text(
            json.dumps(receipt, separators=(',', ':'), ensure_ascii=False) + '\n',
            encoding='utf8', newline='\n')
        checked = current.CurrentOriginalRpPulseSaddleIndexDerivative(
            saddle_owner=owner.saddle, require_checked=True)
        checked_view = checked.evaluate()
        if checked_view['source_family'] != owner.family:
            raise ValueError('Checked derivative family changed')
        print('CURRENT_RP_PULSE_SADDLE_INDEX_DERIVATIVE_CHECK PASS', flush=True)
        return receipt


if __name__ == '__main__':
    run()
