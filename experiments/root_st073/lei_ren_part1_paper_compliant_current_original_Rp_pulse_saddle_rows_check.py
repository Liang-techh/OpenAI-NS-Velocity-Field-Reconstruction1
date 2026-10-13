"""Independent checks for the finite-mu original Rp saddle-row refinement.

The producer keeps the accepted pulse source and only narrows its two positive
scaled full-row intervals.  This checker re-derives the saddle identity and
the finite-mu comparison on an independent moderate-mu fixture, then checks
every recorded directed interval and the source/tail provenance.
"""
import copy
import gzip
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_rows as current


HERE = Path(__file__).parent
NAME = current.NAME
RECEIPT = current.RECEIPT
SOURCE_DPS = current.SOURCE_DPS
ROW_NAMES = current.ROW_NAMES


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _read_raw():
    return json.loads(gzip.decompress((HERE / NAME).read_bytes()))


def _check_hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if _sha(name) != digest:
            raise ValueError('Refined saddle input changed: ' + name)


def _endpoint_bounds(value):
    """Return exact binary endpoints as ordinary mpf values."""
    if isinstance(value, dict):
        lo = tuple(value['lower_exact_mpf_tuple'])
        hi = tuple(value['upper_exact_mpf_tuple'])
    else:
        lo, hi = current.ep(value)
    return mp.make_mpf(tuple(lo)), mp.make_mpf(tuple(hi))


def _record_bounds(value):
    lo = mp.make_mpf(tuple(value['lower_exact_mpf_tuple']))
    hi = mp.make_mpf(tuple(value['upper_exact_mpf_tuple']))
    return lo, hi


def _relative_width(value):
    lo, hi = _endpoint_bounds(value)
    if lo <= 0:
        raise ArithmeticError('Positive interval required for width check')
    return (hi - lo) / lo


def _positive(value):
    lo, _ = _endpoint_bounds(value)
    return lo > 0


def _contains_interval(outer, inner):
    olo, ohi = _record_bounds(outer)
    ilo, ihi = _record_bounds(inner)
    return olo <= ilo and ihi <= ohi


def _overlap(left, right):
    llo, lhi = _record_bounds(left)
    rlo, rhi = _record_bounds(right)
    return max(llo, rlo) <= min(lhi, rhi)


def symbolic_kernel_fixture():
    x, w = sy.symbols('x w', real=True)
    v = 1 + x * w
    phi = x**2 * (2 * v + 1) / (6 * v**2)
    phase_residual = sy.factor(
        phi - x**2 / 2 + x**3 * w * (3 * v + 1) / (6 * v**2))
    if phase_residual != 0:
        raise ArithmeticError('Saddle phase identity failed')

    u, row = sy.symbols('u row', positive=True)
    denominator = 1 + sy.exp(-1 / u**2 + (1 - u)**(-2))
    kernel = ((sy.Rational(1099, 100) - u) *
              sy.exp((1 - u)**(-2) - phi + row * (2 + u)) /
              denominator)
    natural = sy.Rational(1099, 100) * sy.exp(1 + 2 * row - x**2 / 2)
    exact_ratio = ((1 - u / sy.Rational(1099, 100)) *
                   sy.exp((1 - u)**(-2) - 1 + row * u -
                          (phi - x**2 / 2)) / denominator)
    ratio_residual = sy.simplify(kernel / natural - exact_ratio)
    if ratio_residual != 0:
        raise ArithmeticError('Exact kernel ratio identity failed')
    return dict(phase_residual=str(phase_residual),
                ratio_residual=str(ratio_residual),
                exact_ratio='(1-u/10.99)*exp((1-u)^(-2)-1+row*u-(phi-x^2/2))/(1+exp(-1/u^2+(1-u)^(-2)))',
                ratio_factor='(1-u/10.99)*exp((1-u)^(-2)-1+row*u)/(1+exp(-1/u^2+(1-u)^(-2)))')


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
    """Compare direct quadrature to the natural Gaussian for mu=10^-24."""
    with mp.workdps(180):
        mu = mp.mpf('1e-24')
        band = mp.mpf(48)
        u0 = (4 * mu) ** (mp.mpf(1) / 3)
        w = u0 / mp.sqrt(6)
        u_hi = u0 * (1 + band * w)
        if not (u_hi < mp.mpf('.25') and abs(band * w) < mp.mpf('.25')):
            raise ArithmeticError('Moderate-mu saddle guards failed')
        rows = {}
        for row in (1, 2):
            actual = mp.quad(lambda xx: _direct_kernel(mu, row, xx),
                             [-band, -24, 0, 24, band])
            reference = (mp.mpf('10.99') * mp.exp(1 + 2 * row) *
                         mp.sqrt(2 * mp.pi))
            ratio = actual / reference
            log_ratio = mp.log(ratio)
            error = ((mp.mpf(2) / mp.mpf('10.99') + 4 + row) * u_hi
                     + 2 * band**3 * w + mp.exp(-1000))
            if not (mp.exp(-error) <= ratio <= mp.exp(error)):
                raise ArithmeticError('Direct moderate-mu ratio escaped bound')
            rows[str(row)] = dict(actual=mp.nstr(actual, 80),
                                  reference=mp.nstr(reference, 80),
                                  ratio=mp.nstr(ratio, 80),
                                  log_ratio=mp.nstr(log_ratio, 80),
                                  error=mp.nstr(error, 80),
                                  ratio_in_exp_error=True)
        return dict(mu='1e-24', band='48', u0=mp.nstr(u0, 80),
                    w=mp.nstr(w, 80), rows=rows)


def _check_row_diagnostics(view):
    rows = view['row_diagnostics']
    old = view['old_scaled_full_rows']
    refined = view['refined_scaled_full_rows']
    if len(rows) != 2 or len(old) != 2 or len(refined) != 2:
        raise ValueError('Exactly two saddle rows are required')
    if tuple(view['refined_rows_by_name']) != ROW_NAMES:
        raise ValueError('Unexpected row names')
    # Recompute the directed saddle guards, logarithmic error allowance, and
    # convex tangent tails from the recorded source parameters.  This avoids
    # treating producer booleans or copied diagnostics as the proof.
    c = MPIntervalContext()
    c.dps = SOURCE_DPS
    mu = c.make_mpf((tuple(view['original_mu']['lower_exact_mpf_tuple']),
                     tuple(view['original_mu']['upper_exact_mpf_tuple'])))
    u0 = c.make_mpf((tuple(view['original_saddle_u0']['lower_exact_mpf_tuple']),
                     tuple(view['original_saddle_u0']['upper_exact_mpf_tuple'])))
    w = c.make_mpf((tuple(view['original_saddle_w']['lower_exact_mpf_tuple']),
                    tuple(view['original_saddle_w']['upper_exact_mpf_tuple'])))
    band = int(view['band'])
    if current.ep(u0 * (1 + band * w))[1] > mp.mpf('.25'):
        raise ArithmeticError('Recorded saddle u guard is false')
    if abs(band * current.ep(w)[1]) > mp.mpf('.25'):
        raise ArithmeticError('Recorded saddle w guard is false')
    band_u = u0 * (1 + band * w)
    flat_log = -1 / band_u**2 + 1 / (1 - band_u)**2
    if current.ep(flat_log)[1] >= -1000:
        raise ArithmeticError('Independent flat-denominator cap is not below -1000')
    widths = []
    for index, diag in enumerate(rows):
        if diag['name'] != ROW_NAMES[index]:
            raise ValueError('Row order changed')
        if not _positive(refined[index]) or not _positive(old[index]):
            raise ArithmeticError('Pulse rows must remain positive')
        if not _overlap(old[index], refined[index]):
            raise ArithmeticError('Refinement lost overlap with accepted row')
        if _relative_width(refined[index]) >= _relative_width(old[index]):
            raise ArithmeticError('Refinement did not narrow row ' + str(index + 1))
        if not _contains_interval(refined[index], diag['core']):
            raise ArithmeticError('Core interval is not inside refined row')
        if not _positive(diag['old_tail']) or not _positive(diag['tangent_tail']):
            raise ArithmeticError('Positive tail evidence was lost')
        tail_lo, tail_hi = _record_bounds(diag['tail_bound'])
        old_tail_lo, old_tail_hi = _record_bounds(diag['old_tail'])
        support_lo, support_hi = _record_bounds(diag['support_upper'])
        # tail_bound is a directed upper-bound envelope produced by taking
        # endpoint maxima; its lower endpoint need not contain either input.
        if tail_hi < old_tail_hi:
            raise ArithmeticError('Tail bound does not retain old tail')
        if tail_hi < support_hi:
            raise ArithmeticError('Tail bound does not contain tangent support')
        late_lo, late_hi = _record_bounds(diag['late_log'])
        if late_hi >= -1000:
            raise ArithmeticError('Late support cap is not below -1000')
        err_lo, err_hi = _record_bounds(diag['log_ratio_error'])
        allow_lo, _ = _record_bounds(diag['log_ratio_allowance'])
        if not (err_lo > 0 and err_hi < allow_lo):
            raise ArithmeticError('Finite-mu error is not strictly below allowance')
        u_hi = current.ep(u0 * (1 + band * w))[1]
        w_hi = current.ep(w)[1]
        independent_error = ((mp.mpf(2) / mp.mpf('10.99') + 4 + index + 1) *
                             u_hi + 2 * mp.mpf(band)**3 * w_hi +
                             mp.exp(-1000))
        if independent_error > err_hi:
            raise ArithmeticError('Recorded error misses independently recomputed bound')
        flat_record_lo, flat_record_hi = _record_bounds(diag['flat_denominator_log'])
        flat_independent_hi = current.ep(flat_log)[1]
        if flat_record_hi >= -1000 or flat_independent_hi > flat_record_hi:
            raise ArithmeticError('Recorded flat-denominator cap misses independent bound')
        mu_lo, mu_hi = current.ep(mu)
        u0_lo, _ = current.ep(u0)
        independent_late_upper = (mp.log(121) - mp.mpf(1) / (4 * mu_hi) +
                                  mp.mpf(3) / u0_lo**2 +
                                  13 * (index + 1) - 2 * mp.log(u0_lo) +
                                  mp.log(6) / 2)
        late_lo, late_hi = _record_bounds(diag['late_log'])
        if independent_late_upper > late_hi or independent_late_upper >= -1000:
            raise ArithmeticError('Recorded late support cap misses independent bound')
        # Recompute the tangent envelope in a fresh interval context.  The
        # recorded support upper must enclose its upper endpoint plus cap.
        tangent = c.mpf(0)
        for sign in (-1, 1):
            x = c.mpf(sign * band)
            v = 1 + x * w
            phi = x**2 * (2 * v + 1) / (6 * v**2)
            derivative = x * (v**2 + v + 1) / (3 * v**3)
            dlo, dhi = current.ep(derivative)
            if dlo * dhi > 0:
                absolute = c.mpf([min(abs(dlo), abs(dhi)), max(abs(dlo), abs(dhi))])
            else:
                raise ArithmeticError('Independent tangent derivative crossed zero')
            tangent += c.exp(-phi) / absolute
        tangent *= 11 * c.exp(4 + c.mpf('2.5') * (index + 1))
        tangent_hi = current.ep(tangent)[1]
        support_hi = _record_bounds(diag['support_upper'])[1]
        cap_hi = mp.exp(-1000)
        if tangent_hi + cap_hi > support_hi:
            raise ArithmeticError('Support upper misses independent tangent tail')
        widths.append(dict(name=diag['name'],
                           old_relative_width=mp.nstr(_relative_width(old[index]), 30),
                           refined_relative_width=mp.nstr(_relative_width(refined[index]), 30),
                           width_ratio=mp.nstr(_relative_width(refined[index]) /
                                                _relative_width(old[index]), 30)))
    return widths


def _check_source_scope(raw, view, owner):
    if raw.get('candidate_original_Rp_pulse_saddle_rows_constructed') is not True:
        raise ValueError('Producer candidate record missing')
    if raw.get(current.GATE) is not False:
        raise ValueError('Producer gate must remain candidate-only before receipt')
    if view['source_family'] != owner.family:
        raise ValueError('Source-family mismatch')
    if view['source_binding']['selected_branch_unchanged'] is not True:
        raise ValueError('Original selection branch changed')
    if view['original_common_logpref_unchanged'] is not True:
        raise ValueError('Common log pref was not preserved')
    if view['original_selection_branch_unchanged'] is not True:
        raise ValueError('Selection branch provenance missing')
    if view['no_radius_or_N_materialization'] is not True:
        raise ValueError('Forbidden radius/N materialization claim')
    if any(view.get(key) for key in current.OPEN):
        raise ValueError('Out-of-scope gate was enabled')
    if not view['original_convex_tangent_tail_reused']:
        raise ValueError('Original tangent tail not retained')
    if not view['original_late_support_cap_reused']:
        raise ValueError('Original late support cap not retained')
    if view['guards']['all_rows_error_below_allowance'] is not True:
        raise ArithmeticError('Producer did not certify both row errors')
    if view['guards']['nonzero_mu'] is not True:
        raise ArithmeticError('Finite original mu was lost')
    if view['guards']['nonzero_u0'] is not True or view['guards']['nonzero_w'] is not True:
        raise ArithmeticError('Finite saddle parameters were lost')
    if view['phase_identity'] != 'phi-x^2/2=-x^3*w*(3*v+1)/(6*v^2)':
        raise ValueError('Unexpected phase identity')
    return dict(source_family_match=True, original_source_rows_bound=True,
                common_logpref_preserved=True, no_radius_or_N_materialization=True,
                scope_open_flags_false=True)


def _disjoint_intersection_fixture():
    """Ensure a disjoint old/new row cannot be fabricated into an interval."""
    c = MPIntervalContext()
    c.dps = 120
    mu = c.mpf('1e-300')
    u0 = c.exp(c.ln(4 * mu) / 3)
    w = u0 / c.sqrt(6)
    try:
        current._refine_row(c, mu, u0, w, c.mpf(0), c.mpf([0, 1]),
                            c.mpf(0), 1, 48, c.mpf('1e-60'), c.exp(-1000))
    except ArithmeticError as exc:
        if 'does not intersect' not in str(exc):
            raise
        return True
    raise ArithmeticError('Disjoint row fixture was accepted')


def run(owner=None, view=None):
    with mp.workdps(SOURCE_DPS):
        if owner is None:
            owner = current.CurrentOriginalRpPulseSaddleRows(require_checked=False)
        if view is None:
            view = owner.evaluate()
        raw = _read_raw()
        _check_hashes(raw)
        symbolic = symbolic_kernel_fixture()
        moderate = moderate_mu_fixture()
        widths = _check_row_diagnostics(view)
        scope = _check_source_scope(raw, view, owner)
        disjoint = _disjoint_intersection_fixture()
        # The public constructor must also work when the caller has left the
        # process at its ordinary default precision.  The source decorator
        # supplies its own 540-digit context.
        saved_dps = mp.mp.dps
        with mp.workdps(15):
            public = current.CurrentOriginalRpPulseSaddleRows(require_checked=False)
            public_view = public.evaluate()
            public.report(public_view)
            if mp.mp.dps != 15:
                raise ArithmeticError('Public query leaked inner working precision')
            public_widths = _check_row_diagnostics(public_view)
            if [item['name'] for item in public_widths] != [item['name'] for item in widths]:
                raise ArithmeticError('Default-precision public row set differs')
        if mp.mp.dps != saved_dps:
            raise ArithmeticError('Public query leaked outer precision')
        own_hashes = dict(owner.hashes)
        own_hashes[NAME] = _sha(NAME)
        own_hashes[Path(__file__).name] = _sha(Path(__file__).name)
        own_hashes[current.OUTER_RAW] = _sha(current.OUTER_RAW)
        receipt = dict(
            all_passed=True,
            source_family=owner.family,
            checks=dict(symbolic_kernel_identity=True,
                        moderate_mu_direct_quadrature=True,
                        source_scope=scope,
                        every_row_interval_checked=True,
                        every_row_width_narrowed=True,
                        finite_mu_error_strictly_below_allowance=True,
                        positive_tangent_and_late_support_retained=True,
                        disjoint_intersection_rejected=disjoint,
                        default_precision_public_query=True,
                        no_radius_or_N_materialization=True),
            row_widths=widths,
            symbolic=symbolic,
            moderate_mu_fixture=moderate,
            original_mu=view['original_mu'],
            common_logpref=view['original_common_logpref'],
            refined_scaled_full_rows=view['refined_scaled_full_rows'],
            original_selection_branch_unchanged=True,
            positive_unresolved_tail_provenance_retained=True,
            finite_mu_saddle_error_nonzero=True,
            **{current.GATE: True},
            **dict.fromkeys(current.OPEN, False),
            input_hashes=own_hashes,
        )
        data = json.dumps(receipt, separators=(',', ':'), ensure_ascii=False) + '\n'
        (HERE / RECEIPT).write_text(data, encoding='utf-8', newline='\n')
        # Reload through the checked public API.  This catches receipt/hash
        # cycles and exercises the real constructor acceptance guard.
        checked = current.CurrentOriginalRpPulseSaddleRows(require_checked=True)
        checked_view = checked.evaluate()
        _check_row_diagnostics(checked_view)
        if checked_view['source_family'] != owner.family:
            raise ValueError('Checked public family changed')
        print('CURRENT_RP_PULSE_SADDLE_ROWS_CHECK PASS', flush=True)
        return receipt


if __name__ == '__main__':
    run()
