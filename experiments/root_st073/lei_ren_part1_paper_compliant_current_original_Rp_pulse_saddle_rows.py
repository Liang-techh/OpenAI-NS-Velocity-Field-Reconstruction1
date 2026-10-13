"""Analytic refinement of the accepted original Rp pulse saddle rows.

The original pulse kernel and its common logarithmic prefactor remain the
source.  This adapter only narrows the two positive scaled full-row
enclosures.  It binds the accepted outer-pulse source record (and, when
provided, a live CurrentOriginalRpSelectedPulse owner), proves a nonzero
finite-mu saddle error allowance, and retains the original convex tangent
and late-support tail bounds.

No selected amplitude, end coefficient, pressure moment, higher repair, or
global stress claim is made here.
"""
import copy
import gzip
import hashlib
import json
import time
from dataclasses import dataclass
from functools import wraps
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext


HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_original_Rp_pulse_saddle_rows.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_pulse_saddle_rows_check.json'
OUTER_SOURCE = PREFIX + 'outer_pulse_map.py'
OUTER_RAW = PREFIX + 'outer_pulse_map.json'
OUTER_CHECK_SOURCE = PREFIX + 'outer_pulse_map_check.py'
OUTER_CHECK = PREFIX + 'outer_pulse_map_check.json'
GATE = 'current_original_Rp_pulse_saddle_rows_refined_installed'
OPEN = (
    'current_original_Rp_selected_amplitude_installed',
    'current_original_Rp_end_coefficients_installed',
    'current_original_Rp_n1_pressure_moments_installed',
    'current_original_Rp_n1_full_tails_installed',
    'current_original_Rp_five_moment_repair_installed',
    'current_original_Rp_global_stress_installed',
    'current_original_Rp_global_NS_solution_installed',
)
SOURCE_DPS = 540
ROW_NAMES = ('F1', 'F2')


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _source_precision(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        with mp.workdps(SOURCE_DPS):
            return fn(*args, **kwargs)
    return call


def ep(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def record(value):
    if hasattr(value, '_mpi_'):
        lo, hi = ep(value)
        return dict(lower=str(lo), upper=str(hi), width=str(hi - lo),
                    lower_exact_mpf_tuple=list(value._mpi_[0]),
                    upper_exact_mpf_tuple=list(value._mpi_[1]))
    if isinstance(value, dict):
        return {key: record(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [record(item) for item in value]
    return value


def read_interval(c, value):
    if hasattr(value, '_mpi_'):
        return value
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                       tuple(value['upper_exact_mpf_tuple'])))


def _read_json(path, compressed=False):
    data = Path(path).read_bytes()
    if compressed:
        data = gzip.decompress(data)
    return json.loads(data)


def _family(raw):
    return dict(actual_five_defect_family_sha256=raw['actual_five_defect_family_sha256'],
                implicit_source_sha256=raw['implicit_source_sha256'],
                datum_enclosure_sha256=raw['datum_enclosure_sha256'])


def _verify_input_hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if sha(name) != digest:
            raise ValueError('Accepted pulse source changed: ' + name)


def _load_accepted_outer_source():
    raw = _read_json(HERE / OUTER_RAW)
    check = _read_json(HERE / OUTER_CHECK)
    _verify_input_hashes(raw)
    _verify_input_hashes(check)
    if raw.get('actual_ap_selected') or not raw.get('near_equal_rows_recovered_by_analytic_divided_difference'):
        raise ValueError('Original pulse source selection branch changed')
    if not check.get('actual_trial_source_and_interval_inverse_checked'):
        raise ValueError('Accepted original pulse source check missing')
    if not check.get('positive_unresolved_pulse_and_incoming_tails_preserved'):
        raise ValueError('Accepted positive pulse tail provenance missing')
    if check.get('actual_ap_selected') or check.get('complete_O4_corrected_field_built'):
        raise ValueError('Pulse source record has expanded beyond its accepted scope')
    rows = raw.get('pulse_rows', {})
    for key in ('common_logpref', 'scaled_full_rows', 'positive_scaled_tail_bounds',
                'saddle_u0', 'saddle_L', 'band', 'flat_denominator_positive_cap'):
        if key not in rows:
            raise ValueError('Accepted pulse row record missing ' + key)
    if len(rows['scaled_full_rows']) != 2 or len(rows['positive_scaled_tail_bounds']) != 2:
        raise ValueError('Exactly two accepted original pulse rows required')
    return raw, check, _family(raw)


def _same_family(typed_owner, family):
    actual = getattr(typed_owner, 'family_record', None)
    if actual != family:
        raise ValueError('Typed selected pulse family differs from outer pulse source')


def _typed_binding(typed_owner, family, raw):
    if typed_owner is None:
        return dict(kind='accepted_outer_pulse_source_record',
                    source_module=OUTER_SOURCE, raw_record=OUTER_RAW,
                    checker_record=OUTER_CHECK, typed_owner_required=False,
                    selected_branch_unchanged=True)
    if not getattr(typed_owner, 'acceptance_loaded', False):
        raise ValueError('Accepted typed CurrentOriginalRpSelectedPulse required')
    _same_family(typed_owner, family)
    graph = typed_owner.assert_graph()
    if isinstance(graph, dict) and not all(graph.values()):
        raise ValueError('Typed selected pulse graph is not live')
    pulse_facade = getattr(typed_owner, 'pulse', None)
    pulse_map = getattr(pulse_facade, 'pulse', None)
    if pulse_map is None or not hasattr(pulse_map, 'rows') or not hasattr(pulse_map, 'mu'):
        raise ValueError('Typed selected pulse does not expose the original fixed pulse rows')
    old_mu = read_interval(typed_owner.ctx, raw['mu'])
    typed_mu = pulse_map.mu
    alo, ahi = ep(old_mu)
    blo, bhi = ep(typed_mu)
    if bhi < alo or ahi < blo:
        raise ValueError('Typed pulse mu and accepted source mu are disjoint')
    return dict(kind='live_typed_CurrentOriginalRpSelectedPulse',
                selected_owner_type=type(typed_owner).__name__,
                selected_branch_unchanged=True,
                typed_original_pulse_owner=type(pulse_map).__name__,
                source_module=OUTER_SOURCE, raw_record=OUTER_RAW,
                checker_record=OUTER_CHECK, typed_owner_required=True,
                typed_mu_overlaps_accepted_record=True,
                original_rows_path='selected_owner.pulse.pulse.rows',
                original_mu_path='selected_owner.pulse.pulse.mu')


def _max_interval(c, left, right):
    l0, l1 = ep(left)
    r0, r1 = ep(right)
    return c.mpf([max(l0, r0), max(l1, r1)])


def _min_interval(c, left, right):
    l0, l1 = ep(left)
    r0, r1 = ep(right)
    return c.mpf([min(l0, r0), min(l1, r1)])


def _relative_width(c, value):
    lo_tuple = tuple(value._mpi_[0])
    lower = c.make_mpf((lo_tuple, lo_tuple))
    width = value - lower
    return width / lower


def _tangent_tail(c, mu, u0, w, row, band):
    gaussian = c.mpf(0)
    pieces = []
    for sign in (-1, 1):
        x = c.mpf(sign * band)
        v = 1 + x * w
        phi = x**2 * (2 * v + 1) / (6 * v**2)
        derivative = x * (v**2 + v + 1) / (3 * v**3)
        dlo, dhi = ep(derivative)
        if dlo == 0 or dhi == 0 or dlo * dhi > 0:
            absolute = c.mpf([min(abs(dlo), abs(dhi)), max(abs(dlo), abs(dhi))])
        else:
            raise ValueError('Convex tangent derivative crosses zero')
        piece = c.exp(-phi) / absolute
        gaussian += piece
        pieces.append(dict(sign=sign, phi=record(phi), derivative=record(derivative),
                           tangent_integral=record(piece)))
    tail = 11 * c.exp(4 + c.mpf('2.5') * row) * gaussian
    return tail, pieces


def _refine_row(c, mu, u0, w, common, old_row, old_tail, row, band, allowance, cap):
    band_c = c.mpf(band)
    u_hi = ep(u0 * (1 + band_c * w))[1]
    w_hi = ep(w)[1]
    guard_u = u_hi <= c.mpf('.25')
    guard_w = abs(band_c * w_hi) <= c.mpf('.25')
    if not guard_u or not guard_w:
        raise ArithmeticError('Saddle band guards u<=1/4 and |48w|<=1/4 failed')
    band_u = u0 * (1 + band_c * w)
    flat_denominator_log = -1 / band_u**2 + 1 / (1 - band_u)**2
    if ep(flat_denominator_log)[1] >= -1000:
        raise ArithmeticError('Flat denominator cap is not below -1000')
    error = ((c.mpf(2) / c.mpf('10.99') + 4 + row) * u_hi
             + 2 * band_c**3 * w_hi + cap)
    if ep(error)[1] <= 0 or not ep(error)[1] < ep(allowance)[0]:
        raise ArithmeticError('Nonzero finite-mu saddle error is not below allowance')
    amplitude = c.mpf('10.99') * c.exp(1 + 2 * row)
    reference = amplitude * c.sqrt(2 * c.pi)
    gaussian_tail = c.mpf(2) * c.exp(-band_c**2 / 2) / band_c
    core_lower = reference * c.exp(-allowance) - amplitude * gaussian_tail
    core_upper = reference * c.exp(allowance)
    tangent, tangent_pieces = _tangent_tail(c, mu, u0, w, row, band)
    support_upper = tangent + cap
    # The archived tail is a genuine old enclosure.  Retain it if its
    # directed rounding is wider than the reconstructed tangent formula.
    tail_bound = _max_interval(c, old_tail, support_upper)
    full_lower = core_lower
    full_upper = core_upper + tail_bound
    lower_lo, lower_hi = ep(full_lower)
    upper_lo, upper_hi = ep(full_upper)
    old_lo, old_hi = ep(old_row)
    refined_lo = max(old_lo, lower_lo)
    refined_hi = min(old_hi, upper_hi)
    if refined_lo > refined_hi:
        raise ArithmeticError('Refined pulse row does not intersect accepted full row')
    refined = c.mpf([refined_lo, refined_hi])
    late_log = c.ln(121) - c.mpf(1) / (4 * mu) + c.mpf(3) / (u0**2) + 13 * row \
        - 2 * c.ln(u0) + c.ln(6) / 2
    if ep(late_log)[1] >= -1000:
        raise ArithmeticError('Original late support log cap is not below -1000')
    return dict(row=row, name=ROW_NAMES[row - 1], refined=refined,
                old=old_row, core=c.mpf([ep(core_lower)[0], ep(core_upper)[1]]),
                old_tail=old_tail, tangent_tail=tangent, support_upper=support_upper,
                tail_bound=tail_bound, late_log=late_log,
                log_ratio_allowance=allowance, log_ratio_error=error,
                relative_width_old=_relative_width(c, old_row),
                relative_width_new=_relative_width(c, refined),
                flat_denominator_log=flat_denominator_log,
                tangent_pieces=tangent_pieces,
                support_tail_formula='11*exp(4+2.5*i)*sum(exp(-phi_endpoint)/|phi_prime_endpoint|)+exp(-1000)',
                late_support_cap='log(121)-1/(4*mu)+3/u0^2+13*i-2*log(u0)+log(6)/2 < -1000')


@dataclass(frozen=True)
class OriginalRpPulseSaddleRowsPacket:
    marker: str = 'original_Rp_pulse_saddle_rows'


class CurrentOriginalRpPulseSaddleRows:
    @_source_precision
    def __init__(self, selected_owner=None, require_checked=True):
        raw, check, family = _load_accepted_outer_source()
        self.family = family
        self.selected_owner = selected_owner
        self.source_binding = _typed_binding(selected_owner, family, raw)
        self.ctx = MPIntervalContext()
        self.ctx.dps = SOURCE_DPS
        c = self.ctx
        self.raw_outer = raw
        self.outer_check = check
        self.mu = read_interval(c, raw['mu'])
        rows = raw['pulse_rows']
        self.common_logpref = read_interval(c, rows['common_logpref'])
        self.old_rows = [read_interval(c, value) for value in rows['scaled_full_rows']]
        self.old_tails = [read_interval(c, value) for value in rows['positive_scaled_tail_bounds']]
        self.u0 = c.exp(c.ln(4 * self.mu) / 3)
        self.w = self.u0 / c.sqrt(6)
        self.allowance = c.mpf('1e-60')
        self.cap = c.exp(-1000)
        self.band = int(rows['band'])
        self.refined_rows = []
        for row in (1, 2):
            self.refined_rows.append(_refine_row(
                c, self.mu, self.u0, self.w, self.common_logpref,
                self.old_rows[row - 1], self.old_tails[row - 1], row,
                self.band, self.allowance, self.cap))
        self.hashes = dict(raw.get('input_hashes', {}))
        self.hashes.update({name: sha(name) for name in
                            (OUTER_CHECK_SOURCE, OUTER_CHECK, Path(__file__).name)})
        self.acceptance_loaded = False
        self._view = None
        if require_checked:
            receipt = _read_json(HERE / RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or
                    receipt.get('source_family') != self.family or
                    any(receipt.get(key, False) for key in OPEN)):
                raise ValueError('Checked pulse saddle rows required')
            for name, digest in receipt.get('input_hashes', {}).items():
                if sha(name) != digest:
                    raise ValueError('Changed pulse saddle dependency: ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hash_snapshot = copy.deepcopy(self.hashes)

    @_source_precision
    def assert_graph(self):
        if self.hashes != self._hash_snapshot:
            raise ValueError('Pulse saddle source state changed')
        raw, check, family = _load_accepted_outer_source()
        if family != self.family:
            raise ValueError('Pulse saddle source family changed')
        if self.selected_owner is not None:
            _typed_binding(self.selected_owner, self.family, raw)
        return dict(accepted_outer_source=True,
                    selected_branch_unchanged=self.source_binding['selected_branch_unchanged'],
                    common_logpref_source_preserved=True,
                    old_scaled_rows_preserved=True)

    def _view_data(self):
        c = self.ctx
        guard_u = ep(self.u0 * (1 + self.band * self.w))[1] <= c.mpf('.25')
        guard_w = abs(self.band * ep(self.w)[1]) <= c.mpf('.25')
        return dict(
            source_family=self.family,
            source_binding=copy.deepcopy(self.source_binding),
            original_mu=record(self.mu), original_saddle_u0=record(self.u0),
            original_saddle_w=record(self.w), original_common_logpref=record(self.common_logpref),
            old_scaled_full_rows=[record(value) for value in self.old_rows],
            refined_scaled_full_rows=[record(value['refined']) for value in self.refined_rows],
            refined_rows_by_name={value['name']: record(value['refined']) for value in self.refined_rows},
            core_rows={value['name']: record(value['core']) for value in self.refined_rows},
            old_positive_tail_bounds=[record(value) for value in self.old_tails],
            row_diagnostics=[dict(
                row=value['row'], name=value['name'],
                refined=record(value['refined']), old=record(value['old']),
                core=record(value['core']), old_tail=record(value['old_tail']),
                tangent_tail=record(value['tangent_tail']),
                support_upper=record(value['support_upper']),
                tail_bound=record(value['tail_bound']),
                late_log=record(value['late_log']),
                log_ratio_allowance=record(value['log_ratio_allowance']),
                log_ratio_error=record(value['log_ratio_error']),
                flat_denominator_log=record(value['flat_denominator_log']),
                relative_width_old=record(value['relative_width_old']),
                relative_width_new=record(value['relative_width_new']),
                tangent_pieces=value['tangent_pieces'],
                late_support_cap=value['late_support_cap'],
                support_tail_formula=value['support_tail_formula'])
                for value in self.refined_rows],
            band=self.band,
            guards=dict(u_le_one_quarter=guard_u, abs_48w_le_one_quarter=guard_w,
                        all_rows_error_below_allowance=all(
                            ep(value['log_ratio_error'])[1] < ep(self.allowance)[0]
                            for value in self.refined_rows),
                        allowance=record(self.allowance),
                        nonzero_mu=ep(self.mu)[0] > 0,
                        nonzero_u0=ep(self.u0)[0] > 0,
                        nonzero_w=ep(self.w)[0] > 0),
            exact_kernel_ratio_identity=(
                'log((10.99-u)/10.99)+(1-u)^(-2)-1+i*u'
                '-(phi-x^2/2)-log(1+exp(-1/u^2+(1-u)^(-2)))'),
            phase_identity='phi-x^2/2=-x^3*w*(3*v+1)/(6*v^2)',
            natural_reference='10.99*exp(1+2*i)*exp(-x^2/2)',
            original_convex_tangent_tail_reused=True,
            original_late_support_cap_reused=True,
            original_common_logpref_unchanged=True,
            original_selection_branch_unchanged=True,
            finite_mu_error_proved_nonzero=True,
            scaled_rows_refined_only=True,
            no_radius_or_N_materialization=True,
            numerical_evidence_row_enclosures_only=True,
            **{GATE: self.acceptance_loaded},
            **dict.fromkeys(OPEN, False))

    @_source_precision
    def evaluate(self):
        self.assert_graph()
        self._view = self._view_data()
        return self._view

    @_source_precision
    def report(self, view):
        self.assert_graph()
        if self._view is not view:
            raise ValueError('Live pulse saddle view required')
        return copy.deepcopy(view)


@_source_precision
def run(selected_owner=None):
    began = time.monotonic()
    owner = CurrentOriginalRpPulseSaddleRows(selected_owner, require_checked=False)
    view = owner.evaluate()
    raw = dict(candidate_original_Rp_pulse_saddle_rows_constructed=True,
               source_family=owner.family,
               source_binding=owner.source_binding,
               source_records=dict(outer_raw=OUTER_RAW, outer_check=OUTER_CHECK,
                                   outer_source=OUTER_SOURCE),
               view=owner.report(view), input_hashes=owner.hashes,
               **{GATE: False}, **dict.fromkeys(OPEN, False),
               execution_seconds=time.monotonic() - began)
    data = json.dumps(raw, separators=(',', ':')) + '\n'
    (HERE / NAME).write_bytes(gzip.compress(data.encode('utf8'), mtime=0))
    print('CURRENT_RP_PULSE_SADDLE_ROWS_REFINED', 2, flush=True)
    return owner, view


if __name__ == '__main__':
    run()
