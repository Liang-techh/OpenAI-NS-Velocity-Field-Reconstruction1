"""Original scaled F1 index derivative for the accepted Rp saddle source.

The accepted saddle row is an integral of the positive original kernel with
the finite index factor ``exp(i * (2 + u))``.  The common logarithmic
prefactor is independent of the index, so differentiation gives

    F1prime = integral (2 + u) * kernel_1

The core uses the genuine upper bound ``u <= u_hi`` and the accepted positive
tail uses ``u <= 11``.  No radius, N, amplitude, pressure, or full-tail value
is materialized here.
"""

import copy
import gzip
import hashlib
import json
import time
from functools import wraps
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_rows as saddle


HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_original_Rp_pulse_saddle_index_derivative.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_pulse_saddle_index_derivative_check.json'
GATE = 'current_original_Rp_pulse_saddle_index_derivative_installed'
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


def _weighted_bounds(c, saddle_owner):
    """Return directed core/tail bounds for integral ``u * kernel_1``."""
    row = saddle_owner.refined_rows[0]
    core = row['core']
    tail = row['tail_bound']
    band_u = saddle_owner.u0 * (1 + saddle_owner.band * saddle_owner.w)
    u_hi = ep(band_u)[1]
    u_lo = ep(saddle_owner.u0 *
              (1 - saddle_owner.band * saddle_owner.w))[0]
    core_lo, core_hi = ep(core)
    tail_hi = ep(tail)[1]
    weighted_core = c.mpf([u_lo * core_lo, u_hi * core_hi])
    weighted_tail = c.mpf([0, c.mpf(11) * tail_hi])
    weighted_total = c.mpf([
        ep(weighted_core)[0] + ep(weighted_tail)[0],
        ep(weighted_core)[1] + ep(weighted_tail)[1]])
    return dict(
        u_core_lower=c.mpf(u_lo),
        u_core_upper=c.mpf(u_hi),
        core_integral=core,
        accepted_positive_tail=tail,
        weighted_core=weighted_core,
        weighted_tail=weighted_tail,
        weighted_total=weighted_total,
        formulas=dict(
            core='0 <= u*K1 <= u_hi*K1 on |x|<=48',
            tail='0 <= u*K1 <= 11*K1 on the accepted positive support tail',
            total='weighted_core + weighted_tail'))


class CurrentOriginalRpPulseSaddleIndexDerivative:
    @_source_precision
    def __init__(self, saddle_owner=None, require_checked=True):
        if saddle_owner is None:
            saddle_owner = saddle.CurrentOriginalRpPulseSaddleRows(
                require_checked=True)
        if type(saddle_owner) is not saddle.CurrentOriginalRpPulseSaddleRows:
            raise ValueError('Checked original saddle-row owner required')
        if require_checked and not saddle_owner.acceptance_loaded:
            raise ValueError('Checked original saddle-row owner required')
        saddle_owner.assert_graph()
        self.saddle = saddle_owner
        self.ctx = c = saddle_owner.ctx
        self.family = copy.deepcopy(saddle_owner.family)
        self.mu = saddle_owner.mu
        self.u0 = saddle_owner.u0
        self.w = saddle_owner.w
        self.band = saddle_owner.band
        self.common_logpref = saddle_owner.common_logpref
        self.f1 = saddle_owner.refined_rows[0]['refined']
        self.bounds = _weighted_bounds(c, saddle_owner)
        self.weighted_kernel = self.bounds['weighted_total']
        self.index_derivative = self.f1 * 2 + self.weighted_kernel
        self.source_binding = dict(
            kind='accepted_original_Rp_saddle_rows',
            source_module=Path(saddle.__file__).name,
            source_owner_type=type(saddle_owner).__name__,
            source_row_path='saddle_owner.refined_rows[0].refined',
            source_tail_path='saddle_owner.refined_rows[0].tail_bound',
            source_family=self.family,
            common_logpref_preserved=True,
            original_index_derivative='integral (2+u)*kernel_i1 over u in [0,11]',
            selected_branch_unchanged=True)
        self.hashes = dict(saddle_owner.hashes)
        for name in (saddle.NAME, saddle.RECEIPT, Path(saddle.__file__).name,
                     Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded = False
        if require_checked:
            receipt = _read_json(HERE / RECEIPT)
            if (not receipt.get('all_passed') or
                    not receipt.get(GATE) or
                    receipt.get('source_family') != self.family or
                    any(receipt.get(key, False) for key in OPEN)):
                raise ValueError('Checked index-derivative receipt required')
            for name, digest in receipt.get('input_hashes', {}).items():
                if sha(name) != digest:
                    raise ValueError('Changed index-derivative dependency: ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hash_snapshot = copy.deepcopy(self.hashes)
        self._state = record((self.family, self.source_binding, self.bounds,
                              self.f1, self.index_derivative))

    @_source_precision
    def assert_graph(self):
        self.saddle.assert_graph()
        if self.hashes != self._hash_snapshot:
            raise ValueError('Index-derivative source state changed')
        if self._state != record((self.family, self.source_binding, self.bounds,
                                  self.f1, self.index_derivative)):
            raise ValueError('Index-derivative values changed')
        return dict(saddle_source_live=True, source_family=self.family,
                    common_logpref_preserved=True)

    @_source_precision
    def evaluate(self):
        self.assert_graph()
        c = self.ctx
        band_u = self.u0 * (1 + self.band * self.w)
        view = dict(
            source_family=self.family,
            source_binding=copy.deepcopy(self.source_binding),
            saddle_source_hashes=copy.deepcopy(self.saddle.hashes),
            original_mu=record(self.mu),
            original_saddle_u0=record(self.u0),
            original_saddle_w=record(self.w),
            original_common_logpref=record(self.common_logpref),
            original_F1=record(self.f1),
            weighted_core=record(self.bounds['weighted_core']),
            weighted_tail=record(self.bounds['weighted_tail']),
            weighted_u_kernel=record(self.weighted_kernel),
            F1prime=record(self.index_derivative),
            core_integral=record(self.bounds['core_integral']),
            accepted_positive_tail=record(self.bounds['accepted_positive_tail']),
            formulas=copy.deepcopy(self.bounds['formulas']),
            exact_derivative_identity=(
                'd/di exp(i*(2+u))=(2+u)*exp(i*(2+u)); '
                'F1prime=2*F1+integral(u*kernel_1)'),
            support='original u in [0,11], saddle core |x|<=48 plus accepted positive tail',
            guards=dict(
                u_core_lower_positive=ep(self.bounds['u_core_lower'])[0] > 0,
                u_core_upper_le_quarter=ep(band_u)[1] <= c.mpf('.25'),
                abs_48w_le_quarter=abs(self.band * ep(self.w)[1]) <= c.mpf('.25'),
                weighted_core_positive=ep(self.bounds['weighted_core'])[0] > 0,
                weighted_tail_positive=ep(self.bounds['weighted_tail'])[1] > 0,
                weighted_u_remainder_nonzero=ep(self.weighted_kernel)[1] > 0),
            original_source_rows_unchanged=True,
            original_positive_tail_reused=True,
            common_logpref_index_independent=True,
            finite_u_remainder_retained=True,
            no_radius_or_N_materialization=True,
            numerical_evidence_directed_source_enclosure=True,
            **{GATE: self.acceptance_loaded},
            **dict.fromkeys(OPEN, False))
        return view

    @_source_precision
    def report(self, view):
        self.assert_graph()
        if view is not getattr(self, '_view', view):
            # Keep compatibility with callers that pass a freshly evaluated
            # live object while rejecting unrelated dictionaries.
            if view.get('source_family') != self.family:
                raise ValueError('Live index-derivative view required')
        return copy.deepcopy(view)


@_source_precision
def run(saddle_owner=None):
    began = time.monotonic()
    owner = CurrentOriginalRpPulseSaddleIndexDerivative(
        saddle_owner=saddle_owner, require_checked=False)
    view = owner.evaluate()
    owner._view = view
    raw = dict(
        candidate_original_Rp_pulse_saddle_index_derivative_constructed=True,
        source_family=owner.family,
        source_binding=owner.source_binding,
        source_records=dict(
            saddle_rows=saddle.NAME,
            saddle_receipt=saddle.RECEIPT,
            saddle_source=Path(saddle.__file__).name),
        view=owner.report(view),
        input_hashes=owner.hashes,
        **{GATE: False},
        **dict.fromkeys(OPEN, False),
        execution_seconds=time.monotonic() - began)
    data = json.dumps(raw, separators=(',', ':'), ensure_ascii=False) + '\n'
    (HERE / NAME).write_bytes(gzip.compress(data.encode('utf8'), mtime=0))
    print('CURRENT_RP_PULSE_SADDLE_INDEX_DERIVATIVE', flush=True)
    return owner, view


if __name__ == '__main__':
    run()
