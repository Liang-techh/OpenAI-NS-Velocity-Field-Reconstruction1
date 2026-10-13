"""Refined original correction-basis and fixed-pulse energy enclosures.

The accepted pulse source fixes the shape, normalization, finite ``mu`` and
source family.  This adapter narrows only the old bump basis and ``Kpulse``
inputs.  It does not select an amplitude or claim any later pressure, repair,
stress, or Navier--Stokes gate.
"""
import copy
import gzip
import hashlib
import json
import time
from functools import wraps
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_rows as saddle


HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_original_Rp_pulse_basis_energy.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_pulse_basis_energy_check.json'
OUTER_SOURCE = PREFIX + 'outer_pulse_map.py'
OUTER_RAW = PREFIX + 'outer_pulse_map.json'
OUTER_CHECK_SOURCE = PREFIX + 'outer_pulse_map_check.py'
OUTER_CHECK = PREFIX + 'outer_pulse_map_check.json'
SADDLE_SOURCE = saddle.PREFIX + 'current_original_Rp_pulse_saddle_rows.py'
SADDLE_RAW = saddle.NAME
SADDLE_CHECK_SOURCE = saddle.PREFIX + 'current_original_Rp_pulse_saddle_rows_check.py'
SADDLE_RECEIPT = saddle.RECEIPT
NORM_SOURCE = 'lei_ren_part1_paper_bump_integral_enclosures.py'
NORM_CHECK_SOURCE = 'lei_ren_part1_paper_bump_integral_enclosures_check.py'
NORM_RAW = 'lei_ren_part1_paper_bump_integral_enclosures_check.json'
GATE = 'current_original_Rp_pulse_basis_energy_refined_installed'
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
BASIS_CELLS = 2048
ENERGY_CELLS = 4096
ELL_TEXT = '0.15'
ROWS = ('first_row', 'exact_divided_difference_row')


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
    if isinstance(value, (list, tuple)):
        return [record(item) for item in value]
    return value


def read_interval(c, value):
    if hasattr(value, '_mpi_'):
        return value
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                       tuple(value['upper_exact_mpf_tuple'])))


def _read_json(name, compressed=False):
    data = (HERE / name).read_bytes()
    if compressed:
        data = gzip.decompress(data)
    return json.loads(data)


def _family(raw):
    return dict(actual_five_defect_family_sha256=raw['actual_five_defect_family_sha256'],
                implicit_source_sha256=raw['implicit_source_sha256'],
                datum_enclosure_sha256=raw['datum_enclosure_sha256'])


def _verify_hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if sha(name) != digest:
            raise ValueError('Changed accepted source: ' + name)


def _intersect(c, left, right):
    l0, l1 = ep(left)
    r0, r1 = ep(right)
    result = c.mpf([max(l0, r0), min(l1, r1)])
    if ep(result)[0] > ep(result)[1]:
        raise ArithmeticError('Refined interval misses accepted source interval')
    return result


def _abs_upper(value):
    lo, hi = ep(value)
    return max(abs(lo), abs(hi))


def _load_inputs(saddle_owner=None):
    outer = _read_json(OUTER_RAW)
    outer_check = _read_json(OUTER_CHECK)
    _verify_hashes(outer)
    _verify_hashes(outer_check)
    if outer.get('actual_ap_selected') or not outer.get('near_equal_rows_recovered_by_analytic_divided_difference'):
        raise ValueError('Original pulse source selection branch changed')
    if outer_check.get('actual_ap_selected') or outer_check.get('complete_O4_corrected_field_built'):
        raise ValueError('Original pulse source expanded beyond accepted scope')

    if saddle_owner is None:
        saddle_owner = saddle.CurrentOriginalRpPulseSaddleRows(require_checked=True)
    else:
        if not getattr(saddle_owner, 'acceptance_loaded', False):
            raise ValueError('Checked saddle-row owner required')
        saddle_owner.assert_graph()
    family = saddle_owner.family
    if family != _family(outer):
        raise ValueError('Saddle and original pulse family differ')

    norm = _read_json(NORM_RAW)
    if not norm.get('directed_integrals_certified'):
        raise ValueError('Accepted bump normalization is not certified')
    if not norm.get('normalization_contains_one_for_all_centers'):
        raise ValueError('Accepted bump normalization support check missing')
    basis = outer['bump_basis']
    rows = outer['pulse_rows']
    if basis.get('beta_radius', {}).get('lower') != ELL_TEXT or basis.get('beta_radius', {}).get('upper') != ELL_TEXT:
        raise ValueError('Original beta radius is not ell=.15')
    if len(basis.get('first_row', [])) != 2 or len(basis.get('exact_divided_difference_row', [])) != 2:
        raise ValueError('Original 2 by 2 bump basis missing')
    return dict(outer=outer, outer_check=outer_check, norm=norm,
                basis=basis, rows=rows, family=family,
                saddle_owner=saddle_owner,
                mu=saddle_owner.mu)


def _bump_derivative_bounds(c, n):
    """Global analytic bounds from exp(-n/D) D^-k, D=1-r^2."""
    def q(k):
        if k == 0:
            return c.exp(-n)
        return (c.mpf(k) / (c.mpf(n) * c.exp(1)))**k
    b0 = q(0)
    b1 = 2 * n * q(2)
    b2 = 4 * n * n * q(4) + 2 * n * q(2) + 8 * n * q(3)
    return b0, b1, b2


def _midpoint_integral(c, t, ell, n, kind, cells):
    """Composite midpoint enclosure with a global analytic f'' remainder."""
    cells = int(cells)
    h = c.mpf(2) / cells
    total = c.mpf(0)
    b0, b1, b2 = _bump_derivative_bounds(c, n)
    t_lo, t_hi = ep(t)
    t_abs = c.mpf([min(abs(t_lo), abs(t_hi)), max(abs(t_lo), abs(t_hi))])
    exp_t = c.exp(t_abs)
    if kind == 'j':
        sup = exp_t * (t_abs**2 * b0 + 2 * t_abs * b1 + b2)
    elif kind == 'jr':
        g1 = exp_t * (t_abs * b0 + b1)
        g2 = exp_t * (t_abs**2 * b0 + 2 * t_abs * b1 + b2)
        sup = ell * (g2 + 2 * g1)
    elif kind == 'gram':
        sup = c.mpf(b2)
    else:
        raise ValueError('Unknown midpoint integrand')
    for i in range(cells):
        left = c.mpf(-1) + c.mpf(2 * i) / cells
        right = c.mpf(-1) + c.mpf(2 * (i + 1)) / cells
        midpoint = (left + right) / 2
        denominator = 1 - midpoint**2
        b = c.exp(-n / denominator)
        if kind == 'j':
            value = c.exp(t * ell * midpoint) * b
        elif kind == 'jr':
            value = ell * midpoint * c.exp(t * ell * midpoint) * b
        else:
            value = b
        total += h * value
    remainder = c.mpf(2) * h**2 * sup / 24
    lo, hi = ep(total)
    return c.mpf([lo - ep(remainder)[1], hi + ep(remainder)[1]]), dict(
        kind=kind, cells=cells, midpoint_step=record(h),
        global_second_derivative_bound=record(sup),
        bump_derivative_components=dict(b0=record(b0), b1=record(b1), b2=record(b2)),
        midpoint_sum=record(total), remainder=record(remainder),
        remainder_formula='length*h^2*sup_abs_f_second/24',
        derivative_bound_formula='exp(-n/D)*D^(-k)<= (k/(n*e))^k for k>0',
    )


def _finite_factor(c, mu, bound):
    return c.mpf([c.exp(-c.mpf(bound) * mu), c.exp(c.mpf(bound) * mu)])


def _refined_basis(c, inputs):
    mu = inputs['mu']
    ell = c.mpf(ELL_TEXT)
    norm = read_interval(c, inputs['norm']['normalization'])
    # Ordinary mu=0 moments for both later end-pressure exponents.
    moment_data = {}
    for a_text in ('0.5', '1.5'):
        a = c.mpf(a_text)
        j_num, j_diag = _midpoint_integral(c, a, ell, 1, 'j', BASIS_CELLS)
        jr_num, jr_diag = _midpoint_integral(c, a, ell, 1, 'jr', BASIS_CELLS)
        moment_data[a_text] = dict(
            J=j_num / norm, J_r=jr_num / norm,
            J_numerator=j_num, J_r_numerator=jr_num,
            J_diagnostics=j_diag, J_r_diagnostics=jr_diag,
        )

    # For the finite source, A_d uses J(1/2-mu).  The D row is enclosed
    # without subtracting nearly equal boxes, using its positive mu=0 form
    # and the true exp(+-8 mu) factor.
    half_mu = c.mpf('0.5') - mu
    half_2mu = c.mpf('0.5') - 2 * mu
    j1_num, j1_diag = _midpoint_integral(c, half_mu, ell, 1, 'j', BASIS_CELLS)
    j2_num, j2_diag = _midpoint_integral(c, half_2mu, ell, 1, 'j', BASIS_CELLS)
    j1 = j1_num / norm
    j2 = j2_num / norm
    j0 = moment_data['0.5']['J']
    jr0 = moment_data['0.5']['J_r']
    A0 = []
    D0 = []
    A = []
    D = []
    for d in (3, 1):
        a0 = c.exp(-c.mpf(d) / 2) * j0
        d0 = c.exp(-c.mpf(d) / 2) * (d * j0 - jr0)
        if ep(d0)[0] <= 0:
            raise ArithmeticError('Positive divided-difference reference lost')
        A0.append(a0)
        D0.append(d0)
        # Exact A source factor is exp(-d/2 + mu*d) J(1/2-mu); the
        # symmetric factor is retained as a conservative original-mu guard.
        a_actual = c.exp(-c.mpf(d) / 2 + mu * d) * j1
        d_actual = d0 * _finite_factor(c, mu, 8)
        A.append(a_actual)
        D.append(d_actual)
    q_mu = 2 * c.mpf([1, ep(c.exp(2 * mu))[1]])
    det_abs = c.exp(-2 + 5 * mu) * j1 * j2 * q_mu
    determinant = -det_abs
    if ep(determinant)[1] >= 0:
        raise ArithmeticError('Refined divided-difference determinant is not negative')

    gram_num, gram_diag = _midpoint_integral(c, c.mpf(0), ell, 2, 'gram', BASIS_CELLS)
    gram0 = gram_num / (ell * norm**2)
    gram = gram0 * _finite_factor(c, mu, c.mpf('.3'))

    old = inputs['basis']
    old_first = [read_interval(c, value) for value in old['first_row']]
    old_divided = [read_interval(c, value) for value in old['exact_divided_difference_row']]
    first = [_intersect(c, old_first[i], A[i]) for i in range(2)]
    divided = [_intersect(c, old_divided[i], D[i]) for i in range(2)]
    old_det = read_interval(c, old['divided_difference_determinant'])
    determinant = _intersect(c, old_det, determinant)
    old_gram = read_interval(c, old['energy_gram'])
    gram = _intersect(c, old_gram, gram)
    return dict(
        first_row=first, exact_divided_difference_row=divided,
        divided_difference_determinant=determinant, energy_gram=gram,
        end_energy_weights=[c.exp(-26 + 6 * mu) * gram,
                            c.exp(-26 + 2 * mu) * gram],
        beta_radius=ell, cells=BASIS_CELLS,
        original_mu_correction=dict(A_factor='exp(-4*mu)..exp(4*mu)',
                                    D_factor='exp(-8*mu)..exp(8*mu)',
                                    gram_factor='exp(-.3*mu)..exp(.3*mu)',
                                    mu_nonzero=True),
        moments=moment_data,
        finite_moments=dict(J_half_minus_mu=j1, J_half_minus_2mu=j2,
                            J_half_minus_mu_diagnostics=j1_diag,
                            J_half_minus_2mu_diagnostics=j2_diag),
        determinant_identity='det=-exp(-2+5*mu)*J(.5-mu)*J(.5-2*mu)*(exp(2*mu)-1)/mu',
        quotient_identity='(exp(2*mu)-1)/mu=2*integral_0^1 exp(2*mu*v)dv',
        gram_diagnostics=gram_diag,
        identity='row1=A_d; row2-row1=mu*D_d; x=d-s; beta=b(r)/(ell*N)',
    )


def _bulk_pulse_energy(c):
    left = c.mpf('.02')
    right = c.mpf(10)
    z = c.mpf('.01')
    def primitive(x):
        y = x - z
        return -c.exp(-2 * x) * (y**2 / 2 + y / 2 + c.mpf(1) / 4)
    bulk = primitive(right) - primitive(left)
    startup = c.mpf([0, c.mpf('.02')**3 / 3])
    # gp(xi)= (xi-.01)*sigma(11-xi) on [10,11], and 0<=sigma<=1.
    tail_upper = -primitive(10)
    tail = c.mpf([0, ep(tail_upper)[1]])
    return bulk, startup, tail, dict(
        bulk=record(bulk), startup=record(startup), exit=record(tail),
        bulk_primitive='-exp(-2*xi)*((xi-.01)^2/2+(xi-.01)/2+1/4)',
        startup_bound='0<=integral_0^.02 exp(-2*xi)*gp(xi)^2 dxi<=.02^3/3',
        exit_bound='0<=integral_10^11 ... <= integral_10^infinity exp(-2*xi)*(xi-.01)^2 dxi',
        sigma_symmetry_bulk=True, gp_exact_on_bulk=True,
    )


def _refined_energy(c, inputs):
    bulk, startup, tail, diag = _bulk_pulse_energy(c)
    candidate = bulk + startup + tail
    old = read_interval(c, inputs['outer']['fixed_pulse_energy_Kpulse'])
    refined = _intersect(c, old, candidate)
    if ep(refined)[0] <= 0:
        raise ArithmeticError('Refined fixed pulse energy lost positivity')
    diag['old'] = record(old)
    diag['refined'] = record(refined)
    diag['relative_width_old'] = record((old - ep(old)[0]) / ep(old)[0])
    diag['relative_width_refined'] = record((refined - ep(refined)[0]) / ep(refined)[0])
    diag['energy_cells_unused_for_analytic_route'] = ENERGY_CELLS
    return refined, diag


class CurrentOriginalRpPulseBasisEnergy:
    @_source_precision
    def __init__(self, saddle_owner=None, require_checked=True):
        inputs = _load_inputs(saddle_owner)
        self.ctx = MPIntervalContext()
        self.ctx.dps = SOURCE_DPS
        c = self.ctx
        self.inputs = inputs
        self.family = inputs['family']
        self.saddle_owner = inputs['saddle_owner']
        self.mu = read_interval(c, inputs['mu'])
        self.normalization = read_interval(c, inputs['norm']['normalization'])
        self.basis = _refined_basis(c, inputs)
        self.Kpulse, self.energy_diagnostics = _refined_energy(c, inputs)
        self.acceptance_loaded = False
        self._view = None
        self.hashes = dict(inputs['outer'].get('input_hashes', {}))
        self.hashes.update({name: sha(name) for name in (
            OUTER_RAW, OUTER_CHECK_SOURCE, OUTER_CHECK,
            SADDLE_SOURCE, SADDLE_RAW, SADDLE_CHECK_SOURCE, SADDLE_RECEIPT,
            NORM_SOURCE, NORM_CHECK_SOURCE, NORM_RAW, Path(__file__).name)})
        if require_checked:
            receipt = _read_json(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or
                    receipt.get('source_family') != self.family or
                    any(receipt.get(key, False) for key in OPEN)):
                raise ValueError('Checked basis-energy owner required')
            for name, digest in receipt.get('input_hashes', {}).items():
                if sha(name) != digest:
                    raise ValueError('Changed basis-energy dependency: ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hash_snapshot = copy.deepcopy(self.hashes)

    @_source_precision
    def assert_graph(self):
        if self.hashes != self._hash_snapshot:
            raise ValueError('Basis-energy source state changed')
        inputs = _load_inputs(self.saddle_owner)
        if inputs['family'] != self.family:
            raise ValueError('Basis-energy source family changed')
        return dict(accepted_outer_pulse_source=True,
                    accepted_saddle_rows=True,
                    normalization_source_preserved=True,
                    source_shapes_unchanged=True,
                    finite_mu_retained=True)

    def _view_data(self):
        return dict(
            source_family=self.family,
            source_binding=dict(outer_source=OUTER_SOURCE, outer_raw=OUTER_RAW,
                                saddle_source=SADDLE_SOURCE, saddle_raw=SADDLE_RAW,
                                normalization_raw=NORM_RAW,
                                original_shape_and_normalization_preserved=True),
            original_mu=record(self.mu),
            original_normalization=record(self.normalization),
            bump_basis={key: record(value) for key, value in self.basis.items()},
            fixed_pulse_energy_Kpulse=record(self.Kpulse),
            energy_diagnostics=copy.deepcopy(self.energy_diagnostics),
            basis_cells=BASIS_CELLS,
            finite_mu_nonzero=ep(self.mu)[0] > 0,
            normalization_source_retained=True,
            original_sigma_shape_retained=True,
            original_common_logpref_retained=True,
            no_selected_amplitude=True,
            no_radius_or_N_materialization=True,
            no_pressure_or_tail_closure=True,
            exact_bulk_energy_used=True,
            startup_unknown_positive_tail_retained=True,
            exit_unknown_positive_tail_retained=True,
            **{GATE: self.acceptance_loaded},
            **dict.fromkeys(OPEN, False),
        )

    @_source_precision
    def evaluate(self):
        self.assert_graph()
        self._view = self._view_data()
        return self._view

    @_source_precision
    def report(self, view):
        self.assert_graph()
        if self._view is not view:
            raise ValueError('Live basis-energy view required')
        return copy.deepcopy(view)


@_source_precision
def run(saddle_owner=None):
    began = time.monotonic()
    owner = CurrentOriginalRpPulseBasisEnergy(saddle_owner, require_checked=False)
    view = owner.evaluate()
    raw = dict(candidate_original_Rp_pulse_basis_energy_constructed=True,
               source_family=owner.family, source_binding=view['source_binding'],
               view=owner.report(view), input_hashes=owner.hashes,
               **{GATE: False}, **dict.fromkeys(OPEN, False),
               execution_seconds=time.monotonic() - began)
    data = json.dumps(raw, separators=(',', ':'), ensure_ascii=False) + '\n'
    (HERE / NAME).write_bytes(gzip.compress(data.encode('utf8'), mtime=0))
    print('CURRENT_RP_PULSE_BASIS_ENERGY_REFINED', flush=True)
    return owner, view


if __name__ == '__main__':
    run()
