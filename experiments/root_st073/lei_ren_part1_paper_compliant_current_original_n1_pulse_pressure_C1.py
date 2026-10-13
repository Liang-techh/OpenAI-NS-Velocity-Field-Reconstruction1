"""Archive-bound first-Z derivative bounds for the original Rp--Rv pulse.

This is deliberately a small adapter around the accepted pulse owner.  The
linear pressure jets are reused exactly.  The only new numerical work is the
order-three source cap needed by q_ZZ and the derivative of the non-linear
reduced-pressure error.  No source midpoint is selected and no ancestor owner
is rebuilt by ``run``; callers pass the retained checked pulse owner.
"""
import copy
import gzip
import hashlib
import json
import math
import time
from pathlib import Path

import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source as pulse
import lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source as end

HERE, PREFIX = pulse.HERE, pulse.PREFIX
sha, ep = pulse.sha, pulse.ep
record, rebox, source_precision = pulse.record, pulse.rebox, pulse.source_precision
NAME = PREFIX + 'current_original_n1_pulse_pressure_C1.json.gz'
RECEIPT = PREFIX + 'current_original_n1_pulse_pressure_C1_check.json'
GATE = 'current_original_n1_pulse_pressure_C1_source_installed'
LOCAL_GATE = 'current_original_n1_pulse_pressure_C1_local_grid_accuracy_installed'
OPEN = tuple(dict.fromkeys(pulse.OPEN + (
    'complete_signed_pressure_tails_certified',
    'positive_order_five_moment_repair_closed',
    'temporal_recursion',
    'uniform_full_pressure_derivative_accuracy_certified',
)))
CELLS = pulse.CELLS
SAMPLES = ('-1', '-.639', '0', '.211', '.371', '.639', '1')
C1_GRID_COUNT = 4096
SOURCE_DPS = 540


def _read_json_gz(path):
    return json.loads(gzip.decompress(Path(path).read_bytes()))


def _interval(c, value):
    if hasattr(value, '_mpi_'):
        return c.make_mpf(value._mpi_)
    if isinstance(value, dict) and 'lower_exact_mpf_tuple' in value:
        return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                           tuple(value['upper_exact_mpf_tuple'])))
    raise TypeError('Directed interval required')


def _abs_upper(value):
    lo, hi = ep(value)
    return max(abs(lo), abs(hi))


def _scalar_bound(c, value):
    return c.mpf([0, ep(value)[1]])


def _json_safe(value):
    """Recursively encode live interval leaves while preserving metadata."""
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item) for item in value]
    if hasattr(value, '_mpi_'):
        return record(value)
    if hasattr(value, 'coefficients') and hasattr(value, 'order'):
        return record(value)
    if isinstance(value, mp.mpf):
        # A few inherited metadata fields are ordinary mpf scalars rather
        # than directed interval leaves; preserve them textually in JSON.
        return mp.nstr(value, n=max(40, mp.mp.dps))
    return value


def _jet_row(c, value, index):
    if isinstance(value, dict) and 'coefficients' in value:
        return _interval(c, value['coefficients'][index])
    return _interval(c, value[index])


def _source_hashes(owner):
    hashes = dict(owner.hashes)
    parent_raw = _read_json_gz(pulse.HERE / pulse.NAME)
    parent_receipt = json.loads((pulse.HERE / pulse.RECEIPT).read_bytes())
    if parent_raw.get('source_family') != owner.family:
        raise ValueError('Accepted pulse raw has a different source family')
    if parent_receipt.get('source_family') != owner.family or not parent_receipt.get('all_passed'):
        raise ValueError('Accepted pulse receipt is not complete')
    if not parent_receipt.get(pulse.GATE):
        raise ValueError('Accepted pulse gate is not installed')
    for name, digest in parent_raw.get('input_hashes', {}).items():
        if sha(name) != digest:
            raise ValueError('Changed accepted pulse input ' + name)
        hashes[name] = digest
    for name, digest in parent_receipt.get('input_hashes', {}).items():
        if sha(name) != digest:
            raise ValueError('Changed accepted pulse receipt input ' + name)
        hashes[name] = digest
    for name in (pulse.NAME, pulse.RECEIPT, Path(pulse.__file__).name,
                 Path(end.__file__).name, Path(__file__).name):
        hashes[name] = sha(name)
    return hashes, parent_raw


def _order_three_source_caps(owner):
    """Recover B_j and T_j (j=0..3) from the actual C5 source rows."""
    c = owner.c
    if not hasattr(owner, '_cells') or len(owner._cells) != len(CELLS):
        raise ValueError('Retained four-cell pulse source rows are required')
    _, ubar, _, _ = end.polynomials(c, ['-1', '1'], 5, owner.delta)
    normalization = owner.coefficients.basis.normalization
    beta = c.exp(-1) / (c.mpf('.15') * normalization)
    B = [c.mpf(0) for _ in range(4)]
    T = [c.mpf(0) for _ in range(4)]
    cell_records = []
    for cell in owner._cells:
        if tuple(cell['ends']) not in CELLS:
            raise ValueError('Unexpected pulse coefficient cell')
        packet = owner.coefficients.evaluate(list(cell['ends']))
        view = owner.coefficients.report(packet)
        first = end.read_jet(c, view['source_inputs']['incoming_moment_Taylor'][0])
        if first.order < 5:
            raise ValueError('Incoming original mean is not a C5 source')
        initial = ubar * first
        main = ubar * cell['ap'] * 11
        end_terms = [ubar * coefficient * beta for coefficient in cell['C']]
        for j in range(4):
            b = _abs_upper(main[j]) + sum(_abs_upper(value[j]) for value in end_terms)
            t = _abs_upper(initial[j])
            B[j] = c.mpf(max(ep(B[j])[1], math.factorial(j) * b))
            T[j] = c.mpf(max(ep(T[j])[1], math.factorial(j) * t))
        cell_records.append(dict(
            ends=list(cell['ends']),
            source_order=5,
            incoming_mean_order=first.order,
            B=[record(value) for value in B],
            Tinitial=[record(value) for value in T],
        ))
    if any(ep(value)[1] <= 0 for value in B + T):
        raise ValueError('Order-three source cap is empty')
    return B, T, cell_records


def _energy_log_caps(c, owner, B, T):
    mu = owner.mu
    logE = rebox(c, owner.coefficients.selected.amplitude.log_end_scale)
    logA0 = owner.end.logP + owner.end.logUp
    reference = owner.reference
    length = 13 / mu
    common = logA0 - logE - 13 / (2 * mu) + 13 - c.ln(reference)
    rows = []
    for b, t in zip(B, T):
        row = dict(
            unweighted_U=common + 2 * c.ln(b) - c.ln(2 * mu),
            unweighted_initial=common + 2 * c.ln(t),
            weighted_U_with_length_squared=common - 26 + c.ln(2)
                + 2 * c.ln(b) + c.ln(1 + length ** 2),
            weighted_initial_with_length=logA0 - logE - 39 / (2 * mu) + 13
                + 2 * c.ln(t) + c.ln(1 + length) - c.ln(reference),
        )
        if any(ep(value)[1] >= -1000 for value in row.values()):
            raise ValueError('Order-three pulse energy cap is not below -1000')
        rows.append(row)
    epsilon = c.exp(-1000)
    # The q_ZZ coefficient bound is obtained by differentiating the exact
    # radial recovery twice.  With delta<1/2 and L>=1/2:
    # |q_ZZ| <= 56|U|+24|U_Z|+4|U_ZZ|
    #             +28|m|+38|m_Z|+16|m_ZZ|+2|m_ZZZ|.
    q_energy = 288 * epsilon
    qz_energy = 7960 * epsilon
    # Seven-term Cauchy is retained here because the seven source components
    # do not share an accepted correlation enclosure.
    qzz_energy = 200256 * epsilon
    boundary_q_squared = 16 * epsilon
    boundary_qz_squared = 400 * epsilon
    return dict(
        log_caps=rows, length=length, epsilon=epsilon,
        U_derivative_bounds=B, initial_mean_derivative_bounds=T,
        source_energy_bounds_nonzero=True,
        source_order_three_retained=True,
        unweighted_Hardy_mean_factor=10,
        weighted_Cauchy_mean_factor=3,
        conservative_common_mean_factor=10,
        q_energy=q_energy, qZ_energy=qz_energy, qZZ_energy=qzz_energy,
        qZZ_majorant_coefficients=[56, 24, 4, 28, 38, 16, 2],
        qZZ_cauchy_factor=7,
        boundary_q_squared_bound=boundary_q_squared,
        boundary_qZ_squared_bound=boundary_qz_squared,
        exact_units_combined_before_boxing=True,
    )


def _quadratic_C1_bounds(c, owner, energy, qbase, qbase_Z):
    """Differentiate the reduced-S quadratic majorants in fixed units."""
    eps = energy['epsilon']
    U = eps
    UZ = eps
    UZZ = eps
    q = energy['q_energy']
    qZ = energy['qZ_energy']
    qZZ = energy['qZZ_energy']
    qbase2 = energy['boundary_q_squared_bound']
    qbaseZ2 = energy['boundary_qZ_squared_bound']
    cross_boundary = c.sqrt(qbase2 * qbaseZ2)
    cross_Uq = c.sqrt(UZ * qZ) + c.sqrt(U * qZZ)
    cross_UUZ = c.sqrt(U * UZ)
    # These are energy bounds: E(U_Z^2)<=epsilon and
    # E(U*U_ZZ)<=sqrt(E(U)E(U_ZZ))=epsilon.  They are not products of
    # scalar epsilon caps.
    derivative_UUZ = UZ + c.sqrt(U * UZZ)
    envelopes = _reduced_coefficient_envelopes(c)
    coefficient_majorant = c.mpf('256')
    if max(envelopes.values()) > ep(coefficient_majorant)[1]:
        raise ValueError('Reduced-S coefficient envelope exceeded proof majorant')
    # The source-unit values below are the accepted C0 term bounds.  The
    # common 256 majorant covers both each coefficient and its Z derivative;
    # the exact coefficient records are retained for independent checking.
    value0 = c.exp(4) * 8 * eps + cross_Uq + 2 * U + U
    value1 = 72 * eps + cross_Uq + 2 * U + U
    deriv0 = c.exp(4) * cross_boundary + cross_Uq + derivative_UUZ + derivative_UUZ
    deriv1 = cross_boundary + c.mpf('.5') * c.sqrt(q * qZ) \
        + cross_Uq + derivative_UUZ + derivative_UUZ
    S0_Z = coefficient_majorant * (value0 + deriv0)
    S1_Z = coefficient_majorant * (value1 + deriv1)
    reference = owner.reference
    A_cap = owner.end.A_cap
    quadratic = dict(
        m13=_scalar_bound(c, A_cap * S0_Z / (2 * reference)),
        m15=_scalar_bound(c, A_cap * S1_Z / (2 * reference)),
        reduced_S0_Z_absolute_bound=S0_Z,
        reduced_S1_Z_absolute_bound=S1_Z,
        reduced_S_coefficient_majorant=coefficient_majorant,
        reduced_S_coefficient_envelopes=envelopes,
        qbase_bound=qbase,
        qbase_Z_bound=qbase_Z,
        qZZ_energy=qZZ,
        all_quadratic_sources_retained=True,
    )
    viscous_Z = _scalar_bound(
        c, owner.invRv * qbase_Z / owner.reference)
    full = dict(m13=quadratic['m13'], m15=quadratic['m15'] + viscous_Z)
    return quadratic, viscous_Z, full


def _reduced_coefficient_envelopes(c):
    """Directed envelopes for the exact reduced-S coefficients and Z rows."""
    z = c.mpf([-1, 1])
    delta = c.mpf([0, '.5'])
    a = 1 - delta
    d = 1 - z * z
    L = 1 - delta * z * z
    Lz = -2 * delta * z
    dz = -2 * z
    coeff = dict(
        invL=1 / L,
        az_over_2L=a * z / (2 * L),
        d_over_L=d / L,
        two_z2_over_L2=2 * z * z / (L * L),
        two_delta_z2_over_L2=2 * delta * z * z / (L * L),
        two_zd_over_L2=2 * z * d / (L * L),
    )
    deriv = dict(
        invL= -Lz / (L * L),
        az_over_2L=a * (L - z * Lz) / (2 * L * L),
        d_over_L=(dz * L - d * Lz) / (L * L),
        two_z2_over_L2=4 * z / (L * L) - 4 * z * z * Lz / (L * L * L),
        two_delta_z2_over_L2=4 * delta * z / (L * L)
            - 4 * delta * z * z * Lz / (L * L * L),
        two_zd_over_L2=2 * (d + z * dz) / (L * L)
            - 4 * z * d * Lz / (L * L * L),
    )
    return {name: max(_abs_upper(coeff[name]), _abs_upper(deriv[name]))
            for name in coeff}


class CurrentOriginalN1PulsePressureC1Source:
    """Checked C1 source adapter; a retained pulse owner is mandatory."""

    @source_precision
    def __init__(self, pulse_owner=None, require_checked=True):
        if type(pulse_owner) is not pulse.CurrentOriginalN1PulsePressureSource:
            raise ValueError('A retained checked pulse owner is required')
        if not pulse_owner.acceptance_loaded:
            raise ValueError('Checked original pulse owner is required')
        pulse_owner.assert_graph()
        self.pulse = pulse_owner
        self.c = pulse_owner.c
        self.family = copy.deepcopy(pulse_owner.family)
        self.hashes, self.parent_raw = _source_hashes(pulse_owner)
        self.B, self.Tinitial, self.source_cells = _order_three_source_caps(pulse_owner)
        self.energy = _energy_log_caps(self.c, pulse_owner, self.B, self.Tinitial)
        qbase = c_qbase = self.c.mpf(0)
        qbase_Z = self.c.mpf(0)
        for cell in pulse_owner._cells:
            terms = pulse_owner._terms(list(cell['ends']), cell)
            qbase = self.c.mpf(max(ep(qbase)[1], _abs_upper(terms['qinit'][0])))
            qbase_Z = self.c.mpf(max(ep(qbase_Z)[1], _abs_upper(terms['qinit'][1])))
        parent_qbase = pulse_owner.bounds.get('qinit_bound')
        if parent_qbase is not None and ep(parent_qbase)[1] < ep(qbase)[1]:
            raise ValueError('Endpoint qbase source cap is inconsistent')
        self.qbase = qbase
        self.qbase_Z = qbase_Z
        self.quadratic, self.viscous_Z, self.full_errors = _quadratic_C1_bounds(
            self.c, pulse_owner, self.energy, qbase, qbase_Z)
        self.units = copy.deepcopy(pulse_owner.units)
        self.units.update(
            C1_source_order=3,
            qZZ_source_order=3,
            physical_lower_endpoint_derivative='-qinit_Z*invRv/reference for m15',
            unweighted_lower_endpoint_derivative_exact_zero=True,
            full_pressure_derivative_accuracy_not_claimed=True,
        )
        self.source_bindings = dict(
            order_three_source_caps='actual C5 coefficient/incoming moment Taylor rows',
            exact_linear_jets='current_original_n1_pulse_pressure_source._terms',
            qZZ_recovery='exact radial recovery differentiated twice; coefficients [56,24,4,28,38,16,2]',
            quadratic_C1='_quadratic_C1_bounds',
            weighted_viscous_C1='-qinit[1]*invRv/reference',
            local_C1_grid='actual pulse._make_view at every center and whole subcell',
        )
        self.models = self._extend_models(self.parent_raw.get('models', []))
        self.acceptance_loaded = False
        self.local_accuracy_loaded = False
        if require_checked:
            raw = _read_json_gz(HERE / NAME)
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if (not receipt.get('all_passed') or not receipt.get(GATE) or
                    not receipt.get(LOCAL_GATE) or
                    receipt.get('source_family') != self.family or
                    any(receipt.get(key, False) for key in OPEN)):
                raise ValueError('Checked pulse C1 receipt required')
            for name, digest in receipt.get('input_hashes', {}).items():
                if sha(name) != digest:
                    raise ValueError('Changed pulse C1 input ' + name)
                self.hashes[name] = digest
            if raw.get('source_family') != self.family:
                raise ValueError('Pulse C1 raw source family changed')
            if len(raw.get('c1_grid', [])) != C1_GRID_COUNT:
                raise ValueError('Checked local C1 grid has wrong coverage')
            self.c1_grid = raw['c1_grid']
            self.acceptance_loaded = True
            self.local_accuracy_loaded = True
        else:
            self.c1_grid = self._make_c1_grid()
            # The fine batch uses only the retained pulse source's private
            # view builder; certify that its graph stayed unchanged across it.
            self.pulse.assert_graph()
        self._packets = {}
        self._state = record((self.family, self.units, self.energy, self.quadratic,
                              self.full_errors, self.source_cells, self.models,
                              self.c1_grid))
        self._hashes = copy.deepcopy(self.hashes)

    def _extend_models(self, models):
        if len(models) != 256:
            raise ValueError('Accepted pulse model archive must contain 256 cells')
        out = []
        c = self.c
        for model in models:
            row = copy.deepcopy(model)
            half = _interval(c, model['half_width'])
            h = ep(half)[1]
            for key in ('m13', 'm15'):
                fields = row['fields'][key]
                curvature = _interval(c, fields['linear_Z2_half_range'])
                slope = _interval(c, fields['center_slope'])
                remainder = 2 * _abs_upper(curvature) * h
                ambiguity = (ep(slope)[1] - ep(slope)[0]) / 2
                source = ep(self.full_errors[key])[1]
                total = c.mpf([0, remainder + ambiguity + source])
                fields.update(
                    first_derivative_source_error=record(self.full_errors[key]),
                    first_derivative_linear_Taylor_remainder=record(c.mpf(remainder)),
                    first_derivative_coefficient_ambiguity=record(c.mpf(ambiguity)),
                    first_derivative_absolute_error=record(total),
                    first_derivative_accuracy_certified=False,
                    first_derivative_scope='C1 enclosure; uniform derivative gate remains false',
                )
            out.append(row)
        return out

    def _make_c1_grid(self):
        """Evaluate an independent fine C1 grid from the actual pulse views."""
        c = self.c
        count = C1_GRID_COUNT
        grid = []
        for index in range(count):
            left = c.mpf(-1) + c.mpf(2) * index / count
            right = c.mpf(-1) + c.mpf(2) * (index + 1) / count
            center = (left + right) / 2
            half = (right - left) / 2
            point = self.pulse._make_view(center)
            whole = self.pulse._make_view([ep(left)[0], ep(right)[1]])
            fields = {}
            for key in ('m13', 'm15'):
                center_jet = point['linear_pressure_Taylor'][key]
                whole_jet = whole['linear_pressure_Taylor'][key]
                slope = center_jet[1]
                z2 = whole_jet[2]
                # The grid approximates the first derivative by its center
                # value.  Taylor's bound is variation <= 2*|p''|*h across
                # [center-h, center+h], not a value-remainder h^2 bound.
                remainder = c.mpf(2 * _abs_upper(z2) * half)
                ambiguity = (ep(slope)[1] - ep(slope)[0]) / 2
                source = c.mpf(ep(self.full_errors[key])[1])
                total = c.mpf([0, remainder + ambiguity + source])
                fields[key] = dict(
                    center_slope=record(slope),
                    whole_Z2_bound=record(z2),
                    half_width=record(half),
                    source_error=record(self.full_errors[key]),
                    slope_ambiguity=record(c.mpf(ambiguity)),
                    total_error=record(total),
                    absolute_accuracy_certified=ep(total)[1] < ep(c.mpf('.001'))[0],
                )
            grid.append(dict(index=index, left=record(left), right=record(right),
                             center=record(center), half_width=record(half),
                             source_coefficient_cell=point['source_coefficient_cell'],
                             evaluation='pulse._make_view(center), pulse._make_view([left,right])',
                             fields=fields))
        return grid

    @source_precision
    def assert_graph(self):
        self.pulse.assert_graph()
        if self.hashes != self._hashes:
            raise ValueError('Pulse C1 source hashes changed')
        if self._state != record((self.family, self.units, self.energy,
                                  self.quadratic, self.full_errors,
                                  self.source_cells, self.models, self.c1_grid)):
            raise ValueError('Pulse C1 source state changed')

    @source_precision
    def evaluate(self, Z):
        self.assert_graph()
        cell = self.pulse._cell(Z)
        terms = self.pulse._terms(Z, cell)
        linear = {key: terms[key + '_linear'] / self.pulse.reference
                  for key in ('m13', 'm15')}
        qinit = terms['qinit']
        viscous = -qinit[0] * self.pulse.invRv / self.pulse.reference
        viscous_Z = -qinit[1] * self.pulse.invRv / self.pulse.reference
        derivative = {}
        full = {}
        for key in ('m13', 'm15'):
            derivative[key] = linear[key][1]
            full[key] = derivative[key] + self.c.mpf(
                [-ep(self.full_errors[key])[1], ep(self.full_errors[key])[1]])
        full['m15'] += viscous_Z
        packet = object()
        view = dict(
            Z=self.c.mpf(Z), source_family=self.family,
            source_coefficient_cell=list(cell['ends']),
            linear_pressure_Taylor=linear,
            full_pressure_first_derivative=full,
            quadratic_C1_source_errors=self.full_errors,
            source_linear_inputs=terms,
            weighted_viscous_endpoint=viscous,
            weighted_viscous_endpoint_derivative=viscous_Z,
            unweighted_viscous_boundary_exact_zero=True,
            lower_endpoint_derivative_sign='-qinit_Z*invRv/reference',
            order_three_energy_bounds=self.energy,
            qZZ_source_bound=self.energy['qZZ_energy'],
            exact_linear_jets_reused=True,
            finite_mu_endpoint_retained=True,
            full_pressure_derivative_accuracy_not_claimed=True,
            **{GATE: self.acceptance_loaded},
            **{LOCAL_GATE: self.local_accuracy_loaded},
            **dict.fromkeys(OPEN, False),
        )
        self._packets[id(packet)] = (packet, view, record(view))
        return packet

    @source_precision
    def report(self, packet):
        self.assert_graph()
        item = self._packets.get(id(packet))
        if item is None or item[0] is not packet or item[2] != record(item[1]):
            raise ValueError('Live pulse C1 packet required')
        # mpmath interval constants are intentionally not deepcopyable.  Use
        # the accepted end-source rebox walker so public reports remain live,
        # directed views at the owner's source context.
        return end.copy_view(self.c, item[1])


@source_precision
def run(pulse_owner=None):
    began = time.monotonic()
    owner = CurrentOriginalN1PulsePressureC1Source(pulse_owner, require_checked=False)
    samples = {z: owner.report(owner.evaluate(z)) for z in SAMPLES}
    raw = dict(
        source_family=owner.family,
        parent_archive=dict(name=pulse.NAME, source_family=owner.family,
                            parent_input_hashes=owner.parent_raw.get('input_hashes', {})),
        source_cells=_json_safe(owner.source_cells),
        order_three_energy=_json_safe(owner.energy),
        quadratic_C1=_json_safe(owner.quadratic),
        full_errors=_json_safe(owner.full_errors),
        qbase=_json_safe(owner.qbase), qbase_Z=_json_safe(owner.qbase_Z),
        units=_json_safe(owner.units), source_bindings=_json_safe(owner.source_bindings),
        samples=_json_safe(samples), models=_json_safe(owner.models),
        c1_grid=_json_safe(owner.c1_grid), c1_grid_count=C1_GRID_COUNT,
        input_hashes=owner.hashes,
        **{GATE: False, LOCAL_GATE: False}, **dict.fromkeys(OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (HERE / NAME).write_bytes(gzip.compress(
        (json.dumps(raw, separators=(',', ':')) + '\n').encode(), mtime=0))
    print('ORIGINAL_N1_PULSE_PRESSURE_C1_COMPUTED', len(samples), len(owner.models), flush=True)
    return owner, samples, owner.models


if __name__ == '__main__':
    raise SystemExit('Pass a retained checked pulse owner to run()')
