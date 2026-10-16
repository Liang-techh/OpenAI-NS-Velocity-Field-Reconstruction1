"""Original reserved-band P1 as source-bounded physical exponential sums.

Incoming pressure, fixed leading N1p and the new swirl integral have their
original units. Their very different scales are never flattened into a
float or a single midpoint. C1 means ordinary slow-Z differentiation.
"""
import copy
from fractions import Fraction
import gzip
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_compact_n1_velocity import OriginalCompactN1Velocity


def weighted_beta_prefix(engine, endpoint, rate):
    """Positive int_-1^endpoint beta(t) exp(rate*ell*t) dt.

    A degree-three midpoint polynomial and a full-cell fourth-derivative
    interval enclose the defining integral, including the flat endpoints.
    """
    c, iv, J = engine.ctx, engine.iv, engine.J
    endpoint = max(c.mpf(-1), min(c.mpf(1), endpoint))
    rate = iv.mpf(rate)
    key = tuple(endpoint._mpf_), rate._mpi_
    if key in engine._pressure_integrals:
        return engine._pressure_integrals[key]
    if endpoint == -1:
        return iv.mpf(0)
    width = (iv.mpf(endpoint) + 1) / engine.cells
    k = rate * engine.ell
    total = iv.mpf(0)
    for i in range(engine.cells):
        left = -1 + width * i
        right = -1 + width * (i + 1)
        mid = -1 + width * (iv.mpf(i) + iv.mpf('.5'))
        box = iv.mpf([max(c.mpf(-1), _bounds(c, left)[0]),
                      min(c.mpf(1), _bounds(c, right)[1])])

        def jet(q):
            beta = J(iv, engine.beta_jets(iv, q))
            ex = iv.exp(k * q)
            return beta * J(iv, [ex * k**j / math.factorial(j) for j in range(5)])

        at, whole = jet(mid), jet(box)
        total += at[0] * width + at[2] * width**3 / 12 + whole[4] * width**5 / 80
    lo, hi = _bounds(c, total)
    result = iv.mpf([max(c.mpf(0), lo), max(c.mpf(0), hi)])
    engine._pressure_integrals[key] = result
    return result


class OriginalReservedN1PressureNumeric:
    def __init__(self, velocity):
        self.velocity = velocity
        self.background, self.reserved = velocity.numeric, velocity.reserved
        r = self.reserved
        self.ctx, self.iv, self.J = r.ctx, r.iv, r.J
        self.family, self.directory = r.family, r.directory
        self.cells, self.ell, self.beta_jets = r.cells, r.ell, r.beta_jets
        self.hashes, self.code_hashes = dict(velocity.hashes), dict(velocity.code_hashes)
        self._pressure_integrals, self._packets = {}, {}
        specifications = (
            ('current_original_n1_full_five_numeric_bridge', 'current_original_true_complete_five_defect_C1_numeric_bridge_installed'),
            ('current_original_n1_preRp_numeric_bridge', 'current_original_true_preRp_pressure_C1_numeric_bridge_installed'),
            ('current_original_n1_preRp_pressure_C1', 'current_original_n1_preRp_pressure_C1_installed'),
            ('current_original_n1_preRp_C1_inputs', 'current_original_n1_preRp_C1_inputs_installed'),
            ('current_original_n1_repair_band_pressure_source', 'current_original_n1_repair_band_pressure_source_installed'),
            ('current_original_n1_reserved_own_fields', 'current_original_n1_reserved_own_C1_histories_and_C0_QV_installed'))
        raw = {}
        for stem, gate in specifications:
            receipt_path = self.directory / (PREFIX + 'compliant_' + stem + '_check.json')
            receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not receipt[gate] or receipt['source_family'] != self.family:
                raise ValueError('Checked same-source P1 input required: ' + stem)
            for name, digest in receipt['input_hashes'].items():
                if (self.hashes.get(name, digest) != digest or
                        hashlib.sha256((self.directory / name).read_bytes()).hexdigest() != digest):
                    raise ValueError('Changed/conflicting P1 input: ' + name)
                self.hashes[name] = digest
            self.hashes[receipt_path.name] = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            path = self.directory / (PREFIX + 'compliant_' + stem + '.json.gz')
            self.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            raw[stem] = json.loads(gzip.decompress(path.read_bytes()))
            if raw[stem]['source_family'] != self.family:
                raise ValueError('Same original P1 family required')
        self.full = raw['current_original_n1_full_five_numeric_bridge']['samples']['whole_Z']
        self.pre = raw['current_original_n1_preRp_numeric_bridge']['samples']['whole_Z']
        self.c1 = raw['current_original_n1_preRp_pressure_C1']['source_data']
        self.inputs = raw['current_original_n1_preRp_C1_inputs']
        self.band = raw['current_original_n1_repair_band_pressure_source']['samples']
        self.binding = raw['current_original_n1_reserved_own_fields']['own_field_bindings']
        if not self.full['pressure_normalization'] == self.pre['normalization']:
            raise ValueError('Incoming and fixed pressure units differ')
        self.logk = r._decode(self.inputs['quiet_source']['K']['log_coefficient'])
        self.mu = self.background.outer.mui
        name = 'src/openai_ns_reconstruction/original_compact_n1_pressure.py'
        self.code_hashes[name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()

    def _row(self, value):
        return self.reserved._rows([value])[0]

    def _rows(self, values):
        return self.reserved._rows(values)

    def _hull(self, values):
        c, iv = self.ctx, self.iv
        return iv.mpf([min(_bounds(c, v)[0] for v in values),
                       max(_bounds(c, v)[1] for v in values)])

    def _matched(self, Z):
        c, iv = self.ctx, self.iv
        zl, zh = _bounds(c, iv.mpf(Z))
        return [v for v in self.reserved.coefficients['cells']
                if not (_bounds(c, self.reserved._decode(v['clipped_Z']))[1] < zl or
                        _bounds(c, self.reserved._decode(v['clipped_Z']))[0] > zh)]

    def radius_atoms(self, packet):
        return copy.deepcopy(packet['coordinate']['global_y_log_terms'][:4]) + [self._row(self.iv.log(110))]

    def pressure_unit_atoms(self, packet):
        """Exact fixed Phi*A*Rv unit, combined before interval arithmetic."""
        iv = self.iv
        lp = self.background.patch.leading.defect.logP
        lm = self.background.outer.logmui
        phi = iv.mpf('10.99') * iv.sqrt(2 * iv.pi) * (iv.exp(5) - iv.exp(3))
        rest = [iv.log(phi) + lp + self.reserved.logUp - 15,
                iv.mpf('5.5') * iv.exp(-lm),
                -3 * iv.exp(-iv.mpf(2) / 3 * (iv.log(4) + lm)),
                -iv.mpf(184) / 3 * lm, iv.mpf(2) / 3 * iv.log(4) - iv.log(6) / 2]
        return self.radius_atoms(packet) + self._rows(rest)

    def swirl_prefix(self, s, index):
        """int_0^s exp(-(.5+mu)u)*shape_j(u)du, a positive scalar."""
        c, iv = self.ctx, self.iv
        center = self.reserved.centers[index]
        alpha = iv.mpf('.5') + self.mu
        r = (s - center) / self.ell
        lo, hi = _bounds(c, r)
        if hi <= -1:
            return iv.mpf(0)
        lo, hi = max(c.mpf(-1), lo), min(c.mpf(1), hi)
        lower = weighted_beta_prefix(self, lo, -alpha)
        upper = weighted_beta_prefix(self, hi, -alpha)
        integral = iv.mpf([_bounds(c, lower)[0], _bounds(c, upper)[1]])
        return iv.exp(-(alpha + iv.mpf('.5')) * center) * integral / self.reserved.normal

    def evaluate(self, Z, x, packet=None, reserved=None):
        key = Z, x
        if key in self._packets:
            return self._packets[key]
        packet = packet or self.background.evaluate(Z, x)
        reserved = reserved or self.reserved.evaluate(Z, x, packet)
        c, iv = self.ctx, self.iv
        cells = self._matched(Z)
        if not cells:
            raise ValueError('Original P1 axial cover required')
        decode = self.reserved._decode
        parts = self.full['pressure_component_enclosures']['m13']
        unit = self.pressure_unit_atoms(packet)
        terms = []

        def append(name, values, atoms, binding):
            terms.append(dict(name=name, coefficient_C0_C1=self._rows(values),
                              log_scale_terms=copy.deepcopy(atoms), actual_function_binding=binding))

        for name in ('inner', 'lower'):
            append(name, [decode(v) for v in parts[name]['C1_enclosures']], unit,
                   dict(source=PREFIX + 'compliant_current_original_n1_full_five_numeric_bridge.json.gz',
                        component=name, actual_ordinary_signed_source_sectors_retained=True))
        ends = {tuple(v['original_source_Z_cell']) for v in cells}
        active = [v for v in self.pre['cell_piece_enclosures'] if tuple(v['exact_Z_cell']) in ends]
        append('active', [self._hull([decode(v['piece_enclosures']['m13']['active'][k]) for v in active]) for k in (0, 1)],
               unit, dict(source=PREFIX + 'compliant_current_original_n1_preRp_numeric_bridge.json.gz', piece='active'))
        # Both error majorants use absolute endpoint/integral inequalities
        # uniform on 1<=x<=2, so the same bounds apply to every partial prefix.
        logs = [[], []]
        for end in sorted(ends):
            label = '|'.join(end)
            view = self.band[label]['view']
            err = view['perturbation_bounds']['error_m13']
            first = self.c1['original_target_and_control_C3_update']['compact'][label]['normalized_first_Z_upper']['m13']
            ratio = decode(self.inputs['scale_ratios']['m13'])  # log(Rc / fixed pressure unit), Rc cancelled
            for k, row in enumerate((err, first)):
                if row['exact_zero']:
                    raise ValueError('Original finite-N pressure bound must retain nonzero sectors')
                log = decode(row['log_absolute_upper']) + ratio
                if _bounds(c, log)[1] >= -1000:
                    raise ArithmeticError('Partial compact P1 source exceeds fixed-unit budget')
                logs[k].append(log)
        cap = _bounds(c, iv.exp(-1000))[1]
        finite = [iv.mpf(0), iv.mpf(0)] if Fraction(x) == 1 else [iv.mpf([-cap, cap])] * 2
        append('compact_finite_N', finite, unit,
               dict(source=PREFIX + 'compliant_current_original_n1_reserved_own_fields.json.gz',
                    function=self.binding['base_pressure_prefix']['m13'], partial_prefix=True,
                    component='finite_N_increment_to_fixed_leading_prefix',
                    C0_absolute_endpoint_factors=[3, 5], C1_absolute_endpoint_factors=[3, 5],
                    finite_viscous_endpoints_retained=True, Rc_cancelled_before_boxing=True))
        s = iv.log(iv.mpf(x))
        append('compact_exact_leading', [-iv.expm1(-s), iv.mpf(0)],
               self.radius_atoms(packet) + self._rows([2 * self.logk - iv.log(4)]),
               dict(source='original radial-free pressure identity', q0='-k/x', k='physical K/Rc',
                    formula='Rc*k^2/4*(1-exp(-s))', ordinary_Z1_exact_zero=True))
        integrals = [self.swirl_prefix(s, j) for j in (1, 2, 3)]
        patch_cells, scale_logs = [], []
        for cell in cells:
            a = [decode(v) for v in cell['actual_A_rc_C0_C1_coefficient_enclosures']]
            values = [iv.mpf(0), iv.mpf(0)]
            for name, integral in zip(('beta1', 'beta2', 'beta3'), integrals):
                b = [decode(v) for v in cell['coefficients'][name]['normalized_C0_C1_enclosures']]
                values[0] += a[0] * b[0] * integral
                values[1] += (a[1] * b[0] + a[0] * b[1]) * integral
            patch_cells.append(values)
            scale_logs.append(decode(cell['actual_A_rc_shared_fixed_log_scale']))
        # This source-fixed common scale is independent of Z. The derivatives
        # of the actual angular amplitude remain explicitly in a[1] above.
        if any(v._mpi_ != scale_logs[0]._mpi_ for v in scale_logs):
            raise ValueError('One fixed original angular amplitude scale required')
        append('own_swirl', [self._hull([v[k] for v in patch_cells]) for k in (0, 1)],
               copy.deepcopy(reserved['Utheta1_log_unit_atoms']) +
               self._rows([self.background.patch.leading.defect.logP, scale_logs[0]]),
               dict(source=PREFIX + 'compliant_current_original_n1_reserved_own_fields.json.gz',
                    function=self.binding['patch_histories']['m13'], measure='theta0*Utheta1 ds',
                    leading_theta0='Pstar*A_rc*exp(-(.5+mu)s)', ordinary_product_rule=True))
        out = dict(Z=Z, native_compact_x=x, source_family=self.family, P1_C0_C1_terms=terms,
                   fixed_pressure_unit_atoms=unit, swirl_partial_integral_factors=self._rows(integrals),
                   actual_partial_finite_N_normalized_log_bounds=[self._rows(v) for v in logs],
                   cumulative_true_function=self.binding['cumulative']['m13'],
                   incoming_true_function=self.binding['incoming']['m13'],
                   contract=dict(core_datum_lower_active_retained_once=True, fixed_leading_N1p_retained=True,
                                 pressure_not_rebuilt_from_n1_patch=True, physical_Pstar_applied_once=True,
                                 physical_dR_equals_R_ds_applied_once=True, disjoint_original_supports_used=True,
                                 original_positive_mu_retained=True, ordinary_amplitude_Z1_retained=True,
                                 zero_centers_are_proxies_not_signed_source_values=True,
                                 no_hyper_exponential_materialization=True),
                   limits=dict(P1_values_and_ordinary_Z1=True, resolved_signed_function_centers=False,
                               relative_pressure_cancellation_accuracy=False, pressure_radial_derivatives=False,
                               terminal_pressure_zero_at_2Rc=False, global_n1_joins=False,
                               beyond_2Rc=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out


class OriginalCompactN1PressureVelocity(OriginalCompactN1Velocity):
    """Physical orders0/1 velocity plus full reserved-band P1 enclosures."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.pressure = OriginalReservedN1PressureNumeric(self)
        self.hashes.update(self.pressure.hashes)
        self.code_hashes.update(self.pressure.code_hashes)

    def evaluate(self, Z, x, *, log_tau='-1', scaled_time=None, azimuth='0'):
        result = super().evaluate(Z, x, log_tau=log_tau, scaled_time=scaled_time, azimuth=azimuth)
        n1 = self.pressure.evaluate(Z, x, self.numeric.evaluate(Z, x), result['reserved_n1'])
        c = self.ctx
        ell, _, _ = OriginalCoreVelocity._time(self, Z, log_tau, scaled_time)
        time = [self._summary((-2 * v[0], -2 * v[1])) for v in ell]
        terms = [copy.deepcopy(result['physical']['leading_pressure'])]
        for source in n1['P1_C0_C1_terms']:
            terms.append(dict(name='n1_pressure_' + source['name'], temporal_order=1,
                              coefficient=copy.deepcopy(source['coefficient_C0_C1'][0]),
                              log_scale_terms=copy.deepcopy(time + source['log_scale_terms'])))
        result['physical'].update(pressure=terms, full_n1_pressure_included=True, n1_pressure_available=True,
                                  pressure_orders_included=[0, 1], n1_pressure_time_weight='lambda^(-2)',
                                  n1_pressure_extra_Lambda_squared_factor_applied=False,
                                  useful_combined_pressure_relative_error_certified=False)
        result['reserved_n1_pressure'] = n1
        result['limits'].update(n1_pressure=True, physical_pressure_leading_only=False,
                                compact_n1_pressure_C1=True, relative_pressure_cancellation_accuracy=False)
        return result
