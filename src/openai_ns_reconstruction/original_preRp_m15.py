"""Independent repaired weighted m15 transport from Rc through the true Rp.

The signed source functions stay bound to their directed C0/C1 enclosures.
Neither enclosure midpoints nor a terminal-derived inlet drive this flow.
"""
import copy
from fractions import Fraction
import gzip
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_compact_n1_pressure import weighted_beta_prefix
from .original_quiet_numeric import OriginalQuietNumeric


class OriginalPreRpM15Numeric:
    def __init__(self, *, quiet=None, **kwargs):
        if quiet is not None and kwargs:
            raise ValueError('Supply an existing quiet owner or constructor options')
        self.quiet = quiet if quiet is not None else OriginalQuietNumeric(**kwargs)
        self.reader = p = self.quiet.compact.pressure
        self.reserved = p.reserved
        self.ctx, self.iv, self.J = p.ctx, p.iv, p.J
        self.family, self.directory = p.family, p.directory
        self.cells, self.ell, self.beta_jets = p.cells, p.ell, p.beta_jets
        self.mu, self.logk = p.mu, p.logk
        quiet_source = p.inputs['quiet_source']
        if (not quiet_source['source_scaled_Z1_exact_zero']
                or quiet_source['m15']['sign'] != 'negative'
                or quiet_source['m15']['radial_scale_power'] != 2):
            raise ValueError('Original negative p=2 quiet history and fixed K required')
        self.hashes, self.code_hashes = dict(self.quiet.hashes), dict(self.quiet.code_hashes)
        self._pressure_integrals, self._packets = {}, {}
        stem = PREFIX+'compliant_current_original_n1_terminal_integral_identities'
        receipt_path = self.directory/(stem+'_check.json')
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        gates = ('current_original_n1_actual_terminal_integrals_reduce_to_reserved_response',
                 'n1_repaired_own_terminal_five_function_zeros')
        checks = ('actual_integrals', 'transported_partition', 'selected_inverse')
        if (receipt['source_family'] != self.family or not receipt['all_passed']
                or not all(receipt[k] for k in gates)
                or not all(receipt['checks'][k]['passed'] for k in checks)):
            raise ValueError('Same-source signed response and transported partition proofs required')
        for name, digest in receipt['input_hashes'].items():
            if self.hashes.get(name, digest) != digest:
                raise ValueError('Conflicting repaired m15 source: '+name)
            if name not in self.hashes and hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed repaired m15 source: '+name)
            self.hashes[name] = digest
        path = self.directory/(stem+'.json.gz')
        raw = json.loads(gzip.decompress(path.read_bytes()))
        self.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.hashes[receipt_path.name] = hashlib.sha256(receipt_path.read_bytes()).hexdigest()
        terminal = raw['terminal_identity_bindings']
        if (raw['source_family'] != self.family
                or terminal['terminal_sources']['m15'] != p.binding['terminal']['m15']):
            raise ValueError('Actual transported m15 and proved terminal function must agree')
        self.join_identity = dict(source=stem+'.json.gz',
            terminal_function=terminal['terminal_sources']['m15'],
            patch_response=terminal['reduced_patch']['m15'],
            full_defect_plus_selected_response=terminal['reduced_terminal']['m15'],
            checks={k: copy.deepcopy(receipt['checks'][k]) for k in checks},
            exact_source_identity='inner+lower+active+compact+own_swirl+quiet+pulse=0',
            actual_Rp_source_equals_negative_complete_pulse=True,
            enclosure_overlap_is_not_the_identity=True,
            independently_resolved_signed_numeric_seam=False)
        name = 'src/openai_ns_reconstruction/original_preRp_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self, values): return self.reader._rows(values)

    def radius_atoms(self, packet): return self.reader.radius_atoms(packet)

    def fixed_unit_atoms(self, packet):
        """Original log(Phi*A*Rv^2), with shared powers combined first."""
        iv = self.iv
        lp = self.reader.background.patch.leading.defect.logP
        lm = self.reader.background.outer.logmui
        phi = iv.mpf('10.99')*iv.sqrt(2*iv.pi)*(iv.exp(5)-iv.exp(3))
        rest = [iv.log(phi)+lp+self.reserved.logUp-17,
                iv.mpf('18.5')*iv.exp(-lm),
                -3*iv.exp(-iv.mpf(2)/3*(iv.log(4)+lm)),
                -iv.mpf(364)/3*lm, iv.mpf(2)/3*iv.log(4)-iv.log(6)/2]
        return [self.reader._row(self.reserved._source_row(v)*2) for v in self.radius_atoms(packet)]+self._rows(rest)

    def leading_unit_atoms(self, packet):
        return [self.reader._row(self.reserved._source_row(v)*2) for v in self.radius_atoms(packet)]+self._rows([2*self.logk-self.iv.log(4)])

    def swirl_prefix(self, s, index):
        """Positive int_0^s exp((1/2-mu)u)*shape_j(u)du.

        The extra R in -R*theta0*Utheta1 ds changes the m13 rate.
        Shape_j includes exp(-center_j/2); its ell cancels the Jacobian.
        """
        c, iv = self.ctx, self.iv
        center = self.reserved.centers[index]
        rate = iv.mpf('.5')-self.mu
        lo, hi = _bounds(c, (s-center)/self.ell)
        if hi <= -1: return iv.mpf(0)
        lo, hi = max(c.mpf(-1), lo), min(c.mpf(1), hi)
        lower = weighted_beta_prefix(self, lo, rate)
        upper = weighted_beta_prefix(self, hi, rate)
        value = iv.mpf([_bounds(c, lower)[0], _bounds(c, upper)[1]])
        return iv.exp(-self.mu*center)*value/self.reserved.normal

    def own_swirl(self, cells, s):
        iv, decode = self.iv, self.reserved._decode
        integrals = [self.swirl_prefix(s, j) for j in (1, 2, 3)]
        values = []
        for cell in cells:
            a = [decode(v) for v in cell['actual_A_rc_C0_C1_coefficient_enclosures']]
            pair = [iv.mpf(0), iv.mpf(0)]
            for name, integral in zip(('beta1', 'beta2', 'beta3'), integrals):
                b = [decode(v) for v in cell['coefficients'][name]['normalized_C0_C1_enclosures']]
                pair[0] -= a[0]*b[0]*integral
                pair[1] -= (a[1]*b[0]+a[0]*b[1])*integral
            values.append(pair)
        return [self.reader._hull([v[k] for v in values]) for k in (0, 1)], integrals

    def _inlet(self, Z):
        if Z not in self.quiet.inlets:
            self.quiet.inlets[Z] = self.quiet.compact.evaluate(Z, '2')
        return self.quiet.inlets[Z]

    def evaluate(self, Z, chart, coordinate):
        c, iv = self.ctx, self.iv
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:
            raise ValueError('Interior original decimal Z required')
        if chart not in ('compact', 'quiet') or type(coordinate) is not str:
            raise ValueError('Choose compact x in[1,2] or quiet fraction q in[0,1]')
        if not c.isfinite(c.mpf(coordinate)):
            raise ValueError('Finite exact decimal coordinate required')
        q = Fraction(coordinate)
        if not (1<=q<=2 if chart == 'compact' else 0<=q<=1):
            raise ValueError('Original native chart domain required')
        key = Z, chart, coordinate
        if key in self._packets: return self._packets[key]
        packet = self._inlet(Z)
        cells = self.reader._matched(Z)
        if not cells: raise ValueError('Original axial source cover required')
        decode = self.reserved._decode
        unit, terms = self.fixed_unit_atoms(packet), []

        def append(name, values, atoms, binding):
            terms.append(dict(name=name, coefficient_C0_C1=self._rows(values),
                log_scale_terms=copy.deepcopy(atoms), actual_function_binding=binding))

        parts = self.reader.full['pressure_component_enclosures']['m15']
        for name in ('inner', 'lower'):
            append(name, [decode(v) for v in parts[name]['C1_enclosures']], unit,
                dict(source=PREFIX+'compliant_current_original_n1_full_five_numeric_bridge.json.gz',
                     component=name, ordinary_signed_source_function_enclosures=True))
        ends = {tuple(v['original_source_Z_cell']) for v in cells}
        active = [v for v in self.reader.pre['cell_piece_enclosures'] if tuple(v['exact_Z_cell']) in ends]
        append('active', [self.reader._hull([decode(v['piece_enclosures']['m15']['active'][k]) for v in active]) for k in (0, 1)],
            unit, dict(source=PREFIX+'compliant_current_original_n1_preRp_numeric_bridge.json.gz', piece='active',
                       signed_m15_already_applied_once=True))
        # These positive endpoint/integral majorants cover every compact
        # prefix: x<=2, B(x^2*f)<=5 sup|f|, int x dx<=2, and both q_y ends.
        logs = [[], []]
        for end in sorted(ends):
            label = '|'.join(end)
            error = self.reader.band[label]['view']['perturbation_bounds']['error_m15']
            first = self.reader.c1['original_target_and_control_C3_update']['compact'][label]['normalized_first_Z_upper']['m15']
            ratio = decode(self.reader.inputs['scale_ratios']['m15'])
            for k, row in enumerate((error, first)):
                if row['exact_zero']: raise ValueError('Actual compact finite-N sectors must remain nonzero')
                value = decode(row['log_absolute_upper'])+ratio
                if _bounds(c, value)[1] >= -1000:
                    raise ArithmeticError('Compact weighted prefix exceeds original fixed-unit budget')
                logs[k].append(value)
        cap = _bounds(c, iv.exp(-1000))[1]
        empty = chart == 'compact' and q == 1
        finite = [iv.mpf(0)]*2 if empty else [iv.mpf([-cap, cap])]*2
        append('compact_finite_N', finite, unit,
            dict(function=self.reader.binding['base_pressure_prefix']['m15'],
                 component='finite_N_increment_to_exact_compact_weighted_prefix',
                 boundary_x_squared_factor=5, integral_x_dx_factor=2,
                 finite_RQ_y_and_RQ_yZ_endpoints_retained=True, ordinary_Z1_source_bound=True))
        compact_s = iv.log(iv.mpf(coordinate)) if chart == 'compact' else iv.log(2)
        append('compact_exact_leading', [-compact_s, iv.mpf(0)], self.leading_unit_atoms(packet),
            dict(function=self.reader.binding['base_pressure_prefix']['m15'],
                 formula='-K^2*log(R/Rc)/4', K='Rc*k, original fixed Z-independent scalar',
                 ordinary_Z1_exact_zero=True))
        own, integrals = self.own_swirl(cells, compact_s)
        # Original Utheta1 unit + logP + logSigma + logRc equals this
        # fixed p=2 unit exactly. Cancel before boxing huge log atoms.
        append('own_swirl', own, unit,
            dict(function=self.reader.binding['patch_histories']['m15'],
                 measure='-R*theta0*Utheta1 ds', rate='1/2-mu',
                 fixed_unit_identity='Utheta1_unit+logP+logSigma+logRc=pressure_unit_m15',
                 original_A_rc_Z1_and_beta_Z1_product_rule=True))
        quiet_s = iv.mpf(0) if chart == 'compact' else self.quiet.lengthi*iv.mpf(coordinate)
        if chart == 'quiet':
            append('quiet_exact_leading', [-quiet_s, iv.mpf(0)], self.leading_unit_atoms(packet),
                dict(source=PREFIX+'compliant_current_original_n1_quiet_pressure_source.json.gz',
                     function='actual_quiet_source.m15; partial prefix in the original log radius',
                     component='weighted quiet history from unchanged K',
                     formula='-K^2*s/4; s=log(R/(2Rc))', ordinary_Z1_exact_zero=True))
        out = dict(Z=Z, chart=chart, coordinate=coordinate, source_family=self.family,
            selected_integer=packet['selected_integer'], m15_C0_C1_terms=terms,
            fixed_m15_unit_atoms=unit, swirl_partial_integral_factors=self._rows(integrals),
            actual_partial_finite_N_normalized_log_bounds=[self._rows(v) for v in logs],
            compact_cumulative_true_function=self.reader.binding['cumulative']['m15'],
            incoming_true_function=self.reader.binding['incoming']['m15'],
            source_partition_identity=self.join_identity,
            quiet_log_radius=self._rows([quiet_s])[0],
            contract=dict(inner_lower_active_incoming_retained_once=True,
                actual_selected_reserved_coefficient_functions_enclosed=True,
                compact_partial_RdR_measure_and_finite_endpoints_retained=True,
                original_negative_weight_and_extra_R_applied_once=True,
                exact_nonzero_K_squared_compact_and_quiet_histories_retained=True,
                original_fixed_p2_unit_cancellation_before_boxing=True,
                actual_A_rc_and_beta_ordinary_Z1_product_rule=True,
                original_positive_mu_and_selected_N_retained=True,
                independent_of_terminal_derived_numeric_inlet=True,
                enclosure_centers_not_selected_as_signed_source_values=True,
                no_hyper_exponential_materialization=True),
            limits=dict(independent_preRp_m15_C0_C1_enclosure_transport=True,
                preRp_pulse_signed_source_identity=True,
                resolved_signed_point_centers=False, independent_signed_numeric_seam_residual=False,
                relative_cancellation_accuracy=False, global_five_moment_numeric_join=False,
                common_global_N=False, exact_heat_exterior=False, stress_cone=False,
                genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out
