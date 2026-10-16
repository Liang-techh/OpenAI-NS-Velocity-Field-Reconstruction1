"""Original cumulative pressure through constant-U charts to exact R110."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_switch_m15 import OriginalFirstSwitchM15, frozen_coefficients, density_difference
from .original_macro_m15 import add, proxy_modes
from .original_r110_pressure_rows import OriginalR110PressureRows


class OriginalR110M15(OriginalFirstSwitchM15):
    def __init__(self, *, rows=None, first=None, second=None, **kwargs):
        if rows is not None and kwargs:
            raise ValueError('Supply R110 pressure source or options')
        self.rows = rows if rows is not None else OriginalR110PressureRows(**kwargs)
        if not isinstance(self.rows, OriginalR110PressureRows):
            raise ValueError('Genuine original R110 pressure source required')
        self.first = first if first is not None else OriginalFirstSwitchM15(rows=self.rows.first_rows)
        for name in ('t', 'o', 'ctx', 'iv', 'J', 'audit', 'directory', 'family'):
            setattr(self, name, getattr(self.first, name))
        self.hashes, self.code_hashes = dict(self.first.hashes), dict(self.first.code_hashes)
        self.hashes.update(self.rows.hashes); self.code_hashes.update(self.rows.code_hashes)
        path = self.directory/(PREFIX+'compliant_current_original_first_switch_m15_numeric.json')
        raw = path.read_bytes(); data = json.loads(raw)
        if (not data['checks']['passed'] or data['source_family'] != self.family
                or not data['limits']['first_switch_cumulative_signed_m15_with_complete_factored_errors']):
            raise ValueError('Accepted complete original first-switch m15 required')
        for base, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                     (ROOT, data['implementation_hashes'], self.code_hashes)):
            for name, digest in group.items():
                if target.get(name, digest) != digest or hashlib.sha256((base/name).read_bytes()).hexdigest() != digest:
                    raise ValueError('Changed original R110 moment prerequisite: '+name)
                target[name] = digest
        self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        self.first_inlets = {Fraction(p['Z']): p for p in data['samples'] if p['coordinate']['value'] == '1'}
        self.chart, self.second = self.rows.chart, second
        if self.chart == 'post_power' and (second is None or second.chart != 'second_switch'):
            raise ValueError('Explicit cumulative second-switch moment inlet required')
        self._packets = {}; self._leading = {}; self._post = None
        name = 'src/openai_ns_reconstruction/original_r110_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def post_power_field(self):
        if self.chart != 'second_switch':
            raise ValueError('Post-power cumulative source requires second-switch endpoint')
        if self._post is None:
            self._post = OriginalR110M15(rows=self.rows.post_power_field(), first=self.first, second=self)
        return self._post

    def inlet(self, Z):
        if self.chart == 'post_power':
            return self.second.evaluate(Z, '1')
        key = Fraction(Z)
        if key not in self.first_inlets:
            self.first_inlets[key] = self.first.evaluate(Z, '1')
            self.hashes.update(self.first.hashes); self.code_hashes.update(self.first.code_hashes)
        return self.first_inlets[key]

    def leading(self, Z, q):
        key = Z, q._mpi_
        if key in self._leading:
            return self._leading[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        post = self.chart == 'post_power'
        first_inputs = self.first.rows.inputs(Z, iv.mpf(1))
        macro = first_inputs['incoming']; core = macro['incoming']['inlet']
        point, truth, tm, qa = {}, {}, None, None
        for ctx, values, de, axis, target in ((c, core['point']['values'], o.core.delta, o.core.j, point),
                                             (iv, core['truth']['values'], o.core.delta_iv, o.core.j_iv, truth)):
            modes, a = proxy_modes(J, ctx, Z, de, axis, values['Psi'], values['Mpsi'])
            target.update(frozen_coefficients(J, ctx, Z, de, modes))
            if ctx is iv: tm, qa = modes, a
        L = self.rows.micro.post_log
        maximum_eta = L*q+2*o.hcap if post else o.hcap*(1+q)
        decay = iv.mpf([_bounds(self.audit, iv.exp(-maximum_eta))[0], 1])*iv.exp(-self.rows.micro.Y)
        base = dict(U=tm['U'], Q={(0, 0): tm['Q'][(0, 0)],
                    (1, 0): tm['Q'][(1, 0)]+tm['Q'][(1, 1)]*decay}, Uy={}, Qy={(1, 0): -qa*decay})
        bp = self.first.macro.uniform_source_per_hb(Z)
        values = core['truth']['values']
        bm = (macro['incoming']['delta']['Mpsi']+3*bp
              +o._norm(values['Mpsi']-values['Psi'])*2*iv.exp(2*o.hcap))
        P = self.first.source_norms(Z, iv.mpf(1))['primitive']
        z = J.variable(iv, iv.mpf(Z), 3); A, D = 1-z*z*o.core.delta_iv, 1-z*z
        def errors(pnorm, mnorm):
            p, m = (o._norm_error_jet(norm, 3) for norm in (pnorm, mnorm))
            qrow = (p*(2*z)-m*(z*(1-o.core.delta_iv))-o._derivative(m)*D)/A
            qy = (-((p-m)*(z*(1-o.core.delta_iv)))-o._derivative(p-m)*D)/A
            return {name: {(1, 0): row.truncate(2)} for name, row in
                    dict(U=p, Q=qrow, Uy=J.constant(iv, 0, 3), Qy=qy).items()}
        eb = errors(bp, bm)
        mean_source = P*(2*o.hcap+self.rows.geometry(q)['mass']) if post else P*o.hcap*(1+q)
        es = errors(P, mean_source)
        kb = density_difference(J, iv, Z, o.core.delta_iv, base, eb, o.hcap)
        broadened = {name: add(base[name], {k: row*o.hcap for k, row in eb[name].items()}) for name in base}
        ks = density_difference(J, iv, Z, o.core.delta_iv, broadened, es, o.hcap**2)
        weight = L*q*iv.exp(2*(L*q+2*o.hcap)) if post else q*iv.exp(2*o.hcap*(1+q))
        def collect(poly):
            result = {}
            for (degree, mode), row in poly.items():
                if mode != 0:
                    raise ValueError('Error baseline must first enclose all inherited mean modes')
                term = self.symmetric(row)*weight
                result[degree] = result[degree]+term if degree in result else term
            return result
        result = dict(point=point, truth=truth, inlet_errors=collect(kb), source_errors=collect(ks),
                      normalized_inlet_bp=bp, normalized_inlet_bm=bm, normalized_first_switch_P=P)
        self._leading[key] = result
        return result

    def finite_errors(self, Z, q, packet):
        _, base, geo = self.rows.background(Z, q)
        owner = self.rows.second_rows if self.chart == 'post_power' else self.rows.first_rows
        lowpacket = owner.evaluate(Z, '1'); _, lowbase, lowgeo = owner.background(Z, self.iv.mpf(1))
        up = self.correction_densities(Z, packet, base)
        low = self.correction_densities(Z, lowpacket, lowbase)
        volume = self.correction_densities(Z, packet, self.rows.uniform_background(Z, q))['volume']
        result = {}
        for name, power in (('endpoint', 2), ('viscous', 1)):
            result[name] = [self.rows.sum_logs((None if up[name][k] is None else up[name][k]+power*geo['logR'],
                None if low[name][k] is None else low[name][k]+power*lowgeo['logR'])) for k in (0, 1)]
        post, iv = self.chart == 'post_power', self.iv
        logmass = iv.log(self.rows.micro.post_log*q) if post else iv.log(q)
        uppereta = self.rows.micro.post_log*q+2*self.o.hcap if post else self.o.hcap*(1+q)
        result['volume'] = [None if v is None else v+2*iv.log(100)+2*uppereta+logmass for v in volume]
        return result, dict(endpoint=self.rows.log_rows(low['endpoint']), viscous=self.rows.log_rows(low['viscous']),
            log_radius=self.o._summary(lowgeo['logR_point'], lowgeo['logR']))

    def evaluate(self, Z, coordinate):
        packet = self.rows.evaluate(Z, coordinate); q, _, recipe = self.rows.micro._coordinate(coordinate, False)
        key = Z, coordinate
        if key in self._packets: return self._packets[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        post = self.chart == 'post_power'; incoming = self.inlet(Z)
        expected = 69 if post else 51
        if len(incoming['m15_C0_C1_terms']) != expected or incoming['actual_selected_integer'] != packet['selected_integer']:
            raise ValueError('Complete same-N original constant-U inlet required')
        terms = copy.deepcopy(incoming['m15_C0_C1_terms'])
        empty = q._mpi_ == iv.mpf(0)._mpi_; zero, zeroi = J.constant(c, 0, 1), J.constant(iv, 0, 1)
        if not empty: lead = self.leading(Z, q); finite, lower = self.finite_errors(Z, q, packet)
        else: lead = None; finite = {name: [None, None] for name in ('endpoint', 'volume', 'viscous')}; lower = None
        qp = c.mpf(coordinate)
        def width(name):
            atom = o._summary(sum(_bounds(c, o.logh))/2, o.logh)
            atom.update(source_log_atom_only=True, original_native_parameter_name=name,
                        distinct_graph_nodes_with_checked_common_numeric_source=True)
            return atom
        hbB, hbS = width('hbB'), width('hbS')
        def physical(degree, mode=0):
            return o._summary(2*c.log(100)-degree*o.core.logLambda-mode*self.rows.micro.Ypoint,
                              2*iv.log(100)-degree*o.core.ll_iv-mode*self.rows.micro.Y)
        def append(name, p, t, atoms, **binding):
            terms.append(dict(name=self.chart+'_'+name,
                coefficient_C0_C1=o._jet_summary(p.truncate(1), t.truncate(1)), log_scale_terms=copy.deepcopy(atoms),
                actual_function_binding=dict(chart=self.chart, incoming_terms_preserved_once=expected,
                    epsilon_and_mean_modes_separate=True, **binding)))
        Lp, Lt = self.rows.micro.post_log_point, self.rows.micro.post_log
        for d in range(3):
            for n in range(d+1):
                k = 2-n
                if post:
                    jp = Lp*qp if k == 0 else c.expm1(k*Lp*qp)/k
                    jt = Lt*q if k == 0 else iv.expm1(k*Lt*q)/k
                    atoms = [physical(d, n)]
                else: jp, jt, atoms = qp, q, [physical(d, n), hbS]
                p, t = (zero, zeroi) if empty else (lead['point'][(d, n)]*jp, lead['truth'][(d, n)]*jt)
                append(f'half_S1_epsilon_{d}_mean_{n}', p, t, atoms, inherited_exp_minus_nY_atom_retained=True)
                if post or k:
                    if empty: error = zeroi
                    elif post:
                        error = self.symmetric(lead['truth'][(d, n)])*jt*(2*k+2/Lt)*iv.exp(2*k*o.hcap)
                    else:
                        error = self.symmetric(lead['truth'][(d, n)])*k*q*(1+q)*iv.exp(k*o.hcap*(1+q))
                    append(f'geometry_error_epsilon_{d}_mean_{n}', zero, error,
                        [physical(d, n), hbS] if post else [physical(d, n), hbS, hbS], error_enclosure_only=True)
        for d in (1, 2):
            for kind, atoms in (('inlet', [hbB] if post else [hbS, hbB]),
                                ('source', [hbS, hbS] if post else [hbS, hbS, hbS])):
                append(f'{kind}_error_epsilon_{d}', zero, zeroi if empty else lead[kind+'_errors'].get(d, zeroi),
                       [physical(d)]+atoms, error_enclosure_only=True, errors_before_native_measure=True)
        for name in ('endpoint', 'volume', 'viscous'):
            atoms = [hbS] if name == 'volume' and not post else []; truth = zeroi
            bound = self.rows.sum_logs(finite[name])
            if bound is not None:
                top = iv.mpf(_bounds(self.audit, bound)[1]); atom = o._summary(_bounds(c, top)[1], top)
                atom.update(positive_error_bound_atom_only=True, source_value_selected=False)
                atoms += [atom]; truth = J(iv, [iv.mpf([-1, 1]), iv.mpf([-1, 1])])
            append('finite_N_'+name+'_error', zero, truth, atoms, error_enclosure_only=True,
                   both_nonflat_endpoints=True, prefix_background_uniform=True)
        out = dict(Z=Z, coordinate=recipe, chart=self.chart, source_family=self.family,
            actual_selected_integer=copy.deepcopy(packet['selected_integer']), m15_C0_C1_terms=terms,
            incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']), exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']), full_prefix_cover=copy.deepcopy(packet['source_cover']),
            finite_N_error_logs_before_native_measure={name: self.rows.log_rows(v) for name, v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower, constant_U_viscous_increment_exact_zero=True,
            contract=dict(complete_inlet_retained_once=True, original_selected_N_unchanged=True,
                original_constant_U_with_nonzero_mean=True, actual_background_R_Qy_increment_exact_zero=True,
                individual_Qy_not_zeroed=True, both_finite_N_endpoints=True, all_quadratic_error_products=True,
                inherited_core_mean_atoms_not_rounded=True, native_measure_before_width_errors=True,
                original_bridge_and_switch_width_nodes_distinct=True, uniform_volume_background=True,
                source_error_two_width_factors_retained=True, physical_R110_exact_if_post_endpoint=True),
            limits=dict(second_switch_post_power_cumulative_signed_m15_with_complete_factored_errors=True,
                resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
                independent_preRp_signed_seam=False, global_five_moment_numeric_join=False, common_global_N=False,
                exact_heat_exterior=False, stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out
