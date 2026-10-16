"""Long m15 with correlated formal exponential pairs and exact finite T."""
import copy
from fractions import Fraction
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_switch_m15 import frozen_coefficients, density_difference
from .original_macro_m15 import add, proxy_modes
from .original_r110_m15 import OriginalR110M15
from .original_long_pressure_rows import OriginalLongPressureRows


class OriginalLongM15(OriginalR110M15):
    def __init__(self, *, rows=None, **kwargs):
        if rows is not None and kwargs:
            raise ValueError('Supply genuine long pressure rows or constructor options')
        self.rows = rows if rows is not None else OriginalLongPressureRows(**kwargs)
        if not isinstance(self.rows, OriginalLongPressureRows):
            raise ValueError('Actual original long pressure source required')
        self.r110 = OriginalR110M15(rows=self.rows.r110_rows).post_power_field()
        self.first = self.r110.first
        for name in ('t', 'o', 'ctx', 'iv', 'J', 'audit', 'directory', 'family'):
            setattr(self, name, getattr(self.r110, name))
        self.hashes, self.code_hashes = dict(self.r110.hashes), dict(self.r110.code_hashes)
        self.hashes.update(self.rows.hashes); self.code_hashes.update(self.rows.code_hashes)
        path = self.directory/(PREFIX+'compliant_current_original_r110_m15_numeric.json')
        raw = path.read_bytes(); data = json.loads(raw)
        if (not data['checks']['passed'] or data['source_family'] != self.family
                or not data['limits']['second_switch_post_power_cumulative_signed_m15_with_complete_factored_errors']):
            raise ValueError('Accepted complete original R110 pressure moment required')
        self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
        self.inlets = {Fraction(p['Z']): p for p in data['samples']
                       if p['chart'] == 'post_power' and p['coordinate']['value'] == '1'}
        self.chart = 'long_reshape'; self._packets = {}; self._leading = {}
        name = 'src/openai_ns_reconstruction/original_long_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self, Z):
        key = Fraction(Z)
        if key not in self.inlets:
            self.inlets[key] = self.r110.evaluate(Z, '1')
            self.hashes.update(self.r110.hashes); self.code_hashes.update(self.r110.code_hashes)
        return self.inlets[key]

    def leading(self, Z, q):
        key = Z, q._mpi_
        if key in self._leading: return self._leading[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        first_inputs = self.first.rows.inputs(Z, iv.mpf(1))
        macro = first_inputs['incoming']; core = macro['incoming']['inlet']
        point, truth, tm, qa = {}, {}, None, None
        for ctx, values, de, axis, target in ((c, core['point']['values'], o.core.delta, o.core.j, point),
                                             (iv, core['truth']['values'], o.core.delta_iv, o.core.j_iv, truth)):
            modes, a = proxy_modes(J, ctx, Z, de, axis, values['Psi'], values['Mpsi'])
            target.update(frozen_coefficients(J, ctx, Z, de, modes))
            if ctx is iv: tm, qa = modes, a
        E110 = iv.exp(-self.rows.micro.Y-self.rows.long.owner.post_log)
        # Y=logLambda+log25 is about9.4e17, rather than the enormous T.
        # This sparse mpf needs only a small integer exponent. Keep its
        # actual directed endpoints; an ordinary-float underflow is invalid.
        if any(v[1] == 0 or abs(v[2]).bit_length() > 128 for v in E110._mpi_):
            raise ArithmeticError('Inherited E110 mean mode must remain a nonzero sparse mpf')
        decay = E110*iv.mpf([0, 1])
        base = dict(U=tm['U'], Q={(0, 0): tm['Q'][(0, 0)],
            (1, 0): tm['Q'][(1, 0)]+tm['Q'][(1, 1)]*decay}, Uy={}, Qy={(1, 0): -qa*decay})
        bp = self.first.macro.uniform_source_per_hb(Z)
        values = core['truth']['values']
        bm = macro['incoming']['delta']['Mpsi']+4*bp+o._norm(values['Mpsi']-values['Psi'])*2*iv.exp(2*o.hcap)
        P = self.first.source_norms(Z, iv.mpf(1))['primitive']
        z = J.variable(iv, iv.mpf(Z), 3); A, D = 1-z*z*o.core.delta_iv, 1-z*z
        def errors(pnorm, mnorm):
            p, m = (o._norm_error_jet(v, 3) for v in (pnorm, mnorm))
            qrow = (p*(2*z)-m*(z*(1-o.core.delta_iv))-o._derivative(m)*D)/A
            qy = (-(p-m)*(z*(1-o.core.delta_iv))-o._derivative(p-m)*D)/A
            return {name: {(1, 0): row.truncate(2)} for name, row in
                dict(U=p, Q=qrow, Uy=J.constant(iv, 0, 3), Qy=qy).items()}
        eb = errors(bp, bm); es = errors(P, P*(2*o.hcap+2))
        kb = density_difference(J, iv, Z, o.core.delta_iv, base, eb, o.hcap)
        widened = {name: add(base[name], {k: row*o.hcap for k, row in eb[name].items()}) for name in base}
        ks = density_difference(J, iv, Z, o.core.delta_iv, widened, es, o.hcap**2)
        y = self.rows.geometry(q)['y']
        mass = iv.mpf(min(_bounds(self.audit, y)[1], self.audit.mpf('.5')))
        def collect(poly):
            result = {}
            for (d, n), row in poly.items():
                if n != 0: raise ValueError('Uniform error baseline must first enclose all mean modes')
                term = self.symmetric(row)*mass
                result[d] = result[d]+term if d in result else term
            return result
        result = dict(point=point, truth=truth, inlet_errors=collect(kb), source_errors=collect(ks),
                      native_positive_mass_bound=mass)
        self._leading[key] = result
        return result

    def finite_errors(self, Z, q, packet):
        _, base, geo = self.rows.background(Z, q)
        owner = self.rows.inlet_rows
        lowpacket = owner.evaluate(Z, '1'); _, lowbase, _ = owner.background(Z, self.iv.mpf(1))
        up = self.correction_densities(Z, packet, base)
        low = self.correction_densities(Z, lowpacket, lowbase)
        volume = self.correction_densities(Z, packet, self.rows.uniform_background(Z, q))['volume']
        iv = self.iv; result = {}
        for name, power in (('endpoint', 2), ('viscous', 1)):
            result[name] = [self.rows.sum_logs((None if up[name][k] is None else up[name][k]+power*geo['logR'],
                None if low[name][k] is None else low[name][k]+power*iv.log(110))) for k in (0, 1)]
        mass = iv.mpf(min(_bounds(self.audit, geo['y'])[1], self.audit.mpf('.5')))
        result['volume'] = [None if v is None else v+2*geo['logR']+iv.log(mass) for v in volume]
        return result, dict(endpoint=self.rows.log_rows(low['endpoint']), viscous=self.rows.log_rows(low['viscous']),
            log_radius=self.o._summary(self.ctx.log(110), iv.log(110)))

    def at_log_radius(self, Z, y):
        return self.evaluate(Z, 'y:'+y)

    def evaluate(self, Z, coordinate):
        packet = self.rows.evaluate(Z, coordinate); q, _, recipe = self.rows._coordinate(coordinate, False)
        key = Z, coordinate
        if key in self._packets: return self._packets[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        incoming = self.inlet(Z)
        if len(incoming['m15_C0_C1_terms']) != 88 or incoming['actual_selected_integer'] != packet['selected_integer']:
            raise ValueError('Complete same-N original R110 88-term inlet required')
        terms = copy.deepcopy(incoming['m15_C0_C1_terms']); geo = self.rows.geometry(q)
        yp, y = geo['y_point'], geo['y']; empty = y._mpi_ == iv.mpf(0)._mpi_
        zero, zeroi = J.constant(c, 0, 1), J.constant(iv, 0, 1)
        if not empty: lead = self.leading(Z, q); finite, lower = self.finite_errors(Z, q, packet)
        else: lead = None; finite = {name: [None, None] for name in ('endpoint', 'volume', 'viscous')}; lower = None
        def physical(d, n=0, growth=0):
            lp = self.rows.long.owner.post_log_point; lt = self.rows.long.owner.post_log
            return o._summary(2*c.log(110)-d*o.core.logLambda-n*(self.rows.micro.Ypoint+lp)+growth*yp,
                2*iv.log(110)-d*o.core.ll_iv-n*(self.rows.micro.Y+lt)+growth*y)
        def width(name):
            atom = o._summary(sum(_bounds(c, o.logh))/2, o.logh)
            atom.update(source_log_atom_only=True, original_native_parameter_name=name,
                        distinct_graph_nodes_with_checked_common_numeric_source=True)
            return atom
        hbB, hbS = width('hbB'), width('hbS')
        def append(name, p, t, atoms, **binding):
            terms.append(dict(name='long_'+name, coefficient_C0_C1=o._jet_summary(p.truncate(1), t.truncate(1)),
                log_scale_terms=copy.deepcopy(atoms), actual_function_binding=dict(chart='long_reshape',
                    incoming88_terms_preserved_once=True, original_finite_T_not_infinite_window=True, **binding)))
        for d in range(3):
            for n in range(d+1):
                k = 2-n; group = f'long_H_epsilon_{d}_mean_{n}'
                if k:
                    for sign, endpoint in ((1, 'upper'), (-1, 'lower')):
                        p, t = (zero, zeroi) if empty else (lead['point'][(d, n)]*(c.mpf(sign)/k), lead['truth'][(d, n)]*(iv.mpf(sign)/k))
                        append(f'half_S1_epsilon_{d}_mean_{n}_{endpoint}', p, t,
                            [physical(d, n, k if sign == 1 else 0)], correlated_coefficient_id=group,
                            correlated_coefficient_sign=sign, kernel=f'(exp({k}*y)-1)/{k}',
                            exponential_kept_in_log_atom=True, inherited_exp_minus_nY_and_R110_offset_retained=True)
                else:
                    atoms = [physical(d, n)]
                    if not empty:
                        atom = o._summary(c.log(yp), iv.log(y)); atom.update(positive_native_log_measure_only=True,
                            original_measure='y=T*q', huge_y_not_multiplied_into_coefficient=True); atoms.append(atom)
                    append(f'half_S1_epsilon_{d}_mean_{n}_resonance', zero if empty else lead['point'][(d, n)],
                           zeroi if empty else lead['truth'][(d, n)], atoms, exact_resonant_kernel_y=True)
        for d in (1, 2):
            for kind, atoms in (('inlet', [hbB]), ('source', [hbS, hbS])):
                append(f'{kind}_error_epsilon_{d}', zero, zeroi if empty else lead[kind+'_errors'].get(d, zeroi),
                       [physical(d, growth=2)]+atoms, error_enclosure_only=True,
                       finite_positive_mass_min_y_half=True, native_dy_measure_already_integrated=True)
        for name in ('endpoint', 'volume', 'viscous'):
            bound = self.rows.sum_logs(finite[name]); atoms = []; truth = zeroi
            if bound is not None:
                top = iv.mpf(_bounds(self.audit, bound)[1]); atom = o._summary(_bounds(c, top)[1], top)
                atom.update(positive_error_bound_atom_only=True, source_value_selected=False)
                atoms.append(atom); truth = J(iv, [iv.mpf([-1, 1]), iv.mpf([-1, 1])])
            append('finite_N_'+name+'_error', zero, truth, atoms, error_enclosure_only=True,
                   both_nonflat_endpoints=True, finite_native_long_measure_included_once=True)
        out = dict(Z=Z, coordinate=recipe, source_family=self.family, chart='long_reshape',
            m15_C0_C1_terms=terms, incoming_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']),
            actual_selected_integer=copy.deepcopy(packet['selected_integer']), exact_empty_increment=empty,
            geometry=copy.deepcopy(packet['geometry']), full_prefix_cover=copy.deepcopy(packet['source_cover']),
            finite_N_error_logs_in_physical_units={name: self.rows.log_rows(v) for name, v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower, constant_U_viscous_increment_exact_zero=True,
            contract=dict(complete_R110_88term_inlet_once=True, original_T_singleton_preserved=True,
                inherited_mean_modes_exp_minus_nY_plus_log_11_over_10=True,
                huge_exp_y_and_radius_never_materialized=True, signed_exponential_pairs_correlated=True,
                resonance_y_is_separate_formal_atom=True, actual_finite_interval_pressure_primitive=True,
                both_finite_N_pressure_endpoints=True, all_density_quadratic_errors=True,
                uniform_prefix_background_and_correct_T_measure=True, exact_background_RQy_increment_zero=True,
                mean_Z3_QZZ_Uy_Qy_Z1_errors_complete=True),
            limits=dict(long_cumulative_signed_m15_with_complete_factored_errors=True,
                resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
                later_long_Rm_pressure=False, independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
                stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out
