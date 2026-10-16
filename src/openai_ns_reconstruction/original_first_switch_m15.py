"""Original first-switch cumulative pressure with native width atoms.

Core mean modes keep their exp(-nY) atoms instead of subtracting rounded
R100 values. Source-error densities preserve both bridge and switch widths.
"""
import copy
from fractions import Fraction
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_switch_pressure_rows import OriginalFirstSwitchPressureRows
from .original_macro_m15 import (OriginalMacroM15, add, scale, multiply, dz, proxy_modes,
    half_S1_coefficients, densities, combine, difference_densities, difference_products)


def density_kernel(J, ctx, Z, delta, fields):
    U, Q, Uy, Qy = (fields[k] for k in ('U', 'Q', 'Uy', 'Qy'))
    bc, vc = half_S1_coefficients(J, ctx, Z, delta)
    values = densities(U, Q, J)
    return add(scale(combine(bc, values), 2), combine(vc, values), scale(Qy, bc['q']),
               scale(multiply(Q, Qy), bc['q2']*2), scale(multiply(U, Uy), bc['U2']*2))


def density_difference(J, ctx, Z, delta, fields, error, hcap):
    U, Q, Uy, Qy = (fields[k] for k in ('U', 'Q', 'Uy', 'Qy'))
    eu, eq, euy, eqy = (error[k] for k in ('U', 'Q', 'Uy', 'Qy'))
    bc, vc = half_S1_coefficients(J, ctx, Z, delta)
    values = difference_densities(U, Q, eu, eq, J, hcap)
    return add(scale(combine(bc, values), 2), combine(vc, values), scale(eqy, bc['q']),
               scale(difference_products(Q, Qy, eq, eqy, hcap), bc['q2']*2),
               scale(difference_products(U, Uy, eu, euy, hcap), bc['U2']*2))


def frozen_coefficients(J, ctx, Z, delta, modes):
    bc, vc = half_S1_coefficients(J, ctx, Z, delta)
    values = densities(modes['U'], modes['Q'], J)
    F, G = combine(bc, values), combine(vc, values)
    zero = J.constant(ctx, 0, 1)
    return {(d, n): F.get((d, n), zero)*(2-n)+G.get((d, n), zero)
            for d in range(3) for n in range(d+1)}


class OriginalFirstSwitchM15(OriginalMacroM15):
    def __init__(self, *, rows=None, **kwargs):
        if rows is not None and kwargs:
            raise ValueError('Supply retained first-switch source or options')
        self.rows = rows if rows is not None else OriginalFirstSwitchPressureRows(**kwargs)
        if not isinstance(self.rows, OriginalFirstSwitchPressureRows):
            raise ValueError('Actual original first-switch pressure rows required')
        self.macro = OriginalMacroM15(rows=self.rows.macro_rows)
        for name in ('t', 'o', 'ctx', 'iv', 'J', 'audit', 'directory', 'family'):
            setattr(self, name, getattr(self.macro, name))
        self.hashes, self.code_hashes = dict(self.macro.hashes), dict(self.macro.code_hashes)
        self.hashes.update(self.rows.hashes); self.code_hashes.update(self.rows.code_hashes)
        reports = {}
        for suffix in ('macro_m15_numeric', 'first_switch_pressure_rows'):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if not data['checks']['passed'] or data['source_family'] != self.family:
                raise ValueError('Accepted original first-switch moment prerequisite required')
            for root, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                        (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if target.get(name, digest) != digest or hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
                        raise ValueError('Changed original first-switch moment source: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest(); reports[suffix] = data
        if not reports['macro_m15_numeric']['limits']['macro_cumulative_leading_signed_m15_with_complete_factored_errors']:
            raise ValueError('Complete original macro moment inlet required')
        self.inlet_report = PREFIX+'compliant_current_original_macro_m15_numeric.json'
        self.inlets = {Fraction(p['Z']): p for p in reports['macro_m15_numeric']['samples']
                       if p['coordinate'] == dict(kind='exact_decimal_macro_fraction', value='1')}
        self._packets = {}; self._leading = {}; self._source_norms = {}
        name = 'src/openai_ns_reconstruction/original_first_switch_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self, Z):
        key = Fraction(Z)
        if key not in self.inlets:
            self.inlets[key] = self.macro.evaluate(Z, '1')
            self.hashes.update(self.macro.hashes); self.code_hashes.update(self.macro.code_hashes)
        p = self.inlets[key]
        if p['source_family'] != self.family or len(p['m15_C0_C1_terms']) != 28:
            raise ValueError('Same-family original twenty-eight-sector macro inlet required')
        return p

    def source_norms(self, Z, s):
        """Before-width source value, primitive and R-weight variation norms."""
        key = Z, s._mpi_
        if key in self._source_norms:
            return self._source_norms[key]
        o, iv, b = self.o, self.iv, self.o.bridge
        p = o._prepared[Z]
        module = __import__(PREFIX+'compliant_actual_bridge_integrals', fromlist=['jet_norm'])
        phase = b.ctx.mpf(_bounds(self.audit, s)[1])
        with mp.workdps(280):
            mass, _ = b.pulse_masses(phase)
            def modal(name, power, radius_weight=0, differentiated=False):
                total = b.ctx.mpf(0)
                for j, row in enumerate(p['modes'][name]):
                    rate = power-j+radius_weight
                    total += (module.jet_norm(row, b.weight)*module.upper(b.ctx, abs(p['R0']**j*b.ctx.mpf(100)**(power-j)))
                              *b.ctx.exp(abs(rate)*b.cap*phase)*(abs(rate) if differentiated else 1))
                return total
            angular = o._transfer(modal('D_over_R', 1)*phase/2)
            ell = self.rows.micro._switch_increments(Z, s)[0]
            macro_ell = o._transfer(6*b.cap)
            ratio = o._transfer(module.jet_norm(p['quotient'], b.weight))*iv.exp(macro_ell+ell)
            value, weighted, weighted_derivative = iv.mpf(0), iv.mpf(0), iv.mpf(0)
            for name, power, logscale in (('hydro', 1, b.ctx.mpf(0)), ('pressure', 1, b.pressure_log),
                                           ('swirl', 2, b.swirl_log)):
                term = b.source_term(p['data']['phi0']*0, 1, logscale)
                for dest in ('value', 'weighted', 'derivative'):
                    norm = o._transfer(modal('drive_'+name, power, int(dest != 'value'), dest == 'derivative'))*ratio
                    cap = o._abs(o._scaled(norm, term, divide_eps=True))
                    if dest == 'value': value += cap
                    elif dest == 'weighted': weighted += cap
                    else: weighted_derivative += cap
        m = o._transfer(mass)
        # R-weighted source differs from its R100 value by hbS times
        # radial exponent variation, plus hbS^2 angular variation.
        variation = s*m*weighted_derivative+o.hcap*angular*iv.exp(ell)*m*weighted
        result = dict(value=value, primitive=value*m, mass=m,
                      weighted_primitive_remainder_per_hb=variation, angular_per_hb_squared=angular)
        self._source_norms[key] = result
        return result

    def leading(self, Z, s):
        key = Z, s._mpi_
        if key in self._leading:
            return self._leading[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        v = self.rows.inputs(Z, s); macro = v['incoming']
        core = macro['incoming']['inlet']
        point, truth = {}, {}
        raw_modes = {}
        for ctx, values, de, j, target in ((c, core['point']['values'], o.core.delta, o.core.j, point),
                                          (iv, core['truth']['values'], o.core.delta_iv, o.core.j_iv, truth)):
            modes, qa = proxy_modes(J, ctx, Z, de, j, values['Psi'], values['Mpsi'])
            target.update(frozen_coefficients(J, ctx, Z, de, modes)); raw_modes[ctx is iv] = modes, qa
        tm, qa = raw_modes[True]
        E100 = iv.exp(-self.rows.micro.Y)
        Eshort = iv.mpf([_bounds(self.audit, iv.exp(-o.hcap*s))[0], 1])
        base = dict(U=tm['U'], Q=dict(tm['Q']), Uy={}, Qy={(1, 0): -qa*E100*Eshort})
        base['Q'] = {(0, 0): tm['Q'][(0, 0)],
                     (1, 0): tm['Q'][(1, 0)]+tm['Q'][(1, 1)]*E100*Eshort}
        bp = self.macro.uniform_source_per_hb(Z)
        core_values = core['truth']['values']
        bm = (macro['incoming']['delta']['Mpsi']+bp
              +o._norm(core_values['Mpsi']-core_values['Psi'])*2*iv.exp(2*o.hcap)+bp)
        norms = self.source_norms(Z, s)
        z = J.variable(iv, iv.mpf(Z), 3); D, A = 1-z*z, 1-z*z*o.core.delta_iv
        def errors(pnorm, mnorm, ynorm):
            p, m, y = (o._norm_error_jet(norm, 3) for norm in (pnorm, mnorm, ynorm))
            q = (p*(2*z)-m*(z*(1-o.core.delta_iv))-o._derivative(m)*D)/A
            qy = (y*(2*z)-(p-m)*(z*(1-o.core.delta_iv))-o._derivative(p-m)*D)/A
            return {name: {(1, 0): row.truncate(2)} for name, row in dict(U=p, Q=q, Uy=y, Qy=qy).items()}
        eb = errors(bp, bm, iv.mpf(0))
        es = errors(o.hcap*norms['primitive'], o.hcap**2*s*norms['primitive'], norms['value'])
        kb = density_difference(J, iv, Z, o.core.delta_iv, base, eb, o.hcap)
        # Keep Z2 until density_difference takes the Z derivative for C1
        # products; the generic polynomial scale() helper truncates to Z1.
        after_b = {name: add(base[name], {key: row*o.hcap for key, row in eb[name].items()})
                   for name in base}
        ks = density_difference(J, iv, Z, o.core.delta_iv, after_b, es, o.hcap)
        weight = s*iv.exp(2*o.hcap*s)
        kb = {d: self.symmetric(row)*weight for (d, n), row in kb.items()}
        ks = {d: self.symmetric(row)*weight for (d, n), row in ks.items()}
        def coefficient(ctx, value, delta):
            zz = J.variable(ctx, ctx.mpf(Z), 2)
            return (value.truncate(2)*(2*zz)/(1-zz*zz*delta)).truncate(1)
        low = macro['source_Uy']; up = v['source_Uy']
        endpoints = [(coefficient(c, up['point'], o.core.delta)*(-1), coefficient(iv, up['truth'], o.core.delta_iv)*(-1)),
                     (coefficient(c, low['point'], o.core.delta), coefficient(iv, low['truth'], o.core.delta_iv))]
        radius_error = self.symmetric(coefficient(iv, up['truth'], o.core.delta_iv))*s*iv.exp(o.hcap*s)
        def C(ctx, B, delta):
            zz = J.variable(ctx, ctx.mpf(Z), 2)
            return ((B.truncate(2)*(zz*(1-delta))+o._derivative(B.truncate(2))*(1-zz*zz))
                    /(1-zz*zz*delta)).truncate(1)
        vis_point = C(c, low['point'], o.core.delta)*(sum(_bounds(c, norms['mass']))/2)
        vis_truth = C(iv, low['truth'], o.core.delta_iv)*norms['mass']
        vis_remainder = self.symmetric(C(iv, o._norm_error_jet(norms['weighted_primitive_remainder_per_hb'], 3), o.core.delta_iv))
        result = dict(point=point, truth=truth, inlet_errors=kb, source_errors=ks,
                      endpoints=endpoints, radius_error=radius_error, vis_point=vis_point, vis_truth=vis_truth,
                      vis_remainder=vis_remainder, inherited_mean_log_decay=-self.rows.micro.Y,
                      before_width_source_norm=norms['primitive'])
        self._leading[key] = result
        return result

    def finite_errors(self, Z, s, packet):
        _, upbase, geo = self.rows.background(Z, s)
        low_packet = self.rows.macro_rows.evaluate(Z, '1')
        _, lobase, _ = self.rows.macro_rows.background(Z, self.iv.mpf(1))
        up = self.correction_densities(Z, packet, upbase)
        low = self.correction_densities(Z, low_packet, lobase)
        # The switch background boxes cover every prefix0..s, unlike the
        # large-Y macro mean, so the current source box bounds the volume.
        volume = self.correction_densities(Z, packet, upbase)['volume']
        result = {}
        for name, power in (('endpoint', 2), ('viscous', 1)):
            result[name] = [self.rows.sum_logs((None if up[name][k] is None else up[name][k]+power*geo['logR'],
                                               None if low[name][k] is None else low[name][k]+power*self.iv.log(100))) for k in (0, 1)]
        result['volume'] = [None if v is None else v+2*self.iv.log(100)+self.iv.log(s)+2*self.o.hcap*s for v in volume]
        return result, {name: self.rows.log_rows(low[name]) for name in ('endpoint', 'viscous')}

    def evaluate(self, Z, coordinate):
        packet = self.rows.evaluate(Z, coordinate); s, _, recipe = self.rows.micro._coordinate(coordinate, False)
        key = Z, coordinate
        if key in self._packets: return self._packets[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        incoming = self.inlet(Z); terms = copy.deepcopy(incoming['m15_C0_C1_terms'])
        if incoming['actual_selected_integer'] != packet['selected_integer']: raise ValueError('Same selected N required')
        empty = s._mpi_ == iv.mpf(0)._mpi_; zero, zeroi = J.constant(c, 0, 1), J.constant(iv, 0, 1)
        if not empty:
            lead = self.leading(Z, s); finite, lower = self.finite_errors(Z, s, packet)
        else: lead = None; finite = {name: [None, None] for name in ('endpoint', 'volume', 'viscous')}; lower = None
        def width(name):
            atom = o._summary(sum(_bounds(c, o.logh))/2, o.logh)
            atom.update(source_log_atom_only=True, original_native_parameter_name=name,
                        distinct_graph_nodes_with_checked_common_numeric_source=True)
            return atom
        hbB, hbS = width('hbB'), width('hbS')
        def physical(Rpower, degree, mode=0):
            return o._summary(Rpower*c.log(100)-degree*o.core.logLambda-mode*self.rows.micro.Ypoint,
                              Rpower*iv.log(100)-degree*o.core.ll_iv-mode*self.rows.micro.Y)
        binding = dict(interval='first_switch0..s only', unchanged28sector_macro_inlet=True,
                       native_measure='dy=hbS ds', mean_mode='core incoming mode times exp(-Y)*exp(-hbS*s)',
                       hbB_and_hbS_graph_nodes_distinct=True)
        def append(name, p, t, atoms, **extra):
            terms.append(dict(name=name, coefficient_C0_C1=o._jet_summary(p.truncate(1), t.truncate(1)),
                              log_scale_terms=copy.deepcopy(atoms), actual_function_binding=dict(binding, **extra)))
        for d in range(3):
            for n in range(d+1):
                p, t = (zero, zeroi) if empty else (lead['point'][(d, n)]*c.mpf(coordinate), lead['truth'][(d, n)]*s)
                append(f'first_switch_half_S1_epsilon_{d}_mean_{n}', p, t, [physical(2, d, n), hbS],
                       epsilon_degree=d, inherited_mean_degree=n, inherited_exp_minus_nY_atom_retained=True)
                if n < 2:
                    err = zeroi if empty else self.symmetric(lead['truth'][(d, n)])*s*s*iv.exp((2-n)*o.hcap*s)*(iv.mpf(1) if n == 0 else iv.mpf('.5'))
                    append(f'first_switch_short_kernel_error_epsilon_{d}_mean_{n}', zero, err,
                           [physical(2, d, n), hbS, hbS], error_enclosure_only=True)
        for d in (1, 2):
            for kind, atoms in (('inlet', [hbS, hbB]), ('source', [hbS, hbS])):
                error = zeroi if empty else lead[kind+'_errors'].get(d, zeroi)
                append(f'first_switch_{kind}_error_epsilon_{d}', zero, error, [physical(2, d)]+atoms,
                       error_enclosure_only=True, width_retained_by_density_2F_plus_G_plus_Fy=True)
        for k, endpoint in enumerate(('upper', 'lower')):
            p, t = (zero, zeroi) if empty else lead['endpoints'][k]
            append('first_switch_viscous_'+endpoint, p, t, [physical(1, 1), hbS if k == 0 else hbB],
                   actual_Uy_source_endpoint=True, endpoint=endpoint)
        append('first_switch_viscous_upper_radius_error', zero, zeroi if empty else lead['radius_error'],
               [physical(1, 1), hbS, hbS], error_enclosure_only=True)
        append('first_switch_viscous_source_integral', zero if empty else lead['vis_point'], zeroi if empty else lead['vis_truth'],
               [physical(1, 1), hbS, hbS], actual_positive_complement_mass=True, incoming_mean_cancelled_by_exact_ODE_only=True)
        append('first_switch_viscous_source_integral_error', zero, zeroi if empty else lead['vis_remainder'],
               [physical(1, 1), hbS, hbS, hbS], error_enclosure_only=True)
        for name in ('endpoint', 'volume', 'viscous'):
            atoms = [hbS] if name == 'volume' else []; truth = zeroi
            bound = self.rows.sum_logs(finite[name])
            if bound is not None:
                top = iv.mpf(_bounds(self.audit, bound)[1]); atom = o._summary(_bounds(c, top)[1], top)
                atom.update(positive_error_bound_atom_only=True, source_value_selected=False)
                atoms += [atom]; truth = J(iv, [iv.mpf([-1, 1]), iv.mpf([-1, 1])])
            append('first_switch_finite_N_'+name+'_error', zero, truth, atoms,
                   error_enclosure_only=True, both_nonflat_endpoints=True, ordinary_C0_C1_not_cap_derivatives=True)
        out = dict(Z=Z, coordinate=recipe, source_family=self.family, m15_C0_C1_terms=terms,
            incoming_macro_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']), inlet_source_report=self.inlet_report,
            original_sc_tuple=self.t.sc_tuple, actual_selected_integer=copy.deepcopy(packet['selected_integer']),
            exact_empty_increment=empty, first_switch_geometry=copy.deepcopy(packet['geometry']),
            full_prefix_cover=copy.deepcopy(packet['source_cover']),
            finite_N_error_logs_before_volume_width={name: self.rows.log_rows(v) for name, v in finite.items()},
            actual_lower_finite_N_endpoint_budget=lower,
            inherited_core_mean_log_decay=o._summary(-self.rows.micro.Ypoint, -self.rows.micro.Y),
            contract=dict(unchanged28term_inlet_once=True, only_first_switch_moment_added=True,
                inherited_exp_minus_Y_mean_mode_not_rounded_away=True, native_hbS_not_power_or_cap=True,
                short_kernel_width_factored_before_evaluation=True, separate_epsilon_and_mean_degrees=True,
                native_y_derivatives_and_one_phase_measure=True, inlet_error_keeps_hbS_times_hbB=True,
                both_actual_Uy_source_endpoints=True, exact_mean_ODE_used_for_viscous_increment=True,
                signed_control_mass_source_integral_with_error=True, both_finite_N_pressure_endpoints=True,
                complete_prefix_and_quadratic_errors=True, original_selected_N_unchanged=True),
            limits=dict(first_switch_cumulative_signed_m15_with_complete_factored_errors=True,
                resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
                second_switch_post_power_m15=False, independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
                stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out
