"""Cumulative original macro m15 with exponential mean and factored errors.

The Y-geometry proxy is integrated analytically, keeping epsilon powers
separate. Actual L=Y-2hb geometry and changing source errors retain a
formal hb factor. Selected-N errors retain both nonflat endpoints.
"""
import copy
from fractions import Fraction
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_macro_pressure_rows import OriginalMacroPressureRows
from .original_second_micro_m15 import OriginalSecondMicroPressureRows, OriginalSecondMicroM15
from .original_first_micro_m15 import OriginalFirstMicroM15


def add(*polys):
    out = {}
    for poly in polys:
        for key, value in poly.items():
            out[key] = out[key]+value if key in out else value
    return out


def scale(poly, coefficient):
    return {key: (value*coefficient).truncate(1) for key, value in poly.items()}


def multiply(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = i+k, j+l
            value = (x*y).truncate(1)
            out[key] = out[key]+value if key in out else value
    return out


def dz(poly, J):
    return {key: J(value.ctx, [(k+1)*value[k+1] for k in range(value.order)])
            for key, value in poly.items()}


def half_S1_coefficients(J, ctx, Z, delta):
    z = J.variable(ctx, ctx.mpf(Z), 2)
    one = J.constant(ctx, 1, 2)
    A, D = 1-z*z*delta, 1-z*z
    return (dict(q=one/(2*A), q2=one/4, U2=-z*z/(A*A)),
            dict(q=-one/(2*A), qZ=z*(1-delta)/(4*A), q2=-one/4,
                 UqZ=D/(2*A), U2=z*z*(1-delta)/(A*A), UUZ=z*D/(A*A)))


def densities(U, Q, J):
    return dict(q=Q, qZ=dz(Q, J), q2=multiply(Q, Q), U2=multiply(U, U),
                UqZ=multiply(U, dz(Q, J)), UUZ=multiply(U, dz(U, J)))


def combine(coefficients, values):
    return add(*(scale(values[name], coefficient) for name, coefficient in coefficients.items()))


def proxy_modes(J, ctx, Z, delta, axis_j, psi, mean):
    z = J.variable(ctx, ctx.mpf(Z), 3)
    A, D = 1-z*z*delta, 1-z*z
    u0 = 4*z+axis_j
    p, a = psi.truncate(3), (mean-psi).truncate(3)
    derivative = lambda q: J(ctx, [(k+1)*q[k+1] for k in range(q.order)])
    q0 = (u0*(z*(1+delta))-4*D)/A
    qc = (p*(z*(1+delta))-derivative(p)*D)/A
    qa = -(a*(z*(1-delta))+derivative(a)*D)/A
    return (dict(U={(0, 0): u0.truncate(2), (1, 0): p.truncate(2)},
                 Q={(0, 0): q0.truncate(2), (1, 0): qc.truncate(2), (1, 1): qa.truncate(2)}),
            qa.truncate(2))


def kernels(ctx, Y, q):
    """Jn/Rupper^2 for mean modes0,1,2; mode2 is the exact resonance."""
    x = Y*q
    E = ctx.exp(-x)
    return E, [-ctx.expm1(-2*x)/2, E*(-ctx.expm1(-x)), x*ctx.exp(-2*x)]


def proxy_integral(J, ctx, Z, delta, modes, Y, q):
    bc, vc = half_S1_coefficients(J, ctx, Z, delta)
    rho = densities(modes['U'], modes['Q'], J)
    F, G = combine(bc, rho), combine(vc, rho)
    E, W = kernels(ctx, Y, q)
    zero = J.constant(ctx, 0, 1)
    result = []
    for degree in range(3):
        # F(E)-E^2 F(1) cancels the boundary E^2 term identically.
        # This avoids subtracting nearly equal endpoints near q=0.
        result.append((F.get((degree, 0), zero)*2+G.get((degree, 0), zero))*W[0]
                      +(F.get((degree, 1), zero)+G.get((degree, 1), zero))*W[1]
                      +G.get((degree, 2), zero)*W[2])
    return result, E, W


def difference_products(a, b, ea, eb, hcap):
    """Product difference / hb; residual hb in the quadratic error is bounded."""
    return add(multiply(a, eb), multiply(b, ea), scale(multiply(ea, eb), hcap))


def difference_densities(U, Q, EU, EQ, J, hcap):
    return dict(q=EQ, qZ=dz(EQ, J),
        q2=difference_products(Q, Q, EQ, EQ, hcap),
        U2=difference_products(U, U, EU, EU, hcap),
        UqZ=difference_products(U, dz(Q, J), EU, dz(EQ, J), hcap),
        UUZ=difference_products(U, dz(U, J), EU, dz(EU, J), hcap))


class OriginalMacroM15(OriginalFirstMicroM15):
    def __init__(self, *, rows=None, **kwargs):
        if rows is not None and kwargs:
            raise ValueError('Supply retained macro rows or constructor options')
        self.rows = rows if rows is not None else OriginalMacroPressureRows(**kwargs)
        if not isinstance(self.rows, OriginalMacroPressureRows):
            raise ValueError('Actual original macro pressure rows required')
        self.t, self.o = self.rows.t, self.rows.o
        self.ctx, self.iv, self.J, self.audit = self.rows.ctx, self.rows.iv, self.rows.J, self.rows.audit
        self.directory, self.family = self.rows.directory, self.rows.family
        self.hashes, self.code_hashes = dict(self.rows.hashes), dict(self.rows.code_hashes)
        reports = {}
        for suffix, gate in (('second_micro_m15_numeric', 'second_micro_leading_signed_m15_with_complete_factored_finite_N_errors'),
                             ('macro_pressure_rows', 'actual_macro_U_Q_Z2_mean_Z3_Uy_Qy_Z1_with_factored_errors')):
            path = self.directory/(PREFIX+'compliant_current_original_'+suffix+'.json')
            raw = path.read_bytes(); data = json.loads(raw)
            if data['source_family'] != self.family or not data['checks']['passed'] or not data['limits'][gate]:
                raise ValueError('Accepted original macro prerequisite required: '+suffix)
            for root, group, target in ((self.directory, data['input_hashes'], self.hashes),
                                        (ROOT, data['implementation_hashes'], self.code_hashes)):
                for name, digest in group.items():
                    if (target.get(name, digest) != digest
                            or hashlib.sha256((root/name).read_bytes()).hexdigest() != digest):
                        raise ValueError('Changed original macro prerequisite: '+name)
                    target[name] = digest
            self.hashes[path.name] = hashlib.sha256(raw).hexdigest()
            reports[suffix] = data
        self.inlet_report = PREFIX+'compliant_current_original_second_micro_m15_numeric.json'
        self.inlets = {Fraction(p['Z']): p for p in reports['second_micro_m15_numeric']['samples']
                       if p['coordinate'] == dict(kind='exact_decimal_phase', value='2')}
        self._lower_rows = None; self._second_engine = None
        self._uniform = {}; self._packets = {}; self._leading = {}
        name = 'src/openai_ns_reconstruction/original_macro_m15.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def lower_rows(self):
        if self._lower_rows is None:
            self._lower_rows = OriginalSecondMicroPressureRows(precision=self.o.precision, directory=self.directory)
            if self._lower_rows.family != self.family or self._lower_rows.micro.integer != self.rows.micro.integer:
                raise ValueError('Same source family and selected original N at macro inlet required')
            self.hashes.update(self._lower_rows.hashes); self.code_hashes.update(self._lower_rows.code_hashes)
        return self._lower_rows

    def inlet(self, Z):
        key = Fraction(Z)
        if key not in self.inlets:
            if self._second_engine is None:
                self._second_engine = OriginalSecondMicroM15(rows=self.lower_rows())
            self.inlets[key] = self._second_engine.evaluate(Z, '2')
            self.hashes.update(self._second_engine.hashes); self.code_hashes.update(self._second_engine.code_hashes)
        p = self.inlets[key]
        if (p['source_family'] != self.family or Fraction(p['Z']) != key
                or p['original_sc_tuple'] != self.t.sc_tuple or len(p['m15_C0_C1_terms']) != 17):
            raise ValueError('Same-Z complete seventeen-sector second-micro inlet required')
        return p

    def uniform_source_per_hb(self, Z):
        """Absolute original modal prefix norm before its one formal hb."""
        if Z in self._uniform:
            return self._uniform[Z]
        o, b = self.o, self.o.bridge
        p = o._prepared[Z]
        module = __import__(PREFIX+'compliant_actual_bridge_integrals', fromlist=['jet_norm'])
        total = self.iv.mpf(0)
        with mp.workdps(280):
            for name, power, logscale in (('hydro', 1, b.ctx.mpf(0)), ('pressure', 1, b.pressure_log),
                                           ('swirl', 2, b.swirl_log)):
                n = (module.jet_norm(p['micros'][0]['quotient']*p['micros'][0]['direction']['drive_'+name], b.weight)
                     *module.upper(b.ctx, p['R0']**power)*(1+b.cap)
                     +module.jet_norm(p['micros'][1]['quotient']*p['micros'][1]['direction']['drive_'+name], b.weight)
                     *module.upper(b.ctx, p['R0']**power)*b.cap
                     +module.jet_norm(p['quotient'], b.weight)*module.weighted_mode_norm(p['modes']['drive_'+name],
                        module.radius_kernels(b.ctx, p['R0'], b.ctx.mpf(100), p['L'], power), b.weight))*b.ctx.exp(6*b.cap)
                total += o._abs(o._scaled(n, b.source_term(p['data']['phi0']*0, 1, logscale), divide_eps=True))
        self._uniform[Z] = total
        return total

    def at_mean(self, poly, E):
        result = {}
        for (degree, mode), value in poly.items():
            term = value*(E**mode)
            result[degree] = result[degree]+term if degree in result else term
        return result

    def symmetric(self, jet):
        return self.J(self.iv, [self.iv.mpf([-max(abs(v) for v in _bounds(self.audit, jet[k])),
                                          max(abs(v) for v in _bounds(self.audit, jet[k]))]) for k in (0, 1)])

    def leading(self, Z, q):
        key = Z, q._mpi_
        if key in self._leading:
            return self._leading[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        inputs = self.rows.inputs(Z, q); inlet = inputs['incoming']['inlet']
        def modes(ctx, values, de, j):
            return proxy_modes(J, ctx, Z, de, j, values['Psi'], values['Mpsi'])
        pm, pa = modes(c, inlet['point']['values'], o.core.delta, o.core.j)
        tm, ta = modes(iv, inlet['truth']['values'], o.core.delta_iv, o.core.j_iv)
        qp = sum(_bounds(c, q))/2
        pp, _, _ = proxy_integral(J, c, Z, o.core.delta, pm, self.rows.micro.Ypoint, qp)
        tt, E, W = proxy_integral(J, iv, Z, o.core.delta_iv, tm, self.rows.micro.Y, q)
        # Source differences are bounded before hb. The nonzero mean itself
        # is still a full exp(-Yq) mode; only actual-minus-proxy is small.
        Bpsi = self.uniform_source_per_hb(Z)
        Bmean = (inputs['incoming']['delta']['Mpsi']+Bpsi
                 +o._norm(inlet['truth']['values']['Mpsi']-inlet['truth']['values']['Psi'])
                 *(2*q*iv.exp(2*o.hcap*q)))
        epsi, emean = o._norm_error_jet(Bpsi, 3), o._norm_error_jet(Bmean, 3)
        z = J.variable(iv, iv.mpf(Z), 3); D, A = 1-z*z, 1-z*z*o.core.delta_iv
        eq = (epsi*(2*z)-emean*(z*(1-o.core.delta_iv))-o._derivative(emean)*D)/A
        EU, EQ = {(1, 0): epsi.truncate(2)}, {(1, 0): eq.truncate(2)}
        # Uniform mean-mode range is used only in directed error products.
        uniform = {name: {(degree, 0): value for degree, value in self.at_mean(poly, iv.mpf([0, 1])).items()}
                   for name, poly in tm.items()}
        bc, vc = half_S1_coefficients(J, iv, Z, o.core.delta_iv)
        base_density = densities(uniform['U'], uniform['Q'], J)
        diffs = difference_densities(uniform['U'], uniform['Q'], EU, EQ, J, o.hcap)
        ef, eg = combine(bc, diffs), combine(vc, diffs)
        F = combine(bc, densities(tm['U'], tm['Q'], J))
        Fu, Fl = self.at_mean(F, E), self.at_mean(F, iv.mpf(1))
        Gu = self.at_mean(combine(vc, base_density), iv.mpf(1))
        ef, eg = self.at_mean(ef, iv.mpf(1)), self.at_mean(eg, iv.mpf(1))
        # Relative radius/measure errors divided by the original hb.
        gr_up = 4*(1-q)*iv.exp(4*o.hcap*(1-q))
        gr_lo = 4*iv.exp(4*o.hcap)
        gm = gr_lo+(2/self.rows.micro.Y)*(1+o.hcap*gr_lo)
        zero = J.constant(iv, 0, 1)
        bg_errors = []
        for degree in range(3):
            f = ef.get(degree, zero); g = eg.get(degree, zero)
            err = (self.symmetric(Fu.get(degree, zero))*gr_up+self.symmetric(f)*(1+o.hcap*gr_up)
                   +(self.symmetric(Fl.get(degree, zero))*gr_lo+self.symmetric(f)*(1+o.hcap*gr_lo))*E**2
                   +(self.symmetric(Gu.get(degree, zero))*gm+self.symmetric(g)*(1+o.hcap*gm))*W[0])
            bg_errors.append(err)
        # R*Qy for the frozen proxy is constant, so its increment cancels
        # algebraically. Both genuine source endpoints enter the error.
        lower = self.lower_rows().inputs(Z, iv.mpf(2))
        lower_dy = (lower['actual']['Phi']/lower['comparison']['Phi'][0]
                    *lower['comparison']['Psi'][1]*lower['rho']).truncate(2)
        # The independently owned second-chart source has its own interval
        # context. Transfer exact endpoints before combining Taylor jets.
        lower_dy = J(iv, [o._transfer(v) for v in lower_dy.coefficients])
        dy = [inputs['source_Uy']['truth'].truncate(2), lower_dy]
        viscous = []
        for t, decay, source in ((q, E, dy[0]), (iv.mpf(0), iv.mpf(1), dy[1])):
            ey = (source*(2*z)-((epsi-emean)*(z*(1-o.core.delta_iv))
                              +o._derivative(epsi-emean)*D))/A
            gr = 2*(1-t)*iv.exp(2*o.hcap*(1-t))
            viscous.append(self.symmetric(ey.truncate(1))*(1+o.hcap*gr)
                           +self.symmetric((-ta*decay).truncate(1))*gr)
        actual_uniform = dict(U=uniform['U'][(0, 0)]+(uniform['U'][(1, 0)]+epsi.truncate(2)*o.hcap)*o.core.epsilon_iv,
            Q=uniform['Q'][(0, 0)]+(uniform['Q'][(1, 0)]+eq.truncate(2)*o.hcap)*o.core.epsilon_iv)
        result = dict(point=pp, truth=tt, background_error_per_hb=bg_errors,
                    viscous_error_per_epsilon_hb=viscous, actual_uniform=actual_uniform,
                    E=E, W=W, Bpsi=Bpsi, Bmean=Bmean,
                    proxy_R_Qy_constant=True)
        self._leading[key] = result
        return result

    def correction_densities(self, Z, packet, base):
        u, q = self.magnitude(base['U']), self.magnitude(base['Q'])
        errors = {name: [self.rows.read_log(v) for v in packet['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']]
                  for name in ('U', 'Q', 'Q_y')}
        eu, eq = errors['U'], errors['Q']
        values = dict(q=eq[:2], qZ=eq[1:], q2=self.delta_product(q, q, eq, eq),
                      U2=self.delta_product(u, u, eu, eu), UqZ=self.delta_product(u, q[1:], eu, eq[1:]),
                      UUZ=self.delta_product(u, u[1:], eu, eu[1:]))
        bc, vc = half_S1_coefficients(self.J, self.iv, Z, self.o.core.delta_iv)
        return dict(endpoint=self.rows.combine(bc, values, 1), volume=self.rows.combine(vc, values, 1),
                    viscous=errors['Q_y'])

    def finite_errors(self, Z, q, packet, uniform):
        iv = self.iv
        _, upper_base, geo = self.rows.background(Z, q)
        lower_rows = self.lower_rows(); lower_packet = lower_rows.evaluate(Z, '2')
        _, lower_base, _ = lower_rows.background(Z, iv.mpf(2))
        upper = self.correction_densities(Z, packet, upper_base)
        lower = self.correction_densities(Z, lower_packet, lower_base)
        volume = self.correction_densities(Z, packet, uniform)['volume']
        low_radius = self.rows.geometry(iv.mpf(0))['logR']
        def weighted(v, logscale):
            return None if v is None else v+logscale
        result = {}
        for name, power in (('endpoint', 2), ('viscous', 1)):
            result[name] = [self.rows.sum_logs((weighted(upper[name][k], geo['logR']*power),
                                               weighted(lower[name][k], low_radius*power))) for k in (0, 1)]
        # Exact positive physical mass / Rupper^2 is <=min(Lq,1/2),
        # and L<=Y. There is no extra hb on this macro measure.
        mass = iv.mpf(min(_bounds(self.audit, self.rows.micro.Y*q)[1], self.audit.mpf('.5')))
        result['volume'] = [weighted(v, 2*geo['logR']+iv.log(mass)) for v in volume]
        return result, dict(lower_endpoint=self.rows.log_rows(lower['endpoint']),
                            lower_viscous=self.rows.log_rows(lower['viscous']),
                            lower_log_radius=self.o._summary(self.rows.geometry(iv.mpf(0))['logR_point'], low_radius))

    def evaluate(self, Z, coordinate):
        packet = self.rows.evaluate(Z, coordinate)
        q, _, recipe = self.rows.micro._coordinate(coordinate, False)
        key = Z, coordinate
        if key in self._packets:
            return self._packets[key]
        c, iv, J, o = self.ctx, self.iv, self.J, self.o
        incoming = self.inlet(Z); terms = copy.deepcopy(incoming['m15_C0_C1_terms'])
        if incoming['actual_selected_integer'] != packet['selected_integer']:
            raise ValueError('Actual selected integer must be unchanged at macro inlet')
        empty = q._mpi_ == iv.mpf(0)._mpi_
        zero, zeroi = J.constant(c, 0, 1), J.constant(iv, 0, 1)
        logR = iv.log(100)-self.rows.micro.Y*(1-q)
        qp = sum(_bounds(c, q))/2
        logRp = c.log(100)-self.rows.micro.Ypoint*(1-qp)
        h = o._summary(sum(_bounds(c, o.logh))/2, o.logh); h['source_log_atom_only'] = True
        def physical(power, eps_degree):
            return o._summary(power*logRp-eps_degree*o.core.logLambda, power*logR-eps_degree*o.core.ll_iv)
        def append(name, point, truth, atoms, binding):
            terms.append(dict(name=name, coefficient_C0_C1=o._jet_summary(point.truncate(1), truth.truncate(1)),
                              log_scale_terms=copy.deepcopy(atoms), actual_function_binding=binding))
        if empty:
            leading = dict(point=[zero]*3, truth=[zeroi]*3, background_error_per_hb=[zeroi]*3,
                           viscous_error_per_epsilon_hb=[zeroi]*2, proxy_R_Qy_constant=True)
            errors = {name: [None, None] for name in ('endpoint', 'volume', 'viscous')}
            lower_budget = None
        else:
            leading = self.leading(Z, q)
            errors, lower_budget = self.finite_errors(Z, q, packet, leading['actual_uniform'])
        binding = dict(interval='original macro0..q only', incoming_second_m15_retained_once=True,
                       actual_geometry='R=4epsilon*exp(2hb+(Y-2hb)*q)',
                       defining_measure='R*dR=R^2*(Y-2hb)dq', proxy_geometry='Rproxy=4epsilon*exp(Yq)',
                       proxy_discrepancy_charged_with_original_hb=True)
        for degree in range(3):
            append('macro_half_S1_epsilon_'+str(degree), leading['point'][degree], leading['truth'][degree],
                   [physical(2, degree)], dict(binding, epsilon_degree=degree,
                    integrated='(2F0+G0)W0+(F1+G1)W1+G2W2; F2 boundary cancels exactly',
                    nonlinear_and_resonant_mean_modes_retained=True))
            append('macro_source_geometry_error_epsilon_'+str(degree), zero, leading['background_error_per_hb'][degree],
                   [physical(2, degree), h], dict(binding, error_enclosure_only=True,
                    source_and_radius_and_length_errors_before_hb=True, quadratic_source_cross_errors=True))
        for i, name in enumerate(('upper', 'lower')):
            scale_atom = physical(1, 1) if i == 0 else o._summary(c.log(4)-2*o.core.logLambda, iv.log(4)-2*o.core.ll_iv)
            append('macro_leading_viscous_'+name+'_error', zero, leading['viscous_error_per_epsilon_hb'][i],
                   [scale_atom, h], dict(binding, error_enclosure_only=True, endpoint=name,
                    proxy_R_Qy_increment_zero_by_exact_exponential_identity=True,
                    actual_Uy_and_nonzero_mean_endpoint_retained=True))
        for name in ('endpoint', 'volume', 'viscous'):
            bound = self.rows.sum_logs(errors[name]); truth = zeroi; atoms = []
            if bound is not None:
                top = iv.mpf(_bounds(self.audit, bound)[1])
                atom = o._summary(_bounds(c, top)[1], top)
                atom.update(positive_error_bound_atom_only=True, source_value_selected=False)
                atoms = [atom]; truth = J(iv, [iv.mpf([-1, 1]), iv.mpf([-1, 1])])
            append('macro_finite_N_'+name+'_error', zero, truth, atoms,
                   dict(binding, error_enclosure_only=True, both_nonflat_pressure_endpoints_retained=True,
                        complete_first_second_macro_history=True, quadratic_cross_errors_retained=True,
                        ordinary_C0_C1_bounds_not_derivatives_of_caps=True))
        out = dict(Z=Z, coordinate=recipe, source_family=self.family, m15_C0_C1_terms=terms,
            incoming_second_micro_m15_terms=copy.deepcopy(incoming['m15_C0_C1_terms']), inlet_source_report=self.inlet_report,
            original_sc_tuple=self.t.sc_tuple, actual_selected_integer=copy.deepcopy(packet['selected_integer']),
            exact_empty_increment=empty, macro_geometry=copy.deepcopy(packet['geometry']),
            macro_full_prefix_cover=copy.deepcopy(packet['source_cover']),
            integration_method='analytical exponential mean modes with separate epsilon0/1/2; source and geometry errors factored by hb',
            finite_N_error_logs_in_complete_physical_units={k: self.rows.log_rows(v) for k, v in errors.items()},
            incoming_second_endpoint_error_budget=lower_budget,
            contract=dict(complete_seventeen_sector_inlet_retained_once=True, only_macro_increment_added=True,
                nonzero_second_mean_exponential_relaxation=True, epsilon_degrees_kept_separate=True,
                exact_resonant_second_mean_mode_retained=True, near_empty_boundary_subtraction_avoided=True,
                original_hb_factored_source_and_geometry_errors=True, actual_L_not_replaced_by_Y_without_error=True,
                macro_measure_has_no_extra_hb=True, both_actual_pressure_endpoints_retained=True,
                proxy_constant_R_Qy_does_not_zero_actual_endpoint_error=True,
                original_macro_modes_not_core_radial_extrapolation=True,
                mean_ZZZ_and_Q_ZZ_log_errors_retained=True, actual_selected_N_not_candidate_257=True,
                all_source_cells_and_axial_frames_retained=True),
            limits=dict(macro_cumulative_leading_signed_m15_with_complete_factored_errors=True,
                resolved_finite_N_signed_centers=False, useful_total_relative_accuracy=False,
                remaining_active_compact_signed_pressure=False, independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False, common_global_N=False, exact_heat_exterior=False,
                stress_cone=False, genuine_temporal_n_ge_2=False, full_corrected_NS=False))
        self._packets[key] = out
        return out
