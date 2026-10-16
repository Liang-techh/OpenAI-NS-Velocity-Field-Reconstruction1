"""Repaired original P1 through the native selected pulse, as C1 enclosures.

Main uses the actual incoming pressure and radial-free partial primitives.
Gap/end use the checked repaired terminal identity and finite future sources.
Signed coefficient centers and relative pressure cancellation remain open.
"""
import copy
from fractions import Fraction
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_core_velocity import OriginalCoreVelocity
from .original_native_pulse_numeric import OriginalNativePulseVelocity


class OriginalNativePulsePressureNumeric:
    def __init__(self, velocity, *, pressure_cells=128):
        if type(pressure_cells) is not int or not 32 <= pressure_cells <= 1024:
            raise ValueError('Pressure integration cells32..1024 required')
        self.velocity = velocity
        self.native = n = velocity.numeric
        self.ctx, self.iv, self.J = n.ctx, n.iv, n.J
        self.mu, self.delta, self.family = n.mu, n.delta, n.family
        self.directory, self.cells = n.directory, pressure_cells
        self.hashes, self.code_hashes = dict(n.hashes), dict(n.code_hashes)
        self.a, self.b = self.iv.mpf('.5') - self.mu, self.iv.mpf('.5') + self.mu
        self.ell, self.normal = self.iv.mpf('.15'), n.normal
        self._kernels, self._bumps, self._quadratics, self._packets = {}, {}, {}, {}
        namespace = n.field.main.__func__.__globals__
        self.gp, self.gp_energy = namespace['gp'], namespace['gp_energy']
        self.raw_beta = namespace['raw_beta']
        specifications = (
            ('current_original_pulse_radial_pressure_C3_functions',
             'current_original_n1_partial_pulse_pressure_physical_y_FTC_C3_functions_installed'),
            ('current_original_n1_terminal_integral_identities',
             'actual_terminal_integral_functions_zero_not_only_linear_response'),
            ('current_original_n1_repaired_own_C3_terminal_functions',
             'actual_terminal_integral_functions_C3_zero'))
        for stem, gate in specifications:
            p = self.directory / (PREFIX + 'compliant_' + stem + '_check.json')
            receipt = json.loads(p.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not receipt[gate] or receipt['source_family'] != self.family:
                raise ValueError('Checked actual pressure/terminal source required: ' + stem)
            for name, digest in receipt['input_hashes'].items():
                if name in self.hashes:
                    if self.hashes[name] != digest:
                        raise ValueError('Conflicting original pressure source: ' + name)
                elif hashlib.sha256((self.directory / name).read_bytes()).hexdigest() != digest:
                    raise ValueError('Changed original pressure source: ' + name)
                self.hashes[name] = digest
            self.hashes[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        self.source_bindings = dict(fixed_N1p=specifications[0][0], terminal=specifications[2][0],
            radial_free_map=PREFIX+'compliant_current_original_n1_intermediate_own_C3_chart_functions.py::pressure_functions',
            terminal_value_zero_not_only_derivatives=True)
        name = 'src/openai_ns_reconstruction/original_native_pulse_pressure.py'
        self.code_hashes[name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()

    def rat(self, x):
        x = Fraction(x)
        return self.iv.mpf(x.numerator) / x.denominator

    def positive(self, x):
        lo, hi = _bounds(self.ctx, x)
        return self.iv.mpf([max(self.ctx.mpf(0), lo), max(self.ctx.mpf(0), hi)])

    def shift(self, f):
        return self.J(self.iv, [(k+1)*f[k+1] for k in range(f.order)])

    def geometry(self, Z):
        z = self.J.variable(self.iv, self.iv.mpf(Z), 4)
        d, L = 1-z*z, 1-z*z*self.delta
        thetaZ = -2*z/(1+z*z)
        return z, d, L, thetaZ

    def D(self, f, thetaZ):
        return self.shift(f) + f.truncate(f.order-1)*thetaZ

    def kernel_pair(self, xi, t):
        """Positive exp(-a v) and v exp(-a v) gp(xi-mu*v) primitives."""
        iv, c, a = self.iv, self.ctx, self.a
        key = xi._mpi_, t._mpi_
        if key in self._kernels:
            return self._kernels[key]
        lo, hi = _bounds(c, t)
        length = min(lo, c.mpf(1000))
        out = [iv.mpf(0), iv.mpf(0)]
        if length:
            h = iv.mpf(length)/self.cells
            for i in range(self.cells):
                left, right = h*i, h*(i+1)
                box = iv.mpf([_bounds(c, left)[0], _bounds(c, right)[1]])
                g = self.gp(iv, xi-self.mu*box)['value']
                mass = self.native.decay_integral(iv, a, h)
                if _bounds(c, a*h)[1] < c.mpf('1e-20'):
                    first = iv.mpf([_bounds(c, h*h*iv.exp(-a*h)/2)[0], _bounds(c, h*h/2)[1]])
                else:
                    first = self.positive((-iv.expm1(-a*h)-a*h*iv.exp(-a*h))/(a*a))
                decay = iv.exp(-a*left)
                out[0] += decay*mass*g
                out[1] += decay*(left*mass+first)*g
        if hi > length:
            # A positive omitted tail, including uncertain tiny final lengths.
            edge = iv.mpf(length)
            if hi <= 1000:
                tail0 = (iv.mpf(hi)-edge)*iv.exp(-a*edge)
                tail1 = iv.mpf(hi)*tail0
            else:
                tail0 = iv.exp(-a*edge)/a
                tail1 = iv.exp(-a*edge)*(edge/a+1/(a*a))
            for j, tail in enumerate((tail0, tail1)):
                out[j] += iv.mpf([0, _bounds(c, 11*tail)[1]])
        self._kernels[key] = out
        return out

    def bump(self, s, center):
        """Positive future beta mass/first-time moment and local beta rows."""
        s = Fraction(s)
        key = s, center
        if key in self._bumps:
            return self._bumps[key]
        iv, c = self.iv, self.ctx
        ell = Fraction(3, 20)
        start = max(Fraction(-1), min(Fraction(1), (s-center)/ell))
        W, WT, cells = iv.mpf(0), iv.mpf(0), []
        if start < 1:
            width = (1-start)/self.cells
            for i in range(self.cells):
                l, r = start+i*width, start+(i+1)*width
                box = iv.mpf([_bounds(c, self.rat(l))[0], _bounds(c, self.rat(r))[1]])
                v = iv.mpf(center)+self.ell*box
                beta = self.raw_beta(iv, box)/(self.ell*self.normal)
                ds = self.ell*self.rat(width)
                mass = ds*iv.exp(self.a*self.ell*box)*beta
                W += mass
                WT += mass*self.positive(v-self.rat(s))
                cells.append(dict(v=v, beta=beta, ds=ds, mass=mass))
        local = self.rat((s-center)/ell)
        if s <= center-ell or s >= center+ell:
            beta, by = iv.mpf(0), iv.mpf(0)
        else:
            beta = self.raw_beta(iv, local)/(self.ell*self.normal)
            by = beta*(-2*local)/(self.ell*(1-local*local)**2)
        result = dict(W=W, WT=WT, cells=cells, beta=beta, by=by,
            H=-iv.exp(self.a*(iv.mpf(center)-self.rat(s)))*W)
        self._bumps[key] = result
        return result

    def end_quadratic(self, s):
        """Partial finite end integrals; all errors are local positive bounds."""
        s = Fraction(s)
        if s in self._quadratics:
            return self._quadratics[s]
        iv, c = self.iv, self.ctx
        data = [self.bump(s, j) for j in (-3, -1)]
        energies, cross = [iv.mpf(0), iv.mpf(0)], [[iv.mpf(0) for _ in range(2)] for _ in range(2)]
        for j, item in enumerate(data):
            suffix = iv.mpf(0)
            for cell in reversed(item['cells']):
                v, beta, ds = cell['v'], cell['beta'], cell['ds']
                weight = ds*iv.exp(-2*self.mu*v)
                energies[j] += weight*beta*beta
                whole = suffix+cell['mass']
                future = iv.mpf([_bounds(c, suffix)[0], _bounds(c, whole)[1]])
                ownH = -iv.exp(self.a*(iv.mpf((-3, -1)[j])-v))*future
                cross[j][j] += weight*beta*ownH
                for k in range(j+1, 2):
                    later = -iv.exp(self.a*(iv.mpf((-3, -1)[k])-v))*data[k]['W']
                    cross[j][k] += weight*beta*later
                suffix = whole
        result = energies, cross
        self._quadratics[s] = result
        return result

    def unit_atoms(self, packet, kind, extra=()):
        iv = self.iv
        rp = copy.deepcopy(packet['coordinate']['Rp_global_y_log_terms']) + [self.native._row(iv.log(110))]
        # Convert the independent amplitude context without mixing Taylor jets.
        lp = iv.mpf(self.native.amplitude.lp_iv)
        lm = self.native.quiet.numeric.outer.logmui
        Ls = iv.exp(-iv.mpf(2)/3*(iv.log(4)+lm))
        econst = iv.mpf(2)/3*iv.log(4)-iv.log(6)/2
        if kind == 'linear': base = rp + [self.native._row(lp)]
        elif kind == 'quadratic': base = rp + [self.native._row(2*lp)]
        elif kind == 'viscous': base = [self.native._row(lp)]
        elif kind == 'end_linear':
            base = rp + self.native._rows(self.J(iv, [lp, iv.mpf('5.5')/self.mu, -3*Ls, -iv.mpf(4)/3*lm, econst-13]))
        elif kind == 'end_quadratic':
            base = rp + self.native._rows(self.J(iv, [2*lp, -2/self.mu, -6*Ls, -iv.mpf(8)/3*lm, 2*econst-26]))
        elif kind == 'end_quadratic_main_gap':
            base = rp + self.native._rows(self.J(iv, [2*lp, extra[0], -6*Ls, -iv.mpf(8)/3*lm, 2*econst-26]))
            extra = extra[1:]
        elif kind == 'end_viscous':
            base = self.native._rows(self.J(iv, [lp, -iv.mpf('7.5')/self.mu, -3*Ls, -iv.mpf(4)/3*lm, econst-13]))
        else: raise ValueError('Unknown original pressure unit')
        return base + [self.native._row(v) for v in extra]

    def term(self, packet, name, value, kind, extra=()):
        return dict(name=name, coefficient_C0_C1=self.native._rows(value.truncate(1)),
            log_scale_terms=self.unit_atoms(packet, kind, extra),
            actual_function_binding=dict(source=self.source_bindings, partial_integral_not_full_signed_cap=True))

    def main(self, Z, chart, coordinate, packet):
        iv, J = self.iv, self.J
        if Fraction(coordinate) == 0:
            return copy.deepcopy(packet['incoming_P1_at_Rp']['P1_C0_C1_terms'])
        selected, _, u, _, _ = self.native._data(Z)
        z, d, L, thetaZ = self.geometry(Z)
        ap = selected['selected_ap_Taylor'].truncate(4)
        bb, cc = z*ap*2/L, -(z*ap*(1-self.delta)+d*self.D(ap, thetaZ))/L
        Db, Dc, Dap = self.D(bb, thetaZ), self.D(cc, thetaZ), self.D(ap, thetaZ)
        if chart == 'entrance':
            t = iv.mpf(coordinate); xi = self.mu*t
        else:
            xi = iv.mpf(coordinate); t = xi/self.mu
        k, kv = self.kernel_pair(xi, t)
        q = self.native.quiet.numeric._jet(packet['normalized_fields']['Q_over_Utheta'], True)
        g = self.gp(iv, xi)
        B, By = ap*g['value'], ap*(self.mu*g['derivative'])
        A0 = -(1+z*z)*self.native.constants['C1']
        IqZ = Db*k+Dc*kv  # Dtheta(A0)=0 by the canonical original inlet.
        linear = -(q/L+z*IqZ*(1-self.delta)/(2*L))/2
        upper = -q*q/4+z*z*B*B/(L*L)
        G = self.gp_energy(iv, xi, self.native.field.selection.K)
        base = -(d*ap*(Db+Dc/self.a)/L-ap*ap*z*z*(2*self.delta)/(L*L)
                 +z*d*ap*Dap*2/(L*L))*G/2
        hcap = 89*11*self.native.decay_integral(iv, 2, xi)/(self.a*self.a)
        cap = _bounds(self.ctx, hcap)[1]
        herr = iv.mpf([-cap, cap])
        error = -d*ap*Dc*herr/(2*L)
        viscous = (z*B*(1+self.delta)+z*(By-B*self.b)*2-d*self.D(B, thetaZ))/L
        out = copy.deepcopy(packet['incoming_P1_at_Rp']['P1_C0_C1_terms'])
        append = lambda name, value, kind, extra=(): out.append(self.term(packet, name, value, kind, extra))
        append('pulse_main_linear_upper', u*linear, 'linear', (t/2, -xi))
        append('pulse_main_linear_lower', u*A0/(2*L), 'linear')
        append('pulse_main_quadratic_upper', u*u*upper, 'quadratic', (-2*xi,))
        append('pulse_main_quadratic_lower', u*u*A0*A0/4, 'quadratic')
        append('pulse_main_quadratic_integral', u*u*base, 'quadratic', (-self.native.quiet.numeric.outer.logmui,))
        append('pulse_main_quadratic_mu_error', u*u*error, 'quadratic')
        append('pulse_main_finite_viscous_upper', u*viscous, 'viscous', (-t/2, -xi))
        return out

    def end(self, Z, s, packet):
        iv, J = self.iv, self.J
        s = Fraction(s); S = self.rat(s)
        selected, _, u, _, _ = self.native._data(Z)
        z, d, L, thetaZ = self.geometry(Z)
        C = [v.truncate(4) for v in selected['selected_scaled_end_coefficient_Taylor']]
        bb = [z*v*2/L for v in C]
        cc = [-(z*v*(1-self.delta)+d*self.D(v, thetaZ))/L for v in C]
        Db, Dc = [self.D(v, thetaZ) for v in bb], [self.D(v, thetaZ) for v in cc]
        data = [self.bump(s, j) for j in (-3, -1)]
        zero = J.constant(iv, 0, 3)
        B = sum((v*row['beta'] for v, row in zip(C, data)), zero)
        By = sum((v*row['by'] for v, row in zip(C, data)), zero)
        A = sum((v*row['beta']+w*row['H'] for v, w, row in zip(bb, cc, data)), zero)
        IqZ = zero
        for j, center in enumerate((-3, -1)):
            IqZ += (Db[j]*data[j]['W']-Dc[j]*data[j]['WT'])*iv.exp(self.a*center)
        energy, cross = self.end_quadratic(s)
        UU, UUz, Uqz = zero, zero, zero
        for j in range(2):
            UU += C[j]*C[j]*energy[j]
            UUz += C[j]*self.D(C[j], thetaZ)*energy[j]
            Uqz += C[j]*Db[j]*energy[j]
            for k in range(j, 2): Uqz += C[j]*Dc[k]*cross[j][k]
        linear = (-A*iv.exp(self.a*S)/L+z*IqZ*(1-self.delta)/(2*L))/2
        weight = iv.exp(-2*self.mu*S)
        quadratic = (-A*A*(weight/2)+d*Uqz/L+z*z*B*B*(2*weight)/(L*L)
                     -z*z*UU*(2*self.delta)/(L*L)+z*d*UUz*2/(L*L))/2
        viscous = (z*B*(1+self.delta)+z*(By-B*self.b)*2-d*self.D(B, thetaZ))/L
        terms = [self.term(packet, 'pulse_terminal_backward_linear', u*linear, 'end_linear'),
                 self.term(packet, 'pulse_terminal_backward_quadratic', u*u*quadratic, 'end_quadratic'),
                 self.term(packet, 'pulse_terminal_backward_viscous', u*viscous, 'end_viscous', (-S/2, -self.mu*S))]
        return terms, A, self.D(A, thetaZ)

    def gap(self, Z, chart, coordinate, packet):
        iv = self.iv
        out, A, DA = self.end(Z, '-4', packet)
        if chart == 'gap_end' and Fraction(coordinate) == -4:
            return out  # Empty gap increment, before boxing opposite boundary terms.
        selected, _, u, _, _ = self.native._data(Z)
        z, d, L, thetaZ = self.geometry(Z)
        K, DK = A*iv.exp(-4*self.a), DA*iv.exp(-4*self.a)
        f = u*z*DK*(1-self.delta)/(4*L)
        if chart == 'gap_end':
            S = iv.mpf(coordinate)
            out.append(self.term(packet, 'pulse_gap_linear_history', f*(-S-4), 'end_linear'))
            local_extra = (-S,)
        else:
            D = 13-iv.mpf(coordinate)
            out.append(self.term(packet, 'pulse_gap_linear_history_inverse_mu', f*D, 'end_linear', (-self.native.quiet.numeric.outer.logmui,)))
            out.append(self.term(packet, 'pulse_gap_linear_history_finite', f*(-4), 'end_linear'))
            # Combine the inverse-mu atoms before boxing: 2logE-s gives (11-xi)/mu.
            local_extra = ((11-iv.mpf(coordinate))/self.mu,)
        out.append(self.term(packet, 'pulse_gap_quadratic_at_minus4', u*u*K*K*iv.exp(4)/4, 'end_quadratic'))
        local = self.term(packet, 'pulse_gap_quadratic_at_query', -u*u*K*K/4,
                          'end_quadratic' if chart == 'gap_end' else 'end_quadratic_main_gap', local_extra)
        out.append(local)
        return out

    def evaluate(self, Z, chart, coordinate, packet=None):
        key = Z, chart, coordinate
        if key in self._packets: return self._packets[key]
        packet = packet or self.native.evaluate(Z, chart, coordinate)
        with mp.workdps(self.iv.dps+30):
            if chart == 'end': terms = self.end(Z, coordinate, packet)[0]; mode = 'backward_terminal'
            elif chart in ('gap', 'gap_end') or (chart == 'main' and Fraction(coordinate) == 11):
                terms = self.gap(Z, 'gap' if chart == 'main' else chart, coordinate, packet); mode = 'backward_terminal_and_gap_history'
            else:
                terms = self.main(Z, chart, coordinate, packet); mode = 'forward_actual_Rp'
        result = dict(P1_C0_C1_terms=terms, representation=mode, ordinary_Z_order=1,
            full_incoming_P1_at_Rp_retained=copy.deepcopy(packet['incoming_P1_at_Rp']),
            terminal_identity_consumed=True, terminal_pressure_value_zero_certified=True,
            fixed_N1p_and_original_physical_measures=True, finite_viscous_endpoints_retained=True,
            leading_pressure_not_reused_as_P1=True, original_logE_applied_once=True,
            inherited_canonical_q_Z_exact_zero=True, partial_quadratic_errors_not_full_signed_caps=True,
            signed_coefficients_are_function_enclosures=True, relative_pressure_cancellation_accuracy=False)
        self._packets[key] = result
        return result


class OriginalNativePulsePressureVelocity:
    def __init__(self, *, pressure_cells=128, **kwargs):
        self.velocity = OriginalNativePulseVelocity(**kwargs)
        self.pressure = OriginalNativePulsePressureNumeric(self.velocity, pressure_cells=pressure_cells)
        f = self.velocity
        self.ctx, self.iv, self.J, self.amplitude = f.ctx, f.iv, f.J, f.amplitude
        self.family = f.family
        self.hashes, self.code_hashes = dict(self.pressure.hashes), dict(self.pressure.code_hashes)
        self._summary, self._pair, self._restore = f._summary, f._pair, f._restore

    def evaluate(self, Z, chart, coordinate, *, log_tau='-1', scaled_time=None, azimuth='0'):
        result = self.velocity.evaluate(Z, chart, coordinate, log_tau=log_tau, scaled_time=scaled_time, azimuth=azimuth)
        P1 = self.pressure.evaluate(Z, chart, coordinate, result)
        ell, _, _ = OriginalCoreVelocity._time(self, Z, log_tau, scaled_time)
        time = [self._summary((p*(-2), q*(-2))) for p, q in ell]
        for part in P1['P1_C0_C1_terms']:
            result['physical']['pressure'].append(dict(name='n1_pressure_'+part['name'], temporal_order=1,
                coefficient=copy.deepcopy(part['coefficient_C0_C1'][0]),
                log_scale_terms=copy.deepcopy(time)+copy.deepcopy(part['log_scale_terms'])))
        result['reserved_n1_pressure'] = P1
        result['physical'].update(pressure_orders_included=[0, 1], full_n1_pressure_included=True,
            n1_pressure_time_weight='lambda^(-2)', useful_combined_pressure_relative_error_certified=False)
        result['limits'].update(interior_P1=True, P1_ordinary_Z1=True)
        return result
