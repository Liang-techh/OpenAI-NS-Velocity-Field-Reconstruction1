"""Original first collar: evaluated, width-factored leading fields and histories.

Point coefficients approximate the original core inlet and its first width
variation. The unchanged prescribed-shear source encloses the *finite* collar
increment. Neither a width cap nor a comparison history defines the field.
"""
import hashlib
import importlib
import json
import math

import mpmath as mp

from .original_core_numeric import OriginalCoreNumeric, PREFIX, ROOT, _box, _bounds


class OriginalFirstCollarNumeric:
    def __init__(self, *, precision=70, directory=None):
        self.core = OriginalCoreNumeric(precision=precision, directory=directory)
        self.ctx, self.iv = self.core.ctx, self.core.iv
        self.precision, self.family = precision, self.core.family
        self.hashes = dict(self.core.hashes)
        self.code_hashes = {}
        name = PREFIX+'compliant_actual_bridge_integrals_check.json'
        path = self.core.directory/name
        receipt = json.loads(path.read_text(encoding='utf-8'))
        if not receipt['all_passed'] or any(receipt[k] != v for k, v in self.family.items()):
            raise ValueError('Accepted same-source actual bridge required')
        for filename, digest in receipt['input_hashes'].items():
            if hashlib.sha256((self.core.directory/filename).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed actual bridge input: '+filename)
            self.hashes[filename] = digest
        self.hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        module = importlib.import_module(PREFIX+'compliant_actual_bridge_integrals')
        with mp.workdps(280):
            self.bridge = module.CompliantActualBridgeIntegrals()
        if (self.bridge.core.family != self.family['actual_five_defect_family_sha256']
                or self.bridge.core.source != self.family['implicit_source_sha256']
                or self.bridge.core.datum.datum_sha != self.family['datum_enclosure_sha256']):
            raise ValueError('Numerical inlet and actual collar must share all source data')
        self.hashes.update(self.bridge.hashes)
        self.logh = self._transfer(self.bridge.logh)
        self.hcap = self._transfer(self.bridge.cap)
        self.hbox = self.iv.mpf([0, _bounds(self.ctx, self.hcap)[1]])
        if _bounds(self.ctx, self.logh)[1] >= _bounds(self.ctx, self.iv.ln(self.hcap))[0]:
            raise ValueError('Original hb does not admit the remainder cap')
        for relative in ('src/openai_ns_reconstruction/original_core_numeric.py',
                         'src/openai_ns_reconstruction/original_first_collar_numeric.py'):
            self.code_hashes[relative] = hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._inlets, self._prepared = {}, {}
        self.scale_proofs = []

    def _transfer(self, value):
        return self.iv.mpf([mp.make_mpf(t) for t in value._mpi_])

    def _abs(self, value):
        return self.iv.mpf(max(abs(v) for v in _bounds(self.ctx, value)))

    def _symmetric(self, value):
        bound = _bounds(self.ctx, self._abs(value))[1]
        return self.iv.mpf([-bound, bound])

    def _summary(self, point, truth):
        c, iv = self.ctx, self.iv
        text = c.nstr(point, self.precision)
        error = max(abs(v) for v in _bounds(c, truth-iv.mpf(text)))
        return dict(approximation=text, absolute_error_upper=c.nstr(error*(1+c.mpf('1e-30')), 40),
                    raw_directed_error_exact_mpf_tuple=list(error._mpf_),
                    source_enclosure_exact_mpf_tuples=[list(v) for v in truth._mpi_],
                    error_applies_to_serialized_decimal=True)

    def _width_error(self, coefficient):
        """Bound hb times a large norm using the original log before capping."""
        size = self._abs(coefficient)
        if not size:
            return self.iv.mpf(0)
        log_bound = self.logh+self.iv.ln(size)
        tiny = self.iv.mpf('1e-1000')
        if _bounds(self.ctx, log_bound)[1] >= _bounds(self.ctx, self.iv.ln(tiny))[0]:
            raise ArithmeticError('Original width does not control this history error')
        self.scale_proofs.append(dict(source='hb times normalized Volterra error norm', remaining_hb_power=1,
            log_absolute_upper_exact_mpf_tuples=[list(t) for t in log_bound._mpi_], cap='1e-1000',
            source_value_selected=False))
        return self.iv.mpf(min(_bounds(self.ctx, tiny)[1], _bounds(self.ctx, self.hcap*size)[1]))

    def inlet(self, Z, *, degree=22):
        key = Z, degree
        if key in self._inlets:
            return self._inlets[key]
        packet = self.core.coefficients(Z, degree=degree, depth=1)
        outputs = {}
        for directed, label in ((False, 'point'), (True, 'truth')):
            c = self.iv if directed else self.ctx
            eps = self.core.epsilon_iv if directed else self.core.epsilon
            rows = packet['directed' if directed else 'point']
            phi = [row[0] for row in rows['A']]
            psi = [c.mpf(0)]+[rows['Uz'][n][0]/eps for n in range(1, degree+1)]
            r, zero = c.mpf(4), c.mpf(0)
            def linear(a, weight, norm):
                return sum((v*r**(n+weight+1)/(n+weight+1) for n, v in enumerate(a)), zero)/norm
            def product(a, b, weight, norm):
                return sum((v*w*r**(n+m+weight+1)/(n+m+weight+1)
                            for n, v in enumerate(a) for m, w in enumerate(b)), zero)/norm
            values = dict(Phi=sum((v*r**n for n, v in enumerate(phi)), zero),
                          Psi=sum((v*r**n for n, v in enumerate(psi)), zero),
                          H=linear(phi, 1, 8), Mpsi=linear(psi, 0, 4),
                          Kpsi=product(phi, psi, 1, 8), Apsi=product(psi, psi, 0, 4),
                          B=product(phi, phi, 1, 16), C=product(phi, phi, 0, 4))
            if directed:
                correction = (_box(c, self.core.major['scaled_map_size_upper'])
                              /(1-_box(c, self.core.major['scaled_map_Lipschitz_upper'])))
                tail = correction*self.core.tail_factor(c, degree=degree, radial_order=0,
                    axial_order=0, radius=r, h=_box(c, self.core.tube['Xh_parameter']))['tail_per_Xh_norm']
                ep = tail+c.mpf(2)**(degree+1)/math.factorial(degree+1)/math.factorial(degree+2)/(1-c.mpf(2)/((degree+2)*(degree+3)))
                eu = tail
                pn = sum((abs(v)*r**n for n, v in enumerate(phi)), zero)
                un = sum((abs(v)*r**n for n, v in enumerate(psi)), zero)
                mixed = pn*eu+ep*un+ep*eu
                squared_phi = 2*pn*ep+ep*ep
                errors = dict(Phi=ep, Psi=eu, H=ep, Mpsi=eu, Kpsi=mixed,
                              Apsi=2*un*eu+eu*eu, B=squared_phi/2, C=squared_phi)
                values = {name: v+self._symmetric(errors[name]) for name, v in values.items()}
            else:
                values['Phi_rho'] = sum((n*phi[n]*r**(n-1) for n in range(1, degree+1)), zero)
                values['Psi_rho'] = sum((n*psi[n]*r**(n-1) for n in range(1, degree+1)), zero)
            outputs[label] = values
        self._inlets[key] = outputs
        return outputs

    def _scaled(self, row, term, *, extra_power=0, divide_eps=False, symmetric=False):
        """Remove one formal hb, bounding only the *remaining* source powers."""
        iv = self.iv
        value = self._transfer(row)
        if divide_eps:
            value /= self.core.epsilon_iv
        if not self._abs(value):
            return iv.mpf(0)
        power = term['width_power']-1+extra_power
        extra = self._transfer(term['additional_positive_scale_log_enclosure'])
        source = term['additional_positive_scale_source']
        if source.startswith('F0') or power:
            # F0 and hb are never exponentiated. Check the original normalized
            # product log before using a representable positive error cap.
            bound = iv.ln(self._abs(value))+extra+power*self.logh
            cap = iv.mpf('1e-1000')
            if _bounds(self.ctx, bound)[1] >= _bounds(self.ctx, iv.ln(cap))[0]:
                raise ArithmeticError('Unmaterializable source exceeds its normalized cap')
            self.scale_proofs.append(dict(source=source, remaining_hb_power=power, divide_original_epsilon=divide_eps,
                                          log_absolute_upper_exact_mpf_tuples=[list(t) for t in bound._mpi_],
                                          cap='1e-1000', source_value_selected=False))
            lo, hi = _bounds(self.ctx, value)
            return iv.mpf([-_bounds(self.ctx, cap)[1] if lo < 0 or symmetric else 0,
                           _bounds(self.ctx, cap)[1] if hi > 0 or symmetric else 0])
        value *= iv.exp(extra)
        return self._symmetric(value) if symmetric else value

    def _source_increments(self, prepared, s):
        with mp.workdps(280):
            angular, axial = self.bridge.contributions(prepared, 'first', s)
        iv = self.iv
        ell, B_ell = iv.mpf(0), iv.mpf(0)
        for term in angular:
            value = self._scaled(term['normalized_axial_coefficients'][0], term)
            ell += value
            B_ell += self._abs(value)
        delta_psi, B_psi = iv.mpf(0), iv.mpf(0)
        # Sum absolute contributions before signed cancellation. Because wc/ws
        # are nonnegative and monotone, these bound every prefix 0<=u<=s.
        for terms in axial.values():
            for term in terms:
                value = self._scaled(term['normalized_axial_coefficients'][0], term, divide_eps=True)
                error = self._scaled(term['nonlinear_prefix_error_weighted_norm_per_next_hb'],
                                     term, divide_eps=True, extra_power=1, symmetric=True)
                delta_psi += value+error
                B_psi += self._abs(value)+self._abs(error)
        return ell, B_ell, delta_psi, B_psi

    def evaluate(self, Z, s, *, degree=22):
        c, iv = self.ctx, self.iv
        if type(s) is not str or not c.isfinite(c.mpf(s)) or not 0 <= c.mpf(s) <= 1:
            raise ValueError('Use a decimal first-chart phase in [0,1]')
        inlet = self.inlet(Z, degree=degree)
        p, a = inlet['point'], inlet['truth']
        if Z not in self._prepared:
            with mp.workdps(280):
                self._prepared[Z] = self.bridge.prepare(Z)
        sp, si = c.mpf(s), iv.mpf(s)
        def sigma(u):
            if not u:
                return c.mpf(0)
            if u == 1:
                return c.mpf(1)
            exponent = 1/u**2-1/(1-u)**2
            return c.exp(-exponent)/(1+c.exp(-exponent)) if exponent >= 0 else 1/(1+c.exp(exponent))
        mass = (c.mpf('.5') if sp == 1 else c.quad(lambda u: 1-sigma(u), [0, sp/2, sp])) if sp else c.mpf(0)
        ell, B_ell, dpsi, B_psi = self._source_increments(self._prepared[Z], s)
        exponential = iv.exp(self.hcap*B_ell)
        nonlinear = self.hcap*B_ell**2*exponential/2
        dphi = a['Phi']*(ell+self._symmetric(nonlinear))
        B_phi = self._abs(a['Phi'])*B_ell*exponential
        points = dict(Phi=4*p['Phi_rho']*mass, Psi=4*p['Psi_rho']*mass)
        truths = dict(Phi=dphi, Psi=dpsi)
        source0 = dict(H=2*a['Phi'], Mpsi=a['Psi'], Kpsi=2*a['Phi']*a['Psi'],
                       Apsi=a['Psi']**2, B=a['Phi']**2, C=a['Phi']**2)
        source_point = dict(H=2*p['Phi'], Mpsi=p['Psi'], Kpsi=2*p['Phi']*p['Psi'],
                            Apsi=p['Psi']**2, B=p['Phi']**2, C=p['Phi']**2)
        phi_norm, psi_norm = self._abs(a['Phi']), self._abs(a['Psi'])
        phi_square = 2*phi_norm*B_phi+self.hcap*B_phi**2
        variation = dict(H=2*B_phi, Mpsi=B_psi,
            Kpsi=2*(phi_norm*B_psi+psi_norm*B_phi+self.hcap*B_phi*B_psi),
            Apsi=2*psi_norm*B_psi+self.hcap*B_psi**2, B=phi_square, C=phi_square)
        remainders = {}
        for name, rate in dict(H=2, Mpsi=1, Kpsi=2, Apsi=1, B=2, C=1).items():
            base = source0[name]-rate*a[name]
            error = self._width_error(self._abs(base)*rate*si**2/2+si*variation[name])
            truths[name] = si*base+self._symmetric(error)
            points[name] = sp*(source_point[name]-rate*p[name])
            remainders[name] = error
        if not sp:
            points = {name:c.mpf(0) for name in points}
            truths = {name:iv.mpf(0) for name in truths}
        entries = {name:dict(inlet=self._summary(p[name], a[name]),
                            increment_per_original_hb=self._summary(points[name], truths[name])) for name in truths}
        def terms(name, eps_power=0, multiplier=None):
            multiplier = multiplier or (c.mpf(1), iv.mpf(1))
            return [dict(epsilon_power=eps_power, hb_power=k,
                         coefficient=self._summary(multiplier[0]*pp, multiplier[1]*tt))
                    for k, pp, tt in ((0, p[name], a[name]), (1, points[name], truths[name]))]
        z, zi = c.mpf(Z), iv.mpf(Z)
        U0, U0i = 4*z+self.core.j, 4*zi+self.core.j_iv
        constant = lambda value, truth:dict(epsilon_power=0, hb_power=0, coefficient=self._summary(value, truth))
        raw = dict(Phi=terms('Phi'), Uz=[constant(U0, U0i)]+terms('Psi', 1),
            H=terms('H'), M=[constant(U0, U0i)]+terms('Mpsi', 1),
            K=terms('H', multiplier=(U0, U0i))+terms('Kpsi', 1),
            A=[constant(U0**2, U0i**2)]+terms('Mpsi', 1, (2*U0, 2*U0i))+terms('Apsi', 2),
            B=terms('B'), C=terms('C'))
        return dict(Z=Z, phase=s, chart='original_first_micro_flat_collar', source_family=self.family,
            representation='sum coefficient * epsilon^epsilon_power * hb^hb_power; original positive scales remain formal',
            source_scales=dict(log_epsilon=self._summary(-self.core.logLambda, -self.core.ll_iv),
                               log_hb_enclosure_exact_mpf_tuples=[list(t) for t in self.logh._mpi_],
                               rho='4*exp(hb*s)', R='4*epsilon*exp(hb*s)'),
            original_flat_complement_mass=c.nstr(mass, self.precision),
            normalized_fields_and_own_histories=entries, original_fields_and_six_histories=raw,
            uniform_prefix_increment_bounds=dict(Phi=self._summary(c.mpf(0), self._symmetric(B_phi)),
                                                  Psi=self._summary(c.mpf(0), self._symmetric(B_psi))),
            finite_width_Volterra_error_bounds={n:self._summary(c.mpf(0), self._symmetric(e)) for n,e in remainders.items()},
            numerical_approximation_with_complete_absolute_errors=True,
            nominal_coefficients_are_first_width_approximations=True,
            source_width_or_F0_materialized=False, cap_selected_as_source_value=False,
            comparison_histories_substituted=False, all_three_axial_drive_parts_retained_in_directed_branch=True,
            limits=dict(leading_first_collar_values_only=True, normalized_axial_Z_derivatives=False,
                        physical_collar_velocity=False, n1_collar_transport=False, global_moment_repair=False,
                        temporal_n_ge_2_recursion=False, full_corrected_NS=False))
