"""Evaluated original first-collar axial jets, own-mean Q and pressure.

All derivatives are native fixed-phase Z derivatives. Finite-width errors
come from the actual prescribed-shear source, including its comparison
direction, before Q/pressure recovery from the same actual histories.
"""
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import PREFIX, ROOT, _box, _bounds
from .original_first_collar_numeric import OriginalFirstCollarNumeric


def bessel_jet_tail(ctx, chi, degree, radius=4, radial_order=0):
    """Differentiated factorial tail retaining the true chi derivatives.

    Expand chi**n about chi[0], using positive convolutions of the higher
    Taylor coefficients. Each binomial/factorial tail has a geometric ratio.
    This is an axial Taylor-coefficient bound, not an ordinary derivative.
    """
    N, i, order = degree, radial_order, chi.order
    if N < order or not 0 <= i <= N:
        raise ValueError('Tail requires degree>=axial order and 0<=radial order<=degree')
    if _bounds(mp.mp, chi[0])[0] < 0 or _bounds(mp.mp, chi[0])[1] > 1:
        raise ValueError('The actual real chi0 must be enclosed in [0,1]')
    J = type(chi)
    v = J(ctx, [ctx.mpf(0)]+[abs(chi[k]) for k in range(1, order+1)])
    powers = [J.constant(ctx, 1, order)]
    for _ in range(order):
        powers.append(powers[-1]*v)
    n, q = N+1, ctx.mpf(radius)/2
    first = ctx.mpf(radius)**(n-i)/(ctx.mpf(2)**n*math.factorial(n-i)*math.factorial(n+1))
    answer = []
    for k in range(order+1):
        result = ctx.mpf(0)
        for m in range(k+1):
            ratio = q*(N+2)/((N+2-i)*(N+3)*(N+2-m))
            if _bounds(mp.mp, ratio)[1] >= 1:
                raise ValueError('Differentiated model tail does not contract')
            result += powers[m][k]*math.comb(n, m)*first/(1-ratio)
        answer.append(result)
    return answer


class OriginalFirstCollarJets(OriginalFirstCollarNumeric):
    def __init__(self, *, precision=80, directory=None):
        super().__init__(precision=precision, directory=directory)
        path = self.core.directory/(PREFIX+'compliant_current_original_first_collar_numeric.json')
        report = json.loads(path.read_text(encoding='utf-8'))
        if not report['checks']['passed'] or report['source_family'] != self.family:
            raise ValueError('Accepted numerical actual first-collar values required')
        for name, digest in report['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed first-collar value implementation: '+name)
            self.code_hashes[name] = digest
        self.hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        self.Jet = self.core.taylor
        self.weight = self._transfer(self.bridge.weight)
        self._jet_inlets = {}
        name = 'src/openai_ns_reconstruction/original_first_collar_jets.py'
        self.code_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _norm(self, jet):
        return sum((self._abs(jet[k])*self.weight**k for k in range(jet.order+1)), self.iv.mpf(0))

    def _norm_error_jet(self, norm, order):
        return self.Jet(self.iv, [self._symmetric(norm/self.weight**k) for k in range(order+1)])

    def _derivative(self, jet):
        return self.Jet(jet.ctx, [(k+1)*jet[k+1] for k in range(jet.order)])

    def _jet_summary(self, point, truth):
        assert point.order == truth.order
        return [self._summary(point[k], truth[k]) for k in range(point.order+1)]

    def inlet_jets(self, Z, *, degree=24, order=5):
        c = self.ctx
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1 < c.mpf(Z) < 1
                or type(degree) is not int or not 8 <= degree <= 24
                or type(order) is not int or not 1 <= order <= 6):
            raise ValueError('Use interior decimal Z, degree8..24, inlet Taylor order1..6')
        key = Z, degree, order
        if key in self._jet_inlets:
            return self._jet_inlets[key]
        result = {}
        for directed, label in ((False, 'point'), (True, 'truth')):
            c = self.iv if directed else self.ctx
            fixed, rows, delta, eps = self.core._seed(Z, degree, order, directed=directed)
            for n in range(degree):
                self.core.advance(c, fixed, rows, n, c.mpf(Z), delta, eps)
            J, r = self.Jet, c.mpf(4)
            phi = [J(c, row[:order+1]) for row in rows['A']]
            psi = [J.constant(c, 0, order)]+[J(c, rows['Uz'][n][:order+1])/eps for n in range(1, degree+1)]
            zero = J.constant(c, 0, order)
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
            radial = {name:[sum((polys[n]*(math.factorial(n)//math.factorial(n-i))*r**(n-i)
                for n in range(i, degree+1)), zero) for i in range(3)] for name, polys in (('Phi', phi), ('Psi', psi))}
            z = J.variable(c, c.mpf(Z), order)
            j = self.core.j_iv if directed else self.core.j
            H = z*((1-delta)/2)+(1-z*z)*(4*z+j)
            h2 = J(c, [H[0]**2]+list((H*H).coefficients[1:]))
            chi = h2/(h2+(j/500)**2)
            if directed:
                lo, hi = _bounds(self.ctx, chi[0])
                chi = J(c, [c.mpf([max(0, lo), min(1, hi)])]+list(chi.coefficients[1:]))
                correction = _box(c, self.core.major['scaled_map_size_upper'])/(1-_box(c, self.core.major['scaled_map_Lipschitz_upper']))
                tails = {}
                for i in range(3):
                    model = bessel_jet_tail(c, chi, degree, radial_order=i)
                    small = [correction*self.core.tail_factor(c, degree=degree, radial_order=i,
                        axial_order=k, radius=r, h=_box(c, self.core.tube['Xh_parameter']))['tail_per_Xh_norm']
                        /math.factorial(k) for k in range(order+1)]
                    tails[i] = [model[k]+small[k] for k in range(order+1)], small
                    for name, error in (('Phi', tails[i][0]), ('Psi', small)):
                        radial[name][i] += J(c, [self._symmetric(t) for t in error])
                ep, eu = tails[0]
                pn = [sum((abs(v[k])*r**n for n, v in enumerate(phi)), c.mpf(0)) for k in range(order+1)]
                un = [sum((abs(v[k])*r**n for n, v in enumerate(psi)), c.mpf(0)) for k in range(order+1)]
                def pt(a, b, ea, eb, k):
                    return sum((a[l]*eb[k-l]+ea[l]*b[k-l]+ea[l]*eb[k-l] for l in range(k+1)), c.mpf(0))
                errors = dict(Phi=ep, Psi=eu, H=ep, Mpsi=eu,
                    Kpsi=[pt(pn, un, ep, eu, k) for k in range(order+1)],
                    Apsi=[pt(un, un, eu, eu, k) for k in range(order+1)],
                    B=[pt(pn, pn, ep, ep, k)/2 for k in range(order+1)],
                    C=[pt(pn, pn, ep, ep, k) for k in range(order+1)])
                values = {name:v+J(c, [self._symmetric(t) for t in errors[name]]) for name, v in values.items()}
            pscale = c.exp(2*(self.core.lp_iv if directed else self.core.logP)
                           -(self.core.ll_iv if directed else self.core.logLambda))
            pressure = J(c, fixed['P0_Z_taylor'][:order+1])/pscale
            # Relative Taylor jet F0(Z+zeta)^2/F0(Z)^2. The constant is exact;
            # original gradient coefficients dress the primitive only once.
            logarithmic_derivative = [2*v/eps for v in fixed['ell_Z_taylor'][:order]]
            ratios = [c.mpf(1)]
            for k in range(1, order+1):
                ratios.append(sum((logarithmic_derivative[i]*ratios[k-1-i] for i in range(k)), c.mpf(0))/k)
            result[label] = dict(values=values, radial=radial, P0_normalized=pressure,
                                 F0_squared_ratio=J(c, ratios))
        self._jet_inlets[key] = result
        return result

    def _source_jet_increments(self, prepared, s, order):
        with mp.workdps(280):
            angular, axial = self.bridge.contributions(prepared, 'first', s)
        J, iv = self.Jet, self.iv
        ell, psi = J.constant(iv, 0, order), J.constant(iv, 0, order)
        B_ell, B_psi = iv.mpf(0), iv.mpf(0)
        for term in angular:
            row = J(iv, [self._scaled(v, term) for v in term['normalized_axial_coefficients'][:order+1]])
            ell += row
            B_ell += self._norm(row)
        for terms in axial.values():
            for term in terms:
                row = J(iv, [self._scaled(v, term, divide_eps=True) for v in term['normalized_axial_coefficients'][:order+1]])
                error = self._scaled(term['nonlinear_prefix_error_weighted_norm_per_next_hb'],
                    term, extra_power=1, divide_eps=True, symmetric=True)
                error_jet = self._norm_error_jet(error, order)
                psi += row+error_jet
                # A norm bound e gives coefficient bounds e/weight**k. We
                # use its original norm e once, not the sum of those caps.
                B_psi += self._norm(row)+self._abs(error)
        return ell, B_ell, psi, B_psi

    def _refined_increments(self, Z, s, degree, order, B_ell, B_phi):
        """Bound finite-width change from the exact core-join coefficient.

        The known comparison's width polynomials bound its change from the
        *same actual inlet*. Difference norms propagate through the original
        direction after F0 cancellation. This avoids subtracting independent
        enormous F0-relative jets to recover a small normalized derivative.
        """
        iv, J, owner = self.iv, self.Jet, self
        high = self.inlet_jets(Z, degree=degree, order=order+1)['truth']
        a = high['values']
        m = order+1
        zero = iv.mpf(0)
        class Change:
            def __init__(self, base, error=0):
                self.base, self.error = base, iv.mpf(error)
            def coerce(self, other):
                if isinstance(other, Change):
                    return other
                return Change(other if isinstance(other, J) else J.constant(iv, other, self.base.order))
            def __add__(self, other):
                other = self.coerce(other)
                return Change(self.base+other.base, self.error+other.error)
            __radd__ = __add__
            def __neg__(self):
                return Change(-self.base, self.error)
            def __sub__(self, other):
                return self+-self.coerce(other)
            def __rsub__(self, other):
                return self.coerce(other)+-self
            def __mul__(self, other):
                other = self.coerce(other)
                return Change(self.base*other.base, owner._norm(self.base)*other.error
                    +self.error*owner._norm(other.base)+owner._width_error(self.error*other.error))
            __rmul__ = __mul__
            def inverse(self):
                base = self.base.reciprocal()
                norm = owner._norm(base)
                relative = owner._width_error(norm*self.error)
                if _bounds(owner.ctx, relative)[1] >= owner.ctx.mpf('.5'):
                    raise ArithmeticError('Comparison denominator leaves its inlet inverse neighborhood')
                return Change(base, norm**2*self.error/(1-relative))
            def __truediv__(self, other):
                return self*self.coerce(other).inverse()
            def dz(self):
                return Change(owner._derivative(self.base), self.error*self.base.order/owner.weight)
        prepared = self._prepared[Z]
        cover = prepared['micros'][0]['comparison']
        def norm_old(jet):
            return self._norm(J(iv, [self._transfer(v) for v in jet.coefficients[:m+1]]))
        def variation(poly):
            # Exact width polynomial: row0+h row1+h² row2+h³ remainder.
            return (norm_old(poly.rows[1])+self._width_error(norm_old(poly.rows[2]))
                    +self._width_error(self._width_error(self._transfer(poly.remainder))))
        z = Change(J.variable(iv, iv.mpf(Z), m))
        de, eps = self.core.delta_iv, self.core.epsilon_iv
        U0, d, L = 4*z+self.core.j_iv, 1-z*z, 1-z*z*de
        phi = Change(a['Phi'], variation(cover['phi']))
        V = Change(U0.base+a['Psi']*eps, variation(cover['V']))
        raw = dict(H=a['H'], M=U0.base+a['Mpsi']*eps,
            K=U0.base*a['H']+a['Kpsi']*eps,
            A=U0.base*U0.base+U0.base*a['Mpsi']*(2*eps)+a['Apsi']*(eps*eps),
            B=a['B'], C=a['C'])
        moments = {name:Change(base, variation(cover['moments'][name])) for name, base in raw.items()}
        H, M, K, A, B, C = (moments[n] for n in ('H', 'M', 'K', 'A', 'B', 'C'))
        # Original rho=4exp(hs): its difference divided by h is <=4s exp(hs).
        radius = Change(J.constant(iv, 4, m), 4*iv.mpf(s)*iv.exp(self.hcap*iv.mpf(s)))
        # Fixed original gradient and implicit S=eps²F0² (all true Taylor
        # source derivatives retained in its directed Cauchy enclosure).
        fixed, _, _, _ = self.core._seed(Z, degree, m, directed=True)
        g = Change(J(iv, [-v for v in fixed['ell_Z_taylor'][:m+1]]))
        S = Change(J(iv, fixed['S_Z_taylor'][:m+1]))
        P0 = Change(high['P0_normalized'])
        W = 1-z*M*(1-de)-d*M.dz()
        angular = H*(1-de/2)-z*H.dz()*((1-de)/2)-d*K.dz()+z*K*(2*de-1)
        Dbar = radius*((-W+angular/(2*phi))*eps+g*(z*H*((1-de)/2)+d*K)/(2*phi))/L
        hydro = (-W*V+(M-z*M.dz())*((1-de)/2)+z*A*(2*de)-d*A.dz())/(2*L)
        pressure = (z*P0*(2*(1+de))-d*P0.dz())/(2*L)
        SB, SC = S*B, S*C
        swirl = (-z*SB*(2*de)+d*SB.dz()+z*SC*(2*(1+de))-d*SC.dz())/(2*L)
        actual_phi = Change(a['Phi'].truncate(order), B_phi)
        ratio = actual_phi/Change(a['Phi'].truncate(order), variation(cover['phi']))
        drive = ratio*(radius*hydro+radius*pressure*iv.exp(2*self.core.lp_iv)+radius*radius*swirl/eps)
        phase = iv.mpf(s)
        exponent = iv.exp(self._width_error(B_ell))
        error_phi = self._width_error(self._norm(a['Phi'])*(phase*(self._norm(Dbar.base)
            +(1+self.hcap)*Dbar.error)+B_ell**2*exponent)/2)
        error_psi = self._width_error(phase*(self._norm(drive.base)+(1+self.hcap)*drive.error))
        with mp.workdps(280):
            mass = self._transfer(self.bridge.pulse_masses(s)[0])
        actual_inlet = self.inlet_jets(Z, degree=degree, order=order)['truth']
        return (actual_inlet['radial']['Phi'][1]*(4*mass)+self._norm_error_jet(error_phi, order),
                actual_inlet['radial']['Psi'][1]*(4*mass)+self._norm_error_jet(error_psi, order))

    def _calculate(self, Z, s, *, degree=24, order=5):
        c, iv, J = self.ctx, self.iv, self.Jet
        if (type(s) is not str or not c.isfinite(c.mpf(s)) or not 0 <= c.mpf(s) <= 1
                or type(order) is not int or not 1 <= order <= 5):
            raise ValueError('Use decimal first-chart phase in [0,1] and axial order1..5')
        inlet = self.inlet_jets(Z, degree=degree, order=order)
        p, a = inlet['point']['values'], inlet['truth']['values']
        if Z not in self._prepared:
            with mp.workdps(280):
                self._prepared[Z] = self.bridge.prepare(Z)
        # Reuse the existing scalar source as the independent nominal mass.
        scalar = super().evaluate(Z, s, degree=degree)
        mass = c.mpf(scalar['original_flat_complement_mass'])
        ell, B_ell, dpsi, B_psi = self._source_jet_increments(self._prepared[Z], s, order)
        exponent = iv.exp(self._width_error(B_ell))
        nonlinear = self._width_error(B_ell**2*exponent/2)
        dphi = a['Phi']*(ell+self._norm_error_jet(nonlinear, order))
        B_phi = self._norm(a['Phi'])*B_ell*exponent
        dphi, dpsi = self._refined_increments(Z, s, degree, order, B_ell, B_phi)
        points = dict(Phi=inlet['point']['radial']['Phi'][1]*(4*mass),
                      Psi=inlet['point']['radial']['Psi'][1]*(4*mass))
        truths = dict(Phi=dphi, Psi=dpsi)
        source = lambda v:dict(H=2*v['Phi'], Mpsi=v['Psi'], Kpsi=2*v['Phi']*v['Psi'],
                              Apsi=v['Psi']**2, B=v['Phi']**2, C=v['Phi']**2)
        source0, source_point = source(a), source(p)
        phi_norm, psi_norm = self._norm(a['Phi']), self._norm(a['Psi'])
        variation = dict(H=2*B_phi, Mpsi=B_psi,
            Kpsi=2*(phi_norm*B_psi+psi_norm*B_phi+self.hcap*B_phi*B_psi),
            Apsi=2*psi_norm*B_psi+self.hcap*B_psi**2,
            B=2*phi_norm*B_phi+self.hcap*B_phi**2,
            C=2*phi_norm*B_phi+self.hcap*B_phi**2)
        for name, rate in dict(H=2, Mpsi=1, Kpsi=2, Apsi=1, B=2, C=1).items():
            base = source0[name]-rate*a[name]
            error = self._width_error(self._norm(base)*rate*iv.mpf(s)**2/2+iv.mpf(s)*variation[name])
            truths[name] = base*iv.mpf(s)+self._norm_error_jet(error, order)
            points[name] = (source_point[name]-rate*p[name])*c.mpf(s)
        if not c.mpf(s):
            points = {name:J.constant(c, 0, order) for name in points}
            truths = {name:J.constant(iv, 0, order) for name in truths}
        return inlet, points, truths, scalar

    def evaluate(self, Z, s, *, degree=24, order=5):
        inlet, points, truths, scalar = self._calculate(Z, s, degree=degree, order=order)
        p, a = inlet['point']['values'], inlet['truth']['values']
        recovered = {}
        for directed, label in ((False, 'point'), (True, 'truth')):
            c = self.iv if directed else self.ctx
            values, increments = (a, truths) if directed else (p, points)
            z = self.Jet.variable(c, c.mpf(Z), order)
            de = self.core.delta_iv if directed else self.core.delta
            j = self.core.j_iv if directed else self.core.j
            U0, d, L = 4*z+j, 1-z*z, 1-z*z*de
            def Q(v):
                return (2*z*v['Psi']-z*v['Mpsi']*(1-de)-d*self._derivative(v['Mpsi']))/L
            recovered[label] = dict(Qaxis=((z*U0*(1+de)-4*d)/L).truncate(order-1),
                Qpsi_inlet=Q(values), Qpsi_increment_per_hb=Q(increments),
                pressure_C_inlet=values['C']*inlet[label]['F0_squared_ratio'],
                pressure_C_increment_per_hb=increments['C']*inlet[label]['F0_squared_ratio'])
        pp, tt = recovered['point'], recovered['truth']
        return dict(Z=Z, phase=s, source_family=self.family, axial_Taylor_order=order,
            units='coefficient k = native fixed-phase Z derivative / k!; multiply by k! for ordinary derivatives',
            normalized_fields_and_histories={name:dict(inlet=self._jet_summary(p[name], a[name]),
                increment_per_original_hb=self._jet_summary(points[name], truths[name])) for name in points},
            core_inlet_radial_derivatives={name:[self._jet_summary(inlet['point']['radial'][name][i],
                inlet['truth']['radial'][name][i]) for i in range(3)] for name in ('Phi', 'Psi')},
            radial_recovery={name:self._jet_summary(pp[name], tt[name])
                for name in ('Qaxis', 'Qpsi_inlet', 'Qpsi_increment_per_hb')},
            original_Q='Qaxis + epsilon*(Qpsi_inlet + hb*Qpsi_increment_per_hb)',
            pressure=dict(P0_over_Pstar_squared=self._jet_summary(inlet['point']['P0_normalized'], inlet['truth']['P0_normalized']),
                C_inlet_dressed_by_relative_F0_squared=self._jet_summary(pp['pressure_C_inlet'], tt['pressure_C_inlet']),
                C_increment_per_hb_dressed_by_relative_F0_squared=self._jet_summary(pp['pressure_C_increment_per_hb'], tt['pressure_C_increment_per_hb']),
                exact_recovery='P=Pstar^2*p0 + 4*epsilon*F0(Z_base)^2*exp(hb*s)*(dressed_C_inlet+hb*dressed_C_increment)',
                plus_sign=True, original_amplitude_applied_once=True,
                original_radius_exp_hb_s_not_replaced_by_one=True),
            source_scales=scalar['source_scales'],
            axial_jet_weight_exact_mpf_tuples=[list(t) for t in self.weight._mpi_],
            all_axial_drive_parts_and_actual_nonzero_histories_used=True,
            numerical_original_first_collar_axial_jets_with_absolute_errors=True,
            standalone_core_radial_derivative_tails_included=True,
            limits=dict(exact_original_point_oracle=False, nominal_first_width_approximation=True,
                first_second_interface_derivative_join=False, n1_collar_transport=False,
                physical_Cartesian_collar_velocity=False, global_NS_field=False, temporal_n_ge_2_recursion=False))
