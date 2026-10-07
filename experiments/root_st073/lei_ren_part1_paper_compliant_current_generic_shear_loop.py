"""Executable Section 11 loop, before global source assembly and moment repair.

The scalar engine consumes I/F, not T/F, and keeps arbitrary signed b and
nonzero original U_z. It implements the phase change and zero-mean velocity
primitives, rather than reusing the O3-only a=2+mu*g construction. Its
point computations are arbitrary-precision evaluations, not interval proofs.
The current graph attachment below binds known input gates/support only;
it does not supply the missing whole-chart packet or a new common N.
"""
import hashlib
import json
import operator
from pathlib import Path
import mpmath as mp
import sympy as sy

HERE = Path(__file__).resolve().parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX + 'current_generic_shear_loop.json'
RECEIPT = PREFIX + 'current_generic_shear_loop_check.json'
GATE = 'generic_Section11_shear_loop_and_zero_mean_primitives_implemented'
OPEN = ('whole_upstream_current_source_packet_assembled',
        'current_upstream_generic_shear_loop_installed',
        'current_upstream_modified_cumulative_moments_recovered',
        'current_upstream_new_family_moment_repair_certified',
        'current_upstream_new_common_finite_N_certified',
        'global_completed_stress_cone_certified',
        'actual_coefficient_recursion_certified',
        'full_corrected_NS_velocity_field_certified')


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def encoded(value):
    if isinstance(value, dict):
        return {key: encoded(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encoded(item) for item in value]
    if hasattr(value, '_mpf_'):
        return dict(decimal=mp.nstr(value, 75), exact_mpf_tuple=list(value._mpf_))
    return value


def flat_step(c, x):
    """The fixed paper step, evaluated without cancellation of 1-sigma."""
    if x <= 0:
        return c.mpf(0)
    if x >= 1:
        return c.mpf(1)
    odds = 1 / (1 - x) ** 2 - 1 / x ** 2
    if odds <= 0:
        e = c.exp(odds)
        return e / (1 + e)
    return 1 / (1 + c.exp(-odds))


class GenericLoopScales:
    """Conservative whole-input bounds for (11.6)-(11.7).

    Supplied lower/upper bounds must hold on the complete modification box.
    A collection of scalar fixtures is not evidence for that requirement.
    The engine checks each point against them but cannot certify the box.
    """
    def __init__(self, *, a_min, margin_min, boundary_kappa_excess_min,
                 t0_abs_max, p1_abs_max, p2_abs_max, dps=90):
        if isinstance(dps, bool) or operator.index(dps) < 40:
            raise ValueError('At least 40 decimal digits required')
        self.ctx = c = mp.mp.clone()
        c.dps = operator.index(dps)
        names = ('a_min', 'margin_min', 'boundary_kappa_excess_min',
                 't0_abs_max', 'p1_abs_max', 'p2_abs_max')
        for name, value in zip(names, (a_min, margin_min, boundary_kappa_excess_min,
                                      t0_abs_max, p1_abs_max, p2_abs_max)):
            number = c.mpf(value)
            if not c.isfinite(number):
                raise ValueError('Finite whole-input bound required: ' + name)
            setattr(self, name, number)
        if min(self.a_min, self.margin_min, self.boundary_kappa_excess_min) <= 0:
            raise ValueError('Positive a, stronger relaxed margin and strict edge margins required')
        if min(self.t0_abs_max, self.p1_abs_max, self.p2_abs_max) < 0:
            raise ValueError('Nonnegative whole-input upper bounds required')
        self.d_star = self.margin_min / 4
        self.q_star = c.sqrt(3 / (2 * self.a_min))
        self.B_star = (self.t0_abs_max + 2 * self.q_star
                       + 4 * self.p2_abs_max * self.q_star ** 2 / self.d_star)
        self.J_star = self.p1_abs_max * self.B_star + self.p2_abs_max
        self.eta = min(c.mpf(1) / 2, self.margin_min / 8,
                       self.margin_min ** 2 / (8 * (1 + self.J_star ** 2)),
                       self.boundary_kappa_excess_min / 2)
        if self.eta <= 0 or 2 + self.eta == 2:
            raise ArithmeticError('Increase precision to resolve the input-dependent eta')

    def bounds(self):
        return {name: getattr(self, name) for name in (
            'a_min', 'margin_min', 'boundary_kappa_excess_min',
            't0_abs_max', 'p1_abs_max', 'p2_abs_max',
            'd_star', 'q_star', 'B_star', 'J_star', 'eta')}


class GenericShearLoop:
    def __init__(self, scales, *, a, b, p1, p2, Utheta):
        if type(scales) is not GenericLoopScales:
            raise ValueError('Explicit whole-input Section 11 scales required')
        self.scales = scales
        self.ctx = c = scales.ctx
        for name, value in (('a', a), ('b', b), ('p1', p1), ('p2', p2), ('Utheta', Utheta)):
            number = c.mpf(value)
            if not c.isfinite(number):
                raise ValueError('Finite scalar loop input required: ' + name)
            setattr(self, name, number)
        if self.a <= 0 or self.Utheta <= 0:
            raise ValueError('Positive angular velocity and a required')
        self.t0 = -self.b / self.a
        self.kappa = self.a + self.b ** 2 / self.a
        self.H0 = self.p1 + self.p2 * self.t0
        self.J0 = self.p2 - self.p1 * self.t0
        self.input_D = self.H0 - self.kappa
        self.input_Q = 2 * self.input_D ** 2 - (self.kappa - 2) * self.J0 ** 2
        if self.a < scales.a_min or self.H0 - 2 < scales.margin_min:
            raise ValueError('Scalar input violates whole-input lower bounds')
        if (abs(self.t0) > scales.t0_abs_max or abs(self.p1) > scales.p1_abs_max
                or abs(self.p2) > scales.p2_abs_max):
            raise ValueError('Scalar input violates whole-input upper bounds')
        # (3.23): for kappa<=2, H0>2 is essential; for kappa>2,
        # also require the original signed quadratic cone, including q=0.
        if self.kappa > 2 and (self.input_D <= 0 or self.input_Q <= 0):
            raise ValueError('Original kappa>2 input fails the relaxed cone')
        eta = scales.eta
        if self.kappa >= 2 + eta:
            self.q = c.mpf(0)
        else:
            self.q = flat_step(c, (2 + eta - self.kappa) / eta) * c.sqrt(
                (2 + 2 * eta - self.kappa) / (2 * self.a))
        self.v = self.kappa + 2 * self.a * self.q ** 2
        self.u = self.p2 * self.q / scales.d_star
        self.h = c.sqrt(1 + self.u ** 2)
        self.r = self.u / self.h
        # Keep 1-|r| and 1-r² as positive correlated expressions even
        # when u/h rounds to +/-1. Never subtract two rounded units.
        self.one_minus_abs_r = 1 / (self.h * (self.h + abs(self.u)))
        self.one_minus_r2 = 1 / self.h ** 2
        if self.one_minus_abs_r < c.sqrt(c.eps):
            raise ArithmeticError('Increase precision to resolve the Poisson phase condition number')

    @classmethod
    def from_stress(cls, scales, *, a, b, theta_over_F, axial_over_F, Utheta):
        """Convert T/F to the paper I/F using T=I+S, S/F=(-a,b)."""
        c = scales.ctx
        return cls(scales, a=a, b=b, p1=c.mpf(theta_over_F) + c.mpf(a),
                   p2=c.mpf(axial_over_F) - c.mpf(b), Utheta=Utheta)

    def _denominator(self, psi):
        c = self.ctx
        trig = c.sin(psi / 2) if self.r >= 0 else c.cos(psi / 2)
        return self.one_minus_abs_r ** 2 + 4 * abs(self.r) * trig ** 2

    def direction(self, psi):
        c = self.ctx
        if self.q == 0:
            return self.t0
        den = self._denominator(psi)
        if self.r >= 0:
            numerator = self.one_minus_abs_r - 2 * c.sin(psi / 2) ** 2
        else:
            numerator = 2 * c.cos(psi / 2) ** 2 - self.one_minus_abs_r
        return self.t0 + 2 * self.q / self.h * numerator / den

    def _w_integrals(self, psi):
        c = self.ctx
        if self.r == 0:
            return c.sin(psi), psi / 2 + c.sin(2 * psi) / 4
        r = self.r
        if abs(r) < c.mpf('.2'):
            # w=sum_{n>=1} r^(n-1)cos(n psi); integrate its square.
            # The remainder is geometric and does not involve r^-2.
            W1 = c.mpf(0)
            W2 = psi / (2 * self.one_minus_r2)
            power = c.mpf(1)
            tolerance = c.eps / 32
            for k in range(1, 4 * c.prec + 16):
                sine = c.sin(k * psi) / k
                W1 += power * sine
                coeff = power * r / self.one_minus_r2
                if k >= 2:
                    coeff += (k - 1) * (power / r) / 2
                W2 += coeff * sine
                if k > 3 and abs(power / r) * (k + 1) < tolerance:
                    return W1, W2
                power *= r
            raise ArithmeticError('Small-r antiderivative series did not converge')
        minus = self.one_minus_abs_r if r > 0 else 2 - self.one_minus_abs_r
        plus = 2 - self.one_minus_abs_r if r > 0 else self.one_minus_abs_r
        E = 2 * c.atan2(plus * c.sin(psi / 2), minus * c.cos(psi / 2))
        alpha = 2 * self.h ** 2 - 1
        W1 = (E - psi) / (2 * r)
        W2 = ((alpha - 2) * E + psi + 2 * r * c.sin(psi) / self._denominator(psi)) / (4 * r ** 2)
        return W1, W2

    def integrals(self, psi):
        c = self.ctx
        psi = c.mpf(psi)
        if not c.isfinite(psi) or not 0 <= psi <= 2 * c.pi:
            raise ValueError('Angle must lie in the closed period [0,2pi]')
        if psi == 0:
            return c.mpf(0), c.mpf(0)
        if psi == 2 * c.pi or self.q == 0:
            return self.t0 * psi, (self.t0 ** 2 + 2 * self.q ** 2) * psi
        W1, W2 = self._w_integrals(psi)
        T1 = self.t0 * psi + 2 * self.q / self.h * W1
        T2 = self.t0 ** 2 * psi + 4 * self.t0 * self.q / self.h * W1 + 4 * self.q ** 2 / self.h ** 2 * W2
        return T1, T2

    def phase_at_angle(self, psi):
        c = self.ctx
        psi = c.mpf(psi)
        if psi == 2 * c.pi:
            return c.mpf(1)
        _, T2 = self.integrals(psi)
        return self.a * (psi + T2) / (2 * c.pi * self.v)

    def angle_at_phase(self, phase):
        """Monotone inversion, including exact periodic endpoints."""
        c = self.ctx
        phase = c.mpf(phase)
        if not c.isfinite(phase):
            raise ValueError('Finite phase required')
        original_phase = phase
        phase -= c.floor(phase)
        if phase == 0 and original_phase > 0:
            return 2 * c.pi
        if phase == 0 or self.q == 0:
            return 2 * c.pi * phase
        lo, hi = c.mpf(0), 2 * c.pi
        t_cap = abs(self.t0) + 2 * self.q * (self.h + abs(self.u))
        phase_slope_cap = self.a * (1 + t_cap ** 2) / (2 * c.pi * self.v)
        extra_bits = int(c.ceil(c.log(1 + phase_slope_cap, 2)))
        for _ in range(c.prec + extra_bits + 12):
            mid = (lo + hi) / 2
            if mid == lo or mid == hi:
                break
            if self.phase_at_angle(mid) < phase:
                lo = mid
            else:
                hi = mid
        answer = (lo + hi) / 2
        if abs(self.phase_at_angle(answer)-phase) > 64*c.sqrt(c.eps):
            raise ArithmeticError('Increase precision for the requested phase inversion')
        return answer

    def at_angle(self, psi):
        c = self.ctx
        phi = self.phase_at_angle(psi)
        t = self.direction(psi)
        T1, _ = self.integrals(psi)
        aL = self.v / (1 + t ** 2)
        bL = -aL * t
        A = self.a / 2 * (phi - psi / (2 * c.pi))
        B = self.Utheta / 2 * (-self.a * T1 / (2 * c.pi) - self.b * phi)
        if self.q == 0 or psi == 0 or psi == 2 * c.pi:
            A = B = c.mpf(0)
        H = self.p1 + self.p2 * t
        J = self.p2 - self.p1 * t
        D = H - self.v
        Q = 2 * D ** 2 - (self.v - 2) * J ** 2
        return dict(angle=psi, phase=phi, t=t, aL=aL, bL=bL,
                    A=A, B=B, kappa_L=self.v,
                    frozen_signed_direction=D, frozen_signed_quadratic=Q,
                    phase_angle_derivative=self.a * (1 + t ** 2) / (2 * c.pi * self.v))

    def evaluate(self, phase):
        return self.at_angle(self.angle_at_phase(phase))

    def modulate(self, *, logR_offset, N, Uz, slow_A_y, slow_B_y):
        """(11.17)-(11.18), with phase-held source derivatives supplied.

        No spatial derivative is inferred from a chart selector. The caller
        must provide the actual slow partial derivatives at phase=N*y.
        This returns a local velocity candidate, before cumulative recovery.
        """
        if isinstance(N, bool):
            raise ValueError('Finite positive integer N required')
        try:
            N = operator.index(N)
        except TypeError as exc:
            raise ValueError('Finite positive integer N required') from exc
        if N < 1:
            raise ValueError('Finite positive integer N required')
        c = self.ctx
        y, Uz, Ay, By = (c.mpf(x) for x in (logR_offset, Uz, slow_A_y, slow_B_y))
        if not all(c.isfinite(x) for x in (y, Uz, Ay, By)):
            raise ValueError('Finite phase coordinate, velocity and slow derivatives required')
        point = self.evaluate(N * y)
        x = point['A'] / N
        delta_theta = self.Utheta * c.expm1(x)
        delta_z = point['B'] / N
        return dict(loop=point, N=N, actual_phase=N*y,
                    Utheta=self.Utheta*c.exp(x), Uz=Uz+delta_z,
                    delta_theta=delta_theta, delta_z=delta_z,
                    a_N=point['aL']-2*Ay/N,
                    b_N=c.exp(-x)*(point['bL']+2*By/(N*self.Utheta)),
                    slow_A_y=Ay, slow_B_y=By,
                    whole_current_source_or_common_N_admitted=False)


def moment_increment_densities(c, *, R, Utheta, Uz, delta_theta, delta_z):
    """Exact five increments per dR, including the nonzero original Uz.

    They must be integrated from the unchanged inlet. Local zero increments
    never erase incoming histories. Pressure is P0+Mp with the same P0.
    """
    R, E, V, e, u = (c.mpf(x) for x in (R, Utheta, Uz, delta_theta, delta_z))
    if R <= 0 or not all(c.isfinite(x) for x in (R, E, V, e, u)):
        raise ValueError('Positive finite R and finite velocity increments required')
    return dict(Mz=u, Mtheta=c.sqrt(2*R)*e,
                Mztheta=c.sqrt(2*R)*(V*e+E*u+u*e),
                M2=2*V*u+u*u-E*e-e*e/2,
                Mp=(E*e+e*e/2)/R)


def exact_theorem():
    a, b, p1, p2, t, v, q, psi, phi, E, V, e, u, R, N = sy.symbols(
        'a b p1 p2 t v q psi phi Utheta Uz delta_theta delta_z R N', real=True)
    checks = {}
    def zero(name, lhs, rhs):
        if sy.cancel(sy.expand(lhs-rhs)) != 0:
            raise ArithmeticError('Generic shear-loop identity failed: '+name)
        checks[name] = True
    t0 = -b/a
    zero('T_to_I_full_signed_source_mapping',
         (p1-a)-b/a*(p2+b)+(a+b*b/a-2), p1+p2*t0-2)
    aL = v/(1+t*t)
    bL = -aL*t
    zero('loop_kappa_exact', aL+bL*bL/aL, v)
    D = (p1-aL)-bL/aL*(p2+bL)
    J = (p2+bL)+bL/aL*(p1-aL)
    zero('frozen_signed_direction', D, p1+p2*t-v)
    zero('frozen_signed_transverse', J, p2-p1*t)
    zero('weighted_mean_angular_integrand', aL*a*(1+t*t)/(2*sy.pi*v), a/(2*sy.pi))
    zero('weighted_mean_axial_integrand', bL*a*(1+t*t)/(2*sy.pi*v), -a*t/(2*sy.pi))
    zero('phase_period_from_two_Poisson_moments', a*(1+t0*t0+2*q*q), a+b*b/a+2*a*q*q)
    phi_psi = a*(1+t*t)/(2*sy.pi*v)
    zero('periodic_A_derivative', (a/2*(phi_psi-1/(2*sy.pi)))/phi_psi, -(aL-a)/2)
    zero('periodic_B_derivative', (E/2*(-a*t/(2*sy.pi)-b*phi_psi))/phi_psi, E*(bL-b)/2)
    # Closed period reflection implies zero means with respect to phi.
    T1 = sy.Symbol('T1')
    A = a/2*(phi-psi/(2*sy.pi))
    B = E/2*(-a*T1/(2*sy.pi)-b*phi)
    zero('A_reflection_zero_mean', A.xreplace({phi:1-phi,psi:2*sy.pi-psi}), -A)
    zero('B_reflection_zero_mean', B.xreplace({phi:1-phi,T1:2*sy.pi*t0-T1}), -B)
    for name, lhs, rhs in (
        ('Mz', (V+u)-V, u),
        ('Mtheta', sy.sqrt(2*R)*((E+e)-E), sy.sqrt(2*R)*e),
        ('Mztheta', sy.sqrt(2*R)*((V+u)*(E+e)-V*E), sy.sqrt(2*R)*(V*e+E*u+u*e)),
        ('M2', (V+u)**2-(E+e)**2/2-V*V+E*E/2, 2*V*u+u*u-E*e-e*e/2),
        ('Mp', ((E+e)**2-E*E)/(2*R), (E*e+e*e/2)/R)):
        zero('full_nonzero_source_five_increment_'+name,lhs,rhs)
    Ay, By = sy.symbols('slow_A_y slow_B_y')
    zero('finite_N_theta_total_derivative', a-2*(Ay/N-(aL-a)/2), aL-2*Ay/N)
    zero('finite_N_axial_total_derivative',
         sy.exp(-A/N)*(b+2*(By/N+E*(bL-b)/2)/E),
         sy.exp(-A/N)*(bL+2*By/(N*E)))
    return dict(passed=True, exact_identities=checks,
                source='Lei-Ren supplied v2, Section 11, (11.3)-(11.18); five integrands (2.21)',
                strict_active_phase_margins='H(t)-v>=m/2; 2(H(t)-v)^2-(v-2)J(t)^2>=m^2/4',
                exact_means='integral_0^1 (aL,bL)dphi=(a,b)',
                zero_mean_primitives='A(1-phi)=-A(phi); B(1-phi)=-B(phi)',
                automatic_flat_edges='eta<=m_boundary/2 makes q identically zero near strict K edges',
                independent_spatial_multiplier_not_applied=True,
                finite_scalar_evaluation_is_not_interval_admission=True)


def current_source_attachment():
    """Read only checked receipts; no predecessor constructors or re-runs."""
    records, hashes = {}, {Path(__file__).name: sha(Path(__file__).name)}
    for stem in ('current_inner_exit_strict_collar', 'current_O2_reference_slope_relaxed_cone',
                 'current_O2_axial_relaxed_cone'):
        for suffix in ('.json', '_check.json'):
            name = PREFIX+stem+suffix
            row = json.loads((HERE/name).read_bytes())
            if suffix == '_check.json' and not row['all_passed']:
                raise ValueError('Checked current loop prerequisite required: '+stem)
            for dependency, digest in row['input_hashes'].items():
                if sha(dependency) != digest:
                    raise ValueError('Changed current loop prerequisite: '+dependency)
            hashes.update(row['input_hashes'])
            hashes[name] = sha(name)
            records[stem+suffix] = row
    left = records['current_inner_exit_strict_collar.json']
    proof = left['explicit_current_inner_exit_strict_collar']
    family = {key: left[key] for key in ('actual_five_defect_family_sha256',
              'implicit_source_sha256', 'datum_enclosure_sha256')}
    for stem in ('current_O2_reference_slope_relaxed_cone','current_O2_axial_relaxed_cone'):
        if records[stem+'.json']['source_family'] != family:
            raise ValueError('Loop input current family/source/pressure differs')
    return dict(current_source_family=family,
                cross_graph_attachment_sha256=left['cross_graph_attachment_sha256'],
                left_support_fraction=['1/2','1'],
                selected_first_phase_endpoint=proof['selected_first_phase_endpoint'],
                exact_hb_log_enclosure=proof['exact_source_width_log_enclosure'],
                left_support_strict_kappa_lower=proof['full_source_kappa_lower'],
                prospective_left_coordinate='y=log(R/r_minus); r_minus=Ra*exp(hb*s_c/2)',
                prospective_left_taper_argument='2*first_phase/s_c-1',
                ordinary_y_taper_scale_log='log(2)-log(hb)-log(s_c)',
                original_cutoff_mechanism='q=sigma((2+eta-kappa)/eta)*sqrt((2+2eta-kappa)/(2a))',
                extra_spatial_multiplier_not_applied=True,
                whole_chain_remaining=['first/second/macro bridge','switch','reshape','restore','patch',
                                       'common whole-box p1/p2 and stronger margin','right edge strict collar',
                                       'reserved repair interval in the changed source'],
                scalar_engine_not_a_materialized_current_physical_field=True,
                input_hashes=hashes)


def run():
    attached = current_source_attachment()
    scales = GenericLoopScales(a_min='.8',margin_min='1',boundary_kappa_excess_min='.019',
                               t0_abs_max='1',p1_abs_max='8',p2_abs_max='2')
    examples = {}
    for name, a, b, p1, p2 in (
        ('reference_a_four_fifths_p2_zero','.8','0','4','0'),
        ('signed_general_input','.8','.2','5','.2'),
        ('opposite_signed_general_input','.8','-.2','5','-.2'),
        ('strict_edge_loop_exactly_constant','2.019','0','5','0')):
        loop = GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta='1.3')
        examples[name] = dict(q=loop.q,v=loop.v,scalar_fixture_not_current_source_point=True,
                             phases=[loop.evaluate(p) for p in ('0','.137','.5','.863','1')])
    result = dict(exact_Section11_theorem=exact_theorem(),
                  current_source_input_and_left_support_attachment=attached,
                  computational_example_scales=scales.bounds(),scalar_examples=examples,
                  scope='General signed scalar loop, phase inversion, zero-mean primitives, local finite-N velocity/shear formulas and full five increment densities. Current checked input/support binding only; no global packet, modified histories, repair or new N.',
                  **{GATE:True},**dict.fromkeys(OPEN,False),input_hashes=attached['input_hashes'])
    (HERE/NAME).write_text(json.dumps(encoded(result),indent=2)+'\n',encoding='utf8')
    print('Section 11 generic shear loop and zero-mean velocity primitives generated',flush=True)
    return result


if __name__ == '__main__':
    run()
