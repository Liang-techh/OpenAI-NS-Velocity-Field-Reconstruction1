"""Actual incoming functions and positive selected axial branch through C4.

Fixed incoming row factors multiply every coefficient. Positive end scales
remain formal nonzero. Higher orders use the actual scalar implicit equation,
not derivatives of an interval root iteration or a nominal amplitude.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_axial_amplitude_selection import CompliantAxialAmplitude
from lei_ren_part1_paper_compliant_future_energy_high_jets import CompliantFutureEnergyHighJets, copy_jet
from lei_ren_part1_paper_compliant_outer_angular_repair import intersect
from lei_ren_part1_paper_compliant_outer_initial import turnoff_kernels
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def positive_quadratic_jets(c, a0, A2, A1, A0, prior_C1=None):
    """Implicit C4 of A2*a(Z)^2+A1(Z)*a(Z)+A0(Z)=0, A2 constant."""
    if A1.order != A0.order:
        raise ValueError('Quadratic coefficient orders must agree')
    D = 2*A2*a0+A1[0]
    if endpoints(D)[0] <= 0:
        raise ArithmeticError('Selected positive root derivative denominator lost')
    coefficients = [a0]
    for n in range(1, A0.order+1):
        cross = A2*sum((coefficients[i]*coefficients[n-i] for i in range(1, n)), c.mpf(0))
        mixed = sum((A1[i]*coefficients[n-i] for i in range(1, n+1)), c.mpf(0))
        an = -(cross+mixed+A0[n])/D
        if n == 1 and prior_C1 is not None:
            an = intersect(c, an, prior_C1[1])
        coefficients.append(an)
    return IntervalTaylor(c, coefficients), D


class CompliantAxialHighJets:
    def __init__(self):
        self.energy = CompliantFutureEnergyHighJets()
        self.base = CompliantAxialAmplitude(); self.ctx = self.base.ctx
        self.cache = {}; self.hashes = dict(self.energy.hashes); self.hashes.update(self.base.hashes)
        for stem, gate in (('compliant_future_energy_high_jets_check', 'complete_future_corrected_energy_C4_available'),
                           ('compliant_axial_amplitude_selection_check', 'all_passed')):
            name = PREFIX+stem+'.json'; r = json.loads((HERE/name).read_bytes())
            if not r.get(gate) or r['actual_five_defect_family_sha256'] != self.base.future.angular.initial.family:
                raise ValueError('Axial C4 prerequisite/family mismatch: '+name)
            for source, digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                    raise ValueError('Axial C4 source changed: '+source)
                self.hashes[source] = digest
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self._incoming_constants()
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def _incoming_constants(self):
        c = self.ctx; initial = self.base.pulse.initial; buffer = self.base.pulse.buffer
        # At Z=0 the C1 packet isolates U, M, K and E_Q exactly. The
        # positive E_Z is accumulated directly, never recovered by subtracting
        # it from the much larger negative swirl energy.
        with mp.workdps(210):
            p0 = buffer.power(0, 1, cells=128)
            box = lambda value: c.mpf(endpoints(value))
            U = box(p0['Utheta_over_Pstar'][0])
            M = box(p0['Mz_over_R'][1])
            K = box(p0['Mtheta_z_over_sqrt2_R_3half_Pstar'][1])
            EQ = box(p0['Mztheta_over_R_Pstar_squared'][0])
            kernels = turnoff_kernels(initial.ctx, initial.params.yd, initial.params.Md, cells=128)
            td = box(initial.params.yd)-1
            EZ = 16*box(initial.invP2)*(c.exp(-td)+box(kernels['B_squared_mass']))
            EZ *= c.exp(-1)*c.exp(-box(initial.params.Tw))
            invP = c.exp(-box(initial.params.logPstar))
            if any(endpoints(v)[0] <= 0 for v in (U, M, K, EZ)) or endpoints(EQ)[1] >= 0:
                raise ArithmeticError('Actual incoming shape constant signs lost')
            self.constants = dict(U=U, M=M, K=K, E_Q=EQ, E_Z=EZ,
                C1=invP*M/U, C2=invP*K/(U*U), C0=EQ/(U*U), C_E=EZ/(U*U),
                td=td, Tw=box(initial.params.Tw), B_squared_mass=box(kernels['B_squared_mass']),
                retained_far_tail=box(kernels['retained_far_tail']),
                source_shape='u=U/(1+Z^2), m=M*Z, k=K*Z/(1+Z^2), e=E_Z*Z^2+E_Q/(1+Z^2)^2',
                positive_E_Z_definition='16*Pstar^-2*(exp(-td)+K_B2(yd))*exp(-1-Tw)',
                energy_normalization_has_no_extra_invP=True)

    def incoming(self, Z):
        c = self.ctx; Z = c.mpf(Z)
        if endpoints(Z)[0] < -1 or endpoints(Z)[1] > 1:
            raise ValueError('Z in [-1,1] required')
        g = IntervalTaylor(c, [Z+Z**3, 1+3*Z**2, 3*Z, 1, 0])
        h = IntervalTaylor(c, [Z**2+2*Z**4+Z**6, 2*Z+8*Z**3+6*Z**5,
            1+12*Z**2+15*Z**4, 8*Z+20*Z**3, 2+15*Z**2])
        k = self.constants
        return dict(moment_Taylor=[g*k['C1'], g*k['C2']],
            energy_Taylor=h*k['C_E']+k['C0'],
            exact_shapes='mi=Ci*(Z+Z^3); e=C0+C_E*(Z^2+2Z^4+Z^6)',
            Z_independent_constant_definitions=k, ordinary_Taylor_order=4)

    def select(self, Z):
        c = self.ctx; Z = c.mpf(Z); key = tuple(endpoints(Z))
        if key in self.cache: return self.cache[key]
        with mp.workdps(210):
            b = self.base; incoming = self.incoming(Z)
            rows = [jet*factor for jet, factor in zip(incoming['moment_Taylor'], b.incoming_factor_caps)]
            u = b.linear_inverse(rows[0]*(-b.mu), -(rows[1]-rows[0]))
            future = copy_jet(c, self.energy.future(endpoints(Z))['Section7_34_weighted_future_Taylor'])
            incoming_energy = incoming['energy_Taylor']*b.mu
            target = future-incoming_energy+b.base
            A2 = b.K; A1 = target*0; A0 = -target
            for nu, uj, vj in zip(b.nu, u, b.v):
                A2 += nu*vj*vj; A1 += uj*(2*nu*vj); A0 += uj*uj*nu
            old = b.select(Z); prior = old['selected_ap_Taylor']
            selected, denominator = positive_quadratic_jets(c, prior[0], A2, A1, A0, prior)
            scaled = [uj+selected*vj for uj, vj in zip(u, b.v)]
            if not endpoints(scaled[0][0])[1] < 0 < endpoints(scaled[1][0])[0]:
                raise ArithmeticError('Actual C4 end-coefficient signs lost')
            out = dict(Z=Z, selected_ap_Taylor=selected, selected_scaled_end_coefficient_Taylor=scaled,
                actual_affine_incoming_Taylor=u, actual_row_normalized_incoming_Taylor=rows,
                incoming=incoming, energy_target_Taylor=target,
                actual_weighted_incoming_energy_Taylor=incoming_energy,
                actual_weighted_future_energy_Taylor=future,
                quadratic_coefficients=dict(A2=A2, A1=A1, A0=A0),
                positive_root_derivative_denominator=denominator,
                common_log_end_coefficient_scale=b.log_end_scale,
                fixed_incoming_row_factor_definitions=b.incoming_factor_log_definitions,
                exact_selected_energy_equation=old['exact_selected_energy_equation'],
                actual_selected_ap_c1_c2_C4_available=True, ordinary_Taylor_order=4,
                actual_positive_C0_branch_preserved=True,
                all_incoming_orders_use_same_fixed_row_factors=True,
                fifth_derivative_Taylor_remainder_available=False,
                full_pulse_C4_installed=False, full_outer_C4_certified=False,
                whole_outer_cone_certified=False, temporal_recursion=False)
            self.cache[key] = out
            return out

    def report(self):
        with mp.workdps(210):
            return dict(samples=[self.select(z) for z in ('-1', '0', '.5', '1')],
                whole_Z_C4=self.select([-1, 1]), incoming_constants=self.constants,
                actual_five_defect_family_sha256=self.base.future.angular.initial.family,
                implicit_source_sha256=self.base.future.angular.initial.datum.source_sha,
                actual_selected_ap_c1_c2_C4_available=True,
                full_pulse_C4_installed=False, full_outer_C4_certified=False,
                whole_outer_cone_certified=False, temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantAxialHighJets().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)), indent=2)+'\n', encoding='utf-8')
    print('Same-source incoming functions and actual selected axial amplitude/end coefficients C4 generated; full pulse C4/cone pending', flush=True)
    return result


if __name__ == '__main__': run()
