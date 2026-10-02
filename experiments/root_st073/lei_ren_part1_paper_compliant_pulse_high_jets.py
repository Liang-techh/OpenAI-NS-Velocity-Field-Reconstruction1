"""Selected O.4 partial primitives through C4 and radial recovery through C3.

The existing segmented main/gap/end kernels are independent of Z. This
provider transports the admitted actual C4 coefficients through those same
source formulas. Paper (3.9) consumes one Mz derivative, so radial velocity
is explicitly C3 here; no whole-pulse C4 claim or Taylor remainder is made.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_axial_high_jets import CompliantAxialHighJets
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField, decay_integral
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_candidate_pressure_function import q_jets
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def derivative(jet):
    if jet.order < 1: raise ValueError('A retained derivative coefficient is required')
    return IntervalTaylor(jet.ctx, [(n+1)*jet[n+1] for n in range(jet.order)])


class _FutureSource:
    def __init__(self, high): self.high = high
    def __getattr__(self, name): return getattr(self.high.base.future, name)
    def future(self, Z):
        r = dict(self.high.energy.future(endpoints(self.high.ctx.mpf(Z))))
        r['complete_future_energy_Taylor'] = copy_jet(self.high.ctx, r['complete_future_energy_Taylor'])
        return r


class _SelectedSource:
    def __init__(self, high): self.high = high; self.future = _FutureSource(high)
    def __getattr__(self, name): return getattr(self.high.base, name)
    def select(self, Z): return self.high.select(Z)


class CompliantPulseHighJets(CompliantAxialPulseField):
    def __init__(self):
        self.high = CompliantAxialHighJets(); self.ctx = self.high.ctx
        self.selection = _SelectedSource(self.high); self.pulse = self.high.base.pulse
        self.hashes = dict(self.high.hashes); self.data_cache = {}
        for stem, gate in (('compliant_axial_high_jets_check', 'actual_selected_ap_c1_c2_C4_available'),
                           ('compliant_axial_pulse_field_check', 'actual_O4_partial_pulse_field_installed')):
            name = PREFIX+stem+'.json'; r = json.loads((HERE/name).read_bytes())
            if not r.get(gate) or r['actual_five_defect_family_sha256'] != self.selection.future.angular.initial.family:
                raise ValueError('Pulse high-jet prerequisite/family mismatch: '+name)
            for source, digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                    raise ValueError('Pulse high-jet source changed: '+source)
                self.hashes[source] = digest
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        with mp.workdps(210):
            c = self.ctx; box = lambda v:c.mpf(endpoints(v))
            self.mu = box(self.high.base.mu); self.delta = box(self.high.base.future.delta)
            self.rate = 1-self.mu; self.prate = 1+2*self.mu
            self.logE = box(self.high.base.log_end_scale)
            self.logEcap = box(self.high.base.future.params.log_mu)-1000
            self.logE2cap = 2*box(self.high.base.future.params.log_mu)-1000
            if endpoints(self.logEcap-self.logE)[0] <= 0 or endpoints(self.logE2cap-2*self.logE)[0] <= 0:
                raise ArithmeticError('High pulse formal end-factor cap not proved')
            self.Ecap = c.mpf([0,endpoints(c.exp(self.logEcap))[1]])
            self.E2cap = c.mpf([0,endpoints(c.exp(self.logE2cap))[1]])
            p0 = self.pulse.buffer.power(0, 1)
            self.inlet_H = box(p0['Mtheta_over_sqrt2_R_3half_Pstar'][0])
            self.inlet_P = box(p0['Mp_over_Pstar_squared'][0])
            self.Xp = self.inlet_H/box(p0['Utheta_over_Pstar'][0])
            self.beta_derivative_sup = box(self.high.base.future.repair.beta_derivative_sup)
            self.full_end_energy_weights = [c.exp(-2*self.mu*center)*box(self.pulse.basis['energy_gram'])
                                           for center in (-3,-1)]
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def data(self, Z):
        c = self.ctx; Z = c.mpf(Z); key = tuple(endpoints(Z))
        if key in self.data_cache: return self.data_cache[key]
        with mp.workdps(210):
            selected = self.high.select(Z); k = self.high.constants
            z = IntervalTaylor(c, [Z, 1, 0, 0, 0])
            q = IntervalTaylor(c, [1+Z**2, 2*Z, 1, 0, 0]); r = q.reciprocal()
            u = r*k['U']
            z2 = IntervalTaylor(c, [Z**2, 2*Z, 1, 0, 0])
            inlet = dict(Utheta_over_Pstar=list(u.coefficients),
                Mz_over_R=list((z*k['M']).coefficients),
                Mtheta_z_over_sqrt2_R_3half_Pstar=list((z*r*k['K']).coefficients),
                Mtheta_over_sqrt2_R_3half_Pstar=list((r*self.inlet_H).coefficients),
                Mztheta_over_R_Pstar_squared=list((z2*k['E_Z']+r*r*k['E_Q']).coefficients),
                Mp_over_Pstar_squared=[self.inlet_P*v for v in q_jets(c, Z, 4)])
            incoming = selected['incoming']
            result = selected, inlet, u, incoming['energy_Taylor'], incoming['moment_Taylor']
            self.data_cache[key] = result
            return result

    def pressure_moment(self, Z, t, inlet, u, log_parts=None):
        c = self.ctx
        if log_parts is None: log_parts = dict(time_term=-self.prate*t)
        log_decay = sum(log_parts.values(), c.mpf(0)); decay = self.factor(log_decay)
        kernel = (decay_integral(c, self.prate, t) if endpoints(self.prate*t)[1] < mp.mpf('1e-20')
                  else (1-decay)/self.prate)
        p = IntervalTaylor(c, inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)
        rows = self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z), 4)['normalized_pressure_coefficients']
        p0 = IntervalTaylor(c, rows)
        return dict(Mp_over_Pstar_squared=p, P0_over_Pstar_squared=p0, P_over_Pstar_squared=p+p0,
            P_y_over_Pstar_squared=u*u*(decay/2),
            exact_pressure_source_unchanged=True, retained_axial_order=4,
            pressure_enclosure_is_absolute_not_a_relative_cone_bound=True,
            retained_swirl_pressure_decay_log_parts=log_parts,
            retained_swirl_pressure_decay_log_enclosure=log_decay)

    def radial(self, Z, B, m1):
        c = self.ctx; Z = c.mpf(Z)
        if B.order != m1.order or m1.order < 1: raise ValueError('Consistent positive-order pulse jets required')
        order = m1.order-1
        z = IntervalTaylor(c, [Z, 1]+[0]*max(0, order-1))
        q = IntervalTaylor(c, [1+Z**2, 2*Z, 1]+[0]*max(0, order-2)).truncate(order)
        d = IntervalTaylor(c, [1-Z**2, -2*Z, -1]+[0]*max(0, order-2)).truncate(order)
        L = IntervalTaylor(c, [1-self.delta*Z**2, -2*self.delta*Z, -self.delta]+[0]*max(0, order-2)).truncate(order)
        b = B.truncate(order); m = m1.truncate(order); mz = derivative(m1)
        numerator = z*b*2-z*m*(1-self.delta)-d*(mz-z*m/q*2)
        return numerator/L

    def _high_packet(self, point):
        c = self.ctx; Z = point['Z']; B = point['Uz_over_Utheta']; m = point['Mz_over_R_Utheta']
        radial = point['Ur_over_sqrt_R_over_2_Utheta']
        # The differentiated normalized Mz ODE retains the same source.
        By = point['Uz_y_over_Utheta']+B*(c.mpf('.5')+self.mu)
        my = B-m*(c.mpf('.5')-self.mu)
        radial_y = self.radial(Z, By, my)-radial*self.mu
        _, inlet, u, _, _ = self.data(Z)
        theta_z = derivative(u)/u.truncate(3)
        radial_z = derivative(radial)+radial.truncate(2)*theta_z.truncate(2)
        radial_yz = derivative(radial_y)+radial_y.truncate(2)*theta_z.truncate(2)
        point.update(Mtheta_over_sqrt2_R_3half_Utheta=IntervalTaylor.constant(c, point['Mtheta_over_sqrt2_R_3half_Utheta'], 4),
            Ur_y_over_sqrt_R_over_2_Utheta=radial_y,
            Ur_Z_over_sqrt_R_over_2_Utheta=radial_z,
            Ur_yZ_over_sqrt_R_over_2_Utheta=radial_yz,
            inlet_Utheta_over_Pstar_Taylor=u,
            Utheta_Z_over_Utheta_Taylor=theta_z,
            factored_Uz_over_Pstar_Taylor=u*B,
            factored_Ur_over_sqrt_R_over_2_Pstar_Taylor=u.truncate(3)*radial,
            physical_velocity_common_radial_factor='exp(inverse_mu_term+finite_offset); inlet Taylor factor excludes this exact positive radial factor',
            inlet_five_primitive_Taylor=inlet,
            Ur_Z_available=True, retained_axial_order=4, retained_radial_velocity_axial_order=3,
            retained_radial_velocity_logradial_derivative_axial_order=3,
            primitive_C4_and_radial_C3_installed=True,
            radial_velocity_C4_available=False, full_pulse_C4_installed=False,
            fifth_derivative_Taylor_remainder_available=False,
            full_outer_C4_certified=False)
        return point

    def main(self, *args, **kwargs): return self._high_packet(super().main(*args, **kwargs))
    def _gap(self, *args, **kwargs): return self._high_packet(super()._gap(*args, **kwargs))
    def end(self, *args, **kwargs): return self._high_packet(super().end(*args, **kwargs))

    def report(self):
        with mp.workdps(210):
            r = super().report()
            r.update(whole_Z_main=self.main([-1,1], 1), whole_Z_end_support=self.end([-1,1], -3),
                whole_Z_support_crossing=self.end([-1,1], ['-3.16','-3.14']),
                pulse_primitives_axial_C4_installed=True, pulse_radial_velocity_axial_C3_installed=True,
                radial_velocity_Z_and_yZ_available=True,
                full_pulse_C4_installed=False, full_outer_C4_certified=False,
                radial_velocity_C4_available=False, fifth_derivative_Taylor_remainder_available=False,
                next_dependency='Recover fifth-order Mz/ap/end functions for radial C4; higher radial mixed derivatives and complete pulse/post-pulse interface bounds')
            return r


def run():
    r = CompliantPulseHighJets().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(r)), indent=2)+'\n', encoding='utf-8')
    print('Same selected pulse: all charts/five primitives C4, radial velocity and logradial derivative C3, Ur_Z/Ur_yZ generated; full pulse C4 pending', flush=True)
    return r


if __name__ == '__main__': run()
