"""Complete corrected future energy with same-source axial derivatives to4.

All pieces retain Rv*Utheta(Rv,Z)^2 units, q^2 factors, actual angular
corrections, epsilon atoms and the entire exact Gamma exterior. This is
radial swirl energy, not a full physical-domain kinetic energy certificate.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_angular_high_jets import CompliantAngularHighJets, log_taylor
from lei_ren_part1_paper_compliant_future_swirl_energy import CompliantFutureSwirlEnergy
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_compliant_outer_buffer import decay_integral
from lei_ren_part1_paper_compliant_outer_angular_repair import intersect
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def copy_jet(c, jet):
    return IntervalTaylor(c, [c.mpf(endpoints(v)) for v in jet.coefficients])


class CompliantFutureEnergyHighJets:
    def __init__(self):
        self.angular = CompliantAngularHighJets()
        self.base = CompliantFutureSwirlEnergy(); self.ctx = self.base.ctx
        self.order = 4; self.cache = {}; self.hashes = dict(self.base.hashes)
        self.hashes.update(self.angular.hashes)
        for stem, gate in (('compliant_angular_high_jets_check', 'angular_coefficient_C4_available'),
                           ('compliant_future_swirl_energy_check', 'all_passed')):
            name = PREFIX+stem+'.json'; r = json.loads((HERE/name).read_bytes())
            if not r.get(gate) or r['actual_five_defect_family_sha256'] != self.base.angular.initial.family:
                raise ValueError('Future C4 prerequisite/family mismatch: '+name)
            for source, digest in r['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                    raise ValueError('Future C4 source changed: '+source)
                self.hashes[source] = digest
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        if self.angular.repair.angular.initial.datum.source_sha != self.base.angular.initial.datum.source_sha:
            raise ValueError('Future C4 implicit pressure sources disagree')
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def future(self, Z):
        c = self.ctx; Z = c.mpf(Z); key = tuple(endpoints(Z))
        if key in self.cache: return self.cache[key]
        with mp.workdps(210):
            f = self.base; high = self.angular.coefficients(endpoints(Z))
            q = IntervalTaylor(c, [1+Z**2, 2*Z, 1, 0, 0]); logq = log_taylor(q)-c.ln(2)
            flatten = q*0
            for i in range(f.cells):
                left = c.mpf(100)*i/f.cells; right = c.mpf(100)*(i+1)/f.cells
                t = c.mpf([endpoints(left)[0], endpoints(right)[1]])
                sig = stable_sigma(c, t/100)[0]
                flatten += (logq*(2*sig)).exp()*((right-left)*c.exp(-2*f.mu*t))
            Nf = q*q*(c.exp(-200*f.mu)/4)
            Nrel = q*q*(c.exp(-2*f.mu*(100+f.Lrel))/4)
            power = Nf*decay_integral(c, 2*f.mu, f.Lrel)
            angular_change = q*0; w = f.repair.weights
            for jet, center in zip(high['physical_coefficient_Taylor'], (-3, -1)):
                dj = copy_jet(c, jet)
                angular_change += (dj*(2*w['E'])+dj*dj*w['F'])*c.exp(-2*f.mu*center)
            angular_change = angular_change*Nrel
            atoms = f.heat.preheat_atoms['energy']; eps = f.heat.epsilon
            baseline = 1/f.delta; epsatom = -2*eps*atoms['W']; eps2atom = eps**2*atoms['W_squared']
            postrel = Nrel*(f.kernels['steep_in']+f.Ns*f.steep_power
                +f.Nq*f.kernels['steep_out']+f.Nt*f.waiting_energy
                +f.tail_multiplier*(baseline+epsatom+eps2atom))
            heat = copy_jet(c, high['scaled_Gamma_future_defect_Taylor']['energy'])
            heatcap = f.tail_multiplier*(f.delta/2)*f.repair.strong_S_cap
            heatdifference = Nrel*(heat*c.mpf([0, endpoints(heatcap)[1]]))
            total = flatten+power+angular_change+postrel-heatdifference
            old = f.future(Z)
            coeffs = list(total.coefficients)
            for j in (0, 1):
                coeffs[j] = intersect(c, coeffs[j], old['complete_future_energy_Taylor'][j])
            total = IntervalTaylor(c, coeffs)
            if total.order != 4 or endpoints(total[0])[0] <= 0:
                raise ArithmeticError('Positive complete future energy C4 lost')
            out = dict(Z=Z, complete_future_energy_Taylor=total,
                Section7_34_weighted_future_Taylor=total*f.weighted_factor,
                units=old['units'], Section7_34_units=old['Section7_34_units'],
                pieces=dict(flatten=flatten, power_buffer=power, angular_bump_energy_change=angular_change,
                    steep_waiting_preheat_collar_and_infinite_tail=postrel,
                    positive_Gamma_heat_energy_deficit_enclosure=heatdifference),
                separate_preheat_tail_atoms=dict(baseline=baseline, epsilon_atom=epsatom, epsilon_squared_atom=eps2atom),
                Nf=Nf, Nrel=Nrel, ordinary_Taylor_order=4,
                complete_future_corrected_energy_C4_available=True,
                same_C1_source_intersection_used=True,
                exact_heat_energy_deficit_definition=old['exact_heat_energy_deficit_definition'],
                entire_infinite_heat_tail_included=True, signed_angular_changes_retained=True,
                fifth_derivative_Taylor_remainder_available=False,
                full_physical_kinetic_energy_certified=False, full_outer_C4_certified=False,
                whole_outer_cone_certified=False, temporal_recursion=False)
            self.cache[key] = out
            return out

    def report(self):
        with mp.workdps(210):
            return dict(samples=[self.future(z) for z in ('-1', '0', '.5', '1')],
                whole_Z_C4=self.future([-1, 1]), ordinary_Taylor_order=4,
                implicit_source_sha256=self.base.angular.initial.datum.source_sha,
                actual_five_defect_family_sha256=self.base.angular.initial.family,
                Section7_34_weighted_factor=self.base.weighted_factor,
                complete_future_corrected_energy_C4_available=True,
                full_physical_kinetic_energy_certified=False, full_outer_C4_certified=False,
                whole_outer_cone_certified=False, temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantFutureEnergyHighJets().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)), indent=2)+'\n', encoding='utf-8')
    print('Complete same-source corrected future energy C4 generated; actual amplitude high jets/full outer C4/cone pending', flush=True)
    return result


if __name__ == '__main__': run()
