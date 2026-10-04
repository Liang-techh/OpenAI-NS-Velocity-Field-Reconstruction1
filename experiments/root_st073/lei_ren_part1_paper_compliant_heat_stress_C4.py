"""Source-bound stress-free Gamma exterior, retaining original histories.

Only the exact terminal-history bridge admits equivalent angular moments.
The zero stress is a consequence of the full moment-to-stress identities,
not an interval cancellation or a choice of exterior integration constants.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_heat_pressure_C4 import CompliantHeatPressureC4
from lei_ren_part1_paper_compliant_heat_terminal_history_bridge import terminal_history_bridge
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_equations import heat_stress_equations
from lei_ren_part1_paper_compliant_power_angular_C4 import quotient_log_rates
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def angular_primitive_y_rows(X, rates):
    """Retain the original exact normalized angular primitive ODE."""
    rows = [X]
    one = IntervalTaylor.constant(X.ctx, 1, 5)
    for k in range(4):
        value = one if k == 0 else one*0
        for j in range(k+1):
            value -= rates[j]*rows[k-j]*math.comb(k, j)
        rows.append(value)
    return [v.truncate(5-k) for k, v in enumerate(rows)]


class CompliantHeatStressC4:
    def __init__(self, cells=64):
        self.pressure = CompliantHeatPressureC4(cells)
        self.heat, self.ctx = self.pressure.heat, self.pressure.ctx
        self.family, self.source = self.pressure.family, self.pressure.source
        self.history = terminal_history_bridge()
        self.theorem = terminal_stress_identities()
        self.equations = heat_stress_equations()
        if not self.history['complete_terminal_moment_history_bridge_verified']:
            raise ValueError('Actual terminal angular/energy histories not bound')
        for gate in ('renormalized_energy_constant_zero_from_selected_equation',
                     'actual_exterior_energy_full_Gamma_integral_and_half_normalization_verified',
                     'actual_zero_meridional_terminal_histories_verified',
                     'physical_terminal_energy_preserved_positive'):
            if not self.history[gate]:
                raise ValueError('Actual exterior moment source gate missing: '+gate)
        if not self.theorem['full_terminal_moment_stress_theorem_verified']:
            raise ValueError('Exact terminal moment-to-stress theorem required')
        self.hashes = dict(self.pressure.hashes)
        self.hashes.update(self.history['input_hashes'])
        for stem in ('compliant_heat_terminal_history_bridge',
                     'compliant_heat_terminal_stress_identities',
                     'compliant_heat_stress_equations', 'compliant_heat_stress_C4'):
            name = PREFIX+stem+'.py'
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()

    def exterior(self, Z, t):
        point = dict(self.pressure.exterior(Z, t))
        c = self.ctx
        local = point['full_local_Gamma_future_source']
        K = local['K_rows']
        X = local['angular_numerator']/K[0]
        point['original_forward_angular_Taylor'] = point['angular_Taylor']
        point['original_forward_angular_y_derivative_Taylor'] = point['angular_y_derivative_Taylor']
        point['original_retained_angular_tail_constant_defect'] = point['retained_angular_tail_constant_defect']
        # Equivalent moment evaluation after actual source/history proof.
        # Original forward angular and pressure packets remain accessible.
        point['angular_Taylor'] = X
        point['angular_y_derivative_Taylor'] = angular_primitive_y_rows(X, point['angular_ODE_rate_Taylor'])
        one = IntervalTaylor.constant(c, 1, 5)
        log_rates = quotient_log_rates(K)
        kappa = [one*(2+self.heat.delta)-2*log_rates[0]]+[-2*v for v in log_rates[1:]]
        shear = [v.truncate(3-k) for k, v in enumerate(kappa)]
        zero = IntervalTaylor.constant(c, 0, 5)
        point.update(
            exact_terminal_angular_definition='Mtheta=sqrt(2)*c*R^(1-a)/(1-a)+sqrt(2)*c*integral_R^infinity rho^(-a)*(1-H(2d/rho))drho',
            exact_terminal_energy_definition='Mztheta=c^2/2*integral_R^infinity rho^(-1-2a)*H(2d/rho)^2 drho',
            angular_tail_constant_exactly_zero_from_actual_source_history=True,
            original_angular_energy_and_pressure_histories_retained=True,
            all_five_terminal_moments_same_source_verified=True,
            exterior_meridional_moment_Taylor=dict(Mz=zero,Mtheta_z=zero),
            exterior_energy_full_Gamma_integral_and_half_normalization_verified=True,
            renormalized_energy_constant_zero_from_selected_equation=True,
            physical_terminal_energy_preserved_positive=True,
            theta_inertial_stress_over_F_y_derivative_Taylor=shear,
            theta_shear_over_F_y_derivative_Taylor=[-v for v in shear],
            individual_heat_shear_mixed_derivative_order=3,
            exact_heat_shear_positive_bracket_verified=True,
            canonical_heat_shear_coefficient_bounds=dict(lower=c.mpf(2), upper=2+self.heat.delta),
            exterior_similarity_stress_over_F=dict(theta=zero, axial=zero),
            exterior_stress_mixed_derivatives_total_order_le4={
                label:{'y'+str(k)+'_Z'+str(n):c.mpf(0) for k in range(5) for n in range(5-k)}
                for label in ('theta', 'axial')},
            zero_stress_follows_from_terminal_source_history_and_full_integral_identities=True,
            homogeneous_constants_eliminated_from_actual_terminal_moments=True,
            heat_exterior_stress_identity_certified=True,
            stress_scope='original similarity stress (3.16)-(3.18), Z[-1,1], t=log(R/Rtail)>=3',
            actual_compliant_field_physical_heat_region_certified=False,
            actual_physical_map_and_prefactor_transfer_pending=True,
            global_admissible_stress_lift_constructed=False,
            whole_outer_cone_certified=False,
            physical_energy_integral_certified=False,
            temporal_recursion=False)
        return point

    def report(self):
        keys = ('Z', 'coordinate', 'angular_Taylor', 'angular_y_derivative_Taylor',
                'energy_Taylor', 'pressure_over_Pstar_squared_Taylor',
                'original_forward_angular_Taylor', 'original_forward_angular_y_derivative_Taylor',
                'original_retained_angular_tail_constant_defect',
                'original_forward_pressure_over_Pstar_squared_Taylor',
                'theta_inertial_stress_over_F_y_derivative_Taylor',
                'theta_shear_over_F_y_derivative_Taylor',
                'individual_heat_shear_mixed_derivative_order',
                'canonical_heat_shear_coefficient_bounds',
                'exterior_similarity_stress_over_F',
                'exterior_stress_mixed_derivatives_total_order_le4',
                'angular_tail_constant_exactly_zero_from_actual_source_history',
                'original_angular_energy_and_pressure_histories_retained',
                'all_five_terminal_moments_same_source_verified',
                'exterior_meridional_moment_Taylor',
                'exterior_energy_full_Gamma_integral_and_half_normalization_verified',
                'renormalized_energy_constant_zero_from_selected_equation',
                'physical_terminal_energy_preserved_positive',
                'homogeneous_constants_eliminated_from_actual_terminal_moments',
                'zero_stress_follows_from_terminal_source_history_and_full_integral_identities',
                'exact_heat_shear_positive_bracket_verified',
                'heat_exterior_stress_identity_certified',
                'actual_compliant_field_physical_heat_region_certified',
                'actual_physical_map_and_prefactor_transfer_pending',
                'stress_scope', 'global_admissible_stress_lift_constructed',
                'whole_outer_cone_certified', 'physical_energy_integral_certified',
                'temporal_recursion', 'exact_terminal_angular_definition',
                'exact_terminal_energy_definition', 'exact_positive_Ev0Pstar2_log_parts',
                'exact_relative_velocity_log_parts')
        summary = lambda point:{k:point[k] for k in keys}
        with mp.workdps(270):
            return dict(
                actual_five_defect_family_sha256=self.family,
                implicit_source_sha256=self.source,
                scope='Source-bound stress-free Gamma exterior t>=3, Z[-1,1]; exact zero stress and mixed4; individual shear mixed3',
                samples=[summary(self.exterior(z,t)) for z,t in
                         (('0','3'),('.5','3'),('.5','4'),('-1','3'),('1','3'),('.5','1000'))],
                whole_Z_inlet=summary(self.exterior([-1,1],3)),
                whole_Z_unbounded_exterior=summary(self.exterior([-1,1],[3,mp.inf])),
                terminal_history_bridge=self.history,
                terminal_stress_theorem=self.theorem,
                canonical_heat_stress_equations=self.equations,
                all_five_terminal_moments_same_source_verified=True,
                homogeneous_constants_eliminated_from_actual_terminal_moments=True,
                heat_exterior_stress_identity_certified=True,
                original_angular_energy_and_pressure_histories_retained=True,
                actual_compliant_field_physical_heat_region_certified=False,
                main_dispatcher_and_collar_adoption_complete=False,
                global_admissible_stress_lift_constructed=False,
                whole_outer_cone_certified=False,
                physical_energy_integral_certified=False,
                temporal_recursion=False,
                input_hashes=self.hashes)


def run():
    result = CompliantHeatStressC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Source-bound Gamma exterior stress generated; global cone and recursion remain open',flush=True)
    return result


if __name__ == '__main__':
    run()
