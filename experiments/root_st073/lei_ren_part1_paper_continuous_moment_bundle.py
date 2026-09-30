"""Actual five-moment and pressure bundle for one continuous candidate.

Assembly does not prove moment closure, pressure compatibility, finite energy
or admissible stress. Component rows stay available alongside rounded sums.
Provider component checks do not establish coupled moment closure.
"""
import mpmath as mp
from lei_ren_part1_paper_continuous_velocity_radial_jets import velocity_radial_jets
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


class ContinuousMomentBundle:
    def __init__(self,profile,*,pressure_provider=None):
        self.profile=profile;self.precision=profile.precision
        self.pressure_provider=pressure_provider

    def evaluate(self,logR,Z):
        with mp.workdps(self.precision):
            p=self.profile
            angular=p.angular_moments_jet(logR,Z)
            linear=p.linear_axial_moments_jet(logR,Z)
            quadratic=p.quadratic_moments_jet(logR,Z)
            pressure=(p.pressure_moments_jet(logR,Z) if self.pressure_provider is None
                      else self.pressure_provider.moments_jet(logR,Z))
            moments=dict(theta=angular['theta'],z=linear['z'],
                theta_z=linear['theta_z'],z_theta=quadratic['z_theta'],p=pressure['p'])
            jets=dict(theta=angular['theta_Z'],z=linear['z_Z'],
                theta_z=linear['theta_z_Z'],z_theta=quadratic['z_theta_Z'],p=pressure['p_Z'])
            return dict(moments=moments,moments_Z=jets,P=pressure['P'],P_Z=pressure['PZ'],
                components=dict(angular=angular,linear=linear,quadratic=quadratic,pressure=pressure),
                all_five_actual_moments_evaluated=True,shared_candidate=True,
                pressure_propagated_from_actual_inner=True,
                terminal_moments_forced_zero=False,moment_closure_certified=False,
                finite_energy_certified=False,scale_recursion_certified=False)

    def stress(self,logR,Z):
        with mp.workdps(self.precision):
            row=self.evaluate(logR,Z)
            velocity=velocity_radial_jets(self.profile,logR,Z)
            stress=evaluate_mp_stress(logR,Z,self.profile.schedule.delta,
                Utheta=velocity['Utheta'],Uz=velocity['Uz'],
                Utheta_y=velocity['Utheta_y'],Utheta_Z=velocity['Utheta_Z'],
                Uz_y=velocity['Uz_y'],Uz_Z=velocity['Uz_Z'],
                moments=row['moments'],moments_Z=row['moments_Z'],
                P=row['P'],P_Z=row['P_Z'],precision=self.precision,
                include_components=True)
            components=stress.pop('components')
            angular_matching=stress.pop('angular_matching')
            return dict(row,velocity=velocity,stress=stress,
                stress_components=components,angular_matching=angular_matching,
                stress_derivatives_from_shared_continuous_jets=True,
                tiny_stress_terms_separately_resolved=False,
                stress_cone_certified=False,stress_remainder_decomposed=False,
                oscillatory_correction_installed=False,full_NS_residual_certified=False)
