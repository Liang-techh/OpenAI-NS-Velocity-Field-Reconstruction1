"""Two angular-bump targets from the shared actual primitive and heat atoms.

Pressure restoration targets the analytic H=1 preheat datum of v2 (6.10).
It does not use a late pressure gauge shift to repair an incompatible core.
"""
import mpmath as mp
from lei_ren_part1_paper_continuous_angular_target import _tail_terms


class CoupledAngularTargets:
    def __init__(self, profile):
        self.profile = profile
        self.schedule = profile.schedule
        self.precision = profile.precision
        self.angular = profile.angular_moment_provider
        self.heat = self.angular.heat_provider

    def evaluate(self, Z, *, include_inner=True):
        with mp.workdps(self.precision):
            z = mp.mpf(str(Z))
            s = self.schedule
            if abs(z) >= 1:
                raise ValueError('Actual inner targets require |Z|<1')
            tail_log = mp.mpf(str(s.logR_tail))
            rel_log = mp.mpf(str(s.logR_rel))
            key = mp.nstr(tail_log, self.precision)
            z_key = mp.nstr(z, self.precision)
            y_tail = mp.mpf(str(s.y_tail))
            baseline = self.angular._baseline(mp.nstr(y_tail,self.precision),z_key)
            offsets = (self.angular._inner_offsets(z_key) if include_inner else
                       [mp.mpf(0)]*4)
            collar = self.heat.increments(3, z)
            heat_point = self.profile.angular_schedule_provider.heat_jet(s.logR_tail,z)
            terms = _tail_terms(delta=s.delta, log_r_tail=tail_log,
                log_c_inf=s._log_c_inf, xi=heat_point['xi'], xi_z=heat_point['xi_Z'],
                collar_ref_theta=collar[0], collar_heat_correction_theta=collar[2],
                collar_heat_correction_theta_z=collar[4])
            e_rel = self.angular._base(str(s.y_rel))
            theta_scale = mp.sqrt(2)*mp.exp(mp.mpf('1.5')*rel_log+e_rel)
            theta_ratio = terms['scale']/theta_scale
            pressure_scale = mp.exp(2*e_rel)
            # Remove the current bump by using baseline/offset components,
            # rather than subtracting two materialized cumulative moments.
            # Baseline uses the actual tail amplitude Cinf*Rtail^-a*(1-epsilon),
            # whereas heat atoms use Cinf*Rtail^-a. Retain the conversion.
            tail_factor=1-mp.mpf(str(s.epsilon))
            r_parts = dict(
                preheat=-theta_ratio*mp.fsum([tail_factor*baseline[0],terms['normalized_reference']]),
                heat=-theta_ratio*terms['normalized_heat_correction'],
                inner=-offsets[0]/theta_scale)
            direct_preheat=r_parts['preheat']
            coherent=bool(getattr(s,'_coherent_preheat_waiting',False))
            if coherent:
                yt=mp.mpf(str(s.y_t));yr=mp.mpf(str(s.y_rel))
                et=self.angular._base(mp.nstr(yt,self.precision))
                ratio_t=mp.exp(mp.mpf('1.5')*(yt-yr)+et-e_rel)
                delta_mean=self.angular.angular_difference_from_axis(mp.nstr(yt,self.precision),z_key)
                r_parts['preheat']=-ratio_t*delta_mean
            rz_parts = dict(preheat=-theta_ratio*tail_factor*baseline[2],
                heat=-theta_ratio*terms['normalized_heat_correction_Z'],
                inner=-offsets[2]/theta_scale)
            square = self.heat.complete_pressure_heat_integral(z)
            s_target = -square['heat_correction']/pressure_scale
            sz_target = -square['heat_correction_Z']/pressure_scale
            return dict(Z=z,r=mp.fsum(r_parts.values()),rZ=mp.fsum(rz_parts.values()),
                s=s_target,sZ=sz_target,r_components=r_parts,rZ_components=rz_parts,
                theta_scale=theta_scale,pressure_scale=pressure_scale,
                current_angular_bump_excluded=True,
                coherent_zero_axis_waiting_constraint_used=coherent,
                direct_forward_preheat_target=direct_preheat,
                actual_moment_overwritten=False,
                analytic_preheat_pressure_target=True,
                pressure_target_definition='-complete heat pressure correction / Urel^2',
                actual_inner_offsets_included=include_inner,
                heat_pressure_truncation_bound=square['analytic_heat_truncation_absolute_bound']/pressure_scale,
                heat_model='shared quadratic small-x kernel; exact remainder not enclosed in Z',
                actual_core_uses_complete_preheat_datum=False,
                quadrature_error_enclosed=False,terminal_matching_certified=False)
