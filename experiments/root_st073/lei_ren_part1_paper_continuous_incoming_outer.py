"""Shared continuous incoming velocity/mean installed before the axial pulse.

Regenerate reference axial inputs, reapply measured inner offsets once, and
reuse the continuous exterior integral atoms for a new coefficient solve.
Swirl reference energy is regenerated from shared angular propagation;
quadrature, input uncertainty and complete exterior closure remain open.
"""
from copy import deepcopy
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp

from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_tangents import ContinuousAxialAlgebra
from lei_ren_part1_paper_continuous_seeded_outer_jets import (
    ContinuousSeededAxialProfileJets, ContinuousSeededOuterFieldJets)
from lei_ren_part1_paper_seeded_axial_inputs import seeded_axial_inputs
from lei_ren_part1_paper_seeded_outer_field import SeededAxialProfile
from lei_ren_part1_paper_joined_outer import build_joined_field, _mp


def regenerate_incoming(source, seeded, atoms, provider, *, order=96):
    """Replace three reference axial integrals and solve with retained atoms."""
    schedule=source.schedule
    precision=source.precision
    with mp.workdps(precision):
        result=deepcopy(seeded)
        incoming=deepcopy(seeded['incoming'])
        z=mp.mpf(str(incoming['Z']))
        # All normalized rows must use the newly installed amplitude at Rp.
        incoming['log_Ep']=mp.nstr(_mp(source.schedule.at_log_radius(
            str(schedule.logR_p),z)['log_angular_amplitude']),precision)
        from lei_ren_part1_paper_continuous_incoming_angular import ContinuousIncomingAngular
        angular=ContinuousIncomingAngular(str(schedule.logPstar),str(schedule.Md),
                                          precision=provider.precision)
        mixed_factor=angular.mixed_factor(order=order)
        full=provider.full_incoming_rows(z)
        dimensions=incoming['dimensionless_integrals']
        old_dimensions=deepcopy(dimensions)
        from lei_ren_part1_paper_continuous_incoming_swirl import ContinuousIncomingSwirl
        if not hasattr(source,'_continuous_incoming_swirl_reference'):
            source._continuous_incoming_swirl_reference=ContinuousIncomingSwirl(source,order=order)
        swirl=source._continuous_incoming_swirl_reference.reference(z)
        # Serializing at source algebra precision cannot add information to
        # the incoming quadrature atoms. Retain their declared precision.
        dimensions.update(I_z=mp.nstr(full['I_z'],provider.precision),
                          I_uz2=mp.nstr(full['I_uz2'],provider.precision),
                          I_theta_z=mp.nstr(mixed_factor*z/(1+z*z),provider.precision),
                          I_swirl=mp.nstr(swirl['I_swirl'],precision))
        yp=mp.mpf(str(schedule.y_p));ep=mp.mpf(incoming['log_Ep'])
        mu=mp.mpf(str(schedule.mu))
        reference=[mp.mpf(dimensions['I_z'])*mp.exp(-yp-ep),
                   mp.mpf(dimensions['I_theta_z'])*mp.exp(-mp.mpf('1.5')*yp-2*ep)]
        prior=mu*(mp.mpf(dimensions['I_uz2'])-mp.mpf(dimensions['I_swirl']))*mp.exp(-yp-2*ep)
        for key,value in zip(('m1_Mz_over_RpEp','m2_Mtheta_z_over_Rp_sqrt2RpEp2'),reference):
            incoming[key]=signed_log(value,precision)
        incoming['E_prior_mu_Mztheta_over_RpEp2']=signed_log(prior,precision)
        offsets={key:from_signed_log(value) for key,value in
                 seeded['incoming']['inner_seed_transport']['raw_offsets'].items()}
        incoming=seeded_axial_inputs(incoming,log_Rp=str(schedule.logR_p),mu=str(schedule.mu),
                                     offsets=offsets,precision=precision)
        new_prior=from_signed_log(incoming['E_prior_mu_Mztheta_over_RpEp2'])
        future=source.outer.tail.evaluate(float(z),quadrature_order=source.outer.order)
        future_nominal=mp.mpf(future['energy_target_contribution_nominal'])
        target=(1-mp.exp(-26))/4-new_prior+future_nominal
        incoming['continuous_incoming']={
            'precision':provider.precision,'mixed_quadrature_order':order,
            'reference_linear_factors':{
                'I_z_over_Z':mp.nstr(provider.full_incoming_rows(1)['I_z'],provider.precision),
                'I_theta_z_times_1plusZ2_over_Z':mp.nstr(mixed_factor,provider.precision)},
            'old_dimensionless_integrals':old_dimensions,
            'dimensionless_integral_working_precision':{
                'I_z':provider.precision,'I_uz2':provider.precision,
                'I_theta_z':provider.precision,'I_swirl':precision},
            'mixed_axial_factor_working_precision':provider.precision,
            'serialized_digits_are_not_integral_accuracy':True,
            'angular_primitive_float_backed':False,'swirl_energy_inherited':False,
            'swirl_reference_shared_angular_propagation':True,
            'swirl_quadrature_order':order,'swirl_quadrature_error_enclosed':False,
            'angular_heat_kernel_inherited':False,
            'heat_integral_targets_inherited_Taylor':True,
            'future_energy_regenerated_from_live_tail':True,
            'future_energy_nominal':mp.nstr(future_nominal,precision),
            'previous_energy_target':seeded['energy_target'],
            'future_energy_uncertainty_fully_enclosed':False,
            'angular_correction_and_pressure_moments_complete':False,
            'mixed_reference_primitive_float_backed':False,
            'mixed_reference_primitive_definition':'ContinuousIncomingAngular.J with exact J(1)=1/2',
            'mixed_reference_logPstar':str(schedule.logPstar),
            'mixed_reference_Md':str(schedule.Md),
            'quadrature_enclosure_certified':False,'inner_offsets_reapplied_once':True}
        result['incoming']=incoming
        result['energy_target']=mp.nstr(target,precision)
        result['axial']['energy_target']=result['energy_target']
        base=[incoming['row_normalization'][key] for key in ('scaled_base_m1','scaled_base_m2')]
        result['axial']['linear_rhs_inputs']['base']=base
        algebra=ContinuousAxialAlgebra(result,atoms)
        with mp.workdps(algebra.precision):
            a,c=algebra.solve()
            updated=deepcopy(atoms)
            # These diagnostics belong to the previous coefficient solve,
            # not to the retained integration atoms. Do not relabel them.
            for key in ('coefficient_relative_changes','amplitude_relative_change',
                        'old_coefficients_continuous_row_relative_defects',
                        'continuous_component_cumulative_rows_after_support',
                        'end_component_primitive_derivative_errors'):
                updated.pop(key,None)
            updated.update(a_p=mp.nstr(a,algebra.precision),
                           c=[signed_log(v,algebra.precision) for v in c],
                           incoming_rows_regenerated=True,
                           inherited_base_rows=False,inherited_energy_target=False,
                           installed_in_global_profile=True,global_mean_closed=False,
                           energy_target=mp.nstr(algebra.target,algebra.precision))
            errors=[]
            for i in range(2):
                rhs=algebra.base[i]+a*algebra.p[i]
                errors.append(mp.nstr(abs((sum(algebra.M[i][j]*c[j] for j in range(2))+rhs)/rhs),40))
            energy=algebra.Kp*a*a+algebra.mu*sum(k*v*v for k,v in zip(algebra.K,c))
            updated['linear_relative_replay']=errors
            updated['energy_relative_replay']=mp.nstr(abs(energy/algebra.target-1),40)
        return result,updated


class ContinuousIncomingProfile(ContinuousSeededAxialProfileJets):
    def __init__(self,source,*,prepared=None,incoming_precision=100,order=96):
        from lei_ren_part1_paper_continuous_angular_schedule import install_continuous_angular_schedule
        self.angular_schedule_provider=install_continuous_angular_schedule(
            source.schedule,precision=source.precision,primitive_precision=incoming_precision)
        from lei_ren_part1_paper_continuous_angular_correction import install_continuous_angular_correction
        self.angular_correction_provider=install_continuous_angular_correction(
            source.outer.angular,precision=source.precision,basis_precision=incoming_precision)
        # Consumers share the correction object. Discard receipts computed
        # before replacing its atoms, rather than mixing two bump definitions.
        source.outer.tail._coefficient_cache.clear()
        source.outer.pressure._coefficients_cache.clear()
        source.outer.pressure._logU_rel=source.schedule.at_log_radius(
            source.schedule.logR_rel,0)['log_angular_amplitude']
        source.outer.coefficients.cache_clear()
        from lei_ren_part1_paper_continuous_incoming import ContinuousIncomingAxial
        self.incoming_provider=ContinuousIncomingAxial(str(source.schedule.Md),precision=incoming_precision)
        self.incoming_order=order
        revised={z:regenerate_incoming(source,seeded,atoms,self.incoming_provider,order=order)
                 for z,(seeded,atoms) in (prepared or {}).items()}
        super().__init__(source,prepared=revised)
        self.shared_runtimes={}
        from lei_ren_part1_paper_continuous_angular_moments import ContinuousAngularMoments
        self.angular_moment_provider=ContinuousAngularMoments(self,order=order)
        from lei_ren_part1_paper_continuous_mixed_moments import ContinuousMixedMoments
        self.mixed_moment_provider=ContinuousMixedMoments(self)
        from lei_ren_part1_paper_continuous_axial_energy_moments import ContinuousAxialEnergyMoments
        self.axial_energy_moment_provider=ContinuousAxialEnergyMoments(self)

    def angular_moments_jet(self,logR,Z):
        return self.angular_moment_provider.moments_jet(logR,Z)

    def linear_axial_moments_jet(self,logR,Z):
        return self.mixed_moment_provider.moments_jet(logR,Z)

    def pressure_moments_jet(self,logR,Z):
        if not hasattr(self,'pressure_moment_provider'):
            from lei_ren_part1_paper_continuous_pressure_moments import ContinuousPressureMoments
            self.pressure_moment_provider=ContinuousPressureMoments(self,order=self.incoming_order)
        return self.pressure_moment_provider.moments_jet(logR,Z)

    def quadratic_moments_jet(self,logR,Z):
        if self.offset(logR,self.schedule.logR_v)<0:
            if not hasattr(self,'partial_axial_energy_moment_provider'):
                from lei_ren_part1_paper_continuous_partial_axial_energy import ContinuousPartialAxialEnergy
                self.partial_axial_energy_moment_provider=ContinuousPartialAxialEnergy(self,quadrature_order=self.incoming_order)
            return self.partial_axial_energy_moment_provider.moments_jet(logR,Z)
        return self.axial_energy_moment_provider.moments_jet(logR,Z)

    def _ensure_prepared(self,Z):
        z=float(Z)
        if z not in self.prepared:
            from lei_ren_part1_paper_continuous_axial_solve import solve_continuous
            old=SeededAxialProfile.seed_solve(self,z)
            atoms=solve_continuous(old,precision=200)
            regenerated=regenerate_incoming(self.seed_source,old,atoms,self.incoming_provider,order=self.incoming_order)
            self.prepared[z]=regenerated
        return self.prepared[z]

    def runtime(self,Z):
        z=float(Z)
        if z not in self.shared_runtimes:
            from lei_ren_part1_paper_continuous_axial_runtime import SharedContinuousAxialRuntime
            seeded,atoms=self._ensure_prepared(z)
            runtime=SharedContinuousAxialRuntime(seeded,atoms)
            self.shared_runtimes[z]=runtime
            self.prepared[z]=(seeded,runtime.receipt)
        return self.shared_runtimes[z]

    @lru_cache(maxsize=64)
    def seed_solve(self,Z):
        runtime=self.runtime(Z)
        seeded,_=self.prepared[float(Z)]
        return dict(seeded,axial=runtime.receipt)

    def component(self,Z):
        return self.runtime(Z).component

    @lru_cache(maxsize=64)
    def coefficient_tangent(self,Z):
        from lei_ren_part1_paper_seeded_input_tangents import actual_seeded_input_tangents
        runtime=self.runtime(Z)
        seeded,_=self.prepared[float(Z)]
        input_jet=actual_seeded_input_tangents(self.seed_source,seeded,float(Z))
        with mp.workdps(max(self.precision,self.jet_precision,runtime.algebra.precision)):
            tangent=runtime.algebra.tangent(
                [from_signed_log(value) for value in input_jet['base_Z']],
                from_signed_log(input_jet['energy_target_Z']))
        return dict(a_Z=tangent['a_Z'],c_Z=tangent['c_Z'],input=input_jet,
            algebra_precision=runtime.algebra.precision,
            linear_tangent_relative_replay=tangent['linear_tangent_relative_replay'],
            energy_tangent_relative_replay=tangent['energy_tangent_relative_replay'],
            shared_complete_atoms_installed=True)

    def _pulse_integral_and_derivative(self,component,xi,tangent):
        if _mp(xi)>=11:
            balance=component.terminal_balance(1,max(self.precision,self.jet_precision))
            return balance['pulse_integral'],mp.mpf(0),balance['pulse_atom']
        return super()._pulse_integral_and_derivative(component,xi,tangent)

    def _end_primitive_derivative(self,component,c_Z,end):
        if all(_mp(end)-center>=component.basis.ell for center in (-3,-1)):
            # The complete end moment is the very matrix used in the solve.
            return component.runtime.full_end(c_Z,1)
        return super()._end_primitive_derivative(component,c_Z,end)

    def _incoming_y(self,logR):
        return str(self.offset(logR,self.schedule.logRref))

    def values(self,logR,Z):
        result=super().values(logR,Z)
        self._install_angular_jet(result,logR,Z)
        if self.offset(logR,self.schedule.logR_p)<0:
            result['Uz']=self.incoming_provider.values(self._incoming_y(logR),str(Z))['Uz']
            result['continuous_incoming_values_installed']=True
        return result

    def values_with_jets(self,logR,Z):
        result=super().values_with_jets(logR,Z)
        self._install_angular_jet(result,logR,Z)
        if self.offset(logR,self.schedule.logR_p)<0:
            values=self.incoming_provider.values(self._incoming_y(logR),str(Z))
            result.update(Uz=values['Uz'],Uz_Z=values['Uz_Z'],
                          axial_Z_method='continuous_incoming_analytic_Z')
        return result

    def _install_angular_jet(self,result,logR,Z):
        """Differentiate the same multiplier used in the angular velocity."""
        with mp.workdps(max(self.precision,self.jet_precision)):
            base=self.schedule.at_log_radius(logR,Z)
            offset=_mp(str(self.offset(logR,self.schedule.logR_rel)))
            jet=self.angular_correction_provider.value_jet(offset,Z)
            amplitude=mp.exp(_mp(base['log_angular_amplitude']))
            multiplier=1+jet['h']
            LZ=_mp(base['dlogU_dZ'])
            slope=_mp(base['logarithmic_slope'])+jet['h_t']/multiplier
            result.update(Utheta=amplitude*multiplier,
                Utheta_Z=amplitude*(LZ*multiplier+jet['h_Z']),
                logarithmic_slope=slope,logF_slope=slope-mp.mpf('.5'),
                angular_relative_correction=jet['h'],
                angular_relative_correction_Z=jet['h_Z'],
                angular_relative_correction_t=jet['h_t'],
                angular_log1p_correction=mp.nstr(mp.log1p(jet['h']),self.precision),
                angular_Z_method='continuous_schedule_and_implicit_bump_coefficients',
                angular_correction_Z_jet_installed=True,
                angular_Z_jet_complete=False)
            for key in ('heat_deficit','heat_log_amplitude_correction','heat_logarithmic_slope_correction',
                        'heat_kernel_method','heat_truncation_absolute_bounds','heat_integral_targets_regenerated'):
                if key in base:result[key]=base[key]
            if 'heat_log_amplitude_correction' in base:
                result['angular_heat_relative_correction']=mp.expm1(base['heat_log_amplitude_correction'])
            # These local jets differentiate the installed numerical input
            # model. Heat/baseline uncertainty and full moments remain open.

    values_Z=values_with_jets

    def axial_average_jet(self,logR,Z):
        if self.offset(logR,self.schedule.logR_p)>=0:
            result=super().axial_average_jet(logR,Z)
            result['shared_complete_atoms_installed']=True
            return result
        with mp.workdps(self.precision):
            local=self.incoming_provider.mean_jet(self._incoming_y(logR),str(Z))
            offsets=self.seed_source.terminal_offsets(self.seed_source._z_key(_mp(Z)))
            inverseR=mp.exp(_mp(logR))**-1
            return dict(value=local['mean']+offsets['mass_offset']*inverseR,
                        derivative=local['mean_Z']+offsets['mass_offset_Z']*inverseR,
                        method='continuous_incoming_plus_actual_inner_seed',
                        terminal_mean_forced_zero=False,
                        incoming_uncertainty_enclosed=False)


class ContinuousIncomingOuterField(ContinuousSeededOuterFieldJets):
    def __init__(self,source,*,prepared=None,incoming_precision=100,order=96):
        super().__init__(source,prepared=prepared)
        self.outer=ContinuousIncomingProfile(source,prepared=prepared,
            incoming_precision=incoming_precision,order=order)

    def _outer_average(self,log_radius,z):
        jet=self.outer.axial_average_jet(mp.nstr(log_radius,self.precision),z)
        return jet['value'],jet['derivative']


def run():
    folder=Path(__file__).parent
    print('building retained joined field',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    print('regenerating continuous incoming rows and coefficient solve',flush=True)
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    field.outer.runtime(.3)
    revised,solved=field.outer.prepared[.3]
    with mp.workdps(field.precision):
        z=mp.mpf('.3');logRp=_mp(str(source.schedule.logR_p))
        # Incoming mass is constant beyond cutoff. Compare the two formulas
        # at the identical Rp coordinate, avoiding an artificial terminal zero.
        before=field.outer.axial_average_jet(mp.nstr(logRp-1,field.precision),z)
        at=field.outer.axial_average_jet(mp.nstr(logRp,field.precision),z)
        jumps=[abs(at[key]*mp.e/before[key]-1) for key in ('value','derivative')]
        samples=[]
        for y in ('-5','0','1.3','2'):
            radius=_mp(str(source.schedule.logRref))+mp.mpf(y)
            text=mp.nstr(radius,field.precision)
            velocity=field.outer.values_with_jets(text,z)
            mean=field.outer.axial_average_jet(text,z)
            samples.append(dict(y=y,Uz=signed_log(velocity['Uz'],field.precision),
                Uz_Z=signed_log(velocity['Uz_Z'],field.precision),
                mean=signed_log(mean['value'],field.precision),
                mean_Z=signed_log(mean['derivative'],field.precision),method=mean['method']))
        angular_checks=[]
        for t in ('-3','-2.96','-1','-.96','0'):
            radius=_mp(str(source.schedule.logR_rel))+mp.mpf(t)
            text=mp.nstr(radius,field.precision)
            velocity=field.outer.values_with_jets(text,z)
            base=source.schedule.at_log_radius(text,z)
            amp=mp.exp(_mp(base['log_angular_amplitude']))
            correction=field.outer.angular_correction_provider.value_jet(mp.mpf(t),z)
            expected=amp*(_mp(base['dlogU_dZ'])*(1+correction['h'])+correction['h_Z'])
            error=abs((velocity['Utheta_Z']-expected)/expected) if expected else abs(velocity['Utheta_Z'])
            if error>mp.mpf('1e-70'):
                raise ArithmeticError('Actual angular correction Z jet differs from its provider')
            if t=='-3' and not velocity['Utheta_Z']:
                raise ArithmeticError('Actual field discarded the nonzero angular correction Z derivative')
            angular_checks.append(dict(t=t,Utheta_Z=signed_log(velocity['Utheta_Z'],field.precision),
                relative_correction_Z=signed_log(correction['h_Z'],field.precision),
                provider_relative_error=mp.nstr(error,40)))
        pressure_coeff=field.outer.angular.coefficients(z)
        pressure_total,pressure_meta=field.outer.pressure._bump_integrals(float(z),pressure_coeff,order=96)
        pressure_half,half_meta=field.outer.pressure._bump_integrals(float(z),pressure_coeff,order=96,upper_t=-3.)
        if not pressure_half:
            raise ArithmeticError('Partial pressure correction lost the first signed bump')
        pressure_jet=field.outer.pressure.continuous_bump_jet('-2',z)
        dz=mp.mpf('1e-5')
        pressure_fd=sum(weight*field.outer.pressure.continuous_bump_jet('-2',z+shift*dz)['value']
                        for shift,weight in ((-2,1),(-1,-8),(1,8),(2,-1)))/(12*dz)
        pressure_jet_error=abs(pressure_fd/pressure_jet['derivative']-1)
        if pressure_jet_error>mp.mpf('1e-13'):
            raise ArithmeticError('Continuous bump pressure Z derivative failed independent difference')
        logUv=source.schedule.at_log_radius(source.schedule.logR_v,z)['log_angular_amplitude']
        bump_energy,energy_meta=field.outer.tail._angular_bump_correction(float(z),logUv,pressure_coeff,96)
        if not (pressure_meta['continuous_bump_atoms_installed'] and energy_meta['continuous_bump_atoms_installed']):
            raise ArithmeticError('Pressure or energy still uses legacy angular bump quadrature')
        report=dict(Z='.3',incoming=revised['incoming'],continuous_solve=solved,
            Rp_mass_relative_matching=mp.nstr(jumps[0],40),
            Rp_mass_Z_relative_matching=mp.nstr(jumps[1],40),samples=samples,
            installed_in_velocity_and_mean=True,terminal_mean_forced_zero=False,
            shared_complete_atoms_installed=True,
            source_schedule_identity_preserved=field.outer.schedule is source.schedule,
            angular_primitive_float_backed=False,swirl_energy_inherited=True,
            angular_preheat_schedule_installed=True,angular_heat_kernel_inherited=False,
            heat_point_kernel_installed=True,heat_integral_targets_inherited_Taylor=True,
            angular_correction_checks=angular_checks,
            angular_correction_Z_jet_installed=True,
            source_correction_identity_preserved=field.outer.angular is source.outer.angular,
            pressure_correction_same_bump_atoms=hasattr(field.outer.pressure.correction,'_continuous_provider'),
            pressure_full_increment=signed_log(pressure_total,field.precision),
            pressure_first_half_increment=signed_log(pressure_half,field.precision),
            angular_bump_energy=signed_log(bump_energy,field.precision),
            angular_bump_energy_same_atoms=True,
            angular_bump_pressure_Z=signed_log(pressure_jet['derivative'],field.precision),
            angular_bump_pressure_Z_difference_error=mp.nstr(pressure_jet_error,40),
            angular_pressure_and_moments_complete=False,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:report[key] for key in ('Rp_mass_relative_matching',
          'Rp_mass_Z_relative_matching','installed_in_velocity_and_mean')}),flush=True)
    return report


if __name__=='__main__':run()
