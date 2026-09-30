"""Shared continuous incoming velocity/mean installed before the axial pulse.

Regenerate reference axial inputs, reapply measured inner offsets once, and
reuse the continuous exterior integral atoms for a new coefficient solve.
Angular primitives and swirl energy retain their explicitly inherited scope.
"""
from copy import deepcopy
from decimal import Decimal
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
        # The source swirl factor is A(y)/(1+Z^2). Its original schedule
        # primitive remains float-backed; only the axial atom is replaced.
        nodes,weights=mp.gauss_quadrature(order,'legendre')
        def integrate(a,b):
            half=(b-a)/2;center=(a+b)/2
            total=mp.mpf(0)
            for node,weight in zip(nodes,weights):
                y=center+half*node
                logA=mp.mpf(str(schedule._log_A(Decimal(mp.nstr(y,precision)))))
                total+=weight*mp.exp(mp.mpf('1.5')*y+logA)*provider.cutoff(y)
            return half*total
        cutoff_end=mp.exp(mp.mpf(str(schedule.Md)))
        mixed_factor=4*(mp.exp(mp.mpf(str(schedule.logPstar)))/mp.mpf('1.6')+
                        integrate(mp.mpf(0),mp.mpf(1))+integrate(mp.mpf(1),cutoff_end))
        full=provider.full_incoming_rows(z)
        dimensions=incoming['dimensionless_integrals']
        old_dimensions=deepcopy(dimensions)
        # Serializing at source algebra precision cannot add information to
        # the incoming quadrature atoms. Retain their declared precision.
        dimensions.update(I_z=mp.nstr(full['I_z'],provider.precision),
                          I_uz2=mp.nstr(full['I_uz2'],provider.precision),
                          I_theta_z=mp.nstr(mixed_factor*z/(1+z*z),precision))
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
        old_prior=from_signed_log(seeded['incoming']['E_prior_mu_Mztheta_over_RpEp2'])
        new_prior=from_signed_log(incoming['E_prior_mu_Mztheta_over_RpEp2'])
        target=mp.mpf(seeded['energy_target'])+old_prior-new_prior
        incoming['continuous_incoming']={
            'precision':provider.precision,'mixed_quadrature_order':order,
            'reference_linear_factors':{
                'I_z_over_Z':mp.nstr(provider.full_incoming_rows(1)['I_z'],provider.precision),
                'I_theta_z_times_1plusZ2_over_Z':mp.nstr(mixed_factor,precision)},
            'old_dimensionless_integrals':old_dimensions,
            'dimensionless_integral_working_precision':{
                'I_z':provider.precision,'I_uz2':provider.precision,
                'I_theta_z':precision,'I_swirl':seeded['incoming']['precision']},
            'mixed_axial_factor_working_precision':provider.precision,
            'serialized_digits_are_not_integral_accuracy':True,
            'angular_primitive_float_backed':True,'swirl_energy_inherited':True,
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
        from lei_ren_part1_paper_continuous_incoming import ContinuousIncomingAxial
        self.incoming_provider=ContinuousIncomingAxial(str(source.schedule.Md),precision=incoming_precision)
        self.incoming_order=order
        revised={z:regenerate_incoming(source,seeded,atoms,self.incoming_provider,order=order)
                 for z,(seeded,atoms) in (prepared or {}).items()}
        super().__init__(source,prepared=revised)

    @lru_cache(maxsize=64)
    def seed_solve(self,Z):
        z=float(Z)
        if z in self.prepared:
            regenerated=self.prepared[z]
        else:
            from lei_ren_part1_paper_continuous_axial_solve import solve_continuous
            old=SeededAxialProfile.seed_solve(self,z)
            atoms=solve_continuous(old,precision=200)
            regenerated=regenerate_incoming(self.seed_source,old,atoms,self.incoming_provider,order=self.incoming_order)
            self.prepared[z]=regenerated
        return dict(regenerated[0],axial=regenerated[1])

    def _incoming_y(self,logR):
        return str(self.offset(logR,self.schedule.logRref))

    def values(self,logR,Z):
        result=super().values(logR,Z)
        if self.offset(logR,self.schedule.logR_p)<0:
            result['Uz']=self.incoming_provider.values(self._incoming_y(logR),str(Z))['Uz']
            result['continuous_incoming_values_installed']=True
        return result

    def values_with_jets(self,logR,Z):
        result=super().values_with_jets(logR,Z)
        if self.offset(logR,self.schedule.logR_p)<0:
            values=self.incoming_provider.values(self._incoming_y(logR),str(Z))
            result.update(Uz=values['Uz'],Uz_Z=values['Uz_Z'],
                          axial_Z_method='continuous_incoming_analytic_Z')
        return result

    values_Z=values_with_jets

    def axial_average_jet(self,logR,Z):
        if self.offset(logR,self.schedule.logR_p)>=0:
            return super().axial_average_jet(logR,Z)
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
        report=dict(Z='.3',incoming=revised['incoming'],continuous_solve=solved,
            Rp_mass_relative_matching=mp.nstr(jumps[0],40),
            Rp_mass_Z_relative_matching=mp.nstr(jumps[1],40),samples=samples,
            installed_in_velocity_and_mean=True,terminal_mean_forced_zero=False,
            source_schedule_identity_preserved=field.outer.schedule is source.schedule,
            angular_primitive_float_backed=True,swirl_energy_inherited=True,
            finite_energy_certified=False,scale_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:report[key] for key in ('Rp_mass_relative_matching',
          'Rp_mass_Z_relative_matching','installed_in_velocity_and_mean')}),flush=True)
    return report


if __name__=='__main__':run()
