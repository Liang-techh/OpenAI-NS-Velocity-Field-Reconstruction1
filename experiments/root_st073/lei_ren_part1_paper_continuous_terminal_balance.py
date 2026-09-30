"""Actual terminal mass/jet balance and the radial energy tail it implies.

This evaluates the installed continuous incoming/exterior, not the earlier
quadrature graph. Nonzero terminal residuals and signed omitted pieces are
retained. A tiny normalized balance is not a global energy certificate.
"""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_joined_outer import build_joined_field, _mp


def terminal_balance(field,Z,*,nu='.01',tau=1):
    """Return terminal balance terms and physical logarithmic energy density.

    For R beyond axial support: ur = -nu*C(Z)/r, where
    C=((1-delta)*Z*Mz+(1-Z^2)*Mz_Z)/(1-delta*Z^2).
    The kinetic energy density per dZ d(log r) is
    pi*nu^2*C^2*dz/dZ. A nonzero continuous C gives infinite radial energy.
    This function does not set a residual or its derivative to zero.
    """
    with mp.workdps(field.precision):
        z=_mp(Z);viscosity=_mp(nu);time_gap=_mp(tau)
        if abs(z)>=1 or viscosity<=0 or time_gap<=0:
            raise ValueError('Require |Z|<1, nu>0 and tau>0')
        profile=field.outer;logRv=_mp(str(field.schedule.logR_v))
        text=mp.nstr(logRv,field.precision)
        jet=profile.axial_average_jet(text,z)
        component=profile.component(float(z))
        tangent=profile.coefficient_tangent(float(z))
        lam=mp.mpf('.5')-component.mu
        pulse=component.pulse.full_row(component.mu,1)
        pulse_integral=mp.exp(_mp(pulse['log_normalized_pulse_integral']))
        bump_atoms=[mp.exp(lam*center)*component.basis.full(mp.nstr(lam,component.precision))
                    for center in (-3,-1)]
        terms=[component.base[0],component.a*pulse_integral]+[
            c*atom for c,atom in zip(component.c,bump_atoms)]
        baseZ=from_signed_log(tangent['input']['base_Z'][0])
        termsZ=[baseZ,tangent['a_Z']*pulse_integral]+[
            c*atom for c,atom in zip(tangent['c_Z'],bump_atoms)]
        scale=sum(abs(v) for v in terms);scaleZ=sum(abs(v) for v in termsZ)
        mass_over_Rh=jet['value']*mp.exp(logRv-field.logRh)
        massZ_over_Rh=jet['derivative']*mp.exp(logRv-field.logRh)
        delta=field.delta;d=1-z*z;L=1-delta*z*z
        C_over_Rh=((1-delta)*z*mass_over_Rh+d*massZ_over_Rh)/L
        dzdZ=mp.sqrt(viscosity)*time_gap**((1-delta)/2)*L/d**((3-delta)/2)
        log_energy_density=(mp.log(mp.pi*viscosity**2*dzdZ)+
                            2*(field.logRh+mp.log(abs(C_over_Rh)))) if C_over_Rh else None
        bound=jet.get('pulse_omitted_absolute_bound')
        boundZ=jet.get('pulse_omitted_derivative_absolute_bound')
        conditional_C_bound_over_Rh=None
        if bound is not None and boundZ is not None:
            conditional_C_bound_over_Rh=(abs((1-delta)*z)*from_signed_log(bound)+
                d*from_signed_log(boundZ))*mp.exp(logRv-field.logRh)/L
        return dict(Z=mp.nstr(z,40),nu=mp.nstr(viscosity,40),tau=mp.nstr(time_gap,40),
            normalized_mass_terms=[signed_log(v,field.precision) for v in terms],
            normalized_mass_Z_terms=[signed_log(v,field.precision) for v in termsZ],
            normalized_mass_relative_balance=mp.nstr(abs(sum(terms))/scale,60),
            normalized_mass_Z_relative_balance=mp.nstr(abs(sum(termsZ))/scaleZ,60) if scaleZ else None,
            installed_normalized_mass=signed_log(jet['N'],field.precision),
            installed_normalized_mass_Z=signed_log(jet['N_Z'],field.precision),
            installed_mass_relative_balance=mp.nstr(abs(jet['N'])/scale,60),
            installed_mass_Z_relative_balance=mp.nstr(abs(jet['N_Z'])/scaleZ,60) if scaleZ else None,
            primitive_vs_term_replay_difference=signed_log(jet['N']-sum(terms),field.precision),
            primitive_Z_vs_term_replay_difference=signed_log(jet['N_Z']-sum(termsZ),field.precision),
            actual_mass_over_Rh=signed_log(mass_over_Rh,field.precision),
            actual_mass_Z_over_Rh=signed_log(massZ_over_Rh,field.precision),
            radial_tail_coefficient_over_Rh=signed_log(C_over_Rh,field.precision),
            log_energy_per_dZ_dlog_r=mp.nstr(log_energy_density,field.precision) if log_energy_density is not None else None,
            energy_formula='ur=-nu*C(Z)/r; dE/(dZ dlogr)=pi*nu^2*C(Z)^2*dz/dZ',
            pulse_omitted_bound_for_C_over_Rh=(signed_log(conditional_C_bound_over_Rh,field.precision)
                if conditional_C_bound_over_Rh is not None else None),
            omitted_bound_scope='Pulse kernel omitted pieces only, conditional on nominal coefficients; excludes integral quadrature, incoming inputs, and coefficient errors',
            nominal_nonzero_tail=bool(C_over_Rh),terminal_mean_forced_zero=False,
            global_finite_energy_certified=False,scale_recursion_certified=False)


def run():
    folder=Path(__file__).parent
    print('building installed continuous candidate',flush=True)
    source=build_joined_field()
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    print('evaluating actual terminal mass and Z mass',flush=True)
    result=terminal_balance(field,'.3')
    result['source_commit']='9a929d36'
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('normalized_mass_relative_balance',
        'normalized_mass_Z_relative_balance','nominal_nonzero_tail')}),flush=True)
    return result


if __name__=='__main__':run()
