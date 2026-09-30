"""Transport actual inner moment offsets into the Section 7.31/7.34 inputs.

Offsets are raw cumulative moments relative to the temporary reference,
retained unchanged through the intervening exterior. This prepares a new
coefficient solve; it does not assert closure or replace the old profile.
"""
from copy import deepcopy
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


def seeded_axial_inputs(incoming, *, log_Rp, mu, offsets, precision=443):
    """Return updated incoming receipt for actual z/theta_z/z_theta offsets.

Use the SAME schedule's Rp and inherited log_Ep. The energy moment uses
the source convention integral(Uz**2-Utheta**2/2), not kinetic energy.
No default zero is supplied for missing offsets.
"""
    result = deepcopy(incoming)
    with mp.workdps(precision):
        rp = mp.mpf(log_Rp); ep = mp.mpf(incoming['log_Ep'])
        mu = mp.mpf(mu)
        delta = {k: mp.mpf(offsets[k]) for k in ('z','theta_z','z_theta')}
        changes = {
            'm1_Mz_over_RpEp': delta['z']*mp.exp(-rp-ep),
            'm2_Mtheta_z_over_Rp_sqrt2RpEp2': delta['theta_z']*mp.exp(-mp.mpf('1.5')*rp-2*ep)/mp.sqrt(2),
            'E_prior_mu_Mztheta_over_RpEp2': mu*delta['z_theta']*mp.exp(-rp-2*ep),
        }
        for key, change in changes.items():
            result[key] = signed_log(from_signed_log(incoming[key])+change, precision)
        rows = result['row_normalization']
        for i, key in enumerate(('m1_Mz_over_RpEp','m2_Mtheta_z_over_Rp_sqrt2RpEp2'),1):
            scale = mp.exp(mp.mpf(rows['log_scale'+str(i)]))
            base = from_signed_log(result[key])*scale
            rows['scaled_base_m'+str(i)] = signed_log(base, precision)
            rows['scaled_base_rhs'+str(i)+'_minus_m'+str(i)] = signed_log(-base, precision)
        result['inner_seed_transport'] = {
            'raw_offsets': {k:signed_log(v,precision) for k,v in delta.items()},
            'normalized_changes': {k:signed_log(v,precision) for k,v in changes.items()},
            'log_Rp':mp.nstr(rp,precision),
            'scope':'Inputs only; downstream pulse/end-bump solve and primitive replay must use this receipt.',
            'exterior_coefficients_recomputed':False,
            'finite_energy_certified':False,
        }
    return result


def solve_seeded_axial(profile, Z, *, offsets, precision=443):
    """Re-solve the actual exterior equations using a retained inner seed.

    This returns coefficients only. The old profile still uses its original
    coefficients and primitive until an explicit adapter installs both the
    new pulse/end bumps and the propagated incoming moment offsets.
    """
    from lei_ren_part1_paper_axial_correction import (
        solve_actual_axial, independent_end_bump_replay)
    schedule = profile.schedule
    data = profile.coefficients(float(Z))
    with mp.workdps(precision):
        incoming = seeded_axial_inputs(data['incoming'],
            log_Rp=str(schedule.logR_p), mu=str(schedule.mu),
            offsets=offsets, precision=precision)
        prior = from_signed_log(incoming['E_prior_mu_Mztheta_over_RpEp2'])
        future = profile.tail.evaluate(Z, quadrature_order=profile.order)
        target = (1-mp.exp(-26))/4-prior+mp.mpf(
            future['energy_target_contribution_nominal'])
        base = [incoming['row_normalization'][key]
            for key in ('scaled_base_m1','scaled_base_m2')]
        solved = solve_actual_axial(str(schedule.mu), base, profile.pulse,
            mp.nstr(target,precision), precision=precision, order=profile.order)
        replay = independent_end_bump_replay(str(schedule.mu), solved,
            base, profile.pulse, order=solved['quadrature_order'],
            precision=precision)
        return dict(incoming=incoming, axial=solved,
            canonical_quadrature_replay=replay,
            energy_target=mp.nstr(target,precision),
            scope='Actual seeded coefficients; not installed in the velocity/primitive provider.',
            continuous_mean_closure_certified=False,
            finite_energy_certified=False)


def run_regression():
    """Reproducible fixture check, explicitly separate from the shared core."""
    import json
    from pathlib import Path
    from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
    profile = CorrectedSourceProfile(precision=160)
    offsets = dict(z='1e-10',theta_z='-2e-10',z_theta='3e-10')
    result = solve_seeded_axial(profile,.5,offsets=offsets,precision=160)
    rows = result['canonical_quadrature_replay']['row_relative_differences']
    with mp.workdps(160):
        assert all(mp.mpf(v)<mp.mpf('1e-30') for v in rows)
        assert mp.mpf(result['axial']['energy_relative_replay'])<mp.mpf('1e-40')
    report = dict(fixture='Historical default exterior; explicit injected offsets, not measured shared-core offsets',
        Z='.5',raw_offsets=offsets,seed_transport=result['incoming']['inner_seed_transport'],
        a_p=result['axial']['a_p'],linear_relative_replay=rows,
        energy_relative_replay=result['axial']['energy_relative_replay'],
        installed_in_velocity_provider=False,continuous_mean_closure_certified=False,
        finite_energy_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('linear_relative_replay','energy_relative_replay','installed_in_velocity_provider')}))
    return report


def solve_joined_seed(field, Z):
    """Use the measured corrected Rh moments of a JoinedOuterField."""
    from lei_ren_part1_paper_axial_restore import reference_moments
    with mp.workdps(field.precision):
        z = mp.mpf(Z)
        terminal = field.inner.evaluate_x(mp.e,z)
        R = terminal['R']
        u = mp.sqrt(2*R)*terminal['F']
        reference = reference_moments(R,u,z)
        offsets = {k:terminal['moments'][k]-reference[k]
            for k in ('z','theta_z','z_theta')}
        result = solve_seeded_axial(field.outer,z,offsets=offsets,
            precision=field.precision)
        result['seed_origin'] = 'Actual five-bump corrected inner terminal at Rh'
        result['seed_log_Rh'] = mp.nstr(terminal['logR'],field.precision)
        result['unresolved_inner_input_entries'] = terminal['unresolved_input_entries']
        result['subtraction_and_inherited_quadrature_uncertainty_certified'] = False
        return result


if __name__=='__main__':
    import sys
    if '--shared' in sys.argv:
        import json
        from pathlib import Path
        from lei_ren_part1_paper_joined_outer import build_joined_field
        result = solve_joined_seed(build_joined_field(),'.3')
        Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').write_text(
            json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({'linear_replay':result['canonical_quadrature_replay']['row_relative_differences'],
            'energy_replay':result['axial']['energy_relative_replay'],
            'installed_in_velocity_provider':False}))
    else:
        run_regression()
