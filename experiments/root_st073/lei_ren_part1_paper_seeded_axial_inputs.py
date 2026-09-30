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
