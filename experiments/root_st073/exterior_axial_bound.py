"""Screen the axial momentum residual throughout the pure heat exterior.

The result is an analytic spatial supremum with a numerically integrated
radial pressure at the inner edge. It is not a certified interval bound.
"""
import json
import numpy as np

from radial_continuation import ROOT
from heat_exterior import physical
from heated_global_join import HeatedGlobalJoin


def axial_spatial_supremum(field, tau):
    """sup_{r>=r2,z in R} |2 G G_z p_K(r,tau)| at fixed tau.

    p_K=-integral_r^infinity K^2/rho drho, so |p_K| decreases in r.
    The Gaussian factor attains its maximum at |z|=sqrt(nu*age).
    """
    _, r2 = field.ramp_radii(tau)
    p_edge = physical(np.array([[r2, 0., 0.]]), tau,
                      c=field.heat.amplitude, h=field.heat.h,
                      nu=field.nu)['pressure'][0]
    age = field.heat.heat_age + .5 - tau
    supremum = (abs(p_edge) * field.heat.heat_age / age
                 * np.exp(-.5) / np.sqrt(field.nu * age))
    return float(supremum), float(p_edge), float(r2)


def run():
    tau_min, tau_max = .5/64, .5
    rows = []
    for heat_age in (1., 512., 1024., 2048.):
        field = HeatedGlobalJoin(heat_age=heat_age)
        bound, p_edge, r2 = axial_spatial_supremum(field, tau_min)
        age = heat_age + .5 - tau_min
        peak = np.array([[r2, 0., np.sqrt(field.nu*age)]])
        direct_peak = abs(field.heat.analytic_residual(peak, tau_min)[0, 2])
        if not np.isclose(bound, direct_peak, rtol=1e-12, atol=1e-14):
            raise AssertionError('Closed-form spatial supremum disagrees with peak evaluation')
        rows.append(dict(heat_age=heat_age,
                         worst_time_tau=tau_min,
                         exterior_axial_spatial_supremum=bound,
                         direct_peak_evaluation=float(direct_peak),
                         pressure_at_r2=p_edge, r2=r2,
                         axial_energy_factor_at_tau_min=field.heat.axial_energy_factor(tau_min),
                         below_1e_3=bound < 1e-3))
    # r2(tau) is proportional to sqrt(tau), while Z=2tau/s(r2) is
    # constant, so |p_K(r2,tau)| is proportional to tau^(-1-2h).
    # For heat_age>=1 and tau in [1/128,1/2], the log derivative of the
    # bound is -(1+2h)/tau + 3/(2*age) < 0; tau_min is worst.
    report = dict(time_interval=[tau_min, tau_max], rows=rows,
                  spatial_scope='Only the pure exterior r>=r2(tau), all z; radial ramp and inner field excluded.',
                  derivation='Pure swirl has zero radial and angular momentum residual analytically. Axial residual is 2*G*G_z*p_K. Pressure magnitude decreases with radius; Gaussian derivative peaks at |z|=sqrt(nu*age). At r2 proportional to sqrt(tau), |p_K| proportional to tau^(-1-2h), and the bound decreases with tau for heat_age>=1.',
                  pressure_evaluation='Positive quadrature numerical estimate, not an interval-certified bound.',
                  global_pde_validated=False)
    (ROOT/'exterior_axial_bound.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
