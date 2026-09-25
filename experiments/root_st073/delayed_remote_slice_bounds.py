"""Analytic per-slice lower bounds for staged outgoing-moment repair.

At fixed eta, the four radial swirl modes alter only E and the four
radial poloidal modes alter only U. I and J are linear constraints;
Cp and S are quadratic. Minimize each quadratic under its paired
linear equality to detect a genuine radial-shape reachability limit.
The bounds ignore axial coupling and full PDE cost, so they are only
necessary conditions for the forty-mode construction.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import (NEAR5_AXIAL_INTERVALS,
                                          NEAR_RADIAL_INTERVALS, TRAIN_ETAS,
                                          current_mean, defects, make_rows,
                                          quadrature, remote_modes)
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT


def equality_quadratic_minimum(A, B, weights, linear_weight, target_linear):
    """Minimize <weights,(A+B c)^2> with <linear_weight,A+B c> fixed."""
    hessian = B.T@(weights[:, None]*B)
    gradient = B.T@(weights*A)
    constraint = linear_weight@B
    right = target_linear-linear_weight@A
    kkt = np.block([[hessian, constraint[:, None]],
                    [constraint[None, :], np.zeros((1, 1))]])
    solution = np.linalg.solve(kkt, np.r_[-gradient, right])[:B.shape[1]]
    field = A+B@solution
    return solution, field, float(weights@(field**2)), float(
        np.linalg.cond(kkt))


def run():
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16, 2.)
    modes, _ = remote_modes(base, NEAR_RADIAL_INTERVALS,
                            NEAR5_AXIAL_INTERVALS)
    X, weights = quadrature(radial_intervals=NEAR_RADIAL_INTERVALS)
    tau = .5*2**(-5.5)
    rows = make_rows(mean, target, base, modes, X, weights,
                     TRAIN_ETAS, tau)
    prior = json.loads((ROOT/'delayed_remote_moment_repair_near5wide.json').read_text())
    c_prior = np.asarray(prior['coefficients'])
    # Each selected axial bump equals one at its corresponding eta.
    center_axial = (0, 4, 1, 3, 2)
    report_rows = []
    for row, axial in zip(rows, center_axial):
        E = row['E']
        E_prior = defects(row, c_prior)[2]
        e_basis = row['Be'][:, 8*axial:8*axial+4]
        u_basis = row['Bu'][:, 8*axial+4:8*axial+8]
        target_I = float(row['target'][1])
        target_Cp = float(row['target'][4])
        target_J = float(row['target'][2])
        target_S = float(row['target'][3])
        e_coeff, e_min, cp_min, e_cond = equality_quadratic_minimum(
            E, e_basis, weights/(2*X),
            weights*np.sqrt(2*X), target_I)
        # U can be solved after E. Use the near-positive previous E and
        # the unconstrained Cp-minimizing E as two diagnostic scenarios.
        u_scenarios = []
        for label, e_field in (('near_positive_fit', E_prior),
                               ('minimum_Cp_E', e_min)):
            u_coeff, u_min, u_norm_min, u_cond = (
                equality_quadratic_minimum(
                    row['U'], u_basis, weights,
                    weights*np.sqrt(2*X)*e_field, target_J))
            s_min = u_norm_min-float(weights@(e_field**2))/2
            # Optimistic continuum bound: allow arbitrary U throughout
            # the union of the radial supports while holding U fixed
            # elsewhere. Preserve M and J, as any compact streamfunction
            # repair must. Endpoint smoothness is deliberately ignored.
            inside = (X >= 1.04) & (X <= 2.5)
            h = np.sqrt(2*X)*e_field
            w_in = weights[inside]
            h_in = h[inside]
            mass_in = (float(row['target'][0])
                       -float(weights[~inside]@row['U'][~inside]))
            mixed_in = target_J-float(weights[~inside]@(
                row['U'][~inside]*h[~inside]))
            gram = np.array([[np.sum(w_in), w_in@h_in],
                             [w_in@h_in, w_in@(h_in**2)]])
            alpha = np.linalg.solve(gram, [mass_in, mixed_in])
            min_u_inside = float(np.array([mass_in, mixed_in])@alpha)
            min_u_outside = float(weights[~inside]@(
                row['U'][~inside]**2))
            continuum_s_min = (min_u_inside+min_u_outside
                               -float(weights@(e_field**2))/2)
            u_scenarios.append(dict(label=label, S_min=s_min,
                                    target_minus_S_min=target_S-s_min,
                                    optimistic_continuum_S_min=(
                                        continuum_s_min),
                                    target_minus_continuum_S_min=(
                                        target_S-continuum_s_min),
                                    KKT_condition=u_cond,
                                    coefficients=u_coeff.tolist()))
        report_rows.append(dict(eta=row['eta'], target_Cp=target_Cp,
                                Cp_min=cp_min,
                                target_minus_Cp_min=target_Cp-cp_min,
                                E_min_relative_to_target=float(np.min(
                                    e_min/row['E0'])),
                                E_KKT_condition=e_cond,
                                E_coefficients=e_coeff.tolist(),
                                U_scenarios=u_scenarios))
    report = dict(source='delayed_remote_moment_repair_near5wide.json',
                  rows=report_rows,
                  scope='Independent fixed-eta quadratic lower bounds '
                        'for I/Cp and J/S. Negative target-minus-minimum '
                        'would rule out that staged radial subproblem. '
                        'Positive margin does not construct a positive '
                        'global field or prove PDE matching.',
                  accepted=False)
    (ROOT/'delayed_remote_slice_bounds.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps([dict(eta=r['eta'], Cp_margin=r['target_minus_Cp_min'],
                           E_min=r['E_min_relative_to_target'],
                           S_margins=[v['target_minus_S_min']
                                      for v in r['U_scenarios']],
                           continuum_S_margins=[
                               v['target_minus_continuum_S_min']
                               for v in r['U_scenarios']])
                      for r in report_rows]), flush=True)


if __name__ == '__main__':
    run()
