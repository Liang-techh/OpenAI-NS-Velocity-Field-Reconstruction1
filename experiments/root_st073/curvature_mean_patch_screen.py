"""Bounded amplitude-zero mean patch on the curvature-optimized physical lift.

The patch uses twelve compact axisymmetric modes.  Eight modes are generated
from compact cylindrical vector potentials and four are compact pressure
functions.  The vector-potential part has zero divergence analytically; the
screen still checks the complete Cartesian momentum residual with an
independent fourth-order finite-difference evaluator.

This is a local numerical experiment around X = 1.01 and eta in [.2, .3].
It does not evolve a wave amplitude, replace the slice lift, or claim a
moment closure.
"""
import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from joined_field import independent_fd
from joint_collar_fit import kinematics
from radial_continuation import ROOT


TAU_TRAIN = .5 * 2**(-5.5)
TAU_TIME = .5 * 2**(-5.35)
X_BOUNDS = (1.004, 1.065)
ETA_BOUNDS = (.14, .36)

# These sets have no common (X, eta) pair.  The time set is sampled at a
# nearby tau and is also spatially shifted so that it is a separate probe.
TRAIN_X = np.array([1.008, 1.018, 1.028, 1.038])
TRAIN_ETA = np.array([.18, .24, .32])
SPACE_X = np.array([1.012, 1.022, 1.033, 1.045])
SPACE_ETA = np.array([.20, .27, .30])
TIME_X = np.array([1.010, 1.020, 1.031, 1.044])
TIME_ETA = np.array([.195, .255, .315])

# Physical coefficient limits keep the fit local while allowing the pressure
# basis to carry the different units of a scalar pressure correction.
COEFFICIENT_BOUNDS = np.r_[np.full(8, .05), np.full(4, 2000.)]


def _bump(s):
    """Return a C4 compact bump and its derivative with respect to s."""
    s = np.asarray(s, float)
    value = np.zeros_like(s)
    derivative = np.zeros_like(s)
    inside = (s > 0.) & (s < 1.)
    si = s[inside]
    value[inside] = 1024. * si**5 * (1. - si)**5
    derivative[inside] = (5120. * si**4 * (1. - si)**4
                          * (1. - 2. * si))
    return value, derivative


class AxisymmetricMeanPatch:
    """Compact physical basis for two vector potentials and pressure."""

    # Two polynomial degrees in each local coordinate for each of A_theta,
    # A_z, and p: 3 * 2 * 2 = 12 modes.
    names = tuple(f'{family}_{a}{b}'
                  for family in ('A_theta', 'A_z', 'pressure')
                  for a in range(2) for b in range(2))

    def __init__(self, base, tau_values=(TAU_TRAIN, TAU_TIME),
                 x_bounds=X_BOUNDS, eta_bounds=ETA_BOUNDS):
        self.base = base
        self.nu = base.nu
        self.x_bounds = tuple(float(value) for value in x_bounds)
        self.eta_bounds = tuple(float(value) for value in eta_bounds)
        corners = []
        for tau in tau_values:
            X = np.array([self.x_bounds[0], self.x_bounds[1],
                          self.x_bounds[0], self.x_bounds[1]])
            eta = np.array([self.eta_bounds[0], self.eta_bounds[0],
                            self.eta_bounds[1], self.eta_bounds[1]])
            corners.append(base.compact.joined.inner.from_similarity(
                X, eta, tau))
        corners = np.concatenate(corners, axis=0)
        self.r_bounds = (float(np.min(corners[:, 0])),
                         float(np.max(corners[:, 0])))
        self.z_bounds = (float(np.min(corners[:, 2])),
                         float(np.max(corners[:, 2])))
        self.dr = self.r_bounds[1] - self.r_bounds[0]
        self.dz = self.z_bounds[1] - self.z_bounds[0]
        if self.dr <= 0. or self.dz <= 0.:
            raise ValueError('Degenerate physical patch support')

    def _quantities(self, points, a, b):
        """Return q, q_r, q_z for one compact scalar basis function."""
        pts = np.asarray(points, float)
        radius = np.hypot(pts[:, 0], pts[:, 1])
        r0, r1 = self.r_bounds
        z0, z1 = self.z_bounds
        xi = (radius-r0) / self.dr
        eta = (pts[:, 2]-z0) / self.dz
        br, br_s = _bump(xi)
        bz, bz_s = _bump(eta)
        radial_poly = xi**a
        axial_poly = eta**b
        q = br*bz*radial_poly*axial_poly
        q_r = (br_s/self.dr*bz*radial_poly*axial_poly
               + (a/self.dr*xi**(a-1) if a else 0.)
               * br*bz*axial_poly)
        q_z = (br*bz_s/self.dz*radial_poly*axial_poly
               + (b/self.dz*eta**(b-1) if b else 0.)
               * br*bz*radial_poly)
        return q, q_r, q_z

    def mode_fields(self, points, tau, index):
        """Evaluate one correction mode, with no explicit time dependence."""
        del tau  # amplitude-zero screen: the physical patch is time independent
        pts = np.asarray(points, float)
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        family_index = index // 4
        a = (index % 4) // 2
        b = index % 2
        family = ('A_theta', 'A_z', 'pressure')[family_index]
        q, q_r, q_z = self._quantities(pts, a, b)
        radius = np.hypot(pts[:, 0], pts[:, 1])
        safe_radius = np.where(radius > 0., radius, 1.)
        ca = np.divide(pts[:, 0], safe_radius)
        sa = np.divide(pts[:, 1], safe_radius)
        if family == 'A_theta':
            # curl(A_theta e_theta)_r = -q_z and
            # curl(A_theta e_theta)_z = q_r + q/r.
            ur = -q_z
            uz = q_r + q/safe_radius
            velocity[:, 0] = ca*ur
            velocity[:, 1] = sa*ur
            velocity[:, 2] = uz
        elif family == 'A_z':
            # curl(q e_z)_theta = -q_r.
            uth = -q_r
            velocity[:, 0] = -sa*uth
            velocity[:, 1] = ca*uth
        else:
            pressure = q
        return velocity, pressure

    def fields(self, points, tau, coefficients=None):
        """Evaluate the correction represented by ``coefficients``."""
        pts = np.asarray(points, float)
        ts = np.broadcast_to(np.asarray(tau, float), (len(pts),))
        if coefficients is None:
            coefficients = np.ones(len(self.names))
        coefficients = np.asarray(coefficients, float)
        velocity = np.zeros_like(pts)
        pressure = np.zeros(len(pts))
        for index, coefficient in enumerate(coefficients):
            if coefficient == 0.:
                continue
            mode_velocity, mode_pressure = self.mode_fields(pts, ts, index)
            velocity += coefficient*mode_velocity
            pressure += coefficient*mode_pressure
        return velocity, pressure


class PatchMode:
    """Field adapter used to build the linearized momentum columns."""

    def __init__(self, patch, index):
        self.patch = patch
        self.index = index
        self.nu = patch.nu

    def fields(self, points, tau):
        return self.patch.mode_fields(points, tau, self.index)


class PatchedLift:
    """Current physical lift plus a fitted compact correction."""

    def __init__(self, base, patch, coefficients):
        self.base = base
        self.patch = patch
        self.coefficients = np.asarray(coefficients, float)
        self.nu = base.nu

    def fields(self, points, tau):
        velocity, pressure = self.base.fields(points, tau)
        dv, dp = self.patch.fields(points, tau, self.coefficients)
        return velocity + dv, pressure + dp


def _points(field, X_values, eta_values, tau):
    X, eta = np.meshgrid(X_values, eta_values, indexing='ij')
    points = field.compact.joined.inner.from_similarity(
        X.ravel(), eta.ravel(), tau)
    return points, X.ravel(), eta.ravel()


def _kinematic_data(field, points, tau, hs, ht):
    return kinematics(field, points, tau, hs, ht)


def _assemble(base_data, mode_data, coefficients):
    """Return nonlinear momentum residual and its state at coefficients."""
    u0, grad0, part0 = base_data
    mode_u = np.stack([data[0] for data in mode_data], axis=-1)
    mode_grad = np.stack([data[1] for data in mode_data], axis=-1)
    mode_part = np.stack([data[2] for data in mode_data], axis=-1)
    coefficients = np.asarray(coefficients, float)
    u = u0 + np.einsum('nim,m->ni', mode_u, coefficients)
    grad = grad0 + np.einsum('nijm,m->nij', mode_grad, coefficients)
    part = part0 + np.einsum('nim,m->ni', mode_part, coefficients)
    residual = part + np.einsum('nij,nj->ni', grad, u)
    return u, grad, part, residual, mode_u, mode_grad, mode_part


def _fit_coefficients(base, patch, points, tau, hs, ht, iterations=3):
    base_data = _kinematic_data(base, points, tau, hs, ht)
    mode_data = [_kinematic_data(PatchMode(patch, index), points, tau, hs, ht)
                 for index in range(len(patch.names))]
    coefficients = np.zeros(len(patch.names))
    history = []
    for iteration in range(iterations):
        u, grad, part, residual, mode_u, mode_grad, mode_part = _assemble(
            base_data, mode_data, coefficients)
        columns = np.empty((len(points), 3, len(patch.names)))
        for index in range(len(patch.names)):
            columns[:, :, index] = (
                mode_part[:, :, index]
                + np.einsum('nij,nj->ni', mode_grad[:, :, :, index], u)
                + np.einsum('nij,nj->ni', grad, mode_u[:, :, index]))
        flat = columns.reshape(-1, len(patch.names))
        scale = np.maximum(np.linalg.norm(flat, axis=0), 1e-30)
        normalized = flat/scale
        ridge = .02 * np.linalg.norm(normalized, ord=2)
        augmented = np.vstack((normalized, ridge*np.eye(len(patch.names))))
        target = np.r_[-residual.reshape(-1),
                       np.zeros(len(patch.names))]
        step_normalized = np.linalg.lstsq(augmented, target, rcond=None)[0]
        step = step_normalized/scale
        current_score = float(np.sqrt(np.mean(
            np.linalg.norm(residual, axis=1)**2)))
        # The raw Newton step is intentionally damped.  The compact basis has
        # large third derivatives, so accepting it at unit strength would
        # improve a linearized column while worsening the nonlinear residual.
        trial_strengths = np.r_[0., 2.**np.arange(0., -15., -1.)]
        trials = []
        for strength in trial_strengths:
            trial = np.clip(coefficients + strength*step,
                            -COEFFICIENT_BOUNDS, COEFFICIENT_BOUNDS)
            trial_residual = _assemble(base_data, mode_data, trial)[3]
            trial_score = float(np.sqrt(np.mean(
                np.linalg.norm(trial_residual, axis=1)**2)))
            trials.append((trial_score, trial, strength))
        trial_score, trial_coefficients, strength = min(
            trials, key=lambda item: item[0])
        if trial_score > current_score:
            strength = 0.
            trial_coefficients = coefficients.copy()
            trial_score = current_score
        coefficients = trial_coefficients
        history.append(dict(iteration=iteration,
                            residual_max_before=float(np.max(
                                np.linalg.norm(residual, axis=1))),
                            step_norm=float(np.linalg.norm(step)),
                            normalized_step_max=float(np.max(
                                np.abs(step_normalized))),
                            accepted_strength=float(strength),
                            residual_rms_after=trial_score,
                            coefficient_norm=float(np.linalg.norm(
                                coefficients))))
    _, _, _, final_residual, _, _, _ = _assemble(
        base_data, mode_data, coefficients)
    history[-1]['residual_max_after'] = float(np.max(
        np.linalg.norm(final_residual, axis=1)))
    return coefficients, history


def _metrics(label, residual, divergence, X, eta, tau):
    norms = np.linalg.norm(residual, axis=1)
    return dict(
        field=label,
        tau=float(tau),
        node_count=int(len(residual)),
        max_full_momentum=float(np.max(norms)),
        rms_full_momentum=float(np.sqrt(np.mean(norms**2))),
        max_component_abs=np.max(np.abs(residual), axis=0).tolist(),
        max_divergence=float(np.max(np.abs(divergence))),
        nodes=[dict(X=float(x), eta=float(e),
                    residual_norm=float(norm),
                    residual=list(map(float, vector)),
                    divergence=float(div))
               for x, e, norm, vector, div in
               zip(X, eta, norms, residual, divergence)])


def run():
    base = CoupledMomentPhysicalLift(
        slice_filename='wide_taper_curvature_optimize.json')
    patch = AxisymmetricMeanPatch(base)
    hs = .001*np.sqrt(base.nu*TAU_TRAIN)
    ht = .00025*TAU_TRAIN

    train_points, train_X, train_eta = _points(
        base, TRAIN_X, TRAIN_ETA, TAU_TRAIN)
    space_points, space_X, space_eta = _points(
        base, SPACE_X, SPACE_ETA, TAU_TRAIN)
    time_points, time_X, time_eta = _points(
        base, TIME_X, TIME_ETA, TAU_TIME)

    coefficients, fit_history = _fit_coefficients(
        base, patch, train_points, TAU_TRAIN, hs, ht)
    corrected = PatchedLift(base, patch, coefficients)

    all_points = np.concatenate((train_points, space_points, time_points))
    all_tau = np.concatenate((np.full(len(train_points), TAU_TRAIN),
                              np.full(len(space_points), TAU_TRAIN),
                              np.full(len(time_points), TAU_TIME)))
    all_X = np.concatenate((train_X, space_X, time_X))
    all_eta = np.concatenate((train_eta, space_eta, time_eta))
    split_sizes = (len(train_points), len(space_points), len(time_points))
    split_offsets = np.r_[0, np.cumsum(split_sizes)]

    # independent_fd is the validation path.  A single batched call per field
    # keeps the finite-difference stencil identical across the three splits.
    base_residual, base_divergence = independent_fd(
        base, all_points, all_tau, hs, ht)
    corrected_residual, corrected_divergence = independent_fd(
        corrected, all_points, all_tau, hs, ht)

    rows = []
    for name, lo, hi in zip(('train', 'space_holdout', 'time_holdout'),
                            split_offsets[:-1], split_offsets[1:]):
        rows.append(dict(name=name,
                         before=_metrics('base', base_residual[lo:hi],
                                         base_divergence[lo:hi],
                                         all_X[lo:hi], all_eta[lo:hi],
                                         all_tau[lo:hi][0]),
                         after=_metrics('patched', corrected_residual[lo:hi],
                                        corrected_divergence[lo:hi],
                                        all_X[lo:hi], all_eta[lo:hi],
                                        all_tau[lo:hi][0])))
        print(json.dumps({
            'name': name,
            'before_max': rows[-1]['before']['max_full_momentum'],
            'after_max': rows[-1]['after']['max_full_momentum'],
            'before_rms': rows[-1]['before']['rms_full_momentum'],
            'after_rms': rows[-1]['after']['rms_full_momentum'],
            'after_divergence': rows[-1]['after']['max_divergence']}),
              flush=True)

    # Difference of the two divergence traces is the direct FD divergence of
    # the added exact-curl velocity, since divergence is linear.
    patch_fd_divergence = corrected_divergence-base_divergence
    report = dict(
        slice_filename='wide_taper_curvature_optimize.json',
        tau_train=TAU_TRAIN,
        tau_time_holdout=TAU_TIME,
        tau_time_ratio=TAU_TIME/TAU_TRAIN,
        mode_names=list(patch.names),
        mode_count=len(patch.names),
        x_similarity_bounds=list(patch.x_bounds),
        eta_similarity_bounds=list(patch.eta_bounds),
        physical_r_bounds=list(patch.r_bounds),
        physical_z_bounds=list(patch.z_bounds),
        spatial_step=hs,
        time_step=ht,
        coefficients=coefficients.tolist(),
        coefficient_norm=float(np.linalg.norm(coefficients)),
        coefficient_bounds=COEFFICIENT_BOUNDS.tolist(),
        fit_history=fit_history,
        rows=rows,
        exact_curl_analytic_divergence=0.0,
        added_patch_fd_divergence_max=float(np.max(np.abs(
            patch_fd_divergence))),
        scope=('Twelve compact amplitude-zero axisymmetric mean modes: '
               'eight vector-potential velocity modes plus four compact '
               'pressure modes. Full physical Cartesian momentum is fitted '
               'on twelve training nodes and checked by independent fourth-'
               'order finite differences on disjoint spatial and nearby-time '
               'holdouts. Local screen only; no dynamic wave source, '
               'five-moment closure, or global admission.'),
        accepted=False)
    out = ROOT/'curvature_mean_patch_screen.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'output': str(out),
                      'coefficient_norm': report['coefficient_norm']}),
          flush=True)


if __name__ == '__main__':
    run()
