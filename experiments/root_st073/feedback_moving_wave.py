"""Moving-normal exact-curl wave trial on the saved feedback trajectory.

This is a local physical-coordinate principal inverse, not the paper's
normalized supported-pulse construction. Full nonlinear residual decides
whether its centerline cancellation survives spatial reconstruction.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from affine_momentum import jets, momentum
from outer_feedback_evolution import build_current, Trajectory, P_BREAKS
from grouped_outer_cache import outer_cache
from transported_phase_screen import coefficients
from curl_wave_prototype import bump, cylindrical_residual
from midplane_physical_covariance_pairs import support
from radial_continuation import ROOT


class PhasePath:
    def __init__(self, field, mode, initial, interval):
        self.field, self.mode = field, mode
        self.start, self.end = interval
        self.duration = self.end-self.start
        result = solve_ivp(lambda s, x: self.duration*self.rhs(
            self.start+s*self.duration, x), (0., 1.), initial,
            dense_output=True, rtol=1e-9, atol=1e-11, max_step=.125)
        if not result.success:
            raise RuntimeError(result.message)
        self.solution = result
        self._local_cache = {}

    def rhs(self, t, state):
        r, z, kr, kz, phase = state
        ur, uz, F, urr, urz, uzr, uzz, Fr, Fz = coefficients(self.field, r, z, -t)
        return np.array([ur, uz, -urr*kr-uzr*kz-self.mode*Fr,
                         -urz*kr-uzz*kz-self.mode*Fz, -self.mode*F])

    def state(self, t):
        if not self.start <= t <= self.end:
            raise ValueError('Phase path outside interval')
        return self.solution.sol((t-self.start)/self.duration)

    def normal(self, t):
        r, z, kr, kz, phase = self.state(t)
        return np.array([kr, self.mode/r, kz])

    def normal_dot(self, t):
        state = self.state(t)
        ur, uz, F, urr, urz, uzr, uzz, Fr, Fz = self.local_coefficients(t)
        r, z, kr, kz, phase = state
        return np.array([-urr*kr-uzr*kz-self.mode*Fr, -self.mode*ur/r**2,
                         -urz*kr-uzz*kz-self.mode*Fz])

    def local_coefficients(self, t):
        # Only cache along the completed path, never provisional RK stages.
        key = float(t)
        if key not in self._local_cache:
            r, z = self.state(t)[:2]
            self._local_cache[key] = coefficients(self.field, r, z, -t)
        return self._local_cache[key]

    def matrix(self, t):
        r, z = self.state(t)[:2]
        ur, uz, F, urr, urz, uzr, uzz, Fr, Fz = self.local_coefficients(t)
        # Cartesian gradient in cylindrical components plus rotation of basis.
        return np.array([[urr, -2*F, urz], [2*F+r*Fr, ur/r, r*Fz], [uzr, 0., uzz]])


class MovingCurlField:
    def __init__(self, base, paths, amplitudes, pressures, widths):
        self.base, self.paths = base, paths
        self.amplitudes, self.pressures = amplitudes, pressures
        self.dr, self.dz = widths
        self.nu = base.nu

    def fields(self, points, tau):
        points = np.asarray(points)
        times = np.broadcast_to(tau, (len(points),))
        if not np.all(times == times[0]):
            raise ValueError('One time per spatial batch is required')
        t = -float(times[0])
        u, p = self.base.fields(points, tau)
        u, p = u.copy(), p.copy()
        data = []
        for path, amplitude, pressure in zip(self.paths, self.amplitudes, self.pressures):
            state = path.state(t)
            n = path.normal(t)
            potential = 1j*np.cross(n, amplitude(t))/np.dot(n, n)
            data.append((path.mode, state, potential, pressure(t)))
        for i, (x, y, z) in enumerate(points):
            r = np.hypot(x, y)
            if r == 0:
                continue
            theta = np.arctan2(y, x)
            delta = np.zeros(3)
            for m, (rc, zc, kr, kz, phase), c, pressure in data:
                br, br_r = bump(r, rc, self.dr)
                bz, bz_z = bump(z, zc, self.dz)
                if br == 0 or bz == 0:
                    continue
                E, Er, Ez = br*bz, br_r*bz, br*bz_z
                ar, at, az = c
                curl = np.array([1j*m*az*E/r-at*(Ez+1j*kz*E),
                                 ar*(Ez+1j*kz*E)-az*(Er+1j*kr*E),
                                 at*(Er+E/r+1j*kr*E)-1j*m*ar*E/r])
                carrier = np.exp(1j*(m*theta+kr*(r-rc)+kz*(z-zc)+phase))
                delta += (carrier*curl).real
                p[i] += float((carrier*E*pressure).real)
            ca, sa = x/r, y/r
            u[i] += [ca*delta[0]-sa*delta[1], sa*delta[0]+ca*delta[1], delta[2]]
        return u, p


def residual_rows(field, paths, t, offsets, angles, widths):
    dr, dz = widths
    rc, zc = paths[0].state(t)[:2]
    pts = np.array([[(rc+xi*dr)*np.cos(a), (rc+xi*dr)*np.sin(a), zc+eta*dz]
                    for xi, eta in offsets for a in angles])
    tau = -t
    jet = jets(field, pts, tau, .0005*np.sqrt(field.nu*tau), .000005*tau)
    residual = cylindrical_residual(momentum(jet), pts).reshape(len(offsets), len(angles), 3)
    return residual, float(max(abs(np.trace(jet[1], axis1=1, axis2=2))))


def stats(residual):
    norm = np.linalg.norm(residual, axis=-1)
    angular_mean = residual.mean(axis=1)
    oscillatory = residual-angular_mean[:, None, :]
    return dict(max=float(norm.max()), sample_rms=float(np.sqrt(np.mean(norm**2))),
                angular_mean_max=float(np.linalg.norm(angular_mean, axis=-1).max()),
                oscillatory_max=float(np.linalg.norm(oscillatory, axis=-1).max()))


def run():
    # Imported here so the field reconstruction remains independently usable.
    from moving_normal_inverse import solve_path
    inner, _, current = build_current()
    trajectory = json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    mean = Trajectory(current, trajectory['nodes'])
    source = json.loads((ROOT/'compact_potential/midplane_wave_source_k11.json').read_text(encoding='utf-8'))
    initial_k = 11.0003
    interval = (-.5*2.**-initial_k, -.5*2.**-11.0009)
    source['tau'] = -interval[0]
    source['point'] = inner.from_similarity([inner.p.X_max*(1+15*.325)**2],
                                           [-.0125], source['tau'])[0].tolist()
    radius, _, height = source['point']
    args = support(inner, source, mean.nu)
    widths = (args['radial_halfwidth'], 1.4*args['axial_halfwidth'])
    paths, initial = [], []
    for index in (0, 16):
        pulse = source['kelvin_pulses'][index]
        kr, _, kz = pulse['peak_wavevector']
        path = PhasePath(mean, pulse['mode'], [radius, height, kr, kz, 0.], interval)
        n = path.normal(interval[0])
        amplitude = np.array(pulse['peak_amplitude'], complex)
        amplitude -= n*np.dot(n, amplitude)/np.dot(n, n)
        paths.append(path); initial.append(amplitude)
    print(json.dumps(dict(stage='phase_paths_ready')), flush=True)
    # Rebuild the target from the current background with split quadrature.
    cache = outer_cache(mean, [], k=initial_k, order=64,
                        locations=[(-.0125, .325)], radial_breaks=P_BREAKS)
    sl, r, w, R = cache['panels'][0]
    residual = momentum(cache['baseline'])
    target = np.array([-np.dot(w*r*r, residual[sl, 1])/R**2,
                       -np.dot(w*r, residual[sl, 2])/R])
    # Modes are trial carrier directions only; covariance is recomputed here.
    angles = np.arange(24)*2*np.pi/24
    t0 = interval[0]
    points = np.array([[radius*np.cos(a), radius*np.sin(a), height] for a in angles])
    columns = []
    zero = lambda t: 0j
    for j in range(2):
        trial = MovingCurlField(mean, paths,
            [lambda t, i=i, j=j: initial[i] if i == j else np.zeros(3, complex) for i in range(2)],
            [zero, zero], widths)
        delta = trial.fields(points, -t0)[0]-mean.fields(points, -t0)[0]
        cyl = cylindrical_residual(delta, points)
        columns.append(np.mean(cyl[:, 0, None]*cyl[:, 1:], axis=0))
    weights = np.linalg.solve(np.column_stack(columns), target)
    if min(weights) <= 0:
        raise RuntimeError(f'Current physical covariance weights are not positive: {weights}')
    print(json.dumps(dict(stage='current_covariance', weights=weights.tolist(), target=target.tolist())), flush=True)
    initial = [a*np.sqrt(weight) for a, weight in zip(initial, weights)]
    homogeneous = [solve_path(interval, path.normal, path.normal_dot, path.matrix,
                             lambda t: np.zeros(3, complex), mean.nu, initial_amplitude=a)
                   for path, a in zip(paths, initial)]
    moving = MovingCurlField(mean, paths, [v.amplitude for v in homogeneous],
                            [v.pressure for v in homogeneous], widths)
    print(json.dumps(dict(stage='homogeneous_paths_ready')), flush=True)
    sample_k = np.linspace(11.00036, 11.00084, 9)
    sample_t = -.5*2.**-sample_k
    sources = [[], []]
    for index, t in enumerate(sample_t):
        rows, _ = residual_rows(moving, paths, t, [(0., 0.)], angles, widths)
        for j, path in enumerate(paths):
            phase = path.state(t)[4]
            sources[j].append(2*np.mean(rows[0]*np.exp(-1j*(path.mode*angles+phase))[:, None], axis=0))
        if index % 3 == 2:
            print(json.dumps(dict(stage='harmonic_source', completed_times=index+1)), flush=True)
    corrections = []
    for path, values in zip(paths, sources):
        spline = CubicSpline(sample_t, np.array(values), extrapolate=False)
        corrections.append(solve_path((sample_t[0], sample_t[-1]), path.normal,
            path.normal_dot, path.matrix, spline, mean.nu))
    corrected = MovingCurlField(mean, paths,
        [lambda t, i=i: homogeneous[i].amplitude(t)+corrections[i].amplitude(t) for i in range(2)],
        [lambda t, i=i: homogeneous[i].pressure(t)+corrections[i].pressure(t) for i in range(2)], widths)
    report = dict(accepted=False, scale_recursion_established=False, pressure_time_convention='physical t=-tau',
        initial_k=initial_k, source_indices=[0, 16], covariance_target=target.tolist(), weights=weights.tolist(),
        source_sample_k=sample_k.tolist(), widths=list(widths), replay=[],
        scope='Moving center/normal principal harmonic inverse with full compact spatial curl and nonlinear momentum replay. No temporal cutoff, support-wide source inversion, mean repair, volume L2, finite-energy or recursive acceptance.')
    report['initial_carrier_geometry'] = [dict(mode=path.mode,
        radial_phase_halfwidth=float(abs(path.normal(t0)[0])*widths[0]),
        axial_phase_halfwidth=float(abs(path.normal(t0)[2])*widths[1]),
        viscous_damping_over_interval=float(mean.nu*np.dot(path.normal(t0), path.normal(t0))*(interval[1]-interval[0])))
        for path in paths]
    output = ROOT/'feedback_moving_wave.json'
    for k in (11.00045, 11.0006, 11.00075):
        t = -.5*2.**-k
        offsets = [(0., 0.)]+[(x, z) for x in (-.45, .45) for z in (-.45, .45)]
        cases = {}
        for name, field in [('mean', mean), ('moving_homogeneous', moving), ('forced_correction', corrected)]:
            residual, divergence = residual_rows(field, paths, t, offsets, angles, widths)
            cases[name] = dict(center=stats(residual[:1]), spatial_holdout=stats(residual[1:]),
                               sampled_divergence=divergence)
        row = dict(k=k, cases=cases)
        report['replay'].append(row)
        output.write_bytes((json.dumps(report, indent=2)+'\n').encode())
        print(json.dumps(row), flush=True)
    report['holdout_worsening_factors'] = [
        row['cases']['forced_correction']['spatial_holdout']['max'] /
        row['cases']['moving_homogeneous']['spatial_holdout']['max']
        for row in report['replay']]
    report['rejected'] = any(value >= 1 for value in report['holdout_worsening_factors'])
    if report['rejected']:
        report['rejection_reason'] = ('The forced center-path principal correction fails to reduce '
            'full nonlinear momentum maximum at spatial holdouts. Neither field is accepted.')
    output.write_bytes((json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
