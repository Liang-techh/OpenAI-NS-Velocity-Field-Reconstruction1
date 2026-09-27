"""Full spatial Fourier-potential evolution on the current mean trajectory.

Velocity is an exact cylindrical curl. At each state the complete nonlinear
momentum, including cutoff diffusion, is projected jointly onto potential
time derivatives and pressure for mean and harmonics 1..8. Independent
spatial and temporal replay is required; no principal-wave approximation.
"""
import json
from types import SimpleNamespace
import numpy as np
from scipy.interpolate import CubicHermiteSpline, CubicSpline
from numpy.polynomial.legendre import leggauss
from supported_fourier_basis import basis_data
from affine_momentum import jets, momentum
from outer_feedback_evolution import build_current, Trajectory
from curl_wave_prototype import cylindrical_residual
from radial_continuation import ROOT


def basis_jets(points, center, widths, mode, degree, carrier, nu, h):
    velocity, pressure, gp = basis_data(points, center, widths, mode, degree, carrier)
    gradient = np.empty((len(points), 3, 3, velocity.shape[-1]), complex)
    laplacian = np.zeros_like(velocity)
    for axis in range(3):
        e = np.eye(3)[axis]*h
        vm2, vm, vp, vp2 = [basis_data(points+j*e, center, widths, mode, degree, carrier)[0]
                            for j in (-2, -1, 1, 2)]
        gradient[:, :, axis] = (vm2-8*vm+8*vp-vp2)/(12*h)
        laplacian += (-vp2+16*vp-30*velocity+16*vm-vm2)/(12*h*h)
    return velocity, gradient, -nu*laplacian, pressure, gp


class Patch:
    def __init__(self, mean, center, widths, carriers, degree, points, weights, tau, build_solvers=True):
        self.mean, self.center, self.widths = mean, center, widths
        self.carriers, self.degree = carriers, degree
        self.points, self.weights = points, weights
        self.modes = list(range(9))
        self.q = (degree+1)**2
        self.cache = [basis_jets(points, center, widths, m, degree, carriers[m],
                                 mean.nu, .0005*np.sqrt(mean.nu*tau)) for m in self.modes]
        self.angles = np.arange(40)*2*np.pi/40
        self.phase = np.exp(1j*np.outer(self.angles, self.modes))
        self.solvers = []
        if not build_solvers:
            return
        for m, cache in enumerate(self.cache):
            V, _, _, _, Gp = cache
            matrix = np.concatenate((V, Gp), axis=2).reshape(-1, 4*self.q)
            if m == 0:
                matrix = matrix.real
            weight = np.repeat(np.sqrt(weights), 3)
            weighted = matrix*weight[:, None]
            scale = np.maximum(np.linalg.norm(weighted, axis=0), 1e-20)
            normalized = weighted/scale
            U, s, Vh = np.linalg.svd(normalized, full_matrices=False)
            ridge = 1e-6*max(s)
            inverse = (Vh.conj().T*(s/(s*s+ridge*ridge)))@U.conj().T
            self.solvers.append((inverse, scale, weight))

    def background(self, t):
        tau = -t
        return jets(self.mean, self.points, tau, .0005*np.sqrt(self.mean.nu*tau), .000005*tau)

    def residual(self, state, background):
        u, g, linear = [np.repeat(a[:, None], len(self.angles), axis=1) for a in background]
        for m, (V, G, L, _, _) in enumerate(self.cache):
            uv = np.einsum('niq,q->ni', V, state[m])
            gv = np.einsum('nijq,q->nij', G, state[m])
            lv = np.einsum('niq,q->ni', L, state[m])
            u += (uv[:, None]*self.phase[None, :, m, None]).real
            g += (gv[:, None]*self.phase[None, :, m, None, None]).real
            linear += (lv[:, None]*self.phase[None, :, m, None]).real
        return linear+np.einsum('naij,naj->nai', g, u)

    def control(self, state, background):
        residual = self.residual(state, background)
        slopes, pressures = [], []
        for m, (inverse, scale, weight) in enumerate(self.solvers):
            source = np.mean(residual*self.phase[None, :, m, None].conj(), axis=1)
            source = source.real if m == 0 else 2*source
            fit = -(inverse@(source.reshape(-1)*weight))/scale
            slopes.append(fit[:3*self.q]); pressures.append(fit[3*self.q:])
        return np.array(slopes), np.array(pressures), residual

    def corrected_residual(self, residual, slopes, pressures):
        result = residual.copy()
        for m, (V, _, _, _, Gp) in enumerate(self.cache):
            change = np.einsum('niq,q->ni', V, slopes[m])+np.einsum('niq,q->ni', Gp, pressures[m])
            result += (change[:, None]*self.phase[None, :, m, None]).real
        return result


def metrics(residual):
    norm = np.linalg.norm(residual, axis=-1)
    return dict(max=float(norm.max()), sample_rms=float(np.sqrt(np.mean(norm**2))))


def pack(values):
    values = np.asarray(values)
    return np.stack((values.real, values.imag), axis=-1).tolist()


class PotentialField:
    def __init__(self, patch, state, pressure):
        self.patch, self.state, self.pressure = patch, state, pressure
        self.nu = patch.mean.nu
        self.interval = None

    def fields(self, points, tau):
        tau = float(np.asarray(tau).ravel()[0])
        t = -tau
        if self.interval is not None and not self.interval[0] <= t <= self.interval[1]:
            raise ValueError('Potential trajectory evaluation outside its saved interval')
        u, p = self.patch.mean.fields(points, tau)
        u, p = u.copy(), p.copy()
        coefficients, pressures = self.state(t), self.pressure(t)
        for m in self.patch.modes:
            V, P, _ = basis_data(points, self.patch.center, self.patch.widths,
                                  m, self.patch.degree, self.patch.carriers[m])
            u += np.einsum('niq,q->ni', V, coefficients[m]).real
            p += np.einsum('nq,q->n', P, pressures[m]).real
        return u, p


def load_saved_field(path=ROOT/'fourier_patch_evolution.json'):
    """Reconstruct the saved experimental field, regardless of acceptance.

    The caller must inspect report['accepted']; reconstruction is not approval.
    """
    report = json.loads(path.read_text(encoding='utf-8'))
    data = report['state_data']
    decode = lambda value: np.asarray(value)[...,0]+1j*np.asarray(value)[...,1]
    _, _, current = build_current()
    mean_data = json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    mean = Trajectory(current, mean_data['nodes'])
    times = data['physical_times']
    state = CubicHermiteSpline([times[0],times[-1]],
        [decode(data['initial_state']),decode(data['final_state'])],
        [decode(data['initial_slope']),decode(data['final_slope'])],extrapolate=False)
    pressure = CubicSpline(times,decode(data['pressures']),extrapolate=False)
    patch = SimpleNamespace(mean=mean,center=tuple(report['center']),widths=tuple(report['widths']),
        degree=data['degree'],modes=report['modes'],carriers={int(m):v for m,v in report['carriers'].items()})
    field = PotentialField(patch,state,pressure)
    field.interval = (times[0],times[-1])
    return field, report


def run():
    inner, _, current = build_current()
    saved = json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    mean = Trajectory(current, saved['nodes'])
    rejected = json.loads((ROOT/'feedback_moving_wave.json').read_text(encoding='utf-8'))
    source = json.loads((ROOT/'compact_potential/midplane_wave_source_k11.json').read_text(encoding='utf-8'))
    k0 = 11.0003
    tau0 = .5*2.**-k0
    t0, t1 = -tau0, -.5*2.**-11.0008
    point = inner.from_similarity([inner.p.X_max*(1+15*.325)**2], [-.0125], tau0)[0]
    center = (float(point[0]), float(point[2]))
    widths = tuple(rejected['widths'])
    pulse1, pulse4 = [source['kelvin_pulses'][i] for i in rejected['source_indices']]
    n1, n4 = [np.array([p['peak_wavevector'][0], p['peak_wavevector'][2]]) for p in (pulse1, pulse4)]
    carriers = {0: np.zeros(2), 1:n1, 2:2*n1, 3:n4-n1, 4:n4,
                5:n4+n1, 6:n4+2*n1, 7:2*n4-n1, 8:2*n4}
    nodes, gauss = leggauss(7)
    axis = .95*nodes
    train = np.array([[center[0]+x*widths[0], 0., center[1]+z*widths[1]] for x in axis for z in axis])
    weight = np.outer(gauss, gauss).ravel()*train[:, 0]
    weight /= sum(weight)
    held_axis = np.linspace(-.8, .8, 6)
    offsets = [(x,z) for x in held_axis for z in held_axis]+[(x,z) for x in (-.97,.97) for z in (-.97,.97)]
    held = np.array([[center[0]+x*widths[0], 0., center[1]+z*widths[1]] for x,z in offsets])
    baseline_train = jets(mean, train, tau0, .0005*np.sqrt(mean.nu*tau0), .000005*tau0)
    baseline_held = jets(mean, held, tau0, .0005*np.sqrt(mean.nu*tau0), .000005*tau0)
    report = dict(accepted=False, scale_recursion_established=False, projections=[], evolution=[],
        coefficient_convention='u = mean - curl(A); physical t = -tau',
        modes=list(range(9)), angular_samples=40, train_points=len(train), held_points=len(held),
        initial_k=k0, center=list(center), widths=list(widths),
        carriers={m:v.tolist() for m,v in carriers.items()},
        scope='Full spatial Fourier potential/pressure collocation with mean and modes 1..8, nonlinear momentum and compact cutoff derivatives. Sample RMS is not volume L2; no whole-domain, finite-energy, temporal endpoint, or recursion acceptance.')
    output = ROOT/'fourier_patch_evolution.json'
    def save():
        output.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    for degree in (2, 4):
        patch = Patch(mean, center, widths, carriers, degree, train, weight, tau0)
        hold = Patch(mean, center, widths, carriers, degree, held, np.ones(len(held))/len(held), tau0)
        state = np.zeros((9,3*patch.q), complex)
        for pulse, covariance_weight in zip((pulse1,pulse4), rejected['weights']):
            m = pulse['mode']; n = np.array([carriers[m][0], m/center[0], carriers[m][1]])
            a = np.array(pulse['peak_amplitude'], complex)
            a -= n*np.dot(n,a)/np.dot(n,n)
            c = 1j*np.cross(n,a)*np.sqrt(covariance_weight)/np.dot(n,n)
            for component in range(3): state[m,component*patch.q] = c[component]
        slope0, pressure0, residual = patch.control(state, baseline_train)
        held_residual = hold.residual(state, baseline_held)
        held_corrected = hold.corrected_residual(held_residual, slope0, pressure0)
        row = dict(degree=degree, train_before=metrics(residual),
            train_after=metrics(patch.corrected_residual(residual,slope0,pressure0)),
            held_before=metrics(held_residual), held_after=metrics(held_corrected))
        report['projections'].append(row); save(); print(json.dumps(row),flush=True)
        if row['held_after']['max'] >= .95*row['held_before']['max']:
            continue
        # Check the spectral/cylindrical assembly against a callable Cartesian
        # field, with a real time derivative rather than the cached projection.
        check_space = [0,14,27,35]
        check_angles = [1,7,13,19,25,31,37]
        check_points = np.array([[held[i,0]*np.cos(patch.angles[j]),
                                  held[i,0]*np.sin(patch.angles[j]),held[i,2]]
                                 for i in check_space for j in check_angles])
        taylor = PotentialField(patch, lambda t: state+(t-t0)*slope0, lambda t: pressure0)
        direct = cylindrical_residual(momentum(jets(taylor,check_points,tau0,
            .0005*np.sqrt(mean.nu*tau0),.000005*tau0)),check_points)
        expected = held_corrected[check_space][:,check_angles].reshape(-1,3)
        difference = float(max(np.linalg.norm(direct-expected,axis=1)))
        reference = max(float(max(np.linalg.norm(expected,axis=1))),1.)
        report['cartesian_assembly_check'] = dict(max_difference=difference,relative_difference=difference/reference)
        save();print(json.dumps(report['cartesian_assembly_check']),flush=True)
        if difference/reference > 2e-4:
            report['stopped_reason']='Direct Cartesian replay disagrees with the cached spectral assembly.'
            save();return
        # Advance the actual potential state with an explicit midpoint step.
        dt = t1-t0
        midpoint_state = state+dt/2*slope0
        mid_slope, mid_pressure, _ = patch.control(midpoint_state, patch.background(t0+dt/2))
        final_state = state+dt*mid_slope
        final_slope, final_pressure, _ = patch.control(final_state, patch.background(t1))
        state_path = CubicHermiteSpline([t0,t1], np.array([state,final_state]),
                                       np.array([slope0,final_slope]), extrapolate=False)
        pressure_path = CubicSpline([t0,t0+dt/2,t1], np.array([pressure0,mid_pressure,final_pressure]),extrapolate=False)
        field = PotentialField(patch,state_path,pressure_path)
        field.interval = (t0,t1)
        report['state_data'] = dict(degree=degree, physical_times=[t0,t0+dt/2,t1],
            initial_state=pack(state), final_state=pack(final_state),
            initial_slope=pack(slope0), final_slope=pack(final_slope),
            pressures=pack(np.array([pressure0,mid_pressure,final_pressure])))
        # Shift angle samples away from the fitting grid for direct Cartesian replay.
        angles = (np.arange(40)+.37)*2*np.pi/40
        direct_points = np.array([[p[0]*np.cos(a),p[0]*np.sin(a),p[2]] for p in held for a in angles])
        for fraction in (.25,.5,.75):
            t = t0+dt*fraction; tau=-t
            jet = jets(field,direct_points,tau,.0005*np.sqrt(mean.nu*tau),.000005*tau)
            value = momentum(jet)
            row = dict(degree=degree, k=float(-np.log2(2*tau)),fraction=fraction,
                       direct_momentum=metrics(value),sampled_divergence=float(max(abs(np.trace(jet[1],axis1=1,axis2=2)))))
            report['evolution'].append(row);save();print(json.dumps(row),flush=True)
            if row['direct_momentum']['max'] > report['projections'][-1]['held_before']['max']:
                report['rejected']=True
                report['stopped_reason']='Time-evolved full momentum exceeds initial uncorrected wave holdout maximum.'
                save();return
        return
    report['stopped_reason']='Neither polynomial degree improves independent spatial maximum by 5 percent.'
    save()


if __name__ == '__main__':
    run()
