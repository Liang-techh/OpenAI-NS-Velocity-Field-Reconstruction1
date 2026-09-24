"""Direct midpoint residual of an evolved exact-curl patch interval.

The velocity interpolates potential coefficients linearly between the first
two collocation nodes. This tests the actual continuous-time field between
nodes, rather than only the instantaneous projected endpoint residual.
"""
import argparse
import json

import numpy as np

from curl_wave_patch_evolution import FrozenPotentialField, polynomial_data
from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


class FirstIntervalField:
    def __init__(self, base, wave, amplitude, stage,
                 initial_harmonic=None, initial_mean=None):
        self.base, self.wave, self.amplitude = base, wave, amplitude
        self.nu = base.nu
        self.tau0 = stage['tau']
        self.harmonic_slopes = [np.array([complex(re, im) for re, im in row])
                                for row in stage['harmonic_slope_coefficients']]
        self.mean_slope = np.array(stage['mean_slope_coefficients'])
        self.initial_harmonic = (initial_harmonic if initial_harmonic is not None
                                 else [np.zeros(27, complex)
                                       for _ in self.harmonic_slopes])
        self.initial_mean = (initial_mean if initial_mean is not None
                             else np.zeros(18))

    def coefficients(self, t):
        interval = t-self.tau0
        harmonic = [initial+interval*s[:27]
                    for initial, s in zip(self.initial_harmonic,
                                          self.harmonic_slopes)]
        mean = self.initial_mean+interval*self.mean_slope[:18]
        pressures = [s[27:] for s in self.harmonic_slopes]
        mean_pressure = self.mean_slope[18:]
        return harmonic, mean, pressures, mean_pressure

    def fields(self, points, tau):
        pts = np.asarray(points, float)
        ts = np.broadcast_to(tau, (len(pts),))
        if not np.all(ts == ts[0]):
            raise ValueError('Evaluate one time per spatial batch')
        t = float(ts[0])
        harmonic, mean, pressures, mean_pressure = self.coefficients(t)
        frozen = FrozenPotentialField(self.base, self.wave, self.amplitude,
                                      harmonic, mean)
        velocity, pressure = frozen.fields(pts, t)
        for i, (x, y, z) in enumerate(pts):
            r = np.hypot(x, y)
            if r == 0 or abs(r-self.wave.radius) >= self.wave.radial_halfwidth or abs(z-self.wave.zcenter) >= self.wave.axial_halfwidth:
                continue
            theta = np.arctan2(y, x)
            for mode, p_coeff in zip(self.wave.waves, pressures):
                m = mode['m']
                kr, _, kz = mode['normal']
                p = polynomial_data(self.wave, r, z, p_coeff, 1)[0][0]
                carrier = np.exp(1j*(m*theta+kr*(r-self.wave.radius)
                                     +kz*(z-self.wave.zcenter)
                                     +mode['omega']*(t-self.wave.tau0)))
                pressure[i] += (carrier*p).real
            pressure[i] += polynomial_data(self.wave, r, z,
                                           mean_pressure, 1)[0][0]
        return velocity, pressure


class HermiteIntervalField(FirstIntervalField):
    """Endpoint derivative-matched potentials and continuous pressure."""
    def __init__(self, base, wave, amplitude, stage0, stage1, time_step,
                 initial_harmonic=None, initial_mean=None):
        super().__init__(base, wave, amplitude, stage0,
                         initial_harmonic, initial_mean)
        self.time_step = time_step
        self.next_harmonic_slopes = [np.array([complex(re, im)
                                               for re, im in row])
                                     for row in stage1['harmonic_slope_coefficients']]
        self.next_mean_slope = np.array(stage1['mean_slope_coefficients'])

    def coefficients(self, t):
        s = (t-self.tau0)/self.time_step
        h10 = s**3-2*s*s+s
        h01 = -2*s**3+3*s*s
        h11 = s**3-s*s
        a = self.time_step*(h10+h01)
        b = self.time_step*h11
        harmonic = [initial+a*old[:27]+b*new[:27]
                    for initial, old, new in zip(self.initial_harmonic,
                                                 self.harmonic_slopes,
                                                 self.next_harmonic_slopes)]
        mean = (self.initial_mean+a*self.mean_slope[:18]
                +b*self.next_mean_slope[:18])
        pressures = [(1-s)*old[27:]+s*new[27:]
                     for old, new in zip(self.harmonic_slopes,
                                         self.next_harmonic_slopes)]
        mean_pressure = ((1-s)*self.mean_slope[18:]
                         +s*self.next_mean_slope[18:])
        return harmonic, mean, pressures, mean_pressure


def residual_metrics(field, points, tau, time_step):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = min(.0001*tau, time_step/8)
    u, grad, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', grad, u)
    norm = np.linalg.norm(residual, axis=1)
    return {'max_momentum': float(norm.max()),
            'rms_momentum': float(np.sqrt(np.mean(norm**2))),
            'max_fd_divergence': float(np.max(np.abs(np.trace(grad,
                                                             axis1=1,
                                                             axis2=2)))),
            'max_speed': float(np.max(np.linalg.norm(u, axis=1)))}


def run(evolution_name='curl_wave_patch_evolution.json',
        output_name='curl_wave_patch_trajectory.json', interval_index=0):
    evolution = json.loads((ROOT/'compact_potential'/evolution_name).read_text())
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    stages = evolution['stages']
    if interval_index < 0 or interval_index >= len(stages)-1:
        raise ValueError('Interval index requires two adjacent stages')
    dt = evolution['time_step']
    initial_harmonic = [np.zeros(27, complex) for _ in wave.waves]
    initial_mean = np.zeros(18)
    for stage in stages[:interval_index]:
        slopes = [np.array([complex(re, im) for re, im in row])
                  for row in stage['harmonic_slope_coefficients']]
        initial_harmonic = [current+dt*slope[:27]
                            for current, slope in zip(initial_harmonic, slopes)]
        initial_mean += dt*np.array(stage['mean_slope_coefficients'][:18])
    field = FirstIntervalField(base, wave, evolution['amplitude'],
                               stages[interval_index],
                               initial_harmonic, initial_mean)
    hermite = HermiteIntervalField(base, wave, evolution['amplitude'],
                                   stages[interval_index],
                                   stages[interval_index+1], dt,
                                   initial_harmonic, initial_mean)
    tau = stages[interval_index]['tau']+dt/2
    angles = np.arange(8)*2*np.pi/8
    axis = np.array([-.45, -.15, .15, .45])
    points = np.vstack([np.column_stack((r*np.cos(angles), r*np.sin(angles),
                                         np.full(len(angles), z)))
                        for xi in axis for eta in axis
                        for r, z in [(wave.radius+wave.radial_halfwidth*xi,
                                      wave.zcenter+wave.axial_halfwidth*eta)]])
    report = {'tau': tau, 'time_step': dt,
              'interval_index': interval_index,
              'heldout_spatial_nodes': len(axis)**2, 'angles_per_node': len(angles),
              'temporal_fd_step': min(.0001*tau, evolution['time_step']/8),
              'baseline': residual_metrics(base, points, tau,
                                           evolution['time_step']),
              ('evolved_first_interval' if interval_index == 0 else 'evolved_second_interval'):
                  residual_metrics(field, points, tau, dt),
              ('hermite_first_interval' if interval_index == 0 else 'hermite_second_interval'):
                  residual_metrics(hermite, points, tau, dt),
              'scope': 'Direct finite-difference full momentum at the midpoint of one linear and Hermite potential-coefficient interval. The temporal derivative stencil stays within the interval. This is a spatially held-out sample, not an inter-node/global bound or a matched pulse.',
              'accepted': False}
    if interval_index > 0:
        previous_initial_harmonic = [np.zeros(27, complex)
                                     for _ in wave.waves]
        previous_initial_mean = np.zeros(18)
        for stage in stages[:interval_index-1]:
            slopes = [np.array([complex(re, im) for re, im in row])
                      for row in stage['harmonic_slope_coefficients']]
            previous_initial_harmonic = [current+dt*slope[:27]
                                         for current, slope in zip(previous_initial_harmonic,
                                                                   slopes)]
            previous_initial_mean += dt*np.array(stage['mean_slope_coefficients'][:18])
        previous = HermiteIntervalField(base, wave, evolution['amplitude'],
                                        stages[interval_index-1],
                                        stages[interval_index], dt,
                                        previous_initial_harmonic,
                                        previous_initial_mean)
        left = previous.coefficients(stages[interval_index]['tau'])
        right = hermite.coefficients(stages[interval_index]['tau'])
        mismatches = [np.max(np.abs(x-y)) for x, y in zip(left[0], right[0])]
        mismatches.extend([np.max(np.abs(left[1]-right[1]))])
        mismatches.extend([np.max(np.abs(x-y)) for x, y in zip(left[2], right[2])])
        mismatches.extend([np.max(np.abs(left[3]-right[3]))])
        report['interface_coefficient_max_mismatch'] = float(max(mismatches))
    out = ROOT/'compact_potential'/output_name
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--evolution-name', default='curl_wave_patch_evolution.json')
    parser.add_argument('--output-name', default='curl_wave_patch_trajectory.json')
    parser.add_argument('--interval-index', type=int, default=0)
    args = parser.parse_args()
    run(args.evolution_name, args.output_name, args.interval_index)
