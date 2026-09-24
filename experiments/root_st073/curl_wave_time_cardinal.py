"""Three-time derivative-cardinal collocation on the second patch interval.

Each compact time potential has zero value at all three collocation times and
zero value/derivative at both endpoints. Its time derivative is cardinal at
the three interior times. The pressure basis is value-cardinal there.
"""
import json

import numpy as np
from numpy.polynomial import Polynomial

from curl_wave_patch_evolution import fit_slope, metrics, sample_residual
from curl_wave_patch_trajectory import HermiteIntervalField
from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


TIME_NODES = np.array([.25, .5, .75])


def cardinal_polynomials():
    variable = Polynomial([0., 1.])
    envelope = variable**2*(1-variable)**2
    zeros = Polynomial(envelope.coef.copy())
    for node in TIME_NODES:
        zeros = zeros*(variable-node)
    derivative_bubbles, pressure_bubbles = [], []
    for j, node in enumerate(TIME_NODES):
        lagrange = Polynomial([1.])
        for k, other in enumerate(TIME_NODES):
            if k != j:
                lagrange *= (variable-other)/(node-other)
        derivative_bubbles.append(zeros*lagrange/zeros.deriv()(node))
        pressure_bubbles.append(envelope*lagrange/envelope(node))
    return derivative_bubbles, pressure_bubbles


class CardinalIntervalField(HermiteIntervalField):
    def __init__(self, base, wave, amplitude, stage0, stage1, time_step,
                 initial_harmonic, initial_mean, harmonic_fits, mean_fits):
        super().__init__(base, wave, amplitude, stage0, stage1, time_step,
                         initial_harmonic, initial_mean)
        self.harmonic_fits = harmonic_fits
        self.mean_fits = mean_fits
        self.derivative_bubbles, self.pressure_bubbles = cardinal_polynomials()

    def coefficients(self, t):
        harmonic, mean, pressures, mean_pressure = super().coefficients(t)
        s = (t-self.tau0)/self.time_step
        for derivative, pressure, harmonic_fit, mean_fit in zip(
                self.derivative_bubbles, self.pressure_bubbles,
                self.harmonic_fits, self.mean_fits):
            potential_factor = self.time_step*derivative(s)
            pressure_factor = pressure(s)
            harmonic = [current+potential_factor*fit[:27]
                        for current, fit in zip(harmonic, harmonic_fit)]
            mean = mean+potential_factor*mean_fit[:18]
            pressures = [current+pressure_factor*fit[27:]
                         for current, fit in zip(pressures, harmonic_fit)]
            mean_pressure = mean_pressure+pressure_factor*mean_fit[18:]
        return harmonic, mean, pressures, mean_pressure


def direct_metrics(field, points, tau, dt):
    hs = .0005*np.sqrt(field.nu*tau)
    ht = dt/64
    u, gradient, part = kinematics(field, points, tau, hs, ht)
    residual = part+np.einsum('nij,nj->ni', gradient, u)
    norms = np.linalg.norm(residual, axis=1)
    return {'max_momentum': float(norms.max()),
            'rms_momentum': float(np.sqrt(np.mean(norms**2))),
            'max_fd_divergence': float(np.max(np.abs(np.trace(gradient,
                                                             axis1=1,
                                                             axis2=2))))}


def run():
    evolution = json.loads((ROOT/'compact_potential'/'curl_wave_patch_evolution_3stage.json').read_text())
    middle = json.loads((ROOT/'compact_potential'/'curl_wave_temporal_bubble.json').read_text())
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    dt = evolution['time_step']
    stage0, stage1, stage2 = evolution['stages']
    initial_harmonic = [dt*np.array([complex(re, im) for re, im in row])[:27]
                        for row in stage0['harmonic_slope_coefficients']]
    initial_mean = dt*np.array(stage0['mean_slope_coefficients'][:18])
    hermite = HermiteIntervalField(base, wave, evolution['amplitude'],
                                   stage1, stage2, dt,
                                   initial_harmonic, initial_mean)
    angles = np.arange(8)*2*np.pi/8
    train_axis = np.linspace(-.6, .6, 5)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    harmonic_fits, mean_fits, train_projection = [], [], []
    for j, fraction in enumerate(TIME_NODES):
        tau = stage1['tau']+fraction*dt
        if j == 1:
            harmonic_fit = [np.array([complex(re, im) for re, im in row])
                            for row in middle['harmonic_slope_corrections']]
            mean_fit = np.array(middle['mean_slope_correction'])
            projected = middle['midpoint_projected']
        else:
            rows = sample_residual(hermite, wave, tau, train_grid,
                                   angles, dt)
            harmonic_fit, mean_fit = fit_slope(rows, wave, tau, angles)
            projected = metrics(rows, wave, tau, angles,
                                harmonic_fit, mean_fit)
        harmonic_fits.append(harmonic_fit)
        mean_fits.append(mean_fit)
        train_projection.append(projected)
    cardinal = CardinalIntervalField(base, wave, evolution['amplitude'],
                                     stage1, stage2, dt,
                                     initial_harmonic, initial_mean,
                                     harmonic_fits, mean_fits)
    def points_for(axis):
        return np.vstack([np.column_stack((r*np.cos(angles), r*np.sin(angles),
                                           np.full(len(angles), z)))
                          for xi in axis for eta in axis
                          for r, z in [(wave.radius+wave.radial_halfwidth*xi,
                                        wave.zcenter+wave.axial_halfwidth*eta)]])
    heldout = points_for(np.array([-.45, -.15, .15, .45]))
    subset = points_for(np.array([-.15, .15]))
    def screen(fraction, points):
        tau = stage1['tau']+fraction*dt
        return {'fraction': fraction, 'tau': tau,
                'baseline': direct_metrics(base, points, tau, dt),
                'hermite': direct_metrics(hermite, points, tau, dt),
                'cardinal': direct_metrics(cardinal, points, tau, dt)}
    at_nodes = [screen(float(s), heldout) for s in TIME_NODES]
    between_nodes = [screen(s, subset) for s in (.125, .375, .625, .875)]
    endpoint_mismatch = []
    collocation_velocity_mismatch = []
    for fraction in (0., 1.):
        old = hermite.coefficients(stage1['tau']+fraction*dt)
        new = cardinal.coefficients(stage1['tau']+fraction*dt)
        endpoint_mismatch.extend(np.max(np.abs(x-y)) for x, y in zip(
            old[0]+[old[1]]+old[2]+[old[3]],
            new[0]+[new[1]]+new[2]+[new[3]]))
    for fraction in TIME_NODES:
        old = hermite.coefficients(stage1['tau']+fraction*dt)
        new = cardinal.coefficients(stage1['tau']+fraction*dt)
        collocation_velocity_mismatch.extend(np.max(np.abs(x-y))
                                             for x, y in zip(old[0]+[old[1]],
                                                             new[0]+[new[1]]))
    derivative_bubbles, pressure_bubbles = cardinal_polynomials()
    derivative_matrix = [[float(item.deriv()(node)) for item in derivative_bubbles]
                         for node in TIME_NODES]
    pressure_matrix = [[float(item(node)) for item in pressure_bubbles]
                       for node in TIME_NODES]
    report = {'time_step': dt, 'interval_index': 1,
              'time_nodes': TIME_NODES.tolist(),
              'train_spatial_nodes': len(train_grid),
              'heldout_spatial_nodes': 16,
              'train_projection': train_projection,
              'at_nodes': at_nodes, 'between_nodes_subset': between_nodes,
              'endpoint_coefficient_max_mismatch': float(max(endpoint_mismatch)),
              'collocation_potential_max_mismatch': float(max(collocation_velocity_mismatch)),
              'derivative_cardinal_matrix': derivative_matrix,
              'pressure_cardinal_matrix': pressure_matrix,
              'harmonic_fits': [[[[float(v.real), float(v.imag)] for v in row]
                                 for row in fit] for fit in harmonic_fits],
              'mean_fits': [fit.tolist() for fit in mean_fits],
              'scope': 'Three interior-time derivative/pressure collocation conditions on the second Hermite interval. Direct full momentum is checked at all three times on held-out spatial nodes and between times on four-node subsets. No whole-support/time bound.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_time_cardinal.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
