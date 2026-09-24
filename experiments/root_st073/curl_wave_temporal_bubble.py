"""Space-time midpoint correction with zero endpoint potential and pressure.

An odd temporal bubble changes the potential time derivative at the midpoint
without changing endpoint values/slopes or the midpoint velocity. An even
pressure bubble supplies the midpoint pressure fit. Both vanish smoothly at
the time-element endpoints, retaining the Hermite interface conditions.
"""
import json

import numpy as np

from curl_wave_patch_evolution import fit_slope, metrics, sample_residual
from curl_wave_patch_trajectory import HermiteIntervalField, residual_metrics
from curl_wave_prototype import LocalizedCurlWave
from joined_field import ROOT
from radial_peak_cone import current_field


class TemporalBubbleField(HermiteIntervalField):
    def __init__(self, base, wave, amplitude, stage0, stage1, time_step,
                 initial_harmonic, initial_mean,
                 harmonic_corrections, mean_correction):
        super().__init__(base, wave, amplitude, stage0, stage1, time_step,
                         initial_harmonic, initial_mean)
        self.harmonic_corrections = harmonic_corrections
        self.mean_correction = mean_correction

    def coefficients(self, t):
        harmonic, mean, pressures, mean_pressure = super().coefficients(t)
        s = (t-self.tau0)/self.time_step
        even = 16*s*s*(1-s)**2
        odd = even*(2*s-1)
        factor = self.time_step*odd/2
        harmonic = [current+factor*delta[:27]
                    for current, delta in zip(harmonic,
                                              self.harmonic_corrections)]
        mean = mean+factor*self.mean_correction[:18]
        pressures = [current+even*delta[27:]
                     for current, delta in zip(pressures,
                                               self.harmonic_corrections)]
        mean_pressure = mean_pressure+even*self.mean_correction[18:]
        return harmonic, mean, pressures, mean_pressure


def run():
    evolution = json.loads((ROOT/'compact_potential'/'curl_wave_patch_evolution_3stage.json').read_text())
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    dt = evolution['time_step']
    stage0, stage1, stage2 = evolution['stages']
    initial_harmonic = [dt*np.array([complex(re, im) for re, im in row])[:27]
                        for row in stage0['harmonic_slope_coefficients']]
    initial_mean = dt*np.array(stage0['mean_slope_coefficients'][:18])
    original = HermiteIntervalField(base, wave, evolution['amplitude'],
                                    stage1, stage2, dt,
                                    initial_harmonic, initial_mean)
    midpoint = stage1['tau']+dt/2
    angles = np.arange(8)*2*np.pi/8
    train_axis = np.linspace(-.6, .6, 5)
    train_grid = [(x, y) for x in train_axis for y in train_axis]
    rows = sample_residual(original, wave, midpoint, train_grid, angles, dt)
    harmonic_fit, mean_fit = fit_slope(rows, wave, midpoint, angles)
    corrected = TemporalBubbleField(base, wave, evolution['amplitude'],
                                     stage1, stage2, dt,
                                     initial_harmonic, initial_mean,
                                     harmonic_fit, mean_fit)
    test_axis = np.array([-.45, -.15, .15, .45])
    def points_for(axis):
        return np.vstack([np.column_stack((r*np.cos(angles), r*np.sin(angles),
                                           np.full(len(angles), z)))
                          for xi in axis for eta in axis
                          for r, z in [(wave.radius+wave.radial_halfwidth*xi,
                                        wave.zcenter+wave.axial_halfwidth*eta)]])
    heldout = points_for(test_axis)
    subset = points_for(np.array([-.15, .15]))
    def direct_at(tau, points):
        return {'tau': tau, 'baseline': residual_metrics(base, points, tau, dt),
                'hermite': residual_metrics(original, points, tau, dt),
                'with_temporal_bubble': residual_metrics(corrected, points,
                                                         tau, dt)}
    direct_midpoint = direct_at(midpoint, heldout)
    off_midpoints = [direct_at(stage1['tau']+fraction*dt, subset)
                     for fraction in (.25, .75)]
    left_original = original.coefficients(stage1['tau'])
    left_corrected = corrected.coefficients(stage1['tau'])
    right_original = original.coefficients(stage2['tau'])
    right_corrected = corrected.coefficients(stage2['tau'])
    endpoint_mismatch = max(float(np.max(np.abs(x-y)))
                            for pair in ((left_original, left_corrected),
                                         (right_original, right_corrected))
                            for x, y in zip(pair[0][0]+[pair[0][1]]+pair[0][2]+[pair[0][3]],
                                            pair[1][0]+[pair[1][1]]+pair[1][2]+[pair[1][3]]))
    report = {'time_step': dt, 'interval_index': 1,
              'train_spatial_nodes': len(train_grid),
              'heldout_spatial_nodes': len(test_axis)**2,
              'angles_per_node': len(angles),
              'midpoint_projected': metrics(rows, wave, midpoint, angles,
                                            harmonic_fit, mean_fit),
              'direct_midpoint': direct_midpoint,
              'off_midpoint_subset': off_midpoints,
              'endpoint_coefficient_max_mismatch': endpoint_mismatch,
              'harmonic_slope_corrections': [[[float(v.real), float(v.imag)]
                                              for v in row] for row in harmonic_fit],
              'mean_slope_correction': mean_fit.tolist(),
              'scope': 'Midpoint space-time collocation bubble on the second Hermite interval. Direct full momentum is sampled at the midpoint and on four held-out nodes at quarter points. Endpoints remain matched. No whole-interval/full-support bound.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_temporal_bubble.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
