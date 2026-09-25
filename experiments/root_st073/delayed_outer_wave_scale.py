"""Check exact-curl carrier/cutoff and viscosity scales on the new cone.

The envelope fits inside the sampled X=1.012..1.02, eta=.2..3 window
around the selected local Kelvin pair. Even that window has not been
certified continuously; this is a necessary physical scale diagnostic.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from curl_wave_prototype import LocalizedCurlWave
from curl_wave_scale_audit import pieces
from delayed_multimode_cone_fit import BASE_NAME
from radial_continuation import ROOT


def run():
    source = json.loads((ROOT/'delayed_outer_kelvin_source.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    tau = source['tau']
    X = np.array([1.012, 1.016, 1.02, 1.016, 1.016])
    eta = np.array([.25, .25, .25, .2, .3])
    physical = base.compact.joined.inner.from_similarity(X, eta, tau)
    radial_ceiling = float(min(physical[1, 0]-physical[0, 0],
                               physical[2, 0]-physical[1, 0]))
    axial_ceiling = float(min(physical[1, 2]-physical[3, 2],
                              physical[4, 2]-physical[1, 2]))
    radial_halfwidth = .8*radial_ceiling
    axial_halfwidth = .8*axial_ceiling
    time_halfwidth = .05*tau
    wave = LocalizedCurlWave(source, radial_halfwidth,
                             axial_halfwidth,
                             time_halfwidth=time_halfwidth)
    axis = np.linspace(-.6, .6, 5)
    rows = []
    duration = .1*tau
    for mode in wave.waves:
        leading_sq = remainder_sq = 0.
        ratios = []
        for xr in axis:
            for za in axis:
                r = wave.radius+xr*radial_halfwidth
                z = wave.zcenter+za*axial_halfwidth
                leading, remainder = pieces(wave, mode, r, z)
                a, b = np.linalg.norm(leading), np.linalg.norm(remainder)
                ratios.append(float(b/a))
                leading_sq += a*a
                remainder_sq += b*b
        normal_norm = float(np.linalg.norm(mode['normal']))
        minimum_multiplier = int(np.ceil(max(ratios)))
        largest_duration = 1/(wave.nu*(minimum_multiplier*normal_norm)**2)
        rows.append(dict(angular_mode=mode['m'],
                         normal_norm=normal_norm,
                         radial_phase_across_halfwidth=float(
                             abs(mode['normal'][0])*radial_halfwidth),
                         axial_phase_across_halfwidth=float(
                             abs(mode['normal'][2])*axial_halfwidth),
                         carrier_remainder_grid_max=max(ratios),
                         carrier_remainder_grid_rms=float(
                             np.sqrt(remainder_sq/leading_sq)),
                         integer_multiplier_for_remainder_at_most_carrier=(
                             minimum_multiplier),
                         damping_exponent_over_pulse_at_multiplier_one=(
                             wave.nu*normal_norm**2*duration),
                         damping_exponent_over_pulse_at_minimum_multiplier=(
                             wave.nu*(minimum_multiplier*normal_norm)**2
                             *duration),
                         largest_pulse_duration_at_minimum_multiplier=(
                             largest_duration),
                         both_lenient_conditions_feasible=(
                             largest_duration >= duration)))
    report = dict(source='delayed_outer_kelvin_source.json',
                  tau=tau, pulse_duration=duration,
                  sampled_window=dict(X=(1.012, 1.02), eta=(.2, .3)),
                  radial_ceiling=radial_ceiling,
                  axial_ceiling=axial_ceiling,
                  radial_halfwidth=radial_halfwidth,
                  axial_halfwidth=axial_halfwidth,
                  radial_diffusion_time=radial_halfwidth**2/wave.nu,
                  pulse_to_radial_diffusion_time=(
                      duration*wave.nu/radial_halfwidth**2),
                  time_halfwidth=time_halfwidth,
                  covariance_fit_relative_error=wave.covariance_error,
                  modes=rows,
                  scope='One local two-harmonic exact-curl prototype '
                        'inside a sampled cone window. Carrier/remainder '
                        'and viscous-duration conditions are lenient '
                        'necessary diagnostics, not paper asymptotic '
                        'certificates or full momentum acceptance.',
                  accepted=False)
    (ROOT/'delayed_outer_wave_scale.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(radial_halfwidth=radial_halfwidth,
                          axial_halfwidth=axial_halfwidth,
                          modes=rows)), flush=True)


if __name__ == '__main__':
    run()
