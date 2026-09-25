"""Full Cartesian residual of one compact exact-curl wave on the new mean.

The wave uses a local two-Kelvin covariance pair, but its envelope and
viscous terms are retained in the actual physical finite-difference operator.
Only an angular ring at one space/time center is tested here.
"""

import json

import numpy as np

from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from curl_wave_prototype import LocalizedCurlWave, cylindrical_residual
from delayed_multimode_cone_fit import BASE_NAME
from delayed_outer_radial_extension import make_extended_modes
from delayed_similarity_curl_screen import CurlPatchedLift
from delayed_similarity_pressure_screen import SimilarityPressurePatch
from joint_collar_fit import kinematics
from radial_continuation import ROOT


class PressureMeanField:
    def __init__(self, mean, pressure):
        self.mean = mean
        self.pressure = pressure
        self.nu = mean.nu

    def fields(self, points, tau):
        velocity, scalar = self.mean.fields(points, tau)
        values, _ = self.pressure.basis(points, tau)
        return velocity, scalar+values@self.pressure.coefficients


def run():
    source = json.loads((ROOT/'delayed_outer_kelvin_source.json').read_text())
    scale = json.loads((ROOT/'delayed_outer_wave_scale.json').read_text())
    mean_source = json.loads((ROOT/'delayed_outer_admission_search.json').read_text())
    pressure_source = json.loads((ROOT/'delayed_outer_pressure_fit.json').read_text())
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = CurlPatchedLift(base, make_extended_modes(base),
                           mean_source['coefficients'])
    pressure = SimilarityPressurePatch(
        base, pressure_source['pressure_coefficients'])
    field = PressureMeanField(mean, pressure)
    wave = LocalizedCurlWave(
        source, scale['radial_halfwidth'], scale['axial_halfwidth'],
        time_halfwidth=scale['time_halfwidth'])
    tau = source['tau']
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles),
                              radius*np.sin(angles),
                              np.full(len(angles), z)))
    hs = scale['radial_halfwidth']/20
    ht = .0001*tau
    u, grad, part = kinematics(field, points, tau, hs, ht)
    w, wgrad, wpart = kinematics(wave, points, tau, hs, ht)
    baseline = part+np.einsum('pab,pb->pa', grad, u)
    linear = (wpart+np.einsum('pab,pb->pa', grad, w)
              +np.einsum('pab,pb->pa', wgrad, u))
    quadratic = np.einsum('pab,pb->pa', wgrad, w)
    cylindrical_wave = np.column_stack((
        np.cos(angles)*w[:, 0]+np.sin(angles)*w[:, 1],
        -np.sin(angles)*w[:, 0]+np.cos(angles)*w[:, 1],
        w[:, 2]))
    covariance = np.mean(cylindrical_wave[:, 0, None]
                         *cylindrical_wave[:, 1:], axis=0)
    target = np.asarray(source['local_tangential_stress_primitive'])

    def metrics(amplitude):
        R = baseline+amplitude*linear+amplitude**2*quadratic
        norm = np.linalg.norm(R, axis=1)
        return dict(amplitude=float(amplitude),
                    covariance_fraction=float(amplitude**2),
                    maximum=float(np.max(norm)),
                    rms=float(np.sqrt(np.mean(norm**2))),
                    mean_cylindrical=(
                        cylindrical_residual(R, points).mean(axis=0).tolist()))

    samples = [metrics(a) for a in (0., .001, .003, .01, .03,
                                     .1, .3, 1.)]
    dense = [metrics(a) for a in np.r_[0., np.geomspace(1e-5, 1., 101)]]
    report = dict(source='delayed_outer_kelvin_source.json',
                  scale_source='delayed_outer_wave_scale.json',
                  tau=tau, point=source['point'],
                  spatial_step=hs, temporal_step=ht,
                  angular_samples=len(angles),
                  target_covariance=target.tolist(),
                  exact_curl_covariance=covariance.tolist(),
                  covariance_relative_error=float(
                      np.linalg.norm(covariance-target)/np.linalg.norm(target)),
                  samples=samples,
                  best_rms=min(dense, key=lambda row: row['rms']),
                  best_max=min(dense, key=lambda row: row['maximum']),
                  scope='Exact-curl nonaxisymmetric two-mode wave on a '
                        'ten-mode mean plus compact pressure, evaluated '
                        'with complete nonlinear Cartesian momentum on '
                        'one angular ring. No phase/amplitude evolution, '
                        'continuous cone, moment closure, volume L2, '
                        'or PDE acceptance.', accepted=False)
    (ROOT/'delayed_outer_wave_residual.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps(dict(covariance_relative_error=(
                              report['covariance_relative_error']),
                          baseline=samples[0], full_wave=samples[-1],
                          best_rms=report['best_rms'],
                          best_max=report['best_max'])), flush=True)


if __name__ == '__main__':
    run()
