"""Full nonlinear momentum versus covariance fraction for the exact-curl wave."""
import json

import numpy as np

from curl_wave_prototype import LocalizedCurlWave, cylindrical_residual
from joined_field import ROOT
from joint_collar_fit import kinematics
from radial_peak_cone import current_field


def run():
    source = json.loads((ROOT/'compact_potential'/'radial_peak_cone.json').read_text())
    base = current_field()
    wave = LocalizedCurlWave(source)
    tau = source['tau']
    radius, _, z = source['point']
    angles = np.arange(16)*2*np.pi/16
    points = np.column_stack((radius*np.cos(angles), radius*np.sin(angles),
                              np.full(len(angles), z)))
    hs = .0005*np.sqrt(base.nu*tau)
    ht = .0001*tau
    u, grad, part = kinematics(base, points, tau, hs, ht)
    w, wgrad, wpart = kinematics(wave, points, tau, hs, ht)
    baseline = part + np.einsum('nij,nj->ni', grad, u)
    linear = (wpart + np.einsum('nij,nj->ni', grad, w)
              + np.einsum('nij,nj->ni', wgrad, u))
    quadratic = np.einsum('nij,nj->ni', wgrad, w)
    linear_cyl = cylindrical_residual(linear, points)
    quadratic_cyl = cylindrical_residual(quadratic, points)
    def metrics(amplitude):
        cart = baseline + amplitude*linear + amplitude**2*quadratic
        cyl = cylindrical_residual(cart, points)
        norms = np.linalg.norm(cyl, axis=1)
        return {'amplitude': float(amplitude),
                'covariance_fraction': float(amplitude**2),
                'max_momentum': float(norms.max()),
                'rms_momentum': float(np.sqrt(np.mean(norms**2))),
                'mean_cylindrical': cyl.mean(axis=0).tolist()}
    rows = [metrics(x) for x in (0., .001, .003, .01, .03, .1, .3, 1.)]
    candidates = np.r_[0., np.geomspace(1e-5, 1., 101)]
    best = min((metrics(x) for x in candidates), key=lambda row: row['rms_momentum'])
    report = {'tau': tau, 'point': source['point'], 'rows': rows,
              'best_sampled_rms': best,
              'mean_linear_cylindrical': linear_cyl.mean(axis=0).tolist(),
              'mean_quadratic_cylindrical': quadratic_cyl.mean(axis=0).tolist(),
              'max_linear_norm': float(np.max(np.linalg.norm(linear, axis=1))),
              'max_quadratic_norm': float(np.max(np.linalg.norm(quadratic, axis=1))),
              'scope': 'Exact affine residual identity for base+amplitude*wave, retaining quadratic self-transport. Sixteen angles at one physical point; angular covariance fraction is amplitude squared. No global PDE gate.',
              'accepted': False}
    out = ROOT/'compact_potential'/'curl_wave_amplitude.json'
    out.write_bytes((json.dumps(report, indent=2)+'\n').encode())
    print(json.dumps({'rows': rows, 'best_sampled_rms': best}), flush=True)


if __name__ == '__main__':
    run()
