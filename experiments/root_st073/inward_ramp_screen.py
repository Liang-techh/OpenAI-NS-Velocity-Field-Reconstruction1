"""Compare where the axis-regular heat swirl joins the ST073 core."""
import json
import numpy as np

from axial_compact_join import load_compact_candidate
from heated_interior_join import HeatedInteriorJoin
from joined_field import independent_fd
from radial_continuation import ROOT


def run():
    compact = load_compact_candidate()
    tau = .5*2**(-5.5)
    eta, radial_fraction = np.meshgrid(
        [.2725, .345, .4175], [.2, .4, .6, .8], indexing='ij')
    eta, radial_fraction = eta.ravel(), radial_fraction.ravel()
    q = tau/(1-eta**2)
    z = np.sqrt(compact.nu)*q**(.5-compact.joined.inner.h)*eta
    reference = HeatedInteriorJoin(compact=compact)
    r1 = float(reference.attachment_radius(tau))
    points = np.column_stack((r1*radial_fraction, np.zeros(len(eta)), z))
    rows = []
    for fraction in (1., .75, .5, .25):
        field = HeatedInteriorJoin(compact=compact, ramp_fraction=fraction)
        R, div = independent_fd(field, points, tau,
                                .001*np.sqrt(field.nu*tau), .00025*tau)
        norms = np.linalg.norm(R, axis=1)
        row = dict(ramp_fraction=fraction,
                   max_complete=float(np.max(norms)),
                   rms_complete=float(np.sqrt(np.mean(norms**2))),
                   max_angular=float(np.max(np.abs(R[:, 1]))),
                   divergence_max=float(np.max(np.abs(div))),
                   complete_norms=norms.tolist())
        rows.append(row)
        print(json.dumps({k: v for k, v in row.items()
                          if k != 'complete_norms'}), flush=True)
    report = dict(tau=tau, eta=eta.tolist(),
                  radial_fraction=radial_fraction.tolist(), rows=rows,
                  scope='Twelve inner/collar sample points at one time, not a whole-domain gate.',
                  accepted=False)
    (ROOT/'inward_ramp_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
