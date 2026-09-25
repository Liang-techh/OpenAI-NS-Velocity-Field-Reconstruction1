"""Find the smallest tested M-step that keeps the physical cone open."""
import json

from extended_physical_cone_map import physical_cone
from high_frequency_shear_screen import make_field
from radial_continuation import ROOT
from radial_moment_step import RadialMomentStep


def run():
    tau = .5*2**(-5.5)
    base = make_field(16, 2.)
    rows = []
    for amplitude in (-.5, -.45, -.4, -.35, -.3, -.25):
        field = RadialMomentStep(base, amplitude)
        row = dict(moment=amplitude,
                   **physical_cone(field, 1., .2, tau, order=64))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(tau=tau, shear=2., rows=rows,
                  scope='One-point finite-amplitude physical cone refinement. Smaller M-step may reduce nonlinear outgoing moment defects but must pass neighboring eta, radius and time independently.',
                  accepted=False)
    (ROOT/'moment_step_amplitude_refine.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
