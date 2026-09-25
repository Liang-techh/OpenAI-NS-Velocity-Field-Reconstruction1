"""Screen a C2-preserving swirl shape mode in the actual radial bridge."""
import json

from extended_compact_join import load_extended_heated_candidate
from extended_relaxed_cone_screen import cone_point
from radial_continuation import ROOT


def run():
    tau = .5*2**(-5.5)
    rows = []
    for amplitude in (-1., -.5, 0., .5, 1.):
        field = load_extended_heated_candidate(amplitude)
        points = []
        for eta in (.2, .3):
            point = cone_point(field, 1., eta, tau, order=12)
            point['margin'] = (point['upper']-point['vs']
                               if point.get('upper') is not None else None)
            points.append(point)
        row = dict(amplitude=amplitude, points=points)
        rows.append(row)
        print(json.dumps(dict(amplitude=amplitude,
                              summary=[(p['eta'], p.get('a'), p.get('vs'),
                                        p.get('Pc'), p.get('margin'))
                                       for p in points])), flush=True)
    report = dict(tau=tau, X=1., rows=rows,
                  shape='delta u_theta = sqrt(nu) q^(-A) amplitude * 64*y^3*(1-y)^3, y=(r-ri)/(ro-ri). This vanishes with two radial derivatives at each bridge edge.',
                  scope='Two bridge snapshot cone points per amplitude. Full momentum and moment constraints are separate and not accepted.',
                  accepted=False)
    (ROOT/'swirl_bridge_cone_screen.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
