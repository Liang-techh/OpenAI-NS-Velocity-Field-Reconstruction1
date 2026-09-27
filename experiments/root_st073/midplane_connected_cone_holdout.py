"""Independent spatial/edge cone screen for the joint mean repair."""
import argparse
import json
import numpy as np

from adaptive_bridge_recursive_defect import build_fields
from midplane_axial_cone_knob import cone_row
from midplane_physical_covariance_pairs import support as physical_support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run(edge=False):
    inner, fields = build_fields()
    name = 'midplane_connected_cone_edge_repair.json' if edge else 'midplane_connected_cone_repair.json'
    repaired = json.loads((ROOT / name).read_text())
    field = SeparatedMomentModes(fields['two_sided_cone'], repaired['amplitudes'],
                                 windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    supports = json.loads((ROOT / 'midplane_physical_covariance_pairs.json').read_text())
    offsets = (-.95, -.3, .3, .95)
    scales = []
    for k in ((11, 15, 19) if edge else (11, 19)):
        if k in (11, 19):
            source = json.loads((ROOT / 'compact_potential' /
                                 f'midplane_wave_source_k{k}.json').read_text())
            support = next(s['support'] for s in supports['scales'] if s['k'] == k)
        else:
            tau = .5 * 2.**-k
            X = inner.p.X_max * (1 + 15 * .325)**2
            source = dict(tau=tau, point=inner.from_similarity([X], [-.0125], tau)[0])
            support = physical_support(inner, source, field.nu)
        r0, _, z0 = source['point']
        rows = []
        for x in offsets:
            for z in offsets:
                r = r0 + x * support['radial_halfwidth']
                height = z0 + z * 1.4 * support['axial_halfwidth']
                cone = cone_row(field, float(r), float(height), source['tau'], z)
                rows.append(dict(radial_offset=x, axial_offset=z, **cone))
        scales.append(dict(k=k, rows=rows, passing=sum(r['cone_pass'] for r in rows)))
        print(f"k={k}: {scales[-1]['passing']}/16 pass; "
              f"max ratio={max(r.get('cone_ratio', float('inf')) for r in rows):.6g}", flush=True)
    report = dict(offsets=offsets, scales=scales,
                  scope='Independent 4x4 physical support nodes at reported scales. No continuous spatial, temporal, or interscale cone certificate.')
    name = 'midplane_connected_cone_edge_holdout.json' if edge else 'midplane_connected_cone_holdout.json'
    (ROOT / name).write_bytes((json.dumps(report, indent=2) + '\n').encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--edge', action='store_true')
    run(edge=parser.parse_args().edge)
