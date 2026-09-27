"""Check radial quadrature sensitivity of the candidate physical moments."""
import argparse
import json
import numpy as np

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run(split=False):
    inner, fields = build_fields()
    data = json.loads((ROOT / 'midplane_connected_cone_edge_repair.json').read_text())
    field = SeparatedMomentModes(fields['two_sided_cone'], data['amplitudes'],
                                 windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    rows = []
    breaks = sorted({v for lo, hi in RADIAL_WINDOWS_THREE for v in (lo, (lo + hi) / 2, hi)}) if split else None
    for order in ((12, 24, 48) if split else (12, 24, 48, 96)):
        slices = moment_slices(inner, fields['two_sided_cone'], field,
                                orders=(11, 15, 19), n=order, unit_fields=[field], radial_breaks=breaks)
        values = [outer_moments(s, np.zeros(1)).tolist() for s in slices]
        rows.append(dict(gauss_order=order, moments=values))
        print(f"order={order}: maxima={[max(abs(np.array(v))) for v in values]}", flush=True)
    report = dict(scales=[11, 15, 19], radial_breaks=breaks, rows=rows,
                  scope='Direct physical moment quadrature sensitivity of a fitted mean. No momentum-norm or recursive acceptance.')
    name = 'midplane_connected_moment_quadrature_split.json' if split else 'midplane_connected_moment_quadrature.json'
    (ROOT / name).write_bytes((json.dumps(report, indent=2) + '\n').encode())


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--split', action='store_true')
    run(split=parser.parse_args().split)
