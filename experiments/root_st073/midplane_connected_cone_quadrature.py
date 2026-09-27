"""Independent radial-stress quadrature sensitivity of the cone candidate."""
import json
from adaptive_bridge_recursive_defect import build_fields
from midplane_axial_cone_knob import cone_row
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


def run():
    _, fields = build_fields()
    data = json.loads((ROOT / 'midplane_connected_cone_edge_repair.json').read_text())
    field = SeparatedMomentModes(fields['two_sided_cone'], data['amplitudes'],
                                 windows=RADIAL_WINDOWS_THREE, knots=(11., 15., 19.))
    supports = json.loads((ROOT / 'midplane_physical_covariance_pairs.json').read_text())
    results = []
    for k in (11, 19):
        source = json.loads((ROOT / 'compact_potential' / f'midplane_wave_source_k{k}.json').read_text())
        support = next(s['support'] for s in supports['scales'] if s['k'] == k)
        r0, _, z0 = source['point']
        for offset in (-.99, 0.):
            r = r0 + offset * support['radial_halfwidth']
            z = z0 + offset * 1.4 * support['axial_halfwidth']
            rows = [dict(order=order, **cone_row(field, float(r), float(z), source['tau'], offset, stress_order=order))
                    for order in (12, 24, 48, 96)]
            results.append(dict(k=k, diagonal_offset=offset, rows=rows))
            print(k, offset, [row.get('cone_ratio') for row in rows], flush=True)
    (ROOT / 'midplane_connected_cone_quadrature.json').write_bytes((json.dumps(dict(results=results), indent=2) + '\n').encode())


if __name__ == '__main__':
    run()
