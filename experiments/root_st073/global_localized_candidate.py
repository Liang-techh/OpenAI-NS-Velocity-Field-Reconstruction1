"""Globally compact diagnostic field with the constrained inner patch.

load() returns a field whose fields(points, tau) yields velocity and pressure.
Physical time is t=-tau. This is not an accepted Navier--Stokes trajectory.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from global_axial_extension import build_candidate
from global_collar_tangent import CollarCorrection, BOXES
from global_support_energy import support_bounds
from enriched_shape_replay import build_field
from full_wave_tangent import LocalPotentialField, _unpack_full
from localized_constraint_rows import PATCH_CENTER, PATCH_WIDTHS

ROOT = Path(__file__).resolve().parent


def add_patch(base, candidate, snapshot):
    derivatives, pressures = _unpack_full(
        np.asarray(candidate['selected']['local_patch_coefficients']), 9)
    carrier = np.asarray(snapshot['inputs']['wave']['carrier'])
    return LocalPotentialField(base, PATCH_CENTER, PATCH_WIDTHS, 2,
                               {0: np.zeros(2), 1: carrier, 2: 2 * carrier},
                               tuple(np.zeros(27, complex) for _ in range(3)),
                               derivatives, pressures, snapshot['inputs']['mean']['tau'])


def load():
    path = ROOT / 'localized_constrained_tangent.json'
    raw = path.read_bytes()
    candidate = json.loads(raw)
    parent_path = ROOT / candidate['source']
    parent_raw = parent_path.read_bytes()
    if hashlib.sha256(parent_raw).hexdigest() != candidate['source_sha256']:
        raise ValueError('Localized correction parent changed')
    if candidate['status'] != 'completed' or not candidate['assembled_feasible']:
        raise ValueError('Need a frozen assembled-feasible local correction')
    base, localized, _, parent, snapshot, _ = build_candidate(parent_path)
    tau0 = snapshot['inputs']['mean']['tau']
    g = snapshot['inputs']['wave']
    c, w = np.asarray(g['center']), np.asarray(g['widths'])
    if not (np.all(PATCH_CENTER-PATCH_WIDTHS > c-w)
            and np.all(PATCH_CENTER+PATCH_WIDTHS < c+w)):
        raise ValueError('Local patch escaped original compact wave support')
    collar_path = ROOT / 'global_collar_tangent.json'
    collar_raw = collar_path.read_bytes()
    collar = json.loads(collar_raw)
    if collar['status'] != 'completed' or collar['tau'] != tau0:
        raise ValueError('Collar reference is incompatible')
    if not all(b[1] < c[0]-w[0] or b[0] > c[0]+w[0]
               or b[3] < c[1]-w[1] or b[2] > c[1]+w[1] for b in BOXES):
        raise ValueError('Collar intersects inner patch')
    exterior = CollarCorrection(base, collar['control'], tau0)
    # LocalPotentialField needs nu for downstream finite-difference momentum.
    exterior.nu = snapshot['inputs']['mean']['nu']
    full = add_patch(exterior, candidate, snapshot)
    hashes = {path.name: hashlib.sha256(raw).hexdigest(),
              parent_path.name: hashlib.sha256(parent_raw).hexdigest(),
              collar_path.name: hashlib.sha256(collar_raw).hexdigest()}
    return full, exterior, localized, snapshot, candidate, parent, hashes


def run():
    full, exterior, localized, snapshot, candidate, parent, hashes = load()
    inner_base, _ = build_field(parent)
    inner = add_patch(inner_base, candidate, snapshot)
    tau0 = snapshot['inputs']['mean']['tau']
    center = PATCH_CENTER
    inside = np.array([[center[0], 0., center[1]],
                       [0., center[0], center[1]],
                       [center[0], 0., center[1]+0.3*PATCH_WIDTHS[1]]])
    outer = np.array([[(b[0]+b[1])/2, 0., (b[2]+b[3])/2] for b in BOXES])
    report = dict(status='running', accepted=False, pde_validated=False,
                  scale_recursion_established=False, source_hashes=hashes,
                  scope='Global assembly and support probes; no global residual, interval, or recursion certification.',
                  rows=[])
    for dk in (0., 1e-6):
        tau = tau0 * 2**(-dk)
        u, p = full.fields(inside, tau)
        iu, ip = inner.fields(inside, tau)
        ou, op = full.fields(outer, tau)
        eu, ep = exterior.fields(outer, tau)
        bounds = support_bounds(localized, snapshot['inputs']['wave'], tau, extra_boxes=BOXES)
        outside = np.array([[bounds['radius']*1.01, 0., 0.],
                            [0., 0., bounds['z_lower']-1e-5],
                            [0., 0., bounds['z_upper']+1e-5]])
        zu, zp = full.fields(outside, tau)
        errors = [float(np.max(np.abs(a-b))) for a,b in ((u,iu),(p,ip),(ou,eu),(op,ep))]
        if any(e != 0. for e in errors) or np.any(zu != 0.) or np.any(zp != 0.):
            raise ValueError('Global composition changed inner/exterior field or support')
        report['rows'].append(dict(delta_k=dk, tau=tau, support_bounds=bounds,
                                   inner_velocity_error=errors[0], inner_pressure_error=errors[1],
                                   exterior_velocity_error=errors[2], exterior_pressure_error=errors[3],
                                   outside_velocity_max=float(np.abs(zu).max()),
                                   outside_pressure_max=float(np.abs(zp).max())))
    report['status'] = 'completed'
    (ROOT/'global_localized_candidate.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
