"""Actual finite-difference replay of the constrained localized candidate."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from affine_momentum import jets, momentum
from enriched_shape_replay import build_field
from full_wave_tangent import LocalPotentialField, _unpack_full
from localized_constraint_rows import PATCH_CENTER, PATCH_WIDTHS

ROOT = Path(__file__).resolve().parent


def metric(residual, weights):
    square = np.sum(residual**2, axis=1)
    return dict(maximum=float(np.sqrt(square.max())),
                volume_L2=float(np.sqrt(weights @ square)))


def run():
    started = time.perf_counter()
    candidate_path = ROOT / 'localized_constrained_tangent.json'
    candidate_bytes = candidate_path.read_bytes()
    candidate = json.loads(candidate_bytes)
    if candidate['status'] != 'completed' or not candidate['assembled_feasible']:
        raise ValueError('Need a completed feasible localized candidate')
    parent_path = ROOT / candidate['source']
    if hashlib.sha256(parent_path.read_bytes()).hexdigest() != candidate['source_sha256']:
        raise ValueError('Localized candidate parent changed')
    parent, snapshot = build_field(json.loads(parent_path.read_text()))
    tau0 = snapshot['inputs']['mean']['tau']
    carrier = np.asarray(snapshot['inputs']['wave']['carrier'])
    control = np.asarray(candidate['selected']['local_patch_coefficients'])
    derivatives, pressures = _unpack_full(control, 9)
    field = LocalPotentialField(parent, PATCH_CENTER, PATCH_WIDTHS, 2,
                                {0: np.zeros(2), 1: carrier, 2: 2 * carrier},
                                tuple(np.zeros(27, complex) for _ in range(3)),
                                derivatives, pressures, tau0)
    grid_path = ROOT / 'full_wave_dense_tangent.json'
    grid = json.loads(grid_path.read_text())['new_frozen_cache']
    points, weights = np.asarray(grid['points']), np.asarray(grid['weights'])
    output = ROOT / 'localized_actual_replay.json'
    report = dict(status='running', accepted=False, pde_validated=False,
                  scale_recursion_established=False,
                  scope='Independent 18720-point frozen grid; actual Cartesian finite differences at two times. No continuum or time-interval certificate.',
                  source=candidate_path.name,
                  source_sha256=hashlib.sha256(candidate_bytes).hexdigest(),
                  parent_sha256=candidate['source_sha256'],
                  grid_sha256=hashlib.sha256(grid_path.read_bytes()).hexdigest(),
                  point_count=len(points), timesteps=snapshot['timesteps'], rows=[])
    def save():
        output.write_text(json.dumps(report, indent=2) + '\n')
    save()
    for delta_k in (0.0, 1e-6):
        tau = tau0 * 2**(-delta_k)
        row = dict(delta_k=delta_k, tau=tau, physical_time_increment=tau0-tau)
        report['rows'].append(row)
        for name, target in (('parent', parent), ('corrected', field)):
            print(json.dumps(dict(stage='running', delta_k=delta_k, field=name)), flush=True)
            residual = momentum(jets(target, points, tau, snapshot['timesteps']['hspace'],
                                     snapshot['timesteps']['htime']))
            row[name] = metric(residual, weights)
            save()
            print(json.dumps(dict(stage='replayed', delta_k=delta_k, field=name,
                                  metrics=row[name])), flush=True)
        row['both_metrics_improve'] = all(row['corrected'][key] < row['parent'][key]
                                          for key in ('maximum', 'volume_L2'))
    report.update(status='completed', elapsed_seconds=time.perf_counter()-started)
    save()
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    run()
