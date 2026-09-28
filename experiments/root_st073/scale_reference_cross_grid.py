"""Compare cumulative local reference corrections on the original cached grid.

The grid was used to fit the parent, so this is a transfer diagnostic, not a
wholly unseen validation set. No PDE or recursive-step acceptance is implied.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import scale_reference_velocity_step as previous
import scale_reference_trust_fit as trust
from scale_reference_assembly_check import _metric

ROOT = Path(__file__).resolve().parent
PARENT = ROOT/'scale_reference_trust_nonlinear_fit.json'


def run(candidate_path, output_path):
    reports = [json.loads(path.read_text()) for path in (PARENT, candidate_path)]
    coefficients = []
    for report in reports:
        selected = report.get('selected')
        if selected is None:
            raise ValueError('No selected candidate; transfer is undefined')
        c = np.asarray(selected['coefficients'], float)
        if c.shape != (360,) or not np.all(np.isfinite(c)):
            raise ValueError('Expected finite cumulative 360-vector')
        coefficients.append(c)
    data = previous._load_inputs()
    baseline, _, pressure_design, _ = previous._pressure_baseline(data)
    patches, _, _, _ = trust._layout_and_patch_data(data, pressure_design)
    results = []
    for c in coefficients:
        linear = np.zeros(3*len(data['points']))
        for patch in range(2):
            block = c[180*patch:180*(patch+1)]
            linear += patches[patch][0] @ block[:135]
            linear += pressure_design[:, 45*patch:45*(patch+1)] @ block[135:]
        du, dj = previous._joint_velocity_from_coefficients(patches, trust._coefficients_by_patch(c))
        residual = baseline + linear.reshape(-1, 3) + np.einsum('nij,nj->ni', dj, du)
        results.append(dict(residual=_metric(residual, data['weights'], data['points']),
                            correction=_metric(du, data['weights'], data['points'])))
    parent, candidate = [row['residual'] for row in results]
    report = dict(status='completed', accepted=False, pde_validated=False,
                  scale_recursion_established=False,
                  sources={str(path.name): hashlib.sha256(path.read_bytes()).hexdigest()
                           for path in (PARENT, candidate_path, ROOT/'scale_generator_momentum_defect.npz')},
                  parent=results[0], candidate=results[1],
                  improves_L2=candidate['volume_L2'] < parent['volume_L2'],
                  improves_max=candidate['max_norm'] < parent['max_norm'],
                  L2_relative_change=candidate['volume_L2']/parent['volume_L2']-1,
                  max_relative_change=candidate['max_norm']/parent['max_norm']-1,
                  scope='Original 18720-point grid; analytic correction jets added to frozen base jets. Parent training grid, not a new independent holdout. Pressure includes prior projection once.')
    output_path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.candidate, args.output)
