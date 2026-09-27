"""Replay a saved local k-affine correction away from its fitted instant.

This is an audit of the explicit callable, not an integrated NS trajectory.
The unforced residual and volume integral use the same similarity annulus
at every time. No refitting or pressure substitution is performed.
"""
import json
import time
from pathlib import Path

import numpy as np

from broad_meridional_momentum import load_saved_field, quadrature_nodes
from broad_meridional_momentum import evaluate_field, RADIAL_QUADRATURE_BREAKS
from broad_shear_dynamic_control import load_saved_field as load_dynamic_field


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / 'broad_meridional_time_audit.json'


def run(deltas=(0., 1.e-5, 1.e-4, 1.e-3), output=OUTPUT):
    started = time.perf_counter()
    candidate, source = load_saved_field()
    reference, _ = load_dynamic_field()
    k0 = source['k']
    report = dict(
        accepted=False, pde_validated=False, scale_recursion_established=False,
        source='broad_meridional_momentum.json',
        status='running', nu=candidate.nu, k0=k0,
        physical_time='t=-tau, tau=0.5*2**(-k)', forcing='zero',
        domain=dict(eta=[-.5,.5], radial_fraction=[.01,.99], theta=[0.,2*np.pi]),
        quadrature=dict(eta_order=4, radial_order_per_interval=4,
                        radial_breaks=list(RADIAL_QUADRATURE_BREAKS),
                        angle=np.pi/7),
        scope='Off-time replay of the saved affine-k callable without refitting. Not numerical time integration, not constraint preservation, and not a recursion claim.',
        requested_delta_k=list(deltas), rows=[])

    def save():
        report['elapsed_seconds'] = time.perf_counter()-started
        Path(output).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')

    save()
    for delta in deltas:
        k = k0 + delta
        tau = .5*2.**(-k)
        points, weights = quadrature_nodes(candidate.inner,tau,4,4,1,
            radial_bounds=(.01,.99),radial_breaks=RADIAL_QUADRATURE_BREAKS,
            angle_shift=np.pi/7)
        row = dict(k=k, delta_k=delta, tau=tau,
                   elapsed_physical_time=.5*2.**(-k0)-tau,
                   physical_volume=float(weights.sum()))
        for name, field in (('reference_dynamic',reference),('corrected',candidate)):
            metric = evaluate_field(field,points,tau,weights)
            metric.pop('residual_components')
            velocity = field.fields(points,tau)[0]
            metric['velocity_sample_max'] = float(np.max(np.linalg.norm(velocity,axis=1)))
            metric['sampled_domain_kinetic_energy'] = float(.5*np.sum(weights[:,None]*velocity**2))
            row[name] = metric
        row['max_ratio_to_reference'] = row['corrected']['momentum_max']/row['reference_dynamic']['momentum_max']
        row['L2_ratio_to_reference'] = row['corrected']['momentum_volume_L2']/row['reference_dynamic']['momentum_volume_L2']
        report['rows'].append(row)
        save()
        print(json.dumps(row),flush=True)
    report['status']='completed'
    report['all_sampled_times_improve_both_metrics'] = all(
        row['max_ratio_to_reference']<1 and row['L2_ratio_to_reference']<1 for row in report['rows'])
    save()
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--delta-k',type=float,nargs='+',default=[0.,1.e-5,1.e-4,1.e-3])
    parser.add_argument('--output',type=Path,default=OUTPUT)
    args=parser.parse_args()
    run(args.delta_k,args.output)
