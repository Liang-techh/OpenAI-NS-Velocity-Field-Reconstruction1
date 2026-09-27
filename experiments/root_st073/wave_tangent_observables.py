"""Fixed-cylinder shape trends of a saved affine wave tangent, not a trajectory."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from broad_meridional_constrained import load_saved_field
from broad_shear_dynamic_control import load_saved_field as load_dynamic
from full_wave_tangent import LocalPotentialField, _unpack_full
from grouped_joined_field import install_in_field
from meridional_state_cache import value_and_pressure_modes
from vortex_state_observables import fixed_cylinder, observe

ROOT = Path(__file__).resolve().parent


def summarize_intervals(result):
    left, center, right = result['rows']
    keys = ('enstrophy_radial_rms','enstrophy_aspect_ratio',
            'enstrophy_weighted_angular_speed')
    result['one_sided_rates_per_k'] = {
        'backward': {key:(center[key]-left[key])/(center['k']-left['k']) for key in keys},
        'forward': {key:(right[key]-center[key])/(right['k']-center['k']) for key in keys}}
    result['forward_interval_direction_flags'] = dict(
        radial_contraction=right[keys[0]]<center[keys[0]],
        relative_axial_elongation=right[keys[1]]>center[keys[1]],
        angular_speed_magnitude_increase=abs(right[keys[2]])>abs(center[keys[2]]))
    result['interval_scope'] = ('One-sided changes expose curvature or sign reversal hidden by a central derivative. '
        'Even all passing signs are only an affine candidate diagnostic, not an NS time step.')
    return result


def run(candidate_path, output_path, delta_k=1e-6):
    started = time.perf_counter()
    raw = Path(candidate_path).read_bytes()
    candidate = json.loads(raw)
    if candidate['status'] != 'completed':
        raise ValueError('Candidate must be frozen before observation')
    snapshot = json.loads((ROOT/'full_wave_frozen_cache.json').read_text())
    geometry = snapshot['inputs']['wave']
    mean, report = load_saved_field()
    install_in_field(mean)
    if report['coefficients'] != snapshot['inputs']['mean']['coefficients']:
        raise ValueError('Mean differs from frozen candidate geometry')
    dynamic, dynamic_report = load_dynamic()
    install_in_field(dynamic)
    base, *_ = value_and_pressure_modes(dynamic, dynamic_report)
    k0 = dynamic_report['k']
    tau0 = snapshot['inputs']['mean']['tau']
    if not np.isclose(.5*2**(-k0), tau0, rtol=1e-12, atol=0):
        raise ValueError('Reference times differ')
    packed = np.array(candidate['selected']['coefficients_original'])
    wave = packed[:,0]+1j*packed[:,1]
    q = (geometry['degree']+1)**2
    derivatives, pressure = _unpack_full(candidate['selected']['tangent_coefficients'], q)
    # Observation uses velocity only; modal pressure can be omitted exactly.
    pressure = tuple(np.zeros_like(p) for p in pressure)
    carrier = np.array(geometry['carrier'])
    field = LocalPotentialField(mean, geometry['center'], geometry['widths'],
        geometry['degree'], {0:np.zeros(2),1:carrier,2:2*carrier},
        (np.zeros(3*q,complex),wave,np.zeros(3*q,complex)),
        derivatives, pressure, tau0)
    points, weights, domain = fixed_cylinder(base.inner, k0, angles=12)
    result = dict(status='running', accepted=False, pde_validated=False,
        scale_recursion_established=False, source=Path(candidate_path).name,
        source_sha256=hashlib.sha256(raw).hexdigest(), domain=domain,
        delta_k=delta_k, rows=[],
        scope='Fixed physical cylinder and affine tangent extension around one reference time. Sampled enstrophy shape, not an identified global core, resolved trajectory, or scale recursion.')
    output = Path(output_path)
    def save():
        output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    save()
    for offset in (-delta_k,0.,delta_k):
        row = observe(field,points,weights,k0+offset)
        result['rows'].append(row)
        save()
        print(json.dumps(dict(offset=offset, **row)),flush=True)
    left, center, right = result['rows']
    keys = ('enstrophy_radial_rms','enstrophy_aspect_ratio',
            'enstrophy_weighted_angular_speed','cylinder_enstrophy',
            'cylinder_kinetic_energy')
    slopes = {key:(right[key]-left[key])/(2*delta_k) for key in keys}
    result['derivative_per_k'] = slopes
    result['sampled_direction_flags'] = dict(
        radial_contraction=slopes['enstrophy_radial_rms']<0,
        relative_axial_elongation=slopes['enstrophy_aspect_ratio']>0,
        angular_speed_magnitude_increase=(
            center['enstrophy_weighted_angular_speed']*
            slopes['enstrophy_weighted_angular_speed']>0))
    result.update(status='completed', elapsed_seconds=time.perf_counter()-started)
    summarize_intervals(result)
    save()
    print(json.dumps(result['sampled_direction_flags']),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--candidate',type=Path,default=ROOT/'wave_moment_cone_tangent.json')
    parser.add_argument('--output',type=Path,default=ROOT/'wave_tangent_observables.json')
    parser.add_argument('--delta-k',type=float,default=1e-6)
    parser.add_argument('--summarize-existing',action='store_true')
    args = parser.parse_args()
    if args.summarize_existing:
        result = json.loads(args.output.read_text())
        if result['status'] != 'completed' or len(result['rows']) != 3:
            raise ValueError('Only a completed three-point report can be summarized')
        summarize_intervals(result)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(result['forward_interval_direction_flags']),flush=True)
    else:
        run(args.candidate,args.output,args.delta_k)
