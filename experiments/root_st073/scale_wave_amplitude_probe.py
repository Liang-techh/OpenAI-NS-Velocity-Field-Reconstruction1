"""Separate the original oscillatory amplitude from the unchanged mean flow.

An exploratory direction only: pressure is frozen and physical constraints
are not accepted. This scales the original mode-1 wave, not the entire field.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key] = '1'
import hashlib
import json
from pathlib import Path
import numpy as np
from supported_fourier_analytic_jets import basis_jets
from scale_generator_linearization import velocity_correction_response

ROOT = Path(__file__).resolve().parent


def run():
    cp = ROOT/'scale_generator_momentum_defect.npz'
    parent_path = ROOT/'enriched_mean_endpoint_tangent.json'
    geometry_path = ROOT/'full_wave_frozen_cache.json'
    parent = json.loads(parent_path.read_text())
    geometry = json.loads(geometry_path.read_text())['inputs']['wave']
    inputs = json.loads((ROOT/'scale_generator_momentum_defect.json').read_text())['inputs']
    with np.load(cp) as d:
        points, weights, u, j, residual = [d[k] for k in ('points','weights','velocity','gradient','generator_residual')]
    c = np.asarray(parent['selected']['coefficients_original'])
    c = c[:,0]+1j*c[:,1]
    v, dj, visc, _, _ = basis_jets(points,geometry['center'],geometry['widths'],1,geometry['degree'],geometry['carrier'],inputs['viscosity'])
    v = np.einsum('niq,q->ni',v,c).real
    dj = np.einsum('nijq,q->nij',dj,c).real
    visc = np.einsum('niq,q->ni',visc,c).real
    linear, quadratic = velocity_correction_response(
        points,u,j,v,dj,-visc/inputs['viscosity'],inputs['tau0'],
        inputs['similarity_exponent_h'],np.zeros(len(points)),inputs['viscosity'])
    def norm(x):
        return float(np.sqrt(np.sum(weights[:,None]*x*x)))
    rows=[]
    for amplitude in (1,.99,.95,.9,.75,.5,.25,0):
        a=amplitude-1
        r=residual+a*linear+a*a*quadratic
        rows.append(dict(wave_amplitude=amplitude,momentum_volume_L2=norm(r),
                         momentum_max=float(np.linalg.norm(r,axis=1).max()),
                         velocity_change_fraction=norm(a*v)/norm(u),
                         velocity_volume_L2=norm(u+a*v)))
    report=dict(status='completed',accepted=False,pde_validated=False,scale_recursion_established=False,
                sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (cp,parent_path,geometry_path)},
                scope='Reference-time analytic wave amplitude direction on frozen inner grid; original mean and pressure fixed. Mode1 mean swirl is zero. Endpoint shape, moment, cone, independent replay and trajectory unverified. Zero-amplitude row is a diagnostic limit, not a candidate.',
                reference_velocity_L2=norm(u),original_wave_velocity_L2=norm(v),rows=rows)
    (ROOT/'scale_wave_amplitude_probe.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':
    run()
