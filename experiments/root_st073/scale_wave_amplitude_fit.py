"""Optimize original wave amplitude with pressure elimination and a 1% bound.

The mean velocity is fixed. The objective is the exact quadratic generator
residual on the frozen grid; the resulting scalar objective is quartic.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from supported_fourier_analytic_jets import basis_jets
from scale_generator_linearization import velocity_correction_response
from scale_reference_velocity_step import _load_inputs, _pressure_baseline, _metric

ROOT=Path(__file__).resolve().parent


def run():
    d=_load_inputs()
    points,weights=d['points'],d['weights']
    _,_,pressure_design,_=_pressure_baseline(d)
    geometry=d['frozen_report']['inputs']['wave']
    parent_path=ROOT/'enriched_mean_endpoint_tangent.json'
    parent=json.loads(parent_path.read_text())
    coeff=np.asarray(parent['selected']['coefficients_original'])
    coeff=coeff[:,0]+1j*coeff[:,1]
    v,j,visc,_,_=basis_jets(points,geometry['center'],geometry['widths'],1,2,geometry['carrier'],d['nu'])
    v=np.einsum('niq,q->ni',v,coeff).real
    j=np.einsum('nijq,q->nij',j,coeff).real
    visc=np.einsum('niq,q->ni',visc,coeff).real
    linear,quadratic=velocity_correction_response(points,d['velocity'],d['jacobian'],v,j,
        -visc/d['nu'],d['tau0'],d['h'],np.zeros(len(points)),d['nu'])
    row_weight=np.repeat(np.sqrt(weights),3)
    design=row_weight[:,None]*pressure_design
    scales=np.linalg.norm(design,axis=0)
    left,singular,right=np.linalg.svd(design/scales,full_matrices=False)
    keep=singular>singular[0]*1e-12
    q=left[:,keep]
    vectors=[row_weight*x.ravel() for x in (d['generator_residual'],linear,quadratic)]
    projected=[x-q@(q.T@x) for x in vectors]
    a,b,c=projected
    polynomial=np.array([a@a,2*a@b,b@b+2*a@c,2*b@c,c@c])
    bound=.01*np.sqrt(np.sum(weights[:,None]*d['velocity']**2)/np.sum(weights[:,None]*v**2))
    roots=np.polynomial.polynomial.polyroots(polynomial[1:]*np.arange(1,5))
    candidates=[-bound,0.,bound]+[float(x.real) for x in roots if abs(x.imag)<1e-10 and -bound<=x.real<=bound]
    shift=min(candidates,key=lambda x:np.polynomial.polynomial.polyval(x,polynomial))
    raw=d['generator_residual']+shift*linear+shift*shift*quadratic
    rhs=-row_weight*raw.ravel()
    normalized=right[keep].T@((q.T@rhs)/singular[keep])
    pressure_coeff=normalized/scales
    corrected=raw+(pressure_design@pressure_coeff).reshape(-1,3)
    baseline,_,_,_=_pressure_baseline(d)
    before=_metric(baseline,weights,points)
    after=_metric(corrected,weights,points)
    report=dict(status='completed',accepted=False,pde_validated=False,scale_recursion_established=False,
        sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            (parent_path,ROOT/'scale_generator_momentum_defect.npz',ROOT/'scale_generator_pressure_projection.json',
             ROOT/'localized_drift_1800_fit.json',ROOT/'full_wave_frozen_cache.json')},
        scope='Frozen inner-grid scalar amplitude plus pressure fit. Mean velocity fixed; initial mode1 wave scaled, not whole velocity. Pressure coefficients are additions to the original drift1800 pressure, not additions to pressure-projected baseline. Shape, moments, cones and independent replay unverified.',
        pressure_rank=int(keep.sum()),wave_amplitude=1+shift,amplitude_shift_bound=float(bound),
        velocity_change_fraction=float(abs(shift)/bound*.01),pressure_coefficients=pressure_coeff.tolist(),
        quartic_coefficients=polynomial.tolist(),candidate_shifts=candidates,
        baseline=before,selected=after,
        l2_relative_change=after['volume_L2']/before['volume_L2']-1,
        peak_relative_change=after['max_norm']/before['max_norm']-1)
    (ROOT/'scale_wave_amplitude_fit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ('wave_amplitude','velocity_change_fraction','l2_relative_change','peak_relative_change','selected')}))


if __name__=='__main__':run()
