"""Global two-patch diagnostic and direct similarity-profile drift replay."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from global_localized_candidate import load as load_first
from full_wave_tangent import LocalPotentialField, _unpack_full
from localized_similarity_drift import cylindrical

ROOT=Path(__file__).resolve().parent


def load(source=ROOT/'localized_two_patch_constrained.json'):
    raw=Path(source).read_bytes();candidate=json.loads(raw)
    if candidate['status']!='completed' or not candidate['assembled_feasible']:
        raise ValueError('Need a completed assembled-feasible two-patch candidate')
    _,field,localized,snapshot,_,_,hashes=load_first()
    for key in ('candidate','first_fit'):
        entry=candidate['sources'][key]
        if hashlib.sha256((ROOT/entry['path']).read_bytes()).hexdigest()!=entry['sha256']:
            raise ValueError(f'Changed {key} source')
    carrier=np.asarray(snapshot['inputs']['wave']['carrier'])
    original=snapshot['inputs']['wave']
    center,width=np.asarray(original['center']),np.asarray(original['widths'])
    for index,word in ((1,'first'),(2,'second')):
        c=np.asarray(candidate['inputs'][f'{word}_patch_center'])
        w=np.asarray(candidate['inputs'][f'{word}_patch_widths'])
        if not (np.all(c-w>center-width) and np.all(c+w<center+width)):
            raise ValueError('Patch escaped original wave support')
        derivative,pressure=_unpack_full(np.asarray(candidate['selected'][f'patch{index}_coefficients']),9)
        field=LocalPotentialField(field,c,w,2,{0:np.zeros(2),1:carrier,2:2*carrier},
                                  tuple(np.zeros(27,complex) for _ in range(3)),derivative,pressure,
                                  snapshot['inputs']['mean']['tau'])
    hashes[Path(source).name]=hashlib.sha256(raw).hexdigest()
    return field,localized,snapshot,candidate,hashes


def run():
    field,localized,snapshot,candidate,hashes=load()
    grid_path=ROOT/'enriched_endpoint_shape_cache.npz'
    with np.load(grid_path,allow_pickle=False) as cache:
        points,weights=cache['points'],cache['weights']
    tau=snapshot['inputs']['mean']['tau'];h=float(localized.inner.h)
    reference=cylindrical(field.fields(points,tau)[0],points)
    norm=np.sqrt(weights@np.sum(reference**2,axis=1))
    report=dict(status='running',accepted=False,pde_validated=False,scale_recursion_established=False,
                source_hashes=hashes,grid_sha256=hashlib.sha256(grid_path.read_bytes()).hexdigest(),
                scope='Direct assembled global velocity pullback on7776 points; profile stationarity diagnostic, not recursive closure.',rows=[])
    for dk in (1e-6,2e-6,1e-5):
        s=2**(-dk)
        mapped=points*np.array([s**.5,s**.5,s**(.5-h)])
        u=cylindrical(field.fields(mapped,tau*s)[0],mapped)*np.array([s**.5,s**(.5+h),s**(.5+h)])
        difference=u-reference
        relative=float(np.sqrt(weights@np.sum(difference**2,axis=1))/norm)
        row=dict(delta_k=dk,relative_profile_L2=relative,
                 drift_per_abs_log_scale=relative/abs(np.log(s)))
        if dk==1e-6:
            row['quadratic_prediction']=candidate['selected']['drift']
            row['prediction_difference']=row['drift_per_abs_log_scale']-row['quadratic_prediction']
        report['rows'].append(row);print(json.dumps(row),flush=True)
    report['status']='completed'
    (ROOT/'global_two_patch_candidate.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':run()
