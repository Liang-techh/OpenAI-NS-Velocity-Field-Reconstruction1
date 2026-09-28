"""Exact finite-interval profile-drift quadratic for two affine curl patches.

The fixed parent and the two patch supports define this diagnostic. It is
not a PDE certificate or a general criterion for recursive profiles.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from global_localized_candidate import load
from localized_constraint_rows import PATCH_CENTER, PATCH_WIDTHS
from localized_similarity_drift import cylindrical
from supported_fourier_basis import basis_data

ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'localized_drift_constraint.npz'


def response(points, center, widths, carrier):
    columns=[]
    for mode in (0,1,2):
        V,_,_=basis_data(points,center,widths,mode,2,mode*carrier)
        components=('real',) if mode==0 else ('real','imag')
        for index in range(27):
            for component in components:
                columns.append(V[:,:,index].real if component=='real' else -V[:,:,index].imag)
        columns.extend(np.zeros_like(points) for _ in range(9*len(components)))
    return np.stack(columns,axis=-1)


def build():
    _,parent,localized,snapshot,first,_,hashes=load()
    probe_path=ROOT/'localized_next_patch_probe.json'
    probe=json.loads(probe_path.read_text())
    patch2=probe['additional_patch']
    grid_path=ROOT/'enriched_endpoint_shape_cache.npz'
    with np.load(grid_path,allow_pickle=False) as cache:
        points,weights=cache['points'],cache['weights']
    h=float(localized.inner.h);tau0=snapshot['inputs']['mean']['tau']
    dk=1e-6;s=2**(-dk);dt=tau0*(1-s)
    mapped=points*np.array([s**.5,s**.5,s**(.5-h)])
    factors=np.array([s**.5,s**(.5+h),s**(.5+h)])
    u0=cylindrical(parent.fields(points,tau0)[0],points)
    u1=cylindrical(parent.fields(mapped,tau0*s)[0],mapped)*factors
    norm=np.sqrt(weights@np.sum(u0**2,axis=1))
    divisor=norm*abs(np.log(s))
    carrier=np.asarray(snapshot['inputs']['wave']['carrier'])
    V=np.concatenate((response(mapped,PATCH_CENTER,PATCH_WIDTHS,carrier),
                      response(mapped,patch2['center'],patch2['widths'],carrier)),axis=-1)
    theta=np.arctan2(mapped[:,1],mapped[:,0]);c=np.cos(theta)[:,None];sn=np.sin(theta)[:,None]
    cylindrical_V=np.stack((c*V[:,0]+sn*V[:,1],-sn*V[:,0]+c*V[:,1],V[:,2]),axis=1)
    weighted=np.repeat(np.sqrt(weights),3)
    A=(dt*cylindrical_V*factors[None,:,None]).reshape(-1,360)*weighted[:,None]/divisor
    b=(u1-u0).ravel()*weighted/divisor
    gram=A.T@A;linear=A.T@b;constant=float(b@b)
    control=np.r_[first['selected']['local_patch_coefficients'],np.zeros(180)]
    def square(x):return float(x@gram@x+2*linear@x+constant)
    scales=np.maximum(np.linalg.norm(A,axis=0),1e-30)
    direction=np.linspace(-1.,1.,360)/scales
    direction/=np.linalg.norm(A@direction)
    eps=1e-2
    numerical=(square(control+eps*direction)-square(control-eps*direction))/(2*eps)
    analytic=float(2*(gram@control+linear)@direction)
    direct=float(np.linalg.norm(A@control+b))
    np.savez(CACHE,gram=gram,linear=linear,constant=constant,
             delta_k=dk,reference_norm=norm,first_control=control)
    report=dict(status='completed',accepted=False,pde_validated=False,scale_recursion_established=False,
                source_hashes=hashes,probe_sha256=hashlib.sha256(probe_path.read_bytes()).hexdigest(),
                grid_sha256=hashlib.sha256(grid_path.read_bytes()).hexdigest(),
                cache_sha256=hashlib.sha256(CACHE.read_bytes()).hexdigest(),
                scope='7776-point finite-interval pulled-back profile drift per abs(log scale); exact quadratic in 360 physical tangent controls, fixed parent.',
                layout='first patch180 then second patch180; same interleaved real-minus-imag layout as localized fits',
                delta_k=dk,parent_drift=float(np.sqrt(constant)),first_patch_drift=float(np.sqrt(square(control))),
                first_patch_direct_drift=direct,quadratic_vs_direct_error=abs(np.sqrt(square(control))-direct),
                gradient_relative_error=abs(numerical-analytic)/max(abs(analytic),1e-30),
                constraint='c.T@gram@c + 2*linear@c + constant <= chosen_drift_cap**2',
                caveat='Stationarity diagnostic only; a dynamic or log-periodic recursive profile may have nonzero drift.')
    (ROOT/'localized_drift_constraint.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    return gram,linear,constant


if __name__=='__main__':build()
