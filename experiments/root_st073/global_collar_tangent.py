"""Compact exterior tangent fit, preserving the original wave patch exactly."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[name]='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from global_axial_extension import build_candidate
from affine_momentum import jets, momentum
from supported_fourier_basis import basis_data
from wave_higher_harmonic_tangent import _mode_block_columns

ROOT=Path(__file__).resolve().parent
# Three disjoint interiors; all exclude the frozen compact wave support.
BOXES=[(.0001,.0049,.00050,.00094),(.0001,.0049,-.00094,-.00050),
       (.0049,.0105,-.001,.001)]


def geometry(box):
    r0,r1,z0,z1=box
    return np.array([(r0+r1)/2,(z0+z1)/2]),np.array([(r1-r0)/2,(z1-z0)/2])


def grid(order):
    x,w=np.polynomial.legendre.leggauss(order)
    points,weights=[],[]
    for box in BOXES:
        c,d=geometry(box)
        r,z=np.meshgrid(c[0]+d[0]*x,c[1]+d[1]*x,indexing='ij')
        points.append(np.column_stack((r.ravel(),np.zeros(r.size),z.ravel())))
        weights.append((w[:,None]*w[None,:]*d.prod()*2*np.pi*r).ravel())
    return np.concatenate(points),np.concatenate(weights)


def design(points):
    return np.column_stack([_mode_block_columns(points,*geometry(b),0,2,np.zeros(2))[:,::2]
                            for b in BOXES])


def metrics(residual,weights):
    return dict(momentum_max=float(np.max(np.linalg.norm(residual,axis=1))),
                momentum_volume_L2=float(np.sqrt(weights@np.sum(residual**2,axis=1))),
                physical_volume=float(weights.sum()))


class CollarCorrection:
    def __init__(self,base,control,tau0):
        self.base,self.control,self.tau0=base,np.asarray(control).reshape(3,36),tau0
        self.nu=base.nu

    def fields(self,points,tau):
        u,p=self.base.fields(points,tau)
        for box,c in zip(BOXES,self.control):
            V,P,_=basis_data(points,*geometry(box),0,2,np.zeros(2))
            u=u+(self.tau0-float(tau))*np.einsum('niq,q->ni',V,c[:27]).real
            p=p+np.einsum('nq,q->n',P,c[27:]).real
        return u,p


def run():
    full,localized,mean,candidate,snapshot,_=build_candidate()
    tau=snapshot['inputs']['mean']['tau']
    hs,ht=snapshot['timesteps']['hspace'],snapshot['timesteps']['htime']
    points,weights=grid(10)
    residual=momentum(jets(full,points,tau,hs,ht))
    D=design(points)
    sw=np.repeat(np.sqrt(weights),3)
    A=D*sw[:,None]
    scales=np.maximum(np.linalg.norm(A,axis=0),1e-30)
    q,_,rank,_=np.linalg.lstsq(A/scales,-residual.ravel()*sw,rcond=1e-10)
    control=q/scales
    report=dict(status='fitted',accepted=False,pde_validated=False,scale_recursion_established=False,
        scope='Reference-time exterior correction on three disjoint collar rectangles, not the full spatial domain. No temporal matching acceptance.',
        tau=tau,boxes=BOXES,control=control.tolist(),rank=int(rank),
        source_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                       for name in ('global_axial_extension.py','enriched_mean_endpoint_tangent.json')},
        training=dict(before=metrics(residual,weights),after=metrics(residual+(D@control).reshape(-1,3),weights)))
    path=ROOT/'global_collar_tangent.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['training']),flush=True)
    testpoints,testweights=grid(13)
    corrected=CollarCorrection(full,control,tau)
    before=momentum(jets(full,testpoints,tau,hs,ht))
    after=momentum(jets(corrected,testpoints,tau,hs,ht))
    analytic=before+(design(testpoints)@control).reshape(-1,3)
    report['independent_replay']=dict(before=metrics(before,testweights),after=metrics(after,testweights),
        analytic_vs_fd_max=float(np.max(np.linalg.norm(after-analytic,axis=1))))
    wave=snapshot['inputs']['wave'];c,w=np.asarray(wave['center']),np.asarray(wave['widths'])
    report['support_disjoint_from_wave']=all(b[1]<c[0]-w[0] or b[0]>c[0]+w[0] or
        b[3]<c[1]-w[1] or b[2]>c[1]+w[1] for b in BOXES)
    assert report['support_disjoint_from_wave']
    report['status']='completed'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['independent_replay']),flush=True)


if __name__=='__main__':
    run()
