"""Shape observables of the corrected prescribed scale family.

The diagnostic cylinder is pushed forward with the scale map. This checks
kinematic morphology, not NS evolution or a fixed-domain core certificate.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from scale_reference_candidate import load_reference
from scale_transport_generator import _finite_difference_jacobian

ROOT=Path(__file__).resolve().parent


def run(source=ROOT/'scale_reference_velocity_step.json',output=ROOT/'scale_reference_shape_replay.json'):
    source=Path(source)
    reference=load_reference(step=source)
    path=ROOT/'enriched_endpoint_shape_cache.npz'
    with np.load(path) as d:
        points,weights=d['points'],d['weights']
    n=len(points)
    theta=np.arctan2(points[:,1],points[:,0])
    radius=np.linalg.norm(points[:,:2],axis=1)
    if n%12 or np.any(radius<=0):
        raise ValueError('Expected positive-radius twelve-angle rings')
    rings=points.reshape(-1,12,3)
    if not np.allclose(rings[:,:,2],rings[:,:1,2]) or not np.allclose(radius.reshape(-1,12),radius.reshape(-1,12)[:,:1]):
        raise ValueError('Invalid angular ring geometry')
    angles=theta.reshape(-1,12)
    if not np.allclose(np.exp(1j*(angles-angles[:,:1])),np.exp(2j*np.pi*np.arange(12)/12),atol=1e-10):
        raise ValueError('Invalid angular spacing')
    generator_report=json.loads((ROOT/'scale_generator_momentum_defect.json').read_text())
    h=generator_report['inputs']['similarity_exponent_h']
    step=generator_report['inputs']['hspace']
    u=reference.velocity(points)
    j=_finite_difference_jacobian(reference,points,step)
    er=np.column_stack((np.cos(theta),np.sin(theta),np.zeros(n)))
    et=np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(n)))
    def mean(v):return np.repeat(v.reshape(-1,12).mean(axis=1),12)
    m=mean(np.sum(et*u,axis=1))
    mr=mean(np.einsum('ni,nij,nj->n',et,j,er))
    mz=mean(np.einsum('ni,ni->n',et,j[:,:,2]))
    rows=[]
    for s in (1,.5,.25):
        space=np.array([s**.5,s**.5,s**(.5-h)])
        amp=np.array([s**-.5,s**-.5,s**(-.5-h)])
        c=amp[2]-amp[0]
        x=points*space
        velocity=u*amp+c*m[:,None]*et
        dm=mr[:,None]*er/space
        dm[:,2]=mz/space[2]
        gradient=j*amp[None,:,None]/space[None,None,:]
        gradient+=c*(et[:,:,None]*dm[:,None,:]-(m/(radius*space[0]))[:,None,None]*er[:,:,None]*et[:,None,:])
        omega=np.column_stack((gradient[:,2,1]-gradient[:,1,2],gradient[:,0,2]-gradient[:,2,0],gradient[:,1,0]-gradient[:,0,1]))
        w=weights*np.prod(space)
        density=np.sum(omega**2,axis=1)
        total=w@density
        r=radius*space[0]
        zmean=w@(density*x[:,2])/total
        radial=np.sqrt(w@(density*r*r)/total)
        axial=np.sqrt(w@(density*(x[:,2]-zmean)**2)/total)
        rows.append(dict(scale=s,radial_rms=float(radial),axial_rms=float(axial),
                         aspect_ratio=float(axial/radial),
                         signed_angular_speed=float(w@(density*np.sum(velocity*et,axis=1)/r)/total),
                         sampled_kinetic_energy=float(.5*np.sum(w[:,None]*velocity**2)),
                         divergence_max=float(abs(np.trace(gradient,axis1=1,axis2=2)).max())))
    report=dict(status='completed',accepted=False,pde_validated=False,scale_recursion_established=False,
                source=source.name,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                grid_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),point_count=n,spatial_fd_step=step,
                scope=__doc__,rows=rows)
    Path(output).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(rows),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=ROOT/'scale_reference_velocity_step.json')
    parser.add_argument('--output',type=Path,default=ROOT/'scale_reference_shape_replay.json')
    args=parser.parse_args()
    run(args.source,args.output)
