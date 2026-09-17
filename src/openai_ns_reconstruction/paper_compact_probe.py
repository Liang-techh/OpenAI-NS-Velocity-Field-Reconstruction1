"""Export the compact velocity and retain independent full-NS failures."""
import json
from pathlib import Path
import numpy as np
from .paper_compact_field import PaperCompactField
from .constrained_force import RestrictedForce
from .constrained_validation import residual


def run():
    out=Path('artifacts/function_first/compact_field');out.mkdir(parents=True,exist_ok=True)
    f=PaperCompactField()
    force=RestrictedForce(a=9.546541832896667e-11,c=.12148700174414483)
    axis=np.linspace(-2,2,17);times=np.array([.25,.5,.75])
    points=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),axis=-1)
    grid=np.stack([f.at_points(points,t) for t in times])
    np.savez_compressed(out/'grid.npz',x=axis,y=axis,z=axis,times=times,
                        u=grid[...,0],v=grid[...,1],w=grid[...,2])
    rng=np.random.default_rng(926041)
    pts=np.vstack([rng.uniform(-1,1,(48,3)),
                   np.column_stack([np.linspace(.4,.85,16),np.zeros(16),np.linspace(0,.3,16)])])
    pressure=lambda p,t:np.array([f.pressure(*x,t) for x in p])
    rows=[]
    for t in times:
        for step in [.0025,.00125]:
            result=residual(f.at_points,pressure,force,pts,float(t),nu=.01,step=step)
            norms=np.linalg.norm(result['momentum'],axis=1)
            k=int(np.argmax(norms))
            row=dict(time=float(t),step=step,momentum_max=float(norms[k]),
                     momentum_peak=pts[k].tolist(),divergence_max=float(np.abs(result['divergence']).max()))
            rows.append(row);print(json.dumps(row),flush=True)
            (out/'full_residual.json').write_text(json.dumps(dict(
                status='development uniform-plus-targeted samples; not L2 quadrature',
                viscosity=.01,force=dict(a=force.a,c=force.c),points=pts.tolist(),rows=rows),indent=2)+'\n')
    (out/'metadata.json').write_text(json.dumps(dict(field=f.metadata(),
        grid_shape=list(grid.shape),sample_layout='time,x,y,z,component',
        momentum_acceptance='failed',force_source='unchanged coupled_joint parameters'),indent=2)+'\n')


if __name__=='__main__':
    run()
