"""Independent Cartesian energy quadrature and fixed-candidate derivative study."""
import json
from pathlib import Path
import numpy as np
from .quadrature import unit_rule
from .constrained_candidate import CompactCandidate
from .constrained_force import RestrictedForce
from .constrained_validation import residual,sampled_norms


def cartesian_energy(velocity,time,order):
    nodes,weights=unit_rule(order)
    grid=4*nodes-2;weights=4*weights
    # Slice by z to bound memory; this is not the candidate's cylindrical rule.
    xx,yy=np.meshgrid(grid,grid,indexing='ij')
    xy=np.column_stack((xx.ravel(),yy.ravel()))
    wxy=(weights[:,None]*weights[None,:]).ravel()
    total=0.0
    for z,w in zip(grid,weights):
        points=np.column_stack((xy,np.full(len(xy),z)))
        total+=w*np.dot(wxy,np.sum(velocity(points,time)**2,axis=1))/2
    return float(total)


def run():
    folder=Path('artifacts/constrained/optimized_v4')
    c=CompactCandidate.load(folder/'candidate.json')
    training=json.loads((folder/'training.json').read_text())
    f=RestrictedForce(**training['force'])
    cfg=json.loads(Path('configs/constraints_v4.json').read_text())
    rows=[]
    for t in (.25,.5,.75):
        for order in (24,48,96):
            rows.append({'time':t,'order':order,'cartesian_energy':cartesian_energy(c.velocity,t,order)})
    x=np.random.default_rng(914027).uniform(-2,2,(4096,3))
    derivatives=[]
    for step in (.02,.01,.005,.0025,.00125):
        derivatives.append({'time':.75,'step':step,**sampled_norms(
            residual(c.velocity,c.pressure,f,x,.75,step=step),64)})
    boundary=np.array([[2,0,0],[-2,0,0],[0,2,0],[0,-2,0],[0,0,2],[0,0,-2]])
    report={'candidate':'optimized_v4/candidate.json','scope':'fixed candidate; Cartesian full-box quadrature independent of cylindrical normalization; derivative study changes step only',
        'energy':rows,'derivatives':derivatives,
        'boundary_probe_max':float(np.max(np.abs(c.velocity(boundary,.75)))),
        'limits':'Boundary probes alone are not a support proof. Support follows the documented compact bump representation. No PDE acceptance or stability claim.'}
    (folder/'convergence.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':run()
