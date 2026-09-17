"""Pressure-independent obstruction diagnostic for the axisymmetric candidate."""
import json
from pathlib import Path
import numpy as np
from .constrained_candidate import CompactCandidate
from .constrained_force import RestrictedForce
from .constrained_validation import residual


def run():
    root=Path('artifacts/constrained/optimized')
    c=CompactCandidate.load(root/'candidate.json')
    training=json.loads((root/'training.json').read_text())
    f=RestrictedForce(**training['force'])
    cfg=json.loads(Path('configs/constraints.json').read_text())
    rng=np.random.default_rng(cfg['validation']['seed'])
    x=rng.uniform(-2,2,(4096,3));r=np.hypot(x[:,0],x[:,1])
    et=np.column_stack((-x[:,1]/r,x[:,0]/r,np.zeros(len(x))))
    times=cfg['validation']['times'];rows=[]
    for t in times:
        for step in cfg['validation']['derivative_steps']:
            # Axisymmetric pressure contributes identically zero azimuthally.
            result=residual(c.velocity,lambda x,t:np.zeros(len(x)),f,x,t,step=step)
            theta=np.sum(result['momentum']*et,axis=1)
            rows.append({'time':t,'step':step,'azimuthal_sampled_max':float(np.max(np.abs(theta))),
                         'azimuthal_L2_estimate':float(np.sqrt(64*np.mean(theta**2)))})
    report={'scope':'fixed optimized velocity and force; all axisymmetric pressure choices',
            'interpretation':'Azimuthal momentum residual cannot be changed by axisymmetric pressure. This is a lower bound on full residual at these sampled points, not a bound for all velocity parameters.',
            'rows':rows}
    (root/'pressure_independent_residual.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([v for v in rows if v['step']==.005],indent=2))

if __name__=='__main__':run()
