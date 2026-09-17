"""Held-out azimuthal momentum obstruction for the normalized bipolar field."""
import json
from pathlib import Path
import numpy as np
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_force import RestrictedForce
from .constrained_validation import residual


def theta(points):
    radius=np.hypot(points[:,0], points[:,1])
    return np.column_stack((-points[:,1]/radius,points[:,0]/radius,np.zeros(len(points))))


def run():
    cfg=json.loads(Path('configs/constraints.json').read_text())
    field=Eq45SupportedDeliveryField.load_candidate('artifacts/bipolar_energy/normalized_candidate.json')
    zero=lambda x,t:np.zeros(len(x))
    forcezero=lambda x,t:np.zeros_like(x)
    def target(x,t,h):
        m=residual(field.at_points,zero,forcezero,x,t,nu=cfg['nu'],step=h)['momentum']
        return np.sum(m*theta(x),axis=1)
    rng=np.random.default_rng(9172603)
    train=rng.uniform(-2,2,(512,3)); directions=theta(train)
    a=[];b=[]
    for t in (.3125,.5,.6875):
        a.append(np.sum(RestrictedForce(a=0,c=1)(train,t)*directions,axis=1))
        b.append(target(train,t,.005))
    a=np.concatenate(a);b=np.concatenate(b)
    c=float(np.clip(a@b/(a@a),0,10))
    hold=np.random.default_rng(9172604).uniform(-2,2,(2048,3))
    direction=theta(hold)
    rows=[]
    for h in (.02,.01,.005):
        for t in cfg['validation']['times']:
            error=target(hold,t,h)-np.sum(RestrictedForce(a=0,c=c)(hold,t)*direction,axis=1)
            rows.append(dict(step=h,time=t,sampled_max=float(np.max(np.abs(error))),
                volume_L2_estimate=float(np.sqrt(64*np.mean(error**2)))))
    report=dict(candidate_sha256=field.sha256,nu=cfg['nu'],force_c=c,
        force_a='Unidentifiable from azimuthal component; radial/axial force has zero theta projection.',
        train_seed=9172603,holdout_seed=9172604,train_points=512,holdout_points=2048,
        sampling='Uniform Cartesian box [-2,2]^3',rows=rows,
        conclusion='Necessary axisymmetric momentum condition only. Pressure cannot change theta momentum; this does not construct pressure or validate the full PDE.',
        pde_validated=False)
    p=Path('artifacts/bipolar_theta/report.json');p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(force_c=c,finest=[r for r in rows if r['step']==.005]),indent=2))


if __name__=='__main__':run()
