"""Same-sample momentum screen of energy-optimized baseline and odd extension."""
import json
from pathlib import Path
import numpy as np
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_validation import residual
from .constrained_force import RestrictedForce
from .constrained_bipolar_theta import theta

def run(names=("bipolar_joint_baseline96","bipolar_joint_odd3"), seed=9172610, output="artifacts/bipolar_joint_comparison/report.json"):
    x=np.random.default_rng(seed).uniform(-2,2,(2048,3));rows=[];info={}
    for name in names:
        root=Path('artifacts')/name;r=json.loads((root/'report.json').read_text())
        f=Eq45SupportedDeliveryField.load_candidate(root/'candidate.json');force=RestrictedForce(*r['parameters'][-2:])
        info[name]={'sha256':f.sha256,'energy_balance_holdout':r['holdout'],'force':r['parameters'][-2:]}
        for h in (.01,.005):
            for t in (.3125,.5,.6875,.75):
                result=residual(f.at_points,lambda x,t:np.zeros(len(x)),force,x,t,step=h)
                m=result['momentum'];a=np.sum(m*theta(x),axis=1)
                rows.append(dict(candidate=name,time=t,step=h,theta_max=float(np.max(np.abs(a))),theta_L2=float(np.sqrt(64*np.mean(a*a))),zero_pressure_momentum_max=float(np.max(np.linalg.norm(m,axis=1))),divergence_max=float(np.max(np.abs(result['divergence'])))))
    report=dict(seed=seed,points=2048,candidates=info,rows=rows,pde_validated=False,scope='Fresh uniform box sample. Theta residual is pressure-independent for axisymmetric pressure. Full vector uses p=0 and is not minimized over pressure; do not call it a pressure-independent obstruction.')
    p=Path(output);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([r for r in rows if r['step']==.005],indent=2))
if __name__=='__main__':run()
