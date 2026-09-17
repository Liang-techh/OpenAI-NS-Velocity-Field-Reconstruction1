"""Bounded affine-in-time swirl experiment; no automatic acceptance."""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_bipolar_theta import theta
from .constrained_validation import residual
from .constrained_force import RestrictedForce

class TemporalSwirl:
    def __init__(self, base, slope):
        if not np.isfinite(slope) or not -1 <= slope <= 1:
            raise ValueError('slope must be in [-1,1]')
        self.base,self.slope=base,float(slope)
    def at_points(self, points, time):
        x=np.asarray(points,dtype=float); u=self.base.at_points(x,time)
        radius=np.hypot(x[...,0],x[...,1])
        direction=np.stack((-x[...,1],x[...,0],np.zeros_like(radius)),axis=-1)
        direction=np.divide(direction,radius[...,None],out=np.zeros_like(direction),where=radius[...,None]>0)
        swirl=np.sum(u*direction,axis=-1)[...,None]*direction
        return u+self.slope*(np.asarray(time)-.25)[...,None]*swirl

def run():
    cfg=json.loads(Path('configs/constraints.json').read_text())
    base=Eq45SupportedDeliveryField.load_candidate('artifacts/bipolar_energy/normalized_candidate.json')
    def target(field,x,t,h):
        v=residual(field.at_points,lambda x,t:np.zeros(len(x)),lambda x,t:np.zeros_like(x),x,t,nu=cfg['nu'],step=h)['momentum']
        return np.sum(v*theta(x),axis=1)
    x=np.random.default_rng(9172605).uniform(-2,2,(1024,3))
    matrices=[];targets=[]
    for t in (.3125,.5,.6875):
        b=target(base,x,t,.005)
        d=target(TemporalSwirl(base,1),x,t,.005)-b
        f=np.sum(RestrictedForce(a=0,c=1)(x,t)*theta(x),axis=1)
        matrices.append(np.column_stack((d,-f))); targets.append(-b)
    fit=lsq_linear(np.concatenate(matrices),np.concatenate(targets),bounds=([-1,0],[1,10]))
    k,c=map(float,fit.x);child=TemporalSwirl(base,k)
    x=np.random.default_rng(9172606).uniform(-2,2,(2048,3));rows=[]
    for h in (.01,.005):
        for t in cfg['validation']['times']:
            force=np.sum(RestrictedForce(a=0,c=c)(x,t)*theta(x),axis=1)
            for label,field in [('base',base),('child',child)]:
                err=target(field,x,t,h)-force
                rows.append(dict(candidate=label,step=h,time=t,sampled_max=float(np.max(np.abs(err))),L2_estimate=float(np.sqrt(64*np.mean(err**2)))))
    t=np.linspace(.25,.75,21);tau=1-t
    points=np.column_stack((.1*np.sqrt(tau),0*t,.1*tau**.495))
    u=child.at_points(points,t);scaled=u*np.column_stack((tau**.5,tau**.505,tau**.505))
    drift=float(np.max(np.linalg.norm(scaled-scaled[0],axis=1))/np.linalg.norm(scaled[0]))
    report=dict(parent_sha256=base.sha256,slope=k,force_c=c,train_seed=9172605,holdout_seed=9172606,rows=rows,core_drift=drift,core_drift_pass=drift<=.05,pde_validated=False,scope='Bounded affine swirl multiplier 1+k*(t-.25), k in [-1,1]. Initial energy unchanged; no zero-swirl collapse. Fixed force family. Comparator uses same fitted force. Diagnostic, no promotion.')
    p=Path('artifacts/bipolar_temporal/report.json');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))
    print(json.dumps([r for r in rows if r['time']==.75 and r['step']==.005],indent=2))
if __name__=='__main__':run()
