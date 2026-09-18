"""ST047: coordinate-preconditioned continuation; physical field/force unchanged."""
from __future__ import annotations
import argparse, hashlib, json, time, sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2]
for name in ('root_st030','root_st040','root_st043'):
    sys.path.insert(0,str(ROOT/'experiments'/name))
from coupled_fit import JointObjective
from localized_model import LocalizedModel

class Objective(JointObjective):
    def __init__(self,m,acceleration_ratio=.98,axis_weight=.02,**kw):
        super().__init__(m,**kw)
        r,_=m.momentum(np.zeros(m.dim),self.exp)
        self.parent_mz=self.signz*(r[:,2]-self.exp['v']['Qz'])
        self.accscale=max(float(np.sqrt(np.mean(self.parent_mz**2))),.01)
        self.acc_target=self.parent_mz-(1-acceleration_ratio)*np.maximum(self.parent_mz,0)
        R,Z,T=np.meshgrid([0.,.01,.04,.12,.2],[-.2,-.14,-.08,-.02,0.,.02,.08,.14,.2],np.linspace(.25,.75,17),indexing='ij')
        self.axis=m.cache(R*R*(1-T),Z*(1-T)**.495,T)
        self.axis_weight=axis_weight
        self.c_last=None;self.c_value=None;self.f_last=None;self.f_value=None
    def constraints(self,c):
        if self.c_last is not None and np.array_equal(c,self.c_last):return self.c_value
        v,j=super().constraints(c);r,k=self.m.momentum(c,self.exp)
        sl=slice(self.m.nv,self.m.nv+self.m.np)
        pg=self.exp['v']['Qz']+self.exp['M']['Qz']@c[sl]
        mz=self.signz*(r[:,2]-pg)
        km=self.signz[:,None]*k[:,2].copy();km[:,sl]-=self.signz[:,None]*self.exp['M']['Qz']
        self.c_last=c.copy();self.c_value=(np.r_[v,(self.acc_target-mz)/self.accscale],np.vstack((j,-km/self.accscale)))
        return self.c_value
    def fun(self,c):
        if self.f_last is not None and np.array_equal(c,self.f_last):return self.f_value
        value,g=super().fun(c);r,j=self.m.momentum(c,self.axis)
        sq=np.sum(r*r,axis=1);temperature=.0003;ex=sq/temperature;weights=np.exp(ex-logsumexp(ex))
        value+=self.axis_weight*temperature*(logsumexp(ex)-np.log(len(ex)))+.01*float(np.mean(sq))
        g=g+2*np.einsum('n,ni,nik->k',self.axis_weight*weights+.01/len(ex),r,j)
        self.f_last=c.copy();self.f_value=(value,g)
        return self.f_value

def precondition(obj,c,method='whiten',factor=1e4):
    _,J=obj.cached_momentum(c)
    A=(J*np.sqrt(obj.weights)[:,None,None]).reshape(-1,len(c));gram=A.T@A
    d=np.sqrt(np.maximum(np.diag(gram),1e-20));standard=gram/d[:,None]/d[None,:]
    vals,V=eigh(standard,check_finite=False)
    if method=='none':T=np.eye(len(c))
    elif method=='diagonal':T=np.diag(1/np.maximum(np.sqrt(2*factor)*d,1e-6))
    elif method=='whiten':T=(V/np.sqrt(2*factor*np.maximum(vals,1e-3))[None,:])/d[:,None]
    else:raise ValueError('Unknown preconditioner')
    _,K=obj.constraints(c)
    rs=np.clip(1/np.maximum(np.linalg.norm(K@T,axis=1),1e-4),.01,1e4)
    info=dict(method=method,objective_scale=factor,standardized_gram_eigen_min=float(vals[0]),standardized_gram_eigen_max=float(vals[-1]),standardized_gram_rank_rel1e10=int(np.sum(vals>vals[-1]*1e-10)),coordinate_condition=float(np.linalg.cond(T)),coordinate_only=True)
    return T,rs,info

def run(parent,out,method='whiten',maxiter=180,seconds=480.,radial=4,axial=4,td=4,pressure=.99,profile=.995,acc=.98,axis_weight=.02):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Refuse nonempty output')
    config=dict(experiment='ST047',parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),method=method,maxiter=maxiter,wall_seconds=seconds,radial=radial,axial=axial,time_degree=td,pressure_ratio=pressure,profile_ratio=profile,acceleration_ratio=acc,axis_weight=axis_weight,training_space=[32,48],training_time='13Gauss+endpoints',new_validation_seeds=[9174701,9174702],original_gates='unchanged',normalized_training_feasibility_tolerance=1e-7,registered_before_fit=True)
    (out/'registration.json').write_text(json.dumps(config,indent=2))
    start=time.monotonic();m=LocalizedModel(parent,radial,axial,td)
    obj=Objective(m,acc,axis_weight,profile_ratio=profile,pressure_ratio=pressure,space_order=(32,48),peak_weight=0.,softmax_weight=.02,anchor_tolerance=.08,time_cap=1.)
    zero=np.zeros(m.dim);T,rs,info=precondition(obj,zero,method)
    (out/'conditioning.json').write_text(json.dumps(info,indent=2));np.save(out/'coordinate_transform.npy',T)
    np.savez_compressed(out/'model.npz',Ma=m.Ma,Mb=m.Mb,Mp=m.Mp);np.save(out/'row_scale.npy',rs)
    print(json.dumps(info),flush=True)
    factor=1e4;lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    def fun(y):
        v,g=obj.fun(T@y);return v*factor,(T.T@g)*factor
    last_c=[None,None]
    def constraints(y):
        if last_c[0] is not None and np.array_equal(y,last_c[0]):return last_c[1]
        c=T@y;v,j=obj.constraints(c)
        last_c[:]=[y.copy(),(np.r_[v*rs,(c-lo)*bs,(hi-c)*bs],np.vstack(((j@T)*rs[:,None],T*bs[:,None],-T*bs[:,None])))]
        return last_c[1]
    rng=np.random.default_rng(9174790);direction=rng.normal(size=m.dim);direction/=np.linalg.norm(direction)
    h=2e-6;fv,g=fun(zero);cv,cj=constraints(zero)
    cp=constraints(h*direction)[0];cm=constraints(-h*direction)[0]
    check=dict(objective_abs_error=abs((fun(h*direction)[0]-fun(-h*direction)[0])/(2*h)-g@direction),constraint_abs_error=float(np.max(abs((cp-cm)/(2*h)-cj@direction))))
    (out/'derivative_checks.json').write_text(json.dumps(check,indent=2));assert check['constraint_abs_error']<1e-4,check
    print(json.dumps(check),flush=True)
    history=[];best=[np.inf,None];last=[zero.copy()]
    class BudgetStop(Exception):pass
    def callback(y):
        last[0]=y.copy();c=T@y;loss=obj.fun(c)[0];feas=float(obj.constraints(c)[0].min());it=len(history)+1
        inbounds=bool(np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10));original_ok=True
        try:raw=m.candidate(c)
        except ValueError:original_ok=False
        row=dict(iteration=it,loss=loss,min_constraint=feas,bounds_ok=inbounds,original_parameters_ok=original_ok,seconds=time.monotonic()-start)
        history.append(row)
        if feas>=-1e-7 and inbounds and original_ok and loss<best[0]:
            best[:]=[loss,c.copy()];m.f.save(raw,out/'best_feasible.json',row);np.save(out/'best_feasible_modifiers.npy',c)
        if it%5==0:
            np.save(out/'checkpoint.npy',c)
            if original_ok:m.f.save(raw,out/'checkpoint.json',row)
            (out/'history.json').write_text(json.dumps(history,indent=2));print(json.dumps(row),flush=True)
        if time.monotonic()-start>seconds:raise BudgetStop
    try:
        ret=minimize(fun,zero,jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:constraints(y)[0],jac=lambda y:constraints(y)[1])],callback=callback,options=dict(maxiter=maxiter,ftol=2e-9,disp=True))
        final=T@ret.x;success=bool(ret.success);message=str(ret.message);nfev=ret.nfev
    except BudgetStop:
        final=T@last[0];success=False;message='Wall budget checkpoint, no convergence claim';nfev=len(obj.history)
    if best[1] is not None:final=best[1];selection='best_feasible_training_checkpoint'
    else:selection='last_training_iterate_no_feasible_found'
    summary=dict(**config,conditioning=info,success=success,message=message,selection=selection,nfev=nfev,iterations=len(history),loss=obj.fun(final)[0],min_training_constraint=float(obj.constraints(final)[0].min()),elapsed=time.monotonic()-start,variables=m.dim,pde_validated=False,source_correspondence_verified=False)
    m.f.save(m.candidate(final),out/'candidate.json',summary);np.save(out/'modifiers.npy',final)
    (out/'summary.json').write_text(json.dumps(summary,indent=2));(out/'history.json').write_text(json.dumps(history,indent=2));(out/'projection.json').write_text(json.dumps(m.projection,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--method',default='whiten',choices=['whiten','diagonal','none']);p.add_argument('--maxiter',type=int,default=180);p.add_argument('--seconds',type=float,default=480);p.add_argument('--radial',type=int,default=4);p.add_argument('--axial',type=int,default=4);p.add_argument('--td',type=int,default=4);p.add_argument('--pressure',type=float,default=.99);p.add_argument('--profile',type=float,default=.995);p.add_argument('--acc',type=float,default=.98);p.add_argument('--axis-weight',type=float,default=.02)
    a=p.parse_args();run(a.parent,a.out,a.method,a.maxiter,a.seconds,a.radial,a.axial,a.td,a.pressure,a.profile,a.acc,a.axis_weight)
