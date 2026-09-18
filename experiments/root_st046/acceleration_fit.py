"""ST046: joint non-pressure axial momentum, pressure and full-residual continuation.
New local search directions only; physical basis and original gates unchanged.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2]
for name in ('root_st030','root_st040','root_st043'):
    sys.path.insert(0,str(ROOT/'experiments'/name))
from coupled_fit import JointObjective
from localized_model import LocalizedModel

class AccelerationObjective(JointObjective):
    def __init__(self,m,acceleration_ratio=.95,**kwargs):
        super().__init__(m,**kwargs)
        self.acceleration_ratio=acceleration_ratio
        self._cj=None;self._cr=None
        R,_=m.momentum(np.zeros(m.dim),self.exp)
        self.parent_mz=self.signz*(R[:,2]-self.exp['v']['Qz'])
        self.accscale=max(float(np.sqrt(np.mean(self.parent_mz**2))),.01)
        self.acc_target=self.parent_mz-(1-acceleration_ratio)*np.maximum(self.parent_mz,0)
    def constraints(self,c):
        if self._cj is not None and np.array_equal(c,self._cj):return self._cr
        V,J=super().constraints(c)
        R,K=self.m.momentum(c,self.exp)
        pg=self.exp['v']['Qz']+self.exp['M']['Qz']@c[self.m.nv:self.m.nv+self.m.np]
        M=self.signz*(R[:,2]-pg)
        MK=self.signz[:,None]*K[:,2].copy()
        MK[:,self.m.nv:self.m.nv+self.m.np]-=self.signz[:,None]*self.exp['M']['Qz']
        self._cj=c.copy();self._cr=(np.r_[V,(self.acc_target-M)/self.accscale],np.vstack((J,-MK/self.accscale)))
        return self._cr

class AxisPeakObjective(AccelerationObjective):
    def __init__(self,m,acceleration_ratio=.95,axis_weight=.05,**kwargs):
        super().__init__(m,acceleration_ratio,**kwargs)
        R,Z,T=np.meshgrid([0.,.01,.04,.12,.2],np.r_[-.2,-.14,-.08,-.02,0.,.02,.08,.14,.2],np.linspace(.25,.75,17),indexing='ij')
        self.axis=m.cache(R*R*(1-T),Z*(1-T)**.495,T)
        self.axis_weight=axis_weight;self._alast=None;self._aret=None
    def fun(self,c):
        if self._alast is not None and np.array_equal(c,self._alast):return self._aret
        value,g=super().fun(c);R,J=self.m.momentum(c,self.axis);norm=np.sum(R*R,axis=1)
        temp=.0005;lw=norm/temp;prob=np.exp(lw-logsumexp(lw))
        loss=value+self.axis_weight*temp*(logsumexp(lw)-np.log(len(norm)))+.02*float(np.mean(norm))
        grad=g+2*np.einsum('n,ni,nik->k',self.axis_weight*prob+.02/len(norm),R,J)
        self._alast=c.copy();self._aret=(loss,grad)
        return self._aret

def run(parent,out,radial=4,axial=4,td=4,acceleration_ratio=.95,pressure_ratio=.97,profile_ratio=.98,maxiter=600,cap=1.,budget=900.,axis_weight=0.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Refuse nonempty output')
    reg=dict(experiment='ST046',parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),radial=radial,axial=axial,time_degree=td,acceleration_ratio=acceleration_ratio,pressure_ratio=pressure_ratio,profile_ratio=profile_ratio,maxiter=maxiter,time_cap=cap,wall_budget_seconds=budget,axis_peak_weight=axis_weight,training_space=[32,48],training_time='13 Gauss+endpoints',original_physical_gates='UNCHANGED',validation_seeds=[9174601,9174602],written_before_fit=True)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    start=time.monotonic();m=LocalizedModel(parent,radial,axial,td)
    otype=AxisPeakObjective if axis_weight else AccelerationObjective
    extra={"axis_weight":axis_weight} if axis_weight else {}
    obj=otype(m,acceleration_ratio,**extra,profile_ratio=profile_ratio,pressure_ratio=pressure_ratio,space_order=(32,48),peak_weight=0.,softmax_weight=.02,anchor_tolerance=.1,time_cap=cap)
    c=np.zeros(m.dim)
    rng=np.random.default_rng(9174590);v=rng.normal(size=m.dim);v/=np.linalg.norm(v);h=1e-6
    loss,g=obj.fun(c);cv,cj=obj.constraints(c);fp=obj.fun(c+h*v)[0];fm=obj.fun(c-h*v)[0];vp=obj.constraints(c+h*v)[0];vm=obj.constraints(c-h*v)[0]
    check=dict(objective_error=float(abs((fp-fm)/(2*h)-g@v)),constraint_error=float(np.max(abs((vp-vm)/(2*h)-cj@v))),variables=m.dim)
    (out/'derivative_checks.json').write_text(json.dumps(check,indent=2)+'\n');print(check,flush=True)
    if check['constraint_error']>2e-5:raise ValueError('Derivative check failed')
    best=[np.inf,None];last=[c.copy()];trace=[]
    class BudgetStop(Exception):pass
    def callback(x):
        last[0]=x.copy();i=len(trace)+1;value=obj.fun(x)[0];feas=float(np.min(obj.constraints(x)[0]))
        trace.append(dict(iteration=i,loss=value,min_constraint=feas,elapsed=time.monotonic()-start))
        if feas>=-1e-7 and value<best[0]:best[:]=[value,x.copy()]
        if i%10==0:
            np.save(out/'checkpoint.npy',x);m.f.save(m.candidate(x),out/'checkpoint.json',trace[-1]);(out/'iterations.json').write_text(json.dumps(trace,indent=2)+'\n');print(trace[-1],flush=True)
        if time.monotonic()-start>budget:raise BudgetStop
    try:
        opt=minimize(obj.fun,c,jac=True,method='SLSQP',bounds=m.bounds,constraints=[{'type':'ineq','fun':lambda x:obj.constraints(x)[0],'jac':lambda x:obj.constraints(x)[1]}],callback=callback,options=dict(maxiter=maxiter,ftol=2e-12,disp=True))
        chosen=opt.x;success=bool(opt.success);message=str(opt.message);nit=opt.nit;nfev=opt.nfev
    except BudgetStop:
        chosen=last[0];success=False;message='Wall budget reached';nit=len(trace);nfev=len(obj.history)
    if best[1] is not None and obj.constraints(chosen)[0].min()<-1e-7:chosen=best[1];message+='; use best feasible TRAINING checkpoint'
    summary=dict(**reg,success=success,message=message,iterations=nit,nfev=nfev,variables=m.dim,loss=float(obj.fun(chosen)[0]),minimum_constraint=float(obj.constraints(chosen)[0].min()),elapsed=time.monotonic()-start,pde_validated=False,source_correspondence_verified=False)
    m.f.save(m.candidate(chosen),out/'candidate.json',summary)
    np.savez_compressed(out/'recipe.npz',modifiers=chosen,Ma=m.Ma,Mb=m.Mb,Mp=m.Mp)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'iterations.json').write_text(json.dumps(trace,indent=2)+'\n');(out/'history.json').write_text(json.dumps(obj.history,indent=2)+'\n');(out/'projection.json').write_text(json.dumps(m.projection,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--radial',type=int,default=4);p.add_argument('--axial',type=int,default=4);p.add_argument('--td',type=int,default=4);p.add_argument('--acc',type=float,default=.95);p.add_argument('--pressure',type=float,default=.97);p.add_argument('--profile',type=float,default=.98);p.add_argument('--maxiter',type=int,default=600);p.add_argument('--budget',type=float,default=900);p.add_argument('--cap',type=float,default=1.);p.add_argument('--axis-weight',type=float,default=0.)
    a=p.parse_args();run(a.parent,a.out,a.radial,a.axial,a.td,a.acc,a.pressure,a.profile,a.maxiter,a.cap,a.budget,a.axis_weight)
