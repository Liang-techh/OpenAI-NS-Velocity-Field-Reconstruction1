"""ST048: support-edge residual control with signed midplane shear preservation.
All new probe sets and target ratios are autonomous finite-window training choices.
The physical field, original force, normalization and acceptance gates are unchanged.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'experiments/root_st047'))
from continuation import Objective as PriorObjective, LocalizedModel, precondition

class BoundaryShearObjective(PriorObjective):
    def __init__(self, m, *, shear_ratio=.995, edge_weight=.04, edge_temperature=.0001, **kw):
        super().__init__(m, **kw)
        self.shear_ratio = shear_ratio
        self.edge_weight = edge_weight
        self.edge_temperature = edge_temperature
        # Finite physical-radius probes, not merely the earlier scaled radial difference.
        r,t = np.meshgrid([.04,.08,.14,.23,.36,.48,.60], np.linspace(.25,.75,17), indexing='ij')
        self.shear_cache = m.cache(r*r, np.zeros_like(r), t)
        parent,_ = self.shear(np.zeros(m.dim))
        self.shear_parent = parent
        self.shear_sign = np.sign(parent)
        self.shear_scale = np.maximum(abs(parent), 2e-4)
        # Two cylindrical face neighborhoods, with endpoints explicitly included.
        et = np.array([.25,.27,.33,.50,.67,.73,.75])
        rr = np.r_[np.linspace(0,.5,6), np.linspace(.65,1.95,8)]
        zz = np.r_[-np.linspace(1.05,1.96,12)[::-1],np.linspace(1.05,1.96,12)]
        r1,z1,t1 = np.meshgrid(rr,zz,et,indexing='ij')
        r2,z2,t2 = np.meshgrid(np.linspace(1.25,1.96,9),np.linspace(-1.95,1.95,17),et,indexing='ij')
        self.edge_points = np.r_[np.c_[r1.ravel()**2,z1.ravel(),t1.ravel()],np.c_[r2.ravel()**2,z2.ravel(),t2.ravel()]]
        self.edge_cache = m.cache(*self.edge_points.T)
        self.edge_last = None
        self.final_last = None
        self.constraints_last = None
    def shear(self,c):
        m=self.m;D=self.shear_cache
        factor=2*np.sqrt(D['s']);lam,dl=m.norm(c)
        bar=D['v']['Cs']+D['M']['Cs']@c[:m.na]
        val=lam*factor*bar
        J=np.zeros((len(val),m.dim))
        J[:,:m.na]=lam*factor[:,None]*D['M']['Cs']
        J[:,:m.nv]+=(factor*bar)[:,None]*dl[None,:]
        return val,J
    def edge_momentum(self,c):
        if self.edge_last is None or not np.array_equal(c,self.edge_last[0]):
            self.edge_last=(c.copy(),*self.m.momentum(c,self.edge_cache))
        return self.edge_last[1:]
    def constraints(self,c):
        if self.constraints_last is not None and np.array_equal(c,self.constraints_last[0]):
            return self.constraints_last[1:]
        v,J=super().constraints(c)
        shear,S=self.shear(c)
        # Preserve signed shear pointwise to the specified relative tolerance.
        low=(self.shear_sign*shear-self.shear_ratio*abs(self.shear_parent))/self.shear_scale
        LJ=self.shear_sign[:,None]*S/self.shear_scale[:,None]
        # Prevent uncontrolled shear growth as a spurious structural "improvement".
        high=(1.25*abs(self.shear_parent)+1e-5-self.shear_sign*shear)/self.shear_scale
        self.constraints_last=(c.copy(),np.r_[v,low,high],np.vstack((J,LJ,-LJ)))
        return self.constraints_last[1:]
    def fun(self,c):
        if self.final_last is not None and np.array_equal(c,self.final_last[0]):
            return self.final_last[1:]
        value,g=super().fun(c)
        r,J=self.edge_momentum(c);sq=np.sum(r*r,axis=1)
        logw=sq/self.edge_temperature;normalizer=logsumexp(logw)
        weights=np.exp(logw-normalizer)
        value+=self.edge_weight*self.edge_temperature*(normalizer-np.log(len(sq)))+.01*float(np.mean(sq))
        g=g+2*np.einsum('n,ni,nik->k',self.edge_weight*weights+.01/len(sq),r,J)
        self.final_last=(c.copy(),float(value),g)
        return self.final_last[1:]

def run(parent,out,*,radial=5,axial=5,td=4,maxiter=240,seconds=600,pressure=.99,profile=.995,acc=.98,shear=.995,edge_weight=.04,ident='ST048-S'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Output directory must be empty')
    reg=dict(experiment=ident,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),radial=radial,axial=axial,time_degree=td,maxiter=maxiter,wall_seconds=seconds,pressure_ratio=pressure,profile_ratio=profile,acceleration_ratio=acc,signed_shear_ratio=shear,edge_weight=edge_weight,edge_temperature=.0001,training_space=[32,48],training_time='13Gauss plus endpoints',edge_times=[.25,.27,.33,.5,.67,.73,.75],shear_radii=[.04,.08,.14,.23,.36,.48,.6],constraint_times='17uniform',training_feasibility_tolerance=1e-7,validation_seeds=[9174801,9174802],original_physical_gates='UNCHANGED',created_utc=datetime.now(timezone.utc).isoformat(),registered_before_fit=True)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    t0=time.monotonic();m=LocalizedModel(parent,radial,axial,td)
    obj=BoundaryShearObjective(m,shear_ratio=shear,edge_weight=edge_weight,acceleration_ratio=acc,axis_weight=.02,profile_ratio=profile,pressure_ratio=pressure,space_order=(32,48),peak_weight=0.,softmax_weight=.02,anchor_tolerance=.08,time_cap=1.)
    zero=np.zeros(m.dim);T,rs,info=precondition(obj,zero,'whiten');factor=1e4
    (out/'conditioning.json').write_text(json.dumps(info,indent=2)+'\n')
    np.save(out/'coordinate_transform.npy',T);np.save(out/'row_scale.npy',rs)
    np.savez_compressed(out/'model.npz',Ma=m.Ma,Mb=m.Mb,Mp=m.Mp)
    np.save(out/'edge_training_points.npy',obj.edge_points)
    lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    def fun(y):
        v,g=obj.fun(T@y);return factor*v,factor*T.T@g
    cc=[None,None]
    def con(y):
        if cc[0] is not None and np.array_equal(y,cc[0]):return cc[1]
        c=T@y;v,J=obj.constraints(c)
        cc[:]=[y.copy(),(np.r_[v*rs,(c-lo)*bs,(hi-c)*bs],np.vstack(((J@T)*rs[:,None],T*bs[:,None],-T*bs[:,None])))]
        return cc[1]
    rng=np.random.default_rng(9174890);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=2e-6
    v,g=fun(zero);v0,J=con(zero);vp=con(h*d)[0];vm=con(-h*d)[0]
    check=dict(objective_directional_error=abs((fun(h*d)[0]-fun(-h*d)[0])/(2*h)-g@d),constraint_directional_error=float(np.max(abs((vp-vm)/(2*h)-J@d))),parent_shear_range=[float(obj.shear_parent.min()),float(obj.shear_parent.max())],variables=m.dim,edge_probes=len(obj.edge_points))
    (out/'derivative_checks.json').write_text(json.dumps(check,indent=2)+'\n');print(json.dumps(check),flush=True)
    if check['constraint_directional_error']>1e-4 or check['objective_directional_error']>1e-4:raise ValueError('Derivative calibration failed')
    hist=[];best=[np.inf,None];last=[zero.copy()]
    class BudgetStop(Exception):pass
    def callback(y):
        last[0]=y.copy();c=T@y;loss=obj.fun(c)[0];feas=float(obj.constraints(c)[0].min());it=len(hist)+1
        inbounds=bool(np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10));raw=None
        try:raw=m.candidate(c)
        except ValueError:pass
        row=dict(iteration=it,loss=loss,min_constraint=feas,bounds_ok=inbounds,original_parameters_ok=raw is not None,seconds=time.monotonic()-t0)
        hist.append(row)
        if feas>=-1e-7 and inbounds and raw is not None and loss<best[0]:
            best[:]=[loss,c.copy()];m.f.save(raw,out/'best_feasible.json',row);np.save(out/'best_feasible_modifiers.npy',c)
        if it%5==0:
            np.save(out/'checkpoint.npy',c)
            if raw is not None:m.f.save(raw,out/'checkpoint.json',row)
            (out/'history.json').write_text(json.dumps(hist,indent=2)+'\n');print(json.dumps(row),flush=True)
        if time.monotonic()-t0>seconds:raise BudgetStop
    try:
        ret=minimize(fun,zero,jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:con(y)[0],jac=lambda y:con(y)[1])],callback=callback,options=dict(maxiter=maxiter,ftol=2e-9,disp=True))
        final=T@ret.x;success=bool(ret.success);message=str(ret.message);nfev=ret.nfev
    except BudgetStop:
        final=T@last[0];success=False;message='Wall budget checkpoint; no convergence claim';nfev=len(obj.history)
    selected=best[1] if best[1] is not None else final
    summary=dict(**reg,success=success,message=message,iterations=len(hist),nfev=nfev,variables=m.dim,selection='best_feasible_training' if best[1] is not None else 'infeasible_last_iterate',min_constraint=float(obj.constraints(selected)[0].min()),loss=obj.fun(selected)[0],elapsed=time.monotonic()-t0,pde_validated=False,source_correspondence_verified=False)
    m.f.save(m.candidate(selected),out/'candidate.json',summary);np.save(out/'modifiers.npy',selected)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'history.json').write_text(json.dumps(hist,indent=2)+'\n');(out/'projection.json').write_text(json.dumps(m.projection,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--radial',type=int,default=5);p.add_argument('--axial',type=int,default=5);p.add_argument('--td',type=int,default=4);p.add_argument('--maxiter',type=int,default=240);p.add_argument('--seconds',type=float,default=600);p.add_argument('--pressure',type=float,default=.99);p.add_argument('--profile',type=float,default=.995);p.add_argument('--acc',type=float,default=.98);p.add_argument('--shear',type=float,default=.995);p.add_argument('--edge-weight',type=float,default=.04);p.add_argument('--id',default='ST048-S')
    a=p.parse_args();run(a.parent,a.out,radial=a.radial,axial=a.axial,td=a.td,maxiter=a.maxiter,seconds=a.seconds,pressure=a.pressure,profile=a.profile,acc=a.acc,shear=a.shear,edge_weight=a.edge_weight,ident=a.id)
