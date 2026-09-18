"""ST051: jointly reduce momentum and retain pressure alignment in the frozen P field.
New finite training probes/basis directions are autonomous; original NS gates unchanged.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st050r'))
from pressure_morph import PressureMorph, LocalizedModel, precondition, atomic_json
from spacetime import Family

class EdgeModel(LocalizedModel):
    """Extra high-order support-edge columns in the SAME actual spatial field basis."""
    def __init__(self,path,radial=5,axial=5,td=4,edge_modes=True):
        super().__init__(path,radial,axial,td)
        self.edge_modes=edge_modes
        self.edge_pairs=([(i,j) for i in (0,1,2) for j in (7,8)] + [(i,j) for i in (7,8) for j in (0,1)]) if edge_modes else []
        if not edge_modes:return
        f=self.f;old=[self.na,self.nb,self.np];bounds=self.bounds
        groups=[bounds[:old[0]],bounds[old[0]:sum(old[:2])],bounds[sum(old[:2]):-2]]
        mats=[];newbounds=[]
        for block,(M,T,bb) in enumerate(zip((self.Ma,self.Mb,self.Mp),(f.Tp,f.Tw,f.Tq),groups)):
            extra=[]
            for i,j in self.edge_pairs:
                p=np.zeros((f.nr,f.nz));p[i,j]=1
                v=np.linalg.solve(T,p.ravel());v/=np.linalg.norm(v)
                for k in range(td):
                    m=np.zeros((f.ns,f.nt));m[:,k]=v;extra.append(m.ravel())
            mats.append(np.column_stack((M,np.column_stack(extra))))
            newbounds += bb+[(-.02,.02) if block<2 else (-.08,.08)]*len(extra)
        self.Ma,self.Mb,self.Mp=mats;self.na,self.nb,self.np=[M.shape[1] for M in mats]
        self.nv=self.na+self.nb;self.dim=self.nv+self.np+2;self.bounds=newbounds+bounds[-2:]
        self.EM=np.zeros((2*f.ns,self.nv))
        self.EM[:f.ns,:self.na]=np.einsum('stk,t->sk',self.Ma.reshape(f.ns,f.nt,-1),f.q0)
        self.EM[f.ns:,self.na:]=np.einsum('stk,t->sk',self.Mb.reshape(f.ns,f.nt,-1),f.q0)

class AlignedObjective(PressureMorph):
    def __init__(self,m,*,shear_reference=None,robust_grid=True,**kw):
        super().__init__(m,**kw)
        self.robust_grid=robust_grid;self.extra_last=None
        if not robust_grid:return
        # A new training grid, not either old or new validation sample.
        R,Z,T=np.meshgrid([.035,.065,.105,.155,.215],[-.215,-.14,-.085,-.045,-.025,.025,.045,.085,.14,.215],np.linspace(.25,.75,19),indexing='ij')
        self.rob=m.cache(R*R*(1-T),Z*(1-T)**.495,T)
        self.rob_sign=np.sign(Z.ravel())
        self.rob_scales=np.stack(((1-T)**.5,(1-T)**.505,(1-T)**.505),axis=-1).reshape(-1,3)
        self.rob_signs=np.stack((-np.ones(Z.size),np.ones(Z.size),self.rob_sign),axis=1)
        r,t=np.meshgrid(np.linspace(.035,.62,13),np.linspace(.25,.75,19),indexing='ij')
        self.rs=m.cache(r*r,np.zeros_like(r),t)
        sh,_=self.robust_shear(np.zeros(m.dim));self.ref_shear=sh.copy()
        if shear_reference is not None:
            f,raw=Family.load(shear_reference);a=f.coefficients(raw)[0]
            # The reference is in the same actual field basis but separate immutable field.
            vals=[]
            for i in range(0,r.size,128):
                sl=slice(i,i+128);rr=r.ravel()[sl];D=f.bundle(rr*rr,np.zeros(len(rr)),t.ravel()[sl]);vals.append(2*rr*(D['Cs']@a))
            self.ref_shear=np.concatenate(vals)
        self.ref_sign=np.sign(self.ref_shear);self.ref_scale=np.maximum(abs(self.ref_shear),2e-4)
    def robust_shear(self,c):
        m=self.m;D=self.rs;lam,dl=m.norm(c);r=2*np.sqrt(D['s']);bar=D['v']['Cs']+D['M']['Cs']@c[:m.na]
        J=np.zeros((len(r),m.dim));J[:,:m.na]=lam*r[:,None]*D['M']['Cs'];J[:,:m.nv]+=(r*bar)[:,None]*dl
        return lam*r*bar,J
    def constraints(self,c):
        if self.extra_last is not None and np.array_equal(self.extra_last[0],c):return self.extra_last[1:]
        v,J=super().constraints(c)
        if self.robust_grid:
            m=self.m;D=self.rob;sl=slice(m.nv,m.nv+m.np)
            pg=D['v']['Qz']+D['M']['Qz']@c[sl];ps=D['v']['Qs']+D['M']['Qs']@c[sl]
            P=np.zeros((len(pg),m.dim));P[:,sl]=self.rob_sign[:,None]*D['M']['Qz']/self.pgrad_scale
            S=np.zeros_like(P);S[:,sl]=D['M']['Qs']/self.radial_scale
            pv=(self.rob_sign*pg-1e-5)/self.pgrad_scale;sv=(ps-1e-7)/self.radial_scale
            u,K=m.velocity(c,D);scale=max(float(np.linalg.norm(self.coreparent[0])),.01)
            uv=((u*self.rob_scales*self.rob_signs)-1e-7).ravel()/scale
            K=np.pad((K*(self.rob_scales*self.rob_signs)[:,:,None]).reshape(-1,m.nv)/scale,((0,0),(0,m.dim-m.nv)))
            sh,SH=self.robust_shear(c)
            shv=(self.ref_sign*sh-.995*abs(self.ref_shear))/self.ref_scale
            SH=self.ref_sign[:,None]*SH/self.ref_scale[:,None]
            v=np.r_[v,pv,sv,uv,shv];J=np.vstack((J,P,S,K,SH))
        self.extra_last=(c.copy(),v,J);return v,J

def run(parent,out,*,ident='ST051-A',edge=False,maxiter=260,seconds=600,anchor=.12,pp=.003,axis=1.75,reference=None,robust=True,cap=1.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Refusing nonempty experiment output')
    registration=dict(id=ident,parent_file=str(parent),parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),shear_reference_sha256=hashlib.sha256(Path(reference).read_bytes()).hexdigest() if reference else None,edge_modes=edge,radial=5,axial=5,td=4,anchor_tolerance=anchor,poisson_weight=pp,axis_min_swirl_strain_ratio=axis,robust_grid=robust,pressure_target=0.,robust_pressure_margin=1e-5,profile_ratio=.995,acceleration_ratio=1.,morphology_ratio=.9999,shear_parent_ratio=.999,shear_reference_ratio=.995,per_time_mse_cap=cap,training_space=[32,48],training_time='13Gauss plus endpoints',constraint_times=17,robust_times=19,wall_seconds=seconds,maxiter=maxiter,feasibility_tolerance=1e-7,validation_seeds=[9175191,9175192],original_physical_gates='UNCHANGED',registered_utc=datetime.now(timezone.utc).isoformat(),pde_validated=False)
    atomic_json(out/'registration.json',registration)
    start=time.monotonic();m=EdgeModel(parent,edge_modes=edge)
    o=AlignedObjective(m,shear_reference=reference,robust_grid=robust,poisson_weight=pp,morph_ratio=.9999,axis_ratio=axis,shear_ratio=.999,edge_weight=.08,pressure_target=0.,profile_ratio=.995,pressure_ratio=1.,acceleration_ratio=1.,axis_weight=.02,space_order=(32,48),anchor_tolerance=anchor,time_cap=cap,peak_weight=0.,softmax_weight=.03)
    zero=np.zeros(m.dim);T,rs,info=precondition(o,zero,'whiten');atomic_json(out/'conditioning.json',info)
    np.savez_compressed(out/'basis_maps.npz',Ma=m.Ma,Mb=m.Mb,Mp=m.Mp);np.save(out/'solver_transform.npy',T)
    factor=1e4;lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    counts={'fun':0,'constraint':0};cache=[None,None]
    def fun(y):
        counts['fun']+=1;v,g=o.fun(T@y);return factor*v,factor*T.T@g
    def con(y):
        if cache[0] is not None and np.array_equal(cache[0],y):return cache[1]
        counts['constraint']+=1;c=T@y;v,J=o.constraints(c)
        cache[:]=[y.copy(),(np.r_[v*rs,(c-lo)*bs,(hi-c)*bs],np.vstack(((J@T)*rs[:,None],T*bs[:,None],-T*bs[:,None])))];return cache[1]
    rng=np.random.default_rng(9175150);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=1e-6
    v,g=fun(zero);v0,J=con(zero)
    errors={'objective_absolute':abs((fun(h*d)[0]-fun(-h*d)[0])/(2*h)-g@d),'constraint_max':float(np.max(abs((con(h*d)[0]-con(-h*d)[0])/(2*h)-J@d))),'initial_min_constraint':float(o.constraints(zero)[0].min()),'variables':m.dim}
    atomic_json(out/'derivative_check.json',errors);print('CALIBRATION',json.dumps(errors),flush=True)
    if max(errors['objective_absolute'],errors['constraint_max'])>2e-4:raise ValueError('Derivative calibration failed')
    history=[];best=[np.inf,None];last=[zero.copy()]
    class BudgetStop(Exception):pass
    def callback(y):
        last[0]=y.copy();c=T@y;loss=o.fun(c)[0];feas=float(o.constraints(c)[0].min());raw=None
        try:raw=m.candidate(c)
        except ValueError:pass
        good=bool(np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10) and raw is not None)
        row=dict(iteration=len(history)+1,loss=loss,min_constraint=feas,bounds_ok=good,elapsed=time.monotonic()-start);history.append(row)
        if good and feas>=-1e-7 and loss<best[0]:
            best[:]=[loss,c.copy()];m.f.save(raw,out/'best_feasible.json',row);np.save(out/'best_modifiers.npy',c)
        if len(history)%5==0:
            atomic_json(out/'history.json',history);np.save(out/'checkpoint.npy',c)
            if raw is not None:m.f.save(raw,out/'checkpoint.json',row)
            print(json.dumps(row),flush=True)
        if time.monotonic()-start>seconds:raise BudgetStop
    counts={'fun':0,'constraint':0}
    try:
        r=minimize(fun,zero,jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:con(y)[0],jac=lambda y:con(y)[1])],callback=callback,options=dict(maxiter=maxiter,ftol=2e-9,disp=False))
        final=T@r.x;success=bool(r.success);message=str(r.message);nfev=int(r.nfev)
    except BudgetStop:
        final=T@last[0];success=False;message='Registered wall budget; checkpoint retained';nfev=counts['fun']
    selected=best[1] if best[1] is not None else final
    summary=dict(**registration,variables=m.dim,optimizer_success=success,message=message,iterations=len(history),fun_calls=nfev,call_counts=counts,min_constraint=float(o.constraints(selected)[0].min()),loss=o.fun(selected)[0],selection='best_feasible_training' if best[1] is not None else 'infeasible_last_checkpoint',elapsed=time.monotonic()-start)
    raw=m.candidate(selected);m.f.save(raw,out/'candidate.json',summary);np.save(out/'modifiers.npy',selected)
    atomic_json(out/'summary.json',summary);atomic_json(out/'history.json',history);print('SUMMARY',json.dumps(summary),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',default='ST051-A');p.add_argument('--edge',action='store_true');p.add_argument('--maxiter',type=int,default=260);p.add_argument('--seconds',type=float,default=600);p.add_argument('--anchor',type=float,default=.12);p.add_argument('--pp',type=float,default=.003);p.add_argument('--axis',type=float,default=1.75);p.add_argument('--reference');p.add_argument('--no-robust',action='store_true');p.add_argument('--cap',type=float,default=1.)
    a=p.parse_args();run(a.parent,a.out,ident=a.id,edge=a.edge,maxiter=a.maxiter,seconds=a.seconds,anchor=a.anchor,pp=a.pp,axis=a.axis,reference=a.reference,robust=not a.no_robust,cap=None if a.cap<0 else a.cap)
