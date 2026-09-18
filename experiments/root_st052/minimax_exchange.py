"""ST052 sampled minimax exchange, keeping the inherited actual compact NS field.
The epigraph is a finite training bound, never a continuum or scientific certificate.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time, gc
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st051'))
from aligned_continuation import EdgeModel,AlignedObjective,precondition,atomic_json
PARENT_SHA='0071106ef10a5d77b620b942becc00b58c8a6765fd5fb7bc4dcd195ac65c970d'

def pool_points(seed=9175250):
    rng=np.random.default_rng(seed)
    r=np.unique(np.r_[np.linspace(0,1.985,27),np.linspace(.17,.85,12)])
    z=np.unique(np.r_[-np.linspace(1.6,1.985,17),np.linspace(-1.6,1.6,25),np.linspace(1.6,1.985,17)])
    S,Z,T=np.meshgrid(r*r,z,np.linspace(.25,.75,9),indexing='ij')
    grid=np.c_[S.ravel(),Z.ravel(),T.ravel()]
    random=np.c_[rng.uniform(0,4,3072),rng.uniform(-2,2,3072),rng.uniform(.25,.75,3072)]
    return np.r_[grid,random],len(grid)

def scan(m,c,pool):
    raw=m.candidate(c);values=[]
    for i in range(0,len(pool),384):
        s,z,t=pool[i:i+384].T;x=np.c_[np.sqrt(s),np.zeros(len(s)),z]
        values.append(np.linalg.norm(m.f.analytic_residual(raw,x,t),axis=1))
    return np.concatenate(values)

def select_points(pool,values,n_grid,old=()):
    chosen=list(old)
    for t in np.unique(pool[:n_grid,2]):
        ids=np.flatnonzero(pool[:n_grid,2]==t)
        chosen.extend(ids[np.argsort(values[ids])[-48:]])
    ids=np.arange(n_grid,len(pool));chosen.extend(ids[np.argsort(values[ids])[-96:]])
    return np.unique(chosen).astype(int)

class Epigraph:
    def __init__(self,m,o,points,pscale,base_loss):
        self.m=m;self.o=o;self.D=m.cache(*points.T);self.pscale=pscale;self.base_loss=base_loss
        self.memo=None
    def constraints(self,c,q):
        R,J=self.m.momentum(c,self.D);sq=np.sum(R*R,axis=1)
        return q-sq/self.pscale, -2*np.einsum('ni,nik->nk',R,J)/self.pscale
    def objective(self,c,q):
        loss,g=self.o.fun(c)
        return float(q+.10*loss/self.base_loss),.10*g/self.base_loss

def run(parent,out,rounds=2,maxiter=170,seconds=540):
    parent=Path(parent);out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Refuse nonempty output')
    if hashlib.sha256(parent.read_bytes()).hexdigest()!=PARENT_SHA:raise ValueError('Wrong frozen parent')
    out.mkdir(parents=True,exist_ok=True)
    reg=dict(id='ST052-M',parent_sha256=PARENT_SHA,base_commit='4b784f1b8457af2ead49295631d834d4e882000b',rounds=rounds,maxiter_per_round=maxiter,wall_seconds_per_round=seconds,seed=9175250,holdout_seeds=[9175291,9175292],training_pool='39r x57z x9t plus3072uniform cylinder-space-time points, see saved pool',selection='feasible iterate with best full TRAINING-pool maximum at round boundaries; no holdout selection',objective='q + .10 * inherited_structure_aware_loss / parent_loss; q>=|R(x,t)|^2/parent_pool_max^2 at exchanged training points',structure='Inherited ST051 AlignedObjective; anchor .08 relative B, profile ratio1(no worsening), acceleration1, morphology.9999, shear.999relativeB plus.995fixedS, axis ratio1.75; per-training-time L2cap1.0',physical_gates='Original nu/time/support/force/E0/all numerical thresholds unchanged',utc_before_fit=datetime.now(timezone.utc).isoformat(),pde_validated=False)
    atomic_json(out/'registration.json',reg)
    m=EdgeModel(parent,edge_modes=True)
    o=AlignedObjective(m,shear_reference=ROOT/'artifacts/research/ST048-S/candidate.json',robust_grid=True,poisson_weight=.003,morph_ratio=.9999,axis_ratio=1.75,shear_ratio=.999,edge_weight=.08,pressure_target=0.,profile_ratio=1.,pressure_ratio=1.,acceleration_ratio=1.,axis_weight=.02,space_order=(32,48),anchor_tolerance=.08,time_cap=1.,peak_weight=0.,softmax_weight=.03)
    c=np.zeros(m.dim);T,rs,info=precondition(o,c,'whiten');atomic_json(out/'conditioning.json',info)
    np.save(out/'coordinate_map.npy',T)
    pool,ng=pool_points();np.savez_compressed(out/'training_pool.npz',points=pool,n_grid=ng)
    values=scan(m,c,pool);parent_peak=float(values.max());pscale=parent_peak**2
    base_loss=o.fun(c)[0];lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4)
    selected_c=c.copy();selected_peak=parent_peak;report=[];indices=[]
    atomic_json(out/'initial_training.json',dict(peak=parent_peak,loss=base_loss,min_constraint=float(o.constraints(c)[0].min()),points=len(pool),variables=m.dim,epigraph_variables=m.dim+1))
    print('INITIAL',out/'initial_training.json',parent_peak,base_loss,flush=True)
    for rr in range(rounds):
        dest=out/f'round{rr+1}';dest.mkdir()
        indices=select_points(pool,values,ng,indices);np.save(dest/'active_pool_indices.npy',indices)
        e=Epigraph(m,o,pool[indices],pscale,base_loss)
        y0=np.r_[np.linalg.solve(T,c),(float(values[indices].max())**2/pscale)*1.0001]
        memo=[None,None];counts={'fun':0,'constraint':0};hist=[];last=[y0.copy()];start=time.monotonic()
        def fun(y):
            counts['fun']+=1;v,g=e.objective(T@y[:-1],y[-1]);return v,np.r_[T.T@g,1.]
        def con(y):
            if memo[0] is not None and np.array_equal(memo[0],y):return memo[1]
            counts['constraint']+=1;cc=T@y[:-1];v,J=o.constraints(cc);p,K=e.constraints(cc,y[-1])
            values_=np.r_[v*rs,p,(cc-lo)*bs,(hi-cc)*bs,y[-1],4-y[-1]]
            j=np.vstack((np.pad((J@T)*rs[:,None],((0,0),(0,1))),np.c_[K@T,np.ones(len(p))],np.pad(T*bs[:,None],((0,0),(0,1))),np.pad(-T*bs[:,None],((0,0),(0,1))),np.r_[np.zeros(m.dim),1.],np.r_[np.zeros(m.dim),-1.]))
            memo[:]=[y.copy(),(values_,j)];return memo[1]
        if rr==0:
            rng=np.random.default_rng(9175251);d=rng.normal(size=len(y0));d/=np.linalg.norm(d);h=1e-6
            fv,fg=fun(y0);v,J=con(y0)
            check=dict(objective_abs_error=abs((fun(y0+h*d)[0]-fun(y0-h*d)[0])/(2*h)-fg@d),constraint_abs_error=float(np.max(abs((con(y0+h*d)[0]-con(y0-h*d)[0])/(2*h)-J@d))))
            atomic_json(dest/'derivative_check.json',check)
            if max(check.values())>1e-4:raise ValueError(check)
        best=[np.inf,None]
        class BudgetStop(Exception):pass
        def cb(y):
            last[0]=y.copy();cc=T@y[:-1];sv=o.constraints(cc)[0];pv=e.constraints(cc,y[-1])[0]
            feas=float(min(sv.min(),pv.min(),((cc-lo)*bs).min(),((hi-cc)*bs).min(),y[-1],4-y[-1]))
            score=fun(y)[0];row=dict(iteration=len(hist)+1,objective=score,epigraph_bound=float(np.sqrt(max(0,y[-1])*pscale)),min_constraint=feas,elapsed=time.monotonic()-start)
            hist.append(row)
            if feas>=-1e-7:
                try: raw=m.candidate(cc)
                except ValueError:raw=None
                if raw is not None and score<best[0]:
                    best[:]=[score,y.copy()];m.f.save(raw,dest/'best_feasible.json',row);np.save(dest/'best_modifiers.npy',cc)
            if len(hist)%5==0:
                atomic_json(dest/'history.json',hist);np.save(dest/'checkpoint.npy',cc);print(json.dumps(row),flush=True)
            if time.monotonic()-start>seconds:raise BudgetStop
        try:
            res=minimize(fun,y0,jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:con(y)[0],jac=lambda y:con(y)[1])],callback=cb,options=dict(maxiter=maxiter,ftol=2e-9,disp=False))
            fin=res.x;success=bool(res.success);message=str(res.message);nfev=int(res.nfev)
        except BudgetStop:
            fin=last[0];success=False;message='registered wall budget; not converged';nfev=counts['fun']
        y=best[1] if best[1] is not None else fin;cc=T@y[:-1]
        vscan=scan(m,cc,pool);peak=float(vscan.max());feas=float(o.constraints(cc)[0].min())
        summary=dict(round=rr+1,iterations=len(hist),nfev=nfev,counts=counts,success=success,message=message,active_points=len(indices),training_pool_max=peak,epigraph_bound=float(np.sqrt(max(0,y[-1])*pscale)),pool_bound_gap=peak-float(np.sqrt(max(0,y[-1])*pscale)),min_structure_constraint=feas,chosen_feasible=best[1] is not None,elapsed=time.monotonic()-start)
        m.f.save(m.candidate(cc),dest/'candidate.json',summary);np.save(dest/'modifiers.npy',cc);atomic_json(dest/'summary.json',summary);atomic_json(dest/'history.json',hist)
        print('ROUND',json.dumps(summary),flush=True);report.append(summary)
        if best[1] is not None and feas>=-1e-7 and peak<selected_peak:
            selected_c=cc.copy();selected_peak=peak
            m.f.save(m.candidate(cc),out/'candidate.json',dict(selection='best feasible full training-pool peak',**summary));np.save(out/'modifiers.npy',cc)
        c=cc;values=vscan;del e;gc.collect()
    if not (out/'candidate.json').exists():
        m.f.save(m.candidate(selected_c),out/'candidate.json',dict(selection='no qualifying improvement; parent retained'));np.save(out/'modifiers.npy',selected_c)
    result=dict(**reg,actual_pool_points=len(pool),variables=m.dim,parent_pool_peak=parent_peak,selected_pool_peak=selected_peak,round_reports=report,raw_sha256=hashlib.sha256((out/'candidate.json').read_bytes()).hexdigest(),new_improvement_found=selected_peak<parent_peak)
    atomic_json(out/'summary.json',result);print('DONE',json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--rounds',type=int,default=2);p.add_argument('--maxiter',type=int,default=170);p.add_argument('--seconds',type=float,default=540);a=p.parse_args();run(a.parent,a.out,a.rounds,a.maxiter,a.seconds)
