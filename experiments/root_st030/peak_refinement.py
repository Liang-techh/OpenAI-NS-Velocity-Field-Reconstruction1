"""ST033: fresh training-pool residual refinement and a smooth peak penalty.
Full Cartesian validation remains independent and retains the original gates.
"""
from __future__ import annotations
import argparse,json,hashlib,time,gc
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from spacetime import Family,force
from harmonic_fit import Objective

class PeakObjective(Objective):
    def __init__(self,f,points,sigma=.05):
        n=len(points);super().__init__(f,seed=9172904,count=n,moment_weight=10.)
        del self.D;gc.collect()
        s,z,t=points.T;self.D=f.bundle(s,z,t);self.s=s;self.sqrt_s=np.sqrt(s)
        xyz=np.column_stack((np.sqrt(s),np.zeros(n),z))
        self.FA=force(xyz,t,1,0);self.FC=force(xyz,t,0,1);self.sigma=sigma
    def transformed(self,x):
        self.evaluate(x);r=self.cache_r;J=self.cache_j;N=self.N
        R=r[:3*N].reshape(3,N);JG=J[:3*N].reshape(3,N,-1)
        tail=np.sqrt(N)*np.sum(R*R,axis=0)/self.sigma
        JT=(2*np.sqrt(N)/self.sigma)*np.einsum('an,ani->ni',R,JG,optimize=True)
        return np.r_[r,tail],np.vstack((J,JT))


def scan(f,x,pool):
    values=[]
    for lo in range(0,len(pool),512):
        s,z,t=pool[lo:lo+512].T;p=np.column_stack((np.sqrt(s),np.zeros(len(s)),z))
        values.append(np.linalg.norm(f.analytic_residual(x,p,t),axis=1))
    return np.concatenate(values)


def fit(f,x,points,out,iterations=20,sigma=.05):
    o=PeakObjective(f,points,sigma);hist=[];damping=1e-4;start=time.time()
    low=np.r_[np.full(2*f.n,-4.),np.full(f.n,-100.),0.,0.];high=-low;high[-2:]=10.
    for it in range(iterations):
        r,J=o.transformed(x);cost=float(r@r);g=J.T@r;sc=1/np.maximum(np.linalg.norm(J,axis=0),1e-10)
        K=J*sc;H=K.T@K;gg=K.T@r;ok=False;new=cost;dx=np.zeros_like(x)
        del K
        for attempt in range(7):
            M=H.copy();M.flat[::len(x)+1]+=damping
            try:dx=-sc*cho_solve(cho_factor(M,lower=True,check_finite=False),gg,check_finite=False)
            except np.linalg.LinAlgError:damping*=10;continue
            alpha=1.
            for line in range(10):
                y=np.clip(x+alpha*dx,low,high)
                rr,_=o.transformed(y);new=float(rr@rr)
                amp=f.coefficients(y)[4]
                if np.isfinite(new) and new<cost and 1e-4<=amp<=100:
                    ok=True;x=y;break
                alpha*=.5
            if ok:damping=max(1e-10,damping*.4);break
            damping*=10
        row=dict(iteration=it,loss=cost,new_loss=new if ok else cost,accepted=ok,gradient_inf=float(np.max(np.abs(g))),elapsed=time.time()-start,damping=damping)
        hist.append(row);f.save(x,out/'candidate.json',dict(stage='training-only peak refinement',round_iteration=it+1))
        (out/'history.json').write_text(json.dumps(hist,indent=2)+'\n');print(json.dumps(row),flush=True)
        if not ok:break
    return x,hist


def run(warm,out,rounds=2,iterations=20):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    reg=dict(experiment='CR-ROOT-ST033',created_before_fit=True,initial_sha256=hashlib.sha256(Path(warm).read_bytes()).hexdigest(),training_seed=9172904,pool_seed=9172905,validation_seed=9172912,rounds=rounds,iterations_per_round=iterations,sigma=.05,objective='mean(|R|^2)+mean(|R|^4)/sigma^2 plus unchanged core/energy/harmonic/torque penalties',physical_gates='UNCHANGED',selection='Final bounded training iterate; no validation selection',training_count=4096,pde_validated=False)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    f,x=Family.load(warm);rng=np.random.default_rng(9172904)
    fixed=np.column_stack((rng.uniform(0,4,3072),rng.uniform(-2,2,3072),rng.uniform(.25,.75,3072)))
    fixed[:384,0]=rng.uniform(0,.5,384)**2;fixed[:384,1]=rng.uniform(-.5,.5,384)
    fixed[384:640,2]=.25;fixed[640:896,2]=.75
    rp=np.random.default_rng(9172905)
    rnd=np.column_stack((rp.uniform(0,4,4096),rp.uniform(-2,2,4096),rp.uniform(.25,.75,4096)))
    r=np.unique(np.r_[np.linspace(0,1.995,41),np.linspace(1.6,1.995,17)])
    z=np.unique(np.r_[np.linspace(-1.995,1.995,49),np.linspace(-1.995,-1.6,9),np.linspace(1.6,1.995,9)])
    S,Z,T=np.meshgrid(r*r,z,[.25,.34375,.5,.65625,.75],indexing='ij')
    pool=np.vstack((rnd,np.column_stack((S.ravel(),Z.ravel(),T.ravel()))));np.savez_compressed(out/'training_pool.npz',pool=pool,fixed=fixed)
    scores=[]
    def pool_score(stage):
        v=scan(f,x,pool);row=dict(stage=stage,pool_max=float(v.max()),uniform_spacetime_volume_scaled_RMS=float(np.sqrt(16*np.pi*np.mean(v[:4096]**2))),scope='TRAINING pool, not independent validation');scores.append(row);(out/'pool_scores.json').write_text(json.dumps(scores,indent=2)+'\n');print(json.dumps(row),flush=True);return v
    v=pool_score('initial');f.save(x,out/'initial.json')
    for k in range(rounds):
        chosen=np.argsort(v)[-1024:];points=np.vstack((fixed,pool[chosen]));dest=out/f'round{k+1}';dest.mkdir(exist_ok=True);np.save(dest/'training_points.npy',points)
        x,h=fit(f,x,points,dest,iterations=iterations);gc.collect();v=pool_score(f'round{k+1}')
    f.save(x,out/'candidate.json',dict(experiment='CR-ROOT-ST033',selection='last training iterate, no holdout selection'))
    (out/'summary.json').write_text(json.dumps(dict(**reg,pool_points=len(pool),pool_scores=scores,complete=True),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--warm',required=True);p.add_argument('--out',required=True);p.add_argument('--rounds',type=int,default=2);p.add_argument('--iterations',type=int,default=20);a=p.parse_args();run(a.warm,a.out,a.rounds,a.iterations)
