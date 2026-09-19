"""ST054: full existing compact pressure space, fixed velocity and force.
The objective is convex in pressure; grids and optimality are finite dimensional.
No new velocity, force, physical support, or acceptance definition is introduced.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky, cho_solve, solve_triangular
from scipy.optimize import minimize
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st051'))
import aligned_continuation
from spacetime import Family

def save(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)

def pressure_basis(f,s,z,t=None):
    s,z=np.broadcast_arrays(s,z);s,z=s.ravel(),z.ravel()
    B=f.basis(s,z,f.nr,f.nz,False)
    Pr=2*np.sqrt(s)[:,None]*(B[1,0]@f.Tq);Pz=B[0,1]@f.Tq
    if t is None:return Pr,Pz
    T=np.column_stack([Chebyshev.basis(k)(4*np.broadcast_to(t,s.shape)-2) for k in range(f.nt)])
    return tuple((X[:,:,None]*T[:,None,:]).reshape(len(s),f.n) for X in (Pr,Pz))

def residual(f,raw,points,t):
    return np.concatenate([f.analytic_residual(raw,points[i:i+256],np.broadcast_to(t,(len(points),))[i:i+256]) for i in range(0,len(points),256)])

def make_pool(seed):
    rng=np.random.default_rng(seed)
    r=np.unique(np.r_[np.linspace(0,1.99,27),np.linspace(.15,.85,9)])
    z=np.unique(np.r_[np.linspace(-1.99,1.99,43),np.linspace(-1.98,-1.7,9),np.linspace(1.7,1.98,9)])
    R,Z,T=np.meshgrid(r,z,np.linspace(.25,.75,11),indexing='ij')
    grid=np.c_[R.ravel()**2,Z.ravel(),T.ravel()]
    rnd=np.c_[rng.uniform(0,4,3072),rng.uniform(-2,2,3072),rng.uniform(.25,.75,3072)]
    return np.r_[grid,rnd]

def scan(f,raw,pool):
    s,z,t=pool.T;x=np.c_[np.sqrt(s),np.zeros(len(s)),z]
    return residual(f,raw,x,t)

class PressureFit:
    def __init__(self,parent,order=(40,60),ntimes=17):
        self.f,self.raw=Family.load(parent);f=self.f;self.q=self.raw[2*f.n:3*f.n].copy()
        x,wx=leggauss(order[0]);z,wz=leggauss(order[1])
        S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij');s,z=S.ravel(),Z.ravel()
        self.w=(np.outer(wx,wz)/4).ravel();self.s,self.z=s,z
        self.pr,self.pz=pressure_basis(f,s,z)
        g,gw=leggauss(ntimes);self.times=np.r_[.25,.5+.25*g,.75];self.tw=np.r_[.05,.45*gw,.05]
        self.T=np.column_stack([Chebyshev.basis(k)(4*self.times-2) for k in range(f.nt)])
        K=self.pr.T@(self.w[:,None]*self.pr)+self.pz.T@(self.w[:,None]*self.pz)
        self.G=np.kron(K,self.T.T@(self.tw[:,None]*self.T))
        self.L=cholesky(self.G,lower=True);self.W=solve_triangular(self.L.T,np.eye(f.n),lower=False)
        self.R=[];rhs=np.zeros((f.ns,f.nt));cost=0.;theta=0.
        points=np.c_[np.sqrt(s),np.zeros(len(s)),z]
        for j,t in enumerate(self.times):
            R=residual(f,self.raw,points,t);self.R.append(R)
            g=self.pr.T@(self.w*R[:,0])+self.pz.T@(self.w*R[:,2])
            rhs+=self.tw[j]*g[:,None]*self.T[j]
            cost+=self.tw[j]*float(self.w@np.sum(R*R,axis=1));theta+=self.tw[j]*float(self.w@(R[:,1]**2))
        self.g=rhs.ravel();self.cost=cost;self.theta_cost=theta
        self.unconstrained=-cho_solve((self.L,True),self.g)
    def candidate(self,dq):
        r=self.raw.copy();r[2*self.f.n:3*self.f.n]+=dq
        if not np.isfinite(r).all() or np.max(abs(r[2*self.f.n:3*self.f.n]))>100+1e-8:raise ValueError('Pressure bounds')
        return r
    def metrics(self,dq):
        vals=[]
        q=dq.reshape(self.f.ns,self.f.nt)
        for j,R in enumerate(self.R):
            D=q@self.T[j];v=R.copy();v[:,0]+=self.pr@D;v[:,2]+=self.pz@D
            vals.append(float(self.w@np.sum(v*v,axis=1)))
        return dict(weighted_mse=float(self.tw@vals),time_volume_L2=(np.sqrt(16*np.pi*np.array(vals))).tolist(),theta_weighted_mse=self.theta_cost)
    def pressure_constraints(self):
        R,Z,T=np.meshgrid(np.linspace(.03,.22,9),np.r_[-np.linspace(.02,.22,9)[::-1],np.linspace(.02,.22,9)],np.linspace(.25,.75,19),indexing='ij')
        s=(R*R*(1-T)).ravel();z=(Z*(1-T)**.495).ravel();t=T.ravel()
        Pr,Pz=pressure_basis(self.f,s,z,t);Pz*=np.sign(z)[:,None]
        # Zero is the previously declared direction, no new positive margin.
        return np.vstack((Pr,Pz))

def run(parent,out,ident,seconds=300,iterations=150,peak_ratio=1.,margin=0.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Nonempty output directory')
    reg=dict(id=ident,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),original_contract='UNCHANGED',velocity_force='FROZEN',pressure_space='ALL existing spatial/time compact columns',spatial_order=[40,60],temporal_gauss=17,endpoint_weights=.05,training_pool_seed=9175450,validation_seeds=[9175491,9175492],selection='final constrained iterate before all independent validation, no holdout tuning',wall_budget_seconds=seconds,iteration_budget_per_exchange=iterations,peak_cap_ratio=peak_ratio,peak_rule='max(ratio * parent pool max, 1.005 * parent azimuthal max)',pressure_margin=margin,feasibility_tolerance=1e-9,exchange_rounds=2,pressure_directions='radial and signed axial >=0 on new declared core training grid',registered_utc=datetime.now(timezone.utc).isoformat(),pde_validated=False)
    save(out/'registration.json',reg);start=time.monotonic();o=PressureFit(parent);f=o.f
    info=dict(parent=o.metrics(np.zeros(f.n)),unconstrained=o.metrics(o.unconstrained),coefficient_count=f.n,stationarity=float(np.max(abs(o.G@o.unconstrained+o.g))),elapsed=time.monotonic()-start)
    f.save(o.raw,out/'parent.json',dict(scope='frozen parent copied'))
    f.save(o.candidate(o.unconstrained),out/'unconstrained.json',dict(scope='unconstrained pressure diagnostic, NOT claimed structural success'))
    save(out/'diagnostic.json',info);print('DIAGNOSTIC',json.dumps(info),flush=True)
    pool=make_pool(9175450);baseR=scan(f,o.raw,pool);baseNorm=np.linalg.norm(baseR,axis=1);cap=max(peak_ratio*float(baseNorm.max()),1.005*float(np.max(abs(baseR[:,1]))));np.savez_compressed(out/'training_pool.npz',pool=pool,parent_norm=baseNorm)
    C=o.pressure_constraints();C0=C@o.q-margin;CW=C@o.W;scale=np.maximum(np.linalg.norm(CW,axis=1),1e-4)
    CW/=scale[:,None];C0/=scale
    reg.update(training_pool_count=len(pool),training_pool_parent_max=float(baseNorm.max()),training_pool_azimuthal_max=float(np.max(abs(baseR[:,1]))),training_peak_cap=cap,initial_pressure_min=float(C0.min()));save(out/'registration.json',reg)
    active=np.unique(np.r_[np.argsort(baseNorm)[-384:],np.linspace(0,len(pool)-1,64,dtype=int)])
    y=np.zeros(f.n);g=o.W.T@o.g;H=o.W.T@o.G@o.W
    history=[];best=None;bestcost=np.inf
    for exchange in range(2):
        best=None;bestcost=np.inf
        s,z,t=pool[active].T;pr,pz=pressure_basis(f,s,z,t);AR=pr@o.W;AZ=pz@o.W;BR=baseR[active]
        def obj(y):return float(.5*y@H@y+g@y),(H@y+g)
        def constraints(y):
            rr=BR[:,0]+AR@y;rz=BR[:,2]+AZ@y
            peak=(cap*cap-rr*rr-rz*rz-BR[:,1]**2)/(2*cap)
            Jpeak=-(rr[:,None]*AR+rz[:,None]*AZ)/cap
            return np.r_[C0+CW@y,peak],np.vstack((CW,Jpeak))
        def callback(v):
            nonlocal y,best,bestcost
            y=v.copy();cons=constraints(y)[0];row=dict(exchange=exchange,iteration=len(history)+1,objective=obj(y)[0],min_constraint=float(cons.min()),elapsed=time.monotonic()-start)
            history.append(row)
            if cons.min()>=-1e-9 and row['objective']<bestcost:
                best=y.copy();bestcost=row['objective'];f.save(o.candidate(o.W@y),out/'best_feasible.json',row)
            if len(history)%5==0:
                f.save(o.candidate(o.W@y),out/'checkpoint.json',row);save(out/'history.json',history);print('ITER',json.dumps(row),flush=True)
            if time.monotonic()-start>seconds:raise TimeoutError('Registered wall budget')
        try:
            ret=minimize(obj,y,method='SLSQP',jac=True,constraints=[{'type':'ineq','fun':lambda v:constraints(v)[0],'jac':lambda v:constraints(v)[1]}],callback=callback,options=dict(maxiter=iterations,ftol=1e-13))
            y=ret.x;status={'success':bool(ret.success),'message':str(ret.message),'iterations':int(ret.nit)}
        except TimeoutError:
            status={'success':False,'message':'registered wall budget'}
        if best is not None:y=best.copy()
        newR=scan(f,o.candidate(o.W@y),pool);norm=np.linalg.norm(newR,axis=1)
        save(out/f'exchange{exchange}.json',dict(**status,active_count=len(active),pool_max=float(norm.max()),cap=cap,min_core_pressure_constraint=float((C0+CW@y).min()),training=o.metrics(o.W@y)))
        print('EXCHANGE',exchange,status,float(norm.max()),flush=True)
        active=np.union1d(active,np.argsort(norm)[-384:])
        if time.monotonic()-start>seconds:break
    final=o.candidate(o.W@y)
    f.save(final,out/'candidate.json',dict(id=ident,velocity_force_unchanged=True,scientific_acceptance=False))
    np.save(out/'pressure_delta.npy',o.W@y)
    summary=dict(reg | status,total_iterations=len(history),elapsed=time.monotonic()-start,training=o.metrics(o.W@y),pool_max=float(norm.max()),core_pressure_min=float((C0+CW@y).min()),raw_velocity_identical=bool(np.array_equal(final[:2*f.n],o.raw[:2*f.n])),raw_force_identical=bool(np.array_equal(final[-2:],o.raw[-2:])),pressure_max=float(np.max(abs(final[2*f.n:3*f.n]))))
    save(out/'summary.json',summary);save(out/'history.json',history);print('SUMMARY',json.dumps(summary),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',required=True);p.add_argument('--seconds',type=int,default=300);p.add_argument('--iterations',type=int,default=150);p.add_argument('--peak-ratio',type=float,default=1.);p.add_argument('--margin',type=float,default=0.);a=p.parse_args();run(a.parent,a.out,a.id,a.seconds,a.iterations,a.peak_ratio,a.margin)
