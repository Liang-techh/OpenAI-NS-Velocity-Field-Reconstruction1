"""ST055: swirl-only full-momentum correction, with frozen poloidal u, p and f.
The temporal correction vanishes at t=.25. Original energy normalization and
all poloidal kinematics remain unchanged, but radial centrifugal residual is
included exactly. No theta-only loss is reported as full NS acceptance.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh, null_space
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st054'))
from pressure_completion import Family, make_pool, save
from spacetime import NU

def temporal(t, nt):
    t=np.asarray(t).ravel(); k=np.arange(1,nt)
    S=np.column_stack([Chebyshev.basis(int(j))(4*t-2)-(-1.)**j for j in k])
    Sd=np.column_stack([4*Chebyshev.basis(int(j)).deriv()(4*t-2) for j in k])
    return S,Sd

class SwirlModel:
    def __init__(self,parent):
        self.f,self.raw=Family.load(parent);f=self.f
        self.a,self.b,self.p,self.fc,self.amp=f.coefficients(self.raw)
        self.dim=f.ns*(f.nt-1)
        self.q0=f.q0[1:]
    def delta_coefficients(self,c):
        f=self.f;c=np.asarray(c).reshape(f.ns,f.nt-1)
        return np.column_stack((-c@self.q0,c)).ravel()
    def candidate(self,c):
        f=self.f;raw=self.raw.copy()
        raw[f.n:2*f.n]+=self.delta_coefficients(c)/self.amp
        if not np.isfinite(raw).all() or np.max(abs(raw[f.n:2*f.n]))>4+1e-10:
            raise ValueError('Swirl coefficient bound violation')
        return raw
    def cache(self, points):
        # points are s=r^2,z,t; use original exact derivative bundles in chunks.
        points=np.asarray(points,float);n=len(points);f=self.f
        M=np.empty((n,self.dim));L=np.empty_like(M);base=np.empty((n,3));B=np.empty(n)
        for start in range(0,n,128):
            sl=slice(start,start+128);s,z,t=points[sl].T;r=np.sqrt(s)
            D=f.bundle(s,z,t);A=D['A']@self.a;C=D['C']@self.a
            wb=f.basis(s,z,f.nr,f.nz,False)
            W=wb[0,0]@f.Tw;Ws=wb[1,0]@f.Tw;Wz=wb[0,1]@f.Tw
            WL=(8*wb[1,0]+4*s[:,None]*wb[2,0]+wb[0,2])@f.Tw
            S,Sd=temporal(t,f.nt)
            K=2*A[:,None]*W+(2*s*A)[:,None]*Ws+C[:,None]*Wz-NU*WL
            M[sl]=(W[:,:,None]*S[:,None,:]).reshape(-1,self.dim)
            L[sl]=(r[:,None,None]*(W[:,:,None]*Sd[:,None,:]+K[:,:,None]*S[:,None,:])).reshape(-1,self.dim)
            B[sl]=D['B']@self.b
            rr=D['At']@self.a+A*A+2*s*A*(D['As']@self.a)+C*(D['Az']@self.a)-B[sl]**2+2*(D['Qs']@self.p)-NU*(D['AL']@self.a)
            rt=D['Bt']@self.b+2*A*B[sl]+2*s*A*(D['Bs']@self.b)+C*(D['Bz']@self.b)-NU*(D['BL']@self.b)
            rz=D['Ct']@self.a+2*s*A*(D['Cs']@self.a)+C*(D['Cz']@self.a)+(D['Qz']@self.p)-NU*(D['CL']@self.a)
            from spacetime import force
            base[sl]=np.c_[r*rr,r*rt,rz]-force(np.c_[r,np.zeros(len(r)),z],t,*self.fc)
        return dict(points=points,M=M,L=L,B=B,R=base,r=np.sqrt(points[:,0]))
    def residual(self,c,D):
        b=D['M']@c;R=D['R'].copy()
        R[:,0]-=D['r']*(2*D['B']*b+b*b)
        R[:,1]+=D['L']@c
        return R,b
    def pullback(self,c,D,weights):
        R,b=self.residual(c,D)
        g=D['M'].T@(-2*D['r']*(D['B']+b)*weights*R[:,0])+D['L'].T@(weights*R[:,1])
        return R,g

class Objective:
    def __init__(self,m,peak_weight=.025,core_fraction=.03):
        self.m=m;self.peak_weight=peak_weight;self.core_fraction=core_fraction
        x,wx=leggauss(32);z,wz=leggauss(48);tg,tw=leggauss(13)
        times=np.r_[.25,.5+.25*tg,.75];timew=np.r_[.05,.45*tw,.05]
        S,Z,T=np.meshgrid(2*(x+1),2*z,times,indexing='ij')
        q=np.c_[S.ravel(),Z.ravel(),T.ravel()]
        weights=(np.outer(wx,wz)[:,:,None]*timew[None,None,:]/4).ravel()
        pool=make_pool(9175550)
        r,z,t=np.meshgrid(np.linspace(.025,.65,13),np.linspace(-.28,.28,11),np.linspace(.25,.75,13),indexing='ij')
        pool=np.r_[pool,np.c_[(r*r).ravel(),z.ravel(),t.ravel()]]
        self.nq=len(q);self.npool=len(pool);self.D=m.cache(np.r_[q,pool]);self.w=np.r_[weights,np.zeros(len(pool))]
        R,Z,T=np.meshgrid(np.linspace(.035,.215,7),np.r_[-np.linspace(.025,.215,7)[::-1],np.linspace(.025,.215,7)],np.linspace(.25,.75,19),indexing='ij')
        self.C=m.cache(np.c_[(R*R*(1-T)).ravel(),(Z*(1-T)**.495).ravel(),T.ravel()])
        self.core_shape=R.shape
        self.cscale=np.maximum(abs(self.C['B']),.01)
        self.times=np.linspace(.25,.75,25);tau=1-self.times
        self.G=m.cache(np.c_[.01*tau,.1*tau**.495,self.times])
        pts=np.c_[.1*np.sqrt(tau),np.zeros(len(tau)),.1*tau**.495]
        self.gvel=m.f.fields(m.raw,pts,self.times)[0]*np.c_[np.sqrt(tau),tau**.505,tau**.505]
        self.gscale=self.G['r']*tau**.505;self.gnorm=np.linalg.norm(self.gvel[0])
        self.last=None;self.calls=0
    def value(self,c):
        if self.last is not None and np.array_equal(c,self.last[0]):return self.last[1:]
        R,b=self.m.residual(c,self.D);sq=np.sum(R*R,axis=1)
        tau=4e-5;sn=sq[self.nq:]/tau;ls=logsumexp(sn)
        weights=self.w.copy();weights[self.nq:]=self.peak_weight*np.exp(sn-ls)
        loss=.5*float(self.w@sq)+.5*self.peak_weight*tau*(ls-np.log(self.npool))
        g=self.D['M'].T@(-2*self.D['r']*(self.D['B']+b)*weights*R[:,0])+self.D['L'].T@(weights*R[:,1])
        # Preserve core swirl within a declared relative band. Final exact training
        # feasibility is checked separately, not inferred from a penalty value.
        delta=self.C['M']@c;ratio=delta/self.cscale
        v=np.maximum(abs(ratio)-self.core_fraction,0.)
        loss+=.5*.2*np.mean(v*v)
        g+=.2/len(v)*(self.C['M'].T@(v*np.sign(ratio)/self.cscale))
        # Original single-probe scaled-profile condition, tested more often in time.
        gv=self.gvel.copy();gv[:,1]+=self.gscale*(self.G['M']@c)
        diff=gv-gv[0];dn=np.sqrt(np.sum(diff*diff,axis=1)+1e-30);drift=dn/self.gnorm
        excess=np.maximum(drift-.049,0.)
        loss+=.5*np.mean(excess**2)
        g+=(self.G['M'].T@(excess*diff[:,1]/(dn*self.gnorm)*self.gscale))/len(drift)
        # Small physical-swirl change regularization; not an amplitude relaxation.
        loss+=.5*.002*float(np.dot(c,c));g+=.002*c
        # Explicit raw-coefficient limit penalty, independently checked on freezing.
        cb=self.m.raw[self.m.f.n:2*self.m.f.n]+self.m.delta_coefficients(c)/self.m.amp
        ex=np.maximum(abs(cb)-3.999,0.);gx=ex*np.sign(cb)/self.m.amp
        loss+=.5*np.dot(ex,ex)
        g+=(gx.reshape(self.m.f.ns,self.m.f.nt)[:,1:]-gx.reshape(self.m.f.ns,self.m.f.nt)[:,0,None]*self.m.q0).ravel()
        self.calls+=1;self.last=(c.copy(),float(loss),g);return float(loss),g
    def stats(self,c):
        R,b=self.m.residual(c,self.D);sq=np.sum(R*R,axis=1)
        vel=self.gvel.copy();vel[:,1]+=self.gscale*(self.G['M']@c)
        drift=np.linalg.norm(vel-vel[0],axis=1)/self.gnorm
        ratio=(self.C['M']@c)/self.cscale
        bounds=True
        try:self.m.candidate(c)
        except ValueError:bounds=False
        return dict(training_mse=float(self.w@sq),pool_max=float(np.sqrt(sq[self.nq:]).max()),pool_theta_max=float(abs(R[self.nq:,1]).max()),core_relative_swirl_change=float(abs(ratio).max()),core_swirl_min=float((self.C['B']+self.C['M']@c).min()),core_drift_max=float(drift.max()),bounds_ok=bounds)
    def precondition(self):
        # A fixed TRAINING-only subsampled Jacobian metric; no use of validation.
        idx=np.arange(0,self.nq,7);D=self.D;Jr=(-2*D['r'][idx]*D['B'][idx])[:,None]*D['M'][idx]
        w=self.w[idx]*7
        H=D['L'][idx].T@(w[:,None]*D['L'][idx])+Jr.T@(w[:,None]*Jr)+.002*np.eye(self.m.dim)
        sc=1/np.sqrt(np.maximum(np.diag(H),1e-10));H=H*sc[:,None]*sc[None,:]
        vals,V=eigh(H,check_finite=False);floor=1e-5*max(float(vals[-1]),1.)
        T=sc[:,None]*(V/np.sqrt(np.maximum(vals,floor))[None,:])
        return T,dict(min_eigenvalue=float(vals[0]),max_eigenvalue=float(vals[-1]),floor=floor)

def run(parent,out,ident,iterations=250,seconds=360,freeze_core=False):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Nonempty output directory')
    out.mkdir(parents=True,exist_ok=True)
    registration=dict(id=ident,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),physical_contract='UNCHANGED',poloidal_pressure_force='FROZEN',temporal_initial_value='EXACTLY PRESERVED by T_k(tau)-T_k(-1)',training_quad=[32,48,13],pool_seed=9175550,holdout_seeds=[9175591,9175592],core_swirl_band=.03,original_core_drift_training_cap=.049,peak_weight=.025,peak_temperature=4e-5,velocity_regularization=.002,maxiter=iterations,wall_seconds=seconds,selection='lowest feasible TRAINING objective; never select on holdouts',registered_utc=datetime.now(timezone.utc).isoformat(),pde_validated=False,freeze_core_swirl_at_25_times=freeze_core)
    save(out/'registration.json',registration);start=time.monotonic();m=SwirlModel(parent);o=Objective(m)
    T,inf=o.precondition()
    if freeze_core:
        Z=null_space(o.G['M']@T,rcond=1e-10);T=T@Z
        inf.update(free_solver_variables=T.shape[1],core_annihilation_max=float(abs(o.G['M']@T).max()))
    save(out/'conditioning.json',inf)
    c0=np.zeros(m.dim);base=o.stats(c0);basecost=o.value(c0)[0];save(out/'parent_training.json',base)
    rng=np.random.default_rng(9175551);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=1e-8
    v,g=o.value(c0);err=abs((o.value(h*d)[0]-o.value(-h*d)[0])/(2*h)-g@d)
    save(out/'gradient_check.json',dict(absolute_error=err,step=h,note='Core drift hinge is C1 but not C2 at parent; prior separate step refinement confirms derivative convergence.'));print('BASE',base,'GRAD',err,flush=True)
    if err>1e-8:raise ValueError('Objective derivative calibration failed')
    history=[];best=[basecost,c0.copy()];last=[c0.copy()]
    def feasible(st):
        return st['bounds_ok'] and st['core_drift_max']<=.0495 and st['core_relative_swirl_change']<=.03001 and st['core_swirl_min']>0 and st['pool_max']<=base['pool_max']*1.0005 and st['training_mse']<=base['training_mse']
    class BudgetStop(Exception):pass
    def fun(y):
        v,g=o.value(T@y);return v,T.T@g
    def callback(y):
        c=T@y;last[0]=c.copy();v=o.value(c)[0]
        st=o.stats(c);row=dict(iteration=len(history)+1,objective=v,**st,elapsed=time.monotonic()-start)
        history.append(row)
        if feasible(st) and v<best[0]:
            best[:]=[v,c.copy()];m.f.save(m.candidate(c),out/'best_feasible.json',row);np.save(out/'best_modifiers.npy',c)
        if len(history)%5==0:
            save(out/'history.json',history);np.save(out/'checkpoint.npy',c);print('ITER',json.dumps(row),flush=True)
        if time.monotonic()-start>seconds:raise BudgetStop
    try:
        ret=minimize(fun,np.zeros(T.shape[1]),jac=True,method='L-BFGS-B',callback=callback,options=dict(maxiter=iterations,maxls=30,ftol=1e-13,gtol=1e-9,maxcor=20))
        end=T@ret.x;ok=bool(ret.success);msg=str(ret.message)
    except BudgetStop:end=last[0];ok=False;msg='registered wall budget; checkpoint retained'
    # A predeclared training-only recovery along the final direction. Retain all
    # tested values; never inspect a holdout when choosing this scale.
    line=[]
    for alpha in (1.,.8,.6,.4,.2,.1,.05):
        c=alpha*end;st=o.stats(c);v=o.value(c)[0];line.append(dict(alpha=alpha,objective=v,**st,feasible=feasible(st)))
        if feasible(st) and v<best[0]:best[:]=[v,c.copy()]
    c=best[1];raw=m.candidate(c);np.save(out/'swirl_delta.npy',m.delta_coefficients(c)/m.amp);np.save(out/'modifiers.npy',c)
    m.f.save(raw,out/'candidate.json',dict(id=ident,selection='feasible training iterate or registered scaled final direction',poloidal_pressure_force_unchanged=True))
    summary=dict(**registration,variables=m.dim,solver_variables=T.shape[1],iterations=len(history),function_evaluations=o.calls,optimizer_success=ok,message=msg,initial_objective=basecost,selected_objective=best[0],parent_training=base,child_training=o.stats(c),elapsed=time.monotonic()-start,nonzero_change=bool(np.any(c)),initial_B_coefficient_error=float(np.max(abs((raw[m.f.n:2*m.f.n]-m.raw[m.f.n:2*m.f.n]).reshape(m.f.ns,m.f.nt)@m.f.q0))))
    save(out/'summary.json',summary);save(out/'history.json',history);save(out/'final_line.json',line);print('SUMMARY',json.dumps(summary),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',required=True);p.add_argument('--iterations',type=int,default=250);p.add_argument('--seconds',type=float,default=360);p.add_argument('--freeze-core',action='store_true');a=p.parse_args();run(a.parent,a.out,a.id,a.iterations,a.seconds,a.freeze_core)
