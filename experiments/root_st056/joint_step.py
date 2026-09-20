"""ST056: joint full-pressure / initial-flat swirl correction.
New numerical experiment; unchanged original physical and acceptance contract.
"""
from __future__ import annotations
import argparse,hashlib,json,sys,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cholesky,solve_triangular,svd
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st030'))
from spacetime import Family,force,NU

def save(path,obj):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)

def tensor(space,temporal):
    return (space[:,:,None]*temporal[:,None,:]).reshape(len(space),-1)

class JointModel:
    def __init__(self,path):
        self.f,self.raw=Family.load(path);f=self.f
        self.a,self.b,self.p,self.fc,self.amp=f.coefficients(self.raw)
        self.nb=f.ns*(f.nt-1);self.np=f.n;self.dim=self.nb+self.np
    def delta(self,c):
        c=np.asarray(c);f=self.f;sb=c[:self.nb].reshape(f.ns,f.nt-1)
        db=np.column_stack((-sb@f.q0[1:],sb)).ravel()
        return db,c[self.nb:]
    def candidate(self,c):
        f=self.f;db,dp=self.delta(c);raw=self.raw.copy()
        raw[f.n:2*f.n]+=db/self.amp;raw[2*f.n:3*f.n]+=dp
        if not np.isfinite(raw).all():raise ValueError('Nonfinite field')
        if np.max(abs(raw[:2*f.n]))>4+1e-10 or np.max(abs(raw[2*f.n:-2]))>100+1e-10:raise ValueError('Raw coefficient bound')
        return raw
    def cache(self,points):
        p=np.asarray(points,float);n=len(p);f=self.f
        D={'points':p,'M':np.empty((n,self.nb)),'L':np.empty((n,self.nb)),
           'Pr':np.empty((n,self.np)),'Pz':np.empty((n,self.np)),
           'B':np.empty(n),'R':np.empty((n,3)),'u':np.empty((n,3)),'r':np.sqrt(p[:,0])}
        for j in range(0,n,256):
            sl=slice(j,j+256);s,z,t=p[sl].T;r=np.sqrt(s)
            P=f.basis(s,z,f.nr,f.nz,True);W=f.basis(s,z,f.nr,f.nz,False)
            T=np.column_stack([Chebyshev.basis(k)(4*t-2) for k in range(f.nt)])
            Td=np.column_stack([4*Chebyshev.basis(k).deriv()(4*t-2) for k in range(f.nt)])
            S=T[:,1:]-f.q0[1:];St=Td[:,1:]
            ca=self.a.reshape(f.ns,f.nt)@T.T;cat=self.a.reshape(f.ns,f.nt)@Td.T
            cb=self.b.reshape(f.ns,f.nt)@T.T;cbt=self.b.reshape(f.ns,f.nt)@Td.T
            cp=self.p.reshape(f.ns,f.nt)@T.T
            evalmat=lambda X,C:np.einsum('ij,ji->i',X,C)
            PA=-P[0,1]@f.Tp;PC=(2*P[0,0]+2*s[:,None]*P[1,0])@f.Tp
            W0=W[0,0]@f.Tw;Ws=W[1,0]@f.Tw;Wz=W[0,1]@f.Tw
            WL=(8*W[1,0]+4*s[:,None]*W[2,0]+W[0,2])@f.Tw
            A=evalmat(PA,ca);C=evalmat(PC,ca);B=evalmat(W0,cb)
            As=evalmat(-P[1,1]@f.Tp,ca);Az=evalmat(-P[0,2]@f.Tp,ca)
            Cs=evalmat((4*P[1,0]+2*s[:,None]*P[2,0])@f.Tp,ca)
            Cz=evalmat((2*P[0,1]+2*s[:,None]*P[1,1])@f.Tp,ca)
            AL=evalmat((-8*P[1,1]-4*s[:,None]*P[2,1]-P[0,3])@f.Tp,ca)
            CL=evalmat((16*P[1,0]+32*s[:,None]*P[2,0]+8*s[:,None]**2*P[3,0]+2*P[0,2]+2*s[:,None]*P[1,2])@f.Tp,ca)
            Pr=2*r[:,None]*(W[1,0]@f.Tq);Pz=W[0,1]@f.Tq
            D['M'][sl]=tensor(W0,S)
            K=2*A[:,None]*W0+(2*s*A)[:,None]*Ws+C[:,None]*Wz-NU*WL
            D['L'][sl]=r[:,None]*(tensor(W0,St)+tensor(K,S))
            D['Pr'][sl]=tensor(Pr,T);D['Pz'][sl]=tensor(Pz,T);D['B'][sl]=B
            Rr=r*(evalmat(PA,cat)+A*A+2*s*A*As+C*Az-B*B-NU*AL)+evalmat(Pr,cp)
            Rt=r*(evalmat(W0,cbt)+2*A*B+2*s*A*evalmat(Ws,cb)+C*evalmat(Wz,cb)-NU*evalmat(WL,cb))
            Rz=evalmat(PC,cat)+2*s*A*Cs+C*Cz-NU*CL+evalmat(Pz,cp)
            D['R'][sl]=np.c_[Rr,Rt,Rz]-force(np.c_[r,np.zeros(len(s)),z],t,*self.fc)
            D['u'][sl]=np.c_[r*A,r*B,C]
        return D
    def residual(self,c,D):
        b=D['M']@c[:self.nb];dp=c[self.nb:];R=D['R'].copy()
        R[:,0]+=-D['r']*(2*D['B']*b+b*b)+D['Pr']@dp
        R[:,1]+=D['L']@c[:self.nb];R[:,2]+=D['Pz']@dp
        return R,b
    def pullback(self,c,D,W):
        R,b=self.residual(c,D)
        gb=D['M'].T@(-2*D['r']*(D['B']+b)*W[:,0])+D['L'].T@W[:,1]
        gp=D['Pr'].T@W[:,0]+D['Pz'].T@W[:,2]
        return np.r_[gb,gp]

def make_pool(seed):
    rng=np.random.default_rng(seed)
    r=np.unique(np.r_[np.linspace(0,1.99,21),np.linspace(.16,.8,7)])
    z=np.unique(np.r_[np.linspace(-1.99,1.99,35),np.linspace(-1.98,-1.72,7),np.linspace(1.72,1.98,7)])
    R,Z,T=np.meshgrid(r,z,np.linspace(.25,.75,11),indexing='ij')
    grid=np.c_[R.ravel()**2,Z.ravel(),T.ravel()]
    rnd=np.c_[rng.uniform(0,4,2048),rng.uniform(-2,2,2048),rng.uniform(.25,.75,2048)]
    # Include densely sampled near-axis/central theta dynamics in TRAINING only.
    R,Z,T=np.meshgrid(np.linspace(.025,.55,11),np.linspace(-.24,.24,9),np.linspace(.25,.75,15),indexing='ij')
    return np.r_[grid,rnd,np.c_[R.ravel()**2,Z.ravel(),T.ravel()]]

class Objective:
    def __init__(self,m,order=(32,48),nt=13,pool_seed=9205650,peak_weight=.12,time_penalty=.0003,pressure_penalty=.01,moment_penalty=0.):
        self.peak_weight=peak_weight;self.time_penalty=time_penalty;self.pressure_penalty=pressure_penalty;self.moment_penalty=moment_penalty
        self.m=m;x,wx=leggauss(order[0]);z,wz=leggauss(order[1]);g,gw=leggauss(nt)
        self.times=np.r_[.25,.5+.25*g,.75];self.tw=np.r_[.05,.45*gw,.05];self.nt=len(self.times)
        S,Z,T=np.meshgrid(2*(x+1),2*z,self.times,indexing='ij')
        q=np.c_[S.ravel(),Z.ravel(),T.ravel()];self.nq=len(q)
        self.sw=(np.outer(wx,wz)/4).ravel();self.qw=(self.sw[:,None]*self.tw).ravel()
        pool=make_pool(pool_seed);self.npool=len(pool);self.D=m.cache(np.r_[q,pool]);self.last=None
        self.base_time=np.einsum('s,st->t',self.sw,np.sum(self.D['R'][:self.nq]**2,axis=1).reshape(-1,self.nt))
        self.base_peak=float(np.linalg.norm(self.D['R'][self.nq:],axis=1).max())
        R,Z,T=np.meshgrid(np.linspace(.035,.215,7),np.r_[-np.linspace(.025,.215,7)[::-1],np.linspace(.025,.215,7)],np.linspace(.25,.75,19),indexing='ij')
        self.C=m.cache(np.c_[(R*R*(1-T)).ravel(),(Z*(1-T)**.495).ravel(),T.ravel()])
        self.cscale=np.maximum(abs(self.C['B']),.01)
        self.pr0=self.C['Pr']@m.p;self.pz0=np.sign(self.C['points'][:,1])*(self.C['Pz']@m.p)
        self.prscale=np.maximum(abs(self.pr0),.005);self.pzscale=np.maximum(abs(self.pz0),.005)
        self.times_g=np.linspace(.25,.75,25);tau=1-self.times_g
        self.G=m.cache(np.c_[.01*tau,.1*tau**.495,self.times_g])
        self.gscale=np.c_[tau**.5,tau**.505,tau**.505];self.gvel=self.G['u']*self.gscale
        self.gnorm=np.linalg.norm(self.gvel[0]);self.calls=0
        self.reg=np.r_[np.full(m.nb,.0002),np.full(m.np,1e-10)]
    def stats(self,c):
        R,b=self.m.residual(c,self.D);sq=np.sum(R*R,axis=1)
        mse=self.sw@sq[:self.nq].reshape(-1,self.nt)
        ratio=(self.C['M']@c[:self.m.nb])/self.cscale
        g=self.gvel.copy();g[:,1]+=self.G['r']*self.gscale[:,1]*(self.G['M']@c[:self.m.nb])
        drift=np.linalg.norm(g-g[0],axis=1)/self.gnorm
        dp=c[self.m.nb:];pr=self.pr0+self.C['Pr']@dp;pz=self.pz0+np.sign(self.C['points'][:,1])*(self.C['Pz']@dp)
        try:self.m.candidate(c);bounds=True
        except ValueError:bounds=False
        return dict(training_mse=float(mse@self.tw),training_time_mse=mse.tolist(),worst_time_mse_ratio=float(np.max(mse/self.base_time)),
                    pool_max=float(np.sqrt(sq[self.nq:]).max()),core_swirl_relative_change=float(abs(ratio).max()),
                    core_swirl_min=float((self.C['B']+self.C['M']@c[:self.m.nb]).min()),core_profile_drift_max=float(drift.max()),
                    core_radial_gradient_min=float(pr.min()),core_signed_axial_gradient_min=float(pz.min()),bounds_ok=bounds)
    def value(self,c):
        if self.last is not None and np.array_equal(c,self.last[0]):return self.last[1:]
        m=self.m;D=self.D;R,b=m.residual(c,D);sq=np.sum(R*R,axis=1)
        tau=2e-5;ls=logsumexp(sq[self.nq:]/tau);soft=np.exp(sq[self.nq:]/tau-ls)
        weights=np.r_[self.qw,self.peak_weight*soft]
        cost=.5*self.qw@sq[:self.nq]+.5*self.peak_weight*tau*(ls-np.log(self.npool))
        # No-worsening per-time volume-L2 cap: squared differentiable penalty.
        mt=self.sw@sq[:self.nq].reshape(-1,self.nt);ex=np.maximum(mt/self.base_time-.998,0)
        cost+=.5*self.time_penalty*np.mean(ex*ex)
        weights[:self.nq]+=(self.sw[:,None]*(2*self.time_penalty/self.nt*ex/self.base_time)).ravel()
        # Broad pool cap prevents average loss from absorbing a new edge peak.
        cap=.995*self.base_peak;expeak=np.maximum(np.sqrt(sq[self.nq:])/cap-1,0)
        cost+=.5*.0005*np.mean(expeak**2)
        weights[self.nq:]+=.0005/self.npool*expeak/np.maximum(cap*np.sqrt(sq[self.nq:]),1e-30)
        grad=m.pullback(c,D,weights[:,None]*R)
        if self.moment_penalty:
            rt=D['r'][:self.nq];Bt=D['B'][:self.nq]+b[:self.nq]
            v=D['u'][:self.nq]
            mom=16*np.pi*(self.sw@(v[:,2]**2-.5*v[:,0]**2-.5*rt**2*Bt**2).reshape(-1,self.nt))
            cost+=.5*self.moment_penalty*np.mean(mom**2)
            gm=-16*np.pi*self.moment_penalty/self.nt*(self.sw[:,None]*mom).ravel()*rt**2*Bt
            grad[:m.nb]+=D['M'][:self.nq].T@gm
        # Keep actual core swirl and the original probe's scaled-profile gate.
        ratio=(self.C['M']@c[:m.nb])/self.cscale;ex=np.maximum(abs(ratio)-.025,0)
        cost+=.5*.03*np.mean(ex**2);grad[:m.nb]+=.03/len(ex)*(self.C['M'].T@(ex*np.sign(ratio)/self.cscale))
        gv=self.gvel.copy();gv[:,1]+=self.G['r']*self.gscale[:,1]*(self.G['M']@c[:m.nb]);diff=gv-gv[0]
        dn=np.sqrt(np.sum(diff*diff,axis=1)+1e-30);ex=np.maximum(dn/self.gnorm-.049,0)
        cost+=.5*.1*np.mean(ex*ex);grad[:m.nb]+=.1/len(ex)*(self.G['M'].T@(ex*diff[:,1]/(dn*self.gnorm)*self.G['r']*self.gscale[:,1]))
        # Pressure direction constraints are handled jointly, never inferred from a plot.
        dp=c[m.nb:];cz=np.sign(self.C['points'][:,1]);pr=self.pr0+self.C['Pr']@dp;pz=self.pz0+cz*(self.C['Pz']@dp)
        er=np.minimum((pr-1e-6)/self.prscale,0);ez=np.minimum((pz-1e-6)/self.pzscale,0)
        cost+=.5*self.pressure_penalty*(np.mean(er*er)+np.mean(ez*ez))
        grad[m.nb:]+=self.pressure_penalty/len(er)*(self.C['Pr'].T@(er/self.prscale)+self.C['Pz'].T@(ez*cz/self.pzscale))
        cost+=.5*np.dot(self.reg*c,c);grad+=self.reg*c
        # Raw finite-family coefficient bounds, also checked strictly at selection.
        db,dp=m.delta(c);rb=m.raw[m.f.n:2*m.f.n]+db/m.amp;rp=m.p+dp
        eb=np.maximum(abs(rb)-3.999,0);ep=np.maximum(abs(rp)-99.99,0)
        cost+=.5*(eb@eb+ep@ep);xb=(eb*np.sign(rb)/m.amp).reshape(m.f.ns,m.f.nt)
        grad[:m.nb]+=(xb[:,1:]-xb[:,0,None]*m.f.q0[1:]).ravel();grad[m.nb:]+=ep*np.sign(rp)
        self.calls+=1;self.last=(c.copy(),float(cost),grad);return self.last[1:]
    def metric(self):
        m=self.m;D=self.D;n=m.dim
        H=np.zeros((n,n));ids=np.arange(0,self.nq,6)
        for start in range(0,len(ids),256):
            k=ids[start:start+256];w=self.qw[k]*6;Jr=-2*(D['r'][k]*D['B'][k])[:,None]*D['M'][k]
            L=D['L'][k];Pr=D['Pr'][k];Pz=D['Pz'][k]
            H[:m.nb,:m.nb]+=Jr.T@(w[:,None]*Jr)+L.T@(w[:,None]*L)
            H[:m.nb,m.nb:]+=Jr.T@(w[:,None]*Pr)
            H[m.nb:,m.nb:]+=Pr.T@(w[:,None]*Pr)+Pz.T@(w[:,None]*Pz)
        H[m.nb:,:m.nb]=H[:m.nb,m.nb:].T
        H.flat[::n+1]+=self.reg
        sc=1/np.sqrt(np.maximum(np.diag(H),1e-12));Hs=H*sc[:,None]*sc[None,:]
        Hs.flat[::n+1]+=1e-4;L=cholesky(Hs,lower=True)
        return L,sc

def run(parent,out,ident='ST056-J',iterations=240,seconds=360,lock_probe=False,peak_weight=.12,time_penalty=.0003,pressure_penalty=.01,moment_penalty=0.):
    out=Path(out)
    if out.exists() and any(out.iterdir()):raise ValueError('Nonempty output directory')
    out.mkdir(parents=True,exist_ok=True)
    reg=dict(id=ident,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),method='joint initial-flat swirl plus all compact pressure; full centrifugal nonlinearity',
             original_physical_gates='UNCHANGED',frozen='poloidal velocity and restricted original force',swirl_parameters=756,pressure_parameters=864,
             training_space=[32,48],training_time='13Gauss plus endpoints',training_pool_seed=9205650,heldout_seeds=[9205691,9205692],
             core_swirl_band=.03001,core_scaled_probe_cap=.0495,acceptance_peak_training_ratio=1.0,acceptance_all_time_mse_ratio=1.0,
             maxiter=iterations,wall_seconds=seconds,lock_original_probe=lock_probe,peak_weight=peak_weight,time_penalty=time_penalty,pressure_penalty=pressure_penalty,moment_penalty=moment_penalty,registered_utc=datetime.now(timezone.utc).isoformat(),pde_validated=False)
    save(out/'registration.json',reg);started=time.monotonic();m=JointModel(parent);o=Objective(m,peak_weight=peak_weight,time_penalty=time_penalty,pressure_penalty=pressure_penalty,moment_penalty=moment_penalty)
    zero=np.zeros(m.dim);base=o.stats(zero);basecost=o.value(zero)[0];save(out/'parent_training.json',base)
    print('BASE',json.dumps(base),flush=True)
    L,sc=o.metric();np.savez_compressed(out/'preconditioner.npz',L=L,sc=sc)
    U=np.zeros((0,m.dim))
    if lock_probe:
        C=np.c_[o.G['M'],np.zeros((len(o.times_g),m.np))]
        E=solve_triangular(L,(C*sc).T,lower=True,check_finite=False).T
        _,sing,Vt=svd(E,full_matrices=False);rank=int(np.sum(sing>sing[0]*1e-10));U=Vt[:rank]
        save(out/'probe_lock.json',dict(rank=rank,singular_values=sing.tolist(),row_count=len(C),scope='Frozen original core probe at 25 training times; not the complete core.'))
    project=lambda y:y-U.T@(U@y)
    transform=lambda y:sc*solve_triangular(L.T,project(y),lower=False,check_finite=False)
    rng=np.random.default_rng(9205651);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=1e-7
    err=abs((o.value(h*d)[0]-o.value(-h*d)[0])/(2*h)-o.value(zero)[1]@d)
    save(out/'calibration.json',dict(objective_gradient_error=err,step=h))
    if err>1e-7:raise ValueError('Gradient check failed')
    history=[];best=[basecost,zero.copy()];last=[zero.copy()]
    def feasible(st):
        return st['bounds_ok'] and st['pool_max']<=base['pool_max']*(1+1e-8) and st['worst_time_mse_ratio']<=1+1e-8 and st['core_swirl_relative_change']<=.03001 and st['core_swirl_min']>0 and st['core_profile_drift_max']<=.0495 and st['core_radial_gradient_min']>=0 and st['core_signed_axial_gradient_min']>=0
    def fun(y):
        c=transform(y);cost,grad=o.value(c);return cost,project(solve_triangular(L,sc*grad,lower=True,check_finite=False))
    class BudgetStop(Exception):pass
    def callback(y):
        c=transform(y);last[0]=c.copy();st=o.stats(c);cost=o.value(c)[0];ok=feasible(st)
        row=dict(iteration=len(history)+1,objective=cost,feasible=ok,**st,elapsed=time.monotonic()-started);history.append(row)
        if ok and cost<best[0]:
            best[:]=[cost,c.copy()];m.f.save(m.candidate(c),out/'best_feasible.json',row);np.save(out/'best_delta.npy',c)
        if len(history)%5==0:
            np.save(out/'checkpoint.npy',c);save(out/'history.json',history);print('ITER',json.dumps({k:v for k,v in row.items() if k!='training_time_mse'}),flush=True)
        if time.monotonic()-started>seconds:raise BudgetStop
    try:
        ret=minimize(fun,zero,method='L-BFGS-B',jac=True,callback=callback,options=dict(maxiter=iterations,ftol=1e-14,gtol=5e-10,maxcor=25,maxls=35))
        end=transform(ret.x);ok=bool(ret.success);message=str(ret.message);nfev=int(ret.nfev)
    except BudgetStop:end=last[0];ok=False;message='Registered wall budget; numerical convergence not claimed';nfev=o.calls
    line=[]
    for alpha in (1.,.8,.6,.4,.2,.1,.05,.02):
        c=alpha*end;st=o.stats(c);v=o.value(c)[0];line.append(dict(alpha=alpha,objective=v,feasible=feasible(st),**st))
        if feasible(st) and v<best[0]:best[:]=[v,c.copy()]
    c=best[1];raw=m.candidate(c);np.save(out/'joint_delta.npy',c)
    m.f.save(raw,out/'candidate.json',dict(id=ident,source='ST056 joint swirl-pressure correction',parent_sha256=reg['parent_sha256'],poloidal_and_force_frozen=True,scientific_acceptance=False))
    summary=dict(**reg,iterations=len(history),nfev=nfev,optimizer_success=ok,message=message,parent_training=base,selected_training=o.stats(c),
                 nonzero_change=bool(np.any(c)),elapsed=time.monotonic()-started,coefficient_sha256=hashlib.sha256(raw.astype('<f8').tobytes()).hexdigest(),
                 candidate_sha256=hashlib.sha256((out/'candidate.json').read_bytes()).hexdigest())
    save(out/'summary.json',summary);save(out/'history.json',history);save(out/'line_recovery.json',line);print('SUMMARY',json.dumps(summary),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',default='ST056-J');p.add_argument('--seconds',type=int,default=360);p.add_argument('--iterations',type=int,default=240);p.add_argument('--lock-probe',action='store_true');p.add_argument('--peak-weight',type=float,default=.12);p.add_argument('--time-penalty',type=float,default=.0003);p.add_argument('--pressure-penalty',type=float,default=.01);p.add_argument('--moment-penalty',type=float,default=0.);a=p.parse_args();run(a.parent,a.out,a.id,a.iterations,a.seconds,a.lock_probe,a.peak_weight,a.time_penalty,a.pressure_penalty,a.moment_penalty)
