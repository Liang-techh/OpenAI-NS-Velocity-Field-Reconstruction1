"""CR-ROOT-ST030: exact curvature on original mixed collocation, compact p eliminated.
No change to the physical family, force bounds, nonzero energy or acceptance gates.
Pressure columns have a fixed span: orthogonal variable projection is exact.
"""
from __future__ import annotations
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.linalg import qr, solve_triangular, eigh, cho_factor, cho_solve
from spacetime import Family, NU
from harmonic_fit import Objective

class ProjectedExact:
    def __init__(self,f,seed=9172901,count=6144,weight=10.):
        self.f=f;self.o=Objective(f,seed,count,weight);self.n=f.n;self.N=count;self.nv=2*f.n
        D=self.o.D;n=f.n;N=count;r=np.sqrt(D['s'])
        P=np.zeros((3*N,n));P[:N]=2*r[:,None]*D['Qs'];P[2*N:]=D['Qz'];P/=np.sqrt(N)
        sc=1/np.maximum(np.linalg.norm(P,axis=0),1e-30)
        Q,R=qr(P*sc[None,:],mode='economic',check_finite=False)
        if np.min(np.abs(np.diag(R)))<1e-10:raise ValueError('Rank deficient pressure span')
        self.Q=Q;self.R=R;self.psc=sc;self.info={'pressure_scaled_R_condition':float(np.linalg.cond(R)),'pressure_rank':n,'samples':count,'training_seed':seed}
        self.last_y=None;self.calls=0
    def initial(self,raw):return np.r_[raw[:self.nv],raw[-2:]]
    def evaluate(self,y):
        if self.last_y is not None and np.array_equal(y,self.last_y):return
        n=self.n;N=self.N;o=self.o
        raw=np.r_[y[:2*n],np.zeros(n),y[-2:]]
        o.evaluate(raw)
        rr=o.cache_r.copy();J0=o.cache_j;J=np.column_stack((J0[:,:2*n],J0[:,-2:]))
        qp=self.Q.T@rr[:3*N];press=-self.psc*solve_triangular(self.R,qp,lower=False,check_finite=False)
        rr[:3*N]-=self.Q@qp
        J[:3*N]-=self.Q@(self.Q.T@J[:3*N])
        raw[2*n:3*n]=press
        self.raw=raw;self.res=rr;self.J=J;self.last_y=y.copy();self.calls+=1
    def fun(self,y):self.evaluate(y);return self.res
    def linearize(self,y):self.evaluate(y);return self.res,self.J
    def curvature(self,y):
        """H(0.5 ||r||^2)-J^T J, including normalization and all extra penalties."""
        self.evaluate(y);f=self.f;o=self.o;n=self.n;N=self.N;D=o.D
        aa,bb,qq,fc,amp=f.coefficients(self.raw)
        V={k:D[k]@(bb if k.startswith('B') else aa) for k in ('A','As','Az','At','AL','B','Bs','Bz','Bt','BL','C','Cs','Cz','Ct','CL')}
        A,B,C=V['A'],V['B'],V['C'];s=D['s'];r=np.sqrt(s)
        # res already contains 1/sqrt(N), while D matrices have no norm weights.
        wR=self.res[:3*N].reshape(3,N)/np.sqrt(N);wr=r*wR[0];wb=r*wR[1];wz=wR[2]
        def sym(a,b,w):
            K=D[a].T@(w[:,None]*D[b]);return K+K.T
        H=np.zeros((2*n,2*n))
        H[:n,:n]=(sym('A','A',wr)+sym('A','As',2*s*wr)+sym('C','Az',wr)+sym('A','Cs',2*s*wz)+sym('C','Cz',wz))
        H[n:,n:]=-sym('B','B',wr)
        H[:n,n:]=(D['A'].T@((2*wb)[:,None]*D['B'])+D['A'].T@((2*s*wb)[:,None]*D['Bs'])+D['C'].T@(wb[:,None]*D['Bz']))
        H[n:,:n]=H[:n,n:].T
        ga=(D['At'].T@wr+D['A'].T@((2*A+2*s*V['As'])*wr+(2*B+2*s*V['Bs'])*wb+2*s*V['Cs']*wz)
            +D['As'].T@(2*s*A*wr)+D['C'].T@(V['Az']*wr+V['Bz']*wb+V['Cz']*wz)
            +D['Az'].T@(C*wr)-NU*D['AL'].T@wr+D['Ct'].T@wz+D['Cs'].T@(2*s*A*wz)+D['Cz'].T@(C*wz)-NU*D['CL'].T@wz)
        gb=(-D['B'].T@(2*B*wr)+D['Bt'].T@wb+D['B'].T@(2*A*wb)+D['Bs'].T@(2*s*A*wb)+D['Bz'].T@(C*wb)-NU*D['BL'].T@wb)
        gp=np.r_[ga,gb];q=y[:2*n];qmat=q.reshape(2*f.ns,f.nt);ini=qmat@f.q0;den=ini@ini
        v=(ini[:,None]*f.q0).ravel()/den
        Hq=H@q;H=amp**2*(H-np.outer(Hq,v)-np.outer(v,Hq)+(q@Hq)*np.outer(v,v))
        qg=q@gp
        H+=amp*(3*qg*np.outer(v,v)-np.outer(v,gp)-np.outer(gp,v))
        H-=amp*qg/den*np.kron(np.eye(2*f.ns),np.outer(f.q0,f.q0))
        out=np.zeros((2*n+2,2*n+2));out[:2*n,:2*n]=H
        # Extras are small, potentially active piecewise-smooth constraints.
        # Compute their scalar Hessian independently with AD (pressure absent).
        torch=o.torch;tx=torch.tensor(y,requires_grad=True)
        def ext(v):return o.extra(torch.cat((v[:2*n],torch.zeros(n),v[-2:])))
        def loss(v):
            e=ext(v);return .5*torch.sum(e*e)
        HH=torch.autograd.functional.hessian(loss,tx,vectorize=True).detach().numpy()
        ej=self.J[3*N:]
        out+=HH-ej.T@ej
        if not np.all(np.isfinite(out)):raise ValueError('Nonfinite exact Hessian correction')
        return (out+out.T)/2
    def hessian(self,y):
        r,J=self.linearize(y);N=J.T@J;return N+self.curvature(y),N

def trust_step(vals,vec,g,delta):
    c=vec.T@g;low=max(0.,-vals[0]);eps=1e-12
    if vals[0]>eps:
        z=-c/vals
        if np.linalg.norm(z)<=delta:return vec@z,0.
    lo=low+eps
    z=-c/(vals+lo)
    if np.linalg.norm(z)<delta and vals[0]<0:
        z[0]=(-1 if c[0]>=0 else 1)*np.sqrt(max(0.,delta**2-np.sum(z[1:]**2)))
        return vec@z,lo
    hi=max(1.,2*lo)
    while np.linalg.norm(c/(vals+hi))>delta:hi*=2
    for _ in range(65):
        mid=(lo+hi)/2
        if np.linalg.norm(c/(vals+mid))>delta:lo=mid
        else:hi=mid
    return vec@(-c/(vals+hi)),hi

def run(warm,out,iterations=60,seed=9172901,count=6144,weight=10.,gauss_only=False):
    f,x=Family.load(warm);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    reg=dict(experiment='CR-ROOT-ST030',initial_sha256=hashlib.sha256(Path(warm).read_bytes()).hexdigest(),seed=seed,validation_seed=9172910,count=count,iterations=iterations,moment_weight=weight,method='exact mixed-collocation pressure-projected Hessian trust region' if not gauss_only else 'projected Gauss Newton control',physical_gates='UNCHANGED',created_before_fit=True)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n');f.save(x,out/'initial.json')
    t0=time.time();o=ProjectedExact(f,seed,count,weight);y=o.initial(x);delta=.05;history=[];best=np.inf;best_raw=x.copy()
    for it in range(iterations):
        r,J=o.linearize(y);cost=float(r@r);grad=J.T@r;N=J.T@J;H=N if gauss_only else N+o.curvature(y)
        sc=1/np.maximum(np.sqrt(np.maximum(np.diag(N),0)),1e-7);K=H*sc[:,None]*sc[None,:];g=sc*grad
        # Remove the exact amplitude redundancy with an orthogonal Householder.
        gauge=y.copy();gauge[-2:]=0;gauge/=sc;gauge/=np.linalg.norm(gauge)
        v=gauge.copy();v[0]+=np.copysign(1.,v[0]);v/=np.linalg.norm(v)
        Kv=K@v;K=K-2*np.outer(v,Kv)-2*np.outer(Kv,v)+4*(v@Kv)*np.outer(v,v);gt=g-2*v*(v@g)
        vals,vec=eigh(K[1:,1:],check_finite=False);accepted=False;rho=0.;new=cost
        if cost<best:best=cost;best_raw=o.raw.copy();f.save(best_raw,out/'candidate.json',dict(loss=cost,iteration=it))
        for attempt in range(12):
            d,shift=trust_step(vals,vec,gt[1:],delta);d=np.r_[0.,d];d-=2*v*(v@d);d*=sc
            trial=y+d;trial[-2:]=np.clip(trial[-2:],0,10)
            if np.max(np.abs(trial[:2*f.n]))>4:delta*=.25;continue
            try:
                rr=o.fun(trial);new=float(rr@rr);amp=f.coefficients(o.raw)[4]
                if np.max(np.abs(o.raw[2*f.n:-2]))>100 or not 1e-4<=amp<=100:delta*=.25;continue
            except ValueError:delta*=.25;continue
            step=trial-y;pred=-grad@step-.5*step@H@step;rho=.5*(cost-new)/max(pred,1e-30)
            if new<cost and rho>.05:
                y=trial;accepted=True
                if rho>.75:delta=min(delta*2,2.)
                elif rho<.25:delta*=.25
                break
            delta*=.25
        if accepted and new<best:best=new;best_raw=o.raw.copy();f.save(best_raw,out/'candidate.json',dict(loss=best,iteration=it+1))
        row=dict(iteration=it,loss=cost,new_loss=new if accepted else cost,pde_mse=float(r[:3*count]@r[:3*count]),scaled_gradient_max=float(np.max(np.abs(gt[1:]))),lowest_curvature=float(vals[0]),trust_radius=delta,rho=float(rho),accepted=accepted,elapsed=time.time()-t0)
        history.append(row);(out/'history.json').write_text(json.dumps(history,indent=2)+'\n');print(json.dumps(row),flush=True)
        if not accepted or delta<1e-10:break
    summary=dict(**reg,geometry=o.info,elapsed=time.time()-t0,iterations_run=len(history),best_training_loss=best,function_evaluations=o.calls,pde_validated=False,stop='bounded iteration budget or failed trust step; no global optimality claim')
    f.save(best_raw,out/'candidate.json',summary);(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--warm',required=True);p.add_argument('--out',required=True);p.add_argument('--iterations',type=int,default=60);p.add_argument('--seed',type=int,default=9172901);p.add_argument('--count',type=int,default=6144);p.add_argument('--weight',type=float,default=10.);p.add_argument('--gauss-only',action='store_true');a=p.parse_args();run(a.warm,a.out,a.iterations,a.seed,a.count,a.weight,a.gauss_only)
