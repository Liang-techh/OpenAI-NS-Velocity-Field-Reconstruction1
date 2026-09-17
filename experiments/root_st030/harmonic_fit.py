"""ST004 full-momentum continuation with higher compact-pressure moments.
Uses exact-form polynomial weak moments, never replaces momentum acceptance.
"""
from __future__ import annotations
import argparse, json, time, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular, cho_factor, cho_solve
from spacetime import Family, quad, bumped_basis
from gauss_newton import Objective as BaseObjective


def harmonic_tensors(f,order=80,degrees=(2,4,6,8)):
    import sympy as sp
    s,z=sp.symbols('s z',real=True)
    sq,zq,w=quad(order);r=np.sqrt(sq)
    P={k:v@f.Tp for k,v in f.basis(sq,zq,f.nr,f.nz,True).items()}
    B=f.basis(sq,zq,f.nr,f.nz,False)[0,0]@f.Tw
    UR=-r[:,None]*P[0,1];UZ=2*P[0,0]+2*sq[:,None]*P[1,0];UT=r[:,None]*B
    Mp=[];Mb=[];gr=[];gz=[]
    for degree in degrees:
        H=sp.expand(sum((-1)**k*sp.factorial(2*degree-2*k)*z**(degree-2*k)*(s+z*z)**k/(2**degree*sp.factorial(k)*sp.factorial(degree-k)*sp.factorial(degree-2*k)) for k in range(degree//2+1)))
        hs,hz=sp.diff(H,s),sp.diff(H,z)
        hss,hsz,hzz=sp.diff(H,s,2),sp.diff(H,s,z),sp.diff(H,z,2)
        assert sp.expand(4*hs+4*s*hss+hzz)==0
        hs,hz,hss,hsz,hzz=[np.broadcast_to(sp.lambdify((s,z),v,'numpy')(sq,zq),sq.shape) for v in (hs,hz,hss,hsz,hzz)]
        mp=-(UR.T@((w*(2*hs+4*sq*hss))[:,None]*UR)+UZ.T@((w*hzz)[:,None]*UZ)+UR.T@((w*2*r*hsz)[:,None]*UZ)+UZ.T@((w*2*r*hsz)[:,None]*UR))
        mb=-(UT.T@((w*2*hs)[:,None]*UT))
        Mp.append(mp);Mb.append(mb);gr.append(2*r*hs);gz.append(hz)
    gr=np.array(gr).T;gz=np.array(gz).T
    G=gr.T@(w[:,None]*gr)+gz.T@(w[:,None]*gz)
    scale=np.sqrt(np.diag(G));cor=G/scale[:,None]/scale[None,:]
    transform=solve_triangular(cholesky(cor,lower=True),np.diag(1/scale),lower=True)
    return np.einsum('ij,jab->iab',transform,Mp),np.einsum('ij,jab->iab',transform,Mb),G


class Objective(BaseObjective):
    def __init__(self,family,seed=9172627,count=6144,moment_weight=10.):
        self.moment_weight=moment_weight
        super().__init__(family,seed,count)
        mp,mb,self.harmonic_gram=harmonic_tensors(family)
        self.Hp=self.torch.tensor(mp);self.Hb=self.torch.tensor(mb)

    def extra(self,v):
        torch=self.torch;f=self.f;nt=f.nt;n=f.n
        old=super().extra(v)
        # Replace the old degree-2-only moment block; retain energy/core/torque.
        m=len(self.CT);al=v[:n].reshape(f.ns,nt);be=v[n:2*n].reshape(f.ns,nt)
        scale=torch.sqrt(2/(torch.sum((al@self.q0)**2)+torch.sum((be@self.q0)**2)))
        ap=scale*al@self.CT.T;bp=scale*be@self.CT.T
        moments=torch.einsum('at,kab,bt->kt',ap,self.Hp,ap)+torch.einsum('at,kab,bt->kt',bp,self.Hb,bp)
        return torch.cat((old[:-2*m],self.moment_weight*moments.flatten()/np.sqrt(m),old[-m:]))


def normal_lm(warm,out,maxiter=150,count=6144,moment_weight=10.,seed=9172627):
    """Damped Gauss-Newton, explicitly bounded steps and objective line search.
    Normal equations are double precision; damping is never called a proof.
    """
    f,x=Family.load(warm);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    reg=dict(experiment={9172627:'CR-ROOT-ST004',9172631:'CR-ROOT-ST006'}.get(seed,'additional_continuation'),start_sha256=hashlib.sha256(Path(warm).read_bytes()).hexdigest(),seed=seed,training_points=count,moment_degrees=[2,4,6,8],moment_weight=moment_weight,max_iterations=maxiter,method='column-scaled damped normal Gauss-Newton with line search',physical_gates='unchanged from ST003',validation_seed=seed+1,created_before_fit=True)
    (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n');f.save(x,out/'initial.json')
    obj=Objective(f,seed,count,moment_weight);start=time.time();damping=1e-5;log=[]
    lo=np.r_[np.full(2*f.n,-4.),np.full(f.n,-100.),0.,0.];hi=-lo;hi[-2:]=10.
    cost=np.inf;best=x.copy();best_cost=np.inf
    for it in range(maxiter):
        residual=obj.fun(x);J=obj.jac(x);cost=float(residual@residual)
        if cost<best_cost: best_cost=cost;best=x.copy();f.save(best,out/'candidate.json',dict(training_loss=cost,iteration=it,pde_validated=False))
        grad=J.T@residual
        scale=1/np.maximum(np.linalg.norm(J,axis=0),1e-10)
        K=J*scale[None,:];normal=K.T@K;g=K.T@residual
        accepted=False
        for attempt in range(8):
            mat=normal.copy();mat.flat[::len(x)+1]+=damping
            try: dx=-scale*cho_solve(cho_factor(mat,lower=True,check_finite=False),g,check_finite=False)
            except np.linalg.LinAlgError:damping*=10;continue
            # Trust step avoids discontinuous clipping and preserves all bounds.
            limits=np.full_like(dx,np.inf)
            pos=dx>0;neg=dx<0
            limits[pos]=(hi[pos]-x[pos])/dx[pos];limits[neg]=(lo[neg]-x[neg])/dx[neg]
            alpha=min(1.,.995*max(0.,float(np.min(limits))))
            # Active-bound directions require projected clipping (objective still checked).
            if alpha<1e-8: alpha=1.
            for line in range(10):
                trial=np.clip(x+alpha*dx,lo,hi)
                tr=obj.fun(trial);new=float(tr@tr)
                if np.isfinite(new) and new<cost:
                    pred=cost-float(np.linalg.norm(residual+J@(trial-x))**2)
                    rho=(cost-new)/max(pred,1e-30)
                    x=trial;accepted=True
                    if rho>.75:damping=max(damping*.35,1e-12)
                    elif rho<.25:damping=min(damping*4,1e4)
                    break
                alpha*=.5
            if accepted:break
            damping*=10
        row=dict(iteration=it,loss=cost,new_loss=new if accepted else cost,damping=damping,accepted=accepted,gradient_inf=float(np.max(np.abs(grad))),elapsed=time.time()-start,objective_evaluations=len(obj.history))
        log.append(row)
        if it%5==0 or not accepted:print(json.dumps(row),flush=True)
        (out/'optimizer_history.json').write_text(json.dumps(log,indent=2)+'\n')
        if not accepted:break
        if new<best_cost:best_cost=new;best=x.copy();f.save(best,out/'candidate.json',dict(training_loss=new,iteration=it+1,pde_validated=False))
    info=dict(**reg,iterations=len(log),function_evaluations=len(obj.history),best_loss=best_cost,elapsed=time.time()-start,pde_validated=False,stop='iteration budget or unsuccessful damped step, no global optimality claim')
    f.save(best,out/'candidate.json',info);(out/'training_summary.json').write_text(json.dumps(info,indent=2)+'\n');(out/'training_history.json').write_text(json.dumps(obj.history,indent=2)+'\n');print(json.dumps(info,indent=2),flush=True)
    return f,best

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--warm',default='gauss_7x7x8/candidate.json');p.add_argument('--out',default='harmonic_7x7x8');p.add_argument('--maxiter',type=int,default=150);p.add_argument('--count',type=int,default=6144);p.add_argument('--weight',type=float,default=10.);p.add_argument('--seed',type=int,default=9172627);a=p.parse_args();normal_lm(a.warm,a.out,a.maxiter,a.count,a.weight,a.seed)
