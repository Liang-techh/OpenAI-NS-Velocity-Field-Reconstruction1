"""ST050R: restarted pressure-Poisson and vorticity-morphology experiment.
Uses the unchanged compact field, normalization, viscosity and restricted force.
None of the auxiliary diagnostics substitutes for full independent NS acceptance.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st047'))
from continuation import Objective as PriorObjective, LocalizedModel, precondition

def atomic_json(path, obj):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)

class BoundaryShear(PriorObjective):
    """Reimplementation of ST048's same declared edge/shear training formulas."""
    def __init__(self,m,shear_ratio=.995,edge_weight=.04,**kw):
        super().__init__(m,**kw);self.shear_ratio=shear_ratio;self.edge_weight=edge_weight
        r,t=np.meshgrid([.04,.08,.14,.23,.36,.48,.60],np.linspace(.25,.75,17),indexing='ij')
        self.sd=m.cache(r*r,np.zeros_like(r),t)
        self.sp=self.shear(np.zeros(m.dim))[0];self.ss=np.sign(self.sp);self.sscale=np.maximum(abs(self.sp),2e-4)
        et=[.25,.27,.33,.50,.67,.73,.75]
        rr=np.r_[np.linspace(0,.5,6),np.linspace(.65,1.95,8)]
        zz=np.r_[-np.linspace(1.05,1.96,12)[::-1],np.linspace(1.05,1.96,12)]
        r,z,t=np.meshgrid(rr,zz,et,indexing='ij');p1=np.c_[r.ravel()**2,z.ravel(),t.ravel()]
        r,z,t=np.meshgrid(np.linspace(1.25,1.96,9),np.linspace(-1.95,1.95,17),et,indexing='ij')
        self.edge_points=np.r_[p1,np.c_[r.ravel()**2,z.ravel(),t.ravel()]]
        self.edge=m.cache(*self.edge_points.T);self.ec=None;self.cc=None;self.fc=None
    def shear(self,c):
        m=self.m;D=self.sd;lam,dl=m.norm(c);factor=2*np.sqrt(D['s'])
        bar=D['v']['Cs']+D['M']['Cs']@c[:m.na];v=lam*factor*bar
        J=np.zeros((len(v),m.dim));J[:,:m.na]=lam*factor[:,None]*D['M']['Cs'];J[:,:m.nv]+=(factor*bar)[:,None]*dl
        return v,J
    def edge_momentum(self,c):
        if self.ec is None or not np.array_equal(self.ec[0],c):self.ec=(c.copy(),*self.m.momentum(c,self.edge))
        return self.ec[1:]
    def constraints(self,c):
        if self.cc is not None and np.array_equal(self.cc[0],c):return self.cc[1:]
        v,J=super().constraints(c);sh,K=self.shear(c);low=(self.ss*sh-self.shear_ratio*abs(self.sp))/self.sscale
        high=(1.25*abs(self.sp)+1e-5-self.ss*sh)/self.sscale;K=self.ss[:,None]*K/self.sscale[:,None]
        self.cc=(c.copy(),np.r_[v,low,high],np.vstack((J,K,-K)));return self.cc[1:]
    def fun(self,c):
        if self.fc is not None and np.array_equal(self.fc[0],c):return self.fc[1:]
        loss,g=super().fun(c);r,J=self.edge_momentum(c);sq=np.sum(r*r,axis=1);tau=.0001
        normal=logsumexp(sq/tau);weights=np.exp(sq/tau-normal)
        loss+=self.edge_weight*tau*(normal-np.log(len(sq)))+.01*float(np.mean(sq))
        g=g+2*np.einsum('n,ni,nik->k',self.edge_weight*weights+.01/len(sq),r,J)
        self.fc=(c.copy(),float(loss),g);return self.fc[1:]

def pressure_laplacian_cache(m,D):
    f=m.f;s,z,t=D['s'],D['z'],D['t'];base=[];mat=[]
    for k in range(0,len(s),256):
        sl=slice(k,k+256);Q=f.basis(s[sl],z[sl],f.nr,f.nz,False)
        L=(4*Q[1,0]+4*s[sl,None]*Q[2,0]+Q[0,2])@f.Tq
        T=np.column_stack([Chebyshev.basis(j)(4*t[sl]-2) for j in range(f.nt)])
        L=(L[:,:,None]*T[:,None,:]).reshape(-1,f.n)
        base.append(L@m.p);mat.append(L@m.Mp)
    return np.concatenate(base),np.concatenate(mat)

def poisson_residual(m,c,D,lp):
    """div R = tr(grad(u)^2)+Delta(p), for this divergence-free u/f."""
    a=c[:m.na];b=c[m.na:m.nv];p=c[m.nv:m.nv+m.np];s=D['s'];M=D['M'];v=D['v']
    A=v['A']+M['A']@a;As=v['As']+M['As']@a;Az=v['Az']+M['Az']@a
    Cs=v['Cs']+M['Cs']@a;Cz=v['Cz']+M['Cz']@a;B=v['B']+M['B']@b;Bs=v['Bs']+M['Bs']@b
    X=A+2*s*As;K=X*X+A*A+Cz*Cz-2*B*(B+2*s*Bs)+4*s*Az*Cs
    lam,dl=m.norm(c);J=np.zeros((len(s),m.dim))
    J[:,:m.na]=lam**2*(2*X[:,None]*(M['A']+2*s[:,None]*M['As'])+2*A[:,None]*M['A']+2*Cz[:,None]*M['Cz']+4*s[:,None]*(Cs[:,None]*M['Az']+Az[:,None]*M['Cs']))
    J[:,m.na:m.nv]=lam**2*(-4*(B+s*Bs)[:,None]*M['B']-4*(s*B)[:,None]*M['Bs'])
    J[:,:m.nv]+=2*lam*K[:,None]*dl
    J[:,m.nv:m.nv+m.np]=lp[1]
    return lam**2*K+lp[0]+lp[1]@p,J

def vorticity_grams(m,order=(32,48),times=(.25,.375,.5,.625,.75)):
    x,wx=leggauss(order[0]);z,wz=leggauss(order[1]);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij');s=S.ravel();z=Z.ravel();w=(4*np.pi*np.outer(wx,wz)).ravel();r=np.sqrt(s)
    grams=[]
    for t in times:
        D=m.cache(s,z,t);V=D['v'];M=D['M'];U=np.zeros((len(s),3,m.nv+1))
        U[:,0,0]=-r*V['Bz'];U[:,0,1+m.na:]=-r[:,None]*M['Bz']
        U[:,1,0]=r*(V['Az']-2*V['Cs']);U[:,1,1:1+m.na]=r[:,None]*(M['Az']-2*M['Cs'])
        U[:,2,0]=2*(V['B']+s*V['Bs']);U[:,2,1+m.na:]=2*(M['B']+s[:,None]*M['Bs'])
        grams.append([np.einsum('nik,nil,n->kl',U,U,w*factor,optimize=True) for factor in (np.ones(len(s)),z*z,z**4,s)])
    return np.array(grams)

class PressureMorph(BoundaryShear):
    def __init__(self,m,poisson_weight=.003,morph_ratio=.999,axis_ratio=0.,**kw):
        super().__init__(m,**kw);self.ppweight=poisson_weight;self.morph_ratio=morph_ratio
        x,wx=leggauss(18);z,wz=leggauss(26);times=np.linspace(.25,.75,7)
        S,Z,T=np.meshgrid(2*(x+1),2*z,times,indexing='ij')
        self.pd=m.cache(S,Z,T);self.pw=np.broadcast_to((np.outer(wx,wz)/4)[:,:,None]/len(times),S.shape).ravel()
        self.pl=pressure_laplacian_cache(m,self.pd)
        self.grams=vorticity_grams(m);self.gram_ref=self.grams[:,:,0,0].copy()
        self.axis_ratio=axis_ratio
        z,t=np.meshgrid([-.15,-.075,0.,.075,.15],np.linspace(.25,.75,17),indexing="ij")
        self.rot=m.cache(np.zeros_like(z),z*(1-t)**.495,t) if axis_ratio else None
        self.pc=None;self.pfun=None;self.pcon=None
    def morphology(self,c):
        d=np.r_[1,c[:self.m.nv]];Gd=self.grams@d;v=np.einsum('i,tki->tk',d,Gd);J=2*Gd[:,:,1:]
        ref=self.gram_ref;vals=[];rows=[]
        for k,ratio in [(1,self.morph_ratio),(2,self.morph_ratio)]:
            coef=ratio*ref[:,k]/ref[:,0]
            vals.extend((v[:,k]-coef*v[:,0])/ref[:,k]);rows.extend((J[:,k]-coef[:,None]*J[:,0])/ref[:,k,None])
        coef=1.01*ref[:,3]/ref[:,0]
        vals.extend((coef*v[:,0]-v[:,3])/ref[:,3]);rows.extend((coef[:,None]*J[:,0]-J[:,3])/ref[:,3,None])
        return np.array(vals),np.pad(np.array(rows),((0,0),(0,self.m.dim-self.m.nv)))
    def constraints(self,c):
        if self.pcon is not None and np.array_equal(c,self.pcon[0]):return self.pcon[1:]
        v,J=super().constraints(c);mv,MJ=self.morphology(c)
        if self.axis_ratio:
            D=self.rot;m=self.m
            a=D['v']['A']+D['M']['A']@c[:m.na];b=D['v']['B']+D['M']['B']@c[m.na:m.nv]
            ar=self.axis_ratio;K=np.zeros((len(a),m.dim))
            K[:,:m.na]=ar*D['M']['A'];K[:,m.na:m.nv]=D['M']['B']
            # For inward A and positive B, B>=kappa*|A| ensures enough swirl.
            S=max(float(np.sqrt(np.mean(D['v']['A']**2))),.01)
            JA=np.zeros_like(K);JA[:,:m.na]=-D['M']['A']
            v=np.r_[v,(b+ar*a)/S,(-a-1e-5)/S];J=np.vstack((J,K/S,JA/S))
        self.pcon=(c.copy(),np.r_[v,mv],np.vstack((J,MJ)));return self.pcon[1:]
    def fun(self,c):
        if self.pfun is not None and np.array_equal(c,self.pfun[0]):return self.pfun[1:]
        loss,g=super().fun(c);P,J=poisson_residual(self.m,c,self.pd,self.pl)
        loss+=self.ppweight*float(self.pw@(P*P));g=g+2*self.ppweight*((self.pw*P)@J)
        self.pfun=(c.copy(),float(loss),g);return self.pfun[1:]

def run(parent,out,ident='ST050R-C',radial=5,axial=5,td=4,pp=.003,pressure=.99,profile=.995,acc=.98,anchor=.08,cap=1.,maxiter=260,seconds=480,axis_ratio=0.,pressure_target=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()):raise ValueError('Refuse nonempty output directory')
    reg=dict(id=ident,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),radial=radial,axial=axial,td=td,poisson_weight=pp,pressure_ratio=pressure,pressure_target=pressure_target,axis_ratio=axis_ratio,profile_ratio=profile,acceleration_ratio=acc,anchor_tolerance=anchor,time_cap=cap,training_spatial=[32,48],training_times='13Gauss plus endpoints',morphology_spatial=[32,48],morphology_times=[.25,.375,.5,.625,.75],morphology_min_axial_moment_ratio=.999,max_radial_moment_ratio=1.01,shear_min_ratio=.995,edge_count=3423,maxiter=maxiter,seconds=seconds,validation_seeds=[9175091,9175092],feasibility_tolerance=1e-7,original_physical_gates='UNCHANGED',utc_before_fit=datetime.now(timezone.utc).isoformat(),pde_validated=False)
    atomic_json(out/'registration.json',reg);start=time.monotonic();m=LocalizedModel(parent,radial,axial,td)
    o=PressureMorph(m,poisson_weight=pp,axis_ratio=axis_ratio,pressure_target=pressure_target,profile_ratio=profile,pressure_ratio=pressure,acceleration_ratio=acc,axis_weight=.02,space_order=(32,48),anchor_tolerance=anchor,time_cap=cap,peak_weight=0.,softmax_weight=.02)
    c0=np.zeros(m.dim);T,rs,info=precondition(o,c0,'whiten');atomic_json(out/'conditioning.json',info)
    np.savez_compressed(out/'model.npz',Ma=m.Ma,Mb=m.Mb,Mp=m.Mp);np.save(out/'transform.npy',T)
    lo,hi=np.array(m.bounds).T;bs=1/np.maximum(hi-lo,1e-4);factor=1e4
    def fun(y):
        v,g=o.fun(T@y);return factor*v,factor*T.T@g
    lastcon=[None,None]
    def con(y):
        if lastcon[0] is not None and np.array_equal(y,lastcon[0]):return lastcon[1]
        c=T@y;v,J=o.constraints(c);lastcon[:]=[y.copy(),(np.r_[v*rs,(c-lo)*bs,(hi-c)*bs],np.vstack(((J@T)*rs[:,None],T*bs[:,None],-T*bs[:,None])))];return lastcon[1]
    rng=np.random.default_rng(9175050);d=rng.normal(size=m.dim);d/=np.linalg.norm(d);h=2e-6
    v,g=fun(c0);vv,J=con(c0);check=dict(objective_error=abs((fun(h*d)[0]-fun(-h*d)[0])/(2*h)-g@d),constraint_error=float(np.max(abs((con(h*d)[0]-con(-h*d)[0])/(2*h)-J@d))),parent_pp_rms=float(np.sqrt(o.pw@(poisson_residual(m,c0,o.pd,o.pl)[0]**2))))
    atomic_json(out/'calibration.json',check);print('CALIBRATION',json.dumps(check),flush=True)
    if max(check['objective_error'],check['constraint_error'])>1e-4:raise ValueError('Derivative calibration failed')
    history=[];best=[np.inf,None];last=[c0.copy()];calls=[0]
    class BudgetStop(Exception):pass
    def callback(y):
        c=T@y;last[0]=y.copy();val=o.fun(c)[0];feas=float(o.constraints(c)[0].min());raw=None
        try:raw=m.candidate(c)
        except ValueError:pass
        good=bool(np.all(c>=lo-1e-10) and np.all(c<=hi+1e-10) and raw is not None)
        row=dict(iteration=len(history)+1,loss=val,min_constraint=feas,bounds_ok=good,seconds=time.monotonic()-start);history.append(row)
        if feas>=-1e-7 and good and val<best[0]:
            best[:]=[val,c.copy()];m.f.save(raw,out/'best_feasible.json',row);np.save(out/'best_modifiers.npy',c)
        if len(history)%5==0:
            atomic_json(out/'history.json',history);np.save(out/'checkpoint.npy',c)
            if raw is not None:m.f.save(raw,out/'checkpoint.json',row)
            print(json.dumps(row),flush=True)
        if time.monotonic()-start>seconds:raise BudgetStop
    try:
        ret=minimize(fun,c0,jac=True,method='SLSQP',constraints=[dict(type='ineq',fun=lambda y:con(y)[0],jac=lambda y:con(y)[1])],callback=callback,options=dict(maxiter=maxiter,ftol=2e-9,disp=False));final=T@ret.x;ok=bool(ret.success);message=str(ret.message);nfev=ret.nfev
    except BudgetStop:
        final=T@last[0];ok=False;message='Wall budget; saved checkpoint, no convergence claim';nfev=len(o.history)
    selected=best[1] if best[1] is not None else final
    summary=dict(**reg,iterations=len(history),nfev=nfev,optimizer_success=ok,message=message,selection='best_feasible_training' if best[1] is not None else 'infeasible_last',loss=o.fun(selected)[0],min_constraint=float(o.constraints(selected)[0].min()),parent_pp_rms=check['parent_pp_rms'],candidate_pp_rms=float(np.sqrt(o.pw@(poisson_residual(m,selected,o.pd,o.pl)[0]**2))),elapsed=time.monotonic()-start,dim=m.dim)
    m.f.save(m.candidate(selected),out/'candidate.json',summary);np.save(out/'modifiers.npy',selected);atomic_json(out/'summary.json',summary);atomic_json(out/'history.json',history);print('SUMMARY',json.dumps(summary),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',default='ST050R-C');p.add_argument('--radial',type=int,default=5);p.add_argument('--axial',type=int,default=5);p.add_argument('--pp',type=float,default=.003);p.add_argument('--pressure',type=float,default=.99);p.add_argument('--profile',type=float,default=.995);p.add_argument('--acc',type=float,default=.98);p.add_argument('--anchor',type=float,default=.08);p.add_argument('--cap',type=float,default=1.);p.add_argument('--maxiter',type=int,default=260);p.add_argument('--seconds',type=float,default=480);p.add_argument('--axis-ratio',type=float,default=0.);p.add_argument('--pressure-target',type=float)
    a=p.parse_args();run(a.parent,a.out,ident=a.id,radial=a.radial,axial=a.axial,pp=a.pp,pressure=a.pressure,profile=a.profile,acc=a.acc,anchor=a.anchor,cap=None if a.cap<0 else a.cap,maxiter=a.maxiter,seconds=a.seconds,axis_ratio=a.axis_ratio,pressure_target=a.pressure_target)
