"""ST053: explicit new weak-moment and sampled-peak restrictions.
Reuses ST051 solver mechanics unchanged; never modifies its source or NS validator.
"""
from __future__ import annotations
import argparse,sys,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st051'))
import aligned_continuation as inherited
from pressure_morph import atomic_json
NORM2=88*np.pi/3

def moment_matrices(m,times,order=(48,72)):
    x,wx=leggauss(order[0]);z,wz=leggauss(order[1]);S,Z=np.meshgrid(2*(x+1),2*z,indexing='ij')
    s,z=S.ravel(),Z.ravel();w=(4*np.pi*np.outer(wx,wz)).ravel();r=np.sqrt(s);out=[]
    for t in times:
        D=m.cache(s,z,t);parts=[]
        for key,factor,offset in [('A',r,0),('B',r,m.na),('C',np.ones(len(r)),0)]:
            V=np.zeros((len(s),m.nv+1));V[:,0]=factor*D['v'][key]
            n=m.nb if key=='B' else m.na
            V[:,1+offset:1+offset+n]=factor[:,None]*D['M'][key]
            parts.append(V.T@(w[:,None]*V))
        out.append(parts[2]-.5*(parts[0]+parts[1]))
    return np.array(out)

def weak_moment(m,c,G):
    d=np.r_[1.,c[:m.nv]];Gd=G@d;q=np.einsum('i,ti->t',d,Gd);lam,dl=m.norm(c)
    J=np.zeros((len(q),m.dim));J[:,:m.nv]=2*lam**2*Gd[:,1:]+2*lam*q[:,None]*dl
    return lam**2*q,J

def training_pool():
    r=np.unique(np.r_[np.linspace(0,1.99,23),np.linspace(.1,1.,11)])
    z=np.unique(np.r_[np.linspace(-1.98,1.98,31),-np.linspace(1.84,1.98,17),np.linspace(1.84,1.98,17)])
    R,Z,T=np.meshgrid(r,z,np.linspace(.25,.75,11),indexing='ij');p=np.c_[R.ravel()**2,Z.ravel(),T.ravel()]
    rng=np.random.default_rng(9175350)
    return np.r_[p,np.c_[rng.uniform(0,4,1024),rng.uniform(-2,2,1024),rng.uniform(.25,.75,1024)]]

def pool_residual(m,c,p):
    raw=m.candidate(c);out=[]
    for i in range(0,len(p),256):
        s,z,t=p[i:i+256].T;x=np.c_[np.sqrt(s),np.zeros(len(s)),z]
        out.append(np.linalg.norm(m.f.analytic_residual(raw,x,t),axis=1))
    return np.concatenate(out)

class MomentObjective(inherited.AlignedObjective):
    ratio=.9;weight=.2;output=None
    def __init__(self,m,**kw):
        kw['profile_ratio']=1.;kw['acceleration_ratio']=1.
        super().__init__(m,**kw)
        self.G=moment_matrices(m,self.times);c=np.zeros(m.dim);self.parentD=weak_moment(m,c,self.G)[0]
        p=training_pool();res=pool_residual(m,c,p);self.cap=float(res.max());ids=set(np.argsort(res)[-300:].tolist())
        for t in np.linspace(.25,.75,11):
            i=np.flatnonzero(p[:,2]==t);ids.update(i[np.argsort(res[i])[-32:]].tolist())
        ids.update(np.random.default_rng(9175351).choice(len(p),128,replace=False).tolist())
        self.pool=p;self.active=np.array(sorted(ids));self.peakD=m.cache(*p[self.active].T)
        self.last_f=None;self.last_c=None;self.last_r=None
        if self.output:
            np.savez_compressed(self.output/'training_pool.npz',pool=p,parent_residual=res,active=self.active)
            atomic_json(self.output/'moment_setup.json',dict(parentD=self.parentD.tolist(),cap=self.cap,active_points=len(ids),pool_points=len(p)))
        print('MOMENT_SETUP',self.parentD.tolist(),self.cap,len(ids),flush=True)
    def peaks(self,c):
        if self.last_r is None or not np.array_equal(c,self.last_r[0]):self.last_r=(c.copy(),*self.m.momentum(c,self.peakD))
        return self.last_r[1:]
    def fun(self,c):
        if self.last_f is not None and np.array_equal(c,self.last_f[0]):return self.last_f[1:]
        value,g=super().fun(c);D,J=weak_moment(self.m,c,self.G)
        value+=self.weight*np.mean(D*D)/NORM2;g=g+2*self.weight*(D@J)/(len(D)*NORM2)
        self.last_f=(c.copy(),float(value),g);return self.last_f[1:]
    def constraints(self,c):
        if self.last_c is not None and np.array_equal(c,self.last_c[0]):return self.last_c[1:]
        v,J=super().constraints(c);D,K=weak_moment(self.m,c,self.G)
        scale=np.maximum(abs(self.parentD),.01);limit=self.ratio*abs(self.parentD)
        r,RJ=self.peaks(c);pv=(self.cap**2-np.sum(r*r,axis=1))/self.cap**2
        PJ=-2*np.einsum('ni,nik->nk',r,RJ)/self.cap**2
        self.last_c=(c.copy(),np.r_[v,(limit-D)/scale,(limit+D)/scale,pv],np.vstack((J,-K/scale[:,None],K/scale[:,None],PJ)))
        return self.last_c[1:]

def main():
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--reference',required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--ratio',type=float,default=.9);p.add_argument('--weight',type=float,default=.2);p.add_argument('--seconds',type=float,default=420);p.add_argument('--maxiter',type=int,default=160)
    args=p.parse_args();out=args.out
    if out.exists():raise ValueError('Refuse existing output path')
    out.parent.mkdir(parents=True,exist_ok=True)
    atomic_json(out.parent/(out.name+'_ST053_registration.json'),dict(id='ST053-Q',utc_before_fit=datetime.now(timezone.utc).isoformat(),parent_sha256=hashlib.sha256(Path(args.parent).read_bytes()).hexdigest(),
      original_physics_and_validation_gates='UNCHANGED',auxiliary_overrides={'moment_absolute_ratio':args.ratio,'moment_weight':args.weight,'moment_quadrature':[48,72],'profile_ratio':1.,'acceleration_ratio':1.,'core_anchor':.08,'peak_cap':'parent training pool maximum, initially active points only'},
      optimizer_maxiter=args.maxiter,wall_seconds=args.seconds,training_seed=9175350,validation_seeds=[9175391,9175392],structure_seed=9175393,
      note='Inherited run registration describes base class defaults; explicit ST053 overrides here are authoritative for new optimizer only, never for original scientific gates.',pde_validated=False))
    MomentObjective.ratio=args.ratio;MomentObjective.weight=args.weight;MomentObjective.output=out
    inherited.AlignedObjective=MomentObjective
    inherited.run(args.parent,out,ident='ST053-Q',edge=True,maxiter=args.maxiter,seconds=args.seconds,anchor=.08,pp=.003,axis=1.75,reference=args.reference,robust=True,cap=1.)
if __name__=='__main__':main()
