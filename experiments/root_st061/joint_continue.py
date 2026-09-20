"""ST061: actual method comparison with initial-collar sentinels and checkpoints."""
from __future__ import annotations
import sys,json,time,hashlib,argparse
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st060'))
from continuation import Study as PriorStudy,Tangent,PairedGrid,Store,atomic_json,atomic_bytes,filehash,binding as old_binding
from curvature_step import corrected_step,ResidualModel
sys.path.insert(0,str(Path(__file__).parent))
from direct_qcqp import build_subspace,solve_subproblem

class Study(PriorStudy):
    def __init__(self,parent,map):
        super().__init__(parent,True,map)
        r,z,t=np.meshgrid([.16,.22,.28058,.34,.40,.55],[-1.967,-1.93583,-1.905,1.905,1.93583,1.967],[.25,.28,.5,.72,.75],indexing='ij')
        self.sentinels=PairedGrid(self.m.f,r*r,z,t,self.m.fc)
        co,_=self.m.physical(self.zero)
        self.sentinel0=float(np.linalg.norm(self.sentinels.residual(self.sentinels.values(co)),axis=-1).max())
        self.peak0=max(self.peak0,self.sentinel0)
    def stats(self,c):
        s=super().stats(c);co,_=self.m.physical(c)
        p=float(np.linalg.norm(self.sentinels.residual(self.sentinels.values(co)),axis=-1).max())
        s['regular_pool_max']=s['dense_training_max'];s['sentinel_max']=p
        s['dense_training_max']=max(s['dense_training_max'],p);s['dense_peak_ratio']=s['dense_training_max']/self.peak0
        return s
    def local_model(self,c,mw=4.):
        H,g,v,A,p,Ap,meta=super().local_model(c,mw)
        tan=Tangent(self.m,c,self.M);R=self.sentinels.residual(self.sentinels.values(tan.co))[:,0,:];J,_=tan.residual_jac(self.sentinels)
        n=np.linalg.norm(R,axis=1);ids=np.argsort(n)[-36:];P=np.einsum('ni,nik->nk',R[ids],J[ids])/n[ids,None]
        meta['pool_max']=max(meta['pool_max'],float(n.max()))
        return H,g,v,A,np.r_[p,n[ids]],np.vstack((Ap,P)),meta

def direct_step(st,c,H,g,v,A,p,Ap,cap,damping,rank):
    D,info=build_subspace(H,g,A,v,Ap,p,cap,damping,rank)
    co,_=st.m.physical(c);R=st.o.Q.residual(st.o.Q.values(co))
    J,_=Tangent(st.m,c,st.M@D).residual_jac(st.o.Q)
    rm=ResidualModel(st,c);r=st.sentinels.residual(st.sentinels.values(co))[:,0,:]
    jp,_=Tangent(st.m,c,st.M@D).residual_jac(st.sentinels)
    ids=np.argsort(np.linalg.norm(r,axis=1))[-36:]
    Rp=np.r_[rm.Rp,r[ids]];Jp=np.concatenate((np.einsum('nik,kj->nij',rm.Jp,D),jp[ids]))
    weights=np.vstack((st.region_weights,st.o.sw[None,:]/64))
    targets=np.vstack((1.0001*st.region0,.9999*st.o.basetime[None,:]))
    delta,receipt=solve_subproblem(H,g,A,v,R,J,weights,targets,Rp,Jp,cap,D,damping)
    receipt['subspace']=info
    return delta,receipt

def binding():
    b=old_binding()
    for n in ['root_st060/curvature_step.py','root_st061/joint_continue.py','root_st061/direct_qcqp.py']:
        b[n]=filehash(ROOT/'experiments'/n)
    return b

def rawbytes(m,c,ident):
    obj=dict(schema='root_st001_compact_spacetime_v1',nr=m.f.nr,nz=m.f.nz,nt=m.f.nt,basis_kind=m.f.basis_kind,coefficients=m.raw_candidate(c).tolist(),metadata={'id':ident,'parent':'ST060-Q','scope':'Bounded joint continuation; original physical gates unchanged'},pde_validated=False,paper_exact=False,blowup_proved=False,field_identity_claim=False)
    return (json.dumps(obj,indent=2)+'\n').encode()

def run(parent,map,out,ident,method,outer=6,seconds=420,rank=96,resume=False,stop_after=None):
    out=Path(out);store=Store(out);src=binding();start=time.monotonic()
    if (out/'FROZEN.json').exists():raise ValueError('Frozen study; cannot resume fitting')
    if resume:
        b=store.verify();assert b['source_sha256']==src and b['parent_sha256']==filehash(Path(parent));reg=b['registration']
        ident,method,outer,seconds,rank=[reg[k] for k in ['id','method','outer','seconds','rank']]
        st=Study(out/'parent.json',np.load(out/'map.npy',allow_pickle=False));n,c,raw,state=store.load();assert rawbytes(st.m,c,ident)==raw
        history=state['history'];damping=state['damping'];spent=state['elapsed']
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('Existing output; explicit resume required')
        out.mkdir(parents=True,exist_ok=True)
        reg=dict(id=ident,method=method,outer=outer,seconds=seconds,rank=rank,registered_utc=datetime.now(timezone.utc).isoformat(),parent_sha256=filehash(Path(parent)),map_sha256=filehash(Path(map)),original_physical_gates='UNCHANGED',nu=.01,E0=1,time_window=[.25,.75],force_fixed=True,volume_ratios=[.99,1.02],new_holdouts=[9206291,9206292],new_structure_seed=9206293,sentinel_design='180 predetermined paired points; informed by prior peaks, not new holdouts',selection='last committed feasible iterate; both declared training max and worst-time L2 nonincreasing',scientific_target=.001)
        atomic_json(out/'PENDING.json',reg);st=Study(parent,np.load(map,allow_pickle=False));(out/'PENDING.json').unlink()
        store.initialize(dict(registration=reg,source_sha256=src,parent_sha256=filehash(Path(parent))),dict(map=st.M),Path(parent).read_bytes())
        n=0;c=st.zero.copy();history=[];damping=.03;spent=0.
        store.save(n,c,rawbytes(st.m,c,ident),dict(history=history,damping=damping,elapsed=0.,training=st.stats(c)))
    print('START',ident,'generation',n,'dim',st.M.shape,flush=True)
    base=st.stats(st.zero)
    for k in range(n,outer):
        H,g,v,A,p,Ap,meta=st.local_model(c,4.)
        old=st.stats(c);accepted=None;trials=[]
        for ratio in (.985,1.):
            cap=min(st.peak0*.999999,meta['pool_max']*ratio)
            if method=='direct':d,solver=direct_step(st,c,H,g,v,A,p,Ap,cap,damping,rank)
            else:d,solver=corrected_step(st,c,H,g,v,A,p,Ap,cap,damping)
            for alpha in [1.,.7,.5,.3,.1,.03,.01]:
                trial=c+alpha*st.M@d
                try:
                    s=st.stats(trial);good=st.feasible(s) and s['training_volume_L2_max']<old['training_volume_L2_max']-1e-9 and s['dense_training_max']<=old['dense_training_max']+1e-10
                except (ValueError,FloatingPointError):s={'invalid':True};good=False
                trials.append(dict(alpha=alpha,cap=cap,accepted=bool(good),stats=s,solver=solver))
                if good:accepted=trial;break
            if accepted is not None:break
        if accepted is None:damping*=5
        else:c=accepted;damping=max(.001,damping/1.5)
        stats=st.stats(c);history.append(dict(iteration=k+1,accepted=accepted is not None,trials=trials,training=stats))
        state=dict(history=history,damping=damping,elapsed=spent+time.monotonic()-start,training=stats)
        store.save(k+1,c,rawbytes(st.m,c,ident),state)
        print('STEP',ident,k+1,accepted is not None,'L2',stats['training_volume_L2_max'],'max',stats['dense_training_max'],'vol',stats['effective_volume_min_ratio'],flush=True)
        if stop_after is not None and k+1>=stop_after:return
        if state['elapsed']>seconds:break
    n,c,raw,state=store.load();atomic_bytes(out/'candidate.json',raw)
    atomic_json(out/'summary.json',dict(id=ident,iterations=n,elapsed=state['elapsed'],base=base,selected=state['training'],optimizer_converged=False,pde_validated=False,candidate_sha256=hashlib.sha256(raw).hexdigest()))
    print('FINISHED',ident,n,state['elapsed'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--map',required=True);p.add_argument('--out',required=True);p.add_argument('--id',required=True);p.add_argument('--method',choices=['direct','cuts'],required=True);p.add_argument('--outer',type=int,default=6);p.add_argument('--seconds',type=float,default=420);p.add_argument('--rank',type=int,default=96);p.add_argument('--resume',action='store_true');p.add_argument('--stop-after',type=int)
    a=p.parse_args();run(a.parent,a.map,a.out,a.id,a.method,a.outer,a.seconds,a.rank,a.resume,a.stop_after)
