"""Convex supporting cuts for norms of the LINEARIZED full residual.
The cut model includes quadratic residual-energy growth. It remains only a
local model; accepted updates must pass the unchanged nonlinear safeguards.
"""
from __future__ import annotations
import argparse,json,sys,time,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from continuation import Study,Tangent,PairedGrid,ROOT,Store,atomic_json,atomic_bytes,filehash,binding as base_binding,candidate_bytes
from constrained_continuation import quadratic_step

class ResidualModel:
    def __init__(self,st,c):
        self.st=st;self.c=c;self.m=st.m;self.o=st.o
        self.co,self.state=self.m.physical(c);self.q=self.o.Q.values(self.co);self.R=self.o.Q.residual(self.q)
        Rp=st.P.residual(st.P.values(self.co));norm=np.linalg.norm(Rp,axis=-1)
        ids=set(np.argpartition(norm.ravel(),-192)[-192:].tolist())
        for t in range(norm.shape[1]):
            for i in np.argpartition(norm[:,t],-7)[-7:]:ids.add(int(i)*norm.shape[1]+t)
        ids=np.array(sorted(ids));sp,ti=np.unravel_index(ids,norm.shape)
        P=PairedGrid(self.m.f,st.P.s[sp],st.P.z[sp],st.P.t[ti],self.m.fc)
        self.Jp,_=Tangent(self.m,c,st.M).residual_jac(P);self.Rp=Rp.reshape(-1,3)[ids]
    def residual(self,d):
        Jd,_=Tangent(self.m,self.c,(self.st.M@d)[:,None]).residual_jac(self.o.Q,self.q)
        return self.R+Jd[...,0]
    def pull_residual(self,W):
        adj=self.o.Q.residual_adjoints(self.q,W)
        g=self.m.pull_normalization(self.o.Q.pull(adj),self.state)
        return g@self.st.M
    def cuts(self,d,cap,max_energy_cuts=12):
        R=self.residual(d);sq=np.sum(R*R,axis=-1);rows=[];values=[]
        reg=self.st.region_weights@sq;target=1.0001*self.st.region0
        deficits=reg/target-1;ids=np.flatnonzero(deficits.ravel()>2e-6)
        ids=ids[np.argsort(deficits.ravel()[ids])[-max_energy_cuts:]]
        for idx in ids:
            region,t=np.unravel_index(idx,reg.shape);W=np.zeros_like(R);W[:,t,:]=2*self.st.region_weights[region,:,None]*R[:,t,:]
            grad=self.pull_residual(W);rows.append(-grad);values.append(target[region,t]-reg[region,t]+grad@d)
        mse=self.o.sw@sq/64;target_t=.9999*self.o.basetime
        for t in np.flatnonzero(mse>target_t*(1+2e-6)):
            W=np.zeros_like(R);W[:,t,:]=2*self.o.sw[:,None]*R[:,t,:]/64
            grad=self.pull_residual(W);rows.append(-grad);values.append(target_t[t]-mse[t]+grad@d)
        rp=self.Rp+np.einsum('nik,k->ni',self.Jp,d);norm=np.linalg.norm(rp,axis=1)
        ids=np.flatnonzero(norm>cap+1e-7);ids=ids[np.argsort(norm[ids])[-40:]]
        for idx in ids:
            grad=rp[idx]@self.Jp[idx]/norm[idx];rows.append(-grad);values.append(cap-norm[idx]+grad@d)
        info=dict(max_linearized_region_ratio=float((reg/self.st.region0).max()),max_linearized_time_mse_ratio=float((mse/self.o.basetime).max()),max_linearized_peak=float(norm.max()),new_cuts=len(rows))
        return np.array(values),np.array(rows).reshape(-1,len(d)),info

def corrected_step(st,c,H,g,v,A,p,Ap,cap,damping):
    rm=ResidualModel(st,c);rows=A.copy();vals=v.copy();logs=[];d=np.zeros(st.M.shape[1])
    for k in range(5):
        d,qp=quadratic_step(H,g,vals,rows,p,Ap,cap,damping=damping,max_rounds=6)
        vv,jj,info=rm.cuts(d,cap);logs.append(dict(cut_iteration=k,quadratic_model=info,qp=qp))
        if not len(vv):break
        rows=np.vstack((rows,jj));vals=np.r_[vals,vv]
    return d,logs

def run(parent,warm,map_path,out,outer=8,seconds=360,resume=False,stop_after=None):
    out=Path(out);store=Store(out);src=base_binding();src['root_st060/curvature_step.py']=filehash(Path(__file__))
    if (out/'FROZEN.json').exists():raise ValueError('Already frozen for holdouts')
    start=time.monotonic()
    if resume:
        b=store.verify();assert b['source_sha256']==src and b['parent_sha256']==filehash(Path(parent));reg=b['registration'];outer=reg['outer'];seconds=reg['seconds']
        st=Study(out/'parent.json',True,np.load(out/'map.npy',allow_pickle=False));n,c,raw,state=store.load();assert candidate_bytes(st.m,c,'ST060-Q')==raw
        history=state['history'];damping=state['damping'];spent=state['elapsed']
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('Refuse existing run')
        out.mkdir(parents=True,exist_ok=True)
        reg=dict(id='ST060-Q',registered_utc=datetime.now(timezone.utc).isoformat(),outer=outer,seconds=seconds,warm_delta_sha256=filehash(Path(warm)),map_sha256=filehash(Path(map_path)),original_physical_gates='UNCHANGED',holdouts=[9206191,9206192],constraint_screen='same ST060 nonlinear guards including .99..1.02 participation volume',method='supporting-plane outer approximation of convex linearized-residual quadratic norm constraints; nonlinear backtracking',cut_rounds=5)
        atomic_json(out/'PENDING.json',reg);st=Study(parent,True,np.load(map_path,allow_pickle=False));(out/'PENDING.json').unlink()
        store.initialize(dict(registration=reg,source_sha256=src,parent_sha256=filehash(Path(parent))),dict(map=st.M,warm_delta=np.load(warm,allow_pickle=False)),Path(parent).read_bytes())
        c=np.load(warm,allow_pickle=False);n=0;history=[];damping=.03;spent=0.
        assert st.feasible(st.stats(c));store.save(0,c,candidate_bytes(st.m,c,'ST060-Q'),dict(history=[],damping=damping,elapsed=0.,training=st.stats(c)))
    base=st.stats(st.zero);print('BEGIN',n,'dim',st.M.shape,flush=True)
    for k in range(n,outer):
        H,g,v,A,p,Ap,meta=st.local_model(c,4.);old=st.stats(c);accepted=None;trials=[]
        for ratio in (.985,1.):
            cap=min(st.peak0*.999999,meta['pool_max']*ratio)
            d,cuts=corrected_step(st,c,H,g,v,A,p,Ap,cap,damping)
            for alpha in (1.,.7,.5,.3,.1,.03,.01):
                trial=c+alpha*st.M@d
                try:s=st.stats(trial);good=st.feasible(s) and s['training_volume_L2_max']<old['training_volume_L2_max']-1e-9 and s['dense_training_max']<=old['dense_training_max']+1e-10
                except (ValueError,FloatingPointError):s={'invalid':True};good=False
                trials.append(dict(alpha=alpha,cap=cap,accepted=bool(good),stats=s,cuts=cuts))
                if good:accepted=trial;break
            if accepted is not None:break
        if accepted is None:damping*=5
        else:c=accepted;damping=max(.001,damping/1.5)
        s=st.stats(c);history.append(dict(iteration=k+1,accepted=accepted is not None,trials=trials,training=s))
        state=dict(history=history,damping=damping,elapsed=spent+time.monotonic()-start,training=s)
        store.save(k+1,c,candidate_bytes(st.m,c,'ST060-Q'),state);print('STEP',k+1,accepted is not None,s['training_volume_L2_max'],s['dense_training_max'],flush=True)
        if stop_after is not None and k+1>=stop_after:return
        if state['elapsed']>seconds:break
    n,c,raw,state=store.load();atomic_bytes(out/'candidate.json',raw)
    atomic_json(out/'summary.json',dict(id='ST060-Q',iterations=n,elapsed=state['elapsed'],base=base,selected=state['training'],optimizer_converged=False,pde_validated=False,candidate_sha256=hashlib.sha256(raw).hexdigest()))

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--parent',required=True);a.add_argument('--warm',required=True);a.add_argument('--map',required=True);a.add_argument('--out',required=True);a.add_argument('--outer',type=int,default=8);a.add_argument('--seconds',type=float,default=360);a.add_argument('--resume',action='store_true');a.add_argument('--stop-after',type=int)
    p=a.parse_args();run(p.parent,p.warm,p.map,p.out,p.outer,p.seconds,p.resume,p.stop_after)
