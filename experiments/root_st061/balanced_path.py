"""Non-componentwise-monotone merit path, with unchanged parent-based hard screens.
Both final training norms stay no worse than the original parent. An individual
iteration may increase one norm slightly when their normalized sum decreases.
This is NOT an implementation or convergence proof of published SQP-filter methods.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from datetime import datetime,timezone
from joint_continue import ROOT,Study,Store,atomic_json,atomic_bytes,filehash,binding,rawbytes,direct_step

def run(parent,map,warm,out,outer=6,seconds=300,resume=False,stop_after=None):
    out=Path(out);store=Store(out);src=binding();src['root_st061/balanced_path.py']=filehash(Path(__file__));start=time.monotonic()
    if (out/'FROZEN.json').exists():raise ValueError('Frozen study cannot be refitted')
    ident='ST061-P'
    if resume:
        b=store.verify();assert b['source_sha256']==src and b['parent_sha256']==filehash(Path(parent));reg=b['registration'];outer=reg['outer'];seconds=reg['seconds']
        st=Study(out/'parent.json',np.load(out/'map.npy',allow_pickle=False));n,c,raw,state=store.load();assert rawbytes(st.m,c,ident)==raw;history=state['history'];damping=state['damping'];spent=state['elapsed']
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('Refuse existing directory')
        out.mkdir(parents=True,exist_ok=True)
        reg=dict(id=ident,registered_utc=datetime.now(timezone.utc).isoformat(),outer=outer,seconds=seconds,warm_sha256=filehash(Path(warm)),map_sha256=filehash(Path(map)),original_physical_gates='UNCHANGED',auxiliary_nonlinear_screens='UNCHANGED, relative to ST060-Q',merit='sum of parent-normalized worst-time L2 and maximum; strict decrease',individual_iteration_norm_increase_cap=1.005,global_parent_bounds='both <=1, per-time and regional inherited checks retained',subspace_rank=96,holdouts=[9206291,9206292],design_scope='training-only sequential redesign after monotonic-path stalls; no held-out data viewed')
        atomic_json(out/'PENDING.json',reg);st=Study(parent,np.load(map,allow_pickle=False));(out/'PENDING.json').unlink()
        c=np.load(warm,allow_pickle=False);assert st.feasible(st.stats(c));n=0;history=[];damping=.03;spent=0
        store.initialize(dict(registration=reg,source_sha256=src,parent_sha256=filehash(Path(parent))),dict(map=st.M,warm=c),Path(parent).read_bytes());store.save(0,c,rawbytes(st.m,c,ident),dict(history=[],damping=damping,elapsed=0.,training=st.stats(c)))
    base=st.stats(st.zero)
    def score(s):return s['training_volume_L2_max']/base['training_volume_L2_max']+s['dense_training_max']/base['dense_training_max']
    print('START P',n,'initial',st.stats(c)['training_volume_L2_max'],st.stats(c)['dense_training_max'],flush=True)
    for k in range(n,outer):
        H,g,v,A,p,Ap,meta=st.local_model(c,4.);old=st.stats(c);accepted=None;trials=[]
        for ratio in (.985,1.):
            cap=min(st.peak0*.999999,meta['pool_max']*ratio)
            d,solver=direct_step(st,c,H,g,v,A,p,Ap,cap,damping,96)
            for alpha in [1.,.7,.5,.3,.1,.03,.01]:
                trial=c+alpha*st.M@d
                try:
                    s=st.stats(trial);good=st.feasible(s) and score(s)<score(old)-1e-7 and s['training_volume_L2_max']<=1.005*old['training_volume_L2_max'] and s['dense_training_max']<=1.005*old['dense_training_max']
                except (ValueError,FloatingPointError):s={'invalid':True};good=False
                trials.append(dict(alpha=alpha,cap=cap,accepted=bool(good),stats=s,solver=solver))
                if good:accepted=trial;break
            if accepted is not None:break
        if accepted is None:damping*=5
        else:c=accepted;damping=max(.001,damping/1.5)
        s=st.stats(c);history.append(dict(iteration=k+1,accepted=accepted is not None,trials=trials,training=s))
        state=dict(history=history,damping=damping,elapsed=spent+time.monotonic()-start,training=s);store.save(k+1,c,rawbytes(st.m,c,ident),state)
        print('STEP P',k+1,accepted is not None,'L2',s['training_volume_L2_max'],'max',s['dense_training_max'],'merit',score(s),flush=True)
        if stop_after is not None and k+1>=stop_after:return
        if state['elapsed']>seconds:break
    n,c,raw,state=store.load();atomic_bytes(out/'candidate.json',raw);atomic_json(out/'summary.json',dict(id=ident,iterations=n,elapsed=state['elapsed'],base=base,selected=state['training'],optimizer_converged=False,pde_validated=False,candidate_sha256=filehash(out/'candidate.json')))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--map',required=True);p.add_argument('--warm',required=True);p.add_argument('--out',required=True);p.add_argument('--outer',type=int,default=6);p.add_argument('--seconds',type=float,default=300);p.add_argument('--resume',action='store_true');p.add_argument('--stop-after',type=int);a=p.parse_args();run(a.parent,a.map,a.warm,a.out,a.outer,a.seconds,a.resume,a.stop_after)
