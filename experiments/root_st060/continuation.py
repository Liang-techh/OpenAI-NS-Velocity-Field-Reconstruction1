"""ST060 joint temporal enrichment with restartable peak/volume constrained steps.
Same frozen physical family and prescribed force; finite training constraints only.
"""
from __future__ import annotations
import argparse, hashlib, json, sys, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.linalg import null_space, qr, cholesky, solve_triangular
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st059r'))
from bulk_step import BulkStudy,Tangent,PairedGrid,quadratic_step,source_binding as previous_binding
from checkpoint import Store,atomic_json,atomic_bytes,filehash

class Study(BulkStudy):
    def __init__(self,parent,temporal=True,saved_map=None):
        super().__init__(parent,True,saved_map)
        m,o=self.m,self.o
        if temporal and saved_map is None:
            # Add previously unused higher temporal modes of existing spatial columns.
            pairs=[(i,j,k) for i in range(3) for j in (0,1,2,3,9,10,11) for k in (3,4,5)]
            n=len(pairs);E=np.zeros((m.dim,3*n))
            for b,T in enumerate((m.f.Tp,m.f.Tw,m.f.Tq)):
                for q,(i,j,k) in enumerate(pairs):
                    v=np.zeros((m.f.nr,m.f.nz));v[i,j]=1;v=np.linalg.solve(T,v.ravel());v/=np.linalg.norm(v)
                    block=np.zeros((m.f.ns,m.f.nt));block[:,k]=v;E[b*m.f.n:(b+1)*m.f.n,b*n+q]=block.ravel()
            N=len(o.G.s);C=np.zeros((3*N,m.dim));n0=m.f.n
            C[:N,:n0]=o.G.r*o.G.D['A'];C[N:2*N,n0:2*n0]=o.G.r*o.G.D['B'];C[2*N:,:n0]=o.G.D['C']
            E=E@null_space(C@E,rcond=1e-10)
            # Orthonormalize the union, then apply one fixed Jacobian whitening.
            U,R,piv=qr(np.column_stack((self.M,E)),mode='economic',pivoting=True)
            rank=int(np.sum(abs(np.diag(R))>abs(R[0,0])*1e-11));U=U[:,:rank]
            rng=np.random.default_rng(9206151);D=m.grid(rng.uniform(0,4,180),rng.uniform(-2,2,180),np.linspace(.25,.75,11))
            J,_=Tangent(m,self.zero,U).residual_jac(D);J=J.reshape(-1,rank);H=J.T@J/(180*11)
            sc=1/np.sqrt(np.maximum(np.diag(H),1e-12));G=H*sc[:,None]*sc[None,:];G.flat[::rank+1]+=.005
            L=cholesky(G,lower=True);self.M=solve_triangular(L,(U*sc).T,lower=True).T
            self.map_info={'dimension':rank,'added_temporal_modes':[3,4,5],'added_raw_columns':3*n,'probe_null_max':float(abs(C@self.M).max())}
        # Finer collar training, no use of reserved new validation locations.
        r=np.unique(np.r_[np.linspace(0,1.995,37),np.linspace(.021,1.16,43)])
        z=np.unique(np.r_[np.linspace(-1.995,1.995,53),np.linspace(-1.988,-1.755,49),np.linspace(1.755,1.988,49)])
        R,Z=np.meshgrid(r,z,indexing='ij');self.P=m.grid(R.ravel()**2,Z.ravel(),np.linspace(.25,.75,35))
        co,_=m.physical(self.zero);self.peak0=float(np.linalg.norm(self.P.residual(self.P.values(co)),axis=-1).max());self.pool_count=len(self.P.s)*len(self.P.t)
        self.vol0=self.volumes(o.Q.values(co))
        self.blocks=[]
        for k in range(0,len(o.Q.s),64):
            idx=slice(k,k+64);self.blocks.append((idx,m.grid(o.Q.s[idx],o.Q.z[idx],o.times)))
    def densities(self,q):
        s=self.o.Q.s[:,None];e=s*(q['A']**2+q['B']**2)+q['C']**2
        om=s*q['Bz']**2+s*(q['Az']-2*q['Cs'])**2+4*(q['B']+s*q['Bs'])**2
        return e,om
    def volumes(self,q):
        w=self.o.sw;return np.array([(w@v)**2/(w@(v*v)) for v in self.densities(q)])
    def volume_jacobian(self,c,q):
        m,o=self.m,self.o;co,state=m.physical(c);s=o.Q.s[:,None];nt=len(o.times);rows=[]
        for kind,rho in enumerate(self.densities(q)):
            mass=o.sw@rho;square=o.sw@(rho*rho)
            k=2*o.sw[:,None]/mass-2*o.sw[:,None]*rho/square
            jets=({'A':2*s*q['A']*k,'B':2*s*q['B']*k,'C':2*q['C']*k} if kind==0 else
                  {'Bz':2*s*q['Bz']*k,'Az':2*s*(q['Az']-2*q['Cs'])*k,'Cs':-4*s*(q['Az']-2*q['Cs'])*k,
                   'B':8*(q['B']+s*q['Bs'])*k,'Bs':8*s*(q['B']+s*q['Bs'])*k})
            for t in range(nt):
                a={key:np.zeros_like(v) for key,v in jets.items()}
                for key,v in jets.items():a[key][:,t]=v[:,t]
                grad=m.pull_normalization(o.Q.pull(a),state);rows.append(grad@self.M)
        return np.array(rows)
    def stats(self,c):
        st=super().stats(c);co,_=self.m.physical(c);v=self.volumes(self.o.Q.values(co))/self.vol0
        st['effective_volume_min_ratio']=float(v.min());st['effective_volume_max_ratio']=float(v.max());st['effective_volume_ratios']=v.tolist()
        return st
    @staticmethod
    def feasible(st):
        return BulkStudy.feasible(st) and st['effective_volume_min_ratio']>=.99-1e-7 and st['effective_volume_max_ratio']<=1.02+1e-7
    def local_model(self,c,mw=4.):
        H,g,v,A,p,Ap,meta=super().local_model(c,mw);co,_=self.m.physical(c);q=self.o.Q.values(co)
        ratio=self.volumes(q)/self.vol0;J=self.volume_jacobian(c,q)
        v=np.r_[v,np.log(ratio.ravel()/.991),np.log(1.019/ratio.ravel())];A=np.vstack((A,J,-J))
        return H,g,v,A,p,Ap,meta

def binding():
    b=previous_binding();b['root_st060/continuation.py']=filehash(Path(__file__))
    return b

def candidate_bytes(m,c,ident):
    d=dict(schema='root_st001_compact_spacetime_v1',nr=m.f.nr,nz=m.f.nz,nt=m.f.nt,basis_kind=m.f.basis_kind,
           coefficients=m.raw_candidate(c).tolist(),metadata={'id':ident,'parent':'ST059R-V','scope':'Bounded numerical continuation; original physics unchanged'},
           pde_validated=False,paper_exact=False,blowup_proved=False,field_identity_claim=False)
    return (json.dumps(d,indent=2)+'\n').encode()

def run(parent,out,ident='ST060-T',outer=14,seconds=540,temporal=True,mw=4.,peak_ratio=.98,resume=False,stop_after=None):
    out=Path(out);store=Store(out);start=time.monotonic();src=binding()
    if (out/'FROZEN.json').exists():raise ValueError('Frozen for validation; never resume fitting this study')
    if resume:
        b=store.verify();reg=b['registration']
        if b['source_sha256']!=src or b['parent_sha256']!=filehash(Path(parent)):raise ValueError('Resume source or parent changed')
        ident,outer,seconds,temporal,mw,peak_ratio=[reg[k] for k in ('id','outer','seconds','temporal','moment_weight','local_peak_ratio')]
        st=Study(out/'parent.json',temporal,np.load(out/'map.npy',allow_pickle=False));n,c,raw,state=store.load()
        if candidate_bytes(st.m,c,ident)!=raw:raise ValueError('Candidate reconstruction mismatch')
        history=state['history'];damping=state['damping'];spent=state['elapsed']
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('Existing directory; explicit resume required')
        out.mkdir(parents=True,exist_ok=True)
        reg=dict(id=ident,outer=outer,seconds=seconds,temporal=temporal,moment_weight=mw,local_peak_ratio=peak_ratio,
                 registered_utc=datetime.now(timezone.utc).isoformat(),original_physical_gates='UNCHANGED',
                 E0=1.,nu=.01,time_window=[.25,.75],force_fixed=True,initial_velocity_frozen=False,
                 volume_retention=[.99,1.02],region_ratio_cap=1.0005,holdout_seeds=[9206191,9206192],structure_seed=9206193,
                 selection='nonlinear-feasible decreasing sum of normalized max and worst-time volumeL2, with both nonincreasing per step',
                 inherited_auxiliary='ST059R core/pressure/shear/drift/morphology/harmonic screens unchanged; new volume safeguards')
        atomic_json(out/'PENDING.json',reg);st=Study(parent,temporal);(out/'PENDING.json').unlink()
        store.initialize({'registration':reg,'parent_sha256':filehash(Path(parent)),'source_sha256':src,'map_info':st.map_info},{'map':st.M},Path(parent).read_bytes())
        c=st.zero.copy();history=[];damping=.03;n=0;spent=0.
        store.save(0,c,candidate_bytes(st.m,c,ident),dict(history=[],damping=damping,elapsed=0.,training=st.stats(c)))
    base=st.stats(st.zero);print('BASE',json.dumps({k:v for k,v in base.items() if not isinstance(v,list)}),flush=True);print('DIM',st.M.shape,flush=True)
    for k in range(n,outer):
        H,g,v,A,p,Ap,meta=st.local_model(c,mw);old=st.stats(c);accepted=None;trials=[]
        for ratio in (peak_ratio,.995,1.):
            cap=min(st.peak0*.999999,meta['pool_max']*ratio)
            step,inner=quadratic_step(H,g,v,A,p,Ap,cap,damping=damping,max_rounds=7)
            for alpha in (1.,.7,.5,.3,.1,.03,.01):
                trial=c+alpha*st.M@step
                try:s=st.stats(trial);feas=st.feasible(s)
                except (ValueError,FloatingPointError):s={'invalid':True};feas=False
                score=lambda x:x['training_volume_L2_max']/base['training_volume_L2_max']+x['dense_training_max']/base['dense_training_max']
                good=bool(feas and s['training_volume_L2_max']<=old['training_volume_L2_max']+1e-10 and s['dense_training_max']<=old['dense_training_max']+1e-10 and score(s)<score(old)-1e-7)
                trials.append(dict(alpha=alpha,peak_cap=cap,accepted=good,stats=s,inner=inner))
                if good:accepted=trial;break
            if accepted is not None:break
        if accepted is None:damping*=5
        else:c=accepted;damping=max(.001,damping/1.5)
        stats=st.stats(c);history.append(dict(iteration=k+1,accepted=accepted is not None,stats=stats,trials=trials))
        state=dict(history=history,damping=damping,elapsed=spent+time.monotonic()-start,training=stats)
        store.save(k+1,c,candidate_bytes(st.m,c,ident),state)
        print('STEP',k+1,accepted is not None,'L2',stats['training_volume_L2_max'],'max',stats['dense_training_max'],'volume',stats['effective_volume_min_ratio'],flush=True)
        if stop_after is not None and k+1>=stop_after:return {'stopped_after':k+1}
        if state['elapsed']>seconds:break
    n,c,raw,state=store.load();atomic_bytes(out/'candidate.json',raw)
    summary=dict(id=ident,iterations=n,elapsed=state['elapsed'],base=base,selected=state['training'],optimizer_converged=False,pde_validated=False,candidate_sha256=hashlib.sha256(raw).hexdigest())
    atomic_json(out/'summary.json',summary);return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--id',default='ST060-T');p.add_argument('--outer',type=int,default=14);p.add_argument('--seconds',type=float,default=540);p.add_argument('--no-temporal',action='store_true');p.add_argument('--mw',type=float,default=4.);p.add_argument('--peak-ratio',type=float,default=.98);p.add_argument('--resume',action='store_true');p.add_argument('--stop-after',type=int)
    a=p.parse_args();print(json.dumps(run(a.parent,a.out,a.id,a.outer,a.seconds,not a.no_temporal,a.mw,a.peak_ratio,a.resume,a.stop_after)),flush=True)
