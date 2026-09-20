"""ST059R: restartable joint volume minimization in the unchanged NS family.
Rebuild from available ST058-IP, not a fabricated recovery of missing ST059-E.
"""
from __future__ import annotations
import argparse,hashlib,json,os,sys,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import null_space,cholesky,solve_triangular
from scipy.special import softmax
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'experiments/root_st058'))
from constrained_continuation import Study as PriorStudy, Model, Objective, PairedGrid, Tangent, quadratic_step
from harmonic_audit import basis
from checkpoint import Store,atomic_json,atomic_bytes,filehash

SLABS=[(-2.,-1.9),(-1.9,-1.75),(-1.75,1.75),(1.75,1.9),(1.9,2.)]

def mapped_quadrature(nr=32,orders=(14,12,36,12,14)):
    x,wr=leggauss(nr);s=2*(x+1);S=[];Z=[];W=[];labels=[]
    for k,((a,b),n) in enumerate(zip(SLABS,orders)):
        z,w=leggauss(n);z=(a+b)/2+(b-a)*z/2
        ss,zz=np.meshgrid(s,z,indexing='ij')
        S.extend(ss.ravel());Z.extend(zz.ravel());W.extend((np.pi*(b-a)*np.outer(wr,w)).ravel());labels.extend([k]*ss.size)
    return np.array(S),np.array(Z),np.array(W),np.array(labels)

def make_directions(m,o,enriched=True):
    pairs=[(i,j,k) for i in range(3) for j in (0,1,2,3,7,8,9,10,11) for k in range(3)]
    pairs += [(i,j,k) for i in (0,1) for j in (7,8) for k in (3,4,5,6,7)]
    pairs += [(i,j,k) for i in (7,8) for j in (0,1) for k in (0,1,2)]
    if enriched:
        pairs += [(i,j,k) for i in range(3) for j in (4,5,6) for k in range(3)]
        pairs += [(i,j,k) for i in (3,4,5) for j in range(7) for k in range(3)]
    maps=[]
    for transform in (m.f.Tp,m.f.Tw,m.f.Tq):
        columns=[]
        for i,j,k in pairs:
            q=np.zeros((m.f.nr,m.f.nz));q[i,j]=1.;v=np.linalg.solve(transform,q.ravel());v/=np.linalg.norm(v)
            a=np.zeros((m.f.ns,m.f.nt));a[:,k]=v;columns.append(a.ravel())
        maps.append(np.column_stack(columns))
    n=maps[0].shape[1];M=np.zeros((m.dim,3*n))
    for i in range(3):M[i*m.f.n:(i+1)*m.f.n,i*n:(i+1)*n]=maps[i]
    C=np.zeros((3*len(o.G.s),m.dim));N=len(o.G.s);n=m.f.n
    C[:N,:n]=o.G.r*o.G.D['A'];C[N:2*N,n:2*n]=o.G.r*o.G.D['B'];C[2*N:,:n]=o.G.D['C']
    M=M@null_space(C@M,rcond=1e-10)
    return M,{'pairs_per_block':len(pairs),'dimension':M.shape[1],'probe_null_error':float(abs(C@M).max())}

class BulkStudy(PriorStudy):
    def __init__(self,parent,enriched=True,saved_map=None):
        self.m=Model(parent);m=self.m;self.o=Objective(m,moment_weight=4.,safety=100.,pool_seed=9206050);o=self.o
        s,z,w,labels=mapped_quadrature();self.labels=labels
        o.Q=m.grid(s,z,o.times);o.sw=w;o.qw=w[:,None]*o.tw/64
        self.zero=np.zeros(m.dim);co,_=m.physical(self.zero)
        o.baseq=o.Q.values(co);o.baseres=o.Q.residual(o.baseq);o.basetime=w@np.sum(o.baseres**2,axis=-1)/64
        _,Gram,hess=basis();norm=np.sqrt(np.diag(Gram));L=cholesky(Gram/norm[:,None]/norm[None,:],lower=True)
        hh=hess(s,z);HH=np.array([[np.broadcast_to(v,(len(s),)) for v in row] for row in hh]);HH[:,2,:]*=2*np.sqrt(s)
        o.H=solve_triangular(L,HH.reshape(7,-1)/norm[:,None],lower=True).reshape(HH.shape)
        o.baseharm=o.harmonic(o.baseq)[0];o.basemorph=o.morph_values(o.baseq)[0];o.last=None
        self.q0=o.baseq;self.d0=o.baseharm;self.cv0=o.C.values(co);self.sv0=o.S.values(co)
        r=np.unique(np.r_[np.linspace(0,1.992,33),np.linspace(.025,1.15,34)])
        zz=np.unique(np.r_[np.linspace(-1.992,1.992,49),np.linspace(-1.985,-1.76,39),np.linspace(1.76,1.985,39)])
        R,Z=np.meshgrid(r,zz,indexing='ij');self.P=m.grid(R.ravel()**2,Z.ravel(),np.linspace(.25,.75,31))
        self.peak0=float(np.linalg.norm(self.P.residual(self.P.values(co)),axis=-1).max());self.pool_count=len(self.P.s)*len(self.P.t)
        self.region_weights=np.array([w*(labels==k)/64 for k in range(5)])
        self.region0=self.region_weights@np.sum(o.baseres**2,axis=-1)
        if saved_map is None:
            M,self.map_info=make_directions(m,o,enriched)
            rng=np.random.default_rng(9206051);D=m.grid(rng.uniform(0,4,170),rng.uniform(-2,2,170),np.linspace(.25,.75,11))
            J,_=Tangent(m,self.zero,M).residual_jac(D);J=J.reshape(-1,M.shape[1]);H=J.T@J/(170*11)
            sc=1/np.sqrt(np.maximum(np.diag(H),1e-12));G=H*sc[:,None]*sc[None,:];G.flat[::len(G)+1]+=.005
            L=cholesky(G,lower=True);self.M=solve_triangular(L,(M*sc).T,lower=True).T
        else:
            self.M=np.asarray(saved_map).copy();self.map_info={'dimension':self.M.shape[1],'loaded_exact_map':True}
        self.blocks=[]
        for k in range(0,len(s),128):
            idx=slice(k,k+128);self.blocks.append((idx,m.grid(s[idx],z[idx],o.times)))
    def stats(self,c):
        st=super().stats(c);co,_=self.m.physical(c);R=self.o.Q.residual(self.o.Q.values(co))
        sq=np.sum(R*R,axis=-1);reg=self.region_weights@sq
        st['region_ratio_max']=float((reg/self.region0).max());st['region_ratios']=(reg/self.region0).tolist()
        return st
    @staticmethod
    def feasible(st):
        return PriorStudy.feasible(st) and st['region_ratio_max']<=1.0005
    def quadrature_model(self,c,mw):
        m,o,M=self.m,self.o,self.M;k=M.shape[1];tan=Tangent(m,c,M);q=o.Q.values(tan.co);R=o.Q.residual(q)
        mse=o.sw@np.sum(R*R,axis=-1)/64
        weights=.25*o.tw+.75*softmax(8*mse/max(mse.max(),1e-30)+np.log(o.tw))
        H=np.zeros((k,k));g=np.zeros(k);dm=np.zeros((len(o.times),k));dregion=np.zeros((5,len(o.times),k))
        dh=np.zeros((7,len(o.times),k));dn=np.zeros((4,len(o.times),k))
        hd,(ur,ut,uz)=o.harmonic(q);_,(D0,E0,en,nums,factors)=o.morph_values(q)
        for idx,grid in self.blocks:
            vb={key:a[idx] for key,a in q.items()};J,d=tan.residual_jac(grid,vb);w=o.sw[idx];W=w[:,None]*weights/64
            flat=J.reshape(-1,k);H+=flat.T@(J*W[:,:,None,None]).reshape(-1,k)
            g+=np.einsum('ntij,nti,nt->j',J,R[idx],W,optimize=True)
            dsq=2*np.einsum('nti,ntij->ntj',R[idx],J,optimize=True)
            dm+=np.einsum('n,ntj->tj',w/64,dsq)
            dregion+=np.einsum('an,ntj->atj',self.region_weights[:,idx],dsq)
            h=o.H[:,:,idx];dr=grid.r[...,None]*d['A'];dt=grid.r[...,None]*d['B'];dz=d['C']
            dh-=2*np.einsum('n,kn,nt,ntj->ktj',w,h[:,0],ur[idx],dr,optimize=True)
            dh-=2*np.einsum('n,kn,nt,ntj->ktj',w,h[:,1],ut[idx],dt,optimize=True)
            dh-=2*np.einsum('n,kn,nt,ntj->ktj',w,h[:,2],uz[idx],dr,optimize=True)
            dh-=2*np.einsum('n,kn,nt,ntj->ktj',w,h[:,2],ur[idx],dz,optimize=True)
            dh-=2*np.einsum('n,kn,nt,ntj->ktj',w,h[:,3],uz[idx],dz,optimize=True)
            ss=grid.s[:,None,None]
            de=2*ss*q['Bz'][idx,...,None]*d['Bz']+2*ss*D0[idx,...,None]*(d['Az']-2*d['Cs'])+8*E0[idx,...,None]*(d['B']+ss*d['Bs'])
            dn+=np.einsum('an,ntj->atj',factors[:,idx]*w,de,optimize=True)
        H+=mw*np.einsum('ktj,kti,t->ji',dh,dh,o.tw,optimize=True)/64
        g+=mw*np.einsum('kt,ktj,t->j',hd,dh,o.tw,optimize=True)/64
        return tan,q,H,g,mse,dm,hd,dh,dn,dregion,nums
    def local_model(self,c,mw=4.):
        m,o,M=self.m,self.o,self.M;k=M.shape[1]
        tan,q,H,g,mse,dm,hd,dh,dn,dregion,nums=self.quadrature_model(c,mw)
        co=tan.co;rows=[];vals=[]
        def push(v,j):vals.append(np.asarray(v).ravel());rows.append(np.asarray(j).reshape(-1,k))
        push(.98*o.basetime-mse,-dm)
        push(.98**2*np.sum(self.d0**2,axis=0)-np.sum(hd**2,axis=0),-2*np.einsum('kt,ktj->tj',hd,dh))
        Rq=o.Q.residual(q);reg=self.region_weights@np.sum(Rq*Rq,axis=-1)
        push(.999*self.region0-reg,-dregion)
        # Pointwise pressure/velocity/shear guards with exact normalized derivatives.
        cv=o.C.values(co);dc=tan.jets(o.C,('A','B','C','Qs','Qz'))
        for i,key in enumerate(('A','B','C')):
            factor=o.C.r[:,0] if i<2 else np.ones(len(o.C.s));v=factor*cv[key][:,0];jac=factor[:,None]*dc[key]
            base=o.basecore[:,0,i];band=.195*o.cscale[:,0,i]
            push(v-(base-band),jac);push(base+band-v,-jac)
            sign=-np.ones_like(v) if i==0 else np.ones_like(v) if i==1 else o.csign[:,0]
            push(sign*v-1e-7,sign[:,None]*jac)
        for key,sgn,margin in [('Qs',np.ones(len(o.C.s)),1e-6),('Qz',o.csign[:,0],1e-5)]:
            lower=np.minimum(margin,.95*sgn*self.cv0[key][:,0]);push(sgn*cv[key][:,0]-lower,sgn[:,None]*dc[key])
        sv=o.S.values(co);ds=tan.jets(o.S,('Cs','C'));fac=2*o.S.r[:,0]*np.sign(o.baseshear[:,0])
        push(fac*sv['Cs'][:,0]-.997*abs(o.baseshear[:,0]),fac[:,None]*ds['Cs'])
        push(sv['C'][:,0]-.85*o.basebias[:,0],ds['C'])
        U=(o.velocity(cv,o.C)*o.cscales).reshape(o.core_n,o.core_nt,3)
        DUC=np.stack((o.C.r*dc['A'],o.C.r*dc['B'],dc['C']),axis=1)*o.cscales[:,0,:,None]
        DU=DUC.reshape(o.core_n,o.core_nt,3,k);dif=U-U[:,0:1];den=np.sum(U[:,0]**2)
        for t in range(1,o.core_nt):
            cap=o.parent_core_drift[t]+.0003
            v=cap**2*den-np.sum(dif[:,t]**2)
            jac=2*cap**2*np.einsum('ni,nij->j',U[:,0],DU[:,0])-2*np.einsum('ni,nij->j',dif[:,t],DU[:,t]-DU[:,0]);push(v,jac)
        for i in range(3):
            ratio=.995 if i<2 else 1.0105;sign=1 if i<2 else -1;target=ratio*o.basemorph[i]
            push(sign*(nums[i+1]-target*nums[0]),sign*(dn[i+1]-target[:,None]*dn[0]))
        # Current active peaks from every training time, not the published holdouts.
        Rp=self.P.residual(self.P.values(co));norm=np.linalg.norm(Rp,axis=-1)
        ids=set(np.argpartition(norm.ravel(),-128)[-128:].tolist())
        for t in range(norm.shape[1]):
            for i in np.argpartition(norm[:,t],-5)[-6:]:ids.add(int(i)*norm.shape[1]+t)
        ids=np.array(sorted(ids));sp,ti=np.unravel_index(ids,norm.shape)
        pg=PairedGrid(m.f,self.P.s[sp],self.P.z[sp],self.P.t[ti],m.fc);Jp,_=tan.residual_jac(pg);rp=Rp.reshape(-1,3)[ids]
        npv=np.linalg.norm(rp,axis=1);Apeak=np.einsum('ni,nij->nj',rp,Jp)/npv[:,None]
        # A small direct least-squares cost also focuses on active peak directions.
        wt=.002/len(ids);H+=wt*Jp.reshape(-1,k).T@Jp.reshape(-1,k);g+=wt*np.einsum('nij,ni->j',Jp,rp)
        return H,g,np.concatenate(vals),np.vstack(rows),npv,Apeak,dict(pool_max=float(norm.max()),active_peaks=len(ids))

def source_binding():
    paths=['root_st059r/bulk_step.py','root_st059r/checkpoint.py',
           'root_st058/constrained_continuation.py','root_st057/coupled_fit.py',
           'root_st057/tensor_model.py','root_st030/spacetime.py','root_st030/asymmetric_basis.py',
           'root_st030/hybrid_basis.py','root_st056/harmonic_audit.py']
    return {p:filehash(ROOT/'experiments'/p) for p in paths}

def candidate_bytes(m,c,ident):
    raw=m.raw_candidate(c)
    obj=dict(schema='root_st001_compact_spacetime_v1',nr=m.f.nr,nz=m.f.nz,nt=m.f.nt,
             basis_kind=m.f.basis_kind,coefficients=raw.tolist(),metadata={'id':ident,'scope':'Recovery refit from ST058-IP; no historical ST059-E identity claim'},
             pde_validated=False,paper_exact=False,blowup_proved=False,field_identity_claim=False)
    return (json.dumps(obj,indent=2)+'\n').encode()

def run(parent,out,ident='ST059R-V',outer=12,seconds=600,enriched=True,mw=4.,resume=False,stop_after=None):
    out=Path(out);store=Store(out);start=time.monotonic();src=source_binding()
    if resume:
        binding=store.verify()
        if binding['source_sha256']!=src:raise ValueError('Source changed since checkpoint; use a new registered run')
        if filehash(Path(parent))!=binding['parent_sha256']:raise ValueError('Wrong resume parent')
        reg=binding['registration'];outer=reg['outer'];seconds=reg['seconds'];enriched=reg['enriched'];mw=reg['moment_weight'];ident=reg['id']
        st=BulkStudy(out/'parent.json',enriched,np.load(out/'map.npy',allow_pickle=False))
        completed,c,rawbytes,state=store.load();damping=state['damping'];history=state['history'];spent=state['elapsed']
        if candidate_bytes(st.m,c,ident)!=rawbytes:raise ValueError('Stored candidate/delta mismatch')
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('Existing run; use --resume or a new path')
        out.mkdir(parents=True,exist_ok=True)
        reg=dict(id=ident,outer=outer,seconds=seconds,enriched=enriched,moment_weight=mw,
                 registered_utc=datetime.now(timezone.utc).isoformat(),original_physical_gates='UNCHANGED',
                 nu=.01,time_window=[.25,.75],E0=1,force_fixed=True,initial_field_frozen=False,
                 holdout_seeds=[9206091,9206092],structure_seed=9206093,peak_times=31,
                 quadrature_radial=32,quadrature_slab_orders=[14,12,36,12,14],
                 five_signed_slabs=SLABS,slab_max_ratio=1.0005,
                 objective='worst-time reweighted integrated full residual plus harmonic penalty',
                 selection='decreasing training worst-time L2; original physics and finite nonlinear safeguards')
        # Registration is durable before model construction or optimization.
        atomic_json(out/'REGISTRATION_PENDING.json',reg)
        st=BulkStudy(parent,enriched)
        (out/'REGISTRATION_PENDING.json').unlink()
        store.initialize({'parent_sha256':filehash(Path(parent)),'source_sha256':src,'registration':reg}, {'map':st.M},Path(parent).read_bytes())
        c=st.zero.copy();history=[];damping=.03;completed=0;spent=0.
        state=dict(damping=damping,history=history,elapsed=0.,training=st.stats(c))
        store.save(0,c,candidate_bytes(st.m,c,ident),state)
    base=st.stats(st.zero)
    print('RECOVERED',ident,'completed',completed,'dim',st.M.shape[1],flush=True)
    for k in range(completed,outer):
        H,g,v,A,p,Ap,meta=st.local_model(c,mw);accepted=None;trials=[];old=st.stats(c)
        for ratio in (1.,1.002):
            cap=min(st.peak0*.999999,meta['pool_max']*ratio)
            step,inner=quadratic_step(H,g,v,A,p,Ap,cap,damping=damping,max_rounds=7)
            for alpha in (1.,.7,.5,.3,.1,.03,.01):
                trial=c+alpha*st.M@step
                try:stats=st.stats(trial);feasible=st.feasible(stats)
                except (ValueError,FloatingPointError):stats={'invalid':True};feasible=False
                good=bool(feasible and stats['training_volume_L2_max']<old['training_volume_L2_max']-1e-9)
                trials.append(dict(alpha=alpha,cap=cap,accepted=good,stats=stats,inner=inner))
                if good:accepted=trial;break
            if accepted is not None:break
        if accepted is None:damping*=5.
        else:c=accepted;damping=max(.001,damping/1.5)
        state=dict(damping=damping,elapsed=spent+time.monotonic()-start,training=st.stats(c))
        history.append(dict(iteration=k+1,accepted=accepted is not None,training=state['training'],trials=trials))
        state['history']=history
        store.save(k+1,c,candidate_bytes(st.m,c,ident),state)
        print('STEP',k+1,'accepted',accepted is not None,'L2',state['training']['training_volume_L2_max'],'max',state['training']['dense_training_max'],'region',state['training']['region_ratio_max'],flush=True)
        if stop_after is not None and k+1>=stop_after:
            # Real process boundary for recovery regression, not a copied vector.
            return dict(stopped_for_test=True,iteration=k+1)
        if state['elapsed']>seconds:break
    n,c,raw,state=store.load()
    atomic_bytes(out/'candidate.json',raw)
    atomic_json(out/'summary.json',dict(id=ident,iterations=n,selected=state['training'],base=base,
                elapsed=state['elapsed'],optimizer_converged=False,pde_validated=False,
                candidate_sha256=hashlib.sha256(raw).hexdigest(),recovery_scope='New refit; unavailable historical ST059-E not reused'))
    return json.loads((out/'summary.json').read_text())

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--parent',required=True);a.add_argument('--out',required=True)
    a.add_argument('--id',default='ST059R-V');a.add_argument('--outer',type=int,default=12);a.add_argument('--seconds',type=float,default=600)
    a.add_argument('--basic',action='store_true');a.add_argument('--mw',type=float,default=4.);a.add_argument('--resume',action='store_true');a.add_argument('--stop-after',type=int)
    v=a.parse_args();print(json.dumps(run(v.parent,v.out,v.id,v.outer,v.seconds,not v.basic,v.mw,v.resume,v.stop_after)),flush=True)
