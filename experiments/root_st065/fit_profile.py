"""ST065: solve actual similarity-profile momentum instead of rescaling an unfitted snapshot.
The power laws remain imposed. No optimizer or geometry success implies PDE validity.
"""
from pathlib import Path
import sys, json, time, hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import qr, solve_triangular
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'baseline'))
from profile import basis, Profile
from scale_model import ScaleField

def atomic(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    t=path.with_suffix(path.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');t.replace(path)

def quad(n):
    x,w=leggauss(n);S,Z=np.meshgrid(2*(x+1),2*x,indexing='ij')
    return S.ravel(),Z.ravel(),(4*np.pi*np.outer(w,w)).ravel()

class Matrix:
    def __init__(self,s,z,transforms=None):
        self.s=np.asarray(s).ravel();self.z=np.asarray(z).ravel();s=self.s;z=self.z
        r=basis(s,'r');f=basis(z,'z',1);b=basis(z,'z',0)
        def prod(i,j,pol=True):
            return (r[i,:,:,None]*(f if pol else b)[j,:,None,:]).reshape(len(s),108)
        S=s[:,None]
        self.D=dict(A=-prod(0,1),As=-prod(1,1),Az=-prod(0,2),Ass=-prod(2,1),Azz=-prod(0,3),
          C=2*prod(0,0)+2*S*prod(1,0),Cs=4*prod(1,0)+2*S*prod(2,0),Cz=2*prod(0,1)+2*S*prod(1,1),
          Css=6*prod(2,0)+2*S*prod(3,0),Czz=2*prod(0,2)+2*S*prod(1,2),
          B=prod(0,0,False),Bs=prod(1,0,False),Bz=prod(0,1,False),Bss=prod(2,0,False),Bzz=prod(0,2,False),
          P=prod(0,0,False),Ps=prod(1,0,False),Pz=prod(0,1,False))
        self.blocks={k:0 if k.startswith(('A','C')) else (1 if k.startswith('B') else 2) for k in self.D}
        if transforms is not None:
            self.D={k:v@transforms[self.blocks[k]] for k,v in self.D.items()}
    def values(self,x): return {k:v@x[self.blocks[k]] for k,v in self.D.items()}
    def pull(self,g):
        out=np.zeros((3,108))
        for k,v in g.items():out[self.blocks[k]]+=self.D[k].T@v
        return out

def transformations():
    s,z,w=quad(96);M=Matrix(s,z)
    blocks=[np.r_[np.sqrt(w*s)[:,None]*M.D['A'],np.sqrt(w)[:,None]*M.D['C']],
        np.sqrt(w*s)[:,None]*M.D['B'],
        np.r_[2*np.sqrt(w*s)[:,None]*M.D['Ps'],np.sqrt(w)[:,None]*M.D['Pz']]]
    T=[]
    for b in blocks:
        _,R=qr(b,mode='economic');T.append(solve_triangular(R,np.eye(108)))
    return np.array(T)

class Objective:
    def __init__(self,omega_target=4.,core_weight=1.,nquad=40):
        self.parent=ScaleField();self.T=transformations()
        a=np.load(ROOT/'baseline/data/profile.npz');self.x0=np.array([np.linalg.solve(t,a[n].ravel()) for t,n in zip(self.T,['F','G','P'])])
        self.x0[:2]*=np.sqrt(2)/np.linalg.norm(self.x0[:2]);self.x0=self.x0.ravel()
        s,z,w=quad(nquad);self.M=Matrix(s,z,self.T);self.w=w/64
        xx,ww=leggauss(7);self.ks=3*(xx+1);self.kw=ww/2
        self.q=2**(-self.ks)[None,:];self.h=self.parent.h;self.nu=self.parent.nu;self.tau0=self.parent.tau0
        self.ql=self.q**(-self.h);self.q2=self.q**(-2*self.h)
        self.forcing=[]
        for k in self.ks:
            sc=self.parent.scales(k);pts=np.c_[np.sqrt(s)*sc['ar'],np.zeros(len(s)),z*sc['az']]
            self.forcing.append(self.parent.forcing(pts,sc['t'])*sc['q']**1.5)
        self.forcing=np.stack(self.forcing,axis=1)
        # Predetermined local matching samples, separate from held-out audit probes.
        rr,zz=np.meshgrid(np.linspace(.025,.15,7),np.linspace(-.25,.25,11),indexing='ij')
        self.core=Matrix(rr.ravel()**2,zz.ravel(),self.T)
        self.core_weight=core_weight;self.omega_target=omega_target
        self.kappa=(1+self.h)/(2*self.tau0)
        self.cw=np.full(len(rr.ravel()),1/len(rr.ravel()))
        self.count=0;self.last={}
    def physical(self,x):
        x=np.asarray(x).reshape(3,108);n2=np.sum(x[:2]**2);lam=np.sqrt(2/n2)
        y=x.copy();y[:2]*=lam
        return y,(x,n2,lam)
    def pullnorm(self,g,st):
        x,n2,lam=st;out=g.copy();dot=np.sum(g[:2]*x[:2]);out[:2]=lam*(g[:2]-dot*x[:2]/n2)
        return out.ravel()
    def fun(self,x):
        co,state=self.physical(x);v={k:a[:,None] for k,a in self.M.values(co).items()}
        S=self.M.s[:,None];Z=self.M.z[:,None];r=np.sqrt(S);q2=self.q2;ql=self.ql;nu=self.nu;ci=.5-self.h;t0=self.tau0
        A,B,C=v['A'],v['B'],v['C']
        Ra=(A+S*v['As']+ci*Z*v['Az'])/t0+A*A+2*S*A*v['As']+C*v['Az']+q2*(-B*B+2*v['Ps'])-nu*(8*v['As']+4*S*v['Ass']+v['Azz']/q2)
        Rb=((1+self.h)*B+S*v['Bs']+ci*Z*v['Bz'])/t0+2*A*B+2*S*A*v['Bs']+C*v['Bz']-nu*(8*v['Bs']+4*S*v['Bss']+v['Bzz']/q2)
        Rc=((.5+self.h)*C+S*v['Cs']+ci*Z*v['Cz'])/t0+2*S*A*v['Cs']+C*v['Cz']+v['Pz']-nu*(4*v['Cs']+4*S*v['Css']+v['Czz']/q2)
        R=np.stack((r*Ra,r*ql*Rb,ql*Rc),axis=-1)-self.forcing
        self.last_R=R.copy()
        W=self.w[:,None,None]*self.kw[None,:,None]
        loss=.5*np.sum(W*R*R);WR=W*R;ga=WR[:,:,0]*r;gb=WR[:,:,1]*r*ql;gc=WR[:,:,2]*ql
        adj=dict(A=ga*(1/t0+2*A+2*S*v['As'])+gb*(2*B+2*S*v['Bs'])+2*S*gc*v['Cs'],
          As=ga*(S/t0+2*S*A-8*nu),Az=ga*(ci*Z/t0+C),Ass=-4*nu*S*ga,Azz=-nu*ga/q2,
          C=ga*v['Az']+gb*v['Bz']+gc*((.5+self.h)/t0+v['Cz']),
          Cs=gc*(S/t0+2*S*A-4*nu),Cz=gc*(ci*Z/t0+C),Css=-4*nu*S*gc,Czz=-nu*gc/q2,
          B=-2*q2*ga*B+gb*((1+self.h)/t0+2*A),Bs=gb*(S/t0+2*S*A-8*nu),Bz=gb*(ci*Z/t0+C),Bss=-4*nu*S*gb,Bzz=-nu*gb/q2,
          Ps=2*q2*ga,Pz=gc)
        grad=self.M.pull({key:val.sum(axis=1) for key,val in adj.items()})
        cv=self.core.values(co);targets={'A':-self.kappa+np.zeros_like(cv['A']),
          'B':self.omega_target+np.zeros_like(cv['B']),
          'C':2*self.kappa*self.core.z}
        scales={'A':1.,'B':max(1,self.omega_target),'C':.5}
        ga2={};coreloss=0.
        for name,targ in targets.items():
            err=(cv[name]-targ)/scales[name];w=self.core_weight*self.cw
            coreloss+=.5*np.sum(w*err*err);ga2[name]=w*err/scales[name]
        grad+=self.core.pull(ga2);total=loss+coreloss
        self.count+=1;self.last=dict(total_loss=float(total),conditioned_residual_loss=float(loss),core_match_loss=float(coreloss),contraction_ratio_min=float((-cv['A']).min()),contraction_ratio_mean=float((-cv['A']).mean()),axis_rotation_mean=float(cv['B'].mean()),core_swirl_min=float(cv['B'].min()),evaluations=self.count)
        return total,self.pullnorm(grad,state)
    def save(self,x,out):
        co,_=self.physical(x);m=[(t@v).reshape(9,12) for t,v in zip(self.T,co)]
        out=Path(out);out.mkdir(parents=True,exist_ok=True)
        tmp=out/'profile.tmp.npz';np.savez_compressed(tmp,F=m[0],G=m[1],P=m[2]);tmp.replace(out/'profile.npz')
        atomic(out/'stats.json',self.last);np.save(out/'variables.npy',x)
        return m

def run(name,weight=1.,steps=800):
    out=ROOT/'evidence'/name
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    atomic(out/'registration.json',dict(name=name,date=datetime.now(timezone.utc).isoformat(),core_weight=weight,omega_target=4.,maxiter=steps,conditioning='residual multiplied by q^1.5; physical residuals separately audited',velocity_initial_energy=1,core_target='A=-(1+h)/(2 tau0), B=4, C=-2Az, at declared core samples',original_physical_gates='unchanged .001; same ST064 horizon/force',native_scale_exponents='fixed .5 and .5-h, h=.005; not dynamically identified',new_fresh_log_seed=9226591,new_point_seed=9226592))
    obj=Objective(core_weight=weight);x=obj.x0.copy()
    # Directional derivative calibration before any optimizer run.
    rng=np.random.default_rng(9226500);d=rng.normal(size=x.size);d/=np.linalg.norm(d);eps=2e-6
    v,g=obj.fun(x);fd=(obj.fun(x+eps*d)[0]-obj.fun(x-eps*d)[0])/(2*eps)
    atomic(out/'gradient_check.json',dict(analytic=float(g@d),finite_difference=float(fd),absolute_error=float(abs(fd-g@d))))
    assert abs(fd-g@d)<2e-6*max(1.,abs(fd))
    history=[];count=[0];start=time.monotonic();obj.fun(x);obj.save(x,out/'step-000000')
    def callback(x):
        count[0]+=1
        if count[0]%10==0:
            obj.fun(x);record={'iteration':count[0],'elapsed':time.monotonic()-start,**obj.last};history.append(record)
            obj.save(x,out/f'step-{count[0]:06d}');atomic(out/'history.json',history)
            print(name,count[0],record,flush=True)
    result=minimize(obj.fun,x,jac=True,method='L-BFGS-B',callback=callback,options=dict(maxiter=steps,ftol=1e-13,gtol=1e-8,maxls=30,maxcor=30))
    obj.fun(result.x);obj.save(result.x,ROOT/'candidates'/name)
    atomic(out/'summary.json',dict(success=bool(result.success),message=str(result.message),iterations=int(result.nit),elapsed=time.monotonic()-start,**obj.last))
    print('END',name,json.dumps(json.loads((out/'summary.json').read_text())),flush=True)
if __name__=='__main__':
    import argparse;p=argparse.ArgumentParser();p.add_argument('--name',required=True);p.add_argument('--weight',type=float,default=1);p.add_argument('--steps',type=int,default=800);a=p.parse_args();run(a.name,a.weight,a.steps)
