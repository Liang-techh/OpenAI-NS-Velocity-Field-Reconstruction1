"""Separately registered preservation control: same-time parent volumes and bias.
Introduced after the first three frozen trials exposed concentration and pressure losses.
New held-out sample IDs must be used; this is not a retune advertised on old holdouts.
"""
from fit_profile import *

class Guarded(Objective):
    def __init__(self):
        super().__init__(omega_target=4.,core_weight=1.,nquad=40)
        co,_=self.physical(self.x0);self.v0=self.participation(self.M.values(co))[0]
        self.mid=Matrix(np.linspace(.025,.15,13)**2,np.zeros(13),self.T)
        self.bias_target=self.mid.values(co)['C']
    def participation(self,v):
        s=self.M.s[:,None];p=2**(2*self.h*np.array([0.,3.,6.]))[None,:]
        A,B,C=[v[k][:,None] for k in ['A','B','C']]
        W=v['Az'][:,None]/np.sqrt(p)-2*np.sqrt(p)*v['Cs'][:,None]
        V=v['B'][:,None]+s*v['Bs'][:,None]
        densities=[s*(A*A+p*B*B)+p*C*C,s*v['Bz'][:,None]**2+s*W*W+4*p*V*V]
        weights=64*self.w[:,None];vol=[];coeffs=[]
        for rho in densities:
            a=np.sum(weights*rho,axis=0);b=np.sum(weights*rho*rho,axis=0)
            vol.append(a*a/b);coeffs.append(weights*(2/a-2*rho/b))
        return np.array(vol),(densities,coeffs,W,V,p)
    def fun(self,x):
        loss,grad=super().fun(x);co,st=self.physical(x);v=self.M.values(co);vol,details=self.participation(v);ratio=vol/self.v0
        # Narrow interior margins are autonomous and registered in advance.
        lo,hi=.9905,1.0195
        err=np.where(ratio<lo,np.log(ratio/lo),np.where(ratio>hi,np.log(ratio/hi),0.))
        vg=1500*err/6
        loss+=.5*1500*np.mean(err*err)
        rho,fac,W,V,p=details;s=self.M.s[:,None]
        e=fac[0]*vg[0][None,:];o=fac[1]*vg[1][None,:]
        ad=dict(A=2*s*v['A'][:,None]*e,B=2*s*p*v['B'][:,None]*e+8*p*V*o,C=2*p*v['C'][:,None]*e,
            Bz=2*s*v['Bz'][:,None]*o,Az=2*s*W/np.sqrt(p)*o,Cs=-4*s*W*np.sqrt(p)*o,Bs=8*s*p*V*o)
        adj=self.M.pull({k:a.sum(axis=1) for k,a in ad.items()})
        cv=self.core.values(co);sign=np.sign(self.core.z);active=abs(self.core.z)>.03
        diff=np.maximum(.001-sign*cv['Pz'],0)*active;scale=.2
        loss+=.5*10*np.mean((diff/scale)**2)
        pg=-10*diff*sign/(scale**2*len(diff));adj+=self.core.pull({'Pz':pg})
        bias=self.mid.values(co)['C'];e=bias-self.bias_target
        loss+=.5*50*np.mean((e/.02)**2)
        adj+=self.mid.pull({'C':50*e/(.02**2*len(e))})
        grad+=self.pullnorm(adj,st)
        self.last.update(total_loss=float(loss),volume_ratios=ratio.tolist(),core_pressure_correct_fraction=float(np.mean(sign[active]*cv['Pz'][active]>0)),min_bias=float(bias.min()),max_bias_error=float(abs(e).max()))
        return loss,grad

def run():
    out=ROOT/'evidence/ST065-V'
    if out.exists():raise FileExistsError(out)
    out.mkdir()
    atomic(out/'registration.json',dict(utc=datetime.now(timezone.utc).isoformat(),purpose='Post-diagnostic preservation control, separate from frozen R/E/H trials',start='ST064-S',core_targets='same radial/rotation targets as E',new_guards='finite penalty screens for same-time energy/enstrophy effective volumes and midplane bias, axial pressure direction',volume_ratio_goal=[.99,1.02],volume_penalty=1500,volume_interior=[.9905,1.0195],pressure_penalty=10,bias_penalty=50,steps=2400,fresh_log_seed=9226595,fresh_point_seed=9226596,physical_momentum_gates=.001,parameterization='unchanged static 324 coefficients and imposed scale laws'))
    o=Guarded();x=o.x0.copy();rng=np.random.default_rng(9226505);d=rng.normal(size=324);d/=np.linalg.norm(d);v,g=o.fun(x);h=1e-6
    fd=(o.fun(x+h*d)[0]-o.fun(x-h*d)[0])/(2*h)
    atomic(out/'gradient_check.json',dict(analytic=float(g@d),fd=float(fd),error=float(abs(fd-g@d))));assert abs(fd-g@d)<2e-5*max(1,abs(fd))
    n=[0];start=time.monotonic();history=[];o.fun(x);o.save(x,out/'step-000000')
    def cb(x):
        n[0]+=1
        if n[0]%50==0:
            o.fun(x);record={'iteration':n[0],'elapsed':time.monotonic()-start,**o.last};history.append(record);o.save(x,out/f'step-{n[0]:06d}');atomic(out/'history.json',history)
            print(record,flush=True)
    r=minimize(o.fun,x,jac=True,method='L-BFGS-B',callback=cb,options=dict(maxiter=2400,maxcor=40,ftol=1e-13,gtol=1e-8,maxls=40))
    o.fun(r.x);o.save(r.x,ROOT/'candidates/ST065-V')
    atomic(out/'summary.json',dict(success=bool(r.success),message=str(r.message),iterations=int(r.nit),elapsed=time.monotonic()-start,**o.last))
    atomic(out/'FROZEN.json',dict(utc=datetime.now(timezone.utc).isoformat(),profile_sha256=hashlib.sha256((ROOT/'candidates/ST065-V/profile.npz').read_bytes()).hexdigest(),scope='New preservation-control endpoint; no subsequent coefficient changes'))
if __name__=='__main__':run()
