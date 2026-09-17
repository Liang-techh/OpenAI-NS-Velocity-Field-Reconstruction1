"""ST003 bounded Gauss-Newton experiment, not a PDE acceptance shortcut."""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import numpy as np
from numpy.polynomial import Chebyshev
from scipy.optimize import least_squares
from spacetime import Family, force, quad, bump_derivatives, NU

class Objective:
    def __init__(self, family, seed=9172625, count=6144):
        import torch
        torch.set_num_threads(1);torch.set_default_dtype(torch.float64)
        self.torch=torch;self.f=family;self.history=[];self.start=time.time()
        f=family;n=count;rng=np.random.default_rng(seed)
        s=rng.uniform(0,4,n);z=rng.uniform(-2,2,n);t=rng.uniform(.25,.75,n)
        s[:n//6]=rng.uniform(0,.5,n//6)**2;z[:n//6]=rng.uniform(-.5,.5,n//6)
        t[n//3:n//3+n//10]=.25;t[n//2:n//2+n//10]=.75
        self.D=f.bundle(s,z,t);self.s=s;self.sqrt_s=np.sqrt(s);self.N=n
        points=np.column_stack((np.sqrt(s),np.zeros(n),z))
        self.FA=force(points,t,1,0);self.FC=force(points,t,0,1)
        tc=np.linspace(.25,.75,11);tau=1-tc;sc=.01*tau;zc=.1*tau**.495
        core=f.bundle(sc,zc,tc)
        self.core={k:torch.tensor(v) for k,v in core.items() if k in ('A','B','C','s')}
        self.CT=torch.tensor(np.column_stack([Chebyshev.basis(k)(4*tc-2) for k in range(f.nt)]))
        self.CTd=torch.tensor(np.column_stack([4*Chebyshev.basis(k).deriv()(4*tc-2) for k in range(f.nt)]))
        self.scales=torch.tensor(np.column_stack((np.sqrt(tau),tau**.505,tau**.505)))
        self.q0=torch.tensor(f.q0);self.GZ=torch.tensor(f.Gz);self.JB=torch.tensor(f.Jb)
        sq,zq,wq=quad(96)
        Cforce=float((wq*sq)@(bump_derivatives(sq/4)[0]*bump_derivatives(zq*zq/4)[0]))
        self.torque_t=torch.tensor(Cforce*bump_derivatives((2*tc-1)**2)[0])
        seedvel=f.fields(f.initial(),np.column_stack((np.sqrt(sc),np.zeros(11),zc)),tc)[0]*self.scales.numpy()
        self.floor=.5*np.linalg.norm(seedvel[0]);self.cache_x=None
    def extra(self,v):
        torch=self.torch;f=self.f;n=f.n;nt=f.nt;core=self.core
        al=v[:n].reshape(f.ns,nt);be=v[n:2*n].reshape(f.ns,nt)
        scale=torch.sqrt(2/(torch.sum((al@self.q0)**2)+torch.sum((be@self.q0)**2)))
        aa=scale*al.flatten();bb=scale*be.flatten()
        ap=scale*al@self.CT.T;bp=scale*be@self.CT.T
        E=.5*(torch.sum(ap*ap,dim=0)+torch.sum(bp*bp,dim=0))
        cv=torch.stack((torch.sqrt(core['s'])*(core['A']@aa),torch.sqrt(core['s'])*(core['B']@bb),core['C']@aa),dim=1)*self.scales
        ref=cv[0];norm=torch.linalg.norm(ref);drift=torch.linalg.norm(cv-ref,dim=1)/(norm+1e-12)
        D=1.5*torch.sum(ap*(self.GZ@ap),dim=0)-E
        torque=self.JB@(scale*be@self.CTd.T)+v[-1]*self.torque_t
        m=len(E)
        return torch.cat((10*torch.relu(.1-E)/np.sqrt(m),10*torch.relu(E-10)/np.sqrt(m),
            np.sqrt(10/m)*torch.relu(drift-.045),
            (np.sqrt(10)*torch.relu(self.floor-norm)/self.floor).reshape(1),
            np.sqrt(10/m)*torch.relu(cv[:,0]+1e-6),np.sqrt(10/m)*torch.relu(-cv[:,1]+1e-6),np.sqrt(10/m)*torch.relu(-cv[:,2]+1e-6),
            10*D/np.sqrt(m*88*np.pi/3),np.sqrt(10)*torque/np.sqrt(m*32*np.pi)))
    def evaluate(self,x):
        if self.cache_x is not None and np.array_equal(x,self.cache_x):return
        f=self.f;D=self.D;n=f.n;N=self.N;S=self.s;r=self.sqrt_s
        aa,bb,qq,fc,amp=f.coefficients(x)
        vals={k:D[k]@(bb if k.startswith('B') else qq if k.startswith('Q') else aa) for k in D if k not in ('s','z','t','Q')}
        A,B,C=vals['A'],vals['B'],vals['C'];As,Az=vals['As'],vals['Az'];Bs,Bz=vals['Bs'],vals['Bz'];Cs,Cz=vals['Cs'],vals['Cz']
        Rr=r*(vals['At']+A*A+2*S*A*As+C*Az-B*B+2*vals['Qs']-NU*vals['AL'])-fc[0]*self.FA[:,0]-fc[1]*self.FC[:,0]
        Rt=r*(vals['Bt']+2*A*B+2*S*A*Bs+C*Bz-NU*vals['BL'])-fc[0]*self.FA[:,1]-fc[1]*self.FC[:,1]
        Rz=vals['Ct']+2*S*A*Cs+C*Cz+vals['Qz']-NU*vals['CL']-fc[0]*self.FA[:,2]-fc[1]*self.FC[:,2]
        J=np.zeros((3*N,3*n+2))
        J[:N,:n]=r[:,None]*(D['At']+(2*A+2*S*As)[:,None]*D['A']+(2*S*A)[:,None]*D['As']+Az[:,None]*D['C']+C[:,None]*D['Az']-NU*D['AL'])
        J[:N,n:2*n]=(-2*r*B)[:,None]*D['B']
        J[N:2*N,:n]=r[:,None]*((2*B+2*S*Bs)[:,None]*D['A']+Bz[:,None]*D['C'])
        J[N:2*N,n:2*n]=r[:,None]*(D['Bt']+(2*A)[:,None]*D['B']+(2*S*A)[:,None]*D['Bs']+C[:,None]*D['Bz']-NU*D['BL'])
        J[2*N:,:n]=D['Ct']+(2*S*Cs)[:,None]*D['A']+(2*S*A)[:,None]*D['Cs']+Cz[:,None]*D['C']+C[:,None]*D['Cz']-NU*D['CL']
        J[:N,2*n:3*n]=2*r[:,None]*D['Qs'];J[2*N:,2*n:3*n]=D['Qz']
        J[:,-2]=-self.FA.T.ravel();J[:,-1]=-self.FC.T.ravel()
        # Pull back through exact nonzero initial-energy normalization.
        q=x[:2*n];qmat=q.reshape(2*f.ns,f.nt);initial=qmat@f.q0
        v=(initial[:,None]*f.q0[None,:]).ravel();den=initial@initial
        directional=J[:,:2*n]@q
        J[:,:2*n]=amp*(J[:,:2*n]-directional[:,None]*(v/den)[None,:])
        R=np.concatenate((Rr,Rt,Rz))/np.sqrt(N);J/=np.sqrt(N)
        torch=self.torch;tx=torch.tensor(x,requires_grad=True)
        er=self.extra(tx).detach().numpy()
        ej=torch.autograd.functional.jacobian(self.extra,tx,vectorize=True).detach().numpy()
        self.cache_x=x.copy();self.cache_r=np.r_[R,er];self.cache_j=np.vstack((J,ej))
        row={'evaluation':len(self.history)+1,'loss':float(self.cache_r@self.cache_r),'pde_mse':float(R@R),'extra_loss':float(er@er),'elapsed':time.time()-self.start}
        self.history.append(row)
        if len(self.history)%5==0:print(json.dumps(row),flush=True)
    def fun(self,x):self.evaluate(x);return self.cache_r
    def jac(self,x):self.evaluate(x);return self.cache_j


def run(warm,out,maxfun=100):
    f,x=Family.load(warm);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    obj=Objective(f);f.save(x,out/'initial.json')
    low=np.r_[np.full(2*f.n,-4.),np.full(f.n,-100.),0.,0.];high=np.r_[np.full(2*f.n,4.),np.full(f.n,100.),10.,10.]
    result=least_squares(obj.fun,x,jac=obj.jac,bounds=(low,high),method='trf',tr_solver='lsmr',x_scale='jac',max_nfev=maxfun,ftol=1e-12,xtol=1e-12,gtol=1e-10,tr_options={'maxiter':80,'atol':1e-5,'btol':1e-5},verbose=0)
    info={'experiment_id':'CR-ROOT-ST003','optimizer_success':bool(result.success),'message':result.message,'nfev':result.nfev,'njev':result.njev,'loss':float(2*result.cost),'optimality':float(result.optimality),'elapsed':time.time()-obj.start,'training_seed':9172625,'training_points':6144,'parameters':len(x),'pde_validated':False}
    f.save(result.x,out/'candidate.json',info);(out/'training_history.json').write_text(json.dumps(obj.history,indent=2)+'\n');(out/'training_summary.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--warm',default='moment_7x7x8/candidate.json');p.add_argument('--out',default='gauss_7x7x8');p.add_argument('--maxfun',type=int,default=100);a=p.parse_args();run(a.warm,a.out,a.maxfun)
