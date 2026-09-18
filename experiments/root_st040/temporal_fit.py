"""ST041: constrained temporal coefficient repair of the unchanged spatial family."""
import bootstrap
import argparse,json,hashlib,time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
from spacetime import Family
from controlled_fit import Model,embed

class TemporalModel(Model):
 def __init__(self,path):
  self.f,self.raw=Family.load(path);f=self.f
  self.a,self.b,self.p,self.fc,_=f.coefficients(self.raw);cols=[]
  for co in (self.a,self.b,self.p):
   tmp=co.reshape(f.ns,f.nt);M=np.zeros((f.ns,f.nt,f.nt))
   for k in range(f.nt):M[:,k,k]=tmp[:,k]
   cols.append(M.reshape(f.n,f.nt))
  self.Ma,self.Mb,self.Mp=cols;self.na=self.nb=self.np=f.nt;self.nv=2*f.nt;self.dim=3*f.nt+2
  aa=self.a.reshape(f.ns,f.nt)@f.q0;bb=self.b.reshape(f.ns,f.nt)@f.q0;self.e0=np.r_[aa,bb];self.EM=np.zeros((2*f.ns,self.nv));self.EM[:f.ns,:self.na]=np.einsum('stk,t->sk',self.Ma.reshape(f.ns,f.nt,-1),f.q0);self.EM[f.ns:,self.na:]=np.einsum('stk,t->sk',self.Mb.reshape(f.ns,f.nt,-1),f.q0)

class TemporalBiasModel(TemporalModel):
 def __init__(self,path):
  self.f,self.raw=embed(path);f=self.f
  self.a,self.b,self.p,self.fc,_=f.coefficients(self.raw);cols=[]
  for block,co,T in [('A',self.a,f.Tp),('B',self.b,f.Tw),('Q',self.p,f.Tq)]:
   tmp=co.reshape(f.ns,f.nt);M=np.zeros((f.ns,f.nt,f.nt))
   for k in range(f.nt):M[:,k,k]=tmp[:,k]
   columns=[M.reshape(f.n,f.nt)]
   if block!='B':
    extra=[]
    for ir in range(2):
     v=np.zeros((f.nr,f.nz));v[ir,9]=1;vec=np.linalg.solve(T,v.ravel());vec/=np.linalg.norm(vec)
     for k in range(4):
      m=np.zeros((f.ns,f.nt));m[:,k]=vec;extra.append(m.ravel())
    columns.append(np.column_stack(extra))
   cols.append(np.column_stack(columns))
  self.Ma,self.Mb,self.Mp=cols;self.na,self.nb,self.np=[c.shape[1] for c in cols];self.nv=self.na+self.nb;self.dim=self.nv+self.np+2
  aa=self.a.reshape(f.ns,f.nt)@f.q0;bb=self.b.reshape(f.ns,f.nt)@f.q0;self.e0=np.r_[aa,bb];self.EM=np.zeros((2*f.ns,self.nv));self.EM[:f.ns,:self.na]=np.einsum('stk,t->sk',self.Ma.reshape(f.ns,f.nt,-1),f.q0);self.EM[f.ns:,self.na:]=np.einsum('stk,t->sk',self.Mb.reshape(f.ns,f.nt,-1),f.q0)

def run(warm,out,maxiter=300,bias=False):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 reg=dict(id='ST042-time-bias' if bias else 'ST041-time-only',created_before_fit=True,parent_sha256=hashlib.sha256(Path(warm).read_bytes()).hexdigest(),training='Gauss s32,z48,time11 plus endpoints; time-slice constraints at 17uniform slices',validation_seed=9174101,maxiter=maxiter,coefficient_relative_bounds=[-.5,.5],spatial_family='existing opposite-parity extension, 8 poloidal and 8 pressure repair coefficients' if bias else 'UNCHANGED ST006',physical_gates='UNCHANGED',objective='time-averaged full momentum squared plus weak parameter regularization; source-core expanded-profile no worse than parent at training times; original single-probe drift <=.049',scope='source-inspired bias/shear magnitudes are autonomous; no pulses or identity proof',bias_constraints=bool(bias),bias_bounds=[.001,.003] if bias else None)
 (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n');m=(TemporalBiasModel if bias else TemporalModel)(warm)
 x,wx=leggauss(32);z,wz=leggauss(48);tt,wt=leggauss(11);t=np.r_[.25,.5+.25*tt,.75];wt=np.r_[.05,.45*wt,.05];S,Z,T=np.meshgrid(2*(x+1),2*z,t,indexing='ij');ww=(np.outer(wx,wz)[:,:,None]/4*wt).ravel();D=m.cache(S,Z,T)
 tc=np.linspace(.25,.75,17);tau=1-tc;core=m.cache(.01*tau,.1*tau**.495,tc);scales=np.stack((np.sqrt(tau),tau**.505,tau**.505),axis=1)
 # Additional, autonomous broad-core preservation prevents improving just one probe.
 R,Z,T=np.meshgrid([.05,.1,.15,.2],[-.2,-.1,.1,.2],tc,indexing='ij');ts=1-T;expanded=m.cache(R*R*ts,Z*ts**.495,T);es=np.stack((np.sqrt(ts),ts**.505,ts**.505),axis=-1).reshape(-1,3);shape=R.shape
 midR,midT=np.meshgrid([.05,.1,.2,.4],tc,indexing='ij');mid=m.cache(midR**2*(1-midT),np.zeros_like(midR),midT);ms=(1-midT.ravel())**.505
 hist=[];last=None;cache=None;start=time.time()
 def velocity_drift(c,K,scale,shape=None):
  u,J=m.velocity(c,K);u=u*scale;J=J*scale[:,:,None]
  if shape is None:
   ref=u[0];JR=J[0];delta=u-ref;JD=J-JR;den=np.linalg.norm(ref);num=np.sqrt(np.sum(delta*delta,axis=1)+1e-28);dn=np.einsum('ni,nik->nk',delta,JD)/num[:,None];dr=(dn*den-num[:,None]*(ref@JR)/den)/den**2
   return num/den,dr,u,J
  u=u.reshape(shape+(3,));J=J.reshape(shape+(3,m.nv));delta=u-u[:,:,0:1];JD=J-J[:,:,0:1]
  den=np.linalg.norm(u[:,:,0]);num=np.sqrt(np.sum(delta*delta,axis=(0,1,3))+1e-28);dnum=np.einsum('abtc,abtck->tk',delta,JD)/num[:,None];dden=np.einsum('abc,abck->k',u[:,:,0],J[:,:,0])/den
  return num/den,(dnum*den-num[:,None]*dden)/den**2,u,J
 parent_exp=velocity_drift(np.zeros(m.dim),expanded,es,shape)[0]
 def constraints(c):
  drift,dj,u,UJ=velocity_drift(c,core,scales);ex,exj,_,_=velocity_drift(c,expanded,es,shape)
  sign=np.array([-1,1,1]);sgn=u*sign
  vals=np.r_[.049-drift,sgn.ravel()-1e-6,(parent_exp+.005)-ex]
  J=np.vstack((-dj,(UJ*sign[None,:,None]).reshape(-1,m.nv),-exj))
  if bias:
   mu,mj=m.velocity(c,mid);v=(mu[:,2]*ms).reshape(4,-1);vj=(mj[:,2]*ms[:,None]).reshape(4,-1,m.nv)
   vals=np.r_[vals,(v-.001).ravel(),(.003-v).ravel(),v[0]-v[-1]-.0001]
   J=np.vstack((J,vj.reshape(-1,m.nv),-vj.reshape(-1,m.nv),vj[0]-vj[-1]))
  return vals,np.pad(J,((0,0),(0,m.dim-m.nv)))
 def fun(c):
  nonlocal last,cache
  if last is None or not np.array_equal(c,last):
   R,J=m.momentum(c,D);v=float(np.sum(ww[:,None]*R*R));g=2*np.einsum('n,ni,nik->k',ww,R,J);regu=1e-7*np.sum(c*c);v+=regu;g+=2e-7*c;cache=v,g;last=c.copy();hist.append(dict(eval=len(hist)+1,pde_objective=v,elapsed=time.time()-start))
  return cache
 bounds=[(-.5,.5)]*(m.dim-2)+[(-float(v),10-float(v)) for v in m.fc]
 if bias:
  bounds[m.f.nt:m.na]=[(-.005,.005)]*(m.na-m.f.nt)
  bounds[m.nv+m.f.nt:m.nv+m.np]=[(-.05,.05)]*(m.np-m.f.nt)
 c=np.zeros(m.dim);r0=fun(c)[0];ret=minimize(fun,c,jac=True,bounds=bounds,constraints=[{'type':'ineq','fun':lambda c:constraints(c)[0],'jac':lambda c:constraints(c)[1]}],method='SLSQP',options=dict(maxiter=maxiter,ftol=1e-13,disp=True));c=ret.x;raw=m.candidate(c)
 # Verify constraints and bounds before serializing, regardless of solver flag.
 values,jj=constraints(c);summary=dict(**reg,variables=m.dim,optimizer_success=bool(ret.success),message=ret.message,iterations=ret.nit,evaluations=ret.nfev,initial_training_loss=r0,final_training_loss=float(ret.fun),minimum_training_constraint=float(np.min(values)),scaled_core_drift_max=float(np.max(velocity_drift(c,core,scales)[0])),relative_modifiers=c.tolist(),elapsed=time.time()-start,pde_validated=False)
 m.f.save(raw,out/'candidate.json',summary);summary['candidate_sha256']=hashlib.sha256((out/'candidate.json').read_bytes()).hexdigest();(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'history.json').write_text(json.dumps(hist,indent=2)+'\n');np.savez_compressed(out/'cache_parameters.npz',c=c);print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--warm',required=True);p.add_argument('--out',required=True);p.add_argument('--maxiter',type=int,default=300);p.add_argument('--bias',action='store_true');a=p.parse_args();run(a.warm,a.out,a.maxiter,a.bias)
