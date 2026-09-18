"""ST043: bounded joint corrections with genuine expanded-core constraints.
Uses the original physical field and restricted force, not a manufactured source.
Source-inspired diagnostic domains and bounds are autonomous finite-window choices.
"""
from __future__ import annotations
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import Chebyshev
from scipy.optimize import minimize
from scipy.special import logsumexp
import bootstrap
from spacetime import Family
from controlled_fit import Model

class CoupledModel(Model):
 def __init__(self,path,radial=3,axial=3,time_degree=4):
  self.f,self.raw=Family.load(path);f=self.f
  if f.basis_kind!='hybrid9_axial_opposite3_v1':raise ValueError('Expected frozen asymmetric ST042 family')
  self.a,self.b,self.p,self.fc,_=f.coefficients(self.raw);cols=[];bounds=[];self.descriptors=[]
  for block,co,T in [('A',self.a,f.Tp),('B',self.b,f.Tw),('Q',self.p,f.Tq)]:
   tmp=co.reshape(f.ns,f.nt);M=np.zeros((f.ns,f.nt,f.nt))
   for k in range(f.nt):M[:,k,k]=tmp[:,k]
   columns=[M.reshape(f.n,f.nt)];bb=[(-.5,.5)]*f.nt
   self.descriptors.extend([[block,'parent_time',k] for k in range(f.nt)])
   extra=[]
   for ir in range(radial):
    for iz in range(axial):
     poly=np.zeros((f.nr,f.nz));poly[ir,iz]=1.;vec=np.linalg.solve(T,poly.ravel());vec/=np.linalg.norm(vec)
     for k in range(time_degree):
      m=np.zeros((f.ns,f.nt));m[:,k]=vec;extra.append(m.ravel());bb.append((-.05,.05) if block!='Q' else (-.15,.15));self.descriptors.append([block,'lowmode',ir,iz,k])
   for ir in range(2):
    poly=np.zeros((f.nr,f.nz));poly[ir,9]=1.;vec=np.linalg.solve(T,poly.ravel());vec/=np.linalg.norm(vec)
    for k in range(time_degree):
     m=np.zeros((f.ns,f.nt));m[:,k]=vec;extra.append(m.ravel());bb.append((-.005,.005) if block!='Q' else (-.05,.05));self.descriptors.append([block,'opposite',ir,9,k])
   columns.append(np.column_stack(extra));cols.append(np.column_stack(columns));bounds+=bb
  self.Ma,self.Mb,self.Mp=cols;self.na,self.nb,self.np=[c.shape[1] for c in cols];self.nv=self.na+self.nb;self.dim=self.nv+self.np+2
  self.bounds=bounds+[(-float(v),10-float(v)) for v in self.fc]
  aa=self.a.reshape(f.ns,f.nt)@f.q0;bb=self.b.reshape(f.ns,f.nt)@f.q0;self.e0=np.r_[aa,bb]
  self.EM=np.zeros((2*f.ns,self.nv));self.EM[:f.ns,:self.na]=np.einsum('stk,t->sk',self.Ma.reshape(f.ns,f.nt,-1),f.q0);self.EM[f.ns:,self.na:]=np.einsum('stk,t->sk',self.Mb.reshape(f.ns,f.nt,-1),f.q0)

class JointObjective:
 def __init__(self,m,profile_ratio=.75,pressure_ratio=.9,pressure_target=None,space_order=(28,40),peak_weight=.05,core_weight=0.,anchor_tolerance=.1,time_cap=None,softmax_weight=0.):
  self.m=m;self.profile_ratio=profile_ratio;self.pressure_ratio=pressure_ratio;self.pressure_target=pressure_target;self.peak_weight=peak_weight;self.core_weight=core_weight;self.anchor_tolerance=anchor_tolerance;self.time_cap=time_cap;self.softmax_weight=softmax_weight;self.moment_cache=None
  x,wx=leggauss(space_order[0]);z,wz=leggauss(space_order[1]);tt,wt=leggauss(13)
  times=np.r_[.25,.5+.25*tt,.75];wt=np.r_[.075,.425*wt,.075];S,Z,T=np.meshgrid(2*(x+1),2*z,times,indexing='ij')
  self.weights=(np.outer(wx,wz)[:,:,None]/4*wt).ravel();self.D=m.cache(S,Z,T);self.training_shape=S.shape;self.times=times;self.spatial_weights=np.outer(wx,wz)/4
  init_R,_=m.momentum(np.zeros(m.dim),self.D);self.initial_time_mse=np.einsum("ij,ijt->t",self.spatial_weights,np.sum(init_R**2,axis=1).reshape(S.shape))
  tc=np.linspace(.25,.75,17);self.tc=tc;tau=1-tc
  self.core=m.cache(.01*tau,.1*tau**.495,tc);self.scales=np.stack((tau**.5,tau**.505,tau**.505),axis=1)
  R,Z,T=np.meshgrid([.04,.08,.12,.16,.2],[-.2,-.16,-.12,-.08,-.04,.04,.08,.12,.16,.2],tc,indexing='ij');tau=1-T
  self.shape=R.shape;self.exp=m.cache(R*R*tau,Z*tau**.495,T);self.es=np.stack((tau**.5,tau**.505,tau**.505),axis=-1).reshape(-1,3);self.signz=np.sign(Z.ravel());self.z=Z.ravel()
  R,T=np.meshgrid([.05,.1,.2,.4],tc,indexing='ij');self.mid=m.cache(R*R*(1-T),np.zeros_like(R),T);self.ms=(1-T.ravel())**.505
  zero=np.zeros(m.dim);self.coreparent=m.velocity(zero,self.core)[0]*self.scales
  U=m.velocity(zero,self.exp)[0]*self.es;self.up=U.reshape(self.shape+(3,));self.expden=float(np.sum(self.up[:,:,0]**2));self.epdelta=self.up-self.up[:,:,0:1];self.parent_num=np.sum(self.epdelta**2,axis=(0,1,3))
  self.pgrad_parent=self.signz*self.exp['v']['Qz'];self.pgrad_scale=max(float(np.sqrt(np.mean(self.pgrad_parent**2))),.01)
  self.radial_scale=max(float(np.sqrt(np.mean(self.exp['v']['Qs']**2))),.01)
  self.last=None;self.cached=None;self.history=[];self.start=time.time()
 def constraints(self,c):
  m=self.m;u,J=m.velocity(c,self.core);u*=self.scales;J*=self.scales[:,:,None];ref=u[0];JR=J[0];dif=u-ref;DJ=J-JR
  den=float(ref@ref);denj=2*ref@JR;num=np.sum(dif*dif,axis=1);numj=2*np.einsum('ni,nik->nk',dif,DJ)
  vals=[(.049**2*den-num)[1:]/max(float(self.coreparent[0]@self.coreparent[0]),1e-10)]
  jacs=[(.049**2*denj[None,:]-numj)[1:]/max(float(self.coreparent[0]@self.coreparent[0]),1e-10)]
  signs=np.array([-1,1,1]);scale=max(np.linalg.norm(self.coreparent[0]),.01)
  vals.append((u*signs-1e-6).ravel()/scale);jacs.append((J*signs[None,:,None]).reshape(-1,m.nv)/scale)
  U,K=m.velocity(c,self.exp);U*=self.es;K*=self.es[:,:,None];U=U.reshape(self.shape+(3,));K=K.reshape(self.shape+(3,m.nv))
  d=U-U[:,:,0:1];dj=K-K[:,:,0:1];den=float(np.sum(U[:,:,0]**2));denj=2*np.einsum('abq,abqk->k',U[:,:,0],K[:,:,0]);num=np.sum(d*d,axis=(0,1,3));numj=2*np.einsum('abtq,abtqk->tk',d,dj)
  cap2=np.maximum(self.profile_ratio**2*self.parent_num/self.expden,1e-8)
  vals.append(((cap2*den-num)/self.expden)[1:]);jacs.append(((cap2[:,None]*denj-numj)/self.expden)[1:])
  # Prevent artificial relative-drift improvement by inflating/collapsing the initial core.
  delta0=U[:,:,0]-self.up[:,:,0];error0=np.sum(delta0**2)/self.expden
  vals.append(np.array([self.anchor_tolerance**2-error0]));jacs.append((-2*np.einsum('abq,abqk->k',delta0,K[:,:,0])/self.expden)[None,:])
  # Preserve the enlarged-core velocity directions, not just the original single probe.
  sgn=np.stack((-np.ones(len(self.z)),np.ones(len(self.z)),self.signz),axis=1)
  Uf=U.reshape(-1,3);Kf=K.reshape(-1,3,m.nv);vals.append((Uf*sgn-1e-7).ravel()/scale);jacs.append((Kf*sgn[:,:,None]).reshape(-1,m.nv)/scale)
  mu,mj=m.velocity(c,self.mid);mv=(mu[:,2]*self.ms).reshape(4,-1);MJ=(mj[:,2]*self.ms[:,None]).reshape(4,-1,m.nv)
  vals.extend([(mv-.0008).ravel(),(.0032-mv).ravel(),mv[0]-mv[-1]-.00008]);jacs.extend([MJ.reshape(-1,m.nv),-MJ.reshape(-1,m.nv),MJ[0]-MJ[-1]])
  V=np.concatenate(vals);Jac=np.pad(np.vstack(jacs),((0,0),(0,m.dim-m.nv)))
  # Pressure signs are linear and do not use velocity-derived forcing.
  pg=self.signz*(self.exp['v']['Qz']+self.exp['M']['Qz']@c[m.nv:m.nv+m.np]);target=self.pressure_ratio*self.pgrad_parent if self.pressure_target is None else np.full(len(pg),self.pressure_target)
  pressure=(pg-target)/self.pgrad_scale;JP=np.zeros((len(pg),m.dim));JP[:,m.nv:m.nv+m.np]=self.signz[:,None]*self.exp['M']['Qz']/self.pgrad_scale
  radial=(self.exp['v']['Qs']+self.exp['M']['Qs']@c[m.nv:m.nv+m.np]-1e-7)/self.radial_scale;JRp=np.zeros((len(pg),m.dim));JRp[:,m.nv:m.nv+m.np]=self.exp['M']['Qs']/self.radial_scale
  if self.time_cap is not None:
   R,J=self.cached_momentum(c);sr=self.training_shape
   temporal=np.einsum('ij,ijt->t',self.spatial_weights,np.sum(R*R,axis=1).reshape(sr))
   TJ=2*np.einsum('ij,ijtc,ijtck->tk',self.spatial_weights,R.reshape(sr+(3,)),J.reshape(sr+(3,m.dim)))
   V=np.r_[V,(self.time_cap*self.initial_time_mse-temporal)/self.initial_time_mse]
   Jac=np.vstack((Jac,-TJ/self.initial_time_mse[:,None]))
  return np.r_[V,pressure,radial],np.vstack((Jac,JP,JRp))
 def cached_momentum(self,c):
  if self.moment_cache is None or not np.array_equal(self.moment_cache[0],c):
   R,J=self.m.momentum(c,self.D);self.moment_cache=(c.copy(),R,J)
  return self.moment_cache[1:]
 def fun(self,c):
  if self.last is None or not np.array_equal(c,self.last):
   R,J=self.cached_momentum(c);norm=np.sum(R*R,axis=1);weight=self.weights
   pde=float(weight@norm);peak=float(weight@(norm*norm)/(.1**2));loss=pde+self.peak_weight*peak+1e-7*float(c@c)
   gradient=2*np.einsum('n,ni,nik->k',weight*(1+2*self.peak_weight*norm/.1**2),R,J)+2e-7*c
   if self.softmax_weight:
    temperature=.0005;logweights=np.log(weight)+norm/temperature;normalizer=logsumexp(logweights);prob=np.exp(logweights-normalizer)
    loss+=self.softmax_weight*temperature*normalizer;gradient+=self.softmax_weight*2*np.einsum('n,ni,nik->k',prob,R,J)
   if self.core_weight:
    cr,cj=self.m.momentum(c,self.exp);loss+=self.core_weight*float(np.mean(np.sum(cr*cr,axis=1)));gradient+=self.core_weight*2*np.einsum('ni,nik->k',cr,cj)/len(cr)
   self.last=c.copy();self.cached=loss,gradient;self.history.append({'evaluation':len(self.history)+1,'loss':loss,'pde_mean_square':pde,'peak_penalty':peak,'elapsed':time.time()-self.start})
   if len(self.history)%20==0:print(json.dumps(self.history[-1]),flush=True)
  return self.cached

def run(parent,out,profile_ratio=.75,pressure_ratio=.9,pressure_target=None,maxiter=240,radial=3,axial=2,time_degree=4,peak_weight=.05,localized=False,core_weight=0.,anchor_tolerance=.1,time_cap=None,softmax_weight=0.):
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise ValueError('Output must be empty')
 reg=dict(experiment='ST044-localized' if localized else 'ST043',localized=localized,core_weight=core_weight,anchor_tolerance=anchor_tolerance,time_cap=time_cap,softmax_weight=softmax_weight,parent_sha256=hashlib.sha256(Path(parent).read_bytes()).hexdigest(),profile_ratio=profile_ratio,pressure_ratio=pressure_ratio,pressure_target=pressure_target,maxiter=maxiter,radial=radial,axial=axial,time_degree=time_degree,peak_weight=peak_weight,training='Gauss28x40 spatial,13Gauss times+endpoints; expanded grid5x10x17; one-point/midplane17times',independent_seeds=[9174301,9174302],all_original_physical_gates='unchanged',new_structure='autonomous expanded drift/pressure limits, initial core difference controlled by anchor_tolerance; original bias/signs retained',written_before_optimization=True)
 (out/'registration.json').write_text(json.dumps(reg,indent=2)+'\n');t0=time.time()
 if localized:
  from localized_model import LocalizedModel
  m=LocalizedModel(parent,radial,axial,time_degree)
  (out/'projection.json').write_text(json.dumps(m.projection,indent=2)+'\n')
 else:m=CoupledModel(parent,radial,axial,time_degree)
 print('Variables',m.dim,flush=True);obj=JointObjective(m,profile_ratio,pressure_ratio,pressure_target,peak_weight=peak_weight,core_weight=core_weight,anchor_tolerance=anchor_tolerance,time_cap=time_cap,softmax_weight=softmax_weight);c=np.zeros(m.dim)
 rng=np.random.default_rng(9174290);v=rng.normal(size=m.dim);v/=np.linalg.norm(v);h=1e-6;f,g=obj.fun(c);gp=obj.fun(c+h*v)[0];gm=obj.fun(c-h*v)[0];C,J=obj.constraints(c);Cp=obj.constraints(c+h*v)[0];Cm=obj.constraints(c-h*v)[0]
 check=dict(objective_abs_directional_error=abs((gp-gm)/(2*h)-g@v),constraint_max_directional_error=float(np.max(np.abs((Cp-Cm)/(2*h)-J@v))),initial_min_constraint=float(C.min()))
 (out/'derivative_checks.json').write_text(json.dumps(check,indent=2)+'\n');print(check,flush=True)
 if check['constraint_max_directional_error']>1e-5:raise ValueError('Constraint gradient fails')
 iterations=[0]
 def checkpoint(v):
  iterations[0]+=1
  if iterations[0]%10==0:
   try:m.f.save(m.candidate(v),out/'checkpoint_candidate.json',{'iteration':iterations[0],'scope':'unfinished training checkpoint','pde_validated':False})
   except ValueError:pass
   np.save(out/'checkpoint_modifiers.npy',v)
   (out/'history.json').write_text(json.dumps(obj.history,indent=2)+'\n')
 ret=minimize(obj.fun,c,jac=True,bounds=m.bounds,constraints=[dict(type='ineq',fun=lambda c:obj.constraints(c)[0],jac=lambda c:obj.constraints(c)[1])],method='SLSQP',callback=checkpoint,options=dict(maxiter=maxiter,ftol=1e-12,disp=True))
 c=ret.x;vals,j=obj.constraints(c);raw=m.candidate(c);summary=dict(**reg,variables=m.dim,success=bool(ret.success),message=str(ret.message),iterations=ret.nit,nfev=ret.nfev,loss=float(ret.fun),min_constraint=float(vals.min()),elapsed=time.time()-t0,pde_validated=False)
 m.f.save(raw,out/'candidate.json',summary);summary['candidate_sha256']=hashlib.sha256((out/'candidate.json').read_bytes()).hexdigest();(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'history.json').write_text(json.dumps(obj.history,indent=2)+'\n');np.savez_compressed(out/'recipe.npz',c=c,Ma=m.Ma,Mb=m.Mb,Mp=m.Mp);print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--parent',required=True);p.add_argument('--out',required=True);p.add_argument('--profile-ratio',type=float,default=.75);p.add_argument('--pressure-ratio',type=float,default=.9);p.add_argument('--pressure-target',type=float);p.add_argument('--maxiter',type=int,default=240);p.add_argument('--radial',type=int,default=3);p.add_argument('--axial',type=int,default=2);p.add_argument('--time-degree',type=int,default=4);p.add_argument('--peak-weight',type=float,default=.05);p.add_argument('--localized',action='store_true');p.add_argument('--core-weight',type=float,default=0.);p.add_argument('--anchor-tolerance',type=float,default=.1);p.add_argument('--time-cap',type=float);p.add_argument('--softmax-weight',type=float,default=0.);a=p.parse_args();run(a.parent,a.out,a.profile_ratio,a.pressure_ratio,a.pressure_target,a.maxiter,a.radial,a.axial,a.time_degree,a.peak_weight,a.localized,a.core_weight,a.anchor_tolerance,a.time_cap,a.softmax_weight)
