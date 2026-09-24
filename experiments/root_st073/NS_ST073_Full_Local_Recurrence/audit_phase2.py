from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from full_radial import *
from audit import save,ROOT
sys.path.insert(0,str(ROOT/'upstream'))
from general_core import build
from local_field import LocalField,independent_fd

def binding():
 fr=json.loads((ROOT/'evidence/phase2_freeze.json').read_text())
 for p,h in [('full_radial.py',fr['source_sha256']),('data/ST073-V.json',fr['model_sha256'])]:
  if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise ValueError('Frozen input changed')
 return fr

def holdout(seed):
 fr=binding();f=FullRadialField.load(ROOT/'data/ST073-V.json');rng=np.random.default_rng(seed)
 eta=rng.uniform(-.5,.5,32);X=rng.uniform(0,1/64,(32,16));es=np.broadcast_to(eta[:,None],X.shape);ks=np.r_[0,3,6,np.sort(rng.uniform(.02,5.98,5))]
 for j,k in enumerate(ks):
  path=ROOT/'evidence'/f'V_holdout_{seed}_{j}.json'
  if path.exists():continue
  t=time.monotonic();d=f.evaluate_similarity(X,es,.5*2**(-k));rr=np.linalg.norm(d['residual'],axis=-1)
  vel=d['velocity'];pg=d['pressure_gradient'];notplane=abs(es)>.02
  dd=dict(radial_inflow=float(np.mean(vel[...,0]<0)),positive_swirl=float(np.mean(vel[...,1]>0)),bipolar_outflow=float(np.mean((es*vel[...,2])[notplane]>0)),inward_radial_pressure=float(np.mean(pg[...,0]>0)),pressure_toward_midplane=float(np.mean((es*pg[...,2])[notplane]>0)))
  data=ROOT/'data'/f'V_holdout_{seed}_{j}.npz';np.savez_compressed(data,X=X,eta=es,k=k,**d)
  row=dict(seed=seed,index=j,k=float(k),points=X.size,sampling='32 independent eta sites x16 radial points, repeated at8log times. 4096clustered spacetime points per seed, not IID whole-space samples.',full_sampled_max=float(rr.max()),divergence_max=float(abs(d['divergence']).max()),direction_fractions=dd,model_sha256=fr['model_sha256'],data_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),seconds=time.monotonic()-t,global_pde_validated=False)
  save(path,row);print('V HOLDOUT',seed,j,k,rr.max(),dd,flush=True)

def baseline():
 path=ROOT/'data/ST073-V-leading-control.npz'
 if path.exists():
  from general_core import Core
  return Core.load(path)
 c,h=build(order=18,eta_degree=64,swirl=4.)
 c.meta['id']='ST073-V-leading-control';c.meta['scope']='Same axis data as V, leading-only comparison; not a new full solution'
 c.save(path);save(ROOT/'evidence/leading_control_build.json',{'history':h,'source':'unchanged inherited general_core.py','scope':c.meta['scope']})
 return c

def physical(k,n):
 fr=binding();f=FullRadialField.load(ROOT/'data/ST073-V.json');old=LocalField(baseline());tau=.5*2**(-k)
 g,w=leggauss(n);xx,ee=np.meshgrid((g+1)/128,g/2,indexing='ij');q=tau/(1-ee**2);weights=w[:,None]*w[None,:]/256*2*np.pi*.01**1.5*q**(1+f.D)*(1-2*f.h*ee**2)/(1-ee**2)
 d=f.evaluate_similarity(xx,ee,tau);b=old.evaluate(f.from_similarity(xx.ravel(),ee.ravel(),tau),tau)
 row=dict(id='ST073-V',k=k,quadrature_order=n,points=xx.size,volume=float(weights.sum()),scope='Same local physical region and same axis traces; no global energy normalization')
 for name,r in [('full',d),('leading-only same-axis',b)]:
  R=r['residual'].reshape(*xx.shape,3);u=r['velocity'].reshape(*xx.shape,3);omega=r['vorticity'].reshape(*xx.shape,3);en=np.sum(u*u,axis=-1);vo=np.sum(omega*omega,axis=-1)
  row[name]=dict(L2=float(np.sqrt(np.sum(weights*np.sum(R*R,axis=-1)))),max=float(np.linalg.norm(R,axis=-1).max()),energy=float(.5*np.sum(weights*en)),energy_effective_volume=float(np.sum(weights*en)**2/np.sum(weights*en**2)),enstrophy_effective_volume=float(np.sum(weights*vo)**2/np.sum(weights*vo**2)))
 return row

def fd():
 binding();f=FullRadialField.load(ROOT/'data/ST073-V.json');out=[]
 for k in [.4,2.7,5.5]:
  tau=.5*2**(-k);pts=f.from_similarity(np.array([.006,.012]),np.array([-.22,.37]),tau,angle=[.3,.8]);exact=f.evaluate(pts,tau)['residual'];base=np.sqrt(.01*tau);rows=[]
  for kind in ['space','time']:
   for eps in [.004,.002,.001]:
    hs=base*(eps if kind=='space' else .001);ht=tau*(eps if kind=='time' else .0005);t=time.monotonic();R,div=independent_fd(f,pts,tau,hs,ht)
    row=dict(kind=kind,epsilon=eps,hspace=hs,htime=ht,max_error=float(np.linalg.norm(R-exact,axis=-1).max()),max_FD_residual=float(np.linalg.norm(R,axis=-1).max()),divergence_max=float(abs(div).max()),FD=R.tolist(),seconds=time.monotonic()-t);rows.append(row);print('V FD',k,kind,eps,row['max_error'],flush=True)
  out.append(dict(k=k,points=pts.tolist(),analytic=exact.tolist(),levels=rows))
 return dict(rows=out,scope='Independent Cartesian operator; spacings independently varied, six fixed physical points; not a continuum bound')

if __name__=='__main__':
 for seed in [9237395,9237396]:holdout(seed)
 for k in [0,3,6]:
  for n in [8,12,18]:
   p=ROOT/'evidence'/f'V_physical_k{k}_n{n}.json'
   if p.exists():continue
   row=physical(k,n);save(p,row);print('V PHYS',k,n,row,flush=True)
 p=ROOT/'evidence/V_independent_fd.json'
 if not p.exists():save(p,fd())
 print('V AUDITS FINISHED',flush=True)
