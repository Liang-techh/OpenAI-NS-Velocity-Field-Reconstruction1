"""Frozen, scale-resolved independent tests of the local field; no fitting."""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from numpy.polynomial.legendre import leggauss
from full_radial import FullRadialField,Parameters
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'upstream'))
from core_series import Core
from local_field import LocalField,independent_fd

def check():
 fr=json.loads((ROOT/'evidence/freeze.json').read_text())
 for p,h in [('full_radial.py',fr['source_sha256']),('data/ST073-F.json',fr['model_sha256'])]:
  if hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h:raise ValueError('Frozen source/parameters changed')
 return fr

def save(p,obj):
 p=Path(p)
 if p.exists():raise FileExistsError(p)
 p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(p)

def holdout(seed):
 fr=check();f=FullRadialField.load(ROOT/'data/ST073-F.json');rng=np.random.default_rng(seed)
 eta=rng.uniform(-.5,.5,64);X=rng.uniform(0,1/64,(64,8));es=np.broadcast_to(eta[:,None],X.shape)
 kk=np.r_[0,3,6,np.sort(rng.uniform(.02,5.98,5))]
 for j,k in enumerate(kk):
  dest=ROOT/'evidence'/f'holdout_{seed}_{j}.json'
  if dest.exists():continue
  t=time.monotonic();out=f.evaluate_similarity(X,es,.5*2**(-k));norm=np.linalg.norm(out['residual'],axis=-1)
  data=ROOT/'data'/f'holdout_{seed}_{j}.npz';np.savez_compressed(data,X=X,eta=es,k=k,velocity=out['velocity'],pressure=out['pressure'],residual=out['residual'],divergence=out['divergence'])
  row=dict(seed=seed,index=j,k=float(k),points=X.size,sampling='64 independent axial sites, eight new radial points per site; same spatial sites at eight registered/new log-times; 4096 clustered spacetime points per seed',full_sampled_max=float(norm.max()),divergence_max=float(abs(out['divergence']).max()),sampled_local_momentum_pass=bool(norm.max()<.001),model_sha256=fr['model_sha256'],data_sha256=hashlib.sha256(data.read_bytes()).hexdigest(),seconds=time.monotonic()-t,pde_validated=False,scope='Unforced local domain only, not global NS acceptance')
  save(dest,row);print('HOLDOUT',seed,j,k,row['full_sampled_max'],row['seconds'],flush=True)

def physical(k,n):
 fr=check();f=FullRadialField.load(ROOT/'data/ST073-F.json');b=LocalField(Core.load(ROOT/'data/ST068-I.npz'))
 g,w=leggauss(n);X=(g+1)/128;wx=w/128;e=g/2;we=w/2;xx,ee=np.meshgrid(X,e,indexing='ij');weights=wx[:,None]*we[None,:]
 tau=.5*2**(-k);q=tau/(1-ee**2);weights=weights*2*np.pi*.01**1.5*q**(1+f.D)*(1-2*f.h*ee**2)/(1-ee**2)
 d=f.evaluate_similarity(xx,ee,tau);pts=f.from_similarity(xx.ravel(),ee.ravel(),tau);old=b.evaluate(pts,tau)
 row=dict(k=k,quadrature_order=n,point_count=xx.size,region_volume=float(weights.sum()),X_max=1/64,eta_max=.5,nu=.01,force=0,model_sha256=fr['model_sha256'],scope='Identical source-coordinate physical region and ST068 axis traces; no global normalization, force or matching')
 for name,out in [('ST073-F',d),('ST068-I restricted',old)]:
  R=out['residual'].reshape(*xx.shape,3);U=out['velocity'].reshape(*xx.shape,3);sq=np.sum(R*R,axis=-1)
  row[name]=dict(full_local_L2=float(np.sqrt(np.sum(weights*sq))),sampled_max=float(np.sqrt(sq.max())),component_L2=np.sqrt(np.sum(weights[:,:,None]*out['residual_cylindrical'].reshape(*xx.shape,3)**2,axis=(0,1))).tolist(),local_energy=float(.5*np.sum(weights*np.sum(U*U,axis=-1))),local_velocity_L2=float(np.sqrt(np.sum(weights*np.sum(U*U,axis=-1)))),divergence_max=float(abs(out['divergence']).max()))
 # Physical differences are evaluated on the same volume, not an amplitude fit.
 u=d['velocity'];v=old['velocity'].reshape(*xx.shape,3)
 row['relative_velocity_L2_change']=float(np.sqrt(np.sum(weights*np.sum((u-v)**2,axis=-1))/np.sum(weights*np.sum(v*v,axis=-1))))
 return row

def fd_audit():
 check();f=FullRadialField.load(ROOT/'data/ST073-F.json');points=[(.006,-.22,.3),(.012,.37,.8)];out=[]
 for k in [.4,2.7,5.5]:
  tau=.5*2**(-k);pts=f.from_similarity(np.array([p[0] for p in points]),np.array([p[1] for p in points]),tau,np.array([p[2] for p in points]));exact=f.evaluate(pts,tau)['residual']
  base=np.sqrt(.01*tau)
  # Vary space and time separately, retaining fixed physical points.
  rows=[]
  for kind in ['space','time']:
   for eps in [.004,.002,.001]:
    hs=base*(eps if kind=='space' else .001);ht=tau*(eps if kind=='time' else .0005)
    t=time.monotonic();R,div=independent_fd(f,pts,tau,hs,ht)
    rows.append(dict(kind=kind,fraction=eps,space_step=hs,time_step=ht,max_absolute_vector_difference=float(np.linalg.norm(R-exact,axis=-1).max()),fd_sampled_max=float(np.linalg.norm(R,axis=-1).max()),fd_divergence_max=float(abs(div).max()),values=R.tolist(),seconds=time.monotonic()-t));print('FD',k,kind,eps,rows[-1]['max_absolute_vector_difference'],flush=True)
  out.append(dict(k=k,point_count=2,physical_points=pts.tolist(),analytic_values=exact.tolist(),levels=rows))
 return dict(rows=out,scope='Independent fourth-order Cartesian/time derivatives at six physical points, not an error enclosure or whole-domain validation')

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--part',choices=['holdout','physical','fd'],required=True);p.add_argument('--seed',type=int,default=9237391);p.add_argument('--k',type=float,default=6);p.add_argument('--n',type=int,default=18);a=p.parse_args()
 if a.part=='holdout':holdout(a.seed)
 elif a.part=='physical':
  r=physical(a.k,a.n);save(ROOT/'evidence'/f'physical_k{a.k:g}_n{a.n}.json',r);print(json.dumps(r,indent=2))
 else:save(ROOT/'evidence/independent_fd.json',fd_audit())
