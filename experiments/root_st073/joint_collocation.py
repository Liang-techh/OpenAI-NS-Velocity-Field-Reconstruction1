"""Joint swirl, streamfunction and pressure fit with full momentum and moments."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from angular_collocation import CachedBase,AngularModes,sample,metrics
from poloidal_collocation import PoloidalModes
from affine_momentum import jets,momentum,combine
from joined_field import ROOT,independent_fd

class JointModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.a=np.asarray(amplitudes)
  self.field=AngularModes(PoloidalModes(base,self.a[:12]),self.a[12:])
 def fields(self,points,tau):return self.field.fields(points,tau)

def run():
 base=CachedBase();zero=JointModes(base,np.zeros(24));train=[(k,e) for k in (1.,4.) for e in (-.3,0,.3)];data=[]
 for k,e in train:
  tau=.5*2**(-k);g,w=leggauss(12);ip=base.inner.from_similarity([3/64],[e],tau)[0];r=ip[0]*(1+(g+1)/2)
  pts=np.column_stack((r,np.zeros(12),np.full(12,ip[2])));args=(pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau)
  jb=jets(zero,*args);mode=[]
  for j in range(24):
   a=np.zeros(24);a[j]=1;ju=jets(JointModes(base,a),*args);mode.append(tuple(x-y for x,y in zip(ju,jb)))
  modes=tuple(np.stack([v[i] for v in mode]) for i in range(3));scale=np.max(np.linalg.norm(momentum(jb),axis=1))
  wt=w*r*r;wt/=wt.sum();wz=w*r;wz/=wz.sum()
  data.append((jb,modes,scale,wt,wz));print('assembled',k,e,flush=True)
 def objective(a):
  values=[]
  for b,m,s,wt,wz in data:
   res=momentum(combine(b,m,a))/s
   values.append(np.r_[res.ravel(),3*wt@res[:,1],3*wz@res[:,2]])
  return np.r_[np.concatenate(values),.001*a]
 def jacobian(a):
  rows=[]
  for b,m,s,wt,wz in data:
   u,g,_=combine(b,m,a);mu,mg,ml=m
   jac=(ml+np.einsum('knij,nj->kni',mg,u)+np.einsum('nij,knj->kni',g,mu))/s
   rows.append(np.vstack((jac.reshape(24,-1).T,3*np.einsum('n,kn->k',wt,jac[:,:,1]),3*np.einsum('n,kn->k',wz,jac[:,:,2]))))
  return np.vstack((*rows,.001*np.eye(24)))
 start=np.r_[np.array(json.loads((ROOT/'poloidal_collocation/report.json').read_text())['amplitudes']).ravel(),np.zeros(12)]
 fit=least_squares(objective,start,jac=jacobian,bounds=(-4,4),max_nfev=150,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 corrected=JointModes(base,fit.x);previous=JointModes(base,start);rows=[]
 for k,e in [(k,e) for k in (2.5,5.5) for e in (-.2,.2)]:
  row=dict(k=k,eta=e,order=18)
  for label,f in [('baseline',zero),('poloidal_pressure',previous),('joint',corrected)]:
   r,d,mw=sample(f,k,e,18,.0005);row[label]=metrics(r,d,mw)
   radial_y=(leggauss(18)[0]+1)/2
   row[label]['terminal_rz']=float(-np.sum(mw*2/(1+radial_y)*r[:,2]))
  rows.append(row);print(json.dumps(row),flush=True)
 rb,db,_=sample(corrected,4.,.3,12,.0005);predicted=momentum(combine(data[-1][0],data[-1][1],fit.x))
 # Finite-difference check of analytic optimizer Jacobian at fitted coefficients.
 eps=1e-5;direction=np.arange(1,25,dtype=float);direction/=np.linalg.norm(direction)
 jac_error=float(np.max(np.abs((objective(fit.x+eps*direction)-objective(fit.x-eps*direction))/(2*eps)-jacobian(fit.x)@direction)))
 volume=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k);gy,wy=leggauss(12);ge,we=leggauss(8);eta=np.repeat(.4*ge,12);y=np.tile((gy+1)/2,8);q=tau/(1-eta**2)
  ri=np.sqrt(2*base.nu*q*3/64);r=ri*(1+y);z=np.sqrt(base.nu)*q**base.inner.D*eta
  ze=np.sqrt(base.nu)*q**base.inner.D*(1+2*base.inner.D*eta**2/(1-eta**2));weights=np.repeat(.4*we,12)*np.tile(wy/2,8)*2*np.pi*r*ri*ze
  pts=np.column_stack((r,np.zeros_like(r),z))
  for label,f in [('baseline',zero),('poloidal_pressure',previous),('joint',corrected)]:
   res,div=independent_fd(f,pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau)
   item=dict(k=k,field=label,volume=float(weights.sum()),momentum_sampled_max=float(np.max(np.linalg.norm(res,axis=1))),momentum_volume_L2=float(np.sqrt(np.sum(weights*np.sum(res**2,axis=1)))),divergence_sampled_max=float(np.max(np.abs(div))))
   volume.append(item);print(json.dumps(item),flush=True)
 report=dict(amplitudes=fit.x.tolist(),coefficient_order='12 PoloidalModes followed by12 AngularModes',bounds=[-4,4],fit_nfev=fit.nfev,fit_success=bool(fit.success),fit_cost=float(fit.cost),initial_cost=float(.5*np.sum(objective(start)**2)),holdouts=rows,volume_audit=volume,surrogate_full_fd_max_difference=float(np.max(np.abs(rb-predicted))),optimizer_jacobian_directional_error=jac_error,
 objective='Full Cartesian residual normalized per training slice by baseline max, plus weight3 normalized radial theta and axial moment averages, regularization .001.',
 volume_domain='eta in[-.4,.4], ri<r<2ri, all azimuths; 12x8 Gauss quadrature with physical Jacobian. No continuum maximum/integral certificate.',
 scope='24 fixed-scale coefficients, unforced axisymmetric compact corrections. Frozen core and interface jets unchanged. Does not construct global finite energy, axial closure or oscillatory stress realization.',pde_validated=False,global_field_ready=False)
 out=ROOT/'joint_collocation';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
