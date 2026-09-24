"""Matched-data comparison of constant and log-time scale-dependent corrections."""
import json
import numpy as np
from scipy.optimize import least_squares
from numpy.polynomial.legendre import leggauss
from joint_collocation import JointModes,CachedBase,sample,metrics,ROOT,independent_fd
from affine_momentum import jets,momentum,combine

class ScaleModes:
 def __init__(self,base,amplitudes):
  self.base=base;self.inner=base.inner;self.nu=base.nu;self.a=np.asarray(amplitudes).reshape(2,24)
  self.constant=JointModes(base,self.a[0]);self.slope=JointModes(base,self.a[1])
 def fields(self,points,tau):
  u,p=self.constant.fields(points,tau);v,q=self.slope.fields(points,tau);ub,pb=self.base.fields(points,tau)
  s=(-np.log2(2*np.broadcast_to(tau,(len(u),)))-3)/3
  return u+s[:,None]*(v-ub),p+s*(q-pb)

def run():
 base=CachedBase();zero=ScaleModes(base,np.zeros(48));data=[]
 for k in (1.,3.,5.):
  for e in (-.3,0,.3):
   tau=.5*2**(-k);g,w=leggauss(12);ip=base.inner.from_similarity([3/64],[e],tau)[0];r=ip[0]*(1+(g+1)/2)
   pts=np.column_stack((r,np.zeros(12),np.full(12,ip[2])));args=(pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau)
   jb=jets(zero,*args);mode=[]
   # Constant-mode jets; scale-mode derivatives include ds/dt explicitly.
   for j in range(24):
    a=np.zeros(24);a[j]=1;ju=jets(JointModes(base,a),*args);mode.append(tuple(x-y for x,y in zip(ju,jb)))
   s=(k-3)/3;st=1/(3*np.log(2)*tau)
   slope=[(s*u,s*g,s*l+st*u) for u,g,l in mode]
   modes=tuple(np.stack([v[i] for v in mode+slope]) for i in range(3));scale=np.max(np.linalg.norm(momentum(jb),axis=1))
   wt=w*r*r;wt/=wt.sum();wz=w*r;wz/=wz.sum();data.append((jb,modes,scale,wt,wz))
   print('assembled',k,e,flush=True)
 def make_problem(dim):
  def objective(a):
   full=np.r_[a,np.zeros(48-dim)];values=[]
   for b,m,s,wt,wz in data:
    res=momentum(combine(b,m,full))/s
    values.append(np.r_[res.ravel(),3*wt@res[:,1],3*wz@res[:,2]])
   return np.r_[np.concatenate(values),.001*a]
  def jacobian(a):
   full=np.r_[a,np.zeros(48-dim)];rows=[]
   for b,m,s,wt,wz in data:
    u,g,_=combine(b,m,full);mu,mg,ml=(v[:dim] for v in m)
    jac=(ml+np.einsum('knij,nj->kni',mg,u)+np.einsum('nij,knj->kni',g,mu))/s
    rows.append(np.vstack((jac.reshape(dim,-1).T,3*np.einsum('n,kn->k',wt,jac[:,:,1]),3*np.einsum('n,kn->k',wz,jac[:,:,2]))))
   return np.vstack((*rows,.001*np.eye(dim)))
  return objective,jacobian
 prior=np.array(json.loads((ROOT/'joint_collocation/report.json').read_text())['amplitudes']);obj,jac=make_problem(24)
 fixed=least_squares(obj,prior,jac=jac,bounds=(-4,4),max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 obj,jac=make_problem(48);vary=least_squares(obj,np.r_[fixed.x,np.zeros(24)],jac=jac,bounds=(-4,4),max_nfev=120,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 fields=[('baseline',zero),('matched_fixed',ScaleModes(base,np.r_[fixed.x,np.zeros(24)])),('scale_dependent',ScaleModes(base,vary.x))]
 rows=[]
 for k in (2.5,5.5):
  for e in (-.2,.2):
   row=dict(k=k,eta=e)
   for label,f in fields:
    res,div,mw=sample(f,k,e,18,.0005);item=metrics(res,div,mw);y=(leggauss(18)[0]+1)/2
    item['terminal_rz']=float(-np.sum(mw*2/(1+y)*res[:,2]));row[label]=item
   rows.append(row);print(json.dumps(row),flush=True)
 rr,dd,mw=sample(fields[-1][1],5.,.3,12,.0005)
 surrogate_error=float(np.max(np.abs(rr-momentum(combine(data[-1][0],data[-1][1],vary.x)))))
 volume=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k);gy,wy=leggauss(12);ge,we=leggauss(8);eta=np.repeat(.4*ge,12);y=np.tile((gy+1)/2,8);q=tau/(1-eta**2)
  ri=np.sqrt(2*base.nu*q*3/64);r=ri*(1+y);z=np.sqrt(base.nu)*q**base.inner.D*eta
  ze=np.sqrt(base.nu)*q**base.inner.D*(1+2*base.inner.D*eta**2/(1-eta**2));weights=np.repeat(.4*we,12)*np.tile(wy/2,8)*2*np.pi*r*ri*ze
  pts=np.column_stack((r,np.zeros_like(r),z))
  for label,f in fields:
   res,div=independent_fd(f,pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau)
   item=dict(k=k,field=label,volume=float(weights.sum()),momentum_sampled_max=float(np.max(np.linalg.norm(res,axis=1))),momentum_volume_L2=float(np.sqrt(np.sum(weights*np.sum(res**2,axis=1)))),divergence_sampled_max=float(np.max(np.abs(div))))
   volume.append(item);print(json.dumps(item),flush=True)
 report=dict(fixed_amplitudes=fixed.x.tolist(),scale_amplitudes=vary.x.reshape(2,24).tolist(),fits=dict(fixed=dict(cost=float(fixed.cost),nfev=fixed.nfev,success=bool(fixed.success)),scale=dict(cost=float(vary.cost),nfev=vary.nfev,success=bool(vary.success))),holdouts=rows,volume_audit=volume,surrogate_full_fd_max_difference=surrogate_error,
 scale='a(k)=a0+((k-3)/3)*a1, k=-log2(2tau). s_t=1/(3 ln(2) tau) for increasing physical time. Amplitude is spatially constant at fixed time, preserving divergence.',
 bounds='Each intercept and slope in[-4,4]; effective amplitudes need not stay in[-4,4]. Training k=1,3,5, eta=-.3,0,.3; held-out k=2.5,5.5, eta=+/-.2.',
 derivative='Scale-mode jets use product rule with analytic s_t; independent full-field FD includes variable amplitudes at every time stencil. This is time dependence, NOT a proved scale recursion.',
 scope='Same unforced frozen core and interface jets; finite annular/axial slab only. Volume quadrature12x8 on eta in[-.4,.4], ri<r<2ri, all azimuths. No whole-space energy, axial closure, continuum max or converged quadrature certificate.',pde_validated=False,global_field_ready=False)
 out=ROOT/'scale_collocation';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
