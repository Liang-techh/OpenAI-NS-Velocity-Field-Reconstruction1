"""Frozen cylindrical phase/amplitude dynamics on actual ST073 background.
Local diagnostic only: not the paper's normalized leading field or a global wave.
"""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legfit,legint,legval
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.optimize import nnls
from angular_collocation import CachedBase
from poloidal_collocation import PoloidalModes
from affine_momentum import jets,momentum
from joined_field import ROOT

def integrate_pulse(F,g,r,tau,mode,axial_ratio,radial_ratio=-2):
 K=np.array([[0,-2*F,0],[2*F+g[0],0,0],[g[1],0,0]],float)
 nt=np.array([mode/r,axial_ratio*mode/r]);n0=np.r_[radial_ratio*np.linalg.norm(nt),nt]
 nd=np.array([-g@nt,0,0]);basis=null_space(n0[None,:]);duration=.1*tau
 times=np.linspace(0,duration,121)
 def rhs(t,flat):
  n=n0+t*nd;M=flat.reshape(3,2)
  op=-K+np.outer(n,n@K-nd)/(n@n)
  return (op@M).ravel()
 sol=solve_ivp(rhs,(0,duration),basis.ravel(),t_eval=times,rtol=2e-9,atol=1e-11)
 if not sol.success:raise RuntimeError(sol.message)
 undamped=sol.y.T.reshape(-1,3,2)
 damping=np.exp(-.01*((n0@n0)*times+(n0@nd)*times**2+(nd@nd)*times**3/3))
 matrices=undamped*damping[:,None,None]
 singular=np.linalg.svd(matrices,compute_uv=False)[:,0];peak_index=int(np.argmax(singular))
 _,_,vh=np.linalg.svd(matrices[peak_index],full_matrices=False);seed=vh[0];amp=np.einsum('nij,j->ni',matrices,seed)
 normals=n0[None,:]+times[:,None]*nd;norms=np.linalg.norm(amp,axis=1)
 cov=np.trapezoid(.5*amp[:,0,None]*amp[:,1:],times,axis=0)/duration
 raw_amp=np.einsum('nij,j->ni',undamped,seed)
 constraint=np.abs(np.sum(normals*raw_amp,axis=1))/(np.linalg.norm(normals,axis=1)*np.linalg.norm(raw_amp,axis=1))
 return dict(angular_mode=mode,axial_ratio=axial_ratio,radial_ratio=radial_ratio,duration=duration,final_gain=float(singular[-1]),peak_gain=float(singular.max()),transverse_relative_error=float(constraint.max()),mean_radial_tangential_covariance=cov.tolist(),initial_wavevector=n0.tolist(),wavevector_derivative=nd.tolist())

def run():
 base=CachedBase();fit=json.loads((ROOT/'poloidal_collocation/report.json').read_text());fields=[('baseline',base),('poloidal_pressure',PoloidalModes(base,np.array(fit['amplitudes']).ravel()))]
 rows=[];candidates=[];nodes,weights=leggauss(16)
 for label,f in fields:
  for k in (.4,5.5):
   tau=.5*2**(-k)
   for eta in (-.3,0,.3):
    ip=f.inner.from_similarity([3/64],[eta],tau)[0];ri=ip[0];r=ri*(1+(nodes+1)/2);pts=np.column_stack((r,np.zeros(len(r)),np.full(len(r),ip[2])))
    u,J,linear=jets(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau);res=momentum((u,J,linear))
    F=u[:,1]/r;g=np.column_stack((J[:,1,0]-F,J[:,2,0]));gn=np.linalg.norm(g,axis=1);N=g/gn[:,None]
    lam2=-2*F*N[:,0]*(2*F*N[:,0]+gn)
    # Gauss polynomial prefix quadrature at the same physical z.
    stress=[]
    for comp,power in ((1,2),(2,1)):
     anti=legint(legfit(nodes,r**power*res[:,comp],15))
     stress.append(-ri/2*(legval(nodes,anti)-legval(-1,anti))/r**power)
    stress=np.column_stack(stress)
    for i in range(len(r)):
     item=dict(field=label,k=k,eta=eta,y=float((r[i]-ri)/ri),r=float(r[i]),F=float(F[i]),tangential_shear=g[i].tolist(),lambda_squared=float(lam2[i]),radial_inverse_target=stress[i].tolist())
     if lam2[i]>0:
      lam=np.sqrt(lam2[i]);c=lam/(2*F[i]*N[i,0]);perp=np.array([-N[i,1],N[i,0]]);tn=float(stress[i]@N[i]);tk=float(stress[i]@perp)
      item.update(frozen_growth_rate=float(lam),target_dot_N=tn,cone_ratio=float(abs(c*tk/tn)) if tn!=0 else None,cone_like_pass=bool(tn<0 and abs(c*tk)<-tn))
      if label=='poloidal_pressure' and k==5.5:candidates.append((lam,item))
     else:item['cone_like_pass']=False
     rows.append(item)
    print('screened',label,k,eta,flush=True)
 selected=max(candidates,key=lambda x:x[0])[1] if candidates else None;pulses=[];covariance_fit=None
 if selected:
  for mode in (1,2,4):
   for ratio in (-2,-1,0,1,2):
    for radial in (-2,0,2):pulses.append(integrate_pulse(selected['F'],np.array(selected['tangential_shear']),selected['r'],.5*2**(-5.5),mode,ratio,radial))
  C=np.array([p['mean_radial_tangential_covariance'] for p in pulses]).T;target=np.array(selected['radial_inverse_target']);a,error=nnls(C,target)
  covariance_fit=dict(nonnegative_squared_weights=a.tolist(),target=target.tolist(),achieved=(C@a).tolist(),absolute_error=float(error),relative_error=float(error/np.linalg.norm(target)))
 report=dict(rows=rows,selected=selected,pulses=pulses,covariance_fit=covariance_fit,
 equations='Frozen physical cylindrical K=[[0,-2F,0],[2F+gtheta,0,0],[gz,0,0]]. n_dot=(-g dot n_tan,0,0). a_dot=(-K+n(n^T K-n_dot^T)/|n|^2-nu|n|^2 I)a. Includes pressure projection and viscosity.',
 selection='Largest positive frozen growth rate among sampled late poloidal/pressure candidate points. Modes1,2,4, axial ratios -2,-1,0,1,2, radial ratios -2,0,2. Exact scalar viscous damping factored from undamped ODE; seed maximizes sampled peak gain. Nonzero unit initial seeds, NOT zero-start supported pulses. Frozen duration .1*tau.',
 limitations='Paper-inspired diagnostic in physical coordinates, not direct application of the leading normalized theorem or its cone. Omits radial strain, axial gradients and changing background along rays. Prefix stress is restricted radial inverse and not compact at outer boundary. Averaged covariance fit is NOT a full PDE cancellation or realizable global field.',pde_validated=False,global_field_ready=False)
 out=ROOT/'frozen_pulse';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
 print('positive growth',sum(r['lambda_squared']>0 for r in rows),'of',len(rows),flush=True)
 print('cone-like pass',sum(r['cone_like_pass'] for r in rows),flush=True)
 print('peak gain',max((p['peak_gain'] for p in pulses),default=None),'covariance fit',covariance_fit,flush=True)
if __name__=='__main__':run()
