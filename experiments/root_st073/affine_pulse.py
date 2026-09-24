"""Full frozen affine velocity-gradient ray dynamics at directional-screen points."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from scipy.optimize import nnls
from frozen_pulse import integrate_pulse
from affine_momentum import jets
from angular_collocation import CachedBase
from poloidal_collocation import PoloidalModes
from joined_field import ROOT

def affine_pulse(J,nu,n0,duration):
 basis=null_space(np.asarray(n0)[None,:]);times=np.linspace(0,duration,121)
 def rhs(t,state):
  n=state[:3];M=state[3:9].reshape(3,2)
  op=-J+2*np.outer(n,n@J)/(n@n)
  return np.r_[-J.T@n,(op@M).ravel(),nu*(n@n)]
 sol=solve_ivp(rhs,(0,duration),np.r_[n0,basis.ravel(),0.],t_eval=times,rtol=2e-9,atol=1e-11)
 if not sol.success:raise RuntimeError(sol.message)
 normals=sol.y[:3].T;undamped=sol.y[3:9].T.reshape(-1,3,2);M=undamped*np.exp(-sol.y[9])[:,None,None]
 gains=np.linalg.svd(M,compute_uv=False)[:,0];index=int(np.argmax(gains));_,_,vh=np.linalg.svd(M[index],full_matrices=False);seed=vh[0]
 raw=undamped@seed;amp=M@seed;constraint=np.abs(np.sum(normals*raw,axis=1))/(np.linalg.norm(normals,axis=1)*np.linalg.norm(raw,axis=1))
 cov=np.trapezoid(.5*amp[:,0,None]*amp[:,1:],times,axis=0)/duration
 return dict(peak_gain=float(gains.max()),final_gain=float(gains[-1]),transverse_relative_error=float(constraint.max()),mean_radial_tangential_covariance=cov.tolist(),initial_wavevector=np.asarray(n0).tolist(),peak_time=float(times[index]))

def run():
 screen=json.loads((ROOT/'frozen_pulse/report.json').read_text());base=CachedBase();a=json.loads((ROOT/'poloidal_collocation/report.json').read_text())['amplitudes'];f=PoloidalModes(base,np.array(a).ravel())
 selected=[r for r in screen['rows'] if r['field']=='poloidal_pressure' and r['k']==5.5 and r['cone_like_pass']];rows=[]
 for point in selected:
  tau=.5*2**(-point['k']);eta=point['eta'];ri=point['r']/(1+point['y']);z=f.inner.from_similarity([3/64],[eta],tau)[0,2];pts=np.array([[point['r'],0,z]])
  u,J,_=jets(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau);J=J[0];duration=.1*tau;pulses=[]
  for mode in (1,4,16,64,256):
   for axial in (-1,0,1):
    for radial in (-1,0,1):
     nt=np.array([mode/point['r'],axial*mode/point['r']]);n0=np.r_[radial*np.linalg.norm(nt),nt]
     full=affine_pulse(J,f.nu,n0,duration)
     reduced=integrate_pulse(point['F'],np.array(point['tangential_shear']),point['r'],tau,mode,axial,radial)
     pulses.append(dict(mode=mode,axial_ratio=axial,radial_ratio=radial,full=full,reduced_peak_gain=reduced['peak_gain']))
  C=np.array([p['full']['mean_radial_tangential_covariance'] for p in pulses]).T;target=np.array(point['radial_inverse_target']);weights,error=nnls(C,target)
  distance=min(point['r']-ri,2*ri-point['r']);rate=point['frozen_growth_rate']
  row=dict(point=point,velocity_gradient=J.tolist(),gradient_trace=float(np.trace(J)),pulses=pulses,covariance_fit=dict(weights=weights.tolist(),relative_error=float(error/np.linalg.norm(target)),target=target.tolist(),achieved=(C@weights).tolist()),
   localization=dict(distance_to_inner_interface=distance,diffusive_rate_nu_over_distance_squared=float(f.nu/distance**2),frozen_growth_rate=rate,diffusion_to_growth_ratio=float(f.nu/distance**2/rate),note='Dimensional cutoff cost estimate only; no rigorous localized operator bound.'),
   frozen_validity=dict(duration=duration,estimated_advective_displacement=(u[0]*duration).tolist(),radial_displacement_to_interface_distance=float(abs(u[0,0]*duration)/distance)))
  rows.append(row);print('eta',eta,'peak full',max(p['full']['peak_gain'] for p in pulses),'covariance relative error',row['covariance_fit']['relative_error'],'localization',row['localization'],flush=True)
 report=dict(rows=rows,equations='n_dot=-J^T n; a_dot=-J a+2 n(n^T J a)/|n|^2-nu|n|^2 a. Integrate full frozen physical Cartesian J, factor scalar damping using D_dot=nu|n|^2.',
 scope='Actual late candidate points passing previous approximate directional screen. Constant-J local affine rays, not Lagrangian integration in changing ST073. Initial angular-like modes are local wavevectors and do not establish global angular periodicity under Cartesian affine evolution. Nonzero seeds, no compact support or zero-start pulse.',
 interpretation='Positive local growth and nonnegative covariance fit are separate from supported divergence-free wave construction and full residual cancellation. Cutoff and frozen-validity estimates flag omitted terms; not an impossibility proof.',pde_validated=False,global_field_ready=False)
 out=ROOT/'affine_pulse';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
