"""Kelvin phase and viscous amplitude evolution along an actual ST073 trajectory."""
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.linalg import null_space
from affine_pulse import affine_pulse
from wide_modes import CachedWidth,WideJointModes
from joined_field import ROOT

def evolving_pulse(gradient,nu,n0,duration):
 basis=null_space(np.asarray(n0)[None,:]);times=np.linspace(0,duration,161)
 def rhs(t,state):
  J=gradient(t);n=state[:3];M=state[3:9].reshape(3,2)
  projected=-J+2*np.outer(n,n@J)/(n@n)
  return np.r_[-J.T@n,(projected@M).ravel(),nu*(n@n)]
 sol=solve_ivp(rhs,(0,duration),np.r_[n0,basis.ravel(),0.],t_eval=times,rtol=2e-9,atol=1e-11)
 if not sol.success:raise RuntimeError(sol.message)
 n=sol.y[:3].T;raw=sol.y[3:9].T.reshape(-1,3,2);M=raw*np.exp(-sol.y[9])[:,None,None]
 gains=np.linalg.svd(M,compute_uv=False)[:,0];imax=int(np.argmax(gains));_,_,vh=np.linalg.svd(M[imax],full_matrices=False);seed=vh[0]
 unattenuated=raw@seed;err=np.abs(np.sum(n*unattenuated,axis=1))/(np.linalg.norm(n,axis=1)*np.linalg.norm(unattenuated,axis=1))
 return dict(peak_gain=float(gains.max()),final_gain=float(gains[-1]),peak_time_fraction=float(times[imax]/duration),max_transversality_error=float(err.max()),final_wavevector=n[-1].tolist(),damping_exponent_at_exit=float(sol.y[9,-1]))

def run():
 source=json.loads((ROOT/'wide_pulse_screen/report.json').read_text());trace=source['trajectory'];a=json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes'];f=WideJointModes(CachedWidth(6),a)
 times=np.array(trace['times']);path=np.array(trace['path']);tau0=.5*2**(-trace['selected']['k']);h=.0005*np.sqrt(f.nu*tau0)
 gradients=np.empty((len(times),3,3));remaining=tau0-times
 for i in range(3):
  e=np.zeros(3);e[i]=h
  m2=f.fields(path-2*e,remaining)[0];m1=f.fields(path-e,remaining)[0];p1=f.fields(path+e,remaining)[0];p2=f.fields(path+2*e,remaining)[0]
  gradients[:,:,i]=(m2-8*m1+8*p1-p2)/(12*h)
  print('gradient direction',i,flush=True)
 Jof=PchipInterpolator(times,gradients,axis=0);duration=times[-1];r0=np.linalg.norm(path[0,:2]);rows=[]
 for mode in (1,4,16):
  for axial in (-1,0,1):
   for radial in (-1,0,1):
    ntheta=mode/r0;nt=np.array([ntheta,axial*ntheta]);n0=np.r_[radial*np.linalg.norm(nt),nt]
    evolving=evolving_pulse(Jof,f.nu,n0,duration);frozen=affine_pulse(gradients[0],f.nu,n0,duration)
    rows.append(dict(mode=mode,axial_ratio=axial,radial_ratio=radial,evolving=evolving,frozen=dict(peak_gain=frozen['peak_gain'],final_gain=frozen['final_gain'])))
 report=dict(rows=rows,trajectory_source='wide_pulse_screen/report.json',gradient_times=times.tolist(),gradient_matrices=gradients.tolist(),gradient_trace_max=float(np.max(np.abs(np.trace(gradients,axis1=1,axis2=2)))),spatial_fd_step=h,
 equations='Physical Cartesian Kelvin model along actual background particle: n_dot=-J(t)^T n, a_dot=-J(t)a+2n(n^T J(t)a)/|n|²-nu|n|²a. Matrix fundamental solution on n-perpendicular initial plane; damping exponent integrated exactly as separate scalar.',
 scope='J(t) from25-point actual trajectory, fourth-order spatial finite differences and shape-preserving interpolation. 27 sampled initial wavevectors. This is local linearized dynamics, not a global angular-periodic phase or a supported nonlinear wave. Initial approximate stress direction has not been verified along the whole path.',pde_validated=False,global_field_ready=False)
 out=ROOT/'evolving_pulse';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
 print('max evolving gain',max(x['evolving']['peak_gain'] for x in rows),'max frozen gain',max(x['frozen']['peak_gain'] for x in rows),flush=True)
if __name__=='__main__':run()
