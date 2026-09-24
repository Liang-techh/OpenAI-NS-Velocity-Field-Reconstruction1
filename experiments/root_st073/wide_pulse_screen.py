"""Recomputed directional geometry and actual residence for corrected width6 field."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legfit,legint,legval
from scipy.integrate import solve_ivp
from wide_modes import CachedWidth,WideJointModes
from affine_momentum import jets,momentum
from affine_pulse import affine_pulse
from joined_field import ROOT,coordinates

def run():
 a=json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes'];f=WideJointModes(CachedWidth(6),a);g,w=leggauss(20);rows=[];candidates=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k)
  for eta in (-.3,0,.3):
   ip=f.inner.from_similarity([3/64],[eta],tau)[0];ri=ip[0];width=5*ri;r=ri+width*(g+1)/2;pts=np.column_stack((r,np.zeros(len(r)),np.full(len(r),ip[2])))
   u,J,linear=jets(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau);res=momentum((u,J,linear));stress=[]
   for comp,power in ((1,2),(2,1)):
    anti=legint(legfit(g,r**power*res[:,comp],19));stress.append(-width/2*(legval(g,anti)-legval(-1,anti))/r**power)
   stress=np.column_stack(stress);F=u[:,1]/r;shear=np.column_stack((J[:,1,0]-F,J[:,2,0]));gn=np.linalg.norm(shear,axis=1);N=shear/gn[:,None];lam2=-2*F*N[:,0]*(2*F*N[:,0]+gn)
   for i in range(len(r)):
    tn=float(stress[i]@N[i]);tk=float(stress[i]@np.array([-N[i,1],N[i,0]]));c=np.sqrt(max(lam2[i],0))/(2*F[i]*N[i,0]);passed=bool(lam2[i]>0 and tn<0 and abs(c*tk)<-tn)
    row=dict(k=k,eta=eta,y=float((g[i]+1)/2),r=float(r[i]),z=float(ip[2]),inner_radius=float(ri),distance_to_interface=float(min(r[i]-ri,6*ri-r[i])),lambda_squared=float(lam2[i]),target=stress[i].tolist(),target_dot_N=tn,cone_ratio=float(abs(c*tk/tn)) if tn else None,directional_pass=passed)
    rows.append(row)
    if passed and k==5.5:candidates.append((row,J[i]))
   print('screened',k,eta,flush=True)
 selected=max(candidates,key=lambda item:item[0]['distance_to_interface']) if candidates else None;trajectory=None;pulses=[]
 if selected:
  point,J=selected;tau0=.5*2**(-point['k']);x0=np.array([point['r'],0,point['z']])
  def similarity(s,x):
   rad=np.hypot(x[0],x[1]);co=coordinates(rad/np.sqrt(f.nu),x[2]/np.sqrt(f.nu),tau0-s,f.inner.h);ri=np.sqrt(2*f.nu*co['q']*3/64)
   return float((rad-ri)/(5*ri)),float(co['eta'])
  def velocity(s,x):return f.fields(x[None,:],tau0-s)[0][0]
  def inner(s,x):return similarity(s,x)[0]
  def outer(s,x):return 1-similarity(s,x)[0]
  def axial(s,x):return .48-abs(similarity(s,x)[1])
  for event in (inner,outer,axial):event.terminal=True;event.direction=-1
  sol=solve_ivp(velocity,(0,.1*tau0),x0,events=(inner,outer,axial),rtol=2e-9,atol=1e-12,max_step=.005*tau0,dense_output=True)
  if not sol.success:raise RuntimeError(sol.message)
  times=np.linspace(0,sol.t[-1],25);path=sol.sol(times).T;sy=[similarity(s,x) for s,x in zip(times,path)];reason=next((name for name,events in zip(('inner','outer','axial'),sol.t_events) if len(events)),'horizon')
  trajectory=dict(selected=point,event=reason,elapsed=float(sol.t[-1]),elapsed_over_tau=float(sol.t[-1]/tau0),times=times.tolist(),path=path.tolist(),y_eta=sy)
  for mode in (1,4,16):
   for az in (-1,0,1):
    for radial in (-1,0,1):
     nt=np.array([mode/point['r'],az*mode/point['r']]);n0=np.r_[radial*np.linalg.norm(nt),nt]
     pulse=affine_pulse(J,f.nu,n0,sol.t[-1]);pulse.update(mode=mode,axial_ratio=az,radial_ratio=radial);pulses.append(pulse)
  print('selected',point,'residence',trajectory['elapsed_over_tau'],'max gain',max(p['peak_gain'] for p in pulses),flush=True)
 report=dict(rows=rows,trajectory=trajectory,pulses=pulses,scope='Width6 corrected background;20 radial nodes at six scale/axial slices. Same approximate physical frozen directional screen; no paper theorem transfer. Selected late passing point maximizes distance to either radial interface. Actual particle traced to radial/axial probe exit or0.1*tau; passing the initial screen does not certify it along the trajectory. Affine amplitude model still freezes J at the starting point.',pde_validated=False,global_field_ready=False)
 out=ROOT/'wide_pulse_screen';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print('pass count',sum(r['directional_pass'] for r in rows),'/',len(rows),flush=True)
if __name__=='__main__':run()
