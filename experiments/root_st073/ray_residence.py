"""Actual background particle residence near the candidate pulse collar."""
import json
import numpy as np
from scipy.integrate import solve_ivp
from angular_collocation import CachedBase
from poloidal_collocation import PoloidalModes
from joined_field import ROOT,coordinates

def run():
 old=json.loads((ROOT/'affine_pulse/report.json').read_text());base=CachedBase();a=json.loads((ROOT/'poloidal_collocation/report.json').read_text())['amplitudes'];f=PoloidalModes(base,np.array(a).ravel());rows=[]
 for record in old['rows']:
  point=record['point'];tau0=.5*2**(-point['k']);ip=f.inner.from_similarity([3/64],[point['eta']],tau0)[0];x0=np.array([point['r'],0,ip[2]])
  def similarity(s,x):
   r=np.hypot(x[0],x[1]);co=coordinates(r/np.sqrt(f.nu),x[2]/np.sqrt(f.nu),tau0-s,f.inner.h);ri=np.sqrt(2*f.nu*co['q']*3/64)
   return float(r/ri-1),float(co['eta'])
  def velocity(s,x):return f.fields(x[None,:],tau0-s)[0][0]
  def inner(s,x):return similarity(s,x)[0]
  def outer_probe(s,x):return .02-similarity(s,x)[0]
  inner.terminal=True;inner.direction=-1;outer_probe.terminal=True;outer_probe.direction=-1
  sol=solve_ivp(velocity,(0,.1*tau0),x0,events=(inner,outer_probe),rtol=2e-9,atol=1e-12,max_step=.001*tau0,dense_output=True)
  if not sol.success:raise RuntimeError(sol.message)
  times=np.linspace(0,sol.t[-1],25);path=sol.sol(times).T;ys=[similarity(s,x)[0] for s,x in zip(times,path)]
  event='inner_interface' if len(sol.t_events[0]) else 'outer_probe_y=.02' if len(sol.t_events[1]) else 'horizon'
  row=dict(eta0=point['eta'],k=point['k'],event=event,elapsed=float(sol.t[-1]),elapsed_over_tau=float(sol.t[-1]/tau0),times=times.tolist(),path=path.tolist(),relative_radial_y=ys,final_eta=similarity(sol.t[-1],sol.y[:,-1])[1]);rows.append(row);print({k:row[k] for k in ('eta0','event','elapsed','elapsed_over_tau','final_eta')},flush=True)
 report=dict(rows=rows,scope='Lagrangian particles in the actual time-dependent poloidal/pressure candidate. Moving inner-interface event evaluated with q(z,tau). Outer y=.02 is an explicit probe boundary, NOT an established cone boundary. Particle transport only, not a wavepacket or full residual solve.',pde_validated=False)
 out=ROOT/'affine_pulse';(out/'residence.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
