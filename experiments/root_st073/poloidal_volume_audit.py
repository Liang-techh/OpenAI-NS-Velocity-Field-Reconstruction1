"""Physical-volume quadrature on the finite radial transition slab."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from poloidal_collocation import PoloidalModes
from angular_collocation import CachedBase
from joined_field import ROOT,independent_fd

def run():
 fit=json.loads((ROOT/'poloidal_collocation/report.json').read_text());base=CachedBase()
 fields=[('baseline',PoloidalModes(base,np.zeros(12))),('poloidal_pressure_candidate',PoloidalModes(base,np.array(fit['amplitudes']).ravel()))]
 rows=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k);gy,wy=leggauss(12);ge,we=leggauss(8)
  eta=np.repeat(.4*ge,12);y=np.tile((gy+1)/2,8);q=tau/(1-eta**2)
  ri=np.sqrt(2*base.nu*q*3/64);r=ri*(1+y);z=np.sqrt(base.nu)*q**base.inner.D*eta
  ze=np.sqrt(base.nu)*q**base.inner.D*(1+2*base.inner.D*eta**2/(1-eta**2))
  weights=np.repeat(.4*we,12)*np.tile(wy/2,8)*2*np.pi*r*ri*ze
  pts=np.column_stack((r,np.zeros_like(r),z))
  for label,f in fields:
   res,div=independent_fd(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau)
   row=dict(k=k,field=label,volume=float(weights.sum()),momentum_sampled_max=float(np.max(np.linalg.norm(res,axis=1))),momentum_volume_L2=float(np.sqrt(np.sum(weights*np.sum(res**2,axis=1)))) ,divergence_sampled_max=float(np.max(np.abs(div))))
   rows.append(row);print(json.dumps(row),flush=True)
 report=dict(rows=rows,domain='At each specified tau, eta in [-.4,.4], ri(z,t)<r<2ri(z,t), all azimuths. X_inner=3/64. Physical-volume weight=2pi r ri z_eta dy deta.',quadrature=dict(radial=12,axial=8),limitations='Finite quadrature on a subdomain; no continuum supremum, no convergence-certified integral, excludes core, axial ends and exterior. Not a global gate pass.',pde_validated=False)
 (ROOT/'poloidal_collocation/volume_audit.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
