"""Mass flux across the actual moving ST073 local cylinder in similarity coordinates."""
from pathlib import Path
import sys,json
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence'))
from full_radial import FullRadialField

def run():
 f=FullRadialField.load(ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V.json')
 nu=f.nu;D=.5-f.h;xc=1/64;em=.5;rows=[]
 for order in (12,24):
  nodes,weights=leggauss(order);eta=em*nodes;we=em*weights
  for k in (0,3,6):
   tau=.5*2**(-k);q=tau/(1-eta**2)
   r=np.sqrt(2*nu*xc*q);z=np.sqrt(nu)*q**D*eta
   rp=r*eta/(1-eta**2);zp=np.sqrt(nu)*q**D*(1+2*D*eta**2/(1-eta**2))
   u=np.asarray(f.evaluate_similarity(np.full(order,xc),eta,tau)['velocity'],float)
   side=float(np.sum(we*2*np.pi*r*(u[:,0]*zp-u[:,2]*rp)))
   side_boundary=float(np.sum(we*2*np.pi*r*((-r/(2*tau))*zp-(-D*z/tau)*rp)))
   volume=float(np.sum(we*np.pi*r*r*zp));caps=0.;cap_boundary=0.;cap_values=[]
   xx=xc*(nodes+1)/2;wx=xc*weights/2
   for sign in (-1,1):
    e=sign*em;qc=tau/(1-e*e);zc=np.sqrt(nu)*qc**D*e
    uc=np.asarray(f.evaluate_similarity(xx,np.full(order,e),tau)['velocity'],float)
    flux=float(sign*2*np.pi*nu*qc*np.sum(wx*uc[:,2]));caps+=flux;cap_values.append(flux)
    cap_boundary+=float(sign*np.pi*(2*nu*xc*qc)*(-D*zc/tau))
   boundary=side_boundary+cap_boundary;dv=-(1+D)*volume/tau
   rows.append(dict(order=order,k=k,tau=tau,volume=volume,side_fluid_flux=side,lower_cap_fluid_flux=cap_values[0],upper_cap_fluid_flux=cap_values[1],net_fluid_flux=side+caps,boundary_volume_rate=boundary,analytic_volume_rate=dv,relative_outward_flux=side+caps-boundary,reynolds_mass_defect=side+caps-boundary+dv))
 report=dict(rows=rows,geometry='X<=1/64, |eta|<=.5; side r(eta,tau), z(eta,tau), flat caps at eta=+/- .5',normal_convention='outward oriented area: side 2*pi*r*(z_eta,0,-r_eta); caps +/-ez',boundary_motion='at fixed similarity labels: b_r=-r/(2*tau), b_z=-(.5-h)*z/tau; physical t=1-tau',scope='Moving-domain incompressibility/Reynolds mass budget only. No stress matching, exterior construction, full-energy claim or scale-recursion acceptance.',global_field_ready=False)
 out=ROOT/'moving_interface';out.mkdir(exist_ok=True);(out/'mass_flux.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(rows,indent=2))
if __name__=='__main__':run()
