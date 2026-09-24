"""Actual finite-radius boundary traces, not a completed outer matching."""
from pathlib import Path
import json,numpy as np
from full_radial import *
from audit import ROOT,save

rows=[];saved={}
for ident in ['ST073-F','ST073-V']:
 f=FullRadialField.load(ROOT/'data'/f'{ident}.json')
 for k in [0,3,6]:
  tau=.5*2**(-k);es=np.linspace(-.5,.5,17);X=np.ones_like(es)/64;d=f.evaluate_similarity(X,es,tau)
  q=tau/(1-es**2);aa=[];vv=[];coefs=[]
  for e,qq in zip(es,q):
   co=f.coefficients(float(e),float(qq));coefs.append(co);bx=co[1,:,0];cx=co[2,:,0];x=1/64
   F=qq**(1+f.h)*np.polynomial.polynomial.polyval(x,bx);Fx=qq**(1+f.h)*np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(bx));Ux=qq**f.A*np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(cx))
   a=-2*x*Fx/F;v=a+2*x*(Ux/F)**2/a;aa.append(float(a));vv.append(float(v))
  data=ROOT/'data'/f'interface_{ident}_k{k}.npz';np.savez_compressed(data,eta=es,X=X,tau=tau,coefficient_jets=np.asarray(coefs),**d)
  rp=d['pressure_gradient'][:,0];v=d['velocity'];r=np.sqrt(2*f.nu*q/64)
  rows.append(dict(id=ident,k=k,point_count=len(es),a_min=min(aa),v_min=min(vv),v_at_midplane=vv[8],source_v_gt2_fraction=float(np.mean(np.array(vv)>2)),inward_pressure_fraction=float(np.mean(rp>0)),entrainment_mean=float(np.mean(-v[:,0]/(r/(2*tau)))),physical_pressure_gradient=d['pressure_gradient'].tolist(),scope='Actual full corrected field traces, not a leading five-moment match or strict-cone admission'))
  print(ident,k,rows[-1]['v_at_midplane'],rows[-1]['inward_pressure_fraction'],flush=True)
save(ROOT/'evidence/interface_audit.json',dict(rows=rows))
