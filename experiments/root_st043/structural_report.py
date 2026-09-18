"""Independent-grid physical diagnostics; not a calibrated source field comparison."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import bootstrap
from spacetime import Family
from structure_audit import jets

def report(path,reference=None):
 f,raw=Family.load(path)
 if reference is not None:fb,rb=Family.load(reference)
 cases=[]
 for name,rgrid,zgrid,times in [
  ('prior_comparable_grid',np.linspace(.04,.2,11),np.r_[-np.linspace(.04,.2,11)[::-1],np.linspace(.04,.2,11)],np.array([.25,.3125,.4375,.5625,.6875,.75])),
  ('new_offgrid_check',np.linspace(.037,.197,13),np.r_[-np.linspace(.033,.193,13)[::-1],np.linspace(.033,.193,13)],np.r_[.25,np.sort(np.random.default_rng(9174391).uniform(.25,.75,13)),.75])]:
  R,Z=np.meshgrid(rgrid,zgrid,indexing='ij');rows=[];initial=None
  for t in times:
   tau=1-t;r=R*tau**.5;z=Z*tau**.495;v=jets(f,raw,r,z,t);u=v[...,:3];sv=u*np.array([tau**.5,tau**.505,tau**.505]);ax=np.sign(Z)*v[...,7]
   if initial is None:initial=sv.copy()
   row=dict(time=float(t),scaled_profile_drift=float(np.linalg.norm(sv-initial)/np.linalg.norm(initial)),inward_fraction=float(np.mean(u[...,0]<0)),positive_swirl_fraction=float(np.mean(u[...,1]>0)),bipolar_fraction=float(np.mean(Z*u[...,2]>0)),radial_pressure_inward_fraction=float(np.mean(v[...,6]>0)),axial_pressure_toward_fraction=float(np.mean(ax>0)),adverse_axial_pressure_rms=float(np.sqrt(np.mean(np.minimum(ax,0)**2))),aligned_axial_pressure_min=float(ax.min()))
   if reference is not None:
    parent=jets(fb,rb,r,z,t);parent_sv=parent[...,:3]*np.array([tau**.5,tau**.505,tau**.505]);row['velocity_difference_relative_to_parent']=float(np.linalg.norm(sv-parent_sv)/np.linalg.norm(parent_sv));row['axial_pressure_adverse_ratio_to_parent']=float(np.linalg.norm(np.minimum(ax,0))/max(np.linalg.norm(np.minimum(np.sign(Z)*parent[...,7],0)),1e-30))
   rows.append(row)
  cases.append(dict(grid=name,r_range=[float(rgrid[0]),float(rgrid[-1])],z_abs_range=[float(np.min(abs(zgrid))),float(np.max(abs(zgrid)))],points_per_time=int(R.size),rows=rows))
 mid=[]
 for t in [.25,.371,.5,.633,.75]:
  r=np.linspace(.04,.6,65);v=jets(f,raw,r,0,t);mid.append(dict(time=t,min_uz=float(np.min(v[:,2])),max_uz=float(np.max(v[:,2])),max_abs_radial_axial_shear=float(np.max(abs(v[:,8])))))
 return dict(candidate_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),cases=cases,midplane=mid,scope='Autonomous finite-window grids; new off-grid samples not used in fitting. Preserved flow signs, drift and pressure force are separate diagnostics, not a full source identity or continuum certificate.',pde_validated=False,source_correspondence_verified=False)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--reference');p.add_argument('--out',required=True);a=p.parse_args();r=report(a.candidate,a.reference);Path(a.out).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({q['grid']:q['rows'][-1] for q in r['cases']},indent=2))
