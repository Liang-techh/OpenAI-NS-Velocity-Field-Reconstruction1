"""Measured structure of a frozen field. Probe domains are autonomous diagnostics.
No plot-window or source schematic is treated as a calibrated target field.
"""
from __future__ import annotations
import bootstrap
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import label
from spacetime import Family

def jets(f,raw,r,z,t):
 r,z,t=np.broadcast_arrays(r,z,t);shape=r.shape;r,z,t=r.ravel(),z.ravel(),t.ravel();out=[]
 aa,bb,qq,_,_=f.coefficients(raw)
 for start in range(0,len(r),512):
  sl=slice(start,start+512);rr,zz,tt=r[sl],z[sl],t[sl];D=f.bundle(rr**2,zz,tt)
  v={k:D[k]@(bb if k.startswith('B') else qq if k.startswith('Q') else aa) for k in ['A','B','C','Bz','Bs','Az','Cs','Qs','Qz']}
  u=np.column_stack((rr*v['A'],rr*v['B'],v['C']))
  # axisymmetric vorticity: [-r B_z, r A_z-2r C_s, 2B+2r^2 B_s]
  w=np.column_stack((-rr*v['Bz'],rr*(v['Az']-2*v['Cs']),2*v['B']+2*rr**2*v['Bs']))
  out.append(np.column_stack((u,w,2*rr*v['Qs'],v['Qz'],2*rr*v['Cs'])))
 return np.concatenate(out).reshape(shape+(9,))

def weighted_quantile(x,w,p):
 i=np.argsort(x);return float(np.interp(p*np.sum(w),np.cumsum(w[i]),x[i]))

def audit(f,raw,grid=65):
 times=np.array([.25,.3125,.4375,.5625,.6875,.75]);sections=[];cores=[]
 for upper in [.2,.4,.7]:
  R,Z=np.meshgrid(np.linspace(.04,upper,11),np.r_[-np.linspace(.04,upper,11)[::-1],np.linspace(.04,upper,11)],indexing='ij');vals=[];row=[]
  for t in times:
   tau=1-t;v=jets(f,raw,R*np.sqrt(tau),Z*tau**.495,t);u=v[...,:3];scaled=u*np.array([np.sqrt(tau),tau**.505,tau**.505]);vals.append(scaled)
   row.append(dict(time=float(t),inward_fraction=float(np.mean(u[...,0]<0)),positive_swirl_fraction=float(np.mean(u[...,1]>0)),bipolar_fraction=float(np.mean(Z*u[...,2]>0)),inward_pressure_fraction=float(np.mean(v[...,6]>0)),axial_pressure_toward_plane_fraction=float(np.mean(Z*v[...,7]>0)),speed_rms=float(np.sqrt(np.mean(np.sum(u*u,axis=-1))))))
  ref=vals[0];norm=np.linalg.norm(ref)
  for v,d in zip(vals,row):d['scaled_profile_L2_drift']=float(np.linalg.norm(v-ref)/norm)
  cores.append(dict(probe='autonomous scaled grid, not a source-defined core boundary',R_range=[.04,upper],Z_abs_range=[.04,upper],points=R.size,rows=row))
 r=np.linspace(0,2,grid);z=np.linspace(-2,2,2*grid-1);R,Z=np.meshgrid(r,z,indexing='ij');weights=R.copy();weights[[0,-1]]*=.5;weights[:,[0,-1]]*=.5
 for t in [.25,.5,.75]:
  v=jets(f,raw,R,Z,t);u=v[...,:3];w=v[...,3:6];u2=np.sum(u*u,axis=-1);w2=np.sum(w*w,axis=-1);den=np.sum(weights*w2);k=np.unravel_index(np.argmax(w2),w2.shape);ku=np.unravel_index(np.argmax(u2),u2.shape)
  row=dict(time=t,grid=[grid,2*grid-1],speed_max=float(np.sqrt(u2[ku])),speed_peak_rz=[float(R[ku]),float(Z[ku])],vorticity_max=float(np.sqrt(w2[k])),vorticity_peak_rz=[float(R[k]),float(Z[k])],enstrophy_r_rms=float(np.sqrt(np.sum(weights*w2*R**2)/den)),enstrophy_z_rms=float(np.sqrt(np.sum(weights*w2*Z**2)/den)),enstrophy_in_collar=float(np.sum(weights*w2*((R>1.5)|(abs(Z)>1.5)))/den),enstrophy_in_r_lt_half=float(np.sum(weights*w2*(R<.5))/den),energy_in_r_lt_half=float(np.sum(weights*u2*(R<.5))/np.sum(weights*u2)))
  for a in [.25,.5]:
   masks,n=label(w2>a*a*np.max(w2));masses=[float(np.sum(weights*w2*(masks==i))/den) for i in range(1,n+1)]
   row[f'vorticity_{a}_superlevel_components_over_1pct']=sum(m>.01 for m in masses)
  sections.append(row)
 rmid=np.linspace(.04,.6,33);middle=[]
 for t in [.25,.5,.75]:
  v=jets(f,raw,rmid,0,t)
  middle.append(dict(time=t,uz_max_abs=float(np.max(abs(v[:,2]))),uz_min=float(np.min(v[:,2])),radial_axial_shear_max_abs=float(np.max(abs(v[:,8])))))
 return dict(core_grids=cores,whole_support=sections,middle_plane=middle,scope='Sampled diagnostics and exact parity of ansatz must be distinguished; high-vorticity support is not necessarily the source core; no calibrated image similarity score.',pde_validated=False,source_correspondence_verified=False)

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',required=True);p.add_argument('--grid',type=int,default=65);a=p.parse_args();f,raw=Family.load(a.candidate);d=audit(f,raw,a.grid);d['candidate_sha256']=hashlib.sha256(Path(a.candidate).read_bytes()).hexdigest();Path(a.out).write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'whole_support':d['whole_support'],'core_final':[x['rows'][-1] for x in d['core_grids']],'middle_plane':d['middle_plane']},indent=2))
