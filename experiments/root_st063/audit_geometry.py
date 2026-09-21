"""Independent fixed-cylinder and connected-threshold morphology diagnostics.
No candidate fitting. No source-image alignment or continuum claim.
"""
import sys,json,argparse,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import label
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'experiments/root_st057'))
from tensor_model import Model

def fields_grid(m,r,z,t):
 R,Z=np.meshgrid(r,z,indexing='ij');g=m.grid(R.ravel()**2,Z.ravel(),np.atleast_1d(t));co,_=m.physical(np.zeros(m.dim));v=g.values(co)
 shape=(len(r),len(z),len(g.t));out={k:a.reshape(shape) for k,a in v.items()}
 out['ur']=R[...,None]*out['A'];out['utheta']=R[...,None]*out['B'];out['uz']=out['C'];out['omega_z']=2*(out['B']+R[...,None]**2*out['Bs'])
 out['omega_r']=-R[...,None]*out['Bz'];out['omega_theta']=R[...,None]*(out['Az']-2*out['Cs'])
 out['omega_mag']=np.sqrt(out['omega_z']**2+out['omega_r']**2+out['omega_theta']**2)
 out['speed']=np.sqrt(out['ur']**2+out['utheta']**2+out['uz']**2)
 out['residual_mag']=np.linalg.norm(g.residual(v),axis=-1).reshape(shape)
 return out

def audit(candidate,parent,out):
 out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
 m=Model(candidate);mp=Model(parent);times=[.25,.33713,.5,.64137,.75]
 report={'candidate_sha256':hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),'parent_sha256':hashlib.sha256(Path(parent).read_bytes()).hexdigest(),'optimizer_executed':False,'pde_validated':False,'scope':'Finite independent morphology; omega_z moments centered in z and single-transverse variance r^2/2, physical axes; not a detected intrinsic core or source-image match','moments':[],'connected_levels':[],'axis_profiles':[]}
 for nr,nz in [(40,64),(64,96)]:
  a,wa=leggauss(nr);b,wb=leggauss(nz)
  for rmax,zmax in [(.35,.60),(.60,.90)]:
   s=.5*rmax*rmax*(a+1);z=zmax*b;w=(np.pi*.5*rmax*rmax*zmax*np.outer(wa,wb))[:,:,None]
   vals=fields_grid(m,np.sqrt(s),z,times);om=vals['omega_z']**2;M=(w*om).sum((0,1));sz=(w*om*z[None,:,None]).sum((0,1))/M
   varz=(w*om*(z[None,:,None]-sz)**2).sum((0,1))/M;varr=(w*om*s[:,None,None]/2).sum((0,1))/M
   for i,t in enumerate(times):
    report['moments'].append(dict(order=[nr,nz],radius=rmax,z_half=zmax,time=t,chi=float(np.sqrt(varz[i]/varr[i])),sigma_perp=float(np.sqrt(varr[i])),sigma_z=float(np.sqrt(varz[i])),z_centroid=float(sz[i]),omega_z_squared_integral=float(M[i])))
 # Common absolute thresholds; connected COMPONENT CONTAINING AXIS MIDPLANE,
 # rather than merging distant lobes into a misleading single aspect ratio.
 r=np.linspace(0,.65,131);z=np.linspace(-.9,.9,181);v=fields_grid(m,r,z,times)
 for i,t in enumerate(times):
  for scalar,levels in [('B',[.05,.10,.15]),('omega_z',[.15,.25,.35])]:
   a=v[scalar][:,:,i]
   for lev in levels:
    labs,n=label(a>=lev);center=labs[0,len(z)//2];touch=False
    if center:
     reg=labs==center;ri,zi=np.where(reg);height=z[zi.max()]-z[zi.min()];rad=r[ri.max()];touch=bool(reg[-1].any() or reg[:,[0,-1]].any())
     span=float(height);radius=float(rad)
    else:span=0.;radius=0.
    report['connected_levels'].append(dict(time=t,scalar=scalar,level=lev,axis_midplane_component_exists=bool(center),axial_span=span,max_radius=radius,touches_observation_boundary=touch))
  report['axis_profiles'].append(dict(time=t,z=z.tolist(),B_axis=v['B'][0,:,i].tolist(),omega_z_axis=v['omega_z'][0,:,i].tolist()))
 out.write_text(json.dumps(report,indent=2)+'\n');print('MORPH',out,[(x['time'],x['chi']) for x in report['moments'] if x['order']==[64,96] and x['radius']==.35],flush=True)
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--parent',required=True);p.add_argument('--out',required=True);a=p.parse_args();audit(a.candidate,a.parent,a.out)
