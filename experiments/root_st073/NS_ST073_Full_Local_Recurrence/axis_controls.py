"""Preregistered leading-only controls. Never selected as an accepted field."""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'upstream'))
from normalized_core import Parameters,NormalizedCore
from local_field import LocalField
from heat_exterior import profile,tail_moments

def run(sigma,pressure):
 tag=f'sigma{sigma:g}_P{pressure:g}';out=ROOT/'evidence'/f'control_{tag}.json'
 if out.exists():return json.loads(out.read_text())
 p=Parameters(lam=256,sigma=sigma,pressure=pressure,g_peak=.02,order=26,jet_degree=32)
 c=NormalizedCore(p);Hprime=c.D+4-12*c.eta0**2-2*p.bias*c.eta0
 w=p.sigma/np.sqrt(p.lam*Hprime)
 es=np.unique(np.r_[np.linspace(-.5,.5,31),c.eta0+np.array([-3,-2,-1,0,1,2,3])*w]);es=es[abs(es)<=.5]
 xx,ee=np.meshgrid([.003,.009,.015625],es);t=time.monotonic()
 fld=LocalField(c);pts=fld.from_similarity(xx.ravel(),ee.ravel(),.5);v=fld.evaluate(pts,.5)
 seam=c.profiles(np.full(len(es),1/64),es);F=seam['F'][0,0].astype(np.longdouble);fx=seam['F'][1,0].astype(np.longdouble);ux=seam['U'][1,0].astype(np.longdouble)
 aa=-2*(1/64)*fx/F;vv=aa+2*(1/64)*(ux/F)**2/aa
 fc=profile(4.,es,1.)['F'];cmax=np.min(F/fc);cap=F**2*(4-1/64)+cmax**2*tail_moments(4.,es,1.,n=64)[0];needed=-seam['P'][0,0]
 idx=int(np.argmin(abs(es-c.eta0)));res=dict(sigma=sigma,pressure=pressure,g_peak=.02,Lambda=256,eta0=c.eta0,curvature_width=w,axis_curvature=-p.lam*(1-2*p.h*c.eta0**2)*Hprime/sigma**2,sampled_full_max_k0=float(np.linalg.norm(v['residual'],axis=-1).max()),sampled_point_count=len(pts),seam_positive_a_fraction=float(np.mean(aa>0)),seam_v_gt_2_fraction=float(np.mean(vv>2)),pressure_bound_violation_fraction=float(np.mean(needed>cap)),peak_needed=float(needed[idx]),peak_cap=float(cap[idx]),seconds=time.monotonic()-t,scope='Changed autonomous axis regularization and pressure. Diagnostic sample maxima and necessary inequalities only; not comparable peak certificates or corrected full solutions',pde_validated=False)
 out.write_text(json.dumps(res,indent=2)+'\n');c.save(ROOT/'data'/f'control_{tag}.json');print(res,flush=True);return res
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--pressure',type=float,required=True);x=a.parse_args()
 for s in [.008,.04,.12,.25]:run(s,x.pressure)
