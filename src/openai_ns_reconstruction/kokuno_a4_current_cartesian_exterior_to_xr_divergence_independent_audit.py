"""Independent public-velocity Cartesian divergence audit for A1 #980 through X_R."""
from __future__ import annotations
import argparse, copy, hashlib, json, math, tempfile
from pathlib import Path
from typing import Any, Mapping
import numpy as np
from .kokuno_pa16_current_cartesian_exterior_to_xr import KokunoPA16CurrentCartesianExteriorToXR

SCHEMA='kokuno-a4-current-cartesian-exterior-to-xr-divergence-audit-v1'
UPSTREAM_PR=980
UPSTREAM_HEAD='d3c971f2c62e272333e124e532212d23cca4908d'
UPSTREAM_SOURCE_BLOB='2ba28636a56b252fa485719e7b3e8da754d7e588'
SEED=9173691
STEPS=(0.02,0.01,0.005)
TIMES=(0.31,0.47,0.63,0.71)
ETA=(-0.60,0.60)
ZONES=(('near_seam_exterior',-4.80,-4.10),('middle_exterior',-3.45,-1.85),('near_XR_exterior',-1.35,-0.35))
N_PER_ZONE_TIME=6
SEAM_OFFSETS=(-0.06,-0.02,0.02,0.06)
DIV_GATE=1e-5
STABILITY_FACTOR=1.25
STABILITY_FLOOR=2e-8
SPEED_FLOOR=1e-10
MUT_EPS=1e-3
MUT_DETECT=5e-4
ORDER_TOL=2e-13
TRUTH={
 'current_cartesian_leading_through_XR_consumed':True,
 'candidate_save_reload_required':True,
 'independent_cartesian_fd_operator_used':True,
 'exterior_through_XR_divergence_scoped_assessed':True,
 'seam_near_divergence_scoped_assessed':True,
 'axis_regression_scoped_assessed':True,
 'leading_only_ns_residual_assessed':False,
 'leading_plus_oscillatory_ns_residual_assessed':False,
 'after_correction_ns_residual_assessed':False,
 'post_XR_global_velocity_materialized':False,
 'matched_global_pressure_materialized':False,
 'restricted_forcing_materialized':False,
 'canonical_24_48_96_volume_admission_assessed':False,
 'same_protocol_comparable_to_st006':False,
 'pde_validated':False,
}

def _canon(x:Mapping[str,Any])->str:return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False)
def _lh(rng,n):
 a=(np.arange(n)+rng.random(n))/n;rng.shuffle(a);return a

def _xyz(field,X,eta,theta,t):
 X,eta,theta,t=np.broadcast_arrays(*[np.asarray(v,float) for v in (X,eta,theta,t)])
 q=(1-t)/(1-eta*eta);r=np.sqrt(2*q*X);z=q**float(field.D)*eta
 return np.stack((r*np.cos(theta),r*np.sin(theta),z),axis=-1)

def _jac_volume(field,eta,t):
 eta,t=np.broadcast_arrays(np.asarray(eta,float),np.asarray(t,float));D=float(field.D)
 q=(1-t)/(1-eta*eta);return q**(D+1)*(1+2*D*eta*eta/(1-eta*eta))

def _velocity(field,p,t):
 p=np.asarray(p,float);t=np.asarray(t,float)
 u=np.asarray(field.velocity(p[:,0],p[:,1],p[:,2],t),float)
 if u.shape!=(len(p),3) or not np.all(np.isfinite(u)):raise RuntimeError('bad public velocity')
 return u

def fd2_jacobian(field,p,t,h):
 p=np.asarray(p,float);t=np.asarray(t,float);J=np.empty((len(p),3,3))
 for k in range(3):
  pp=p.copy();pm=p.copy();pp[:,k]+=h;pm[:,k]-=h
  J[:,:,k]=(_velocity(field,pp,t)-_velocity(field,pm,t))/(2*h)
 return J

def make_probes(field):
 rng=np.random.default_rng(SEED);ps=[];ts=[];ws=[];zi=[];ti=[]
 for it,t in enumerate(TIMES):
  for iz,(_,s0,s1) in enumerate(ZONES):
   s=s0+(s1-s0)*_lh(rng,N_PER_ZONE_TIME);X=float(field.X_R)*np.exp(s)
   eta=ETA[0]+(ETA[1]-ETA[0])*_lh(rng,N_PER_ZONE_TIME);th=2*math.pi*_lh(rng,N_PER_ZONE_TIME)
   tv=np.full(N_PER_ZONE_TIME,t);p=_xyz(field,X,eta,th,tv)
   w=(s1-s0)*(ETA[1]-ETA[0])*2*math.pi*X*_jac_volume(field,eta,tv)/N_PER_ZONE_TIME
   ps.append(p);ts.append(tv);ws.append(w);zi.append(np.full(N_PER_ZONE_TIME,iz));ti.append(np.full(N_PER_ZONE_TIME,it))
 seam=[];seamt=[];labels=[]
 for it,t in enumerate(TIMES):
  for j,o in enumerate(SEAM_OFFSETS):
   X=float(field.X_h)*math.exp(o);eta=-.45+.30*((it+j)%4);th=.37+.61*it+.29*j
   seam.append(_xyz(field,[X],[eta],[th],[t])[0]);seamt.append(t);labels.append('inner' if o<0 else 'exterior')
 axX=np.array([0.,0.,1e-10,1e-8,1e-6]);axe=np.array([0.,.45,-.45,.25,0.]);axth=np.array([0.,1.,.7,1.3,2.1]);axt=np.array([.50,.47,.63,.31,.71])
 return {'p':np.concatenate(ps),'t':np.concatenate(ts),'w':np.concatenate(ws),'zi':np.concatenate(zi).astype(int),'ti':np.concatenate(ti).astype(int),'seam':np.asarray(seam),'seamt':np.asarray(seamt),'labels':labels,'axis':_xyz(field,axX,axe,axth,axt),'axist':axt}

def _point_stats(field,p,t,h):
 J=fd2_jacobian(field,p,t,h);d=np.trace(J,axis1=1,axis2=2);f=np.linalg.norm(J,axis=(1,2));i=int(np.argmax(abs(d)))
 return {'sampled_max':float(np.max(abs(d))),'rms':float(np.sqrt(np.mean(d*d))),'normalized_sampled_max':float(np.max(abs(d)/np.maximum(1.,f))),'worst':{'point':[float(x) for x in p[i]],'time':float(t[i]),'divergence':float(d[i])}}

def metrics(field,P,h):
 J=fd2_jacobian(field,P['p'],P['t'],h);d=np.trace(J,axis1=1,axis2=2);f=np.linalg.norm(J,axis=(1,2));vol=float(np.sum(P['w']));n=abs(d)/np.maximum(1.,f);speed=np.linalg.norm(_velocity(field,P['p'],P['t']),axis=1);i=int(np.argmax(abs(d)))
 zones=[]
 for iz,(name,_,_) in enumerate(ZONES):
  m=P['zi']==iz;w=P['w'][m];dd=d[m];v=float(np.sum(w));zones.append({'zone':name,'sampled_max':float(np.max(abs(dd))),'weighted_rms':float(np.sqrt(np.sum(w*dd*dd)/v)),'volume_l2_estimate':float(np.sqrt(np.sum(w*dd*dd)))})
 return {'step':h,'sampled_max':float(np.max(abs(d))),'pooled_weighted_rms':float(np.sqrt(np.sum(P['w']*d*d)/vol)),'pooled_volume_l2_estimate':float(np.sqrt(np.sum(P['w']*d*d))),'normalized_sampled_max':float(np.max(n)),'normalized_weighted_rms':float(np.sqrt(np.sum(P['w']*n*n)/vol)),'speed_rms':float(np.sqrt(np.mean(speed*speed))),'seam':_point_stats(field,P['seam'],P['seamt'],h),'axis':_point_stats(field,P['axis'],P['axist'],h),'per_zone':zones,'worst':{'point':[float(x) for x in P['p'][i]],'time':float(P['t'][i]),'divergence':float(d[i])}}

class Mut:
 def __init__(self,b):self.b=b;self.D=b.D;self.X_h=b.X_h;self.X_R=b.X_R
 def velocity(self,x,y,z,t):
  u=np.asarray(self.b.velocity(x,y,z,t),float).copy();u[...,0]+=MUT_EPS*np.asarray(x,float);return u

def _stable(m,f):return f<=STABILITY_FACTOR*m+STABILITY_FLOOR
def _param_mut(field):
 p=copy.deepcopy(field.configuration());p['current']['eta_fd_step']*=1.001
 return KokunoPA16CurrentCartesianExteriorToXR.from_configuration(p).semantic_sha256!=field.semantic_sha256

def _post_xr(field):
 p=_xyz(field,[float(field.X_R)*1.001],[.2],[.7],[.5])[0]
 try:field.velocity(*p,.5)
 except ValueError:return True
 return False

def _pass(r):
 m,f=r['resolutions'][1:];c=r['checks']
 return bool(f['sampled_max']<=DIV_GATE and f['pooled_weighted_rms']<=DIV_GATE and f['seam']['sampled_max']<=DIV_GATE and f['axis']['sampled_max']<=DIV_GATE and f['speed_rms']>=SPEED_FLOOR and _stable(m['sampled_max'],f['sampled_max']) and _stable(m['pooled_weighted_rms'],f['pooled_weighted_rms']) and _stable(m['seam']['sampled_max'],f['seam']['sampled_max']) and _stable(m['axis']['sampled_max'],f['axis']['sampled_max']) and all(c.values()))

def materialize_receipt():
 orig=KokunoPA16CurrentCartesianExteriorToXR()
 with tempfile.TemporaryDirectory() as d:
  p=Path(d)/'candidate.json';orig.save_configuration(p);field=KokunoPA16CurrentCartesianExteriorToXR.load_configuration(p)
 P=make_probes(field);res=[metrics(field,P,h) for h in STEPS]
 rev={k:v for k,v in P.items()};idx=np.arange(len(P['p'])-1,-1,-1)
 for k in ('p','t','w','zi','ti'):rev[k]=P[k][idx]
 rr=metrics(field,rev,STEPS[-1]);mut=metrics(Mut(field),P,STEPS[-1])
 order=abs(rr['sampled_max']-res[-1]['sampled_max'])<=ORDER_TOL and abs(rr['pooled_weighted_rms']-res[-1]['pooled_weighted_rms'])<=ORDER_TOL
 r={'schema':SCHEMA,'upstream':{'pr':UPSTREAM_PR,'head':UPSTREAM_HEAD,'source_blob':UPSTREAM_SOURCE_BLOB,'semantic_sha256':field.semantic_sha256,'X_h':float(field.X_h),'X_R':float(field.X_R)},'protocol':{'seed':SEED,'steps':list(STEPS),'times':list(TIMES),'eta_interval':list(ETA),'zones':[list(x) for x in ZONES],'points_per_zone_time':N_PER_ZONE_TIME,'seam_offsets':list(SEAM_OFFSETS),'integration_probe_count':len(P['p']),'seam_probe_count':len(P['seam']),'axis_probe_count':len(P['axis']),'divergence_gate':DIV_GATE,'stability_factor':STABILITY_FACTOR,'stability_floor':STABILITY_FLOOR,'speed_floor':SPEED_FLOOR,'mutation_epsilon':MUT_EPS,'mutation_detection_floor':MUT_DETECT,'operator':'centered Cartesian FD2 on save/reloaded public velocity only','volume_note':'scoped physical-volume exterior estimate, not canonical 24/48/96 whole-domain admission'},'resolutions':res,'mutation':{'kind':'u_x += 1e-3*x','finest_sampled_max':mut['sampled_max'],'finest_weighted_rms':mut['pooled_weighted_rms']},'checks':{'save_reload_semantic_identity':field.semantic_sha256==orig.semantic_sha256,'config_parameter_mutation_changes_semantic_identity':_param_mut(field),'velocity_mutation_detected':mut['sampled_max']>=MUT_DETECT,'offgrid_order_invariant':order,'post_XR_fail_closed':_post_xr(field)},'truth_boundary':copy.deepcopy(TRUTH)}
 r['audit_pass']=_pass(r);r['receipt_sha256']=hashlib.sha256(_canon(r).encode()).hexdigest();return r

def enforce_receipt(r):
 if r.get('schema')!=SCHEMA:raise ValueError('schema drift')
 p=r.get('protocol',{});exp={'seed':SEED,'steps':list(STEPS),'times':list(TIMES),'eta_interval':list(ETA),'zones':[list(x) for x in ZONES],'points_per_zone_time':N_PER_ZONE_TIME,'seam_offsets':list(SEAM_OFFSETS),'divergence_gate':DIV_GATE,'stability_factor':STABILITY_FACTOR,'stability_floor':STABILITY_FLOOR,'speed_floor':SPEED_FLOOR,'mutation_epsilon':MUT_EPS,'mutation_detection_floor':MUT_DETECT}
 for k,v in exp.items():
  if p.get(k)!=v:raise ValueError(f'protocol drift: {k}')
 if r.get('upstream',{}).get('head')!=UPSTREAM_HEAD or r.get('upstream',{}).get('source_blob')!=UPSTREAM_SOURCE_BLOB:raise ValueError('upstream identity drift')
 for k,v in TRUTH.items():
  if r.get('truth_boundary',{}).get(k) is not v:raise ValueError(f'truth drift: {k}')
 q=copy.deepcopy(dict(r));obs=q.pop('receipt_sha256',None)
 if obs!=hashlib.sha256(_canon(q).encode()).hexdigest():raise ValueError('receipt checksum mismatch')
 if bool(r.get('audit_pass'))!=_pass(r):raise ValueError('audit_pass laundering')
 if not _pass(r):
  f=r['resolutions'][-1];raise RuntimeError(f"XR divergence audit failed: max={f['sampled_max']:.6e}, rms={f['pooled_weighted_rms']:.6e}, seam={f['seam']['sampled_max']:.6e}, axis={f['axis']['sampled_max']:.6e}")

def _main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);ap.add_argument('--verify-receipt',type=Path);a=ap.parse_args()
 if a.verify_receipt:
  enforce_receipt(json.loads(a.verify_receipt.read_text()));print('A4 XR divergence receipt passes frozen gates');return
 r=materialize_receipt();s=json.dumps(r,indent=2,sort_keys=True,allow_nan=False)
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s+'\n')
 print(s)
if __name__=='__main__':_main()
