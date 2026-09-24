"""Trade off full momentum against interior pulse-growth and stress directions."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legvander,legint,legval
from scipy.optimize import least_squares
from wide_modes import CachedWidth,WideJointModes,sample
from affine_momentum import jets,momentum,combine
from joined_field import ROOT,independent_fd

def prefix_matrix(g):
 n=len(g);v=legvander(g,n-1)
 a=np.column_stack([legval(g,legint([0]*j+[1]))-legval(-1,legint([0]*j+[1])) for j in range(n)])
 return a@np.linalg.inv(v)

def stress_and_cone(u,J,res,r,width,Q,nu):
 stress=np.column_stack((-width/2*(Q@(r*r*res[:,1]))/(r*r),-width/2*(Q@(r*res[:,2]))/r))
 F=u[:,1]/r;shear=np.column_stack((J[:,1,0]-F,J[:,2,0]));sn=np.linalg.norm(shear,axis=1);N=shear/np.maximum(sn[:,None],1e-15)
 lam2=-2*F*N[:,0]*(2*F*N[:,0]+sn);lam=np.sqrt(np.maximum(lam2,0));tn=np.sum(stress*N,axis=1);tk=stress[:,1]*N[:,0]-stress[:,0]*N[:,1]
 q=2*F*N[:,0];direction=(tn<0)&(lam2>0)&(lam*np.abs(tk)<np.abs(q)*(-tn));growth=lam>nu/r**2
 return stress,direction,growth,lam,tn,tk,q,sn,F

def run():
 base=CachedWidth(6.);prior=np.array(json.loads((ROOT/'wide_collocation/report.json').read_text())['amplitudes']);training=[];nodes,weights=leggauss(16);Q=prefix_matrix(nodes);interior=np.where(((nodes+1)/2>=.12)&((nodes+1)/2<=.88))[0]
 for k in (1.,4.,5.5):
  tau=.5*2**(-k)
  for eta in (-.3,0,.3):
   ip=base.inner.from_similarity([3/64],[eta],tau)[0];ri=ip[0];width=5*ri;r=ri+width*(nodes+1)/2;pts=np.column_stack((r,np.zeros(16),np.full(16,ip[2])))
   args=(pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau);f0=WideJointModes(base,np.zeros(24));jb=jets(f0,*args);modes=[]
   for j in range(24):
    a=np.zeros(24);a[j]=1;ju=jets(WideJointModes(base,a),*args);modes.append(tuple(x-y for x,y in zip(ju,jb)))
   modes=tuple(np.stack([m[i] for m in modes]) for i in range(3))
   up,Jp,Lp=combine(jb,modes,prior);Rp=momentum((up,Jp,Lp));sp=np.max(np.linalg.norm(Rp,axis=1));Tp=stress_and_cone(up,Jp,Rp,r,width,Q,base.nu)[0];ts=max(1e-6,float(np.max(np.linalg.norm(Tp,axis=1))))
   wr=weights*r*r;wr/=wr.sum();wz=weights*r;wz/=wz.sum();training.append(dict(k=k,eta=eta,base=jb,modes=modes,r=r,width=width,res_scale=sp,stress_scale=ts,wr=wr,wz=wz))
   print('assembled',k,eta,flush=True)
 def objective(a):
  values=[]
  for item in training:
   u,J,L=combine(item['base'],item['modes'],a);res=momentum((u,J,L));s=item['res_scale'];R=res/s
   parts=[R.ravel(),np.array([3*item['wr']@R[:,1],3*item['wz']@R[:,2]])]
   if item['k']>=4:
    _,direction,growth,lam,tn,tk,q,sn,F=stress_and_cone(u,J,res,item['r'],item['width'],Q,base.nu)
    ix=interior;damp=base.nu/item['r'][ix]**2;ts=item['stress_scale']
    parts.extend([2*np.maximum(tn[ix],0)/ts,2*np.maximum(lam[ix]*np.abs(tk[ix])+np.abs(q[ix])*tn[ix],0)/((sn[ix]+np.abs(F[ix])+1)*ts),1.5*np.maximum(damp-lam[ix],0)/damp])
   values.append(np.concatenate(parts))
  return np.r_[np.concatenate(values),.01*(a-prior)]
 fit=least_squares(objective,prior,bounds=(-8,8),max_nfev=100,ftol=1e-8,xtol=1e-8,gtol=1e-8,verbose=0)
 f=WideJointModes(base,fit.x);old=WideJointModes(base,prior);rows=[]
 for k in (2.5,5.5):
  for eta in (-.2,.2):
   row=dict(k=k,eta=eta)
   for name,field in [('prior',old),('cone_fit',f)]:
    R,div,mw=sample(field,k,eta,18,.0005);row[name]=dict(momentum_max=float(np.max(np.linalg.norm(R,axis=1))),angular_terminal=float(-mw@R[:,1]),divergence_max=float(np.max(np.abs(div))))
   rows.append(row);print(json.dumps(row),flush=True)
 volume=[]
 for k in (2.5,5.5):
  tau=.5*2**(-k);gy,wy=leggauss(12);ge,we=leggauss(8);eta=np.repeat(.4*ge,12);y=np.tile((gy+1)/2,8);q=tau/(1-eta**2)
  ri=np.sqrt(2*base.nu*q*3/64);r=ri*(1+5*y);z=np.sqrt(base.nu)*q**base.inner.D*eta
  ze=np.sqrt(base.nu)*q**base.inner.D*(1+2*base.inner.D*eta**2/(1-eta**2));wv=np.repeat(.4*we,12)*np.tile(wy/2,8)*2*np.pi*r*5*ri*ze;pts=np.column_stack((r,np.zeros_like(r),z))
  for name,field in [('prior',old),('cone_fit',f)]:
   R,div=independent_fd(field,pts,tau,.0005*np.sqrt(base.nu*tau),.00025*tau);row=dict(k=k,field=name,momentum_max=float(np.max(np.linalg.norm(R,axis=1))),physical_volume_L2=float(np.sqrt(np.sum(wv*np.sum(R*R,axis=1)))) ,divergence_max=float(np.max(np.abs(div))))
   volume.append(row);print(json.dumps(row),flush=True)
 cone=[]
 for item in training:
  for name,a in [('prior',prior),('cone_fit',fit.x)]:
   u,J,L=combine(item['base'],item['modes'],a);R=momentum((u,J,L));T,d,g,lam,tn,tk,q,sn,F=stress_and_cone(u,J,R,item['r'],item['width'],Q,base.nu)
   cone.append(dict(k=item['k'],eta=item['eta'],field=name,interior_nodes=len(interior),direction_pass=int(np.count_nonzero(d[interior])),growth_pass=int(np.count_nonzero(g[interior])),joint_pass=int(np.count_nonzero((d&g)[interior])),all_nodes_joint_pass=int(np.count_nonzero(d&g)),residual_max=float(np.max(np.linalg.norm(R,axis=1)))))
 report=dict(amplitudes=fit.x.tolist(),bounds=[-8,8],fit_nfev=fit.nfev,fit_success=bool(fit.success),fit_cost=float(fit.cost),prior_objective_cost=float(.5*np.sum(objective(prior)**2)),holdouts=rows,volume=volume,cone_training=cone,interior_y=[float((nodes[interior[0]]+1)/2),float((nodes[interior[-1]]+1)/2)],
 objective='Full momentum and radial theta/axial moments, plus soft penalty for interior approximate direction and growth rate exceeding nu/r² for azimuthal mode1; no hard cone guarantee. Approximate physical frozen tangential formula, not direct normalized paper cone.',
 limitations='Exploratory bounded fit. Fitting cone diagnostics can worsen full momentum and physical L2. Training k1,4,5.5 eta=-.3,0,.3; holdout k2.5,5.5 eta=+/-.2. Finite grids only; no supported wave, global closure or proof.',pde_validated=False,global_field_ready=False)
 out=ROOT/'cone_fit';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print('cone',cone,flush=True)
if __name__=='__main__':run()
