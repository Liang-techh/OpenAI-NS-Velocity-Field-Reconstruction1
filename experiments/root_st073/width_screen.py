"""Matched physical-volume and directional screening of variable-width bridges."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss,legfit,legint,legval
from width_field import WidthField,ROOT
from joined_field import JoinedField
from affine_momentum import jets,momentum

def run():
 old=JoinedField();new=WidthField(2);pts=old.inner.from_similarity(np.array([3/64*.9,3/64*1.4**2,3/64*2.2**2]),np.array([-.2,0,.2]),.1)
 # from_similarity only maps coordinates; no evaluation outside registered inner.
 u,p=old.fields(pts,.1);v,q=new.fields(pts,.1);equivalence=dict(velocity=float(np.max(np.abs(u-v))),pressure=float(np.max(np.abs(p-q))))
 gr,wr=leggauss(16);ge,we=leggauss(6);eta=np.repeat(.4*ge,16);y=np.tile((gr+1)/2,6);rows=[]
 for ratio in (2.,3.,4.,6.):
  f=WidthField(ratio)
  for k in (2.5,5.5):
   tau=.5*2**(-k);qq=tau/(1-eta**2);ri=np.sqrt(2*f.nu*qq*3/64);width=(ratio-1)*ri;r=ri+width*y;z=np.sqrt(f.nu)*qq**f.inner.D*eta
   ze=np.sqrt(f.nu)*qq**f.inner.D*(1+2*f.inner.D*eta**2/(1-eta**2));weights=np.repeat(.4*we,16)*np.tile(wr/2,6)*2*np.pi*r*width*ze
   pts=np.column_stack((r,np.zeros_like(r),z));u,J,linear=jets(f,pts,tau,.0005*np.sqrt(f.nu*tau),.00025*tau);res=momentum((u,J,linear));stress=np.zeros((len(r),2))
   for axial in range(6):
    sl=slice(axial*16,(axial+1)*16)
    for col,comp,power in ((0,1,2),(1,2,1)):
     anti=legint(legfit(gr,r[sl]**power*res[sl,comp],15));stress[sl,col]=-width[sl]/2*(legval(gr,anti)-legval(-1,anti))/r[sl]**power
   F=u[:,1]/r;g=np.column_stack((J[:,1,0]-F,J[:,2,0]));N=g/np.linalg.norm(g,axis=1)[:,None];lam2=-2*F*N[:,0]*(2*F*N[:,0]+np.linalg.norm(g,axis=1));tn=np.sum(stress*N,axis=1);tk=stress[:,1]*N[:,0]-stress[:,0]*N[:,1]
   c=np.sqrt(np.maximum(lam2,0))/(2*F*N[:,0]);passed=(lam2>0)&(tn<0)&(np.abs(c*tk)<-tn)
   item=dict(ratio=ratio,k=k,volume=float(weights.sum()),momentum_sampled_max=float(np.max(np.linalg.norm(res,axis=1))),momentum_volume_L2=float(np.sqrt(np.sum(weights*np.sum(res**2,axis=1)))),divergence_sampled_max=float(np.max(np.abs(np.trace(J,axis1=1,axis2=2)))) ,directional_pass_count=int(passed.sum()),sample_count=len(r),directional_pass_y=y[passed].tolist(),directional_pass_eta=eta[passed].tolist(),velocity_max=float(np.max(np.linalg.norm(u,axis=1))))
   rows.append(item);print(json.dumps(item),flush=True)
 report=dict(ratio2_equivalence=equivalence,rows=rows,quadrature=dict(radial=16,axial=6),domain='eta in[-.4,.4], ri<r<ratio*ri, all azimuths. Volume grows with ratio and L2 uses actual physical weights, not volume normalization.',
 construction='Same corrected inner jets and fixed heat amplitude; quintic swirl/pressure and septic streamfunction rebuilt for physical width=(ratio-1)*ri. Analytic psi_z includes width and endpoint motion. No fitted force.',
 limitations='Exploratory geometry comparison, no independent holdout optimization or continuum certificate. Directional screen uses physical frozen tangential approximation; it is not the paper leading-cone theorem. No axial closure, global energy or supported wave realization.',pde_validated=False,global_field_ready=False)
 out=ROOT/'width_screen';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
if __name__=='__main__':run()
