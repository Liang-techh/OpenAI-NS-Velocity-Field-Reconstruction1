"""Finite-energy solenoidal localization control, NOT an NS-accepted transition."""
from pathlib import Path
import sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from source_coordinates import coordinates
from local_field import independent_fd

def cutoff(s):
    s=np.asarray(s);a=np.clip(s,0,1)
    value=1-126*a**5+420*a**6-540*a**7+315*a**8-70*a**9
    derivative=-630*a**4*(1-a)**4
    return np.where(s<=0,1,np.where(s>=1,0,value)),np.where((s>0)&(s<1),derivative,0)

class CompactControl:
    def __init__(self):
        self.inner=FullRadialField.load(ROOT/'radial_continuation/candidate.json');self.nu=self.inner.nu
    def fields(self,points,tau):
        if np.any(np.asarray(tau)<.5/64) or np.any(np.asarray(tau)>.5):raise ValueError('Registered finite time window only')
        pts=np.asarray(points,float);src=pts/np.sqrt(self.nu);r=np.hypot(src[:,0],src[:,1]);c=coordinates(r,src[:,2],tau,self.inner.h)
        X,e=c['X'],c['eta'];active=(X<3/64)&(abs(e)<.5);velocity=np.zeros_like(pts);pressure=np.zeros(len(pts))
        if not np.any(active):return velocity,pressure
        idx=np.flatnonzero(active);xx=X[idx];ee=e[idx];rr=r[idx]
        d=self.inner.evaluate_similarity(xx,ee,np.broadcast_to(tau,(len(pts),))[idx])
        R,Rd=cutoff((xx-1/64)/(2/64));Z,Zd=cutoff((abs(ee)-.3)/.2)
        H=R*Z;Hx=Rd/(2/64)*Z;He=R*Zd/.2*np.sign(ee)
        hz=Hx*c['X_z'][idx]+He*c['eta_z'][idx];hr=Hx*c['X_r'][idx]
        psi=[]
        for x,eta,q in zip(xx,ee,c['q'][idx]):
            co=self.inner.coefficients(float(eta),float(q))[2,:,0]
            psi.append(float(q*np.polynomial.polynomial.polyval(x,np.r_[0,co/np.arange(1,len(co)+1)])))
        pr=np.divide(psi,rr,out=np.zeros(len(rr)),where=rr>0)
        # u_phys=sqrt(nu)*u_source; streamfunction product differentiated analytically.
        ur=H*d['velocity'][:,0]-np.sqrt(self.nu)*pr*hz
        uz=H*d['velocity'][:,2]+np.sqrt(self.nu)*pr*hr
        ut=H*d['velocity'][:,1];ca=np.divide(src[idx,0],rr,out=np.ones(len(rr)),where=rr>0);sa=np.divide(src[idx,1],rr,out=np.zeros(len(rr)),where=rr>0)
        velocity[idx]=np.column_stack((ur*ca-ut*sa,ur*sa+ut*ca,uz));pressure[idx]=H*d['pressure']
        return velocity,pressure

def run():
 f=CompactControl();rows=[]
 for k in (.4,5.5):
  tau=.5*2**(-k);pts=f.inner.from_similarity([.008,.025,.008,.03],[.1,.1,.4,.4],tau,angle=[.2,.3,.4,.5])
  u,p=f.fields(pts,tau);core=f.inner.evaluate(pts[:1],tau)
  for factor in (.002,.001):
   res,div=independent_fd(f,pts,tau,factor*np.sqrt(f.nu*tau),.00025*tau)
   rows.append(dict(k=k,step_factor=factor,regions=['core','radial_collar','axial_collar','corner'],momentum_norms=np.linalg.norm(res,axis=1).tolist(),divergence=div.tolist(),core_velocity_replay_error=float(np.max(abs(u[0]-core['velocity'][0])))))
 out=ROOT/'compact_control';out.mkdir(exist_ok=True)
 report=dict(rows=rows,scope='Streamfunction C4 taper: preserves only X<=1/64 and |eta|<=.3; NOT the entire frozen ST073 domain. Compact support X<3/64, |eta|<.5, finite energy at each registered time by bounded compact velocity. No force added; pressure tapered explicitly; transition momentum NOT assumed small.',global_field_ready=False,pde_validated=False)
 (out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows,indent=2))
if __name__=='__main__':run()
