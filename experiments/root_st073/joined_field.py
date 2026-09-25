"""Actual C2 velocity/pressure bridge; no dynamic or global admission."""
import math,sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
from poloidal_bridge import coefficients,septic
from swirl_bridge import traces,quintic
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from source_coordinates import coordinates
from heat_exterior import physical
from local_field import independent_fd

class JoinedField:
 def __init__(self,inner=None,join_X=3/64,heat_amplitude=None,outer_ratio=2.,swirl_bubble_amplitude=0.,outer_swirl_bubble_amplitude=0.,slope_swirl_bubble_amplitude=0.,shear_swirl_amplitude=0.,shear_swirl_cycles=16):
  self.inner=inner if inner is not None else FullRadialField.load(ROOT/'radial_continuation/candidate.json')
  self.nu=self.inner.nu;self.join_X=float(join_X);self.outer_ratio=float(outer_ratio)
  self.swirl_bubble_amplitude=float(swirl_bubble_amplitude)
  self.outer_swirl_bubble_amplitude=float(outer_swirl_bubble_amplitude)
  self.slope_swirl_bubble_amplitude=float(slope_swirl_bubble_amplitude)
  self.shear_swirl_amplitude=float(shear_swirl_amplitude)
  self.shear_swirl_cycles=int(shear_swirl_cycles)
  if self.join_X<=0 or self.join_X>self.inner.p.X_max:raise ValueError('Join must be inside the inner field')
  if self.outer_ratio<=1:raise ValueError('Outer radius must exceed inner radius')
  if self.slope_swirl_bubble_amplitude and not self.join_X<1<self.join_X*self.outer_ratio**2:raise ValueError('Slope reference X=1 must lie inside the bridge')
  if self.shear_swirl_amplitude and not self.join_X<1<self.join_X*self.outer_ratio**2:raise ValueError('Shear reference X=1 must lie inside the bridge')
  if self.shear_swirl_cycles<1 or self.shear_swirl_cycles!=shear_swirl_cycles:raise ValueError('Shear cycles must be a positive integer')
  self.c=(float(heat_amplitude) if heat_amplitude is not None else
          json.loads((ROOT/'heat_join/screen.json').read_text())['heat_amplitude'])
 def fields(self,points,tau):
  pts=np.asarray(points,float);ts=np.broadcast_to(tau,(len(pts),));uv=[];pv=[];sn=np.sqrt(self.nu)
  for point,t in zip(pts,ts):
   r=np.hypot(*point[:2]);z=point[2];coord=coordinates(r/sn,z/sn,t,self.inner.h);e=float(coord['eta']);q=float(coord['q'])
   if abs(e)>self.inner.p.eta_max or not .5/64<=t<=.5:raise ValueError('Only registered axial/time slab supported')
   ri,ro,L,a=coefficients(self.inner,e,t,self.join_X,self.outer_ratio);w=ro-ri
   if r<=ri:
    d=self.inner.evaluate(point[None,:],t);uv.append(d['velocity'][0]);pv.append(d['pressure'][0]);continue
   if r>=ro:
    d=physical(point[None,:],t,c=self.c);uv.append(d['velocity'][0]);pv.append(d['pressure'][0]);continue
   X=self.join_X;rho=ri/sn;co=self.inner.coefficients(e,q)
   def val(comp,j=0,col=0):return float(np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(co[comp,:,col],j))/(2*q)**j)
   A=val(0);Cs,Css,Csss=val(2,1),val(2,2),val(2,3)
   Cz,Csz,Cssz=val(2,0,1),val(2,1,1),val(2,2,1)
   fourth=(6*Cs+24*rho*rho*Css+8*rho**4*Csss)/sn
   fixed_z=np.array([-ri*sn*rho*A,ri*Cz,Cz+2*rho*rho*Csz,(6*rho*Csz+4*rho**3*Cssz)/sn])
   riz=ri*float(coord['q_z'])/(2*q*sn)
   total_z=fixed_z+np.r_[L[1:],fourth]*riz
   # Differentiate the linear Hermite solve at fixed dimensionless radius.
   modified=np.array([total_z[j]+j*L[j]*riz/ri for j in range(4)])
   az=septic(modified,w)
   y=(r-ri)/w;yp=-(1+(self.outer_ratio-1)*y)*riz/w
   psi_r=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(a))/w
   psi_z=np.polynomial.polynomial.polyval(y,az)+np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(a))*yp
   _,_,ls,rs=traces(self.inner,e,t,self.c,self.join_X,self.outer_ratio);sw=quintic(ls,rs,w);uth=np.polynomial.polynomial.polyval(y,sw)
   uth+=np.sqrt(self.nu)*q**(-self.inner.A)*self.swirl_bubble_amplitude*64*y**3*(1-y)**3
   uth+=np.sqrt(self.nu)*q**(-self.inner.A)*self.outer_swirl_bubble_amplitude*(64*y**3*(1-y)**3/.421875)*(y/.75)**8
   if self.slope_swirl_bubble_amplitude:
    y0=(np.sqrt(1/self.join_X)-1)/(self.outer_ratio-1)
    uth+=np.sqrt(self.nu)*q**(-self.inner.A)*self.slope_swirl_bubble_amplitude*64*y**3*(1-y)**3*(y-y0)
   if self.shear_swirl_amplitude:
    y0=(np.sqrt(1/self.join_X)-1)/(self.outer_ratio-1)
    phase=2*np.pi*self.shear_swirl_cycles*(y-y0)
    uth+=np.sqrt(self.nu)*q**(-self.inner.A)*self.shear_swirl_amplitude*64*y**3*(1-y)**3*np.sin(phase)/(2*np.pi*self.shear_swirl_cycles)
   outpoint=np.array([[ro,0,z]]);outer=physical(outpoint,t,c=self.c)
   pressure_left=[self.nu*val(3),2*sn*rho*val(3,1),2*val(3,1)+4*rho*rho*val(3,2)]
   pressure_right=[float(outer['pressure'][0]),rs[0]**2/ro,2*rs[0]*rs[1]/ro-rs[0]**2/ro**2]
   pressure=np.polynomial.polynomial.polyval(y,quintic(pressure_left,pressure_right,w))
   ur=-psi_z/r;uz=psi_r/r;ca,sa=point[0]/r,point[1]/r
   uv.append([ur*ca-uth*sa,ur*sa+uth*ca,uz]);pv.append(pressure)
  return np.array(uv),np.array(pv)

def run():
 f=JoinedField();rows=[]
 for k in (.4,5.5):
  tau=.5*2**(-k);eta=np.array([-.2,.2]);X=3/64*np.array([1.3,1.7])**2;pts=f.inner.from_similarity(X,eta,tau,angle=[.2,.5])
  for factor in (.002,.001):
   residual,div=independent_fd(f,pts,tau,factor*np.sqrt(f.nu*tau),.00025*tau)
   rows.append(dict(k=k,step_factor=factor,residual_norms=np.linalg.norm(residual,axis=1).tolist(),divergence=div.tolist()))
 report=dict(rows=rows,scope='C2 radial joining with analytic moving-interface psi_z and pressure interpolation. Only axial slab |eta|<=.5; no axial closure or finite global energy claim. Unforced full Cartesian momentum FD.',pde_validated=False,global_field_ready=False)
 out=ROOT/'joined_field';out.mkdir(exist_ok=True);(out/'report.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows,indent=2))
if __name__=='__main__':run()
