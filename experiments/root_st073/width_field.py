"""Variable-width C2 joining with exact moving-interface streamfunction derivatives."""
import math,sys,json,numpy as np
from radial_continuation import ROOT,FullRadialField
from poloidal_bridge import coefficients,septic
from swirl_bridge import traces,quintic
sys.path.insert(0,str(ROOT/'NS_ST073_Full_Local_Recurrence/upstream'))
from source_coordinates import coordinates
from heat_exterior import physical,heat_factor
from local_field import independent_fd

class WidthField:
 def __init__(self,ratio=2.):
  if ratio<=1:raise ValueError("Outer/inner radius ratio must exceed1")
  self.ratio=float(ratio)
  self.experimental_time_extension=False
  self.inner=FullRadialField.load(ROOT/'radial_continuation/candidate.json');self.nu=self.inner.nu
  self.c=json.loads((ROOT/'heat_join/screen.json').read_text())['heat_amplitude']
 def fields(self,points,tau):
  pts=np.asarray(points,float);ts=np.broadcast_to(tau,(len(pts),));uv=[];pv=[];sn=np.sqrt(self.nu)
  for point,t in zip(pts,ts):
   r=np.hypot(*point[:2]);z=point[2];coord=coordinates(r/sn,z/sn,t,self.inner.h);e=float(coord['eta']);q=float(coord['q'])
   if abs(e)>.5 or t<=0 or (not self.experimental_time_extension and not .5/64<=t<=.5):
    raise ValueError('Only registered axial/time slab supported')
   ri,_,L,_=coefficients(self.inner,e,t);ro=self.ratio*ri;width=ro-ri;a=septic(L,width)
   if r<=ri:
    d=self.inner.evaluate(point[None,:],t);uv.append(d['velocity'][0]);pv.append(d['pressure'][0]);continue
   if r>=ro:
    d=physical(point[None,:],t,c=self.c);uv.append(d['velocity'][0]);pv.append(d['pressure'][0]);continue
   X=3/64;rho=ri/sn;co=self.inner.coefficients(e,q)
   def val(comp,j=0,col=0):return float(np.polynomial.polynomial.polyval(X,np.polynomial.polynomial.polyder(co[comp,:,col],j))/(2*q)**j)
   A=val(0);Cs,Css,Csss=val(2,1),val(2,2),val(2,3)
   Cz,Csz,Cssz=val(2,0,1),val(2,1,1),val(2,2,1)
   fourth=(6*Cs+24*rho*rho*Css+8*rho**4*Csss)/sn
   fixed_z=np.array([-ri*sn*rho*A,ri*Cz,Cz+2*rho*rho*Csz,(6*rho*Csz+4*rho**3*Cssz)/sn])
   riz=ri*float(coord['q_z'])/(2*q*sn)
   total_z=fixed_z+np.r_[L[1:],fourth]*riz
   # Differentiate the linear Hermite solve at fixed dimensionless radius.
   modified=np.array([total_z[j]+j*L[j]*riz/ri for j in range(4)])
   az=septic(modified,width)
   y=(r-ri)/width;yp=-(1+(self.ratio-1)*y)*riz/width
   psi_r=np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(a))/width
   psi_z=np.polynomial.polynomial.polyval(y,az)+np.polynomial.polynomial.polyval(y,np.polynomial.polynomial.polyder(a))*yp
   _,_,ls,_=traces(self.inner,e,t,self.c)
   rr=ro/sn;ss=rr**2/2;exponent=.5+self.inner.h;Z=2*t/ss
   H,H1,H2=[float(heat_factor(Z,self.inner.h,j)) for j in range(3)]
   K=self.c*ss**(-exponent)*H;Ks=self.c*ss**(-exponent-1)*(-exponent*H-Z*H1);Kss=self.c*ss**(-exponent-2)*(exponent*(exponent+1)*H+2*(exponent+1)*Z*H1+Z**2*H2)
   rs=np.array([sn*K,rr*Ks,(Ks+rr**2*Kss)/sn]);sw=quintic(ls,rs,width);uth=np.polynomial.polynomial.polyval(y,sw)
   outpoint=np.array([[ro,0,z]]);outer=physical(outpoint,t,c=self.c)
   pressure_left=[self.nu*val(3),2*sn*rho*val(3,1),2*val(3,1)+4*rho*rho*val(3,2)]
   pressure_right=[float(outer['pressure'][0]),rs[0]**2/ro,2*rs[0]*rs[1]/ro-rs[0]**2/ro**2]
   pressure=np.polynomial.polynomial.polyval(y,quintic(pressure_left,pressure_right,width))
   ur=-psi_z/r;uz=psi_r/r;ca,sa=point[0]/r,point[1]/r
   uv.append([ur*ca-uth*sa,ur*sa+uth*ca,uz]);pv.append(pressure)
  return np.array(uv),np.array(pv)
