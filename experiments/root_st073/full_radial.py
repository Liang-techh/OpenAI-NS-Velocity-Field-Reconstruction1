"""Local all-momentum radial Taylor construction (ST073).

This is an autonomous local solver, not the source's all-order global scheme.
Keep all axial diffusion, time derivatives and radial pressure in the recurrence.
Axis data are the smooth ST068 traces. No global matching, force, or energy
normalization is silently added. All dimensional residuals are physical.
"""
from __future__ import annotations
from dataclasses import dataclass,asdict
from functools import lru_cache
from pathlib import Path
import json,sys
import numpy as np
from scipy.signal import convolve2d
sys.path.insert(0,str(Path(__file__).resolve().parent/'upstream'))
from source_coordinates import coordinates

LD=np.longdouble

def multiply(a,b,shape):
    a=a[:shape[0],:shape[1]]; b=b[:shape[0],:shape[1]]
    return convolve2d(a,b,mode='full')[:shape[0],:shape[1]]

def derivative(a,axis,scale):
    o=np.zeros_like(a)
    if axis==0:o[:-1,:]=np.arange(1,len(a),dtype=LD)[:,None]*a[1:,:]/scale
    else:o[:,:-1]=np.arange(1,a.shape[1],dtype=LD)[None,:]*a[:,1:]/scale
    return o

def binom_jet(b,n):
    out=np.ones(n,dtype=LD)
    for i in range(1,n):out[i]=out[i-1]*(b-i+1)/i
    return out

def q_power_jet(beta,eta,shape,h=.005,terms=120):
    """Taylor coefficients of Q**beta about (eta,1-eta**2).

    Q-(eta+zeta)^2 Q**(2h) = 1-eta**2+theta; Q(0,0)=1.
    Lagrange inversion: Q**beta=sum c_m (eta+zeta)**(2m)
      *(1-eta**2+theta)**(beta+(2h-1)*m).
    Numerical series in m is separately refined; not a rigorous enclosure.
    """
    b=LD(beta);e=LD(eta);eps=2*LD(h);d=1-e*e
    if abs(e)>.5000001:raise ValueError('Axis-jet implementation limited to |eta|<=.5')
    nz,nt=shape;o=np.zeros(shape,dtype=LD)
    for m in range(terms):
        if m==0:cm=LD(1)
        else:
            cm=b/LD(m)
            for j in range(1,m):cm*= (b+eps*m-j)/j
        gamma=b+(eps-1)*m
        z=np.zeros(nz,dtype=LD)
        z[:min(2*m+1,nz)]=binom_jet(LD(2*m),min(2*m+1,nz))*np.array([e**(2*m-i) for i in range(min(2*m+1,nz))],dtype=LD)
        tt=binom_jet(gamma,nt)*d**(gamma-np.arange(nt,dtype=LD))
        o+=cm*z[:,None]*tt[None,:]
    return o

def recur(B0,C0,P0,N,S,wz,wt):
    """Return coefficient Taylor jets A/B/C/P for the full source-nu=1 system."""
    shape=B0.shape;dz=lambda a:derivative(a,0,wz);dt=lambda a:-derivative(a,1,wt)
    mul=lambda a,b:multiply(a,b,shape)
    B=[B0];C=[C0];P=[P0];A=[-dz(C[0])/2]
    for n in range(N):
        rb=dt(B[n])-dz(dz(B[n]));rc=dt(C[n])+dz(P[n])-dz(dz(C[n]))
        ra=dt(A[n])-dz(dz(A[n]))
        for i in range(n+1):
            j=n-i
            rb+=2*(j+1)*mul(A[i],B[j])+mul(C[i],dz(B[j]))
            rc+=2*j*mul(A[i],C[j])+mul(C[i],dz(C[j]))
            ra+=(1+2*j)*mul(A[i],A[j])+mul(C[i],dz(A[j]))-mul(B[i],B[j])
        B.append(S*rb/(4*(n+1)*(n+2)))
        C.append(S*rc/(4*(n+1)**2))
        A.append(-dz(C[-1])/(2*(n+2)))
        P.append(2*(n+2)*A[-1]-S*ra/(2*(n+1)))
    return A,B,C,P

@dataclass(frozen=True)
class Parameters:
    order:int=8
    jet_extra:int=3
    axis_terms:int=120
    h:float=.005
    nu:float=.01
    axis_swirl:float=1.
    axial_slope:float=4.
    axial_bias:float=.02
    pressure:float=1.
    X_max:float=1/64
    eta_max:float=.5
    k_max:float=6.
    def check(self):
        if not(1<=self.order<=18 and self.jet_extra>=3 and self.axis_terms>=60):raise ValueError('Unsupported numerical settings')
        if not(0<self.h<.01 and self.nu>0 and self.axis_swirl>0 and 0<self.eta_max<=.5 and 0<self.X_max<=.25):raise ValueError('Invalid physical settings')

class FullRadialField:
    def __init__(self,p:Parameters=Parameters()):
        p.check();self.p=p;self.nu=p.nu;self.h=p.h;self.A=.5+p.h;self.D=.5-p.h
    @lru_cache(maxsize=3000)
    def axis_jets(self,e):
        p=self.p;shape=(2*p.order+p.jet_extra+2,p.order+p.jet_extra+1)
        powers={b:q_power_jet(b,e,shape,p.h,p.axis_terms) for b in [-1-p.h,-1.,-self.A,-2*self.A,-2.]}
        z=np.zeros(shape,dtype=LD);z[0,0]=LD(e);z[1,0]=1
        B=LD(p.axis_swirl)*powers[-1-p.h]
        C=LD(p.axial_slope)*multiply(z,powers[-1.],shape)+LD(p.axial_bias)*powers[-self.A]
        P=-LD(p.pressure)*powers[-2*self.A]+LD(.5*p.pressure)*multiply(multiply(z,z,shape),powers[-2.],shape)
        return B,C,P
    @lru_cache(maxsize=1500)
    def coefficients(self,e,q):
        p=self.p;q=LD(q);N=p.order;S=2*q;wz=q**LD(self.D);wt=q
        B0,C0,P0=self.axis_jets(float(e));B=[q**LD(-1-self.h)*B0];C=[q**LD(-self.A)*C0];P=[q**LD(-2*self.A)*P0]
        A,B,C,P=recur(B[0],C[0],P[0],N,S,wz,wt)
        # Store just derivative jets needed to evaluate complete physical R.
        jets=[]
        for seq in [A,B,C,P]:
            ar=np.stack(seq)
            jets.append(np.stack([ar[:,0,0],ar[:,1,0]/wz,2*ar[:,2,0]/wz**2,-ar[:,0,1]/wt],axis=1))
        return np.array(jets,dtype=LD)
    def evaluate_similarity(self,X,eta,tau,angle=0):
        X,e,tau,angle=np.broadcast_arrays(np.asarray(X,float),np.asarray(eta,float),np.asarray(tau,float),np.asarray(angle,float))
        if not all(np.isfinite(v).all() for v in [X,e,tau,angle]):raise ValueError('Finite points required')
        if np.any(X<0) or np.any(X>self.p.X_max+1e-13) or np.any(abs(e)>self.p.eta_max+1e-13):raise ValueError('Outside registered local field')
        if np.any(tau<.5*2**(-self.p.k_max)-1e-14) or np.any(tau>.5+1e-14):raise ValueError('Outside time window')
        flat=[]
        for x0,e0,t0,ang in zip(X.ravel(),e.ravel(),tau.ravel(),angle.ravel()):
            q=LD(t0)/(1-LD(e0)**2);x=LD(x0);S=2*q;r=np.sqrt(S*x);co=self.coefficients(float(e0),float(q))
            v=np.empty((4,6),dtype=LD)
            for n in range(4):
                coeff=co[n]
                v[n,:4]=[np.polynomial.polynomial.polyval(x,coeff[:,k]) for k in range(4)]
                v[n,4]=np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(coeff[:,0]))/S
                v[n,5]=np.polynomial.polynomial.polyval(x,np.polynomial.polynomial.polyder(coeff[:,0],2))/S**2
            a,b,c,pv=v[:,0];az,bz,cz,pz=v[:,1];azz,bzz,czz,pzz=v[:,2];at,bt,ct,pt=v[:,3];ass,bs,cs,ps=v[:,4];ass2,bss,css,pss=v[:,5]
            ut=np.array([r*at,r*bt,ct]);adv=np.array([r*(a*a+2*S*x*a*ass+c*az-b*b),r*(2*a*b+2*S*x*a*bs+c*bz),2*S*x*a*cs+c*cz])
            grad=np.array([2*r*ps,0,pz]);lap=np.array([r*(8*ass+4*S*x*ass2+azz),r*(8*bs+4*S*x*bss+bzz),4*cs+4*S*x*css+czz])
            u=np.array([r*a,r*b,c]);res=ut+adv+grad-lap;div=2*a+2*S*x*ass+cz
            omega=np.array([-r*bz,r*(az-2*cs),2*b+2*S*x*bs])
            # cylindrical-to-Cartesian orthogonal change; nu scaling is explicit
            ca,sa=np.cos(ang),np.sin(ang);Q=np.array([[ca,-sa,0],[sa,ca,0],[0,0,1]])
            rn=np.sqrt(self.nu)
            flat.append(dict(velocity=rn*(Q@u),pressure=self.nu*pv,residual=rn*(Q@res),residual_cylindrical=rn*res,divergence=div,vorticity=Q@omega,time_derivative=rn*(Q@ut),advection=rn*(Q@adv),pressure_gradient=rn*(Q@grad),viscous_term=rn*(Q@lap)))
        out={k:np.array([p[k] for p in flat],float).reshape(X.shape+np.asarray(flat[0][k]).shape) for k in flat[0]} if flat else {}
        return out
    def from_similarity(self,X,e,tau,angle=0):
        X,e,tau,angle=np.broadcast_arrays(X,e,tau,angle);q=tau/(1-e*e);r=np.sqrt(2*self.nu*q*X)
        return np.stack([r*np.cos(angle),r*np.sin(angle),np.sqrt(self.nu)*q**self.D*e],axis=-1)
    def evaluate(self,points,tau):
        pts=np.asarray(points,float)
        if pts.ndim!=2 or pts.shape[1]!=3 or not np.isfinite(pts).all():raise ValueError('Finite (n,3) input required')
        t=np.broadcast_to(np.asarray(tau,float),(len(pts),));src=pts/np.sqrt(self.nu)
        r=np.hypot(src[:,0],src[:,1]);coord=coordinates(r,src[:,2],t,self.h)
        return self.evaluate_similarity(coord['X'],coord['eta'],t,np.arctan2(src[:,1],src[:,0]))
    def fields(self,points,tau):
        o=self.evaluate(points,tau);return o['velocity'],o['pressure']
    def save(self,path,ident='ST073-F'):
        p=Path(path)
        if p.exists():raise FileExistsError(p)
        p.write_text(json.dumps(dict(id=ident,parameters=asdict(self.p),scope='Finite local radial series for full unforced momentum; no global matching or energy normalization',pde_validated=False,global_field_ready=False),indent=2)+'\n')
    @classmethod
    def load(cls,path):return cls(Parameters(**json.loads(Path(path).read_text())['parameters']))
