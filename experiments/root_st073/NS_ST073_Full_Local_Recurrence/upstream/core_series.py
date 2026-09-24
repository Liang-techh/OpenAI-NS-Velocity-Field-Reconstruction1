"""Regular-axis radial recurrence for the leading equations (4.7),(4.13).

Independent implementation in bivariate power coefficients (X,eta).
This is a LOCAL, leading-order profile, not a global compact NS candidate.
Axis pressure is an explicitly autonomous datum awaiting outer matching.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import numpy as np
from numpy.polynomial import polynomial as pp
from scipy.signal import convolve2d


def deriv(c: np.ndarray, x: int=0, eta: int=0) -> np.ndarray:
    return pp.polyder(pp.polyder(c,m=x,axis=0),m=eta,axis=1)


def evaluate(c, x, eta):
    x,eta=np.broadcast_arrays(np.asarray(x),np.asarray(eta))
    return pp.polyval2d(x,eta,c)


def product(a,b, eta_degree=None):
    out=convolve2d(np.asarray(a),np.asarray(b),mode='full')
    if eta_degree is not None:out=out[:,:eta_degree+1]
    return out


def add(*arrays):
    nr=max(a.shape[0] for a in arrays);ne=max(a.shape[1] for a in arrays)
    out=np.zeros((nr,ne),dtype=np.result_type(*arrays))
    for a in arrays:out[:a.shape[0],:a.shape[1]]+=a
    return out


def shift(c,x=0,eta=0):
    return np.pad(c,((x,0),(eta,0)))


def integrate_x(c):
    return np.vstack((np.zeros((1,c.shape[1])),c/np.arange(1,len(c)+1)[:,None]))


def radial_average(c):return c/np.arange(1,len(c)+1)[:,None]


def invL(h,degree):
    out=np.zeros(degree+1);out[::2]=(2*h)**np.arange(len(out[::2]));return out


def pder(a):
    out=np.zeros_like(a);out[:-1]=np.arange(1,len(a))*a[1:];return out


def pmul(a,b,n):
    return np.convolve(a,b)[:n]


def eta_mult(a,k=1):
    b=np.zeros_like(a)
    if k==0:return a.copy()
    b[k:]=a[:-k];return b


def build(order=18,eta_degree=64,h=.005,swirl=1.,slope=4.,bias=.02,
          pressure=1.,pressure_shape=.5, checkpoint_dir=None):
    """Analytic Cauchy data: F0=swirl,U0=slope*eta+bias,Pi0=-P+P*c*eta².
    No wall or truncated-eta inflow boundary conditions are imposed.
    Truncation degree in eta is independent of radial truncation order.
    """
    if not(0<h<.01 and 0<=order<=48 and eta_degree>=8):raise ValueError('Invalid truncation/h')
    m=eta_degree+1;F=np.zeros((order+1,m));U=np.zeros_like(F);P=np.zeros((order+2,m))
    F[0,0]=swirl;U[0,:2]=[bias,slope];P[0,[0,2]]=[-pressure,pressure*pressure_shape]
    one=np.zeros(m);one[0]=1;d=one.copy();d[2]=-1
    Li=invL(h,eta_degree);A=.5+h;D=.5-h
    hist=[]
    for n in range(order):
        W=np.zeros((n+1,m));H=np.zeros_like(W)
        for i in range(n+1):
            W[i]=-(2*D*eta_mult(U[i])+pmul(d,pder(U[i]),m))/(i+1)
            H[i]=pmul(d,U[i],m)
        W[0]+=one;H[0,1]+=D
        Rf=np.zeros(m);Ru=np.zeros(m);sqF=np.zeros(m)
        for i in range(n+1):
            j=n-i
            Rf+=(j+1)*pmul(W[i],F[j],m)+pmul(H[i],pder(F[j]),m)
            Rf-=2*h*eta_mult(pmul(U[i],F[j],m))
            Ru+=j*pmul(W[i],U[j],m)+pmul(H[i],pder(U[j]),m)
            Ru-=2*A*eta_mult(pmul(U[i],U[j],m))
            sqF+=pmul(F[i],F[j],m)
        Rf+=h*F[n]
        Ru+=A*U[n]+pmul(d,pder(P[n]),m)-(4*A+2*n)*eta_mult(P[n])
        F[n+1]=pmul(Li,Rf,m)/(2*(n+1)*(n+2))
        U[n+1]=pmul(Li,Ru,m)/(2*(n+1)**2)
        P[n+1]=sqF/(n+1)
        hist.append({'radial_order':n+1,'F_coefficient_max':float(abs(F[n+1]).max()),'U_coefficient_max':float(abs(U[n+1]).max())})
        if checkpoint_dir is not None:
            dest=Path(checkpoint_dir);dest.mkdir(parents=True,exist_ok=True)
            tmp=dest/f'order_{n+1:02d}.tmp'
            with tmp.open('wb') as f:np.savez_compressed(f,F=F[:n+2],U=U[:n+2],P_recursion=P[:n+2])
            tmp.replace(dest/f'order_{n+1:02d}.npz')
    # Keep ALL squared-polynomial coefficients for the pressure identity.
    P=integrate_x(product(F,F));P[0,:3]=[-pressure,0,pressure*pressure_shape]
    meta=dict(schema='local_axis_series_v1',id='ST068-I',order=order,eta_degree=eta_degree,h=h,
              axis=dict(swirl=swirl,slope=slope,bias=bias,pressure=pressure,pressure_shape=pressure_shape),
              X_max=.25,eta_max=1.,nu=.01,T=1.,tau_range=[.5/64,.5],
              pde_validated=False,global_field_ready=False,matched_exterior=False,
              source_correspondence_verified=False,blowup_proved=False)
    return Core(F,U,P,meta),hist


@dataclass
class Core:
    F:np.ndarray
    U:np.ndarray
    P:np.ndarray
    meta:dict
    def __post_init__(self):
        self.h=float(self.meta['h']);self.A=.5+self.h;self.D=.5-self.h
        self.cache={}
        avg=radial_average(self.U)
        # v0 = numerator / L; V0=X*v0. Numerator is a polynomial.
        self.vnum=add(2*shift(self.U,eta=1),-2*self.D*shift(avg,eta=1),-deriv(avg,eta=1),shift(deriv(avg,eta=1),eta=2))
    def dval(self,name,X,eta,nx=0,ne=0):
        key=(name,nx,ne)
        if key not in self.cache:self.cache[key]=deriv(getattr(self,name),nx,ne)
        return evaluate(self.cache[key],X,eta)
    def jets(self,name,X,eta):
        out={(i,j):self.dval(name,X,eta,i,j) for i,j in [(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]}
        if name=='vnum':
            L=1-2*self.h*eta**2;inv=1/L;inv1=4*self.h*eta/L**2;inv2=4*self.h/L**2+32*self.h**2*eta**2/L**3
            a=out.copy()
            out={(0,0):a[0,0]*inv,(1,0):a[1,0]*inv,(2,0):a[2,0]*inv,
                 (0,1):a[0,1]*inv+a[0,0]*inv1,(1,1):a[1,1]*inv+a[1,0]*inv1,
                 (0,2):a[0,2]*inv+2*a[0,1]*inv1+a[0,0]*inv2}
        return out
    def profiles(self,X,eta,check=True):
        X,eta=np.broadcast_arrays(np.asarray(X),np.asarray(eta))
        if check and (np.any(~np.isfinite(X)) or np.any(~np.isfinite(eta)) or np.any(X<0) or np.any(X>self.meta['X_max']+1e-13) or np.any(abs(eta)>1+1e-13)):
            raise ValueError('Outside declared LOCAL profile; no global extrapolation')
        F=self.jets('F',X,eta);U=self.jets('U',X,eta);P=self.jets('P',X,eta);v=self.jets('vnum',X,eta)
        avg=self.dval('U',X,eta)*0+evaluate(radial_average(self.U),X,eta)
        avge=evaluate(deriv(radial_average(self.U),eta=1),X,eta)
        L=1-2*self.h*eta**2;d=1-eta**2
        W=1-2*self.D*eta*avg-d*avge;Hc=self.D*eta+d*U[0,0]
        # Positive viscosity-minus-transport defect; exact root of eq(4.13).
        ef=2*L*(X*F[2,0]+2*F[1,0])-(W+self.h*(1-2*eta*U[0,0]))*F[0,0]-W*X*F[1,0]-Hc*F[0,1]
        eu=2*L*(X*U[2,0]+U[1,0])-W*X*U[1,0]-self.A*(1-2*eta*U[0,0])*U[0,0]-Hc*U[0,1]-d*P[0,1]+4*self.A*eta*P[0,0]+2*eta*X*P[1,0]
        return dict(F=F,U=U,P=P,v=v,leading_F_defect=ef,leading_U_defect=eu,
                    pressure_defect=P[1,0]-F[0,0]**2,divergence_defect=v[0,0]+X*v[1,0]-(2*self.A*eta*U[0,0]-d*U[0,1]+2*eta*X*U[1,0])/L)
    def save(self,path):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        if path.exists():raise FileExistsError(path)
        with path.open('wb') as f:np.savez_compressed(f,F=self.F,U=self.U,P=self.P,metadata=np.array(json.dumps(self.meta)))
    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as d:return cls(d['F'],d['U'],d['P'],json.loads(str(d['metadata'])))
