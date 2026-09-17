"""Finite annular complete-curl modes. Autonomous basis, not paper wave data.
Pure Cartesian polynomial potentials, smooth annular bump, no axis division.
Mode values and first/second velocity derivatives use analytic product jets.
"""
from __future__ import annotations
from functools import lru_cache
from itertools import product
from math import factorial, comb
import numpy as np
from spacetime import bump_derivatives

ZERO=(0,0,0)
UNITS=((1,0,0),(0,1,0),(0,0,1))

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def mul_poly(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():out[add(i,j)]=out.get(add(i,j),0.)+u*v
    return {k:v for k,v in out.items() if v!=0}
def scale_poly(a,c):return {k:c*v for k,v in a.items()}
def sum_poly(a,b):
    out=a.copy()
    for k,v in b.items():out[k]=out.get(k,0.)+v
    return {k:v for k,v in out.items() if v!=0}

def solid(m,phase):
    out={}
    for k in range(m+1):
        a=comb(m,k)*(1j)**k;c=a.real if phase=='cos' else a.imag
        if c:out[(m-k,k,0)]=c
    return out

def pderiv(poly,alpha,x):
    out=np.zeros(len(x))
    for powers,c in poly.items():
        if any(p<a for p,a in zip(powers,alpha)):continue
        term=np.full(len(x),float(c))
        for j,(p,a) in enumerate(zip(powers,alpha)):
            term*=factorial(p)/factorial(p-a)*x[:,j]**(p-a)
        out+=term
    return out

def envelope_jets(points,rin=.20,rout=1.80,zout=1.80):
    x,y,z=np.asarray(points).T;s=x*x+y*y;cen=(rout*rout+rin*rin)/2;h=(rout*rout-rin*rin)/2
    q=((s-cen)/h)**2;v=bump_derivatives(q);qs=2*(s-cen)/h**2;qss=2/h**2
    R=[v[0],v[1]*qs,v[2]*qs**2+v[1]*qss,v[3]*qs**3+3*v[2]*qs*qss]
    v=bump_derivatives(z*z/zout**2);qz=2*z/zout**2;qzz=2/zout**2
    Z=[v[0],v[1]*qz,v[2]*qz**2+v[1]*qzz,v[3]*qz**3+3*v[2]*qz*qzz]
    out={}
    for a in range(4):
        for b in range(4-a):
            for c in range(4-a-b):
                row=np.zeros(len(x))
                for ka in range(a//2+1):
                    for kb in range(b//2+1):
                        co=factorial(a)/factorial(ka)/factorial(a-2*ka)*factorial(b)/factorial(kb)/factorial(b-2*kb)
                        row+=co*(2*x)**(a-2*ka)*(2*y)**(b-2*kb)*R[a+b-ka-kb]
                out[(a,b,c)]=row*Z[c]
    return out

class CurlModes:
    def __init__(self,ms=(1,2),radial_degrees=(0,1),axial_degrees=(0,1)):
        self.potentials=[];self.meta=[];self.pressures=[];self.pressure_meta=[]
        s={ (2,0,0):1.,(0,2,0):1.}
        centered=sum_poly(scale_poly(s,1/2),{ZERO:-.8})
        for m in ms:
            for phase in ('cos','sin'):
                for ir in radial_degrees:
                    for iz in axial_degrees:
                        P=solid(m,phase)
                        if ir==1:P=mul_poly(P,centered)
                        elif ir!=0:raise ValueError('Unsupported radial mode')
                        P=mul_poly(P,{(0,0,iz):1.})
                        self.pressures.append(P);self.pressure_meta.append((m,phase,ir,iz))
                        self.potentials.append(({}, {}, P));self.meta.append(('Az',m,phase,ir,iz))
                        self.potentials.append((scale_poly(mul_poly(P,{(0,1,0):1}),-1),mul_poly(P,{(1,0,0):1}),{}));self.meta.append(('Atheta',m,phase,ir,iz))
        self.size=len(self.meta)
    def jets(self,points):
        x=np.asarray(points,float);env=envelope_jets(x);n=len(x)
        def scalar(poly,alpha):
            result=np.zeros(n)
            if not poly:return result
            for gamma in product(*(range(a+1) for a in alpha)):
                co=np.prod([comb(a,g) for a,g in zip(alpha,gamma)])
                result+=co*pderiv(poly,gamma,x)*env[sub(alpha,gamma)]
            return result
        vals=np.empty((n,self.size,3));jac=np.empty((n,self.size,3,3));lap=np.empty_like(vals)
        # curl(A)_i = partial_j A_k-partial_k A_j for cyclic (i,j,k)
        for k,pot in enumerate(self.potentials):
            cache={}
            def vel(i,alpha):
                j=(i+1)%3;l=(i+2)%3
                return scalar(pot[l],add(alpha,UNITS[j]))-scalar(pot[j],add(alpha,UNITS[l]))
            for i in range(3):
                vals[:,k,i]=vel(i,ZERO)
                for j in range(3):jac[:,k,i,j]=vel(i,UNITS[j])
                lap[:,k,i]=sum(vel(i,add(e,e)) for e in UNITS)
        pg=np.empty((n,len(self.pressures),3))
        for k,p in enumerate(self.pressures):
            for i in range(3):pg[:,k,i]=scalar(p,UNITS[i])
        return vals,jac,lap,pg
