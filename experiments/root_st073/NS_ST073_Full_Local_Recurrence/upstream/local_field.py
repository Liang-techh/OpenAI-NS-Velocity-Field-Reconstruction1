"""Physical evaluator of a SOURCE LEADING inner field (not a global NS solution).
All time chain rules, radial inertia, and axial viscosity are retained in R.
Input coordinates are physical at nu=.01 (or the explicitly requested nu).
The local audit uses no force; it does not change the old fixed-force benchmark.
"""
from __future__ import annotations
import numpy as np
from core_series import Core
from source_coordinates import coordinates


def operator(j,X,e,b,h):
    """T_b, Z_b, Z_{b-D} Z_b from a profile's second-order jets."""
    D=.5-h;d=1-e*e;L=1-2*h*e*e
    f,fx,fxx,fe,fxe,fee=[j[k] for k in [(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]]
    tb=(-b*f+D*e*fe+X*fx)/L
    num=2*b*e*f+d*fe-2*e*X*fx
    zb=num/L
    gx=(2*b*e*fx+d*fxe-2*e*(fx+X*fxx))/L
    ge=(2*b*f+2*b*e*fe-2*e*fe+d*fee-2*X*fx-2*e*X*fxe)/L+num*4*h*e/L**2
    zz=(2*(b-D)*e*zb+d*ge-2*e*X*gx)/L
    return tb,zb,zz


def scaled(j,X,e,q,b,h,factor=1.):
    tb,zb,zz=operator(j,X,e,b,h);D=.5-h
    return dict(val=factor*q**b*j[0,0],s=factor*.5*q**(b-1)*j[1,0],
                ss=factor*.25*q**(b-2)*j[2,0],t=factor*q**(b-1)*tb,
                z=factor*q**(b-D)*zb,zz=factor*q**(b-2*D)*zz)


class LocalField:
    def __init__(self,core:Core,nu=.01):
        if not np.isfinite(nu) or nu<=0:raise ValueError('Positive viscosity')
        self.core=core;self.nu=nu;self.rootnu=np.sqrt(nu)
    def evaluate(self,points,tau,check=True):
        pts=np.asarray(points,float)
        if pts.ndim!=2 or pts.shape[1]!=3 or not np.isfinite(pts).all():raise ValueError('Finite (n,3) points required')
        tau=np.broadcast_to(np.asarray(tau,float),(len(pts),))
        if np.any(tau<self.core.meta['tau_range'][0]-1e-14) or np.any(tau>self.core.meta['tau_range'][1]+1e-14):raise ValueError('Outside registered finite tau interval')
        src=pts/self.rootnu;x,y,z=src.T;s=x*x+y*y;r=np.sqrt(s)
        coord=coordinates(r,z,tau,self.core.h);q,e,X=[coord[k] for k in ['q','eta','X']]
        p=self.core.profiles(X,e,check)
        aa=scaled(p['v'],X,e,q,-1,self.core.h,.5)
        bb=scaled(p['F'],X,e,q,-1-self.core.h,self.core.h)
        cc=scaled(p['U'],X,e,q,-self.core.A,self.core.h)
        pp=scaled(p['P'],X,e,q,-2*self.core.A,self.core.h)
        a,b,c=aa['val'],bb['val'],cc['val']
        u=np.c_[x*a-y*b,y*a+x*b,c]
        ut=np.c_[x*aa['t']-y*bb['t'],y*aa['t']+x*bb['t'],cc['t']]
        conv_a=a*a+2*s*a*aa['s']+c*aa['z']-b*b
        conv_b=2*a*b+2*s*a*bb['s']+c*bb['z']
        conv_c=2*s*a*cc['s']+c*cc['z']
        conv=np.c_[x*conv_a-y*conv_b,y*conv_a+x*conv_b,conv_c]
        la=8*aa['s']+4*s*aa['ss']+aa['zz'];lb=8*bb['s']+4*s*bb['ss']+bb['zz']
        lc=4*cc['s']+4*s*cc['ss']+cc['zz']
        lap=np.c_[x*la-y*lb,y*la+x*lb,lc]
        gradp=np.c_[2*x*pp['s'],2*y*pp['s'],pp['z']]
        R=ut+conv+gradp-lap
        ezero=np.divide(x,r,out=np.ones_like(x),where=r>0);eone=np.divide(y,r,out=np.zeros_like(y),where=r>0)
        Rc=np.c_[R[:,0]*ezero+R[:,1]*eone,-R[:,0]*eone+R[:,1]*ezero,R[:,2]]
        angular_axial_diff=r*bb['zz'];axial_axial_diff=cc['zz']
        lead=np.c_[Rc[:,1]+angular_axial_diff,Rc[:,2]+axial_axial_diff]
        algebraic_lead=-q[:,None]**(-self.core.A-1)/ (1-2*self.core.h*e[:,None]**2)*np.c_[np.sqrt(2*X)*p['leading_F_defect'],p['leading_U_defect']]
        # Source-space divergence equals physical divergence under u_nu=sqrt(nu)u.
        div=2*a+2*s*aa['s']+cc['z']
        omega=np.c_[-y*(aa['z']-2*cc['s'])-x*bb['z'],x*(aa['z']-2*cc['s'])-y*bb['z'],2*b+2*s*bb['s']]
        return dict(velocity=self.rootnu*u,pressure=self.nu*pp['val'],residual=self.rootnu*R,
                    residual_cylindrical=self.rootnu*Rc,divergence=div,vorticity=omega,
                    time_derivative=self.rootnu*ut,advection=self.rootnu*conv,pressure_gradient=self.rootnu*gradp,
                    viscous_term=self.rootnu*lap,axial_diffusion_tangential=self.rootnu*np.c_[angular_axial_diff,axial_axial_diff],
                    leading_tangential=self.rootnu*lead,leading_from_profile=self.rootnu*algebraic_lead,
                    X=X,eta=e,q=q,scope='Local leading field; full unforced residual, no exterior/support/global-energy completion')
    def fields(self,points,tau):
        d=self.evaluate(points,tau);return d['velocity'],d['pressure']
    def from_similarity(self,X,e,tau,angle=0):
        X,e,tau,angle=np.broadcast_arrays(np.asarray(X),np.asarray(e),np.asarray(tau),np.asarray(angle))
        if np.any(abs(e)>=1):raise ValueError('Finite physical points require |eta|<1')
        q=tau/(1-e*e);r=np.sqrt(2*self.nu*q*X);z=self.rootnu*q**self.core.D*e
        return np.stack((r*np.cos(angle),r*np.sin(angle),z),axis=-1)


def independent_fd(field,pts,tau,hspace,htime):
    """Fourth-order Cartesian derivatives at FIXED physical x; tau=T-t."""
    pts=np.asarray(pts);tau=np.broadcast_to(tau,(len(pts),));u,p=field.fields(pts,tau)
    grad=np.empty((len(pts),3,3));gp=np.empty((len(pts),3));lap=np.zeros_like(u)
    for i in range(3):
        e=np.zeros(3);e[i]=hspace
        um2,pm2=field.fields(pts-2*e,tau);um,pm=field.fields(pts-e,tau)
        up,pp=field.fields(pts+e,tau);up2,pp2=field.fields(pts+2*e,tau)
        grad[:,:,i]=(um2-8*um+8*up-up2)/(12*hspace)
        gp[:,i]=(pm2-8*pm+8*pp-pp2)/(12*hspace)
        lap+=(-up2+16*up-30*u+16*um-um2)/(12*hspace*hspace)
    # In increasing physical time, tau decreases.
    utau=(field.fields(pts,tau-2*htime)[0]-8*field.fields(pts,tau-htime)[0]+8*field.fields(pts,tau+htime)[0]-field.fields(pts,tau+2*htime)[0])/(12*htime)
    residual=-utau+np.einsum('nij,nj->ni',grad,u)+gp-field.nu*lap
    return residual,np.trace(grad,axis1=1,axis2=2)
