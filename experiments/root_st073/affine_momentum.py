"""Reusable Cartesian FD jets for affine-field nonlinear momentum fitting."""
import numpy as np

def jets(field,pts,tau,hspace,htime):
 pts=np.asarray(pts);u,p=field.fields(pts,tau)
 grad=np.empty((len(pts),3,3));gp=np.empty((len(pts),3));lap=np.zeros_like(u)
 for i in range(3):
  e=np.zeros(3);e[i]=hspace
  um2,pm2=field.fields(pts-2*e,tau);um,pm=field.fields(pts-e,tau)
  up,pp=field.fields(pts+e,tau);up2,pp2=field.fields(pts+2*e,tau)
  grad[:,:,i]=(um2-8*um+8*up-up2)/(12*hspace)
  gp[:,i]=(pm2-8*pm+8*pp-pp2)/(12*hspace)
  lap+=(-up2+16*up-30*u+16*um-um2)/(12*hspace*hspace)
 ut=-(field.fields(pts,tau-2*htime)[0]-8*field.fields(pts,tau-htime)[0]+8*field.fields(pts,tau+htime)[0]-field.fields(pts,tau+2*htime)[0])/(12*htime)
 return u,grad,ut+gp-field.nu*lap

def momentum(jet):
 u,g,linear=jet
 return linear+np.einsum('nij,nj->ni',g,u)

def combine(base,modes,a):
 return tuple(b+np.tensordot(a,m,axes=1) for b,m in zip(base,modes))
