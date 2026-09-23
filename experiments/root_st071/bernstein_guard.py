"""Sufficient positivity tests for finite axial polynomials, with subdivisions.
All coefficients are computed in floating point; this is not interval certification.
Unlike a sample grid, positive Bernstein coefficients imply polynomial positivity
mathematically. Numerical roundoff must be separately bounded for a strict proof.
"""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial import chebyshev as ch

def transform(degree,a,b):
    """Chebyshev x coefficients -> degree-n Bernstein coefficients on [a,b]."""
    n=degree;polys=np.zeros((n+1,n+1));polys[0,0]=1
    if n:
        polys[1,:2]=[a,b-a]
    for j in range(1,n):
        polys[j+1]+=2*a*polys[j]-polys[j-1]
        polys[j+1,1:]+=2*(b-a)*polys[j,:-1]
    elevation=np.zeros((n+1,n+1))
    for k in range(n+1):
        for i in range(k+1):elevation[k,i]=math.comb(k,i)/math.comb(n,i)
    return elevation@polys.T

def product_rule(n):
    i,j=np.meshgrid(np.arange(n+1),np.arange(n+1),indexing='ij');i=i.ravel();j=j.ravel();W=np.zeros((len(i),2*n+1))
    for row,(ii,jj) in enumerate(zip(i,j)):W[row,ii+jj]=math.comb(n,int(ii))*math.comb(n,int(jj))/math.comb(2*n,int(ii+jj))
    return i,j,W

def certificate(F,U,subdivisions=32,target=2.12,amin=.05):
    n=F.shape[1]-1;nodes=np.linspace(-1,1,subdivisions+1)
    trans=np.array([transform(n,a,b) for a,b in zip(nodes[:-1],nodes[1:])])
    fc=ch.chebval(1,F);fx=8*ch.chebval(1,ch.chebder(F,axis=0));ux=8*ch.chebval(1,ch.chebder(U,axis=0));f=trans@fc;x=trans@fx;u=trans@ux
    i,j,W=product_rule(n)
    G=(.25*x[:,i]*x[:,j]+.5*u[:,i]*u[:,j]+.5*target*x[:,i]*f[:,j])@W
    A=-.5*x-amin*f
    return dict(G=G,A=A,F=f,nodes=nodes/2,min_G=float(G.min()),min_A=float(A.min()),min_F=float(f.min()),scope='Floating sufficient Bernstein condition; no interval rounding bound')
