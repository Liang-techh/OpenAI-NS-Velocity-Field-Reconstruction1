"""Read-only directional quadrature, precision, and boundary diagnostics."""
import sys,time,json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import full_radial as fr
from audit import save,ROOT
from audit_phase2 import baseline,binding
from local_field import LocalField

def run_quad():
 binding();f=fr.FullRadialField.load(ROOT/'data/ST073-V.json');b=LocalField(baseline());rows=[]
 # Separate radial refinement at fixed axial n=32 and axial at fixed radial n=18.
 for k,nr,ne in [(0,18,32),(3,18,32),(6,18,24),(6,18,32),(6,12,32),(6,24,32)]:
  t=time.monotonic();tau=.5*2**(-k);r,w=leggauss(nr);e,we=leggauss(ne);xx,ee=np.meshgrid((r+1)/128,e/2,indexing='ij');q=tau/(1-ee*ee)
  weights=w[:,None]*we[None,:]/256*2*np.pi*.01**1.5*q**(1+f.D)*(1-2*f.h*ee*ee)/(1-ee*ee)
  d=f.evaluate_similarity(xx,ee,tau);old=b.evaluate(f.from_similarity(xx.ravel(),ee.ravel(),tau),tau);v={}
  for name,out in [('full',d),('leading_same_axis',old)]:
   R=out['residual'].reshape(*xx.shape,3);u=out['velocity'].reshape(*xx.shape,3);en=np.sum(u*u,axis=-1)
   v[name]=dict(L2=float(np.sqrt(np.sum(weights*np.sum(R*R,axis=-1)))),max=float(np.linalg.norm(R,axis=-1).max()),energy=float(.5*np.sum(weights*en)))
  row=dict(k=k,nr=nr,ne=ne,results=v,seconds=time.monotonic()-t);rows.append(row);save(ROOT/'evidence'/f'extra_quad_k{k}_r{nr}_e{ne}.json',row);print(row,flush=True)
 save(ROOT/'evidence/directional_quadrature.json',dict(rows=rows,scope='Same frozen field; radial and axial quadrature varied separately. No fit or continuous bound.'))

def precision():
 f=fr.FullRadialField(fr.Parameters(order=8,axis_swirl=4));X=np.array([.002,.01,1/64]);es=np.array([-.5,0,.43]);a=f.evaluate_similarity(X,es,.0078125)
 oldLD=fr.LD
 try:
  fr.LD=np.float64;g=fr.FullRadialField(fr.Parameters(order=8,axis_swirl=4));b=g.evaluate_similarity(X,es,.0078125)
 finally:fr.LD=oldLD
 d={k:float(np.max(abs(a[k]-b[k]))) for k in ['velocity','pressure','residual']}
 save(ROOT/'evidence/arithmetic_sensitivity.json',dict(max_absolute_difference=d,scope='Three fixed points, binary64 versus long-double implementation, not interval arithmetic'))
 print(d,flush=True)
if __name__=='__main__':run_quad();precision()
