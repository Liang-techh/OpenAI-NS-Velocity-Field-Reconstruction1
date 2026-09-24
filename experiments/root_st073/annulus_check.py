"""Finite annulus quadrature for the explicit order10 ST073 radial continuation."""
from dataclasses import replace
from pathlib import Path
import json,numpy as np
from numpy.polynomial.legendre import leggauss
from radial_continuation import ROOT,FullRadialField

def run():
    original=json.loads((ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V.json').read_text())
    base=FullRadialField.load(ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V.json')
    f=FullRadialField(replace(base.p,order=10,X_max=3/64));rows=[]
    for nr,ne in ((8,12),(12,18)):
        a,wa=leggauss(nr);b,wb=leggauss(ne)
        xs=1/64+(a+1)/64;es=.5*b
        X,e=np.meshgrid(xs,es,indexing='ij');weights=np.outer(wa/64,wb*.5)
        for k in (0,3,6):
            tau=.5*2**(-k);q=tau/(1-e*e);jac=2*np.pi*f.nu**1.5*q**(1+f.D)*(1+2*f.D*e*e/(1-e*e))
            d=f.evaluate_similarity(X,e,tau);rr=np.sum(d['residual']**2,axis=-1)
            rows.append(dict(k=k,radial_order=nr,axial_order=ne,volume=float(np.sum(weights*jac)),momentum_max=float(np.sqrt(rr.max())),momentum_volume_L2=float(np.sqrt(np.sum(weights*jac*rr))),divergence_max=float(np.max(np.abs(d['divergence'])))))
    out=ROOT/'radial_continuation'
    candidate=dict(original);candidate['id']='ST073-V-order10-radial-extension';candidate['parameters']=dict(original['parameters'],order=10,X_max=3/64);candidate['scope']='Explicit autonomous radial extension; small finite domain only, no global energy or exterior matching';(out/'candidate.json').write_bytes((json.dumps(candidate,indent=2)+'\n').encode())
    report=dict(rows=rows,annulus_X=[1/64,3/64],eta=[-.5,.5],scope='Full local unforced momentum integrated with physical volume Jacobian; quadrature nodes only, no continuum bound, no independent FD repeat, no exterior or global admission',pde_validated=False,global_field_ready=False)
    (out/'annulus.json').write_bytes((json.dumps(report,indent=2)+'\n').encode());print(json.dumps(rows[-3:],indent=2))
if __name__=='__main__':run()
