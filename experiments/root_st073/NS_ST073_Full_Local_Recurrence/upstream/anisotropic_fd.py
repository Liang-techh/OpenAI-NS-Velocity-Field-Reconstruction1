"""Independent physical Cartesian finite differences with directional scales.
Coordinates are fixed in time differentiation. Step choices use preexisting
profile scales, not the residual at a sample. Field coefficients never change.
"""
import sys,json
from pathlib import Path
import numpy as np
from normalized_core import NormalizedCore
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'upstream'))
from local_field import LocalField

def fd(field,pts,tau,steps,time_step):
    u,p=field.fields(pts,tau);grad=np.zeros((len(pts),3,3));lap=np.zeros_like(u);gp=np.zeros_like(u)
    for j in range(3):
        hh=steps[:,j];dis=np.zeros_like(pts);dis[:,j]=hh
        m2,pm2=field.fields(pts-2*dis,tau);m,pm=field.fields(pts-dis,tau)
        p1,pp1=field.fields(pts+dis,tau);p2,pp2=field.fields(pts+2*dis,tau)
        grad[:,:,j]=(m2-8*m+8*p1-p2)/(12*hh[:,None]);gp[:,j]=(pm2-8*pm+8*pp1-pp2)/(12*hh)
        lap+=(-p2+16*p1-30*u+16*m-m2)/(12*hh[:,None]**2)
    h=time_step;ut=(field.fields(pts,tau+2*h)[0]-8*field.fields(pts,tau+h)[0]+8*field.fields(pts,tau-h)[0]-field.fields(pts,tau-2*h)[0])/(12*h)
    return ut+np.einsum('nij,nj->ni',grad,u)+gp-field.nu*lap,np.trace(grad,axis1=1,axis2=2)

def run():
    c=NormalizedCore.load(ROOT/'data/ST072-N.json');f=LocalField(c);p=c.params;w=p.sigma/np.sqrt(p.lam*4.495)
    eta=np.array([c.eta0,c.eta0-1.1*w,c.eta0+.85*w,-.18,.22]);Y=np.array([2.,1.8,2.3,1.5,2.7]);X=Y/p.lam;res=[]
    for k in [.4,3,5.5]:
        tau=.5*2**(-k);q=tau/(1-eta**2);pts=f.from_similarity(X,eta,tau);ref=f.evaluate(pts,tau)['residual']
        # x/y scales follow radial geometry, z follows local axial sharpness.
        radial=np.sqrt(2*.01*q/p.lam)
        axial=np.sqrt(.01)*q**c.D*np.maximum(w,.1*abs(eta-c.eta0))
        rows=[]
        for fac in [.06,.03,.015]:
            steps=np.column_stack((radial,radial,axial))*fac
            R,div=fd(f,pts,tau,steps,tau*2e-5);err=np.linalg.norm(R-ref,axis=1)/np.maximum(1,np.linalg.norm(ref,axis=1))
            rows.append({'step_factor':fac,'steps':steps.tolist(),'absolute_errors':np.linalg.norm(R-ref,axis=1).tolist(),'relative_errors':err.tolist(),'relative_error_max':float(err.max()),'divergence_max':float(abs(div).max())})
        res.append({'k':k,'rows':rows})
    out={'scope':'Directional Cartesian FD at fixed points, separately scaled radial and axial steps; complements rather than overwrites initial isotropic FD records','results':res}
    pth=ROOT/'evidence/anisotropic_fd.json'
    if pth.exists():raise FileExistsError(pth)
    pth.write_text(json.dumps(out,indent=2));return out
if __name__=='__main__':print(json.dumps(run(),indent=2))
