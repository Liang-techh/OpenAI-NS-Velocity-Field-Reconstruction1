"""Necessary global energy balance under the fixed two-parameter force."""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import lsq_linear
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_force import RestrictedForce

def run():
    f=Eq45SupportedDeliveryField.load_candidate('artifacts/bipolar_energy/normalized_candidate.json')
    rows=[]
    for n,h in [(48,.01),(96,.005)]:
        g,w=leggauss(n);r,z=np.meshgrid(g+1,2*g,indexing='ij')
        x=np.column_stack((r.ravel(),np.zeros(r.size),z.ravel()))
        weights=(2*np.pi*r*w[:,None]*2*w[None,:]).ravel()
        for t in (.3125,.4375,.5625,.6875):
            u=f.at_points(x,t)
            ut=(f.at_points(x,t+h)-f.at_points(x,t-h))/(2*h)
            diss=0.
            for j in range(3):
                d=np.eye(3)[j]*h
                grad=(f.at_points(x+d,t)-f.at_points(x-d,t))/(2*h)
                diss+=.01*np.sum(weights*np.sum(grad**2,axis=1))
            rate=float(np.sum(weights*np.sum(u*ut,axis=1)))
            work=[float(np.sum(weights*np.sum(u*RestrictedForce(a=a,c=c)(x,t),axis=1))) for a,c in [(1,0),(0,1)]]
            rows.append(dict(order=n,step=h,time=t,energy_rate=rate,dissipation=float(diss),required_work=rate+float(diss),force_work_columns=work))
    train=[a for a in rows if a['order']==48]
    fit=lsq_linear(np.array([a['force_work_columns'] for a in train]),np.array([a['required_work'] for a in train]),bounds=(0,10))
    for row in rows:
        row['fitted_work']=float(np.dot(row['force_work_columns'],fit.x))
        row['balance_defect']=row['required_work']-row['fitted_work']
    report=dict(candidate_sha256=f.sha256,force_a=float(fit.x[0]),force_c=float(fit.x[1]),rows=rows,scope='Necessary global energy identity Eprime + nu integral |grad u|² = integral u dot f, assuming divergence-free compact field. Fit coarse quadrature only; finer quadrature reuses force. No full PDE acceptance. Spatial and derivative resolutions varied together; not separate convergence proof.',pde_validated=False)
    p=Path('artifacts/bipolar_energy_balance/report.json');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':run()
