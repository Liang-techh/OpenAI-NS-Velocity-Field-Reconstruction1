"""Fixed-poloidal spatial swirl optimization with exact quadratic energy constraint."""
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_bipolar_theta import theta
from .constrained_validation import residual
from .constrained_force import RestrictedForce

def run():
    base=Eq45SupportedDeliveryField.load_candidate('artifacts/bipolar_energy/normalized_candidate.json')
    basis=base.candidate.parent.profile_basis
    indices=[i for i,(_,j) in enumerate(basis.mode_indices) if j%2==0]
    initial=np.array(basis.swirl_coefficients)[indices]
    def field(v):
        coeff=np.zeros(basis.mode_count);coeff[indices]=v
        return Eq45SupportedDeliveryField(replace(base.candidate,parent=replace(base.candidate.parent,profile_basis=replace(basis,swirl_coefficients=tuple(coeff)))))
    fields=[field(v) for v in np.eye(len(indices))]
    pol=field(np.zeros(len(indices)))
    g,w=leggauss(96);r,z=np.meshgrid(g+1,2*g,indexing='ij')
    q=np.column_stack((r.ravel(),0*r.ravel(),z.ravel()));weights=(2*np.pi*r*w[:,None]*2*w[None,:]).ravel()
    up=pol.at_points(q,.25);ep=float(.5*np.sum(weights*np.sum(up**2,axis=1)))
    swirl=np.stack([f.at_points(q,.25)-up for f in fields],axis=-1)
    gram=.5*np.einsum('n,nci,ncj->ij',weights,swirl,swirl)
    probe=np.array([[.1*np.sqrt(.75),0,.1*.75**.495]])
    core=np.array([f.at_points(probe,.25)[0,1] for f in fields]);core0=core@initial
    def target(f,x,t,h):
        return np.sum(residual(f.at_points,lambda x,t:np.zeros(len(x)),lambda x,t:np.zeros_like(x),x,t,nu=.01,step=h)['momentum']*theta(x),axis=1)
    train=np.random.default_rng(9172607).uniform(-2,2,(1024,3)); matrices=[]
    for t in (.3125,.5,.6875):
        columns=[target(f,train,t,.005) for f in fields]
        columns.append(-np.sum(RestrictedForce(a=0,c=1)(train,t)*theta(train),axis=1))
        matrices.append(np.column_stack(columns))
    matrix=np.concatenate(matrices);hess=matrix.T@matrix/len(matrix)
    constraints=[{'type':'eq','fun':lambda a:ep+a[:-1]@gram@a[:-1]-1,'jac':lambda a:np.r_[2*gram@a[:-1],0]},
      {'type':'ineq','fun':lambda a:core@a[:-1]-.95*core0,'jac':lambda a:np.r_[core,0]}]
    fit=minimize(lambda a:a@hess@a,np.r_[initial,0.],jac=lambda a:2*hess@a,method='SLSQP',bounds=[(-4,4)]*len(indices)+[(0,10)],constraints=constraints,options={'maxiter':300,'ftol':1e-12})
    child=field(fit.x[:-1]);c=float(fit.x[-1])
    out=Path('artifacts/bipolar_spatial');out.mkdir(parents=True,exist_ok=True);child.save_candidate(out/'candidate.json')
    x=np.random.default_rng(9172608).uniform(-2,2,(2048,3));rows=[]
    for h in (.01,.005):
        for t in (.25,.3125,.4375,.5625,.6875,.75):
            for name,f in [('base',base),('child',child)]:
                err=target(f,x,t,h)-np.sum(RestrictedForce(a=0,c=c)(x,t)*theta(x),axis=1)
                rows.append(dict(candidate=name,step=h,time=t,sampled_max=float(np.max(np.abs(err))),L2_estimate=float(np.sqrt(64*np.mean(err**2)))))
    report=dict(parent_sha256=base.sha256,candidate_sha256=child.sha256,optimizer_success=bool(fit.success),message=fit.message,coefficients=fit.x[:-1].tolist(),force_c=c,initial_energy=ep+float(fit.x[:-1]@gram@fit.x[:-1]),core_swirl_ratio=float(core@fit.x[:-1]/core0),train_seed=9172607,holdout_seed=9172608,rows=rows,pde_validated=False,scope='Four existing even-eta swirl modes; fixed poloidal field. Initial energy and 95% central swirl floor constrained. Theta-only loss; other momentum/energy balance must be checked. Same fitted force for comparator. No default promotion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2));print(json.dumps(rows[-2:],indent=2))
if __name__=='__main__':run()
