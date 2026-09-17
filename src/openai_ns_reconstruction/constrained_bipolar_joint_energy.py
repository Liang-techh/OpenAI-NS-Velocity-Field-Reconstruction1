"""Joint existing-mode feasibility search for global energy balance, not NS acceptance."""
from dataclasses import replace
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize
from .eq45_supported_delivery import Eq45SupportedDeliveryField
from .constrained_force import RestrictedForce

def run(multistart=False):
    base=Eq45SupportedDeliveryField.load_candidate('artifacts/bipolar_energy/normalized_candidate.json')
    b=base.candidate.parent.profile_basis
    modes=[('phi',i) for i,(_,j) in enumerate(b.mode_indices) if j%2==1]+[('swirl',i) for i,(_,j) in enumerate(b.mode_indices) if j%2==0]
    initial=np.array([getattr(b,k+'_coefficients')[i] for k,i in modes])
    def field(a):
        phi=np.zeros(b.mode_count);swirl=np.zeros(b.mode_count)
        for v,(k,i) in zip(a,modes):
            (phi if k=='phi' else swirl)[i]=v
        return Eq45SupportedDeliveryField(replace(base.candidate,parent=replace(base.candidate.parent,profile_basis=replace(b,phi_coefficients=tuple(phi),swirl_coefficients=tuple(swirl)))))
    fields=[field(v) for v in np.eye(len(modes))]
    def matrices(n,h,t):
        g,w=leggauss(n);r,z=np.meshgrid(g+1,2*g,indexing='ij')
        x=np.column_stack((r.ravel(),0*r.ravel(),z.ravel()));weights=(2*np.pi*r*w[:,None]*2*w[None,:]).ravel()
        u=np.stack([f.at_points(x,t) for f in fields],axis=-1)
        ut=np.stack([(f.at_points(x,t+h)-f.at_points(x,t-h))/(2*h) for f in fields],axis=-1)
        energy=.5*np.einsum('n,nci,ncj->ij',weights,u,u)
        work=np.einsum('n,nci,ncj->ij',weights,u,ut)
        for j in range(3):
            d=np.eye(3)[j]*h
            grad=np.stack([(f.at_points(x+d,t)-f.at_points(x-d,t))/(2*h) for f in fields],axis=-1)
            work+=.01*np.einsum('n,nci,ncj->ij',weights,grad,grad)
        force=np.column_stack([np.einsum('n,nci,nc->i',weights,u,RestrictedForce(a=a,c=c)(x,t)) for a,c in [(1,0),(0,1)]])
        return energy,(work+work.T)/2,force
    times=(.3125,.5,.6875)
    cache=[matrices(32,.005,t) for t in times]
    # Initial energy uses t=.25 exactly; derivative evaluations stay inside window.
    g,w=leggauss(64);r,z=np.meshgrid(g+1,2*g,indexing='ij');x=np.column_stack((r.ravel(),0*r.ravel(),z.ravel()));weights=(2*np.pi*r*w[:,None]*2*w[None,:]).ravel()
    u=np.stack([f.at_points(x,.25) for f in fields],axis=-1);e0=.5*np.einsum('n,nci,ncj->ij',weights,u,u)
    probe=np.array([[.1*np.sqrt(.75),0,.1*.75**.495]])
    core=np.column_stack([f.at_points(probe,.25)[0] for f in fields]);original=core@initial
    def eq(a):
        v=a[:-2];force=a[-2:]
        return np.r_[v@e0@v-1,[v@q@v-v@f@force for _,q,f in cache]]
    def coreineq(a):
        ratio=(core@a[:-2])/original
        return np.r_[ratio-.95,1.05-ratio]
    start=np.r_[initial,0.,0.]
    fit=minimize(lambda a:np.sum((a[:-2]-initial)**2)+.01*np.sum(a[-2:]**2),start,method='SLSQP',bounds=[(-4,4)]*len(modes)+[(0,10)]*2,constraints=[{'type':'eq','fun':eq},{'type':'ineq','fun':coreineq}],options={'maxiter':500,'ftol':1e-10})
    trials=[]
    if multistart:
        rng=np.random.default_rng(9172609)
        starts=[start,fit.x]+[np.r_[np.clip(initial+rng.normal(0,1.5,len(initial)),-4,4),0.,0.] for _ in range(10)]
        feasible=[]
        for j,seed in enumerate(starts):
            result=minimize(lambda a:np.sum(eq(a)[1:]**2),seed,method='SLSQP',bounds=[(-4,4)]*len(modes)+[(0,10)]*2,constraints=[{'type':'eq','fun':lambda a:eq(a)[0]},{'type':'ineq','fun':coreineq}],options={'maxiter':500,'ftol':1e-11})
            physical=abs(eq(result.x)[0])<1e-6 and np.min(coreineq(result.x))>=-1e-6
            trials.append(dict(start=j,optimizer_success=bool(result.success),physical_constraints_satisfied=bool(physical),objective=float(result.fun),initial_energy_error=float(eq(result.x)[0]),parameters=result.x.tolist()))
            if physical:feasible.append(result)
        if not feasible:raise RuntimeError('No energy/core-feasible multistart result; inspect optimization')
        fit=min(feasible,key=lambda result:result.fun)
    a=fit.x;child=field(a[:-2]);rows=[]
    for t in (.34375,.46875,.59375,.71875):
        e,q,f=matrices(96,.0025,t);v=a[:-2]
        rows.append(dict(time=t,energy=float(v@e@v),required_work=float(v@q@v),force_work=float(v@f@a[-2:]),balance_defect=float(v@q@v-v@f@a[-2:])))
    out=Path('artifacts/bipolar_joint_energy_multistart' if multistart else 'artifacts/bipolar_joint_energy');out.mkdir(parents=True,exist_ok=True);child.save_candidate(out/'candidate.json')
    report=dict(multistart_trials=trials,energy_balance_enforced_as_equality=not multistart,parent_sha256=base.sha256,candidate_sha256=child.sha256,modes=modes,parameters=a.tolist(),optimizer_success=bool(fit.success),message=str(fit.message),training_equality_defects=eq(a).tolist(),core_ratios=((core@a[:-2])/original).tolist(),holdout=rows,pde_validated=False,scope='Existing six velocity modes and original bounded two-parameter force. Initial energy1; central component ratios constrained [.95,1.05]. Energy identity only, full momentum unoptimized. Failed optimizer results retained without promotion.')
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':
    import sys
    run(multistart='--multistart' in sys.argv)
