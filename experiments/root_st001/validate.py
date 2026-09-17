"""Independent Cartesian FD validation; does not call training derivative bundles."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from spacetime import Family, force, quad, TIMES, NU


def time_derivative(field,points,t,h):
    if not .25 <= t <= .75 or h<=0 or 4*h>.5:
        raise ValueError('Invalid time stencil')
    if t-2*h >= .25 and t+2*h <= .75:
        offsets=(-2,-1,1,2);coefs=(1,-8,8,-1)
    elif t+4*h <= .75:
        offsets=(0,1,2,3,4);coefs=(-25,48,-36,16,-3)
    elif t-4*h >= .25:
        offsets=(0,-1,-2,-3,-4);coefs=(25,-48,36,-16,3)
    else:
        raise ValueError('Stencil exits declared time interval')
    return sum(c*field(points,t+o*h)[0] for o,c in zip(offsets,coefs))/(12*h)


def cartesian_residual(field,forcing,points,t,space_step=.01,time_step=.0025):
    x=np.asarray(points,float);h=space_step
    if x.ndim!=2 or x.shape[1]!=3 or not np.all(np.isfinite(x)) or h<=0:
        raise ValueError('Invalid Cartesian sample')
    u,p=field(x,t);jac=np.empty((len(x),3,3));lap=np.zeros_like(u);gp=np.zeros_like(u)
    for j in range(3):
        e=np.zeros(3);e[j]=h
        um2,pm2=field(x-2*e,t);um1,pm1=field(x-e,t)
        up1,pp1=field(x+e,t);up2,pp2=field(x+2*e,t)
        jac[:,:,j]=(um2-8*um1+8*up1-up2)/(12*h)
        lap+=(-um2+16*um1-30*u+16*up1-up2)/(12*h*h)
        gp[:,j]=(pm2-8*pm1+8*pp1-pp2)/(12*h)
    ut=time_derivative(field,x,t,time_step)
    residual=ut+np.einsum('nij,nj->ni',jac,u)+gp-NU*lap-forcing(x,t)
    return residual,np.trace(jac,axis1=1,axis2=2)


def norms(R,div):
    v=np.linalg.norm(R,axis=1)
    return dict(momentum_max=float(v.max()),momentum_L2=float(np.sqrt(64*np.mean(v*v))),
      divergence_max=float(np.abs(div).max()),divergence_L2=float(np.sqrt(64*np.mean(div*div))))


def validate(candidate,out,seed=9172622):
    family,raw=Family.load(candidate);fc=family.coefficients(raw)[3]
    field=lambda x,t:family.fields(raw,x,t)
    forcing=lambda x,t:force(x,t,*fc)
    rng=np.random.default_rng(seed);x=rng.uniform(-2,2,(4096,3))
    rows=[];temporal=[]
    for h in (.02,.01,.005):
        for t in TIMES:
            R,div=cartesian_residual(field,forcing,x,t,h,.0025)
            row=dict(time=float(t),space_step=h,time_step=.0025,**norms(R,div))
            # Independent formulation cross-check; this is NOT used by the FD operator.
            exact=family.analytic_residual(raw,x,t)
            row['analytic_FD_max_difference']=float(np.max(np.linalg.norm(exact-R,axis=1)))
            rows.append(row)
        print(json.dumps({'finished_spatial_step':h,'worst_max':max(a['momentum_max'] for a in rows if a['space_step']==h)}),flush=True)
    # Temporal refinement on independent subset to vary one parameter at a time.
    for h in (.02,.01,.005):
        for t in (.25,.4375,.75):
            R,div=cartesian_residual(field,forcing,x[:256],t,.005,h)
            temporal.append(dict(time=t,space_step=.005,time_step=h,**norms(R,div)))
    energies=[];moments=[]
    for order in (24,48,96):
        s,z,w=quad(order);pts=np.column_stack((np.sqrt(s),np.zeros(len(s)),z))
        for t in TIMES:
            u,p=field(pts,t);E=.5*float(w@np.sum(u*u,axis=1))
            energies.append(dict(order=order,time=float(t),energy=E))
            # Exact weak identity under compact u,p and divergence-free prescribed f:
            # <(x/2,y/2,-z),R> = integral(uz^2-(ux^2+uy^2)/2).
            D=float(w@(u[:,2]**2-.5*(u[:,0]**2+u[:,1]**2)))
            moments.append(dict(order=order,time=float(t),anisotropy=D,
              momentum_L2_lower_bound_estimate=abs(D)/np.sqrt(88*np.pi/3),
              scope='Quadrature estimate of a pressure-independent continuum lower bound; not a rigorous interval bound or sampled max'))
    tau=1-TIMES;points=np.column_stack((.1*np.sqrt(tau),np.zeros(len(tau)),.1*tau**.495))
    u=field(points,TIMES)[0];scaled=u*np.column_stack((np.sqrt(tau),tau**.505,tau**.505))
    drifts=np.linalg.norm(scaled-scaled[0],axis=1)/np.linalg.norm(scaled[0])
    boundary=np.array([[2,0,0],[0,2,0],[0,0,2],[0,0,-2],[2.01,0,0],[0,0,2.01],[3,3,3]],float)
    bound=max(float(np.max(np.abs(v))) for t in TIMES for v in field(boundary,t))
    fine=[r for r in rows if r['space_step']==.005]
    ef=[r['energy'] for r in energies if r['order']==96]
    e48=[r['energy'] for r in energies if r['order']==48]
    econv=max(abs(a-b)/max(abs(b),1e-30) for a,b in zip(e48,ef))
    amp=family.coefficients(raw)[4]
    gates=dict(momentum_max=all(r['momentum_max']<.001 for r in fine),momentum_L2=all(r['momentum_L2']<.001 for r in fine),divergence_max=all(r['divergence_max']<1e-5 for r in fine),divergence_L2=all(r['divergence_L2']<1e-5 for r in fine),initial_energy=abs(ef[0]-1)<=.001,energy_range=all(.1<=v<=10 for v in ef),energy_quadrature=econv<.001,boundary=bound<=1e-10,core_signs=bool(np.all(u[:,0]<0)&np.all(u[:,1]>0)&np.all(u[:,2]>0)),scaled_core_drift=bool(max(drifts)<.05),derived_amplitude=bool(1e-4<=amp<=100))
    result=dict(experiment='CR-ROOT-ST002' if seed==9172624 else 'CR-ROOT-ST001',candidate_sha256=hashlib.sha256(Path(candidate).read_bytes()).hexdigest(),validation_seed=seed,uniform_box_points=4096,force={'a':float(fc[0]),'c':float(fc[1])},spatial_refinement=rows,temporal_refinement_subset=temporal,energy_quadrature=energies,pressure_independent_moments=moments,core={'scaled_vectors':scaled.tolist(),'relative_drifts':drifts.tolist()},boundary_max=bound,energy_relative_change_48_to_96=econv,derived_amplitude=float(amp),gates=gates,pde_validated=False,all_numeric_gates_pass=all(gates.values()),scope='Finite-window alternative-family research. Sampled numerical gates are not a continuum supremum proof, exact paper reconstruction or blow-up result. No default candidate promoted.')
    Path(out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'gates':gates,'worst_momentum_max':max(r['momentum_max'] for r in fine),'worst_momentum_L2':max(r['momentum_L2'] for r in fine),'core_drift':max(drifts)},indent=2),flush=True)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('candidate');p.add_argument('--out',default='validation.json');p.add_argument('--seed',type=int,default=9172622);a=p.parse_args();validate(a.candidate,a.out,a.seed)
