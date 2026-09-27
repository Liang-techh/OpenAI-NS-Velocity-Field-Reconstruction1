"""Broaden outer modes and include midplane viscous cost in moment repair."""
import argparse
import json
import numpy as np
from scipy.optimize import minimize, least_squares, linprog
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices, outer_moments
from affine_momentum import combine, momentum, jets
from midplane_resolved_feasibility import ZeroBackground, RADIAL_BREAKS, _integrated_moment_coefficients, _evaluate, _jacobian, _cone_replay
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE
from radial_continuation import ROOT
from outer_pressure_modes import OuterPressure,outer_cache
from midplane_connected_cone_repair import cached_cones
from midplane_outer_residual_source import outer_cones
WINDOWS=((.12,.38),(.40,.72),(.58,.92))
BREAKS=sorted({v for lo,hi in WINDOWS for v in (lo,(lo+hi)/2,hi)})


def run(k=11):
    inner,fields=build_fields();base=fields['two_sided_cone']
    old=np.array(json.loads((ROOT/'midplane_connected_cone_edge_repair.json').read_text())['amplitudes']).reshape(3,2,3,2)
    start=np.zeros((3,2,3,3));start[:,:,:,:2]=old;start=start.ravel()
    free=np.array([3,4,5,6,7,8,12,13,14,15,16,17]) + 18 * ((11,15,19).index(k))
    def field(a,background=base):
        return SeparatedMomentModes(background,a,windows=WINDOWS,knots=(11.,15.,19.),axial_powers=(0,1,2))
    current=field(start)
    units=[field(np.eye(54)[i],ZeroBackground(base)) for i in free]
    units += [OuterPressure(ZeroBackground(base),np.eye(6)[i]) for i in range(6)]
    data=moment_slices(inner,base,current,orders=(float(k),),n=48,unit_fields=units,
        unit_fields_are_deltas=True,radial_breaks=BREAKS)[0]
    c=_integrated_moment_coefficients(data)
    scale=np.maximum(abs(c[0]),1.)
    constraints=[dict(type='eq',fun=lambda x:_evaluate(c,x)/scale,
        jac=lambda x:_jacobian(c,x)/scale[:,None])]
    before=momentum(data['baseline']);rscale=float(max(np.linalg.norm(before,axis=1)))
    weight=np.concatenate([p['weights']*p['radii'] for p in data['panels']]);weight/=sum(weight)
    tau=.5*2.**-k
    midpoints=inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.35,.95,101))**2,np.zeros(101),tau)
    args=(midpoints,tau,.0005*np.sqrt(inner.nu*tau),.0001*tau)
    midbase=jets(current,*args);midunits=[jets(unit,*args) for unit in units]
    cost_base=tuple(np.concatenate((old,new)) for old,new in zip(data['baseline'],midbase))
    cost_modes=tuple(np.concatenate((data['modes'][j],np.stack([u[j] for u in midunits])),axis=1) for j in range(3))
    midweight=midpoints[:,0]/sum(midpoints[:,0]);weight=np.r_[.5*weight,.5*midweight]
    def objective(x):
        jet=combine(cost_base,cost_modes,x)
        r=momentum(jet)/rscale
        mu,mg,ml=cost_modes
        dr=(ml+np.einsum('pnij,nj->pni',mg,jet[0])+np.einsum('nij,pnj->pni',jet[1],mu))/rscale
        norm2=np.sum(r*r,axis=1)
        value=np.dot(weight,norm2+.1*norm2**2)
        grad=2*np.einsum('n,ni,pni->p',weight*(1+.2*norm2),r,dr)
        return float(value),grad
    previous=json.loads((ROOT/f'midplane_wide_outer_axial_repair_k{k}.json').read_text())
    x0=np.r_[previous['correction'],np.zeros(6)]
    cache=outer_cache(current,units,k,order=48)
    init_cone=cached_cones(cache,x0)
    def inequalities(x):
        rows=cached_cones(cache,x)
        return np.r_[rows[:,0]/np.maximum(abs(init_cone[:,0]),1e-20)-.02,
            -rows[:,1]/np.maximum(abs(init_cone[:,1]),1e-20)-.02,.95-rows[:,2]]
    constraints.append(dict(type='ineq',fun=lambda x:inequalities(x)[:6]))
    print('Joint velocity-pressure moment and cone caches ready',flush=True)
    direction=np.random.default_rng(73).normal(size=len(units));direction/=np.linalg.norm(direction)
    error=abs((objective(x0+1e-5*direction)[0]-objective(x0-1e-5*direction)[0])/2e-5-objective(x0)[1]@direction)
    assert error < 1e-5*max(1,abs(objective(x0)[1]@direction)),error
    fit=minimize(objective,x0,jac=True,method='SLSQP',bounds=[(-40.,40.)]*len(free)+[(-100.,100.)]*6,
        constraints=constraints,options=dict(maxiter=400,ftol=1e-11))
    shear_fit=fit
    shear_success=bool(fit.success)
    print(json.dumps(dict(stage='shear',success=bool(fit.success),moment_max=float(max(abs(_evaluate(c,fit.x)))),lambda_constraint=float(min(inequalities(fit.x)[:6])))),flush=True)
    fixed=fit.x.copy()
    fixedjet=combine(cache['baseline'],cache['modes'],fixed)
    geometry=[]
    for i in range(6):
        u,J=fixedjet[0][i],fixedjet[1][i];R=cache['panels'][i][3];F=u[1]/R
        shear=np.array([J[1,0]-F,J[2,0]]);N=shear/max(np.linalg.norm(shear),1e-30);K=np.array([-N[1],N[0]])
        lam=-2*F*N[0]*(2*F*N[0]+np.linalg.norm(shear))
        geometry.append((N,K,np.sqrt(max(lam,0)),.95*abs(2*F*N[0])))
    def pressure_vector(p):return np.r_[fixed[:12],p]
    def pressure_inequalities(p):
        rr=momentum(combine(cache['baseline'],cache['modes'],pressure_vector(p)));out=[]
        for i,(sl,r,w,R) in enumerate(cache['panels']):
            T=np.array([-np.dot(w*r*r,rr[sl,1])/R**2,-np.dot(w*r,rr[sl,2])/R])
            N,K,L,G=geometry[i];dn,dk=T@N,T@K
            dnscale=max(abs(init_cone[i,1]),1e-20);conescale=max(G*dnscale,1e-20)
            out.extend([-dn/dnscale-.02,(-G*dn-L*dk)/conescale,(-G*dn+L*dk)/conescale])
        return np.array(out)
    pzero=np.zeros(6);b=pressure_inequalities(pzero)
    A=np.column_stack([pressure_inequalities(v)-b for v in np.eye(6)])
    mzero=_evaluate(c,pressure_vector(pzero))[[1,3]]/scale[[1,3]]
    E=_jacobian(c,pressure_vector(pzero))[[1,3],12:]/scale[[1,3],None]
    lp=linprog(np.zeros(6),A_ub=-A,b_ub=b,A_eq=E,b_eq=-mzero,bounds=[(-100.,100.)]*6,method='highs')
    print(json.dumps(dict(stage='pressure_linear_feasibility',success=bool(lp.success),message=lp.message)),flush=True)
    pressure_success=False
    if lp.success and min(inequalities(fixed)[:6])>=-1e-6:
        def pobj(p):
            value,grad=objective(pressure_vector(p));return value,grad[12:]
        pf=minimize(pobj,lp.x,jac=True,method='SLSQP',bounds=[(-100.,100.)]*6,
            constraints=[dict(type='eq',fun=lambda p:E@p+mzero,jac=lambda p:E),dict(type='ineq',fun=lambda p:A@p+b,jac=lambda p:A)],options=dict(maxiter=300,ftol=1e-11))
        fit.x=pressure_vector(pf.x);pressure_success=bool(pf.success)
    fit.success=bool(shear_success and pressure_success)
    stages=[]
    for label,x,success in [('previous',x0,False),('joint_pressure' if pressure_success else 'shear_only',fit.x,fit.success)]:
        rr=momentum(combine(data['baseline'],data['modes'],x))
        stages.append(dict(label=label,success=bool(success),moment_max=float(max(abs(_evaluate(c,x)))),
            correction_norm=float(np.linalg.norm(x)),momentum_peak=float(max(np.linalg.norm(rr,axis=1))),objective=objective(x)[0],minimum_cone_constraint=float(min(inequalities(x)))))
    print(json.dumps(stages),flush=True)
    candidate=start.copy();candidate[free]+=fit.x[:12];f=OuterPressure(field(candidate),fit.x[12:])
    replay=[]
    for n in (96,128):
        d=moment_slices(inner,base,f,orders=(float(k),),n=n,unit_fields=[f],radial_breaks=BREAKS)[0]
        values=outer_moments(d,np.zeros(1));rr=momentum(d['baseline'])
        replay.append(dict(order=n,moments=values.tolist(),absolute_max=float(max(abs(values))),momentum_peak=float(max(np.linalg.norm(rr,axis=1)))))
        print(json.dumps(replay[-1]),flush=True)
    tau=.5*2.**-k
    points=np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.01,.99,61))**2,np.full(61,eta),tau) for eta in (-.3,-.1,0.,.1,.3)])
    args=(points,tau,.0005*np.sqrt(inner.nu*tau),.0001*tau)
    heldout={}
    components={}
    for label,fld in [('baseline',current),('candidate',f)]:
        rr=momentum(jets(fld,*args));heldout[label]=float(max(np.linalg.norm(rr,axis=1)));components[label]=np.max(abs(rr),axis=0).tolist()
    supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    cone=_cone_replay(inner,f,supports,k,order=96)
    baseline_cone=_cone_replay(inner,current,supports,k,order=96)
    outer=outer_cones(f,k,order=64,radial_breaks=BREAKS)
    report=dict(k=k,shear_optimizer_success=shear_success,pressure_linear_model=dict(A=A.tolist(),b=b.tolist(),E=E.tolist(),mzero=mzero.tolist()),pressure_linear_feasible=bool(lp.success),pressure_linear_message=lp.message,pressure_bounds=[-100,100],windows=WINDOWS,outer_cones=outer,fit_order=48,stages=stages,optimizer_message=str(fit.message),correction=fit.x[:12].tolist(),pressure_coefficients=fit.x[12:].tolist(),amplitudes=candidate.tolist(),minimum_cone_constraint=float(min(inequalities(fit.x))),
        gradient_directional_error=float(error),replay=replay,heldout_momentum_peaks=heldout,heldout_components=components,cone=cone,baseline_cone=baseline_cone,free_indices=free.tolist(),
        scope='Single-scale moment-constrained optimization of sampled full momentum; heldout radial/axial grid. No global maximum, spatial-volume L2, finite energy, or scale recursion established.',accepted=False,scale_recursion_established=False)
    (ROOT/f'midplane_outer_pressure_staged_k{k}.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(heldout=heldout,cone_pass=cone['pass_count'],outer_cone_pass=sum(r['cone_pass'] for r in outer))),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--k',type=int,choices=(11,15,19),default=11)
    run(parser.parse_args().k)
