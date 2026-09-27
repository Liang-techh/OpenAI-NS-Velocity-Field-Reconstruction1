"""Separate outer cone feasibility from pressure moment compensation and bounds."""
import json
import numpy as np
from scipy.optimize import linprog,minimize
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices,outer_moments
from affine_momentum import combine,momentum,jets
from midplane_resolved_feasibility import ZeroBackground,_integrated_moment_coefficients,_cone_replay
from outer_pressure_modes import OuterPressure,outer_cache,WINDOWS,BREAKS
from separated_moment_modes import SeparatedMomentModes
from midplane_outer_residual_source import outer_cones
from radial_continuation import ROOT
from outer_swirl_slope import OuterSwirlSlope
P_WINDOWS=((*WINDOWS[1:],(.93,.99)))
P_BREAKS=sorted(set(BREAKS)|{.93,.96,.99})


def run():
    inner,fields=build_fields();base=fields['two_sided_cone'];k=11
    seed=json.loads((ROOT/'midplane_outer_pressure_staged_k11.json').read_text())
    velocity=SeparatedMomentModes(base,seed['amplitudes'],windows=WINDOWS,knots=(11.,15.,19.),axial_powers=(0,1,2))
    current=OuterPressure(velocity,seed['pressure_coefficients'])
    units=[OuterPressure(ZeroBackground(base),np.eye(9)[i],windows=P_WINDOWS) for i in range(9)]
    units += [OuterSwirlSlope(ZeroBackground(base),np.eye(9)[i],k,P_WINDOWS) for i in range(9)]
    data=moment_slices(inner,base,current,orders=(11.,),n=48,unit_fields=units,unit_fields_are_deltas=True,radial_breaks=P_BREAKS)[0]
    c=_integrated_moment_coefficients(data)
    cache=outer_cache(current,units,k,order=48);res=momentum(cache['baseline']);A=[];b=[];node_cases=[]
    for i,(sl,r,w,R) in enumerate(cache['panels']):
        u,J=cache['baseline'][0][i],cache['baseline'][1][i];F=u[1]/R;s=np.array([J[1,0]-F,J[2,0]]);N=s/np.linalg.norm(s);K=np.array([-N[1],N[0]])
        lam=-2*F*N[0]*(2*F*N[0]+np.linalg.norm(s));assert lam>0
        L=np.sqrt(lam);G=.95*abs(2*F*N[0])
        T=np.array([-np.dot(w*r*r,res[sl,1])/R**2,-np.dot(w*r,res[sl,2])/R])
        P=np.stack([-np.einsum('n,pn->p',w*r*r/R**2,cache['modes'][2][:,sl,1]),-np.einsum('n,pn->p',w*r/R,cache['modes'][2][:,sl,2])])
        scale=max(abs(T@N),1.)
        H=np.stack([-N/scale,(-G*N-L*K)/(G*scale),(-G*N+L*K)/(G*scale)])
        margin=np.array([.02,0.,0.]);A.extend(H@P);b.extend(H@T-margin)
        # Permit an arbitrary axial stress, beyond any chosen pressure basis.
        aa=H[:,1,None];bb=H[:,0]*T[0]-margin
        local=linprog([0.],A_ub=-aa,b_ub=bb,bounds=[(None,None)],method='highs')
        node_cases.append(dict(index=i,eta=(-.2,0.,.2)[i//2],radial_fraction=(.5,.75)[i%2],theta_stress=float(T[0]),normal=N.tolist(),lambda_squared=float(lam),arbitrary_axial_stress_feasible=bool(local.success),status=local.message))
    A=np.array(A);b=np.array(b);E=c[1];m=c[0]
    escale=np.maximum(np.max(abs(E),axis=1),1e-12);E=E/escale[:,None];m=m/escale
    cases=[];chosen=None
    for count,bound,with_moments in ((18,100.,True),(18,None,True)):
        lp=linprog(np.zeros(count),A_ub=-A[:,:count],b_ub=b,A_eq=E[:,:count] if with_moments else None,b_eq=-m if with_moments else None,bounds=[(-bound,bound) if bound is not None else (None,None)]*count,method='highs')
        row=dict(mode_count=count,bound=bound,with_moments=with_moments,success=bool(lp.success),status=lp.message)
        cases.append(row);print(json.dumps(row),flush=True)
        if chosen is None and count==18 and with_moments and lp.success:chosen=lp.x
    report=dict(k=k,pressure_windows=P_WINDOWS,cases=cases,nodewise=node_cases,pressure_linear_model=dict(A=A.tolist(),b=b.tolist(),E=E.tolist(),m=m.tolist()),
        scope='Local swirl-time-slope and pressure correction at fixed instantaneous velocity; moment/cone conditions do not establish interval evolution or recursive contraction.',accepted=False,scale_recursion_established=False)
    if chosen is not None:
        weight=np.concatenate([p['weights']*p['radii'] for p in data['panels']]);weight/=sum(weight)
        r0=momentum(data['baseline']);scale=float(max(np.linalg.norm(r0,axis=1)));D=data['modes'][2]
        def objective(x):
            r=(r0+np.einsum('p,pni->ni',x,D))/scale;n2=np.sum(r*r,axis=1)
            return float(weight@(n2+.1*n2*n2)),2*np.einsum('n,ni,pni->p',weight*(1+.2*n2),r,D/scale)
        normfit=minimize(lambda x:(float(x@x),2*x),chosen,jac=True,method='SLSQP',constraints=[dict(type='eq',fun=lambda x:E@x+m,jac=lambda x:E),dict(type='ineq',fun=lambda x:A@x+b,jac=lambda x:A)],options=dict(maxiter=300,ftol=1e-11))
        fit=minimize(objective,normfit.x,jac=True,method='SLSQP',constraints=[dict(type='eq',fun=lambda x:E@x+m,jac=lambda x:E),dict(type='ineq',fun=lambda x:A@x+b,jac=lambda x:A)],options=dict(maxiter=300,ftol=1e-11))
        f=OuterSwirlSlope(OuterPressure(current,fit.x[:9],windows=P_WINDOWS),fit.x[9:],k,P_WINDOWS);replay=[]
        for n in (96,128):
            d=moment_slices(inner,base,f,orders=(11.,),n=n,unit_fields=[f],radial_breaks=P_BREAKS)[0];v=outer_moments(d,np.zeros(1));replay.append(dict(order=n,moments=v.tolist(),absolute_max=float(max(abs(v)))))
            print(json.dumps(replay[-1]),flush=True)
        supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
        inner_cone=_cone_replay(inner,f,supports,k,order=96)
        tau=.5*2.**-k
        pts=np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.01,.99,61))**2,np.full(61,eta),tau) for eta in (-.3,-.1,0.,.1,.3)])
        rr=momentum(jets(f,pts,tau,.0005*np.sqrt(inner.nu*tau),.0001*tau))
        report.update(inner_cone=inner_cone,heldout_momentum_peak=float(max(np.linalg.norm(rr,axis=1))),coefficients=fit.x.tolist(),optimizer_success=bool(fit.success),replay=replay,outer_cones=outer_cones(f,k,order=64,radial_breaks=P_BREAKS))
    (ROOT/'midplane_outer_slope_pressure_repair.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(nodewise_feasible=[x['arbitrary_axial_stress_feasible'] for x in node_cases])),flush=True)

if __name__=='__main__':run()
