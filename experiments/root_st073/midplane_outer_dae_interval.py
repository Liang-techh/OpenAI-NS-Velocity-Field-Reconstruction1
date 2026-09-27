"""Integrate a short outer-coefficient DAE and replay one callable field."""
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
from scipy.optimize import root
from adaptive_bridge_recursive_defect import build_fields
from adaptive_bridge_moment_fit import moment_slices
from affine_momentum import jets,momentum
from midplane_outer_moment_dae import FREE
from midplane_resolved_feasibility import ZeroBackground,RADIAL_BREAKS,_integrated_moment_coefficients,_evaluate,_jacobian,_cone_replay
from midplane_integrated_moment_balance import evaluate
from separated_moment_modes import SeparatedMomentModes,RADIAL_WINDOWS_THREE
from radial_continuation import ROOT


def run():
    inner,fields=build_fields();base=fields['two_sided_cone']
    a=np.array(json.loads((ROOT/'midplane_outer_axial_interscale.json').read_text())['amplitudes'])
    current=SeparatedMomentModes(base,a,windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
    units=[]
    for index in FREE:
        b=np.zeros(18);b[index]=1.
        units.append(SeparatedMomentModes(ZeroBackground(base),np.tile(b,3),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2)))
    nodes=np.linspace(13.,13.1,5);coeff=[];responses=[]
    for k in nodes:
        d=moment_slices(inner,base,current,orders=(k,),n=48,unit_fields=units,unit_fields_are_deltas=True,radial_breaks=RADIAL_BREAKS)[0]
        coeff.append(_integrated_moment_coefficients(d));response=d['modes'][0]/(d['tau']*np.log(2));rows=[]
        for p in d['panels']:
            sl=p['slice'];r=p['radii'];w=p['weights'];R=p['outer_radius']
            rows.extend([-np.einsum('n,pn->p',w*r*r/R**2,response[:,sl,1]),-np.einsum('n,pn->p',w*r/R,response[:,sl,2])])
        responses.append(np.array(rows));print(json.dumps(dict(stage='cache',k=float(k))),flush=True)
    curves=[CubicSpline(nodes,np.stack([c[i] for c in coeff]),axis=0) for i in range(3)]
    response_curve=CubicSpline(nodes,np.stack(responses),axis=0)
    seed=np.array(json.loads((ROOT/'midplane_outer_moment_dae_k13.json').read_text())['coefficient_values'])
    def coefficients(k):return tuple(curve(k) for curve in curves)
    def project(k,swirl):
        values=seed.copy();values[:6]=swirl;c=coefficients(k)
        def fun(p):
            values[6:8]=p
            return _evaluate(c,values)[[1,3]]
        def jac(p):
            values[6:8]=p
            return _jacobian(c,values)[[1,3],6:8]
        sol=root(fun,seed[6:8],jac=jac,tol=1e-10)
        if max(abs(fun(sol.x)))>1e-7:raise RuntimeError('Axial algebraic projection failed')
        values[6:8]=sol.x
        return values.copy()
    def rhs(k,swirl):
        v=project(k,swirl);m=_evaluate(coefficients(k),v)
        return np.linalg.lstsq(response_curve(k)[[0,2],:6],-m[[0,2]],rcond=None)[0]
    solution=solve_ivp(rhs,(nodes[0],nodes[-1]),seed[:6],rtol=1e-9,atol=1e-11,dense_output=True,max_step=.01)
    if not solution.success:raise RuntimeError(solution.message)
    class PathField:
        def __init__(self):
            for name in ('inner','nu','join_X','ratio'):setattr(self,name,getattr(current,name))
        def fields(self,points,tau):
            k=-np.log2(2*float(np.asarray(tau).ravel()[0]))
            if k<nodes[0] or k>nodes[-1]:raise ValueError('Outside integrated coefficient interval')
            values=project(k,solution.sol(k));b=np.zeros(18);b[FREE]=values
            correction=SeparatedMomentModes(ZeroBackground(base),np.tile(b,3),windows=RADIAL_WINDOWS_THREE,knots=(11.,15.,19.),axial_powers=(0,1,2))
            u,p=current.fields(points,tau);v,q=correction.fields(points,tau)
            return u+v,p+q
    path=PathField();supports=json.loads((ROOT/'midplane_physical_covariance_pairs.json').read_text())
    report=dict(interval=nodes.tolist(),solver_success=bool(solution.success),solver_message=solution.message,nfev=solution.nfev,
        integration_nodes=solution.t.tolist(),swirl_values=solution.y.T.tolist(),initial_values=project(13.,seed[:6]).tolist(),
        final_values=project(13.1,solution.sol(13.1)).tolist(),coefficient_samples=[dict(k=float(k),offset=c[0].tolist(),linear=c[1].tolist(),quadratic=c[2].tolist(),slope_response=D.tolist()) for k,c,D in zip(nodes,coeff,responses)],
        heldout=[],scope='One integrated coefficient trajectory on k in [13,13.1], using interpolated moment maps and axial projection. Heldout field replay is independent of map interpolation. No interval supremum, global energy, full residual acceptance, or scale recursion.',accepted=False,scale_recursion_established=False)
    print(json.dumps(dict(stage='integrated',nfev=solution.nfev,steps=len(solution.t))),flush=True)
    for k in (13.0125,13.0875):
        rows=[]
        for order,zstep in ((48,.002),(96,.001)):
            m=evaluate(path,inner,k,order,zstep)
            rows.append(dict(order=order,z_step_factor=zstep,moments=m.tolist(),maximum=float(max(abs(m)))))
            print(json.dumps(dict(k=k,**rows[-1])),flush=True)
        tau=.5*2.**-k
        points=np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*np.linspace(.01,.99,41))**2,np.full(41,eta),tau) for eta in (-.3,-.1,0.,.1,.3)])
        rr=momentum(jets(path,points,tau,.0005*np.sqrt(inner.nu*tau),.0001*tau))
        cone=_cone_replay(inner,path,supports,k,order=96)
        cone["k"]=float(k)
        report['heldout'].append(dict(k=k,replay=rows,momentum_peak=float(max(np.linalg.norm(rr,axis=1))),cone=cone))
        (ROOT/'midplane_outer_dae_interval.json').write_bytes((json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(cone_passes=[r['cone']['pass_count'] for r in report['heldout']])),flush=True)

if __name__=='__main__':
    run()
