"""Independent finite-difference PDE diagnostics; never imports training loss."""
import numpy as np


def residual(velocity, pressure, force, points, time, *, nu=0.01, step=0.005,
             time_bounds=(0.25, 0.75)):
    """Fourth-order Cartesian spatial stencils, second-order time stencils.

    time is scalar. Endpoints use inward one-sided stencils. Returned L2
    estimates require callers to supply appropriate sampling/volume weights.
    """
    x = np.asarray(points, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
        raise ValueError('points must be a finite (N,3) array')
    lo, hi = time_bounds
    if not np.isfinite([time, nu, step, lo, hi]).all() or nu <= 0 or step <= 0:
        raise ValueError('finite positive viscosity and step required')
    if not lo <= time <= hi or 4*step > hi-lo:
        raise ValueError('invalid time domain or stencil width')
    u = velocity(x,time)
    jac = np.empty((len(x),3,3))
    lap = np.zeros_like(u)
    gradp = np.empty_like(u)
    for j in range(3):
        d = np.eye(3)[j]*step
        um2,um1,up1,up2 = [velocity(x+k*d,time) for k in (-2,-1,1,2)]
        jac[:,:,j] = (um2-8*um1+8*up1-up2)/(12*step)
        lap += (-up2+16*up1-30*u+16*um1-um2)/(12*step**2)
        pm2,pm1,pp1,pp2 = [pressure(x+k*d,time) for k in (-2,-1,1,2)]
        gradp[:,j] = (pm2-8*pm1+8*pp1-pp2)/(12*step)
    if time-step < lo:
        ut=(-3*u+4*velocity(x,time+step)-velocity(x,time+2*step))/(2*step)
    elif time+step > hi:
        ut=(3*u-4*velocity(x,time-step)+velocity(x,time-2*step))/(2*step)
    else:
        ut=(velocity(x,time+step)-velocity(x,time-step))/(2*step)
    momentum=ut+np.einsum('nij,nj->ni',jac,u)+gradp-nu*lap-force(x,time)
    return {'momentum':momentum,'divergence':np.trace(jac,axis1=1,axis2=2)}


def sampled_norms(result, volume):
    """Uniform spatial sampling only: Monte Carlo L2, not a certified bound."""
    if not np.isfinite(volume) or volume<=0:
        raise ValueError('positive finite volume required')
    r=np.linalg.norm(result['momentum'],axis=-1)
    d=np.abs(result['divergence'])
    return {'residual_sampled_max':float(r.max()),
            'residual_L2_estimate':float(np.sqrt(volume*np.mean(r*r))),
            'divergence_sampled_max':float(d.max()),
            'divergence_L2_estimate':float(np.sqrt(volume*np.mean(d*d)))}

def structure_metrics(candidate, config):
    """Twenty-one time probes, higher-order energy quadrature; sampled evidence."""
    times=np.linspace(*config['domain']['time_interval'],21)
    tau=1-times
    points=np.column_stack((.1*np.sqrt(tau),np.zeros(len(times)),.1*tau**.495))
    u=candidate.velocity(points,times)
    scaled=u*np.column_stack((np.sqrt(tau),tau**.505,tau**.505))
    reference=np.linalg.norm(scaled[0])
    drift=float(np.max(np.linalg.norm(scaled-scaled[0],axis=1)/reference)) if reference>0 else float('inf')
    energy=[candidate.energy(float(t),96) for t in times]
    n=config['nontriviality']
    signs=bool(np.all(u[:,0]<0) and np.all(u[:,1:]>0))
    passed=(drift<=config['validation']['thresholds']['scaled_core_profile_relative_drift']
            and min(energy)>=n['minimum_energy_each_validation_time']
            and max(energy)<=n['maximum_energy_each_validation_time']
            and abs(energy[0]-n['reference_energy'])<=n['reference_energy_abs_tolerance'] and signs)
    return {'sampled_constraints_pass':bool(passed),'core_drift':drift,
            'energy_range':[min(energy),max(energy)],'core_signs_pass':signs,
            'time_points':21,'energy_quadrature_order':96}


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    from .constrained_candidate import CompactCandidate
    from .constrained_force import RestrictedForce
    parser=argparse.ArgumentParser(description='Independent held-out candidate PDE diagnostics')
    parser.add_argument('--candidate',default='artifacts/constrained/initial_candidate.json')
    parser.add_argument('--config',default='configs/constraints.json')
    parser.add_argument('--training',help='training artifact containing the fitted force coefficients')
    parser.add_argument('--output',default='artifacts/constrained/initial_pde_validation.json')
    args=parser.parse_args()
    cfg=json.loads(Path(args.config).read_text())
    family=json.loads(Path(args.candidate).read_text())['family']
    if family == 'temporal_swirl_v1':
        from .constrained_temporal_swirl import TemporalSwirlCandidate
        c=TemporalSwirlCandidate.load(args.candidate)
    elif family == 'angular_momentum_outer_v1':
        from .constrained_outer_momentum import AngularMomentumCandidate
        c=AngularMomentumCandidate.load(args.candidate)
    elif family == 'compact_axisymmetric_tensor_v1':
        from .constrained_tensor_candidate import TensorCandidate
        c=TensorCandidate.load(args.candidate)
    else:
        c=CompactCandidate.load(args.candidate)
    force_parameters=(json.loads(Path(args.training).read_text())['force']
                      if args.training else cfg['forcing']['initial_parameters'])
    f=RestrictedForce(**force_parameters)
    v=cfg['validation']
    box=np.asarray(cfg['domain']['evaluation_box'],dtype=float)
    volume=float(np.prod(box[:,1]-box[:,0]))
    x=np.random.default_rng(v['seed']).uniform(box[:,0],box[:,1],(v['held_out_points'],3))
    rows=[]
    for t in v['times']:
        for h in v['derivative_steps']:
            rows.append({'time':t,'step':h,**sampled_norms(
                residual(c.velocity,c.pressure,f,x,t,nu=cfg['nu'],step=h,
                         time_bounds=cfg['domain']['time_interval']),volume)})
    fine=[r for r in rows if r['step']==min(v['derivative_steps'])]
    limits=v['thresholds']
    passed=all(r['residual_sampled_max']<=limits['pde_residual_max']
        and r['residual_L2_estimate']<=limits['pde_residual_L2']
        and r['divergence_sampled_max']<=limits['divergence_max']
        and r['divergence_L2_estimate']<=limits['divergence_L2'] for r in fine)
    structure=structure_metrics(c,cfg)
    passed=passed and structure['sampled_constraints_pass']
    result={'structure':structure,'status':'sampled_pde_thresholds_passed' if passed else 'failed_validation',
        'candidate':args.candidate,'force':force_parameters,
        'seed':v['seed'],'points':len(x),'scope':'uniform held-out spatial samples; sampled maxima and Monte Carlo L2; no full acceptance claim','rows':rows}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
