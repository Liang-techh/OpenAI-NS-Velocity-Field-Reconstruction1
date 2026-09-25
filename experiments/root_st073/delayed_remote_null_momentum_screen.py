"""Select a J/S-preserving U null direction by sampled full momentum.

This finite screen uses the smoother X=1.01, width=.02, degree=11
physical lift and compares its 24 canonical quadratic null roots.
"""

import json

import numpy as np

import delayed_remote_u_physical as physical
from coupled_moment_physical_lift import CoupledMomentPhysicalLift
from delayed_momentum_tangent_screen import nodes, residual, stats
from delayed_multimode_cone_fit import BASE_NAME
from delayed_remote_moment_repair import current_mean
from extended_relaxed_cone_screen import profile
from high_frequency_shear_screen import make_field
from moment_shear_slice_repair import moment_vector
from radial_continuation import ROOT


def run():
    physical.START,physical.WIDTH,physical.DEGREE = 1.01,.02,11
    e_source = json.loads((ROOT/'delayed_remote_positive_e.json').read_text())
    floor_source = json.loads((ROOT/'delayed_remote_swirl_floor_capacity.json').read_text())
    e_row = next(entry['result'] for entry in floor_source['rows']
                 if entry['relative_E_floor']==.05)
    base = CoupledMomentPhysicalLift(slice_filename=BASE_NAME)
    mean = current_mean(base)
    target = make_field(16,2.)
    field, solved = physical.solve_coefficients(mean,base,target,
                                                e_source,e_row)
    tau = .5*2**(-5.5)
    train_points,Xtrain,etatrain = nodes(base,(1.015,1.03),(.3,),tau)
    before,_ = residual(mean,train_points,tau)
    options = []
    for option in solved['null_options']:
        field.u_coefficients = np.asarray(option['coefficients'])
        R,div = residual(field,train_points,tau)
        options.append(dict(index=option['index'],sign=option['sign'],
                            coefficients=option['coefficients'],
                            max_abs_corrected_U=(
                                option['max_abs_corrected_U']),
                            norms=np.linalg.norm(R,axis=1).tolist(),
                            momentum=stats(R),
                            max_abs_fd_divergence=float(np.max(np.abs(div)))))
    selected = min(options,key=lambda row:row['momentum']['max'])
    field.u_coefficients = np.asarray(selected['coefficients'])
    screens = []
    for name,xs,etas,t in (
        ('train',(1.015,1.03),(.3,),tau),
        ('space_holdout',(1.02,1.04,1.3),(.28,.32),tau),
        ('time_holdout',(1.02,1.035,1.5),(.29,.31),.5*2**(-5.4))):
        points,X,eta = nodes(base,xs,etas,t)
        old,_ = residual(mean,points,t)
        new,div = residual(field,points,t)
        screens.append(dict(name=name,X=X.tolist(),eta=eta.tolist(),
                            before=stats(old),after=stats(new),
                            max_abs_fd_divergence=float(np.max(np.abs(div)))))
    Xq,w = physical.resolved_quadrature(order=48)
    U,E = profile(field,Xq,.3,tau)
    U0,E0 = profile(target,Xq,.3,tau)
    moment_defect = (moment_vector(U,E,Xq,w)
                     -moment_vector(U0,E0,Xq,w))
    cone_points,_,_ = nodes(base,(1.008,1.016,1.02,1.03),
                            (.2,.25,.3),tau)
    old_local,_ = mean.fields(cone_points,tau)
    new_local,_ = field.fields(cone_points,tau)
    report = dict(source='delayed_remote_u_physical_smooth.json',
                  options=options,selected=dict(index=selected['index'],
                                                sign=selected['sign'],
                                                coefficients=selected['coefficients']),
                  training_before=stats(before),screens=screens,
                  independent_moment_defect=moment_defect.tolist(),
                  max_abs_local_velocity_change=float(np.max(np.abs(
                      new_local-old_local))),
                  scope='Sampled full-momentum selection among 24 '
                        'fixed-slice M/J/S-preserving null roots. '
                        'Moment and PDE checks are finite grids; no '
                        'continuous cone or global acceptance.',
                  accepted=False)
    (ROOT/'delayed_remote_null_momentum_screen.json').write_bytes(
        (json.dumps(report,indent=2)+'\n').encode())
    print(json.dumps(dict(selected=report['selected']['index'],
                          sign=report['selected']['sign'],
                          selected_train=next(s for s in screens
                                              if s['name']=='train'),
                          independent_moment_defect=(
                              report['independent_moment_defect']),
                          max_abs_local_velocity_change=(
                              report['max_abs_local_velocity_change']))),
          flush=True)


if __name__=='__main__':
    run()
