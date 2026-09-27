"""Replay the saved constrained matrices using the equivalent grouped field.

Assembly is never repeated here. The complete residual uses Cartesian FD;
moment/cone replay uses different radial orders. Results remain instantaneous.
"""
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from broad_meridional_constrained import constrained_volume_fit
from broad_meridional_momentum import (build_compact_baseline,build_modes_with_parameters,
    CorrectedField,quadrature_nodes,evaluate_field,RADIAL_QUADRATURE_BREAKS)
from broad_shear_dynamic_control import load_saved_field as load_dynamic_field
from grouped_joined_field import install_in_field
from midplane_integrated_moment_balance import evaluate as replay_moments
from midplane_outer_residual_source import outer_cones

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'broad_meridional_constrained.json'
OUTPUT=ROOT/'broad_meridional_constrained_replay.json'


def load_saved_field(path=OUTPUT):
    report=json.loads(Path(path).read_text())
    dynamic,source=load_dynamic_field()
    install_in_field(dynamic)
    base=build_compact_baseline(dynamic,source)
    modes,_=build_modes_with_parameters(dynamic,source,source['k'],compact_pressure=base)
    return CorrectedField(base,modes,report['coefficients']),report


def run():
    started=time.perf_counter()
    source=json.loads(SOURCE.read_text())
    assembled=source['assembled_problem']
    volume={'columns':np.asarray(assembled['volume_columns']),
            'baseline_residual':np.asarray(assembled['volume_baseline_residual']),
            'weights':np.asarray(assembled['volume_weights'])}
    coefficients,solve=constrained_volume_fit(volume,source['moment_problem'],source['cone_problem'])
    R=volume['baseline_residual']+np.einsum('ncp,p->nc',volume['columns'],coefficients)
    payload={key:source[key] for key in ('assembled_problem','moment_problem','cone_problem')}
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    report=dict(accepted=False,pde_validated=False,scale_recursion_established=False,
        status='selected_cached_candidate',source_matrices=SOURCE.name,
        matrix_content_sha256=digest,k=source['k'],tau=source['tau'],
        coefficients=coefficients.tolist(),
        solver={key:solve[key] for key in ('success','status','message','iterations','used_solver_point',
                'fallback_to_known_feasible','equality_residual_max','minimum_cone_margin')},
        predicted_training=dict(momentum_max=float(np.linalg.norm(R,axis=1).max()),
            momentum_volume_L2=float(np.sqrt(np.sum(volume['weights'][:,None]*R**2)))),
        scope='Cached finite constrained problem, independent full-FD spatial holdout and higher-order integrated moment/cone replay. Grouped JoinedField is an equivalent evaluator, not a reduced physical model. No time integration, continuum, finite-energy, PDE or recursion acceptance.')
    def save():
        report['elapsed_seconds']=time.perf_counter()-started
        OUTPUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    save()
    dynamic,ds=load_dynamic_field()
    report['grouped_backend_replacements']=install_in_field(dynamic)
    baseline=build_compact_baseline(dynamic,ds)
    modes,_=build_modes_with_parameters(dynamic,ds,ds['k'],compact_pressure=baseline)
    candidate=CorrectedField(baseline,modes,coefficients)
    tau=source['tau'];k=source['k']
    points,weights=quadrature_nodes(dynamic.inner,tau,4,4,1,
        radial_bounds=(.01,.99),radial_breaks=RADIAL_QUADRATURE_BREAKS,angle_shift=np.pi/7)
    report['holdout_domain']=dict(eta=[-.5,.5],radial_fraction=[.01,.99],theta=[0.,2*np.pi],
        eta_order=4,radial_order_per_interval=4,point_count=len(points),physical_volume=float(weights.sum()))
    report['holdout']={}
    for label,field in (('raw',dynamic),('corrected',candidate)):
        result=evaluate_field(field,points,tau,weights)
        result.pop('residual_components')
        report['holdout'][label]=result
        save()
    report['status']='full_momentum_replayed';save()
    print(json.dumps({'stage':report['status'],'holdout':report['holdout']}),flush=True)
    breaks=source['quadrature']['radial_split_breaks']
    mom=replay_moments(candidate,dynamic.inner,k,96,.002,radial_breaks=breaks)
    report['moment_replay']=dict(order=96,values=mom.tolist(),maximum_absolute=float(abs(mom).max()))
    report['status']='moments_replayed';save()
    print(json.dumps({'stage':report['status'],'moments':report['moment_replay']}),flush=True)
    locations=[(row['label']['eta'],row['label']['y']) for row in source['cone_problem']['diagnostics']]
    cones=outer_cones(candidate,k,order=64,radial_breaks=breaks,locations=locations)
    report['cone_replay']=dict(order=64,count=len(cones),passes=sum(row['cone_pass'] for row in cones),rows=cones)
    report['status']='completed'
    report['both_holdout_metrics_improve']=all(report['holdout']['corrected'][key]<report['holdout']['raw'][key]
        for key in ('momentum_max','momentum_volume_L2'))
    report['finite_compatibility_pass']=bool(abs(mom).max()<1e-3 and all(row['cone_pass'] for row in cones))
    save()
    print(json.dumps({'stage':'completed','both_metrics_improve':report['both_holdout_metrics_improve'],
        'moment_max':float(abs(mom).max()),'cone_passes':report['cone_replay']['passes'],
        'elapsed_seconds':report['elapsed_seconds']}),flush=True)
    return report


if __name__=='__main__':run()
