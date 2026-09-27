"""Repair instantaneous mean moment/stress constraints for a growth seed.

The seed's positive energy direction is not a time-evolved NS solution.
Pressure and swirl time slopes chosen here leave instantaneous velocity
and its spatial energy matrix unchanged. Direct moment replay is required.
"""
import json
import time
import argparse
import numpy as np

from outer_feedback_evolution import build_current, Trajectory, SwirlValue, solve_control, choose_control
from outer_pressure_modes import OuterPressure
from midplane_outer_slope_patch_repair import P_WINDOWS, P_BREAKS
from midplane_integrated_moment_balance import evaluate
from midplane_outer_residual_source import outer_cones
from affine_momentum import jets, momentum
from radial_continuation import ROOT
from broad_annular_shear import BroadAnnularShear


def apply_mean_increment(current, descriptor):
    if descriptor is None:
        return current
    if descriptor['type'] != 'broad_annular_shear':
        raise ValueError('Unknown mean increment type.')
    return BroadAnnularShear(current, **descriptor['parameters'])


class ControlledMean:
    """Local affine state path in k, preserving the physical-time convention."""
    def __init__(self, current, state, control, k0):
        self.current = current
        self.state = np.asarray(state)
        self.control = np.asarray(control)
        self.k0 = float(k0)
        for name in ('inner', 'nu', 'join_X', 'ratio'):
            setattr(self, name, getattr(current, name))

    def fields(self, points, tau):
        k = -np.log2(2*float(np.asarray(tau).ravel()[0]))
        value = self.state + (k-self.k0)*self.control[9:]
        return OuterPressure(SwirlValue(self.current, value), self.control[:9],
                             windows=P_WINDOWS).fields(points, tau)


def load_saved_field(path=ROOT/'fourier_shear_feasibility.json'):
    """Reload a diagnostic local mean; callers must inspect acceptance data."""
    report = json.loads(path.read_text(encoding='utf-8'))
    if 'control' not in report:
        raise ValueError('Mean-control solve did not produce a field.')
    _, _, current = build_current()
    current = apply_mean_increment(current, report.get('mean_increment'))
    return ControlledMean(current, report['state'], report['control'], report['k']), report


def run(seed_path=ROOT/'fourier_shear_codesign.json', reuse_cache=False):
    seed = json.loads(seed_path.read_text(encoding='utf-8'))
    candidate = seed.get('selected_candidate')
    if candidate is None:
        raise ValueError('No selected growth candidate; do not run mean feasibility.')
    k = float(seed['initial_k'])
    inner, base, current = build_current()
    saved = json.loads((ROOT/'outer_feedback_evolution.json').read_text(encoding='utf-8'))
    original = Trajectory(current, saved['nodes'])
    descriptor = candidate.get('mean_increment')
    current = apply_mean_increment(current, descriptor)
    breaks = list(P_BREAKS)
    if descriptor is not None:
        for ends in (current.onset, current.offset):
            breaks.extend([ends[0],sum(ends)/2,ends[1]])
    breaks = sorted(set(breaks))
    delta = np.asarray(candidate['delta_state'], dtype=float)
    if delta.shape != (9,) or not np.isfinite(delta).all():
        raise ValueError('Expected nine finite swirl value increments.')
    state = original.state(k) + delta
    reference = np.r_[original.pressure(k), original.state(k, 1)]
    report = dict(accepted=False, scale_recursion_established=False,
                  source=seed_path.name, k=k, delta_state=delta.tolist(),
                  state=state.tolist(), mean_increment=descriptor, radial_breaks=breaks,
                  status='solving_controls',
                  scope='Instantaneous mean repair for a supported-wave energy seed. No evolved wave, stress realization, finite-energy or full-domain acceptance.')
    path = ROOT/'fourier_shear_feasibility.json'
    started = time.perf_counter()
    def save():
        report['elapsed_seconds'] = time.perf_counter()-started
        path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    save()
    print(json.dumps(dict(stage='mean_repair_started',k=k)), flush=True)
    cache_path=ROOT/'fourier_shear_control_problem.json'
    cache_identity=dict(k=k,state=state.tolist(),mean_increment=descriptor,
                        radial_breaks=breaks,reference=reference.tolist(),order=48)
    def cache_problem(problem):
        cache_path.write_text(json.dumps(dict(identity=cache_identity,problem=problem),indent=2)+'\n',encoding='utf-8')
    try:
        cache=json.loads(cache_path.read_text(encoding='utf-8')) if reuse_cache and cache_path.exists() else None
        if cache is not None and cache['identity']==cache_identity:
            p=cache['problem']
            En,mn,A,b,ref,feasible=[np.asarray(p[key]) for key in ('En','mn','A','b','reference','feasible')]
            control,selection=choose_control(En,mn,A,b,ref,feasible)
            solve=dict(k=k,state=state.tolist(),control=control.tolist(),
                       moment_max=float(max(abs(np.asarray(p['E'])@control+np.asarray(p['m'])))),
                       minimum_constraint=float(min(A@control+b)),
                       minimum_lambda_squared=min(p['lambdas']),selection_method=selection,
                       problem_cache_reused=True)
        else:
            control, solve = solve_control(inner, base, current, k, state, reference,
                            radial_breaks=breaks,problem_callback=cache_problem)
    except RuntimeError as error:
        report.update(status='candidate_rejected_by_mean_control', reason=str(error))
        save()
        print(json.dumps(dict(stage=report['status'],reason=str(error))), flush=True)
        return report
    report.update(control=control.tolist(), control_solve=solve, status='replaying')
    save()
    field = ControlledMean(current, state, control, k)
    moments = evaluate(field, inner, k, 96, .002, radial_breaks=breaks)
    report['moment_replay'] = dict(order=96, values=moments.tolist(),
                                   max=float(np.max(np.abs(moments))))
    save()
    print(json.dumps(dict(stage='moment_replay',**report['moment_replay'])), flush=True)
    cones = outer_cones(field, k, order=64, radial_breaks=breaks)
    report['outer_cones'] = cones
    report['cone_pass_count'] = sum(row['cone_pass'] for row in cones)
    save()
    tau = .5*2.**-k
    radial = np.linspace(.01,.99,61)
    points = np.concatenate([inner.from_similarity(inner.p.X_max*(1+15*radial)**2,
                             np.full(len(radial),eta), tau) for eta in (-.3,-.1,0.,.1,.3)])
    residual = momentum(jets(field,points,tau,.0005*np.sqrt(field.nu*tau),.0001*tau))
    report['sampled_mean_momentum'] = dict(max=float(np.max(np.linalg.norm(residual,axis=1))),
                                         sample_rms=float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),
                                         point_count=len(points))
    report['status'] = 'instantaneous_repair_replayed'
    report['sampled_outer_mean_constraints_pass'] = bool(report['moment_replay']['max']<1e-3
                                                   and report['cone_pass_count']==len(cones))
    report['constraint_scope']='Four integrated moments and six outer cone samples; excludes the growing-wave-region stress cone.'
    save()
    print(json.dumps({key: report[key] for key in ('status','sampled_outer_mean_constraints_pass',
                      'cone_pass_count','sampled_mean_momentum','elapsed_seconds')}), flush=True)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seed',default='fourier_shear_codesign.json')
    parser.add_argument('--reuse-cache',action='store_true',help='Reuse matching controls only if field/assembly sources are unchanged.')
    args=parser.parse_args()
    run(ROOT/args.seed,reuse_cache=args.reuse_cache)
