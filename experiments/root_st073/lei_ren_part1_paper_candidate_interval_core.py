"""Resumable explicit-center Lambda120 coupled core and analytic tail budget.

Checkpoints hold exact directed interval endpoints. Source identities and
the center interval are checked before resume. This is radial profile
coefficient generation, not temporal n-dependent background recursion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_general_center_factory import GeneralCenterCoreFactory, PARAMETERS
from lei_ren_part1_paper_candidate_gauge_core import _pack_fixed, _unpack_fixed, _pack_rows, _unpack_rows, _atomic_write_json
from lei_ren_part1_paper_functional_core_step import advance_one
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode
from lei_ren_part1_paper_candidate_pressure_function import ACCEPTED_SHA

HERE = Path(__file__).resolve().parent
NORM = 'lei_ren_part1_paper_candidate_core_tail_budget_Lambda120.json'
STEM = 'lei_ren_part1_paper_candidate_interval_core_Z049_Z051'


def input_hashes(factory):
    names = ('lei_ren_part1_paper_candidate_general_center_factory.py',
             'lei_ren_part1_paper_functional_core_step.py',
             'lei_ren_part1_paper_functional_core_recursion.py',
             'lei_ren_part1_paper_candidate_gauge_core.py',
             'lei_ren_part1_paper_candidate_exact_amplitude.py',
             'lei_ren_part1_paper_uniform_axis_jets.py',
             'lei_ren_part1_paper_candidate_pressure_function.py',
             'lei_ren_part1_paper_factored_core_positivity.py',
             'lei_ren_part1_paper_analytic_radial_tail.py', NORM)
    hashes = dict(factory.pressure_hashes)
    hashes.update({name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names})
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return hashes


def tail_budget(ctx, degree):
    norm = json.loads((HERE/NORM).read_bytes())
    if mp.mpf(norm['Lambda']) != mp.mpf(PARAMETERS['Lambda']) or norm['accepted_schedule_sha256'] != ACCEPTED_SHA:
        raise ValueError('Analytic tail datum does not match candidate Lambda/pressure')
    for name, digest in norm['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Analytic norm input changed: '+name)
    def restore(name):
        packet = norm[name]
        return ctx.mpf([mp.make_mpf(tuple(packet['lower_exact_mpf_tuple'])),
                        mp.make_mpf(tuple(packet['upper_exact_mpf_tuple']))])
    h = restore('Xh_parameter')
    phi, psi = restore('analytic_Phi_Xh_norm_upper'), restore('analytic_Psi_Xh_norm_upper')
    rows = []
    worst = ctx.mpf(0)
    for total in range(4):
        for i in range(total+1):
            k = total-i
            if i > degree+1:
                continue
            factor = tail_factor(ctx, degree=degree, radial_order=i, axial_order=k,
                                 radius='4.1', h=h)['tail_per_Xh_norm']
            pb, ub = phi*factor, psi*factor
            upper = max(endpoints(pb)[1], endpoints(ub)[1])
            worst = ctx.mpf([0, max(endpoints(worst)[1], upper)])
            rows.append(dict(scaled_radial_order=i, axial_order=k, Phi_tail=pb, Psi_tail=ub))
    return dict(completed_degree=degree, domain_scaled_R=['0','4.1'], domain_Z=['-1','1'],
                exact_candidate_Taylor_tail_relative_to_accepted_data=True,
                rows=rows, maximum_normalized_mixed_C3_tail=worst,
                target_met=bool(endpoints(worst)[1] <= mp.mpf('1e-12')),
                finite_interval_width_not_an_adapter_error_bound=True,
                physical_derivative_conversion='radial order i multiplies by Lambda**i; Uz correction also by 1/Lambda',
                full_NS_residual_bound=False)


def run(seconds=30, max_steps=8, center=('.49','.51'), target_degree=124, stem=STEM):
    if seconds < 0 or max_steps < 0 or target_degree < 3:
        raise ValueError('Nonnegative time/step budgets and target degree >=3 required')
    factory = GeneralCenterCoreFactory()
    ctx = factory.ctx
    center = list(center) if isinstance(center, tuple) else center
    identity = dict(center=center, target_radial_degree=target_degree, required_axial_depth=3,
                    initial_axis_length=target_degree+4, parameters=PARAMETERS,
                    precision=factory.precision, accepted_schedule_sha256=factory.datum['accepted_schedule_sha256'])
    hashes = input_hashes(factory)
    state_path = HERE/(stem+'_state.json')
    receipt_path = HERE/(stem+'.json')
    with mp.workdps(factory.precision+40):
        if state_path.exists():
            state = json.loads(state_path.read_bytes())
            if state['identity'] != identity or state['input_hashes'] != hashes:
                raise ValueError('Resume center, construction datum or source identity changed')
            fixed = _unpack_fixed(ctx, state['fixed'])
            rows = {name: _unpack_rows(ctx, state['rows'][name]) for name in ('A','Uz','P')}
            completed = state['completed_radial_order']
            if not 0 <= completed <= target_degree or any(
                len(rowset) != completed+1 or any(len(row) != target_degree+4-n for n,row in enumerate(rowset))
                for rowset in rows.values()):
                raise ValueError('Resume row count or retained axial depth is inconsistent')
        else:
            seed = factory.build(center, degree=0, required_depth=target_degree+3)
            fixed = seed['fixed']
            rows = {name: seed['core_rows'][name] for name in ('A','Uz','P')}
            completed = 0
            state = dict(identity=identity, input_hashes=hashes, completed_radial_order=0,
                         fixed=_pack_fixed(fixed), rows={name:_pack_rows(rowset) for name,rowset in rows.items()},
                         step_timings_seconds=[], elapsed_compute_seconds=0,
                         temporal_recursion=False, comparison_transition_generated=False)
            _atomic_write_json(state_path, state)
        print('Explicit-center interval core starting at',completed,'of',target_degree,flush=True)
        started = time.monotonic()
        updates = 0
        while completed < target_degree and updates < max_steps and time.monotonic()-started < seconds:
            before = time.monotonic()
            advance_one(ctx, fixed, rows, completed, ctx.mpf(center), ctx.mpf(PARAMETERS['delta']))
            duration = time.monotonic()-before
            completed += 1
            updates += 1
            state['completed_radial_order'] = completed
            state['rows'] = {name:_pack_rows(rowset) for name,rowset in rows.items()}
            state['step_timings_seconds'].append(duration)
            state['elapsed_compute_seconds'] += duration
            _atomic_write_json(state_path, state)
            print('Explicit-center radial order',completed,'of',target_degree,'seconds',round(duration,3),flush=True)
        budget = tail_budget(ctx, completed) if completed >= 3 else None
        report = dict(identity=identity, completed_radial_order=completed,
                      completed_target=completed==target_degree,
                      input_hashes=hashes, state_file=state_path.name,
                      state_sha256=hashlib.sha256(state_path.read_bytes()).hexdigest(),
                      exact_interval_checkpoint=True, original_center_state_reused=False,
                      last_run_updates=updates, analytic_tail_budget=budget,
                      pressure_datum_preserved=True, interval_family_coefficients=True,
                      infinite_tail_and_finite_interval_width_are_separate=True,
                      whole_axis_atlas=False, comparison_transition_generated=False,
                      temporal_recursion=False, independent_Cartesian_residual_validated=False)
        receipt_path.write_text(json.dumps(_encode(report),indent=2)+'\n',encoding='utf-8')
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seconds', type=float, default=30)
    parser.add_argument('--max-steps', type=int, default=8)
    args = parser.parse_args()
    run(seconds=args.seconds, max_steps=args.max_steps)
