"""Nominal actual-candidate exit bridge using analytic comparison derivatives.

Directed input intervals are projected to MP midpoints for the existing
finite exit integrator. This experiment is not a certified bridge/cone.
"""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp

from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets, STATE, HERE, restore
from lei_ren_part1_paper_core_adapter import CorePolynomial
from lei_ren_part1_paper_exit_tangents import ExitTangents
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets


def midpoint(value):
    lo,hi = endpoints(value)
    return (lo+hi)/2


def candidate_core(calc):
    state = json.loads((HERE/STATE).read_bytes())
    c = calc.ctx
    coefficients = [calc.F0]
    for n in range(1,4):
        coefficients.append(sum((calc.ell[k]*coefficients[n-1-k]
                                 for k in range(n)),c.mpf(0))/n)
    amplitude = IntervalTaylor(c,coefficients)
    frows = [[midpoint(v) for v in (p*amplitude*calc.lam**n).coefficients]
             for n,p in enumerate(calc.phi)]
    urows = [[midpoint(v) for v in (p*calc.lam**n).coefficients]
             for n,p in enumerate(calc.u)]
    prows = [[midpoint(restore(c,v)) for v in row[:4]] for row in state['P_rows']]
    packet = dict(F=frows,Uz=urows,P=prows,Z=calc.center,
                  delta=mp.mpf('1e-200'),precision=calc.precision)
    def factory(Z):
        if abs(mp.mpf(str(Z))-calc.center) > mp.mpf(10)**(-calc.precision+5):
            raise ValueError('Candidate tensor only supplies Z=.3')
        return packet
    return CorePolynomial(factory,Lambda='1e120',delta='1e-200',precision=calc.precision)


def build_exit(calc,steps=8):
    core = candidate_core(calc)
    comparison = SimpleNamespace(core=core,precision=calc.precision,h_b=calc.hb)
    bounds = uniform_axis_jets(calc.ctx,radius=1,j='1e-14',Lambda='1e120',
                              logC='5e151',delta='1e-200',length=4)
    upper = endpoints(bounds['G_absolute_upper'])[1]
    logeps = -mp.mpf('1e120')*(2*upper+1)
    bridge = ExitTangents(comparison,steps=steps,epsilon=mp.exp(logeps),
        analytic_driver_provider=lambda y,z:provider(y,z))
    def provider(y,z):
        packet = calc.exit_driver_jets(y,z,bridge.multiplier(y))
        return dict(A=midpoint(packet['A'][0]),B=midpoint(packet['B'][0]),
                    A_Z=midpoint(packet['A_Z']),B_Z=midpoint(packet['B_Z']))
    return bridge, logeps


def run():
    calc = CandidateComparisonJets(32)
    with mp.workdps(calc.precision+40):
        bridge,logeps = build_exit(calc,8)
        out = bridge.evaluate('.01','.3')
        refined,_ = build_exit(calc,16)
        fine = refined.evaluate('.01','.3')
        if not out['analytic_driver_used'] or out['F'] <= 0:
            raise AssertionError('Actual candidate bridge failed analytic path')
        def nstr(v): return mp.nstr(v,60)
        report = dict(state_sha256=calc.state_hash,actual_candidate_core_used=True,
            midpoint_projection_of_directed_inputs=True,
            comparison_steps=32,exit_steps=8,refinement_exit_steps=16,
            epsilon_policy='conservative uniform analytic G bound; cone not certified',
            log_epsilon=nstr(logeps),center='.3',endpoint_y='.01',
            analytic_driver_used=True,axial_finite_difference_used=False,
            F_nonzero=True,log_abs_F=nstr(mp.log(out['F'])),
            endpoint_Uz=nstr(out['Uz']),
            moment_signed_logs={name:dict(sign=int(mp.sign(v)),
                log_abs=nstr(mp.log(abs(v))) if v else None)
                for name,v in out['moments'].items()},
            exit_step_refinement={name:nstr(abs(out[name]-fine[name]))
                                 for name in ('log_F_over_Fa','Uz','F_R_over_F')},
            infinite_radial_remainder_propagated=False,
            ODE_discretization_error_enclosed=False,stress_cone_certified=False,
            full_transition_R110_generated=False,terminal_five_moment_closure=False,
            input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                for name in (STATE,'lei_ren_part1_paper_candidate_comparison_jets.py',
                    'lei_ren_part1_paper_exit_tangents.py','lei_ren_part1_paper_exit_bridge.py',
                    'lei_ren_part1_paper_core_adapter.py','lei_ren_part1_paper_uniform_axis_jets.py')},
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Actual candidate initial exit bridge generated with analytic derivatives; numerical center result only',flush=True)
        return report


if __name__ == '__main__':
    run()
