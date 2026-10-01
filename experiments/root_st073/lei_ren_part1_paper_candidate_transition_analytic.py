"""Actual Lambda120 exit continuation and switches through physical R=110.

Center numerical experiment using analytic comparison tangents. Neither
pointwise cone tests nor finite RK refinement proves a uniform cone/closure.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_candidate_comparison_jets import CandidateComparisonJets,HERE,STATE
from lei_ren_part1_paper_candidate_exit_analytic import build_exit
from lei_ren_part1_paper_exit_continuation import ExitContinuation
from lei_ren_part1_paper_exit_switches import ExitSwitches
from lei_ren_part1_paper_exit_field import cone_receipt


def build_transition(calc,exit_steps=16,switch_steps=16):
    bridge,logeps = build_exit(calc,exit_steps)
    continuation = ExitContinuation(bridge,max_R=100)
    def driver(x,z,blend):
        with mp.workdps(calc.precision):
            R = switches.R0*mp.exp(x)
            bar = bridge.comparison.evaluate(switches._global_y(R),z)
            start = switches._start(mp.nstr(z,calc.precision))
            return dict(D=bar['D'],D_Z=bar['D_Z'],I_z=bar['I_z'],I_z_Z=bar['I_z_Z'],
                        Fbar=bar['F'],Fbar_Z=bar['F_Z'],F0=start['F0'],F0_Z=start['F0_Z'],bar=bar)
    switches = ExitSwitches(continuation,steps=switch_steps,analytic_driver_provider=driver)
    return switches,logeps


def signed_log(value):
    return dict(sign=int(mp.sign(value)),log_abs=mp.nstr(mp.log(abs(value)),70) if value else None)


def exact_packet(sample):
    """Persist MP tuples: short decimal logs cannot recover this amplitude."""
    def pack(value):
        if isinstance(value,dict):
            return {k:pack(v) for k,v in value.items()}
        return dict(exact_mpf_tuple=list(mp.mpf(value)._mpf_))
    names = ('R','F','FZ','Uz','UZ','Ur','P','PZ','P0','P0_Z',
             'moments','momentsZ','raw_quadratic_integrals','stress')
    return {name:pack(sample[name]) for name in names if name in sample}


def sample_record(sample):
    F = sample['F']
    def number(v): return mp.nstr(v,70)
    m = sample['moments']
    mz = sample.get('momentsZ',sample.get('moments_Z'))
    stress = sample['stress']
    # Factoring F keeps physically tiny swirl moments and their derivatives
    # readable without discarding their actual nonzero physical values.
    theta = m['theta']/F
    mixed = m['theta_z']/F
    pressure_moment = m['p']/(F*F)
    return dict(R=number(sample['R']),region=sample.get('region','initial_exit'),
                F_nonzero=bool(F>0),F=signed_log(F),Uz=number(sample['Uz']),
                Ur=signed_log(sample.get('Ur',stress['U_r'])),
                moments={name:signed_log(v) for name,v in m.items()},
                moments_Z={name:signed_log(v) for name,v in mz.items()},
                factored_moments=dict(theta_over_F=number(theta),theta_z_over_F=number(mixed),
                                     p_over_F_squared=number(pressure_moment)),
                pressure=signed_log(sample.get('P')),
                stress={name:signed_log(stress[name]) for name in ('S_theta','S_z','T_theta','T_z')},
                pointwise_cone=cone_receipt(sample),
                analytic_driver_used=bool(sample.get('analytic_driver_used',False)))


def run():
    calc = CandidateComparisonJets(32)
    with mp.workdps(calc.precision+40):
        switches,logeps = build_transition(calc)
        z = mp.mpf('.3')
        rows = []
        for radius in ('1','100','101','110'):
            sample = switches.evaluate_R(mp.mpf(radius),z)
            rows.append(sample_record(sample))
            print('Actual candidate transition reached R='+radius,flush=True)
        coarse = switches.evaluate_R(mp.mpf(110),z)
        # Hold the actual exit and comparison fixed; refine only switches.
        fine = ExitSwitches(switches.provider,steps=32,
                            analytic_driver_provider=switches.analytic_driver_provider)
        better = fine.evaluate_R(mp.mpf(110),z)
        refinement = dict(log_F=mp.nstr(abs(mp.log(better['F']/coarse['F'])),70),
                          Uz=mp.nstr(abs(better['Uz']-coarse['Uz']),70),
                          moments={name:signed_log(abs(better['moments'][name]-coarse['moments'][name]))
                                   for name in coarse['moments']})
        if not all(row['F_nonzero'] for row in rows):
            raise AssertionError('Actual candidate transition lost nonzero swirl')
        if not coarse['analytic_driver_used'] and coarse['region'] != 'power':
            # The final power branch has explicit derivatives and needs no hook.
            if not all(row['analytic_driver_used'] for row in rows if row['region'].startswith('switch')):
                raise AssertionError('Switch region used a finite-difference driver')
        names = (STATE,'lei_ren_part1_paper_candidate_comparison_jets.py',
                 'lei_ren_part1_paper_candidate_comparison_adapter.py',
                 'lei_ren_part1_paper_candidate_exit_analytic.py',
                 'lei_ren_part1_paper_exit_tangents.py','lei_ren_part1_paper_exit_continuation.py',
                 'lei_ren_part1_paper_exit_switches.py','lei_ren_part1_paper_mp_stress.py')
        report = dict(center='.3',Lambda='1e120',state_sha256=calc.state_hash,
                      comparison_steps=32,exit_steps=16,switch_steps=16,refinement_switch_steps=32,
                      log_epsilon=mp.nstr(logeps,70),samples=rows,switch_refinement=refinement,
                      R110_refined_endpoint_exact_mpf=exact_packet(better),
                      initial_exit_endpoint_exact_mpf=exact_packet(switches.tangent.evaluate(2*switches.hb,z)),
                      actual_exit_and_switch_moments_generated=True,
                      amplitude_reset_or_pressure_fit=False,axial_finite_difference_used=False,
                      finite_MP_midpoint_input_projection=True,
                      original_parameter_errors_enclosed=False,
                      infinite_core_radial_error_propagated=False,ODE_error_enclosed=False,
                      whole_axis_transition=False,uniform_cone_certified=False,
                      terminal_functional_five_moment_closure=False,temporal_recursion=False,
                      input_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('Actual candidate transition R110 complete; center numerical fields, no terminal matching claim',flush=True)
        return report


if __name__ == '__main__':
    run()
