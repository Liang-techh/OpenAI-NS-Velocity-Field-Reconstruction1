"""Actual coupled-target replay; never reset pressure or claim a closed field."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_joined_outer import build_joined_field
from lei_ren_part1_paper_continuous_incoming_outer import ContinuousIncomingOuterField
from lei_ren_part1_paper_coupled_angular_targets import CoupledAngularTargets
from lei_ren_part1_paper_axial_correction import signed_log
from lei_ren_part1_paper_coupled_angular_solve import solve_coupled_angular
from lei_ren_part1_paper_coupled_angular_components import solve_coupled_angular_components


def run(coherent_waiting=False):
    folder=Path(__file__).parent
    print('building actual source for coupled angular targets',flush=True)
    source=build_joined_field(continuous_pressure=True,coherent_waiting=coherent_waiting)
    seeded=json.loads((folder/'lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    atoms=json.loads((folder/'lei_ren_part1_paper_continuous_axial_solve.json').read_text())
    field=ContinuousIncomingOuterField(source,prepared={.3:(seeded,atoms)})
    p=field.outer
    with mp.workdps(field.precision):
        target=CoupledAngularTargets(p).evaluate('.3')
        print('actual angular/pressure targets extracted',flush=True)
        p.pressure_moments_jet(p.schedule.logR_tail,'.3')
        actual=p.pressure_moment_provider.terminal_pressure_jet('.3')
        normalized_actual=actual['required_bump']/target['pressure_scale']
        mu=mp.mpf(str(p.schedule.mu))
        solution=solve_coupled_angular(p.angular_correction_provider,target,smallness_limit=mu**29,require_positive=True)
        components=solve_coupled_angular_components(p.angular_correction_provider,target,smallness_limit=mu**29,require_positive=True)
        encode=lambda v:signed_log(v,80)
        report=dict(component_coupled_solution=components.as_dict(80),coupled_solution=solution.as_dict(80),coupled_rows_closed_at_target_scale=all(abs(v)<mp.mpf('1e-30') for v in solution.target_normalized_row_residuals) if solution.target_normalized_row_residuals is not None else False,direct_forward_preheat_target=encode(target['direct_forward_preheat_target']),waiting_match_residual=encode(source.waiting_match_receipt['normalized_match_residual']) if coherent_waiting else None,Z='.3',coherent_waiting=coherent_waiting,waiting_length=str(p.schedule.waiting_length),actual_targets={k:encode(target[k]) for k in ('r','rZ','s','sZ')},
            angular_target_components={k:encode(v) for k,v in target['r_components'].items()},
            angular_target_Z_components={k:encode(v) for k,v in target['rZ_components'].items()},
            angular_target_over_mu29=encode(abs(target['r'])/mu**29),
            canonical_pressure_target=encode(target['s']),
            fixed_existing_core_pressure_target=encode(normalized_actual),
            canonical_and_existing_core_pressure_targets_equal=normalized_actual==target['s'],
            current_P_infinity=encode(actual['P_infinity']),
            shared_heat_truncation_bound=encode(target['heat_pressure_truncation_bound']),
            current_bump_excluded_from_targets=True,actual_inner_offsets_retained=True,
            pressure_gauge_shift_applied=False,actual_core_uses_complete_preheat_datum=False,
            targets_installed=False,quadrature_error_enclosed=False,
            five_terminal_moments_closed=False,finite_energy_certified=False,
            stress_cone_certified=False,scale_recursion_certified=False)
    output=Path(__file__).with_name('lei_ren_part1_paper_coupled_angular_coherent_waiting.json') if coherent_waiting else Path(__file__).with_suffix('.json')
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('actual_targets','angular_target_over_mu29',
        'canonical_and_existing_core_pressure_targets_equal')},indent=2),flush=True)
    return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--coherent-waiting',action='store_true')
    run(parser.parse_args().coherent_waiting)
