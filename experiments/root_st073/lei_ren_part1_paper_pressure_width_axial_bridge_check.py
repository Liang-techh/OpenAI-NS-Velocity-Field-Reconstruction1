"""Actual-source first-Z jets of the prescribed pressure/width exit."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_adapter import build_source_core
from lei_ren_part1_paper_pressure_width_axial_comparison import AxialPressureWidthComparison
from lei_ren_part1_paper_pressure_width_axial_bridge import AxialPressureWidthExitBridge
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    print('building actual complete-pressure core with second Z jets', flush=True)
    bundle = build_source_core(precision=260, degree=18, j='1e-14', Lambda='1e36',
        logC='5e151', logPstar='14', delta='1e-200', continuous_pressure=True,
        coherent_waiting=True, complete_preheat_components=True, component_Z_jet_depth=2)
    with mp.workdps(bundle['precision']):
        hb = mp.exp(-100-100*mp.mpf('1e152'))
        comparison = AxialPressureWidthComparison(bundle, h_b=hb,
            pressure_order=9, width_order=2, transition_steps=4)
        bridge = AxialPressureWidthExitBridge(comparison, steps=8)
        start = bridge.evaluate(0, '.3')
        print('integrating actual differentiated prescribed exit', flush=True)
        end = bridge.evaluate(2, '.3')
        assert end['Uz_Z'].component(1, 0) != 0
        assert end['P_Z'].component(1, 0) != 0
        increments = {k: end['moments_Z'][k]-start['moments_Z'][k]
                      for k in end['moments_Z']}
        assert all(v.component(0, 1) != 0 for v in increments.values())
        core = comparison.core_snapshot(0, '.3')
        ring_error = lambda a, b: max((abs(a.component(*k)-b.component(*k))/max(
            mp.mpf(1), abs(b.component(*k))) for k in set(a.atoms)|set(b.atoms)), default=mp.mpf(0))
        boundary_errors = {k: ring_error(start[k], core[v].value)
            for k, v in [('F', 'F'), ('Uz', 'Uz'), ('P', 'P'), ('Ur', 'Ur')]}
        assert max(boundary_errors.values()) < mp.mpf('1e-200'), boundary_errors
        encode = lambda jet: {str(k): signed_log(v, 70) for k, v in jet.atoms.items()}
        report = dict(Z='.3', s_end=2, actual_h_b=signed_log(hb, 80),
            metadata=bridge.metadata(), comparison_metadata=comparison.metadata(),
            field_components={k: encode(end[k]) for k in
                ('F', 'F_Z', 'Uz', 'Uz_Z', 'P', 'P_Z', 'Ur', 'Ur_R')},
            first_Z_moment_increments={k: encode(v) for k, v in increments.items()},
            start_core_max_scaled_errors={k: mp.nstr(v, 40) for k, v in boundary_errors.items()},
            ODE_divergence_numerator=encode(end['ODE_divergence_numerator']),
            stress_components={k: encode(end['stress'][k]) for k in
                ('I_theta', 'I_z', 'S_theta', 'S_z', 'T_theta', 'T_z')},
            pressure_tail_Z_atoms_retained=True,
            all_five_first_width_Z_moment_increments_nonzero=True,
            divergence_check_scope='algebraic ODE jet identity, not finite-RK radial differentiation',
            independent_Cartesian_divergence_verified=False,
            functional_terminal_moments_closed=False, finite_energy_certified=False,
            stress_cone_certified=False, temporal_recursion_certified=False)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('actual first-Z exit and radial velocity receipt saved', flush=True)
    return report


if __name__ == '__main__':
    run()
