"""Measured exterior-pressure handoff into the finite nonlinear core.

Each iteration derives pressure from the joined swirl, reconstructs the
finite core under it, then measures the new joined pressure. This is a
finite numerical compatibility experiment, not the paper's collar/moments.
"""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
import numpy as np
from lei_ren_part1_outer_pressure import JoinedOuterPressure, build_reference_joined_profile
from openai_ns_reconstruction.paper_core_series import PaperCoreSeries, ChebyshevEtaGrid
from openai_ns_reconstruction.paper_core_reference import PaperCoreReference
from openai_ns_reconstruction.heat_exterior import HeatExterior


def run(eta_nodes=257):
    joined, reference = build_reference_joined_profile()
    grid = ChebyshevEtaGrid(eta_nodes)
    pressure = np.array([joined.P0(float(z)) for z in grid.eta])
    history = []
    for iteration in range(8):
        core = PaperCoreSeries(reference, maxdegree=4, eta_nodes=eta_nodes,
                               axis_pressure_values=pressure)
        joined = JoinedOuterPressure(lambda R, Z, c=core: float(c.F(R, Z)),
                                     joined.heat)
        updated = np.array([joined.P0(float(z)) for z in grid.eta])
        defect = float(np.max(np.abs(updated-pressure)))
        history.append(dict(iteration=iteration+1, axis_pressure_max_defect=defect))
        if defect < 1e-9:
            break
        pressure = updated
    # Off-grid holdout measures interpolation and reassembly together.
    ztest = np.linspace(-.48, .48, 12)
    holdout = np.array([joined.P0(float(z)) for z in ztest])
    core_axis = core.Pi(0., ztest)
    holdout_defect = float(np.max(np.abs(holdout-core_axis)))
    report = dict(
        candidate_id='lr1-pressure-core-prefix-001',
        source='https://arxiv.org/html/2609.35406v1',
        parameters=reference.metadata()['parameters'], c_inf=joined.heat.c_inf,
        R_core=joined.R_core, R_join=joined.R_join,
        eta_nodes=eta_nodes, radial_degree=4, quadrature_order=64,
        core_axis_pressure_values=core.axis_pressure_values.tolist(),
        history=history, last_collocation_pressure_defect=defect,
        holdout_Z=ztest.tolist(), holdout_outer_P0=holdout.tolist(),
        holdout_core_P0=core_axis.tolist(), holdout_max_pressure_defect=holdout_defect,
        finite_pressure_iteration_converged=bool(defect < 1e-9),
        shared_pressure_holdout_satisfied=bool(holdout_defect < 1e-7),
        pressure_holdout_tolerance=1e-7,
        core_metadata=core.metadata(), pde_validated=False,
        scale_recursion_established=False, full_five_moment_matching=False,
        scope='Finite radial prefix and numerical shared-pressure iteration; independent flat blend, no paper collar, cone, global energy, forcing or oscillatory closure. Holdout defect must be retained.',
    )
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('history', 'holdout_max_pressure_defect', 'finite_pressure_iteration_converged')}))
    return report


def load_core(path=None):
    """Rebuild the saved finite core and joined swirl without re-optimizing."""
    path = Path(path) if path else Path(__file__).with_suffix('.json')
    report = json.loads(path.read_text(encoding='utf-8'))
    core = PaperCoreSeries(PaperCoreReference(**report['parameters']),
                           maxdegree=report['radial_degree'], eta_nodes=report['eta_nodes'],
                           axis_pressure_values=report['core_axis_pressure_values'])
    joined = JoinedOuterPressure(lambda R, Z: float(core.F(R, Z)),
                                 HeatExterior(h=core.reference.h, c_inf=report['c_inf']),
                                 R_core=report['R_core'], R_join=report['R_join'])
    return core, joined


if __name__ == '__main__':
    run()
