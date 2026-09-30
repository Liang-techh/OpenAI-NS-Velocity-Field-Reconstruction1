"""Shared-pressure iteration for the actual angular-matched exterior."""
import json
from pathlib import Path

import numpy as np

from lei_ren_part1_extended_swirl import ExtendedSwirl
from lei_ren_part1_pressure_core import load_core
from lei_ren_part1_heat_collar import HeatCollar
from openai_ns_reconstruction.paper_core_series import PaperCoreSeries, ChebyshevEtaGrid
from openai_ns_reconstruction.paper_core_reference import PaperCoreReference


def run(eta_nodes=257,max_iterations=8):
    saved,_=load_core()
    grid=ChebyshevEtaGrid(eta_nodes)
    pressure=np.array(saved.Pi(0,grid.eta))
    history=[]
    for iteration in range(max_iterations):
        core=PaperCoreSeries(saved.reference,maxdegree=4,eta_nodes=eta_nodes,
                             axis_pressure_values=pressure)
        profile=ExtendedSwirl(core=core)
        updated=np.array([profile.axis_pressure(float(z)) for z in grid.eta])
        defect=float(np.max(np.abs(updated-pressure)))
        row={"iteration":iteration+1,"axis_pressure_max_defect":defect,
             "axis_pressure_min":float(updated.min()),"axis_pressure_max":float(updated.max())}
        history.append(row)
        print(json.dumps(row),flush=True)
        if defect<1e-9: break
        pressure=updated
    ztest=np.linspace(-.975,.975,16)
    holdout=np.array([profile.axis_pressure(float(z),n=256) for z in ztest])
    actual_core=np.array(core.Pi(0,ztest))
    holdout_defect=float(np.max(np.abs(holdout-actual_core)))
    report={"source":"https://arxiv.org/html/2609.35406v1",
            "candidate_id":"lr1-extended-angular-pressure-prefix-001",
            "eta_nodes":eta_nodes,"radial_degree":4,"parameters":saved.reference.metadata()['parameters'],
            "core_axis_pressure_values":core.axis_pressure_values.tolist(),
            "history":history,"last_collocation_pressure_defect":defect,
            "holdout_Z":ztest.tolist(),"holdout_outer_P0":holdout.tolist(),
            "holdout_core_P0":actual_core.tolist(),"holdout_max_pressure_defect":holdout_defect,
            "finite_pressure_iteration_converged":bool(defect<1e-9),
            "shared_pressure_holdout_satisfied":bool(holdout_defect<1e-7),
            "collar":profile.collar.metadata(),"R_core":profile.R_core,
            "R_anchor":profile.R_anchor,"anchor_power":profile.anchor_power,
            "whole_cone_validated":False,"five_moments_closed":False,"recursion_validated":False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({"holdout_max_pressure_defect":holdout_defect,
                     "finite_pressure_iteration_converged":report['finite_pressure_iteration_converged']}),flush=True)
    return report


def load_extended_profile(path=None, *, require_pressure_holdout=True):
    """Replay the saved actual core/outer profile without pressure iteration."""
    path=Path(path) if path else Path(__file__).with_suffix('.json')
    report=json.loads(path.read_text(encoding='utf-8'))
    if require_pressure_holdout and not (report['finite_pressure_iteration_converged']
                                        and report['shared_pressure_holdout_satisfied']):
        raise ValueError('saved extended pressure profile has not passed its holdout')
    core=PaperCoreSeries(PaperCoreReference(**report['parameters']),
                         maxdegree=report['radial_degree'],eta_nodes=report['eta_nodes'],
                         axis_pressure_values=report['core_axis_pressure_values'])
    params=report['collar']
    collar=HeatCollar(h=params['h'],c_inf=params['c_inf'],R_b=params['R_b'],
                     ell=params['ell'],epsilon=params['epsilon'])
    return ExtendedSwirl(core=core,collar=collar,R_core=report['R_core'],
                         R_anchor=report['R_anchor'],anchor_power=report['anchor_power'])


if __name__=='__main__': run()
