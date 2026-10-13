"""Load checked original Gamma future functions without ancestor rebuilds."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1GammaFutureFunctions()
    owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        complete_Gamma_future_C1=current.parent.encode(owner.functions),
        exact_nonzero_S=owner.schedule['S'].node,
        complete_future_energy_C1_installed=False,actual_pulse_coefficients_selected=False),indent=2))
