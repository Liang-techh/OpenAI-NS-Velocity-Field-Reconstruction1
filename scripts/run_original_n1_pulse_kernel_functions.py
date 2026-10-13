"""Load checked exact pulse kernels without reconstructing field ancestors."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1PulseKernelFunctions()
    owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        exact_kernels=current.parent.encode(owner.functions),
        actual_native_incoming_C1=current.parent.encode(owner.incoming),
        complete_future_C1_function_installed=False,
        complete_pulse_pressure_C1_functions_installed=False),indent=2))
