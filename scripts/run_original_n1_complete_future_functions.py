"""Load checked complete future C1 functions from saved original ancestors."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1CompleteFutureFunctions();owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        complete_future_C1=current.parent.encode(owner.functions),
        original_small_angular_C1=current.parent.encode(owner.angular),
        actual_pulse_C1=current.parent.encode(owner.pulse_functions),
        actual_positive_pulse_branch_proved=True),indent=2))
