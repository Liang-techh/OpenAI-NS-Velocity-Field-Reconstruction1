"""Restore checked actual lower-collar fields and signed pressure C1."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1LowerCollarFunctions();owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,core_checked=owner.core.acceptance_loaded,
        graph_nodes=len(owner.g.nodes),signed_pressure_C1=current.own.encode(owner.functions),
        leading_fields=current.own.encode({key:owner.fields[key] for key in
            ('Utheta','U','U_y','m','m_y','Q','Q_y','Q_yy','Ur','absolute_swirl_base_pressure')}),
        exact_source_parameters=current.own.encode(owner.parameters),
        radial_scope='Rin=4/Lambda*exp(hb*sc/4) <= R <= r_minus=4/Lambda*exp(hb*sc/2)',
        phase_endpoints_Z_independent=True,all_histories_start_at_true_core_rho4=True,
        core_comparison_not_selected_as_collar=True,numerical_original_collar_point_accuracy_claimed=False,
        complete_inner_cutoff_five_defect_adapter_installed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
