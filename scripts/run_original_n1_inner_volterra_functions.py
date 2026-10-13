"""Restore checked raw n1 regular operator and complete-series definitions."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1InnerVolterraFunctions();owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,collar_checked=owner.collar.acceptance_loaded,
        graph_nodes=len(owner.g.nodes),six_raw_component_C3=current.own.encode(owner.functions),
        n1_mean_primitive_Q_V=current.own.encode(owner.recovered),
        operator_rows={chart:{key:len(system[key]) for key in ('B0','B1','C0','g')} for chart,system in owner.systems.items()},
        saved_exact_finite_word_indices=[0,1],full_series_retains_all_words=True,
        common_inner_scope='core 0<=rho<=4; first collar 0<=s<=sc/4; Rin=4/Lambda*exp(hb*sc/4)',
        n1_radial_recovery_weight_is_one_plus_delta=True,
        original_numeric_n1_solution_oracle_executed=False,numerical_original_n1_point_accuracy_claimed=False,
        full_inner_cutoff_and_five_defect_adapter_installed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
