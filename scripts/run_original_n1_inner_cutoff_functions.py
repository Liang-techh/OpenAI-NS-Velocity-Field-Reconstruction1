"""Restore true n1 radial and cutoff functions with signed inner endpoints."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1InnerCutoffFunctions();owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,radial_checked=owner.radial.acceptance_loaded,
        Volterra_checked=owner.inner.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        original_inner_Rin_endpoint_functions=current.own.encode(owner.endpoint),
        even_F_U_mean_K_rho_orders=[0,1,2],even_axial_orders=[0,1,2,3],
        Q_V_rho_orders=[0,1,2],Q_V_axial_orders=[0,1,2],pressure_rho_order=1,
        cutoff_keep_fraction='1/2',cutoff_zero_fraction='3/4',
        only_F1_Uz1_cut=True,mean_K_Q_V_pressure_recovered_from_own_history=True,
        leading_pressure_forcing_unchanged=True,n1_temporal_mean_weight='1+delta',
        common_inner_scope='core 0<=rho<=4; first 0<=s<=sc/4; endpoint Rin',
        original_numeric_n1_solution_oracle_executed=False,
        full_five_defects_and_linear_repair_installed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
