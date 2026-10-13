"""Restore complete signed n1 defects and actual C1 linear repair functions."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1FullDefectLinearFunctions();owner.assert_graph()
    receipt=current.parent.linear.read(current.RECEIPT)
    print(json.dumps(dict(checked=owner.acceptance_loaded,cutoff_checked=owner.cutoff.acceptance_loaded,
        linear_operator_checked=owner.linear.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        full_five_signed_C1_functions=current.own.encode(owner.moments),
        selected_actual_C1_coefficients=current.own.encode(owner.coefficients),
        reserved_repair_velocity_C1=current.own.encode(owner.repair_velocity),
        algebraic_linear_five_moment_response_identities_checked=receipt['algebraic_linear_five_moment_response_identities_checked'],
        pressure_partition=owner.partition,coefficient_norms_and_higher_derivative_bounds_transferred=False,
        repaired_own_M_V_P_installed=False,positive_order_five_moment_repair_closed=False,
        original_numeric_n1_solution_oracle_executed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
