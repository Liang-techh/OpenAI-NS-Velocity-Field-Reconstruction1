"""Inspect the checked original leading Rm patch without rebuilding sources."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN)
        or not receipt.get('same_exact_original_unique_branch_theorem_bound')
        or pins!=raw['input_hashes'] or archive!=current.sha(current.NAME)
        or checker!=current.sha(checker_name) or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked true Rm leading implicit C6/patch source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed true Rm patch input '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Exact source/weight/conjugacy and scoped calculus checks required')
    print(json.dumps(dict(checked=True,true_prescribed_Rm_inlet_C6=True,
        same_original_unique_leading_five_control_functions_C6=True,
        exact_original_beta_integrals_and_raw_normalized_map_conjugacy=True,
        true_Rm_to_Rh_background_fields_and_Ey_Vy_C6=True,full_signed_generic_C5=True,
        independent_P0_and_original_physical_factors_retained=True,
        finite_N_Rm_patch_transport_C5=False,six_outer_chart_C6_continuation=False,
        complete_preRp_pressure_C3=False,original_scale_C6_numeric_containment=False,
        all17_targets_controls=False,global_higher_derivative_frequency=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
