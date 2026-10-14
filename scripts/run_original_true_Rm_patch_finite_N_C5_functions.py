"""Inspect checked finite-N patch functions and pressure without rebuilding."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual Rm patch finite-N C5/pressure C3 source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed finite-N patch source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Actual source/log transport/native bump/pressure checks required')
    print(json.dumps(dict(checked=True,original_Rm_patch_loop_density_correction_C5=True,
        nonzero_Rm_inlet_and_dx_over_x_and_pressure_rate_zero_memory=True,
        total_modified_radial_C5_and_signed_generic_C4=True,
        physical_U_Uy_mean_C5_Q_Qy_C4=True,complete_rminus_to_Rh_signed_pressure_C3=True,
        original_N_global_phase_and_independent_P0_retained=True,
        six_outer_chart_C5_continuation=False,complete_preRp_pressure_C3=False,
        direct_Qyy_and_required_mixed_norms=False,original_scale_C5_numeric_containment=False,
        all17_targets_controls=False,global_higher_derivative_frequency=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
