"""Inspect the checked six-chart Rc continuation without rebuilding sources."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked original six-chart C6/C5 and Rc pressure C3 required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed outer source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Original source theorem, geometry and scoped calculus checks required')
    print(json.dumps(dict(checked=True,six_outer_prescribed_C6_generic_C5=True,
        all17_prescribed_C6_generic_C5_and_finite_N_C5_functions=True,
        same_original_Rh_leading_function_identity_and_moment_ODE=True,
        nonzero_leading_correction_and_pressure_memories_retained=True,
        original_O2_exponential_Jacobian_and_mu_N_global_phase_P0=True,
        total_modified_radial_C5_and_signed_generic_C4=True,
        signed_rminus_to_Rc_pressure_C3=True,
        compact_Rc_to_2Rc_C5=False,quiet_2Rc_to_Rp_extension=False,
        complete_preRp_pressure_C3=False,original_scale_high_order_numeric_containment=False,
        direct_Qyy_required_mixed_norms=False,all17_normalized_targets_and_compact_controls=False,
        global_higher_derivative_frequency=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
