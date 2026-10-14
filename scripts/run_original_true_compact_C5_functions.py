"""Inspect the checked actual compact C5 source without replaying ancestors."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked original compact C5 and 2Rc pressure C3 required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed compact source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Compact source/normalization/implicit/profile/endpoint checks required')
    print(json.dumps(dict(checked=True,actual_all17_normalized_targets_C5=True,
        same_original_unique_compact_branch_C5=True,accepted_C3_target_prefix_and_mu_N_identity=True,
        actual_Rc_to_2Rc_compact_fields_histories_radial_C5=True,
        compact_first_second_velocity_radial_C5=True,modified_signed_generic_C4=True,
        relative_five_terminal_C5_zero_identities=True,signed_rminus_to_2Rc_pressure_C3=True,
        absolute_leading_histories_and_P0_not_zeroed=True,
        original_scale_C5_numeric_containment=False,quiet_2Rc_to_Rp_extension=False,
        complete_preRp_pressure_C3=False,global_higher_derivative_frequency=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
