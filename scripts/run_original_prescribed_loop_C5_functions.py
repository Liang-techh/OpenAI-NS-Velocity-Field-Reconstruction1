"""Inspect saved true prescribed-loop C5 sources without rebuilding ancestors."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked true prescribed-loop C5 source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prescribed-loop source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Scoped source/implicit calculus/density/transport checks required')
    print(json.dumps(dict(checked=True,true_source_charts=list(current.CHARTS),
        signed_q_inverse_A_B_C5=True,phase_primitives_C5=True,
        finite_N_velocity_and_five_signed_densities_C5=True,
        original_rminus_to_Rm_complete_correction_histories_C5=True,
        modified_background_E_V_m_h_k_e_p_P0_pressure_C5=True,
        original_common_repair_only_N_and_pressure_memory_retained=True,
        numerical_C5_containment=False,mixed_y_C5_bounds=False,
        all17_finite_N_C5_sources=False,compact_C5_targets_controls=False,
        global_higher_derivative_frequency_admission=False,global_n1=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
