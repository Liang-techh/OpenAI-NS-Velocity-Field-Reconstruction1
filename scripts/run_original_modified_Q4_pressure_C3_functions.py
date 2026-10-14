"""Inspect saved modified generic/radial/pressure functions without rebuilds."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked true modified Q4/radial/pressure C3 source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed modified pressure source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Scoped true source/native/implicit/fast-phase/generic/pressure checks required')
    print(json.dumps(dict(checked=True,true_modified_source_charts=list(current.CHARTS),
        complete_modified_Q_and_signed_inertial_C4=True,
        total_modified_Ey_Vy_mean_y_C5=True,physical_U_mean_Uy_C5_Q_Qy_C4=True,
        original_rminus_to_Rm_signed_m13_m15_pressure_C3=True,
        original_fast_phase_terms_Pstar_measures_and_finite_endpoints_retained=True,
        direct_Qyy_and_N1p_C3=False,complete_preRp_pressure_C3=False,
        numerical_C3_C5_containment=False,all17_finite_N_sources=False,
        compact_C5_targets_controls=False,global_higher_derivative_frequency=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
