"""Inspect the source-pinned genuine pulse C3 archive without rebuilding."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual pulse pressure C3 source archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed higher pulse source '+name)
    if not all(row['passed'] for row in receipt['checks'].values()):
        raise ValueError('Complete scoped higher pulse checks required')
    print(json.dumps(dict(checked=True,native_and_Gamma_C5=True,
        leading_pulse_U_mean_C5_auxiliary_Q_C4=True,pulse_pressure_C3=True,
        actual_pressure_C3_functions=receipt['checks']['pulse']['genuine_pressure_C3_functions'],
        original_C1_C2_C3_handles_retained=True,original_Gamma_and_pulse_histories_retained=True,
        pressure_derivative_numerical_containment=False,complete_five_defect_C3=False,
        reserved_coefficient_C3=False,required_C3_C5_mixed_norms=False,
        global_n1_closure=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
