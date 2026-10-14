"""Inspect exact transported n1 terminal identities without rebuilding sources."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_terminal_integral_identities as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker=pins.pop(Path(current.__file__).stem+'_check.py',None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(Path(current.__file__).stem+'_check.py')
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual terminal-integral source archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed actual terminal source '+name)
    checks=receipt['checks']
    print(json.dumps(dict(checked=True,actual_terminal_integral_C1_function_zeros=list(current.MOMENTS),
        physical_basis_integral_identities=checks['actual_integrals']['basis_rows'],
        transported_value_Z1_zero_rows=len(checks['transported_partition']['actual_transported_terminal_rows']),
        tail_recovery_zero_rows=len(checks['tail_recovery']['actual_mass_mean_K_QV_and_pressure_tail_rows']),
        mass_mean_K_C1_and_QV_values_zero_after='2Rc',n1_pressure_C1_zero_after='Rv',
        original_axis_pressure_and_actual_quiet_pulse_histories_retained=True,
        leading_angular_and_heat_pressure_not_zeroed=True,
        interior_QV_Z1=False,required_C3_C5_bounds=False,global_n1_closure=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
