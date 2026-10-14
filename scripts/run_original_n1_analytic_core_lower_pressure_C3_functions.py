"""Inspect source-pinned core C5/lower pressure C3 without rebuilding."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual higher core/lower source archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed higher core/lower source '+name)
    if not all(row['passed'] for row in receipt['checks'].values()):
        raise ValueError('Complete scoped core/lower checks required')
    print(json.dumps(dict(checked=True,actual_analytic_core_required_C5=True,
        actual_lower_U_mean_C5_Q_C4=True,actual_signed_lower_pressure_C3=True,
        source_component_rows=len(receipt['checks']['analytic_source']['actual_admitted_core_component_rows']),
        pressure_rows=receipt['checks']['direct_physical_fixture']['independent_physical_pressure_rows'],
        original_histories_physical_measures_and_m15_sign_retained=True,
        numerical_higher_mixed_core_and_lower_C3_containment=False,
        complete_preRp_pressure_C3=False,complete_five_defect_C3=False,
        reserved_coefficient_C3=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
