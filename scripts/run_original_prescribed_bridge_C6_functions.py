"""Inspect source-pinned prescribed bridge C6 functions without rebuilding."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or receipt['source_family']!=raw['source_family']):
        raise ValueError('Checked true prescribed bridge C6 source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed prescribed bridge source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Complete scoped source/calculus checks required')
    print(json.dumps(dict(checked=True,true_comparison_C7=True,true_prescribed_background_micro_macro_C6=True,
        true_Q_inertial_p2_C5=True,R100_endpoint_functions=True,
        nonzero_all_six_history_memory_retained=True,
        actual_finite_N_modified_bridge_C6=False,all17_chart_axial6_source=False,
        quantitative_bridge_C6_containment=False,compact_C5_targets_controls=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
