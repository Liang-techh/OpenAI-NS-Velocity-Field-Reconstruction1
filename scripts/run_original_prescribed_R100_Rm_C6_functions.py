"""Inspect the sealed original prescribed downstream C6 source chain."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual prescribed R100-Rm source required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed downstream prescribed source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Scoped source/width/geometry/join/control checks required')
    print(json.dumps(dict(checked=True,prescribed_source_charts=list(current.CHARTS),
        genuine_C6_fields_and_six_own_histories=True,full_Q_inertial_p2_C5=True,
        actual_original_T_and_positive_geometry=True,distinct_width_tags_same_actual_K1_source=True,
        Rm_endpoint_functions_and_physical_y_derivatives=True,
        finite_N_modified_source_C6=False,all17_source_C6=False,
        quantitative_downstream_C6_containment=False,compact_C5_targets_controls=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
