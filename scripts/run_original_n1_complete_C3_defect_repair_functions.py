"""Inspect the saved complete C3 partition without rebuilding ancestors."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_complete_C3_defect_repair_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=rec['source_family'] or not all(c['passed'] for c in rec['checks'].values())):
        raise ValueError('Checked complete signed C3 partition and selected inverse required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed complete C3 source '+name)
    print(json.dumps(dict(checked=True,existing_inner_C3_endpoints_reused=True,
        complete_signed_partition_C3=True,selected_reserved_coefficients_C3=True,
        algebraic_physical_five_response_C3_zero=True,repair_Uz_Utheta_F_C3=True,
        lower_active_field_join=False,selected_native_full_C5_pulse_bridge=False,
        C3_numeric_norm_admission=False,repaired_own_terminal_C3=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
