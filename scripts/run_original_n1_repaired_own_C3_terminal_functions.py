"""Inspect the checked repaired own C3 integral functions without replay."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_repaired_own_C3_terminal_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=rec['source_family'] or not all(c['passed'] for c in rec['checks'].values())):
        raise ValueError('Checked actual repaired own C3 terminal integral identities required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed actual repaired own source '+name)
    print(json.dumps(dict(checked=True,actual_signed_patch_histories_C3=True,
        mean_mass_K_pressure_C3=True,own_QV_C2_from_real_mean_Z3=True,
        actual_compact_direct_Qyy_C4=True,total_reserved_pressure_y_C3=True,
        actual_terminal_integral_five_C3_function_zeros=True,
        mass_tail_C3_QV_C2_zero_after_2Rc=True,n1_pressure_tail_C3_zero_after_Rv=True,
        leading_angular_and_heat_pressure_retained=True,
        lower_active_field_join=False,selected_native_full_C5_pulse_bridge=False,
        C3_numeric_norm_admission=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
