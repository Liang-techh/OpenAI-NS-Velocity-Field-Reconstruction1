"""Inspect the checked quiet Rp source without replaying accepted owners."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=rec['source_family'] or not all(c['passed'] for c in rec['checks'].values())):
        raise ValueError('Checked actual quiet Rp C5 and pressure C3 required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed quiet source '+name)
    print(json.dumps(dict(checked=True,actual_2Rc_to_Rp_quiet_C5=True,
        quiet_profile_history_y0_to_y4_C5=True,quiet_direct_Qyy_C4_N1p_C3=True,
        actual_K_equals_original_Rp_Pstar_kappa=True,
        general_pressure_equals_original_closed_quiet_primitive_C3=True,
        complete_rminus_to_Rp_pressure_component_C3=True,
        partitioned_Rin_to_Rp_pressure_sum_C3=True,
        inner_0_Rin_pressure_C3=False,lower_active_rminus_full_field_seam=False,
        complete_0_to_Rp_pressure_C3=False,selected_native_pulse_full_C5_bridge=False,
        original_scale_high_order_numeric_containment=False,global_frequency=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
