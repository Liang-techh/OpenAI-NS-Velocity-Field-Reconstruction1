"""Inspect pinned intermediate own functions without replaying ancestors."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions as current


def main():
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=rec['source_family'] or not all(c['passed'] for c in rec['checks'].values())
        or len(raw['intermediate_own_C3_chart_bindings']['charts'])!=22):
        raise ValueError('Checked actual intermediate own C3 chart functions required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed intermediate own source '+name)
    print(json.dumps(dict(checked=True,actual_regional_charts=22,actual_preRc_charts=17,
        own_mass_mean_K_pressure_C3=True,own_QV_C2=True,
        original_nonzero_incoming_mass_and_pressure_retained=True,
        actual_partial_pulse_integral_bodies_C3=True,
        fixed_pressure_source_is_leading_N1p=True,
        complete_regional_exit_uses_true_terminal_integral_zeros=True,
        lower_active_field_join=False,selected_native_full_C5_pulse_bridge=False,
        C3_numeric_norm_admission=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
