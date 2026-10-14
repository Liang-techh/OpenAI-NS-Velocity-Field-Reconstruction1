"""Inspect checked true reserved own histories without rebuilding ancestors."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields as current


def main():
    raw=current.parent.linear.read(current.NAME)
    receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker=pins.pop(Path(current.__file__).stem+'_check.py',None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(Path(current.__file__).stem+'_check.py')
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked source-bound own-field archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed reserved own-field source '+name)
    print(json.dumps(dict(checked=True,source_bound_own_fields=list(raw['own_field_bindings']['fields']),
        C1_mean_M_K_P_and_five_cumulative_histories=True,Q_V_genuine_value_only=True,
        chart='Rc <= R <= 2Rc; exact s/log(2) in [0,1]',
        fixed_leading_pressure_source_and_actual_core_datum_retained=True,
        patch_pressure_gradient_is_not_total_pressure_gradient=True,
        checks={key:value.get('passed') for key,value in receipt['checks'].items()},
        terminal_function_zeros=False,required_C3_C5_bounds=False,
        positive_order_five_moment_repair_closed=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
