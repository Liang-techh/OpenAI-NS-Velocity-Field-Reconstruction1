"""Read the accepted original native-to-selected source identity checkpoint."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_native_numeric_identity as current


def main():
    raw=current.true.parent.linear.read(current.NAME)
    receipt=current.true.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or receipt.get('source_identity')!=raw['source_identity']
        or receipt.get('source_family')!=raw['source_family']):
        raise ValueError('Checked actual native numeric source identity required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed native pulse source '+name)
    print(json.dumps(dict(checked=True,exact_native_scalar_function_identities=4,
        source_Taylor_rows=receipt['source_identity']['production_Taylor_proof'],
        exact_true_and_numeric_factor_recipe_identities=2,
        source_identity_not_inferred_from_interval_overlap=True,
        pulse_pressure_enclosure_transfer=False,positive_order_five_moment_repair_closed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
