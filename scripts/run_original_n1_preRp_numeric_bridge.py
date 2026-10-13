"""Restore checked source-bound pre-Rp pressure C1 enclosures."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_preRp_numeric_bridge as current


def main():
    raw=current.lower.full.parent.linear.read(current.NAME)
    receipt=current.lower.full.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE) or any(receipt.get(k) for k in current.OPEN)
        or pins!=raw['input_hashes'] or archive!=current.sha(current.NAME)
        or receipt.get('source_family')!=raw['source_family']
        or receipt.get('source_identity')!=raw['source_identity']):
        raise ValueError('Checked true pre-Rp numerical archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed original pre-Rp numerical source '+name)
    print(json.dumps(dict(checked=True,radial_scope='original rminus to Rp',
        actual_original_queries=receipt['query_results'],
        actual_live_field_function_replay_preserved=True,
        signed_defining_functions=['m13','m15'],numeric_rows='source-normalized magnitude enclosures for Z0/Z1',
        complete_external_pressure_and_reserved_coefficients=False,
        positive_order_five_moment_repair_closed=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
