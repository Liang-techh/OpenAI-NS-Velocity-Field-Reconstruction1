"""Inspect the checked true complete-pulse C1 numerical checkpoint."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_numeric_bridge as current


def main():
    raw=current.true.parent.linear.read(current.NAME)
    receipt=current.true.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE) or any(receipt.get(k) for k in current.OPEN)
        or pins!=raw['input_hashes'] or archive!=current.sha(current.NAME)
        or receipt.get('source_identity')!=raw['source_identity']
        or receipt.get('source_family')!=raw['source_family']):
        raise ValueError('Checked actual complete pulse C1 enclosure archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed true pulse numerical source '+name)
    print(json.dumps(dict(checked=True,radial_scope='original Rp to Rv, main/gap/end',
        original_queries=receipt['query_results'],nonlinear_C1_double_scaling_corrected=True,
        units='reference*A*Rv^p; p=1 for m13,2 for m15; no physical exponent materialization',
        combined_inner_to_infinity_pressure=False,positive_order_five_moment_repair_closed=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
