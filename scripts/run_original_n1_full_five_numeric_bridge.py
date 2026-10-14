"""Inspect the saved checked true complete signed five C1 moment enclosures."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge as current


def main():
    raw=current.true.parent.linear.read(current.NAME)
    receipt=current.true.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE) or any(receipt.get(k) for k in current.OPEN)
        or pins!=raw['input_hashes'] or archive!=current.sha(current.NAME)
        or receipt.get('source_identity')!=raw['source_identity']
        or receipt.get('source_family')!=raw['source_family']):
        raise ValueError('Checked actual complete true five C1 enclosure archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed complete true moment source '+name)
    queries=[dict(name=row['name'],Z=[row['Z']['lower'],row['Z']['upper']],
        velocity_C0_C1_rows=row['actual_true_velocity_value_Z1_rows'],
        pressure_C0_C1_rows=row['actual_true_pressure_value_Z1_rows'],
        nonzero_source_rows_retained=row['nonzero_sector_rows_retained']) for row in receipt['query_results']]
    print(json.dumps(dict(checked=True,original_queries=queries,
        complete_signed_moments=['m11','m12','m13','m14','m15'],
        pressure_support='inner + Rin-rminus + rminus-Rp + Rp-Rv; zero N1p beyond Rv',
        pressure_units='reference*A*Rv^p; p=1 for m13, p=2 for m15; fixed in Z',
        velocity_units='actual signed physical log sectors; full support at Rin',
        positive_order_five_moment_repair_closed=False,actual_coefficient_norm_bounds_transferred=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
