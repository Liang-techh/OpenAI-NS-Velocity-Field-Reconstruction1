"""Inspect source-pinned actual two axial/three angular C0/C1 coefficient bounds."""
import json
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_reserved_coefficient_enclosures as current


def main():
    raw=current.full.true.parent.linear.read(current.NAME)
    receipt=current.full.true.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    if (not receipt.get('all_passed') or not receipt.get(current.GATE) or any(receipt.get(k) for k in current.OPEN)
        or pins!=raw['input_hashes'] or archive!=current.sha(current.NAME)
        or receipt.get('source_family')!=raw['source_family']
        or receipt.get('source_identity')!=raw['source_identity']):
        raise ValueError('Checked actual reserved coefficient C1 numerical archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed actual coefficient source '+name)
    def short(value,rounding):
        with localcontext() as context:
            context.prec=10;context.rounding=rounding
            return str(+Decimal(value))
    ranges={name:{key:[[short(v['lower'],ROUND_FLOOR),short(v['upper'],ROUND_CEILING)]
        for v in row['normalized_C0_C1_enclosures']] for key,row in sample['coefficients'].items()}
        for name,sample in raw['samples'].items()}
    print(json.dumps(dict(checked=True,actual_coefficients=list(current.COEFFICIENTS),
        queries=receipt['query_results'],normalized_physical_coefficient_C0_C1_ranges=ranges,
        units=dict(alpha=receipt['source_identity']['alpha_fixed_unit'],beta=receipt['source_identity']['beta_fixed_unit']),
        physical_absolute_log_bounds_retained=True,actual_nonzero_mu_and_Astar_Z_retained=True,
        repaired_own_terminal_zeros=False,higher_derivative_C3_C5_bounds=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
