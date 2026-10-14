"""Inspect, restore or rebuild the finite actual own combined outer source."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_combined_outer_mixed_source_functions as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_n1_combined_outer_mixed_source_functions_check as checker
        owner,_=current.run()
        checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['contracted_rows_per_seam']==246 and rec['actual_contracted_mixed_join_rows']==1230
    assert len(rec['checks'])==6 and all(proof['passed'] for proof in rec['checks'].values())
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalN1CombinedOuterMixedSourceFunctions()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        for name in ('reserved','quiet','main','gap','end'):owner.at_fraction(name,'1/2')
        owner.at_tail_log_offset('1/2')
        owner.assert_graph();restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        finite_combined_outer_source_accepted=True,contracted_seams=list(current.SEAMS),
        state_physical_y_order=4,state_ordinary_Z_order=3,QV_ordinary_Z_order=2,
        pressure_and_cumulative_physical_y_order=3,pressure_and_cumulative_ordinary_Z_order=3,
        contracted_rows_per_seam=246,actual_contracted_mixed_join_rows=1230,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        all_order_smoothness_quantitative_bounds=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
