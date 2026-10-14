"""Inspect or independently restore the actual cutoff/lower Rin source join."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_n1_Rin_source_join_check as checker
        owner,_=current.run();checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_Rin_contracted_join_rows']==138 and len(rec['checks'])==5
    assert all(proof['passed'] for proof in rec['checks'].values())
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed Rin source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalN1RinSourceJoin()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert set(owner.join)=={'left','right'} and owner.contract==current.CONTRACT
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_Rin_common_leading_pressure_source=True,actual_Rin_signed_FTC_and_empty_lower_inlet=True,
        actual_Rin_state_y2_pressure_moment_y1_source_join=True,actual_Rin_contracted_join_rows=138,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        Rin_full_y4_pressure_y3_mixed_contract=False,own_rminus_joined=False,global_n1=False,
        temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
