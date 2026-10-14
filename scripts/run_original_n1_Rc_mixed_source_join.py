"""Inspect or independently restore the actual own Rc mixed source join."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_Rc_mixed_source_join as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_n1_Rc_mixed_source_join_check as checker
        owner,_=current.run();checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_Rc_contracted_mixed_join_rows']==246
    assert len(rec['checks'])==6 and all(proof['passed'] for proof in rec['checks'].values())
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalN1RcMixedSourceJoin()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        owner.at_left_fraction('1/2');owner.combined.at_fraction('reserved','1/2');owner.assert_graph()
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_O3_power_flat_source_and_full_nonzero_histories=True,
        actual_O3_power_y4_C5_Q_y4_C4_pressure_y3_C3=True,
        actual_Rc_compact_empty_inlet_and_nonzero_pressure_identity=True,
        actual_Rc_own_mixed_source_join=True,actual_Rc_contracted_mixed_join_rows=246,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        own_Rin_rminus_joined=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
