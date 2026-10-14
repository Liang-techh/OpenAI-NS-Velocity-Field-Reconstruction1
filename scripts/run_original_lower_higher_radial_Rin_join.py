"""Inspect or independently restore genuine lower radial4 and full Rin join."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_lower_higher_radial_Rin_join as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_lower_higher_radial_Rin_join_check as checker
        owner,_=current.run();checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_Rin_contracted_mixed_join_rows']==246 and len(rec['checks'])==6
    assert all(group['passed'] for group in rec['checks'].values())
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed higher lower source',name)
    restored=False;queries=0
    if args.restore_source:
        owner=current.CurrentOriginalLowerHigherRadialRinJoin()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        g=owner.g;lo,hi=owner.rin.intermediate['domain']
        for fraction in ('0','1/2','1'):
            x=g.add(owner.ref(lo),g.mul(g.constant(Fraction(fraction)),g.sub(owner.ref(hi),owner.ref(lo))))
            for jets in list(owner.radial.values())+[owner.N,owner.pressure_y,owner.m15_y]:
                for jet in jets:current.outer.Algebra(g,len(current.outer.rows(jet))-1).at(jet,owner.coordinate,x)
            queries+=1
        owner.assert_graph();restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,actual_lower_source_fraction_queries=queries,
        actual_lower_radial4_C5_Q_radial4_C4=True,actual_pressure_N_rho2_C3=True,
        actual_Rin_full_y4_pressure_y3_mixed_source_join=True,actual_Rin_contracted_mixed_join_rows=246,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        own_rminus_joined=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
