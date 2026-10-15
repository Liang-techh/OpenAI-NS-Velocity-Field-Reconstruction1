"""Inspect or restore the actual micro2 to macro0 mixed source interface."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join_check as checker
        owner,_=current.run();checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert len(rec['checks'])==7 and all(check['passed'] for check in rec['checks'].values())
    assert rec['actual_second_micro_frozen_macro_contracted_mixed_join_rows']==246
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed micro/macro source',name)
    restored=False;queries=0
    if args.restore_source:
        owner=current.CurrentOriginalMicroMacroMixedSourceJoin()
        assert owner.acceptance_loaded and not owner.numerical_oracle_available
        g=owner.g
        assert owner.jacobian.node==owner.bridge['length']!=owner.base.jacobian.node
        assert owner.loop_input['charts']['second_micro']['domain']==[g.one.node,g.constant(2).node]
        for key in ('S0','S1'):
            for node in owner.modified_input['charts']['second_micro']['local_integrals'][key]:
                integral=g.nodes[node]
                assert integral['lower']==g.one.node and integral['upper']==owner.base.coordinate.node
        for ids in owner.modified_input['endpoints']['second_micro']['integrals'].values():
            for node in ids:
                endpoint=g.nodes[node]
                assert endpoint['operation']=='function_substitution'
                assert endpoint['variable']==owner.base.coordinate.node and endpoint['value']==g.constant(2).node
        assert owner.modified['frozen_macro']['incoming_integrals']==owner.modified_input['endpoints']['second_micro']['integrals']
        for name in current.CHARTS:
            chart=owner.charts[name];lo,hi=chart['domain'];coordinate=owner.ref(owner.loop[name]['coordinate'])
            for fraction in ('0','1/2','1'):
                value=g.add(lo,g.mul(g.constant(Fraction(fraction)),g.sub(hi,lo)))
                for jets in list(chart['physical'].values())+[chart['N_y'],chart['pressure_y'],chart['m15_y']]:
                    for jet in jets:current.outer.Algebra(g,len(current.outer.rows(jet))-1).at(jet,coordinate,value)
                owner.at(name,chart['own'],value);queries+=1
        owner.assert_graph();restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,actual_two_sided_source_queries=queries,
        actual_macro_U_mean_y4_C5_Q_y4_C4=True,actual_macro_pressure_N_y2_C3=True,
        separate_hb_length_and_periodic_phase_lifts=True,
        actual_second_micro_frozen_macro_own_mixed_rows=246,exact_graph_nodes=rec['exact_graph_nodes'],
        numerical_oracle_available=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
