"""Inspect or restore actual active micro radial4 and own mixed seam sources."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_micro_active_radial4_pressure_C3_functions as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_micro_active_radial4_pressure_C3_functions_check as checker
        owner,_=current.run();checker.run(owner)
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert len(rec['checks'])==6 and all(check['passed'] for check in rec['checks'].values())
    assert rec['actual_first_second_micro_contracted_mixed_join_rows']==246
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed active micro source',name)
    restored=False;queries=0
    if args.restore_source:
        owner=current.CurrentOriginalMicroActiveRadial4PressureC3Functions()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert owner.background['fields']==owner.base.bridge['micro']['fields']
        intermediate=current.previous.common_source.previous
        assert current.parent.linear.read(intermediate.RECEIPT)['checks']['source']['fixed_pulse_pressure_is_original_leading_N1p_forcing']
        g=owner.g
        for name in current.CHARTS:
            chart=owner.charts[name];lo,hi=chart['domain']
            for key in ('E','V'):assert owner.loop[name]['roots'][key]==owner.background['fields'][key][:6]
            for fraction in ('0','1/2','1'):
                x=g.add(lo,g.mul(g.constant(Fraction(fraction)),g.sub(hi,lo)))
                for jets in list(chart['physical'].values())+[chart['N_y'],chart['pressure_y'],chart['m15_y']]:
                    for jet in jets:current.outer.Algebra(g,len(current.outer.rows(jet))-1).at(jet,owner.coordinate,x)
                owner.at(chart['own'],x);queries+=1
        owner.assert_graph();restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,actual_two_sided_active_micro_queries=queries,
        actual_active_U_mean_y4_C5_Q_y4_C4=True,actual_pressure_N_y2_C3=True,
        actual_periodic_lift_not_raw_fractional_part_derivative=True,
        actual_first_second_micro_own_mixed_rows=246,exact_graph_nodes=rec['exact_graph_nodes'],
        numerical_oracle_available=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
