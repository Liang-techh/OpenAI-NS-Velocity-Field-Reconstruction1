"""Inspect/restore all ten full original high-source charts and R100 seam."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    parser.add_argument('--rebuild-source',action='store_true')
    args=parser.parse_args()
    if args.rebuild_source:
        import lei_ren_part1_paper_compliant_current_original_prescribed_high_radial_source_functions_check as checker
        owner,_=current.run();checker.run(owner)
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    assert rec['all_passed'] and len(rec['checks'])==7 and all(proof['passed'] for proof in rec['checks'].values())
    assert all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert rec['actual_full_original_charts']==list(current.CHARTS)
    assert rec['actual_native_to_physical_derivative_requests']==180
    assert rec['actual_macro_first_switch_own_mixed_rows']==246
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    restored=False;queries=0
    if args.restore_source:
        owner=current.CurrentOriginalPrescribedHighRadialSourceFunctions()
        assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert tuple(owner.charts)==current.CHARTS
        for name in current.CHARTS:
            assert current.outer.encode(owner.charts[name]['domain'])==owner.loop[name]['domain']
            for fraction in ('0','1/2','1'):
                packet=owner.at_fraction(name,fraction)
                assert not packet['numerical_oracle_available']
                assert len(packet['physical']['U'])==5 and len(packet['physical']['Q'])==5
                assert len(packet['N_y'])==3 and len(packet['own']['fields']['P1'])==4
                queries+=1
        owner.assert_graph();restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_complete_original_charts=list(current.CHARTS),actual_native_to_physical_derivative_requests=180,
        actual_full_chart_fraction_queries=queries,actual_macro_first_switch_own_mixed_rows=246,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        other_intermediate_source_joins_complete=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
