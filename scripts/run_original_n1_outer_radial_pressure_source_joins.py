"""Inspect or restore accepted actual Rp/Rv radial-pressure source joins."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_outer_radial_pressure_source_joins as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES)
    assert not any(rec[key] for key in current.OPEN)
    assert rec['actual_Rp_Rv_total_pressure_y_C3_join_rows']==8
    assert rec['actual_Rp_Rv_N1p_C3_join_rows']==8
    assert rec['actual_signed_pressure_C3_source_map_rows']==24
    assert rec['accepted_ordinary_outer_trace_rows_reused']==204
    assert set(rec['checks'])=={'Rp_endpoint_recovery','Rv_endpoint_recovery','actual_chart_radius_sources',
        'partial_provider_source_identity','own_pressure_transport'}
    assert all(proof['passed'] for proof in rec['checks'].values())
    transport=rec['checks']['own_pressure_transport']
    assert len(transport['actual_Rp_Rv_m13_m15_endpoint_rows_by_whole_function_identity_composition'])==16
    assert len(transport['actual_current_own_pressure_radial_rows_by_accepted_FTC_composition'])==24
    pins=dict(rec['input_hashes']);pins.pop(current.NAME)
    pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and rec['source_family']==raw['source_family']
    for name,digest in rec['input_hashes'].items():
        assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalN1OuterRadialPressureSourceJoins()
        owner.assert_graph()
        assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert set(owner.current_own)=={'main','gap','end'}
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_Rp_Rv_total_pressure_y_C3_join_rows=8,actual_Rp_Rv_N1p_C3_join_rows=8,
        actual_partial_and_complete_C0_integral_source_rows=rec['actual_partial_and_complete_C0_integral_source_rows'],
        complete_higher_rows_by_pinned_derivative_commutation=rec['complete_higher_rows_by_pinned_derivative_commutation'],
        actual_signed_pressure_C3_source_map_rows=24,accepted_ordinary_outer_trace_rows_reused=204,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        combined_outer_seam_accepted=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
