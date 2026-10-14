"""Inspect or restore actual compact/own higher source and 2Rc joins."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_compact_higher_radial_pressure_C3_functions as current

CHECKS={'actual_profile_calculus','actual_history_and_Q_calculus','actual_compact_pressure_calculus',
    'actual_own_velocity_calculus','actual_own_recovery_calculus','actual_compact_pressure_FTC_source',
    'actual_total_pressure_source','actual_outside_own_radial_source','actual_2Rc_source_joins'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert set(rec['checks'])==CHECKS and all(proof['passed'] for proof in rec['checks'].values())
    assert rec['leading_2Rc_higher_pressure_C3_rows']==24
    assert rec['own_2Rc_y4_C3_QV_C2_rows']==150 and rec['total_own_2Rc_higher_pressure_C3_rows']==16
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalCompactHigherRadialPressureC3Functions()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert len(owner.profiles)==5 and len(owner.histories)==5 and set(owner.own_outside)=={'quiet','main','gap','end'}
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        compact_modified_profiles_and_histories_y4_C5=True,compact_Q_y4_C4_and_pressure_y3_C3=True,
        own_reserved_and_outside_y4_C3_QV_C2=True,leading_2Rc_higher_pressure_C3_rows=24,
        own_2Rc_y4_C3_QV_C2_rows=150,total_own_2Rc_higher_pressure_C3_rows=16,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        combined_outer_seam_accepted=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
