"""Inspect or restore actual pulse higher radial pressure source functions."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_pulse_higher_radial_pressure_C3_functions as current

CHECKS={'actual_scalar_and_axial_source','actual_higher_mass_Q_source','actual_pressure_radial_source',
    'actual_quiet_higher_source','actual_internal_higher_pressure_joins','actual_outer_higher_pressure_joins'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert set(rec['checks'])==CHECKS and all(proof['passed'] for proof in rec['checks'].values())
    assert rec['actual_internal_higher_pressure_C3_rows']==48 and rec['actual_Rp_Rv_higher_pressure_C3_rows']==48
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalPulseHigherRadialPressureC3Functions()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert set(owner.charts)=={'main','gap','end'} and set(owner.outer_pressure_traces)=={'Rp','Rv'}
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_Uyyy_Uyyyy_C5_and_Qyyy_Qyyyy_C4=True,actual_pressure_yy_yyy_C3=True,
        internal_higher_pressure_C3_rows=48,Rp_Rv_higher_pressure_C3_rows=48,
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        combined_outer_seam_accepted=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
