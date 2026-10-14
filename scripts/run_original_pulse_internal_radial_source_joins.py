"""Inspect or restore actual internal pulse radial source joins."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_pulse_internal_radial_source_joins as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[key] for key in current.GATES) and not any(rec[key] for key in current.OPEN)
    assert set(rec['checks'])=={'fixed_kernel_source','selected_linear_source','leading_internal_sources',
        'internal_primitive_composition','own_second_radial_source'}
    assert all(proof['passed'] for proof in rec['checks'].values())
    assert rec['actual_internal_radial_pressure_and_density_rows']==24
    assert rec['actual_internal_pressure_primitive_rows']==16 and rec['actual_own_second_radial_rows']==90
    pins=dict(rec['input_hashes']);pins.pop(current.NAME);pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and raw['source_family']==rec['source_family']
    for name,digest in rec['input_hashes'].items():assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalPulseInternalRadialSourceJoins()
        owner.assert_graph();assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert set(owner.leading_traces)=={'main_gap','gap_end'}
        restored=True
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        internal_radial_pressure_and_density_C3_rows=24,retained_pressure_primitive_C3_rows=16,
        actual_own_second_radial_rows=90,exact_graph_nodes=rec['exact_graph_nodes'],
        numerical_oracle_available=False,pressure_second_radial_source=False,
        combined_outer_seam_accepted=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
