"""Inspect or restore the actual current pulse radial/pressure C3 source."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true')
    args=parser.parse_args()
    rec=current.parent.linear.read(current.RECEIPT);raw=current.parent.linear.read(current.NAME)
    assert rec['all_passed'] and all(rec[k] for k in current.GATES) and not any(rec[k] for k in current.OPEN)
    assert rec['actual_second_radial_source_rows']==54 and rec['actual_partial_integrand_FTC_rows']==120
    assert set(rec['checks'])=={'second_radial_source','partial_primitive_FTC','generic_recovery_FTC'}
    assert all(v['passed'] for v in rec['checks'].values())
    pins=dict(rec['input_hashes']);pins.pop(current.NAME)
    pins.pop(Path(current.__file__).stem+'_check.py')
    assert pins==raw['input_hashes'] and rec['source_family']==raw['source_family']
    for name,digest in rec['input_hashes'].items():
        assert current.sha(name)==digest,('Changed source',name)
    restored=False
    if args.restore_source:
        owner=current.CurrentOriginalPulseRadialPressureC3Functions()
        assert owner.acceptance_loaded and not owner.numerical_oracle_available
        assert owner.at_fraction('main','1/2')['exact_fraction']=='1/2'
        restored=True
    proof=rec['checks']['generic_recovery_FTC']
    print(json.dumps(dict(checked=True,source_only_restored=restored,
        actual_second_radial_source_rows=54,actual_partial_integrand_FTC_rows=120,
        actual_base_pressure_FTC_whole_function_identities=len(proof['actual_signed_pressure_FTC_rows']),
        pressure_FTC_higher_rows_by_ordinary_derivative_commutation=len(proof['pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation']),
        exact_graph_nodes=rec['exact_graph_nodes'],numerical_oracle_available=False,
        outer_total_radial_pressure_joins=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
