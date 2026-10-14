"""Inspect or restore the accepted eight-function current Rp pulse source."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true',
        help='Reconstruct the checked source callback and consumer without numerical ancestor replay.')
    args=parser.parse_args()
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker=Path(current.__file__).stem+'_check.py';checksha=pins.pop(checker,None)
    checks=rec.get('checks',{})
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checksha!=current.sha(checker)
        or raw['source_family']!=rec['source_family']
        or set(checks)!={'source_inputs','outer_common_sources','compact_quiet','callback_and_consumer'}
        or not all(check['passed'] for check in checks.values())
        or len(checks['compact_quiet']['actual_eight_Rp_frame_C5_rows'])!=48
        or len(checks['compact_quiet']['actual_terminal_integral_and_unique_branch_source_proof']['actual_integral_to_residual_C5_rows'])!=30
        or len(checks['callback_and_consumer']['actual_absolute_background_pressure_C3_rows'])!=12):
        raise ValueError('Checked actual current Rp native C5 source required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed current Rp pulse source '+name)
    restored=None
    if args.restore_source:
        owner=current.CurrentOriginalRpNativeC5PulseSource()
        fields=owner.pulse_source_functions()
        assert len(owner.incoming_source_rows())==8 and not fields['numerical_oracle_available']
        restored=dict(checked_callback_and_consumer_restored=True,
            actual_incoming_C5_rows=sum(len(current.previous.rows(jet)) for jet in fields['native'].values()),
            actual_absolute_pressure_C3_rows=sum(len(current.previous.rows(jet)) for jet in fields['absolute_pressure'].values()),
            numerical_oracle_available=False,historical_python_session_required=False)
    print(json.dumps(dict(checked=True,actual_current_Rp_native_functions=list(current.FRAME_KEYS),
        actual_C5_incoming_rows=48,
        actual_selected_coefficient_field_and_pressure_rows=checks['callback_and_consumer']['actual_selected_coefficient_field_and_pressure_source_rows'],
        actual_terminal_integral_projection_rows=30,
        actual_absolute_background_pressure_C3_rows=12,
        exact_graph_nodes=len(raw['exact_graph_nodes']),archive_source_recovery=restored,
        source_function_callback=True,interval_numeric_selected_callback=False,
        original_scale_numeric_bounds=False,other_global_own_joins=False,
        global_n1=False,temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
