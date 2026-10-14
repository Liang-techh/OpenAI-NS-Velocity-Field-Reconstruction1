"""Inspect the checked actual lower/active source join without replay."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_lower_active_common_source_join_functions as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore-source',action='store_true',
        help='Restore the checked exact graph and seam views; no numerical oracle or ancestor replay.')
    args=parser.parse_args()
    raw=current.parent.linear.read(current.NAME);rec=current.parent.linear.read(current.RECEIPT)
    pins=dict(rec.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not rec.get('all_passed') or not all(rec.get(k) for k in current.GATES)
        or any(rec.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=rec['source_family'] or len(rec['checks'])!=10
        or not all(c['passed'] for c in rec['checks'].values())
        or len(rec['checks']['trace_transfer']['actual_registered_two_sided_trace_rows'])!=44):
        raise ValueError('Checked actual lower/active common-source join required')
    for name,digest in rec['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed lower/active source '+name)
    restored=None
    if args.restore_source:
        import lei_ren_part1_paper_compliant_current_original_lower_active_archive_recovery as recovery
        owner=recovery.restore(True)
        assert owner.at_seam()['actual_C5_Q4_pressure_C3_source_join']
        assert owner.resolve_sc(owner.loop.sc)==owner.lower.lower.parameters['sc']
        restored=dict(exact_graph_nodes=len(raw['exact_graph_nodes']),
            checked_seam_views_restored=True,numerical_oracle_available=False,
            historical_python_session_required=False)
    print(json.dumps(dict(checked=True,actual_lower_active_leading_field_join=True,
        U_Uy_mean_mean_y_swirl_C5=True,Q_Qy_C4=True,absolute_pressure_C3=True,
        exact_original_sc_source_projection_bound=True,actual_six_true_inlet_histories=True,
        archive_source_recovery=restored,
        original_pressure_Pstar_squared_normalization_once=True,
        actual_two_sided_finite_N_flat_collar=True,
        actual_source_trace_rows=44,leading_recovery='1-delta',
        positive_order_recovery='1+delta, separate provider',
        selected_native_full_C5_pulse_bridge=False,quantitative_mixed_join_norms=False,
        other_global_own_joins=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
