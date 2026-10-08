"""Focused firstbridge wiring and new full-route derivative replay.

The centered mathematical backend is already independently checked. This
check verifies its actual firstbridge source binding, all twelve transport
steps, unchanged C0/P0 and new five targets, without querying any live owner.
"""
import json,time
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_firstbridge_centered_C1 as source
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,encode=source.HERE,source.sha,source.encode


@source.slope.rc.native.inlet.source_precision
def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    assert all(manifest[k] is False for k in source.common.current.FLAGS)
    wm=json.loads((HERE/source.weighted.NAME).read_bytes());cm=json.loads((HERE/source.slope.source.NAME).read_bytes())
    sigma_manifest=json.loads((HERE/(source.PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
    c=MPIntervalContext();c.dps=240
    sigma={int(k):source.centered.previous.read_cap(c,v) for k,v in sigma_manifest['sigma_global_derivative_log_caps'].items()}
    tiles=0;primitive_reductions=0;density_reductions=0;target_reductions=0
    for archive in manifest['actual_firstbridge_centered_C1_archives']:
        data=source.terminal.load_archive(archive);payload=source.terminal.load_archive(data['accepted_original_pre_slope_archive'])
        sd=source.terminal.load_archive(data['accepted_original_slope_archive']);Z=data['exact_Z_range']
        assert data['source_family']==payload['source_family']==sd['source_family']==manifest['source_family']
        assert data['candidate_N']==source.N and Z==payload['exact_Z_range']==sd['exact_Z_range']
        velocity=data['firstbridge_centered_source']['original_live_velocity_source']
        assert velocity['all_original_root_rows_exactly_match_accepted_archive']
        got=source.tightened_source(c,payload,velocity,sigma)
        exact_replay_equal(encode(got['record']),data['firstbridge_centered_source'],'actual-centered-firstbridge-source')
        assert encode(got['selected_Z_contributions'])==data['selected_firstbridge_Z_contributions']
        primitive_reductions+=sum(row['strict_upper_reduction'] for row in got['record']['primitive_bound_comparisons'].values())
        density_reductions+=sum(row['strict_upper_reduction'] for row in got['record']['firstbridge_contribution_comparisons'].values())
        pressure=next(row for row in wm['actual_weighted_pressure_tiles'] if source.weighted.same_Z(row['exact_Z_range'],Z))
        upstream=source.replay_upstream(c,payload,got,pressure)
        exact_replay_equal(encode(upstream['trace']),data['actual_updated_pre_slope_C1_attribution'],'new-upstream-attribution')
        exact_replay_equal(encode(upstream['record']),data['actual_twelve_chart_transport'],'all-twelve-chart-Z-transport')
        assert upstream['jets']==data['actual_updated_pre_slope_Z']
        primary=next(row for row in wm['genuine_original_O2_weighted_pressure_replays'] if row['ordered_source_cells']==2048 and source.weighted.same_Z(row['exact_Z_range'],Z))
        report=next(row for row in cm['original_complete_same_N_source_incoming_O2_integral_refinements'] if row['ordered_source_cells']==2048 and source.weighted.same_Z(row['exact_Z_range'],Z))
        terminal=source.updated_terminal(c,sd,upstream,primary,report)
        exact_replay_equal(encode(terminal),data['actual_updated_Rc_terminal_defects'],'new-five-target-C1-composition')
        previous=sd['actual_updated_Rc_terminal_defects_and_source_ledger']
        assert encode(terminal['actual_signed_five_terminal_defect_C0'])==previous['actual_signed_five_terminal_defect_C0']
        for key,row in terminal['actual_five_target_diagnostics'].items():
            before=source.ep(source.iv(c,previous['actual_five_target_diagnostics'][key]['Z_log_absolute_upper']))[1]
            after=source.ep(source.iv(c,encode(row['Z_log_absolute_upper'])))[1]
            assert after<=before;target_reductions+=after<before
        assert terminal['pressure_rate_zero_memory_and_separate_original_P0_retained']
        assert terminal['actual_other_four_Z_rows_recomposed_from_new_upstream_and_unchanged_sources']
        assert not terminal['actual_nonlinear_controls_or_terminal_function_identity_admitted']
        assert all(data[k] is False for k in source.common.current.FLAGS);tiles+=1
        print('Firstbridge centered native jets, twelve chart transport and new Rc C1 targets:',Z,'PASS',flush=True)
    assert primitive_reductions==4 and density_reductions>0 and target_reductions>0
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        actual_two_tile_firstbridge_source_and_twelve_chart_and_terminal_replays=tiles,strict_native_primitive_Z_cap_reductions=primitive_reductions,
        strict_firstbridge_density_Z_integral_cap_reductions=density_reductions,strict_five_terminal_Z_cap_reductions=target_reductions,
        original_all_terminal_C0_covers_exactly_unchanged=True,original_P0_and_pressure_zero_rate_retained=True,
        independent_centered_backend_math_reused_from_hash_bound_accepted_receipt=True,
        live_owners_or_other11_upstream_charts_slope_producers_or_inverse_solvers_not_requeried=True,
        actual_full_function_controls_global_N_heat_stress_recursion_admitted=False,**dict.fromkeys(source.common.current.FLAGS,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
