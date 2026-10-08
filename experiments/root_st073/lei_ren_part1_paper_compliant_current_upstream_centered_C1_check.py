"""Replay new per-chart centered jets, complete transport and terminal targets."""
import json,time
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_upstream_centered_C1 as source
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal


@source.slope.rc.native.inlet.source_precision
def run():
    began=time.monotonic();m=json.loads((source.HERE/source.NAME).read_bytes())
    assert m[source.GATE] and m['candidate_N']==source.N
    for name,digest in m['input_hashes'].items():assert source.sha(name)==digest,name
    assert all(m[k] is False for k in source.common.current.FLAGS)
    wm=json.loads((source.HERE/source.weighted.NAME).read_bytes());cm=json.loads((source.HERE/source.slope.source.NAME).read_bytes())
    sig=json.loads((source.HERE/(source.PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes());c=MPIntervalContext();c.dps=240
    sigma={int(k):source.centered.previous.read_cap(c,v) for k,v in sig['sigma_global_derivative_log_caps'].items()};tiles=0;queries=0;density=0;primitive=0;targets=0
    for archive in m['actual_all_upstream_centered_C1_archives']:
        data=source.terminal.load_archive(archive);old=source.terminal.load_archive(data['accepted_firstbridge_centered_archive'])
        payload=source.terminal.load_archive(old['accepted_original_pre_slope_archive']);sd=source.terminal.load_archive(old['accepted_original_slope_archive']);Z=data['exact_Z_range']
        assert data['source_family']==old['source_family']==payload['source_family']==m['source_family'] and data['candidate_N']==source.N
        assert Z==old['exact_Z_range']==payload['exact_Z_range']
        charts=data['actual_centered_upstream_chart_sources'];assert list(charts)==list(source.LABELS)
        assert charts['active_first_bridge']['selected_Z_contributions']==old['selected_firstbridge_Z_contributions']
        for index,label in enumerate(source.LABELS[1:],start=1):
            saved=charts[label];live=saved['record']['original_live_source'];row=source.weighted.accepted.route_rows(payload['complete_actual_pre_O2_tile_route'])[1][index][1]
            assert live['label']==label and live['original_source_provenance']==row[source.weighted.accepted.row_fields(label)[1]]['source_provenance']
            got=source.chart_cover(c,payload,live,index,sigma);exact_replay_equal(source.encode(got),saved,'original-centered-upstream-chart')
            primitive+=sum(v['strict_upper_reduction'] for v in got['record']['primitive_comparisons'].values())
            density+=sum(v['strict_upper_reduction'] for v in got['record']['contribution_comparisons'].values());queries+=1
        pressure=next(r for r in wm['actual_weighted_pressure_tiles'] if source.weighted.same_Z(r['exact_Z_range'],Z));upstream=source.transport(c,payload,charts,pressure)
        exact_replay_equal(source.encode(upstream['trace']),data['actual_updated_upstream_attribution'],'new-upstream-attribution')
        exact_replay_equal(source.encode(upstream['record']),data['actual_twelve_chart_transport'],'new-upstream-transport');assert upstream['jets']==data['actual_updated_upstream_Z']
        primary=next(r for r in wm['genuine_original_O2_weighted_pressure_replays'] if r['ordered_source_cells']==2048 and source.weighted.same_Z(r['exact_Z_range'],Z))
        report=next(r for r in cm['original_complete_same_N_source_incoming_O2_integral_refinements'] if r['ordered_source_cells']==2048 and source.weighted.same_Z(r['exact_Z_range'],Z))
        output=source.previous.updated_terminal(c,sd,upstream,primary,report)
        exact_replay_equal(source.encode(output),data['actual_updated_Rc_terminal_defects'],'new-all-upstream-Rc-targets')
        assert source.encode(output['actual_signed_five_terminal_defect_C0'])==old['actual_updated_Rc_terminal_defects']['actual_signed_five_terminal_defect_C0']
        for key,row in output['actual_five_target_diagnostics'].items():
            before=source.ep(source.iv(c,old['actual_updated_Rc_terminal_defects']['actual_five_target_diagnostics'][key]['Z_log_absolute_upper']))[1]
            after=source.ep(source.iv(c,source.encode(row['Z_log_absolute_upper'])))[1];assert after<=before;targets+=after<before
        assert output['pressure_rate_zero_memory_and_separate_original_P0_retained'] and not output['actual_nonlinear_controls_or_terminal_function_identity_admitted']
        assert all(data[k] is False for k in source.common.current.FLAGS);tiles+=1
        print('All original centered upstream charts, twelve steps and five Rc targets:',Z,'PASS',flush=True)
    assert queries==22 and primitive>0 and density>0
    hashes=dict(m['input_hashes']);hashes[source.NAME]=source.sha(source.NAME);hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=m['source_family'],candidate_N=source.N,
        actual_new_native_chart_sources_replayed=queries,actual_complete_upstream_and_terminal_replays=tiles,
        strict_primitive_Z_cap_reductions=primitive,strict_density_Z_integral_cap_reductions=density,strict_terminal_Z_cap_reductions=targets,
        original_C0_targets_exactly_unchanged=True,original_geometry_P0_and_pressure_zero_rate_retained=True,
        accepted_centered_math_firstbridge_slope_and_post_slope_evidence_reused=True,live_owners_and_inverse_solvers_not_rerun_by_checker=True,
        actual_global_functions_controls_N_heat_stress_recursion_admitted=False,**dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (source.HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
