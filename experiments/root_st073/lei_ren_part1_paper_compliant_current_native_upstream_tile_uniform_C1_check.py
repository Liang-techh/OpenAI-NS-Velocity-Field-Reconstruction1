"""Check refined native tile covers and replay accepted O2 affine integrals.

This reads saved complete native source routes. It does not reconstruct the
native owner, rerun original O2 integrals, or rerun accepted scalar fixtures.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_O2_tile_upstream_replay as replay
import lei_ren_part1_paper_compliant_current_original_O2_source_incoming_common_N_check as accepted

source=replay.upstream;common=source.common;ep=common.ep;iv=common.interval
KEYS=source.KEYS;HERE=source.HERE;sha=source.sha;same=accepted.same_factored


def primitive_rows(proof,uniform=False):
    key='original_whole_period_C0_Z_covers' if uniform else 'original_whole_period_primitive_C0_Z_covers'
    if key in proof:return proof[key]
    assert proof['active_support_empty_by_original_kappa_ge_a']
    return None


def scalar_log(value):
    if value is None:return None
    return mp.make_mpf(tuple(value['exact_mpf_tuple']))


def primitive_check(trace):
    assert trace['same_original_source_root_objects']
    assert trace['all_original_cutoff_branches_and_q_Z_retained']
    assert trace['no_derivative_of_interval_selector_or_C0_cap']
    assert trace['original_periodic_parameter_C1_theorem']==source.serial.periodic_parameter_theorem()
    old=primitive_rows(trace['old_whole_period_primitive_cover'])
    new=primitive_rows(trace['parameter_uniform_primitive_cover'],True)
    reductions=0
    for key,decision in trace['same_source_bound_comparisons'].items():
        assert decision['complete_certified_source_function_cover_selected_not_field_value']
        a=scalar_log(decision['old_absolute_upper_log']);b=scalar_log(decision['uniform_absolute_upper_log'])
        tighter=b is None or (a is not None and b<a)
        assert decision['selected_backend']==('parameter_uniform' if tighter else 'inherited_whole_period')
        strict=tighter and a!=b
        assert decision['strict_absolute_upper_reduction']==strict
        selected=trace['selected_primitive_C0_Z_covers'][key]
        candidates=new if tighter else old
        if candidates is None:assert selected['exact_zero']
        else:assert selected==candidates[key]
        assert selected['encloses_original_source_function'] and not selected['point_value_selected']
        reductions+=strict
    return reductions


def route_rows(route):
    switch=route['known_actual_inlet_to_R110_C1_history']
    second=switch['known_actual_inlet_to_phase2_C1_history']
    first=second['known_original_first_bridge_phase1_history']
    rows=[('active_first_bridge',first),('second_bridge',second)]
    rows.extend(switch['actual_serial_chart_C1_history_records'].items())
    rows.extend(route['actual_serial_chart_C1_history_records'].items())
    assert tuple(k for k,r in rows)==source.LABELS
    return first['known_original_initial_collar_history'],rows


def row_fields(label):
    if label=='active_first_bridge':
        return ('actual_entire_remaining_first_bridge_geometry','original_whole_cell_source_and_cutoff_C1',
            'actual_first_bridge_exit_correction_','original_phase1_background_and_separate_P0_Z',
            'actual_first_bridge_exit_own_history_',None)
    if label=='second_bridge':
        return ('actual_whole_second_bridge_geometry','original_whole_second_bridge_source_C1',
            'actual_phase2_correction_','original_phase2_background_and_separate_P0_Z',
            'actual_phase2_own_history_','actual_inherited_phase1_correction_')
    return ('actual_whole_native_chart_geometry','original_full_box_signed_C1_source',
        'actual_right_correction_','actual_right_original_background_and_separate_P0_Z',
        'actual_right_own_history_','actual_inherited_correction_')


def own_checks(c,coordinates,background,values,jets,own,ownZ):
    assert background['original_units']==common.previous.five.UNITS
    assert background['P0_not_merged_into_pressure_history'] and background['ordinary_Z_factorial_conversion']
    assert background['original_histories_not_zeroed_or_reset']
    for kind,correction,saved in (('C0',values,own),('Z',jets,ownZ)):
        for key in KEYS:
            original=common.restore_common_source(background['original_normalized_history_'+kind+'_enclosures'][key],coordinates)
            same(c,saved[key],original+correction[key])


def tile_checks(c,payload,archive):
    coordinates,final,finalZ,P0,P0_Z=replay.inputs(c,payload)
    route=payload['complete_actual_pre_O2_tile_route'];flat,rows=route_rows(route)
    family=payload['source_family'];Z=ep(c.mpf(payload['exact_Z_range']))
    assert ep(iv(c,route['Z_box']))==Z
    assert payload['exact_Z_range']==archive['exact_Z_range']
    assert payload['full_source_bounds_requeried_on_this_tile_not_relabelled_whole_Z']
    assert payload['original_source_normalization_densities_masses_and_incoming_equations_unchanged']
    assert flat['exact_collar_fraction_endpoints']==['1/2','3/4']
    assert flat['original_correction_initial_condition_is_the_checked_sc_half_inlet']
    assert flat['original_background_histories_and_P0_not_reset']
    assert ep(iv(c,flat['original_eta_log_cover']))[1]<=ep(c.ln(iv(c,flat['original_kappa_minus2_lower']))-c.ln(2))[0]
    values={k:common.restore_common_source(v,coordinates) for k,v in flat['actual_correction_C0_enclosures'].items()}
    jets={k:common.restore_common_source(v,coordinates) for k,v in flat['actual_correction_Z_enclosures'].items()}
    assert all(v.zero for v in (*values.values(),*jets.values()))
    reductions=0;attribution=[]
    for (label,row),trace in zip(rows,payload['twelve_same_source_uniform_primitive_comparisons'],strict=True):
        assert row['source_family']==family and row['candidate_N']==source.N and ep(iv(c,row['Z_box']))==Z
        geometry_key,source_key,out,background_key,own,inherited=row_fields(label)
        assert row['original_whole_period_C1_cover']==trace
        reductions+=primitive_check(trace)
        provenance=row[source_key]['source_provenance']
        assert provenance['source_family']==family and ep(iv(c,provenance['Z_box']))==Z
        assert provenance['one_original_seed_graph_asserted'] and not provenance['cache_cover']
        assert provenance['selector_derivatives_not_substituted'] and provenance['numerical_caps_not_used_as_defining_function_values']
        g=row[geometry_key]
        assert ep(iv(c,provenance['coordinate_box']))==ep(iv(c,g['native_coordinate_box']))
        geometry=dict(width=common.restore_common_source(g['positive_true_log_radius_width'],coordinates),
            regular=iv(c,g['regular_true_log_radius_width_cover']),
            scalar_cover=iv(c,g['scalar_width_cover_used_only_for_directed_kernel_bounds']))
        assert not geometry['width'].zero and ep(geometry['width'].coefficient)[0]>0
        adds={k:common.restore_common_source(v,coordinates) for k,v in row['true_log_radius_signed_C0_contributions'].items()}
        addsZ={k:common.restore_common_source(v,coordinates) for k,v in row['true_log_radius_signed_Z_contributions'].items()}
        recorded=payload['per_chart_actual_true_width_contribution_attribution'][label]
        assert recorded==dict(C0=row['true_log_radius_signed_C0_contributions'],Z=row['true_log_radius_signed_Z_contributions'])
        output={};outputZ={}
        for key,rate in source.middle.RATES.items():
            if inherited:
                same(c,row[inherited+'C0'][key],values[key]);same(c,row[inherited+'Z'][key],jets[key])
            factors=source.middle.transfer.true_width_kernel(coordinates,geometry,rate)
            assert not factors['mass'].zero and ep(factors['mass'].coefficient)[0]>0
            if rate==0:
                assert factors['decay'].scale.powers==(0,0,0,0,0) and ep(factors['decay'].coefficient)==(1,1)
            output[key]=factors['decay']*values[key]+adds[key]
            outputZ[key]=factors['decay']*jets[key]+addsZ[key]
            same(c,row[out+'C0'][key],output[key]);same(c,row[out+'Z'][key],outputZ[key])
        values,jets=output,outputZ
        own_checks(c,coordinates,row[background_key],values,jets,row[own+'C0'],row[own+'Z'])
        attribution.append(dict(label=label,C0_log_absolute_uppers={k:common.absolute_log_upper(v) for k,v in adds.items()},
            Z_log_absolute_uppers={k:common.absolute_log_upper(v) for k,v in addsZ.items()}))
    for key in KEYS:
        same(c,route['actual_original_inlet_to_O2_inlet_correction_C0'][key],values[key])
        same(c,route['actual_original_inlet_to_O2_inlet_correction_Z'][key],jets[key])
    own_checks(c,coordinates,route['actual_O2_inlet_original_background_and_separate_P0_Z'],final,finalZ,
        route['actual_O2_inlet_own_history_C0'],route['actual_O2_inlet_own_history_Z'])
    assert reductions==payload['same_source_primitive_strict_upper_reductions']==archive['same_source_primitive_strict_upper_reductions']==34
    return dict(exact_Z_range=payload['exact_Z_range'],source_cells=13,active_charts=12,
        same_source_primitive_strict_upper_reductions=reductions,
        full_factored_serial_C0_Z_affine_rows=120,original_background_added_once=True,
        original_initial_zero_is_source_defined_not_downstream_reset=True,
        positive_true_width_masses_and_rate_zero_pressure_memory_preserved=True,
        per_chart_actual_contribution_log_attribution=attribution)


def replay_checks(c,payload,record,archive):
    result=replay.replay(c,payload,record,archive)
    coordinates,values,jets,P0,P0_Z=replay.inputs(c,payload)
    operator=source.middle.history.C1DuhamelOperator(coordinates)
    original={k:coordinates.scalar(iv(c,v)) for k,v in record['source_defined_original_histories_at_y1'].items()}
    originalZ={k:coordinates.scalar(iv(c,v)) for k,v in record['source_defined_original_history_Z_at_y1'].items()}
    operator.append(c.mpf(1),{k:coordinates.scalar(iv(c,v)) for k,v in record[accepted.TOTALS[0]].items()},
        {k:coordinates.scalar(iv(c,v)) for k,v in record[accepted.TOTALS[1]].items()},payload['source_family'])
    propagated=operator.apply(values,jets,payload['source_family']);comparisons={}
    before=record['actual_source_incoming_and_own_history_transport']
    for kind,field,output,background in (('C0','values',propagated['values'],original),('Z','Z_derivatives',propagated['Z_derivatives'],originalZ)):
        for key in KEYS:
            same(c,result['propagated_actual_correction_'+kind][key],output[key])
            same(c,result['actual_original_plus_propagated_correction_'+kind][key],background[key]+output[key])
            old=common.restore_common_source(before['propagated_actual_correction_'+kind][key],coordinates)
            a=ep(common.absolute_log_upper(old))[1];b=ep(common.absolute_log_upper(output[key]))[1]
            row=dict(before_absolute_upper_log=a,after_absolute_upper_log=b,strict_upper_reduction=b<a,
                logarithmic_envelope_comparison_not_physical_residual=True)
            # Ratios are meaningful only when subtraction of large logs is
            # resolved and safely materializable. Giant derivative logs stay formal.
            if max(abs(a),abs(b))<mp.mpf('1e100') and abs(b-a)<100:
                row['magnitude_upper_ratio']=mp.exp(b-a)
            comparisons[kind+'_'+key]=row
    same(c,result['actual_absolute_pressure'],P0+(original['p']+propagated['values']['p']))
    same(c,result['actual_absolute_pressure_Z'],P0_Z+(originalZ['p']+propagated['Z_derivatives']['p']))
    for bad in (dict(record,explicit_candidate_N=7),dict(record,source_family={}),
            dict(record,exact_Z_range=['.37','.38'])):
        try:replay.replay(c,payload,bad,archive)
        except ValueError:pass
        else:raise AssertionError('Changed original source/N/tile was admitted')
    result['downstream_actual_correction_bound_comparisons']=comparisons
    result['independent_accepted_affine_C1_operator_matches']=True
    return result


def run():
    begin=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    old=json.loads((HERE/common.NAME).read_bytes());receipt=json.loads((HERE/common.RECEIPT).read_bytes())
    assert receipt['all_passed'] and receipt[common.GATE] and old['source_family']==manifest['source_family']
    assert manifest['exact_upstream_source_function_incoming_binding']==old['exact_upstream_source_function_incoming_binding']
    hashes={**manifest['input_hashes'],source.NAME:sha(source.NAME)};checks=[];replays=[]
    c=MPIntervalContext();c.dps=240
    with mp.workdps(300):
        for archive in manifest['actual_strict_sign_upstream_tile_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw)
            assert payload['source_family']==manifest['source_family'] and payload['candidate_N']==source.N
            assert payload['exact_upstream_source_function_incoming_binding']==manifest['exact_upstream_source_function_incoming_binding']
            checks.append(tile_checks(c,payload,archive));hashes[archive['filename']]=sha(archive['filename'])
            for record in old['original_complete_same_N_source_incoming_O2_integral_refinements']:
                if tuple(replay.Fraction(v) for v in record['exact_Z_range'])==tuple(replay.Fraction(v) for v in payload['exact_Z_range']):
                    replays.append(replay_checks(c,payload,record,archive))
            print('Actual refined upstream tile and unchanged N1024 O2 replay PASS',payload['exact_Z_range'],flush=True)
    assert len(checks)==2 and len(replays)==4
    hashes[Path(replay.__file__).name]=sha(Path(replay.__file__).name)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    report=dict(**{replay.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        genuine_refined_upstream_original_O2_replays=replays,accepted_O2_integrals_reused_without_recomputation=True,
        refined_upstream_manifest=source.NAME,refined_upstream_manifest_sha256=sha(source.NAME),
        all_original_P0_and_pressure_correction_memory_preserved=True,
        global_terminal_closure_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False),
        input_hashes=hashes,execution_seconds=time.monotonic()-begin)
    encode=lambda v:source.middle.packets.encode(common.base.encoded(v))
    (HERE/replay.NAME).write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf8')
    result=dict(all_passed=True,**{source.GATE:True,replay.GATE:True},source_family=manifest['source_family'],
        candidate_N=source.N,actual_strict_sign_upstream_tile_checks=checks,accepted_original_O2_replays=4,
        accepted_parameter_uniform_scalar_theorem_reused_not_rerun=True,
        input_hashes={**hashes,replay.NAME:sha(replay.NAME)},execution_seconds=time.monotonic()-begin,
        global_terminal_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.RECEIPT).write_text(json.dumps(encode(result),indent=2)+'\n',encoding='utf8')
    print('Genuine upstream tile refinement and actual O2 replay PASS',flush=True)
    return result


if __name__=='__main__':run()
