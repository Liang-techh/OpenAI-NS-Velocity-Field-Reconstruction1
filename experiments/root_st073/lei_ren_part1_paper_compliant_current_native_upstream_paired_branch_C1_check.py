"""Focused saved-branch checks and actual original O2 incoming replay."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_upstream_paired_branch_C1 as source
import lei_ren_part1_paper_compliant_current_native_upstream_tile_uniform_C1_check as accepted

common=source.common;HERE=source.HERE;sha=source.sha;ep=source.ep;iv=common.interval
same=accepted.same;KEYS=source.KEYS;replay=accepted.replay


def native_restore(record,bases,ledger):
    c=bases[0].ctx;s=record['formal_positive_scale']
    value=source.prior.ScaledEnclosure(source.prior.FormalScale(bases,
        tuple(s['source_exponents'])+(s['radius_power'],),iv(c,s['additional_log_interval'])),
        iv(c,record['coefficient_interval']),ledger)
    same(c,record,value)
    return value


def firstbranch_check(c,trace):
    accepted.primitive_check(trace['accepted_uniform_baseline_cover'])
    saved=trace['genuine_same_source_conditional_paired_cover'];bases=tuple(iv(c,v) for v in saved['native_source_log_bases'])
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    roots={name:{(int(k[1]),int(k[-1])):native_restore(v,bases,ledger) for k,v in rows.items()}
        for name,rows in saved['original_root_rows'].items()}
    computed=source.branch_paired_cover(roots,iv(c,saved['original_eta_log']),
        iv(c,saved['original_log_a_positive_lower']),iv(c,saved['original_dstar_log']))
    actual=computed['record']['conditional_cutoff_paired_branches']
    assert saved['branches_proved_empty']==computed['record']['branches_proved_empty']
    assert len(actual)==len(saved['conditional_cutoff_paired_branches'])
    comparisons=0
    for original,got in zip(saved['conditional_cutoff_paired_branches'],actual,strict=True):
        assert original['name']==got['name'] and original['condition']==got['condition']
        assert original['only_Delta_C0_restricted_all_original_derivative_objects_retained']
        for name,records in (('original_q_jet',got['original_q_jet']),
                ('q2',got['original_direct_q_squared_jet']['ordinary_q_squared_rows'])):
            before=original[name] if name!='q2' else original['original_direct_q_squared_jet']['ordinary_q_squared_rows']
            for k,r in records.items():same(c,before[k],native_restore(source.encode(r),bases,ledger));comparisons+=1
        for k,r in got['paired_derivative_covers'].items():
            same(c,original['paired_derivative_covers'][k],native_restore(source.encode(r),bases,ledger));comparisons+=1
        if original['name']=='flat':
            assert all(v['exact_zero'] for v in original['original_q_jet'].values())
            assert all(v['exact_zero'] for v in original['paired_derivative_covers'].values())
    selected=trace['selected_primitive_C0_Z_covers'];old=trace['accepted_uniform_baseline_cover']['selected_primitive_C0_Z_covers']
    for k in ('A','B_over_Pstar'):assert selected[k]==old[k]
    for k,value in computed['values'].items():
        same(c,saved['whole_source_paired_derivative_hulls'][k],value)
        before=source.baseline.magnitude_log(native_restore(old[k],bases,ledger));after=source.baseline.magnitude_log(value)
        decision=trace['paired_derivative_comparisons'][k]
        assert accepted.scalar_log(decision['accepted_uniform_absolute_upper_log'])==before
        assert accepted.scalar_log(decision['paired_branch_union_absolute_upper_log'])==after
        tighter=after is None or before is not None and after<before
        assert decision['strict_absolute_upper_reduction']==(tighter and before!=after)
        chosen=value if tighter else native_restore(old[k],bases,ledger)
        same(c,selected[k],chosen)
    assert trace['original_C0_A_B_objects_and_native_density_owner_unchanged']
    return dict(nonempty_cutoff_branches=[row['name'] for row in actual],
        saved_original_branch_q_q_squared_and_paired_jet_rows_recomputed=comparisons,
        paired_parameter_scalar_proof_reused=True,all_conditional_branches_hulled_not_added=True)


def tile_check(c,payload,before):
    coordinates,final,finalZ,P0,P0_Z=replay.inputs(c,payload)
    route=payload['complete_actual_pre_O2_tile_route'];flat,rows=accepted.route_rows(route)
    oldflat,oldrows=accepted.route_rows(before['complete_actual_pre_O2_tile_route'])
    assert payload['source_family']==before['source_family'] and payload['exact_Z_range']==before['exact_Z_range']
    assert payload['candidate_N']==before['candidate_N']==source.N
    assert payload['exact_upstream_source_function_incoming_binding']==before['exact_upstream_source_function_incoming_binding']
    assert flat==oldflat
    values={k:common.restore_common_source(v,coordinates) for k,v in flat['actual_correction_C0_enclosures'].items()}
    jets={k:common.restore_common_source(v,coordinates) for k,v in flat['actual_correction_Z_enclosures'].items()}
    assert all(v.zero for v in (*values.values(),*jets.values()))
    branch=firstbranch_check(c,payload['twelve_source_primitive_cover_records'][0])
    Z=ep(c.mpf(payload['exact_Z_range']))
    for i,((label,row),(oldlabel,oldrow),trace) in enumerate(zip(rows,oldrows,payload['twelve_source_primitive_cover_records'],strict=True)):
        assert oldlabel==label and ep(iv(c,row['Z_box']))==Z
        assert row['source_family']==payload['source_family'] and row['candidate_N']==source.N
        gkey,skey,out,bkey,own,inherited=accepted.row_fields(label)
        assert row[skey]==oldrow[skey] and row[gkey]==oldrow[gkey] and row[bkey]==oldrow[bkey]
        assert row['original_whole_period_C1_cover']==trace
        if i:assert trace==oldrow['original_whole_period_C1_cover']
        assert row['true_log_radius_signed_C0_contributions']==oldrow['true_log_radius_signed_C0_contributions']
        g=row[gkey];geometry=dict(width=common.restore_common_source(g['positive_true_log_radius_width'],coordinates),
            regular=iv(c,g['regular_true_log_radius_width_cover']),scalar_cover=iv(c,g['scalar_width_cover_used_only_for_directed_kernel_bounds']))
        output={};outputZ={}
        for k,rate in source.middle.RATES.items():
            if inherited:
                same(c,row[inherited+'C0'][k],values[k]);same(c,row[inherited+'Z'][k],jets[k])
            decay=source.middle.transfer.true_width_kernel(coordinates,geometry,rate)['decay']
            a=common.restore_common_source(row['true_log_radius_signed_C0_contributions'][k],coordinates)
            z=common.restore_common_source(row['true_log_radius_signed_Z_contributions'][k],coordinates)
            output[k]=decay*values[k]+a;outputZ[k]=decay*jets[k]+z
            same(c,row[out+'C0'][k],output[k]);same(c,row[out+'Z'][k],outputZ[k])
        values,jets=output,outputZ
        accepted.own_checks(c,coordinates,row[bkey],values,jets,row[own+'C0'],row[own+'Z'])
    for k in KEYS:
        same(c,route['actual_original_inlet_to_O2_inlet_correction_C0'][k],values[k])
        same(c,route['actual_original_inlet_to_O2_inlet_correction_Z'][k],jets[k])
        assert route['actual_original_inlet_to_O2_inlet_correction_C0'][k]==before['complete_actual_pre_O2_tile_route']['actual_original_inlet_to_O2_inlet_correction_C0'][k]
    return dict(exact_Z_range=payload['exact_Z_range'],original_source_and_C0_contributions_unchanged_on_all12_active_cells=True,
        complete_original13_cell_C0_Z_affine_transport_checked=True,firstbridge_conditional_paired_check=branch)


def run():
    begin=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    baseline=json.loads((HERE/source.baseline.NAME).read_bytes())
    old_O2=json.loads((HERE/common.NAME).read_bytes())
    old_replay=json.loads((HERE/replay.NAME).read_bytes())
    c=MPIntervalContext();c.dps=240;checks=[];replays=[];hashes=dict(manifest['input_hashes'])
    hashes[source.NAME]=sha(source.NAME)
    with mp.workdps(300):
        for archive in manifest['genuine_firstbridge_paired_branch_upstream_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw);hashes[archive['filename']]=sha(archive['filename'])
            previous=next(a for a in baseline['actual_strict_sign_upstream_tile_archives'] if a['exact_Z_range']==payload['exact_Z_range'])
            before=json.loads(gzip.decompress((HERE/previous['filename']).read_bytes()))
            checks.append(tile_check(c,payload,before))
            for record in old_O2['original_complete_same_N_source_incoming_O2_integral_refinements']:
                if tuple(replay.Fraction(v) for v in record['exact_Z_range'])!=tuple(replay.Fraction(v) for v in payload['exact_Z_range']):continue
                row=accepted.replay_checks(c,payload,record,archive)
                row['actual_upstream_binding']['manifest']=source.NAME
                row['actual_upstream_binding']['manifest_sha256']=sha(source.NAME)
                row['actual_upstream_binding']['firstbridge_conditional_paired_derivative_cover']=True
                old=next(r for r in old_replay['genuine_refined_upstream_original_O2_replays'] if r['exact_Z_range']==row['exact_Z_range'] and r['ordered_source_cells']==row['ordered_source_cells'])
                comparisons={}
                coordinates,*unused=replay.inputs(c,payload)
                for kind in ('C0','Z'):
                    for k in KEYS:
                        field='propagated_actual_correction_'+kind
                        a=common.restore_common_source(old[field][k],coordinates);b=common.restore_common_source(row[field][k],coordinates)
                        la=source.baseline.magnitude_log(a);lb=source.baseline.magnitude_log(b)
                        comparisons[kind+'_'+k]=dict(before_absolute_upper_log=la,after_absolute_upper_log=lb,
                            strict_upper_reduction=lb<la,enclosure_not_physical_amplitude_or_residual=True)
                        if kind=='C0':same(c,old[field][k],b)
                row['bound_comparisons_against_accepted_uniform_tile_replay']=comparisons
                replays.append(row)
            print('Firstbridge branch-correlated source/C1 history/O2 replay PASS',payload['exact_Z_range'],flush=True)
    assert len(checks)==2 and len(replays)==4
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    report=dict(**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        genuine_paired_branch_original_O2_replays=replays,original_C0_contributions_and_backgrounds_unchanged=True,
        accepted_O2_integrals_not_recomputed=True,input_hashes=hashes,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.REPLAY).write_text(json.dumps(source.encode(report),indent=2)+'\n',encoding='utf8')
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        genuine_firstbridge_paired_branch_tile_checks=checks,accepted_original_O2_replays=4,
        input_hashes={**hashes,source.REPLAY:sha(source.REPLAY)},execution_seconds=time.monotonic()-begin,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Genuine firstbridge branch-correlated derivative refinement PASS',flush=True)
    return result


if __name__=='__main__':run()
