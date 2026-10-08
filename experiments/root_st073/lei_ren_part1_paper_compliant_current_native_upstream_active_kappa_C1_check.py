"""Check saved active kappa source refinements and replay unchanged O2."""
from contextlib import contextmanager
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_upstream_active_kappa_C1 as source
import lei_ren_part1_paper_compliant_current_native_upstream_paired_branch_C1_check as history

common=source.common;accepted=history.accepted;replay=history.replay
iv=common.interval;ep=source.ep;same=history.same;HERE=source.HERE;sha=source.sha;KEYS=source.KEYS


def active_firstbranch_check(c,trace):
    saved=trace['genuine_same_source_conditional_paired_cover']
    bases=tuple(iv(c,v) for v in saved['native_source_log_bases'])
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    restore=lambda r:history.native_restore(source.encode(r),bases,ledger)
    roots={name:{(int(k[1]),int(k[-1])):restore(v) for k,v in rows.items()} for name,rows in saved['original_root_rows'].items()}
    got=source.active_branch_cover(roots,iv(c,saved['original_eta_log']),
        iv(c,saved['original_log_a_positive_lower']),iv(c,saved['original_dstar_log']))
    assert saved['original_root_rows']==saved['accepted_paired_cover_on_same_original_roots']['original_root_rows']
    assert saved['branches_proved_empty']==got['record']['branches_proved_empty']
    rows=got['record']['conditional_cutoff_paired_branches'];assert len(rows)==len(saved['conditional_cutoff_paired_branches'])
    reduced=0;counts=0
    for before,actual in zip(saved['conditional_cutoff_paired_branches'],rows,strict=True):
        assert before['name']==actual['name'] and before['condition']==actual['condition']
        proof=before['same_source_active_root_refinement'];computed=actual['same_source_active_root_refinement']
        if before['name']=='flat':
            assert proof==computed==dict(flat_branch_no_active_constraint_or_body_evaluation=True)
            assert all(v['exact_zero'] for v in before['original_q_jet'].values())
        else:
            assert proof['original_active_kappa_theorem']==source.active_kappa_theorem()
            expected_K=2 if before['name']=='negative' else mp.mpf('2.5')
            assert ep(iv(c,proof['conditional_K']))==(expected_K,expected_K)
            for field in ('original_a_C0','conditional_same_source_positive_a_C0','original_a_Z','original_b_Z',
                    'K_over_a','sqrt_K_over_a','inherited_original_Delta_Z','same_original_derivative_conditional_majorant','selected_Delta_Z_cover'):
                same(c,proof[field],restore(computed[field]));counts+=1
            assert proof['original_a_Z']==saved['original_root_rows']['a']['y0_Z1']
            assert proof['original_b_Z']==saved['original_root_rows']['b']['y0_Z1']
            assert proof['strict_Delta_Z_upper_reduction']==computed['strict_Delta_Z_upper_reduction']
            assert proof['root_map_copied_original_roots_not_mutated']
            reduced+=proof['strict_Delta_Z_upper_reduction']
        for field in ('original_q_jet','paired_derivative_covers'):
            for k,value in actual[field].items():same(c,before[field][k],restore(value));counts+=1
        for k,value in actual['original_direct_q_squared_jet']['ordinary_q_squared_rows'].items():
            same(c,before['original_direct_q_squared_jet']['ordinary_q_squared_rows'][k],restore(value));counts+=1
    assert reduced==2
    for k,value in got['values'].items():same(c,saved['whole_source_paired_derivative_hulls'][k],value)
    old=trace['accepted_uniform_baseline_cover']['selected_primitive_C0_Z_covers']
    selected=trace['selected_primitive_C0_Z_covers']
    for k in ('A','B_over_Pstar'):assert old[k]==selected[k]
    for k,value in got['values'].items():
        a=source.baseline.magnitude_log(restore(old[k]));b=source.baseline.magnitude_log(value)
        choose=b is None or a is not None and b<a
        same(c,selected[k],value if choose else restore(old[k]))
    return dict(original_active_kappa_source_identity_verified=True,
        nonempty_branches=[r['name'] for r in rows],strict_same_source_Delta_Z_reductions=reduced,
        reconstructed_branch_root_q_q_squared_and_primitive_rows=counts,
        original_a_Z_b_Z_and_remaining_mixed_source_rows_retained=True,
        a_C0_and_Delta_Z_are_conditional_function_covers_not_selector_derivatives=True)


@contextmanager
def new_source_check_binding():
    original=history.firstbranch_check;history.firstbranch_check=active_firstbranch_check
    try:yield
    finally:history.firstbranch_check=original


def run():
    begin=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['original_active_kappa_theorem']==source.active_kappa_theorem()
    assert all(manifest[k] is False for k in common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    previous=json.loads((HERE/source.previous.NAME).read_bytes())
    old_O2=json.loads((HERE/common.NAME).read_bytes());old_replay=json.loads((HERE/source.previous.REPLAY).read_bytes())
    c=MPIntervalContext();c.dps=240;checks=[];replays=[];hashes={**manifest['input_hashes'],source.NAME:sha(source.NAME)}
    with mp.workdps(300),new_source_check_binding():
        for archive in manifest['genuine_active_kappa_upstream_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            payload=json.loads(raw);hashes[archive['filename']]=sha(archive['filename'])
            original=next(a for a in previous['genuine_firstbridge_paired_branch_upstream_archives'] if a['exact_Z_range']==payload['exact_Z_range'])
            before=json.loads(gzip.decompress((HERE/original['filename']).read_bytes()))
            checks.append(history.tile_check(c,payload,before))
            for record in old_O2['original_complete_same_N_source_incoming_O2_integral_refinements']:
                if tuple(replay.Fraction(v) for v in record['exact_Z_range'])!=tuple(replay.Fraction(v) for v in payload['exact_Z_range']):continue
                row=accepted.replay_checks(c,payload,record,archive)
                row['actual_upstream_binding'].update(manifest=source.NAME,manifest_sha256=sha(source.NAME),
                    original_active_kappa_Delta_Z_source_cover_refinement=True)
                old=next(r for r in old_replay['genuine_paired_branch_original_O2_replays'] if r['exact_Z_range']==row['exact_Z_range'] and r['ordered_source_cells']==row['ordered_source_cells'])
                coordinates,*unused=replay.inputs(c,payload);comparisons={}
                for kind in ('C0','Z'):
                    for k in KEYS:
                        field='propagated_actual_correction_'+kind
                        a=common.restore_common_source(old[field][k],coordinates);b=common.restore_common_source(row[field][k],coordinates)
                        la=source.baseline.magnitude_log(a);lb=source.baseline.magnitude_log(b)
                        comparisons[kind+'_'+k]=dict(before_absolute_upper_log=la,after_absolute_upper_log=lb,
                            strict_upper_reduction=lb<la,formal_enclosure_not_physical_amplitude_or_residual=True)
                        if kind=='C0':same(c,old[field][k],b)
                        else:assert lb<la
                row['bound_comparisons_against_accepted_paired_branch_replay']=comparisons;replays.append(row)
            print('Original active kappa source/serial C1/O2 replay PASS',payload['exact_Z_range'],flush=True)
    assert len(checks)==2 and len(replays)==4
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    report=dict(**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        genuine_active_kappa_original_O2_replays=replays,original_C0_contributions_backgrounds_and_P0_unchanged=True,
        accepted_O2_integrals_reused_without_recomputation=True,input_hashes=hashes,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.REPLAY).write_text(json.dumps(source.encode(report),indent=2)+'\n',encoding='utf8')
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        genuine_active_kappa_tile_checks=checks,accepted_original_O2_replays=4,strict_downstream_ordinary_Z_reductions=20,
        input_hashes={**hashes,source.REPLAY:sha(source.REPLAY)},execution_seconds=time.monotonic()-begin,
        full_Rc_closure_or_global_N_or_stress_or_recursion_claimed=False,**dict.fromkeys(common.current.FLAGS,False))
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original active kappa conditional derivative refinement PASS',flush=True)
    return result


if __name__=='__main__':run()
