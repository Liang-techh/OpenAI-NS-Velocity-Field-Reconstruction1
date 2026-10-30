"""Link saved native collar theorem to the independently checked Q operator."""
import hashlib
import json
from pathlib import Path
import traceback

ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='collar_native_operator_link_v1'
DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/name for name in ('source.json','entered.json','result.json'))
NATIVE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_collar_transport.json'
Q='experiments/root_st073/collar_radial_q_operator_v1_receipt.json'
REGULAR='experiments/root_st073/collar_pressure_regular_singular_v6_receipt.json'
LOWER='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py'
RUNNER='scripts/run_collar_native_operator_link_once.py'
REQUIRED=frozenset((NATIVE,Q,REGULAR,LOWER,RUNNER))
FIXED={NATIVE:'7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96',
       Q:'c8fc8ff7c58c54f3ce6438b3bb661b6bb75c8fe680d844ebd72b619cbada0fe0',
       REGULAR:'c4c56da7f0a06a776f87498505950a5f2f077f710eb7fb2d270081ce845ec8fd'}
def sha(raw):return hashlib.sha256(raw).hexdigest()
def require(value,message):
    if not value:raise ValueError(message)
def exclusive(path,value):
    raw=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as f:f.write(raw)
    return sha(raw)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    braw=BRIEF.read_bytes();b=json.loads(braw)
    require(b.get('component')==COMPONENT and b.get('phase')==COMPONENT+'_pending'
        and b.get('goal_active') is True and b.get('full_goal_complete') is False
        and b.get('do_not_reenter') is True and b.get('one_entry_only') is True,'pending guard')
    cp=b['source_checkpoint'];require(Path(cp['path']).resolve()==SOURCE.resolve(),'checkpoint path')
    sraw=SOURCE.read_bytes();require(sha(sraw)==cp['sha256'],'checkpoint hash')
    src=json.loads(sraw);require(src.get('component')==COMPONENT and src.get('one_entry_only') is True,'scope')
    hashes={}
    for item in src['files']:
        path=item['path'];require(path in REQUIRED and path not in hashes,'unexpected input')
        h=sha((ROOT/path).read_bytes());require(h==item['sha256'] and (path not in FIXED or h==FIXED[path]),'changed '+path)
        hashes[path]=h
    require(set(hashes)==REQUIRED,'missing input')
    return hashes,sha(sraw),sha(braw)
def consume(hashes):
    raw={}
    for name in REQUIRED:
        data=(ROOT/name).read_bytes();require(sha(data)==hashes[name],'changed after entry '+name);raw[name]=data
    native,q,regular=(json.loads(raw[name]) for name in (NATIVE,Q,REGULAR))
    require(native['checks']['passed'] is True and q['status']=='returned' and regular['status']=='returned','accepted reports')
    t=native['actual_native_collar_transport'];contract=native['contract']
    true=('native_defining_leading_function_equals_canonical_collar','actual_quotient_equals_ratio',
          'comparison_core_primitive_identity','actual_own_mean_and_pressure_equal_canonical_actual_histories',
          'first_window_only','all_old_error_caps_retained')
    false=('native_actual_six_histories_equal_comparison','old_finite_numeric_producer_is_exact_function',
           'raw_n1_cut_response_identity','full_C_inner_m15_identity','six_signed_source_equality',
           'saved_six_terms_or_error_ledger_contracted','global_five_moment_closure')
    require(t['passed'] is True and all(t[k] is True for k in true) and all(t[k] is False for k in false), 'native source scope')
    require(t['ordinary_Z_orders']==[0,1,2,3] and contract['old_finite_rows_not_promoted_to_exact_function_values'] is True,'native derivative or finite scope')
    lower_name=Path(LOWER).name
    require(native['input_hashes'][lower_name]==hashes[LOWER]
        and regular['input_hashes'][LOWER]==hashes[LOWER] and q['input_hashes'][LOWER]==hashes[LOWER], 'canonical lower source mismatch')
    qp=q['proof'];op=qp['operator_proof']
    require(all(qp['actual_alias_chain_checks'].values()) and op['passed'] is True
        and op['differences_by_order']==['0','0','0'],'actual radial Q operator')
    require(qp['prior_conditional_comparison_result_sha256']==regular['result_sha256'],'different conditional proof')
    require(all(v['residual']=='0' for v in regular['proof']['identities'].values())
        and regular['proof']['ordinary_Z1']['all_residuals']=='0','conditional decomposition missing')
    return {'status':'saved_native_leading_and_operator_source_link',
        'native_source_report_sha256':hashes[NATIVE], 'canonical_lower_source_sha256':hashes[LOWER],
        'first_window_native_leading_equals_canonical_from_accepted_report':True,
        'native_actual_mean_and_pressure_histories_retained':True,
        'actual_histories_not_equated_to_comparison_histories':True,
        'ordinary_Z_orders_from_native_theorem':[0,1,2,3],
        'radial_Q_operator_rows_equal':[0,1,2],
        'conditional_regular_singular_source_formula_reused':True,
        'scope':'same original leading function on first-window common collar; not old finite source identity',
        'old_finite_numeric_function_equality_proved':False,'raw_n1_cut_response_identity_proved':False,
        'forward_six_source_equality_proved':False,'old_errors_removed':False,
        'numeric_value_computed':False,'global_terminal_closure_proved':False}
def main():
    hashes,source_sha,brief_sha=preflight()
    entry_sha=exclusive(ENTRY,{'component':COMPONENT,'entered':True,'source_sha256':source_sha,'brief_sha256':brief_sha,'input_hashes':hashes,'consumer_calls':0})
    out={'component':COMPONENT,'source_sha256':source_sha,'entry_sha256':entry_sha,'status':'failed','consumer_calls':0}
    try:
        out['consumer_calls']=1;out['proof']=consume(hashes);out['status']='returned';out['error']=None
    except BaseException as exc:
        out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    result_sha=exclusive(RESULT,out)
    print(json.dumps({'component':COMPONENT,'status':out['status'],'source_sha256':source_sha,'entry_sha256':entry_sha,'result_sha256':result_sha}))
    return 0 if out['status']=='returned' else 1
if __name__=='__main__':raise SystemExit(main())
