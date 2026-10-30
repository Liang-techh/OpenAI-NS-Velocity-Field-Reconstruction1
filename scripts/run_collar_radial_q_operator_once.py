"""One-entry source-pinned audit of the actual collar radial Q operators."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback

ROOT=Path(__file__).resolve().parents[1]
BRIEF=ROOT.parent.parent/'RESUME_BRIEF_2026_10_16.json'
COMPONENT='collar_radial_q_operator_v1'
DIR=ROOT/'tmp'/COMPONENT
SOURCE,ENTRY,RESULT=(DIR/name for name in ('source.json','entered.json','result.json'))
BASE='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_'
LOWER=BASE+'n1_lower_collar_functions.py'
CORE=BASE+'n1_core_fixed_point_functions.py'
PULSE=BASE+'n1_pulse_pressure_functions.py'
OWN=BASE+'n1_pulse_own_fields.py'
C2=BASE+'C2_phase_primitives.py'
FORWARD='src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
LAURENT='src/openai_ns_reconstruction/original_first_collar_transport.py'
PRIOR='experiments/root_st073/collar_pressure_regular_singular_v6_receipt.json'
LEAF='src/openai_ns_reconstruction/collar_radial_q_operator_proof.py'
RUNNER='scripts/run_collar_radial_q_operator_once.py'
REQUIRED=frozenset((LOWER,CORE,PULSE,OWN,C2,FORWARD,LAURENT,PRIOR,LEAF,RUNNER))
FIXED={PRIOR:'c4c56da7f0a06a776f87498505950a5f2f077f710eb7fb2d270081ce845ec8fd'}

def sha(b):return hashlib.sha256(b).hexdigest()
def require(value,message):
    if not value:raise ValueError(message)
def exclusive(path,value):
    raw=(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n').encode()
    with path.open('xb') as stream:stream.write(raw)
    return sha(raw)
def preflight():
    require(not ENTRY.exists() and not RESULT.exists(),'already entered')
    braw=BRIEF.read_bytes();b=json.loads(braw)
    require(b.get('component')==COMPONENT and b.get('phase')==COMPONENT+'_pending'
        and b.get('goal_active') is True and b.get('full_goal_complete') is False
        and b.get('do_not_reenter') is True and b.get('one_entry_only') is True,'pending guard')
    cp=b['source_checkpoint'];require(Path(cp['path']).resolve()==SOURCE.resolve(),'checkpoint path')
    sraw=SOURCE.read_bytes();require(sha(sraw)==cp['sha256'],'checkpoint SHA')
    src=json.loads(sraw);require(src.get('component')==COMPONENT and src.get('one_entry_only') is True,'scope')
    hashes={}
    for item in src['files']:
        path=item['path'];require(path in REQUIRED and path not in hashes,'unexpected input')
        h=sha((ROOT/path).read_bytes());require(h==item['sha256'] and (path not in FIXED or FIXED[path]==h),'changed '+path)
        hashes[path]=h
    require(set(hashes)==REQUIRED,'missing input')
    return hashes,sha(sraw),sha(braw)
def _contains_statement(text, statement):
    wanted=ast.dump(ast.parse(statement).body[0],include_attributes=False)
    return any(ast.dump(n,include_attributes=False)==wanted for n in ast.parse(text).body)
def consume(hashes):
    raw={}
    for path in REQUIRED:
        data=(ROOT/path).read_bytes();require(sha(data)==hashes[path],'changed after entry '+path);raw[path]=data
    old=json.loads(raw[PRIOR]);require(old['status']=='returned' and old['component']=='collar_pressure_regular_singular_v6','prior receipt')
    require(old['input_hashes'][FORWARD]==hashes[FORWARD] and old['input_hashes'][LOWER]==hashes[LOWER], 'prior source pins')
    alias_checks={
        'lower_core_import':_contains_statement(raw[LOWER].decode(), 'import lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions as core'),
        'lower_own_alias':_contains_statement(raw[LOWER].decode(), 'own=core.pulse.own'),
        'core_pulse_import':_contains_statement(raw[CORE].decode(), 'import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions as pulse'),
        'pulse_own_import':_contains_statement(raw[PULSE].decode(), 'import lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields as own'),
        'own_C2_import':_contains_statement(raw[OWN].decode(), 'import lei_ren_part1_paper_compliant_current_original_C2_phase_primitives as c2')}
    require(all(alias_checks.values()),'actual own/C2 import chain changed')
    spec=importlib.util.spec_from_file_location('collar_radial_q_operator_v1_leaf',ROOT/LEAF)
    require(spec is not None and spec.loader is not None,'leaf import')
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    proof=module.prove(*(raw[path].decode('utf-8') for path in (OWN,FORWARD,C2,LAURENT)))
    return {'actual_alias_chain_checks':alias_checks,'operator_proof':proof,
        'prior_conditional_comparison_result_sha256':old['result_sha256'],
        'prior_other_assumptions_unchanged':True,'leading_field_and_history_same_owner_proved':False,
        'forward_finite_rows_equal_canonical_function_proved':False,
        'error_removed':False,'numeric_value_computed':False,'terminal_matching_proved':False}
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
