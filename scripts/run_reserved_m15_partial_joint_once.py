"""One-entry partial-support source correlation using saved interval records."""
from __future__ import annotations
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
COMPONENT = 'reserved_m15_partial_joint_v1'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / n for n in ('source.json', 'entered.json', 'result.json'))
BASE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_'
PARTIAL = 'tmp/reserved_m15_actual_cell_v1/result.json'
REPAIRED = BASE + 'repaired_Rp_m15_numeric.json'
JOINT = BASE + 'joint_reserved_m15_numeric.json'
OLD_SOURCE = 'src/openai_ns_reconstruction/original_reserved_m15_response.py'
OLD_RUNNER = 'scripts/run_original_repaired_Rp_m15.py'
LEAF = 'src/openai_ns_reconstruction/reserved_m15_partial_joint.py'
RUNNER = 'scripts/run_reserved_m15_partial_joint_once.py'
REQUIRED = frozenset((PARTIAL, REPAIRED, JOINT, OLD_SOURCE, OLD_RUNNER, LEAF, RUNNER))
FIXED = {PARTIAL: 'aaeed53501a5f8b33e9c052f17dad20574745e20c0956208ae5685164eb8f888',
         REPAIRED: 'bd9fb2fc23072c743eeeadd300041ba187aff908a3bff4a3399f7222b72e6871',
         JOINT: 'a46413fd755aa96c385aa7b81e1919fe0d22838868b8b5a7d4348708dd16dad6'}

def sha(raw): return hashlib.sha256(raw).hexdigest()
def exclusive(path, data):
    raw = (json.dumps(data, sort_keys=True, separators=(',', ':')) + '\n').encode()
    with path.open('xb') as f: f.write(raw)
    return sha(raw)
def require(test, message):
    if not test: raise ValueError(message)

def preflight():
    require(not ENTRY.exists() and not RESULT.exists(), 'already entered')
    braw = BRIEF.read_bytes(); b = json.loads(braw)
    require(b.get('component') == COMPONENT and b.get('phase') == COMPONENT+'_pending'
            and b.get('goal_active') is True and b.get('full_goal_complete') is False
            and b.get('one_entry_only') is True and b.get('do_not_reenter') is True, 'brief guard')
    cp = b['source_checkpoint']; require(Path(cp['path']).resolve() == SOURCE.resolve(), 'checkpoint path')
    sraw = SOURCE.read_bytes(); require(sha(sraw) == cp['sha256'], 'checkpoint hash')
    src = json.loads(sraw)
    require(src.get('component') == COMPONENT and src.get('one_entry_only') is True, 'checkpoint scope')
    hashes = {}
    for row in src['files']:
        name = row['path']; require(name in REQUIRED and name not in hashes, 'unexpected frozen input')
        h = sha((ROOT/name).read_bytes())
        require(h == row['sha256'] and (name not in FIXED or h == FIXED[name]), 'input hash '+name)
        hashes[name] = h
    require(set(hashes) == REQUIRED, 'incomplete manifest')
    return hashes, sha(sraw), sha(braw)

def consume(hashes):
    raw = {}
    for name in REQUIRED:
        data = (ROOT/name).read_bytes(); require(sha(data) == hashes[name], 'changed after entry '+name)
        raw[name] = data
    prior, report, joint = (json.loads(raw[n]) for n in (PARTIAL, REPAIRED, JOINT))
    require(prior['status'] == 'returned' and prior['assembly_invocations'] == 1, 'prior saved result')
    partial = prior['result']; scope = partial['selected_partial_row_scope']
    require(scope['x'] == '1.27' and scope['query_Z'] == '0', 'partial scope')
    family = report['source_family']
    require(partial['source_family'] == family == joint['source_family'], 'source family')
    require(report['checks']['passed'] is True and joint['checks']['passed'] is True, 'accepted reports')
    require(report['implementation_hashes'][OLD_SOURCE] == hashes[OLD_SOURCE]
            and report['implementation_hashes'][OLD_RUNNER] == hashes[OLD_RUNNER], 'old source pin')
    full_name = 'lei_ren_part1_paper_compliant_current_original_n1_full_five_numeric_bridge.json.gz'
    full_hash = report['input_hashes'][full_name]
    require(prior['input_hashes']['experiments/root_st073/'+full_name] == full_hash
            and joint['cached_component_source']['report_sha256'] == full_hash, 'different complete source')
    j = joint['joint_reserved_minus_target_source']
    require(j['source_identity_before_boxing'] is True and j['independent_target_consumed_once'] is True
            and j['ordinary_Z_orders'] == [0,1] and j['exact_response'] == [0,0,1]
            and j['selected_compact_x'] == '2' and j['common_normalization']['U_Z_exact_zero'] is True,
            'saved full-support source theorem')
    samples = [s for s in report['samples'] if s['Z'] == '0']
    views = [s for s in joint['joint_discrepancy_views'] if s['Z'] == '0']
    require(len(samples) == len(views) == 1, 'one Z=0 saved packet/view')
    target = samples[0]['independent_target']; old_target = views[0]['independent_target']
    require(target['actual_function_binding']['sign'] == -1
            and target['coefficient_C0_C1'] == old_target['coefficient_C0_C1']
            and target['log_scale_terms'] == old_target['log_scale_terms'], 'independent target/unit mismatch')
    tree = ast.parse(raw[OLD_SOURCE].decode('utf-8'))
    classes = [x for x in tree.body if isinstance(x, ast.ClassDef) and x.name == 'OriginalReservedM15Response']
    methods = [x for x in classes[0].body if isinstance(x, ast.FunctionDef) and x.name == 'partial_row']
    require(len(classes) == len(methods) == 1 and [a.arg for a in methods[0].args.args] == ['self','x']
            and not any(isinstance(x, ast.Name) and x.id == 'Z' for x in ast.walk(methods[0])), 'partial row scope')
    from mpmath.ctx_iv import MPIntervalContext
    import mpmath as mp
    iv = MPIntervalContext(); iv.dps = 130
    def decode(value):
        pairs = value.get('exact_mpf_tuples', value.get('source_enclosure_exact_mpf_tuples'))
        require(isinstance(pairs, list) and len(pairs) == 2, 'exact interval endpoints')
        return iv.make_mpf(tuple(tuple(x) for x in pairs))
    def packet(v):
        return {'source_enclosure_exact_mpf_tuples':[list(x) for x in v._mpi_],
                'display':[mp.nstr(mp.make_mpf(x), 20) for x in v._mpi_]}
    def hull(values):
        require(bool(values), 'empty hull')
        lo = min((x._mpi_[0] for x in values),key=mp.make_mpf)
        hi = max((x._mpi_[1] for x in values),key=mp.make_mpf)
        return iv.make_mpf((lo,hi))
    response = partial['authenticated_input_intervals']['angular_inverse_response']
    require(len(response) == 3, 'response length')
    cells = [{'cell_id':c['cell_id'],'own_C0_C1':[decode(x) for x in c['own_swirl_C0_C1']]}
             for c in partial['assembly']['cell_results']]
    require(len(cells) == 2, 'expected two saved Z=0 cells')
    spec = importlib.util.spec_from_file_location('reserved_m15_partial_joint_v1_leaf',ROOT/LEAF)
    require(spec is not None and spec.loader is not None, 'leaf import')
    module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)
    value = module.assemble_partial_joint(r2=decode(response[2]), source_cells=cells,
        target_C0_C1=[decode(x) for x in target['coefficient_C0_C1']],
        joint_C0_C1=[decode(x) for x in j['coefficient_C0_C1']], hull=hull,
        source_identity='saved old source F=inner+lower+preRp+pulse', fixed_unit=True,
        r2_z_constant=True, error_ledger={'prior_result_sha256':hashes[PARTIAL],
        'repaired_report_sha256':hashes[REPAIRED], 'joint_report_sha256':hashes[JOINT]})
    return {'query':{'Z':'0','compact_x':'1.27'},'source_family':family,
        'formula':'own+r2*joint-(1-r2)*target',
        'r2':packet(decode(response[2])), 'target_C0_C1':[packet(decode(x)) for x in target['coefficient_C0_C1']],
        'joint_C0_C1':[packet(decode(x)) for x in j['coefficient_C0_C1']],
        'cell_results':[{'cell_id':c['cell_id'], 'own_C0_C1':[packet(x) for x in c['own_C0_C1']],
          'joint_term_C0_C1':[packet(x) for x in c['joint_term_C0_C1']],
          'target_term_C0_C1':[packet(x) for x in c['target_term_C0_C1']],
          'residual_C0_C1':[packet(x) for x in c['residual_C0_C1']]} for c in value['cell_results']],
        'hull_C0_C1':[packet(x) for x in value['hull_C0_C1']], 'error_ledger':value['error_ledger'],
        'unit_and_source_semantics':'conditional on accepted old common-source theorem and fixed-unit interpretation',
        'current_physical_F_owner_identity_proved':False, 'forward_212_compared':False,
        'independent_terminal_matching_closed':False, 'old_errors_removed':False}

def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {'component':COMPONENT,'entered':True,'source_sha256':source_sha,
                                  'brief_sha256':brief_sha,'input_hashes':hashes,'consumer_calls':0})
    out = {'component':COMPONENT,'entry_sha256':entry_sha,'source_sha256':source_sha,
           'status':'failed','consumer_calls':0}
    try:
        out['consumer_calls']=1; out['result']=consume(hashes); out['status']='returned'; out['error']=None
    except BaseException as exc:
        out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    result_sha=exclusive(RESULT,out)
    print(json.dumps({'component':COMPONENT,'status':out['status'],'source_sha256':source_sha,
                      'entry_sha256':entry_sha,'result_sha256':result_sha}))
    return 0 if out['status']=='returned' else 1
if __name__ == '__main__': raise SystemExit(main())
