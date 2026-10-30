"""One-entry exact metadata join for the saved Cauchy tube and gap."""
from __future__ import annotations
import gzip
import hashlib
import json
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
COMPONENT = 'n1_seam_gap_binding_v1'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR/name for name in ('source.json', 'entered.json', 'result.json'))
BASE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_'
REPORT = BASE + 'selected_words.json'
ARCHIVE = BASE + 'coupled_first_step.json.gz'
PRIOR = 'tmp/n1_seam_z1_cauchy_v1/result.json'
RUNNER = 'scripts/run_n1_seam_gap_binding_once.py'
REQUIRED = (REPORT, ARCHIVE, PRIOR, RUNNER)
FIXED = {REPORT:'0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',
         ARCHIVE:'30b807af4147658ec0978ebd253cc1b2b31bed67bfed430a965e66461450add1',
         PRIOR:'bf6ad485216d237166af8ea55261650925fbd3b64fee267af6c9bdd2910caaac'}


def sha(raw): return hashlib.sha256(raw).hexdigest()
def require(ok, message):
    if not ok: raise ValueError(message)
def exclusive(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':'))+'\n').encode()
    with path.open('xb') as handle: handle.write(raw)
    return sha(raw)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(), 'already entered')
    brief_raw = BRIEF.read_bytes(); brief = json.loads(brief_raw)
    require(brief.get('component') == COMPONENT and brief.get('phase') == COMPONENT+'_pending'
        and brief.get('goal_active') is True and brief.get('full_goal_complete') is False
        and brief.get('one_entry_only') is True and brief.get('do_not_reenter') is True, 'brief guard')
    checkpoint=brief['source_checkpoint']; raw=SOURCE.read_bytes()
    require(Path(checkpoint['path']).resolve() == SOURCE.resolve() and sha(raw) == checkpoint['sha256'], 'checkpoint')
    source=json.loads(raw); require(source['component'] == COMPONENT and source['one_entry_only'] is True, 'manifest')
    hashes={}
    for row in source['files']:
        path=row['path']; require(path in REQUIRED and path not in hashes, 'unexpected input')
        digest=sha((ROOT/path).read_bytes())
        require(digest == row['sha256'] and (path not in FIXED or digest == FIXED[path]), 'changed '+path)
        hashes[path]=digest
    require(set(hashes) == set(REQUIRED), 'missing input')
    return hashes,sha(raw),sha(brief_raw)


def endpoints(record):
    a,b=tuple(record['lower_exact_mpf_tuple']),tuple(record['upper_exact_mpf_tuple'])
    require(len(a)==len(b)==4, 'four-entry endpoint tuples required')
    return a,b


def consume(hashes):
    raw={}
    for path in REQUIRED:
        data=(ROOT/path).read_bytes(); require(sha(data)==hashes[path], 'changed after entry '+path); raw[path]=data
    report=json.loads(raw[REPORT]); archive=json.loads(gzip.decompress(raw[ARCHIVE])); prior=json.loads(raw[PRIOR])
    require(report['checks']['passed'] is True and prior['status']=='returned', 'accepted inputs')
    tail=report['actual_selected_raw_word_tail']; sample=archive['samples']['actual_first_keep']
    actual=sample['analytic_preconditioned_tail']; derivative=sample['derivative_tail']; arithmetic=tail['arithmetic']
    require(tail['passed'] is True and tail['source_family']==sample['source_family']==archive['source_family']==prior['output']['source_family'], 'source family')
    require(tail['source_pins']['coupled_archive']['expected_sha256']==hashes[ARCHIVE]
        and tail['source_pins']['coupled_archive']['sha256_match'] is True, 'coupled archive pin')
    checks={
        'Mbase_same': endpoints(arithmetic['Mbase_upper'])==endpoints(actual['initial_base_common_tube_bound']),
        'beta1_same': endpoints(arithmetic['beta_N_upper'])==endpoints(actual['geometric_tail_ratio_upper']),
        'tail_Z0_same': endpoints(tail['raw_tail_enclosures']['0']['tail_upper'])==endpoints(actual['omitted_sum_bound']),
        'tail_Z1_same': endpoints(tail['raw_tail_enclosures']['1']['tail_upper'])==endpoints(derivative['remaining_first_Z_derivative_norm_upper']),
        'gap_same': endpoints(tail['raw_tail_enclosures']['0']['cauchy_gap'])==endpoints(derivative['Cauchy_gap'])==endpoints(tail['raw_tail_enclosures']['1']['cauchy_gap']),
        'scope': actual['scope']=='uniform analytic series tail for canonical core/actual-first-collar problem; later completed-leading supports unbound',
    }
    c=MPIntervalContext(); c.prec=1800
    radius=c.make_mpf(endpoints(actual['solution_Z_radius'])); gap=c.make_mpf(endpoints(derivative['Cauchy_gap']))
    half_radius=radius/2
    checks['gap_is_half_solution_radius']=gap._mpi_==half_radius._mpi_
    checks['source_gap_positive']=bool(gap.a>0)
    checks['prior_used_gap_lower']=tuple(prior['output']['directed_bounds']['gap_lower']['lower_exact_mpf_tuple'])==endpoints(derivative['Cauchy_gap'])[0]
    require(all(checks.values()), 'exact Cauchy metadata join failed')
    return {'checks':checks,'source_family':tail['source_family'],
        'solution_radius':actual['solution_Z_radius'],'cauchy_gap':derivative['Cauchy_gap'],
        'gap_identity':'saved Cauchy_gap = solution_Z_radius/2 = eta/16',
        'same_real_base_interval':'[-1,1]', 'prior_result_sha256':hashes[PRIOR],
        'old_math_reentered':False, 'new_field_value_computed':False}


def main():
    hashes,source_sha,brief_sha=preflight()
    entry=exclusive(ENTRY,{'component':COMPONENT,'source_sha256':source_sha,'brief_sha256':brief_sha,'input_hashes':hashes,'entered':True})
    out={'component':COMPONENT,'source_sha256':source_sha,'entry_sha256':entry,'status':'failed'}
    try:
        out['output']=consume(hashes); out['status']='returned'; out['error']=None
    except BaseException as exc:
        out['error']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
    digest=exclusive(RESULT,out); print(json.dumps({'component':COMPONENT,'status':out['status'],'result_sha256':digest}))
    return 0 if out['status']=='returned' else 1


if __name__=='__main__': raise SystemExit(main())
