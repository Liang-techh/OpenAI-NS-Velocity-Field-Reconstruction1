"""One-entry source-bound Cauchy and raw-operator error propagation."""
from __future__ import annotations
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / 'RESUME_BRIEF_2026_10_16.json'
COMPONENT = 'n1_seam_z1_cauchy_v1'
DIR = ROOT / 'tmp' / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ('source.json', 'entered.json', 'result.json'))
BASE = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_'
REPORT = BASE + 'selected_words.json'
COUPLED = BASE + 'coupled_first_step.py'
PRIOR = 'experiments/root_st073/n1_inlet_residual_join_v1_receipt.json.gz'
DOMAIN_DOC = 'docs/CURRENT_ORIGINAL_N1_ACTUAL_COLLAR_DOMAIN_2026_10_10.md'
COUPLED_DOC = 'docs/CURRENT_ORIGINAL_N1_COUPLED_FIRST_STEP_2026_10_11.md'
LEAF = 'src/openai_ns_reconstruction/n1_seam_z1_error_transport.py'
RUNNER = 'scripts/run_n1_seam_z1_cauchy_once.py'
REQUIRED = (REPORT, COUPLED, PRIOR, DOMAIN_DOC, COUPLED_DOC, LEAF, RUNNER)
FIXED = {
    REPORT: '0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e',
    COUPLED: '97913b82b8570b220921cd89918569a4019ff3495bdbea0548b3c53e9d8c93e7',
    PRIOR: '6eb7552ea8253c6ee4ab8a0f0ed507fc8c24b5e083013b0a4d4dac865ce10c5c',
}


def sha(raw): return hashlib.sha256(raw).hexdigest()
def require(ok, message):
    if not ok: raise ValueError(message)


def exclusive(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()
    with path.open('xb') as handle: handle.write(raw)
    return sha(raw)


def preflight():
    require(not ENTRY.exists() and not RESULT.exists(), 'already entered')
    brief_raw = BRIEF.read_bytes(); brief = json.loads(brief_raw)
    require(brief.get('component') == COMPONENT and brief.get('phase') == COMPONENT+'_pending'
            and brief.get('goal_active') is True and brief.get('full_goal_complete') is False
            and brief.get('one_entry_only') is True and brief.get('do_not_reenter') is True, 'brief guard')
    checkpoint = brief['source_checkpoint']; raw = SOURCE.read_bytes()
    require(Path(checkpoint['path']).resolve() == SOURCE.resolve() and sha(raw) == checkpoint['sha256'], 'checkpoint')
    manifest = json.loads(raw); require(manifest['component'] == COMPONENT and manifest['one_entry_only'] is True, 'manifest')
    hashes = {}
    for item in manifest['files']:
        path = item['path']; require(path in REQUIRED and path not in hashes, 'unexpected input')
        digest = sha((ROOT/path).read_bytes())
        require(digest == item['sha256'] and (path not in FIXED or digest == FIXED[path]), 'changed input '+path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), 'missing input')
    return hashes, sha(raw), sha(brief_raw)


def interval(c, record):
    lower = tuple(record['lower_exact_mpf_tuple']); upper = tuple(record['upper_exact_mpf_tuple'])
    require(len(lower) == len(upper) == 4, 'four-entry endpoints required')
    x = c.make_mpf((lower, upper)); require(x.a <= x.b, 'reversed interval')
    return x


def saved_pair(c, pair):
    # Prior runner serialized lower/upper endpoint intervals independently.
    lower = tuple(pair['lower']['lower_exact_mpf_tuple'])
    upper = tuple(pair['upper']['upper_exact_mpf_tuple'])
    x = c.make_mpf((lower, upper)); require(x.a <= x.b, 'reversed sector coefficient')
    return x


def encode(value):
    if hasattr(value, '_mpi_'):
        return {'lower_exact_mpf_tuple': list(value._mpi_[0]), 'upper_exact_mpf_tuple': list(value._mpi_[1]), 'directed': True}
    if isinstance(value, dict): return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)): return [encode(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None: return value
    raise TypeError('unsupported output '+type(value).__name__)


def consume(hashes):
    raw = {}
    for path in REQUIRED:
        data = (ROOT/path).read_bytes(); require(sha(data) == hashes[path], 'changed after entry '+path); raw[path] = data
    report = json.loads(raw[REPORT]); previous_packet = json.loads(gzip.decompress(raw[PRIOR]))
    require(report['checks']['passed'] is True and previous_packet['result']['status'] == 'returned', 'accepted inputs')
    require(previous_packet['frozen_byte_sha256']['result'] == '662ec434484c37ec87e6abfe20acfe57ab0b2527d5a0ade31adf5faa99a1ce69', 'prior result hash')
    previous_output = previous_packet['result']['output']; previous = previous_output['residual_enclosure']
    tail = report['actual_selected_raw_word_tail']; arithmetic = tail['arithmetic']
    require(tail['passed'] is True and tail['tail_is_a_directed_majorant_not_a_field_value'] is True
            and tail['scope'] == 'uniform N=1 raw core/first-collar tail majorant'
            and tail['natural_word_depth'] == 1 and tail['ordinary_Z_orders'] == [0, 1]
            and tail['normalized_error_mapping']['raw_tail_only'] is True
            and tail['normalized_error_mapping']['applied'] is False, 'tail scope')
    require(tail['source_family'] == previous_output['source_family'] == report['source_family']
            and previous_output['selected_Z_interval'] == report['actual_selected_first_word_enclosure']['selected_Z_interval'], 'same family and Z')
    require(previous_output['input_binding']['same_four_log_anchors'] is True
            and previous_output['input_binding']['rho_is_4_and_phase_zero'] is True, 'previous physical inlet binding')
    require(arithmetic['formula'] == 'Mbase_upper * beta_N_upper^(N+1) / (1-beta_N_upper)'
            and arithmetic['retained_preconditioned_indices'] == [0, 1]
            and arithmetic['positive_Mbase_and_beta_below_one'] is True
            and arithmetic['Cauchy_gap_positive'] is True, 'preconditioned norm metadata')
    pin = tail['source_pins']['coupled_source']
    require(pin['path'] == COUPLED and pin['sha256_match'] is True and pin['expected_sha256'] == hashes[COUPLED], 'coupled pin')
    coupled = raw[COUPLED].decode('utf-8'); domain_doc = raw[DOMAIN_DOC].decode('utf-8'); coupled_doc = raw[COUPLED_DOC].decode('utf-8')
    require(all(x in coupled for x in (
        "arbitrary_k_term_formula='Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!'",
        'beta=domain.upper(c,c.e*a*S*(Cp/(N+1)+C1/loss))',
        'initial_function_Z_radius=eta/4,solution_Z_radius=eta/8,matrix_Z_radius=eta/2',
        "derivative_margin=tail['solution_Z_radius']/2")), 'word and gap source formulas')
    require('Radius around the real Z interval' in domain_doc and 'Recursive solution/tail estimate | eta/8' in domain_doc
            and '||H^k Wbase||_(eta/8)' in coupled_doc, 'common tube and W1 norm proof')
    c = MPIntervalContext(); c.prec = 1800
    z = interval(c, previous_output['selected_Z_interval'])
    require(z.a >= -1 and z.b <= 1, 'selected real Z tube outside admitted real interval')
    mbase = interval(c, arithmetic['Mbase_upper']).b
    beta = interval(c, arithmetic['beta_N_upper']).b
    gap = interval(c, tail['raw_tail_enclosures']['0']['cauchy_gap']).a
    gap1 = interval(c, tail['raw_tail_enclosures']['1']['cauchy_gap']).a
    t0 = interval(c, tail['raw_tail_enclosures']['0']['tail_upper']).b
    t1 = interval(c, tail['raw_tail_enclosures']['1']['tail_upper']).b
    require(mbase >= 0 and beta >= 0 and beta < 1 and gap > 0 and gap == gap1 and t0 >= 0 and t1 >= 0, 'positive bounds')
    # k=1 term: a(1+a)(Cp+C1/L) <= beta1, since e/2>1.
    # For each real Z in [-1,1], the disk of radius gap=eta/16 is
    # contained in the eta/8 tube around that same real interval.
    word_norm = (mbase*beta).b
    word_zz = (2*word_norm/(gap*gap)).b
    tail_zz = (2*t0/(gap*gap)).b
    require(word_norm >= 0 and word_zz >= 0 and tail_zz >= 0, 'negative derived radius')

    def absolute_product_upper(coefficient, bound):
        value = saved_pair(c, coefficient)
        result = (abs(value).b * bound).b
        require(result >= 0, 'negative operator-sector radius')
        return result

    spec = importlib.util.spec_from_file_location('n1_seam_z1_transport_leaf', ROOT/LEAF)
    require(spec is not None and spec.loader is not None, 'leaf import')
    module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)
    propagated = module.propagate(previous, word_zz_upper=word_zz, tail_upper=t0, tail_z_upper=t1,
        tail_zz_upper=tail_zz, absolute_product_upper=absolute_product_upper)
    return {'source_family': report['source_family'], 'selected_Z_interval': previous_output['selected_Z_interval'],
        'cauchy_scope': {'base_real_interval': '[-1,1]', 'solution_tube_radius': 'eta/8',
          'real_center_disk_radius': 'eta/16', 'tube_inclusion': 'each real center belongs to the base interval; its eta/16 disk lies inside its eta/8 neighborhood',
          'W1_norm_formula': 'Mbase*beta1; e*(Cp/2+C1/L) >= Cp+C1/L',
          'second_derivative_formula': '2*norm/(eta/16)^2', 'positive_bound_is_not_signed_value': True},
        'directed_bounds': encode({'Mbase_upper': mbase, 'beta1_upper': beta, 'gap_lower': gap,
          'W1_norm_upper': word_norm, 'W1_Z2_upper': word_zz, 'tail_Z0_upper': t0,
          'tail_Z1_upper': t1, 'tail_Z2_upper': tail_zz}),
        'error_transport': encode(propagated), 'Z1_is_enclosure_not_a_signed_point': True,
        'physical_NS_residual_or_terminal_closure': False}


def main():
    hashes, source_sha, brief_sha = preflight()
    entry_sha = exclusive(ENTRY, {'component': COMPONENT, 'source_sha256': source_sha, 'brief_sha256': brief_sha,
        'input_hashes': hashes, 'entered': True})
    result = {'component': COMPONENT, 'source_sha256': source_sha, 'entry_sha256': entry_sha, 'status': 'failed'}
    try:
        result['output'] = consume(hashes); result['status'] = 'returned'; result['error'] = None
    except BaseException as exc:
        result['error'] = {'type': type(exc).__name__, 'message': str(exc), 'traceback': traceback.format_exc()}
    result_sha = exclusive(RESULT, result)
    print(json.dumps({'component': COMPONENT, 'status': result['status'], 'result_sha256': result_sha}))
    return 0 if result['status'] == 'returned' else 1


if __name__ == '__main__': raise SystemExit(main())
