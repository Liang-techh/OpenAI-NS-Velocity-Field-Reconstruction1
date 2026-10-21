"""Fresh source-only J3 successor consuming the saved V2 failure certificate.

The V2 attempt is treated as immutable JSON evidence.  This module checks the
saved Boolean facts, repairs the one producer/consumer key mismatch, and
rechecks the bridge fingerprint and the bounded endpoint rows.  It does not
import or call an earlier kernel, an archive producer, or a runtime operator.
"""

import ast as _ast
import hashlib as _hashlib
import json as _json
from collections.abc import Mapping as _Mapping


SCHEMA = 'original_native_J3_endpoint_transport_v3'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_B_source_closure.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
MOD = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
MOD_CREATOR = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
C6_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py'
BRIDGE_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py'
BACKGROUND = 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz'
WEIGHTED = 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json'
WEIGHTED_CREATOR = 'src/openai_ns_reconstruction/original_native_active_S1_weighted_ibp.py'
FAILED_REPORT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport_v2.json'
FAILURE_EVIDENCE = 'experiments/root_st073/original_native_J3_endpoint_transport_failure_evidence.json'
GAMMA_CREATOR = 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py'

# These six paths are the lossless V2 source, entry, result, and diagnostic
# chain.  The two implementation files are included so the saved chain cannot
# be detached from the source that produced it.
V2_SOURCE_CODE = 'src/openai_ns_reconstruction/original_native_J3_endpoint_transport_v2.py'
V2_RUNNER = 'scripts/run_original_native_J3_endpoint_transport_v2.py'
V2_SOURCE_PACKET = 'tmp/original_native_J3_endpoint_transport_v2/new_endpoint_transport_v2_source.json'
V2_ENTERED = 'tmp/original_native_J3_endpoint_transport_v2/new_endpoint_transport_v2_entered.json'
V2_RESULT = 'tmp/original_native_J3_endpoint_transport_v2/new_endpoint_transport_v2_result.json'
V2_DIAGNOSTICS = 'tmp/original_native_J3_endpoint_transport_v2/new_endpoint_transport_v2_diagnostics.jsonl'

RECORD_NAMES = [
    PARENT, LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, C6_CREATOR,
    BRIDGE_CREATOR, BACKGROUND, WEIGHTED, WEIGHTED_CREATOR,
    GAMMA_CREATOR, FAILED_REPORT, FAILURE_EVIDENCE,
    V2_SOURCE_CODE, V2_RUNNER, V2_SOURCE_PACKET, V2_ENTERED,
    V2_RESULT, V2_DIAGNOSTICS,
]
RECORD_HASHES = {
    PARENT: '37f20f105dcb58691b32deb7f5fca59f1a25ff7641086614aa61885ddc090655',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP_CREATOR: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    MOD_CREATOR: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
    C6_CREATOR: '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32',
    BRIDGE_CREATOR: '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605',
    BACKGROUND: '5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a',
    WEIGHTED: 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28',
    WEIGHTED_CREATOR: 'cdf33563c34e211a00128e3901f49875e50bdea3b4cd9224cc00f81e3288ed8d',
    GAMMA_CREATOR: '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01',
    # The V2 report is published after this source is frozen.  The publisher
    # supplies this one literal before entry; all other hashes are immutable.
    FAILED_REPORT: 'a05f11e84353e4a1665d56f4560c3a5745b2bfc7c5ee37eed85b62147a5fbe83',
    FAILURE_EVIDENCE: 'f0b3869589559344ca2fa2513f81104e50654d25e6560bf9474c682fdc82319d',
    V2_SOURCE_CODE: 'fac9283a9a4e87ddaabd7fc0e6563e97e9c15e6e38457b5389503a75e08ab5bc',
    V2_RUNNER: 'eeea5a1f306d1e7f67f27961ee89d5ccf01d6144c6f4dfb7f0f7c9a93614b1b1',
    V2_SOURCE_PACKET: 'e142d66b66d17c1e83d490fe08a2d413806b6413f133680b04d53c9869fa1f79',
    V2_ENTERED: '9319a7dd6fc5bbfd8f6e5e5f6629d63309013a08a7b4ce4f26f16ee5d798a8b9',
    V2_RESULT: '35b0966181cbc5d0cc992b484bf4b0390ca48a9028b17121ab2cbca5686f0900',
    V2_DIAGNOSTICS: '1b8a9af92e65cc362e1894f4147d755630144f5746d8251051f6b471e6c1465d',
}

STAGE_FLAGS = [
    'actual_J3_source_ports_bound',
    'original_five_moments_Z2_endpoint_identity_proved',
    'finite_N_m_Z2_endpoint_identity_proved',
    'J3_R_C1_endpoint_identity_proved',
    'J3_U_C1_endpoint_identity_proved',
    'J3_Q_C1_endpoint_identity_proved',
    'J3_B_C1_endpoint_identity_proved',
    'J3_native_weighted_jump_zero_proved',
    'native_IN_J3_source_reduction_proved',
    'all_old_terms_retained',
    'accepted_ancestors_not_reexecuted',
]
FALSE_FLAGS = [
    'all_native_B_jumps_zero_proved', 'J4_through_J15_B_identities_proved',
    'quantitative_packet_error_reduction_proved', 'numeric_values_computed',
    'J3_Qy_identity_proved', 'source_C5_installation_claimed',
    'numeric_containment_claimed', 'new_loop_domain_admission_claimed',
    'admissible_cone_admission_claimed', 'physical_layer_completed',
    'old_error_terms_removed', 'global_five_moment_closure',
    'time_order_recursion_claimed', 'full_goal_complete',
    'uninterpreted_graph_parameter_identity_proved', 'old_operators_replayed',
]
FAMILY = {
    'actual_five_defect_family_sha256': '3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894',
    'implicit_source_sha256': '5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae',
    'datum_enclosure_sha256': 'dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010',
}

V2_SOURCE_DIGEST = '24034acbfa96218982867140fb8c37d9b5f3b9c82c1c82055a4c637bb2674efc'
V2_FROZEN_INVENTORY_SHA256 = '3eca9d5d9126d67f924ef5f258db8ea4c69c8698a9dd168711726e618b80b02a'
V2_FROZEN_INVENTORY_PATH = 'work/native_J3_v2_frozen_failure_inventory_2026_10_21.json'
V2_FROZEN_INVENTORY_SCHEMA = 'native_J3_v2_frozen_failure_inventory_v1'
V2_FROZEN_INVENTORY_BYTES = 1602
V2_SOURCE_RECORD_HASHES = {
    PARENT: RECORD_HASHES[PARENT], LOOP: RECORD_HASHES[LOOP], MOD: RECORD_HASHES[MOD],
    LOOP_CREATOR: RECORD_HASHES[LOOP_CREATOR], MOD_CREATOR: RECORD_HASHES[MOD_CREATOR],
    C6_CREATOR: RECORD_HASHES[C6_CREATOR], BRIDGE_CREATOR: RECORD_HASHES[BRIDGE_CREATOR],
    BACKGROUND: RECORD_HASHES[BACKGROUND], WEIGHTED: RECORD_HASHES[WEIGHTED],
    WEIGHTED_CREATOR: RECORD_HASHES[WEIGHTED_CREATOR],
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport.json':
        'cc276ae3ab77c8cd6730e4f44b3ea60e95ea4af4b1e50d2153915781b26b16d0',
    FAILURE_EVIDENCE: RECORD_HASHES[FAILURE_EVIDENCE], GAMMA_CREATOR: RECORD_HASHES[GAMMA_CREATOR],
}
V2_SOURCE_DECODED_RECORD_HASHES = {
    PARENT: RECORD_HASHES[PARENT], LOOP: 'f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2',
    MOD: '0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75',
    LOOP_CREATOR: RECORD_HASHES[LOOP_CREATOR], MOD_CREATOR: RECORD_HASHES[MOD_CREATOR],
    C6_CREATOR: RECORD_HASHES[C6_CREATOR], BRIDGE_CREATOR: RECORD_HASHES[BRIDGE_CREATOR],
    BACKGROUND: '18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a',
    WEIGHTED: RECORD_HASHES[WEIGHTED], WEIGHTED_CREATOR: RECORD_HASHES[WEIGHTED_CREATOR],
    'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport.json':
        'cc276ae3ab77c8cd6730e4f44b3ea60e95ea4af4b1e50d2153915781b26b16d0',
    FAILURE_EVIDENCE: RECORD_HASHES[FAILURE_EVIDENCE], GAMMA_CREATOR: RECORD_HASHES[GAMMA_CREATOR],
}
INHERITED_SOURCE_COUNT = 1864
INHERITED_CHECKPOINT_COUNT = 140
FAILED_ARTIFACT_UNION_COUNT = 26
NEW_SOURCE_UNION_COUNT = 1867
NEW_CHECKPOINT_COUNT = 4
V2_SOURCE_IMPLEMENTATIONS = {
    V2_RUNNER: RECORD_HASHES[V2_RUNNER],
    V2_SOURCE_CODE: RECORD_HASHES[V2_SOURCE_CODE],
}
V2_FALSE_RU_KEYS = {
    'bridge.body_sha256.method.CurrentOriginalPrescribedBridgeC6Functions.__init__',
    'MOD.finite_N_memory_rows',
    'R_U.full_unexecuted_suffix_guards',
}
BRIDGE_BODY_SHA256 = 'df1b4b43d63a07bc8e44aa8c5db8bed8c09d73ac31626b3befbd69e750f93e5b'
BRIDGE_STALE_BODY_SHA256 = 'df1b4b43d63a07bc8e44fa945662ea1535de41a0d8bc2a53a2aacf2954d5b6040'
GAMMA_SIGMA_BODY_SHA256 = '90710db1c1ee4bdec06cb4515c082198c1f3fc006d7f1ec58daebd0cddc2476a'
C6_SIGMA_BODY_SHA256 = '1f72db601deda38dc3a82db0b52df5d6ee8c7ecad4d6bb983b718c2906217ea4'
SIGMA_DEFINITION = ('0 for x<=0; 1 for x>=1; '
                    'exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside')
FLAT_SIGMA_SHA256 = '32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7'

# The six entries are deliberately literal so the publication runner can
# compare both raw digests and decoded byte counts before parsing any JSON.
FAILURE_ARTIFACT_INVENTORY = {
    V2_SOURCE_CODE: {'sha256': RECORD_HASHES[V2_SOURCE_CODE], 'bytes': 154544, 'encoding': 'utf-8'},
    V2_RUNNER: {'sha256': RECORD_HASHES[V2_RUNNER], 'bytes': 16051, 'encoding': 'utf-8'},
    V2_SOURCE_PACKET: {'sha256': RECORD_HASHES[V2_SOURCE_PACKET], 'bytes': 278311, 'encoding': 'utf-8'},
    V2_ENTERED: {'sha256': RECORD_HASHES[V2_ENTERED], 'bytes': 668, 'encoding': 'utf-8'},
    V2_RESULT: {'sha256': RECORD_HASHES[V2_RESULT], 'bytes': 31804, 'encoding': 'utf-8'},
    V2_DIAGNOSTICS: {'sha256': RECORD_HASHES[V2_DIAGNOSTICS], 'bytes': 24157, 'encoding': 'utf-8'},
}

RU_ENDPOINT_PAIRS = {
    'c6_R1_and_second_lower': [[88437, 88438]],
    'c5_first_V_N': [[106932, 106933]],
    'c5_second_V_N': [[109340, 109341]],
    'c3_first_physical_U': [[148085, 148086]],
    'c3_second_physical_U': [[154247, 154248]],
    'common_Pstar': 33572,
    'common_L': 25272,
    'common_N': 18,
}
MQ_SOURCE_PAIRS = {
    'physical_m_through_Z2': [[148097, 154259], [148098, 154260], [148099, 154261]],
    'physical_U_premise': [[148085, 154247], [148086, 154248]],
    'physical_Q': [[148180, 154342], [148184, 154346]],
    'physical_Q_numerators': [[148175, 154337], [148183, 154345]],
}
SOURCE_HASHES = {
    'parent_report_sha256': RECORD_HASHES[PARENT],
    'c6_creator_sha256': RECORD_HASHES[C6_CREATOR],
    'bridge_creator_sha256': RECORD_HASHES[BRIDGE_CREATOR],
    'loop_creator_sha256': RECORD_HASHES[LOOP_CREATOR],
    'mod_creator_sha256': RECORD_HASHES[MOD_CREATOR],
    'c6_raw_archive_sha256': RECORD_HASHES[BACKGROUND],
    'loop_raw_archive_sha256': RECORD_HASHES[LOOP],
    'mod_raw_archive_sha256': RECORD_HASHES[MOD],
}
ENDPOINT_PROJECTION = {
    'first_switch_upper': 'coordinate node 73190 -> source value 1; radius node 85905',
    'second_switch_lower': 'coordinate node 75949 -> source value 0; radius node 88600',
    'radius_rule': '100*exp(9*1) = 100*exp(9*(1+0))',
    'phase_rule': 'same source fractional_part phase after equal-radius substitution',
}

CORRECTION_SPECS = {
    'm': {'rate': '1', 'first_integral': 106958, 'first_history': 106970,
          'first_decay': 106964, 'first_memory': 106940,
          'second_integral': 109367, 'second_history': 109379,
          'second_decay': 109373, 'second_memory': 109349,
          'first_modified': [107095, 107096, 107097],
          'first_base': [87591, 87592, 87593],
          'second_modified': [109504, 109505, 109506],
          'second_base': [90109, 90110, 90111],
          'first_incoming': [104503, 104504, 104505],
          'second_incoming': [107131, 107132, 107133]},
    'h': {'rate': '3/2', 'first_integral': 106993, 'first_history': 107005,
          'first_decay': 106999, 'first_memory': 106977,
          'second_integral': 109402, 'second_history': 109414,
          'second_decay': 109408, 'second_memory': 109386,
          'first_modified': [107101, 107102, 107103],
          'first_base': [87601, 87602, 87603],
          'second_modified': [109510, 109511, 109512],
          'second_base': [90119, 90120, 90121],
          'first_incoming': [104509, 104510, 104511],
          'second_incoming': [107137, 107138, 107139]},
    'k': {'rate': '3/2', 'first_integral': 107023, 'first_history': 107035,
          'first_decay': 107029, 'first_memory': 106977,
          'second_integral': 109432, 'second_history': 109444,
          'second_decay': 109438, 'second_memory': 109386,
          'first_modified': [107107, 107108, 107109],
          'first_base': [87611, 87612, 87613],
          'second_modified': [109516, 109517, 109518],
          'second_base': [90129, 90130, 90131],
          'first_incoming': [104515, 104516, 104517],
          'second_incoming': [107143, 107144, 107145]},
    'e': {'rate': '1', 'first_integral': 107053, 'first_history': 107065,
          'first_decay': 107059, 'first_memory': 106940,
          'second_integral': 109462, 'second_history': 109474,
          'second_decay': 109468, 'second_memory': 109349,
          'first_modified': [107113, 107114, 107115],
          'first_base': [87621, 87622, 87623],
          'second_modified': [109522, 109523, 109524],
          'second_base': [90139, 90140, 90141],
          'first_incoming': [104521, 104522, 104523],
          'second_incoming': [107149, 107150, 107151]},
    'p': {'rate': '0', 'first_integral': 107083, 'first_history': 107089,
          'first_decay': 107089, 'first_memory': 1,
          'second_integral': 109492, 'second_history': 109498,
          'second_decay': 109498, 'second_memory': 1,
          'first_modified': [107119, 107120, 107121],
          'first_base': [87631, 87632, 87633],
          'second_modified': [109528, 109529, 109530],
          'second_base': [90149, 90150, 90151],
          'first_incoming': [104527, 104528, 104529],
          'second_incoming': [107155, 107156, 107157]},
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return _hashlib.sha256(_json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def _sha_text(value):
    return _hashlib.sha256(value.encode('utf-8')).hexdigest() if isinstance(value, str) else None


def _row(rows, number):
    return rows[number] if (isinstance(rows, list) and isinstance(number, int)
                            and 0 <= number < len(rows) and isinstance(rows[number], dict)) else None


def _row_is(rows, number, **expected):
    item = _row(rows, number)
    return item is not None and all(item.get(key) == value for key, value in expected.items())


def _prefix(value, expected):
    return isinstance(value, list) and value[:len(expected)] == list(expected)


def _source_path(name):
    if name.startswith(('experiments/', 'src/', 'scripts/', 'tmp/')):
        return name
    return 'experiments/root_st073/' + name


def _raw_hash_ok(name, item):
    expected = RECORD_HASHES[name]
    actual = item.get('raw_sha256') if isinstance(item, _Mapping) else None
    return actual == expected


def _record_ok(name, item):
    text = item.get('source_text') if isinstance(item, _Mapping) else None
    if not (isinstance(item, _Mapping) and isinstance(text, str)
            and len(text) > 0):
        return False
    return (_raw_hash_ok(name, item)
            and item.get('raw_authenticated_before_decode') is True
            and item.get('source_path') == _source_path(name)
            and item.get('is_gzip') is name.endswith('.gz')
            and item.get('decoded_sha256') == _sha_text(text)
            and _sha_text(text) == V2_SOURCE_DECODED_RECORD_HASHES.get(name, RECORD_HASHES[name]))


def _json_object(text, label):
    value = _json.loads(text)
    if not isinstance(value, dict):
        raise ValueError(label + ': object required')
    return value


def _definition(tree, class_name, function_name):
    if class_name is None:
        return next((node for node in tree.body
                     if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))
                     and node.name == function_name), None)
    cls = next((node for node in tree.body
                if isinstance(node, _ast.ClassDef) and node.name == class_name), None)
    if cls is None:
        return None
    return next((node for node in cls.body
                 if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))
                 and node.name == function_name), None)


def _body_sha(node):
    if node is None:
        return None
    module = _ast.Module(body=list(node.body), type_ignores=[])
    dumped = _ast.dump(module, annotate_fields=True, include_attributes=False)
    return _hashlib.sha256(dumped.encode('utf-8')).hexdigest()


def _bridge_source_guard(text):
    result = {
        'source_sha256': _sha_text(text) == RECORD_HASHES[BRIDGE_CREATOR],
        'ast_parse': False,
        'body_present': False,
        'body_sha256': False,
        'body_sha256_value': None,
        'body_is_not_stale': False,
        'bridge_function_needles': False,
    }
    try:
        tree = _ast.parse(text, filename=BRIDGE_CREATOR)
        node = _definition(tree, 'CurrentOriginalPrescribedBridgeC6Functions', '__init__')
        segment = _ast.get_source_segment(text, node) if node is not None else ''
        histories = _definition(tree, None, 'histories')
        generic = _definition(tree, None, 'generic_fields')
        history_segment = _ast.get_source_segment(text, histories) if histories is not None else ''
        generic_segment = _ast.get_source_segment(text, generic) if generic is not None else ''
        actual = _body_sha(node)
        result.update(ast_parse=True, body_present=node is not None,
                      body_sha256=actual == BRIDGE_BODY_SHA256,
                      body_sha256_value=actual,
                      body_is_not_stale=actual != BRIDGE_STALE_BODY_SHA256,
                      bridge_function_needles=(
                          all(needle in history_segment for needle in (
                              'def histories(a,phi,U,initial,phase,jacobian,label):',
                              'bodies=densities(a,phi,U)',
                              'for name,rate in RATES.items():'))
                          and "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)" in generic_segment))
    except (SyntaxError, TypeError, ValueError):
        pass
    result['passed'] = all(result[key] is True for key in (
        'source_sha256', 'ast_parse', 'body_present', 'body_sha256',
        'body_is_not_stale', 'bridge_function_needles'))
    return result


def _creator_formula_guard(text, name, expected_sha, needles):
    result = {'source_sha256': _sha_text(text) == expected_sha,
              'ast_parse': False, 'body_present': False, 'body_sha256': False,
              'formula_needles': False}
    try:
        tree = _ast.parse(text, filename=name)
        node = _definition(tree, 'CurrentOriginalModifiedQ4PressureC3Functions', '__init__')
        segment = _ast.get_source_segment(text, node) if node is not None else ''
        result.update(ast_parse=True, body_present=node is not None,
                      body_sha256=_body_sha(node) ==
                      '8040fa64b5da2de7fa1e2e29cd2583bb803b5b0735d36357050200ad87d645d4',
                      formula_needles=all(needle in text or needle in segment for needle in needles))
    except (SyntaxError, TypeError, ValueError):
        pass
    result['passed'] = all(result[key] is True for key in (
        'source_sha256', 'ast_parse', 'body_present', 'body_sha256', 'formula_needles'))
    return result


def _sigma_source_guard(gamma_text, c6_text, background, parent_pins):
    predicates = {
        'gamma_source_sha256': _sha_text(gamma_text) == RECORD_HASHES[GAMMA_CREATOR],
        'c6_source_sha256': _sha_text(c6_text) == RECORD_HASHES[C6_CREATOR],
        'gamma_pin': isinstance(parent_pins, _Mapping)
            and parent_pins.get(GAMMA_CREATOR) == RECORD_HASHES[GAMMA_CREATOR],
        'flat_sigma_pin': isinstance(parent_pins, _Mapping)
            and parent_pins.get('lei_ren_part1_paper_compliant_flat_pulse_derivatives.py') == FLAT_SIGMA_SHA256,
        'definition': SIGMA_DEFINITION ==
            '0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside',
        'fixed_cutoff_rows': False,
        'gamma_body_sha256': False,
        'c6_delegate_body_sha256': False,
    }
    try:
        gamma_tree, c6_tree = _ast.parse(gamma_text), _ast.parse(c6_text)
        predicates['gamma_body_sha256'] = _body_sha(_definition(gamma_tree, None, 'sigma')) == GAMMA_SIGMA_BODY_SHA256
        predicates['c6_delegate_body_sha256'] = _body_sha(_definition(c6_tree, None, 'sigma')) == C6_SIGMA_BODY_SHA256
        rows = background['exact_graph_nodes']
        binding = {
            'AST_sha256': '4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c',
            'function': 'stable_sigma', 'module': 'lei_ren_part1_paper_compliant_outer_initial.py',
            'module_sha256': 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657',
        }
        predicates['fixed_cutoff_rows'] = all(_row_is(
            rows, number, operation='exact_original_source_flat_sigma', argument=argument,
            definition=SIGMA_DEFINITION, source_binding=binding,
            exact_function_not_its_interval_cap=True)
            for number, argument in ((74428, 73190), (76947, 75949)))
    except (KeyError, IndexError, TypeError, SyntaxError, ValueError):
        pass
    passed = all(value is True for value in predicates.values())
    return {
        'schema': 'native_J3_v3_exact_cutoff_endpoint_source_laws_v1',
        'passed': passed, 'source_only': True, 'predicates': predicates,
        'source_sha256': RECORD_HASHES[GAMMA_CREATOR],
        'cutoff_definition': SIGMA_DEFINITION, 'fixed_cutoff_rows': [74428, 76947],
        'endpoint_source_values': {'at_0': 0, 'at_1': 1} if passed else None,
        'ordinary_Z_orders': [0, 1], 'endpoint_inputs_Z_independent': True,
        'source_body_hashes': {
            'Gamma.sigma': GAMMA_SIGMA_BODY_SHA256,
            'C6.sigma': C6_SIGMA_BODY_SHA256,
        },
        'derivation': [
            'The exact piecewise cutoff gives sigma(1)=1 and sigma(0)=0.',
            'The two source arguments are endpoint constants, so their ordinary Z derivatives are zero.',
        ],
        'numeric_values_computed': False, 'source_constructor_calls': 0,
        'old_proof_functions_called': 0, 'generic_DAG_walks': 0,
    }


def _saved_chain(records):
    """Authenticate the six saved V2 bytes and return their JSON facts."""
    source_text = records[V2_SOURCE_PACKET]['source_text']
    entered_text = records[V2_ENTERED]['source_text']
    result_text = records[V2_RESULT]['source_text']
    diagnostics_text = records[V2_DIAGNOSTICS]['source_text']
    source = _json_object(source_text, 'saved source')
    entered = _json_object(entered_text, 'saved entered')
    result = _json_object(result_text, 'saved result')
    lines = [_json.loads(line) for line in diagnostics_text.splitlines() if line.strip()]
    stages = result.get('stage_diagnostics')
    stage_map = {event.get('stage'): event for event in stages} if isinstance(stages, list) else {}
    prefix = stage_map.get('boolean_prefix', {}).get('payload', {})
    cutoff = stage_map.get('cutoff_source_laws', {}).get('payload', {})
    ru = stage_map.get('R_U', {}).get('payload', {})
    ru_preds = ru.get('predicates') if isinstance(ru, _Mapping) else None
    false_keys = {name for name, value in ru_preds.items()
                  if value is False} if isinstance(ru_preds, _Mapping) else set()
    prefix_preds = prefix.get('predicates') if isinstance(prefix, _Mapping) else None
    cutoff_preds = cutoff.get('predicates') if isinstance(cutoff, _Mapping) else None
    inventory_ok = all(
        len(records[path]['source_text'].encode('utf-8')) == meta['bytes']
        and records[path].get('raw_sha256') == meta['sha256']
        for path, meta in FAILURE_ARTIFACT_INVENTORY.items())
    predicates = {
        'inventory_exact': inventory_ok,
        'source_schema': source.get('schema') == 'native_J3_endpoint_transport_inputs_v2',
        'source_input_digest': source.get('input_digest') == V2_SOURCE_DIGEST,
        'source_record_hashes': source.get('record_hashes') == V2_SOURCE_RECORD_HASHES,
        'source_decoded_record_hashes': source.get('decoded_record_hashes') == V2_SOURCE_DECODED_RECORD_HASHES,
        'source_implementations': source.get('implementation_hashes') == V2_SOURCE_IMPLEMENTATIONS,
        'entered_component': entered.get('component') == 'endpoint_transport',
        'entered_source_chain': entered.get('source_sha256') == RECORD_HASHES[V2_SOURCE_PACKET],
        'entered_input_chain': entered.get('input_digest') == V2_SOURCE_DIGEST,
        'entered_kernel_chain': entered.get('kernel_sha256') == RECORD_HASHES[V2_SOURCE_CODE],
        'result_failure': result.get('passed') is False and result.get('failed_attempt_preserved') is True,
        'result_error': result.get('error_type') == 'AssertionError'
            and result.get('error') == 'Source prerequisite: new_R_U_C1',
        'result_checkpoints': result.get('source_checkpoint_sha256') == RECORD_HASHES[V2_SOURCE_PACKET]
            and result.get('entry_checkpoint_sha256') == RECORD_HASHES[V2_ENTERED]
            and result.get('diagnostic_checkpoint_sha256') == RECORD_HASHES[V2_DIAGNOSTICS],
        'result_event_count': result.get('diagnostic_events_saved_before_validation') == 3,
        'diagnostic_chain': isinstance(stages, list) and len(stages) == 3 and lines == stages
            and [event.get('stage') for event in stages] == ['boolean_prefix', 'cutoff_source_laws', 'R_U'],
        'prefix_positive': isinstance(prefix_preds, _Mapping) and bool(prefix_preds)
            and all(type(value) is bool and value is True for value in prefix_preds.values()),
        'cutoff_positive': isinstance(cutoff_preds, _Mapping) and bool(cutoff_preds)
            and all(type(value) is bool and value is True for value in cutoff_preds.values()),
        'ru_false_set': false_keys == V2_FALSE_RU_KEYS,
        'ru_failure_preserved': isinstance(ru_preds, _Mapping)
            and ru.get('J3_R_C1_proved') is False and ru.get('J3_U_C1_proved') is False,
    }
    saved_ok = all(value is True for value in predicates.values())
    return {
        'passed': saved_ok, 'predicates': predicates, 'source': source,
        'entered': entered, 'result': result, 'stages': stages,
        'prefix': prefix, 'cutoff': cutoff, 'ru': ru,
        'ru_predicates': dict(ru_preds) if isinstance(ru_preds, _Mapping) else {},
        'saved_failure_flags': {key: False for key in sorted(false_keys)},
        'saved_positive_predicates': {
            key: value for key, value in (ru_preds.items() if isinstance(ru_preds, _Mapping) else ())
            if value is True},
    }


def _inventory_shape(report):
    inventory = report.get('artifact_inventory') if isinstance(report, _Mapping) else None
    return _pin_map_ok(inventory)


def _pin_map_ok(value):
    """Require the exact six-key saved-attempt inventory shape."""
    if isinstance(value, _Mapping):
        value = value.get('artifacts', value.get('items', value))
    if isinstance(value, _Mapping):
        value = [dict(path=path, **(entry if isinstance(entry, _Mapping) else {}))
                 for path, entry in value.items()]
    if not isinstance(value, list):
        return False
    by_path = {item.get('path'): item for item in value if isinstance(item, _Mapping)}
    return (set(by_path) == set(FAILURE_ARTIFACT_INVENTORY)
            and all(by_path[path].get('sha256') == item['sha256']
                    and by_path[path].get('bytes') == item['bytes']
                    for path, item in FAILURE_ARTIFACT_INVENTORY.items()))


def _frozen_inventory_ok(report):
    frozen = report.get('frozen_failure_inventory') if isinstance(report, _Mapping) else None
    artifacts = frozen.get('artifacts') if isinstance(frozen, _Mapping) else None
    expected = {path: {'sha256': item['sha256'], 'bytes': item['bytes']}
                for path, item in FAILURE_ARTIFACT_INVENTORY.items()}
    return (
        isinstance(frozen, _Mapping)
        and frozen.get('path') == V2_FROZEN_INVENTORY_PATH
        and frozen.get('sha256') == V2_FROZEN_INVENTORY_SHA256
        and frozen.get('bytes') == V2_FROZEN_INVENTORY_BYTES
        and frozen.get('schema') == V2_FROZEN_INVENTORY_SCHEMA
        and frozen.get('attempt_count') == 1
        and frozen.get('do_not_reenter') is True
        and frozen.get('result_passed') is False
        and frozen.get('accepted_J3_identity') is False
        and frozen.get('remaining_B_nodes') == list(range(523, 536))
        and frozen.get('ancestor_operators_replayed') is False
        and frozen.get('all_old_errors_retained') is True
        and isinstance(artifacts, _Mapping)
        and set(artifacts) == set(expected)
        and all(artifacts.get(path) == item for path, item in expected.items()))


def _accepted_parent_ok(parent):
    old = parent.get('component_result', {}) if isinstance(parent, _Mapping) else {}
    retained = parent.get('retained_construction_data', {}) if isinstance(parent, _Mapping) else {}
    scope = old.get('B_identity_scope', {}) if isinstance(old, _Mapping) else {}
    pins = parent.get('source_pins', {}) if isinstance(parent, _Mapping) else {}
    return {
        'schema': parent.get('source_only') is True and parent.get('full_goal_complete') is False,
        'component': old.get('passed') is True and old.get('J2_B_identity_proved') is True,
        'pins': isinstance(pins, _Mapping) and len(pins) == 1850,
        'ordinary_Z': old.get('ordinary_Z_orders') == [0, 1]
            and old.get('ordinary_Z_not_temporal_n') is True,
        'B_scope': scope.get('cumulative_B_identity_nodes') == [520, 521, 522]
            and scope.get('unresolved_B_identity_nodes') == list(range(523, 536)),
        'IN_scope': old.get('native_IN_reduced_source_identity') ==
            'I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)',
        'ledger': all(retained.get(key) == value for key, value in {
            'cumulative_terms': 212, 'raw_signed_terms': 213, 'reduced_terms': 183,
            'six_signed_channels_retained': True, 'all_old_error_data_retained': True,
            'reserved_contribution_registered_once': True, 'inherited_errors_deleted': False,
        }.items()),
        'source_family': old.get('source_family') == FAMILY,
    }


def _failed_parent_ok(report, saved_result):
    component = report.get('component_result', {}) if isinstance(report, _Mapping) else {}
    counters = {
        'source_pin_count': report.get('source_pin_count') == INHERITED_SOURCE_COUNT,
        'checkpoint_count': report.get('checkpoint_count') == NEW_CHECKPOINT_COUNT,
        'total_checkpoint_count': report.get('total_checkpoint_count') == INHERITED_CHECKPOINT_COUNT,
        'failed_union_count': report.get('preserved_failed_artifact_count_union') == FAILED_ARTIFACT_UNION_COUNT,
        'artifact_count': report.get('preserved_failed_artifact_count') == 6,
        'consumer_count': report.get('source_consumer_attempt_count') == 1,
        'ancestor_calls': report.get('accepted_ancestor_operator_calls') == 0,
        'old_calls': report.get('old_proof_functions_called') == 0,
    }
    return {
        'schema': report.get('schema') == 'original_native_J3_endpoint_transport_v2_failure_stage_v1',
        'scope': report.get('component') == 'endpoint_transport'
            and report.get('source_only') is True
            and report.get('accepted_identity') is False
            and report.get('full_goal_complete') is False,
        'stage': report.get('failure_stage') == 'new_R_U_C1',
        'counters': all(counters.values()),
        'component_failure': component.get('passed') is False
            and component.get('failed_attempt_preserved') is True
            and component.get('error') == 'Source prerequisite: new_R_U_C1'
            and component.get('diagnostic_events_saved_before_validation') == 3,
        'inventory': _inventory_shape(report),
        'artifact_evidence': (report.get('frozen_artifacts') == report.get('artifact_inventory')
            and _inventory_shape(report)),
        'inventory_pointer': isinstance(report.get('artifact_inventory'), list)
            and len(report['artifact_inventory']) == 6,
        'frozen_inventory_pointer': _frozen_inventory_ok(report),
        'saved_result_chain': component == saved_result,
        'false_flags': report.get('false_flags') == dict.fromkeys(FALSE_FLAGS, False),
    }


def _failure_evidence_ok(evidence, saved_result):
    full = evidence.get('full_result', {}) if isinstance(evidence, _Mapping) else {}
    freeze = evidence.get('freeze_metadata', {}) if isinstance(evidence, _Mapping) else {}
    items = evidence.get('artifact_inventory') if isinstance(evidence, _Mapping) else None
    return {
        'schema': evidence.get('schema') == 'original_native_J3_endpoint_transport_failure_evidence_v1',
        'scope': evidence.get('component') == 'endpoint_transport'
            and evidence.get('metadata_only') is True
            and evidence.get('failed_attempt_preserved') is True
            and evidence.get('raw_bytes_preserved') is True
            and evidence.get('runtime_objects_recovered') is False
            and evidence.get('numeric_values_recovered') is False
            and evidence.get('accepted_identity') is False
            and evidence.get('source_prefix_objects_recovered') is False
            and evidence.get('do_not_reenter') is True
            and evidence.get('error_stage') == 'new_R_U_C1'
            and evidence.get('attempt_count') == 1,
        'full_result': full.get('passed') is False
            and full.get('failed_attempt_preserved') is True
            and full.get('error') == 'Source prerequisite: new_R_U_C1',
        'freeze': freeze.get('accepted_ancestor_operator_calls') == 0
            and freeze.get('accepted_source_consumer_calls') == 0
            and freeze.get('do_not_reenter') is True
            and freeze.get('full_result_saved') is True
            and freeze.get('numeric_evaluations') == 0
            and freeze.get('old_operators_replayed') is False
            and freeze.get('old_proof_functions_called') == 0
            and freeze.get('source_consumer_attempts') == 1
            and freeze.get('source_execution_occurred') is True
            and freeze.get('runner_caught_exception') is True,
        'saved_result_error': saved_result.get('error') == full.get('error')
            and saved_result.get('passed') is False,
        'inventory': isinstance(items, list) and len(items) == 5
            and all(isinstance(item, _Mapping) and item.get('sha256')
                    and item.get('bytes') and item.get('encoding') == 'base64'
                    for item in items),
    }


def _archive_shape(archive, wrapper, count):
    return (isinstance(archive, _Mapping)
            and isinstance(archive.get('exact_graph_nodes'), list)
            and len(archive['exact_graph_nodes']) == count
            and isinstance(archive.get(wrapper), _Mapping)
            and archive.get('source_family') == FAMILY)


def _finite_rows(loop_rows):
    channels = {}
    endpoint_aliases = {}
    for field, spec in CORRECTION_SPECS.items():
        ok = True
        endpoint_aliases[field] = list(spec['second_incoming'])
        for order in range(3):
            first = _row(loop_rows, spec['first_integral'] + order)
            second = _row(loop_rows, spec['second_integral'] + order)
            base_integral = {
                'operation': 'definite_integral', 'lower': 0,
                'exact_function_integral': True,
                'numerical_value_not_installed': True,
                'endpoint_independent_of_Z': True,
                'true_core_inlet_complete_history': True,
                'numerical_cap_not_defining_integrand': True,
                'actual_finite_N_signed_correction': True,
                'physical_Jacobian_applied_once': True,
                'zero_correction_inlet_only_at_original_rminus': True,
                'predecessor_memory_never_reset': True,
            }
            ok = ok and isinstance(first, dict) and isinstance(second, dict)
            if isinstance(first, dict) and isinstance(second, dict):
                ok = ok and all(first.get(key) == value for key, value in base_integral.items())
                ok = ok and all(second.get(key) == value for key, value in base_integral.items())
                ok = ok and first.get('variable') == 'true_loop_C5_transport_first_switch'
                ok = ok and second.get('variable') == 'true_loop_C5_transport_second_switch'
                ok = ok and first.get('upper') == 73190 and second.get('upper') == 75949
                ok = ok and first.get('ordinary_slow_Z_derivative_order') == order
                ok = ok and second.get('ordinary_slow_Z_derivative_order') == order
                ok = ok and first.get('original_rate') == spec['rate']
                ok = ok and second.get('original_rate') == spec['rate']
                ok = ok and first.get('original_chart') == 'first_switch'
                ok = ok and second.get('original_chart') == 'second_switch'
            first_head = spec['first_decay'] + order if spec['rate'] != '0' else spec['first_incoming'][order]
            second_head = spec['second_decay'] + order if spec['rate'] != '0' else spec['second_incoming'][order]
            ok = ok and _row_is(loop_rows, spec['first_history'] + order,
                                operation='sum', arguments=[first_head, spec['first_integral'] + order])
            ok = ok and _row_is(loop_rows, spec['second_history'] + order,
                                operation='sum', arguments=[second_head, spec['second_integral'] + order])
            if spec['rate'] != '0':
                ok = ok and _row_is(loop_rows, spec['first_decay'] + order,
                                    operation='product', arguments=[spec['first_incoming'][order], spec['first_memory']])
                ok = ok and _row_is(loop_rows, spec['second_decay'] + order,
                                    operation='product', arguments=[spec['second_incoming'][order], spec['second_memory']])
            else:
                ok = ok and spec['first_decay'] == spec['first_history']
                ok = ok and spec['second_decay'] == spec['second_history']
        for order, node in enumerate(endpoint_aliases[field]):
            ok = ok and _row_is(loop_rows, node, operation='function_substitution',
                                expression=spec['first_history'] + order,
                                variable=73190, value=1,
                                Z_independent_substitution=True)
        channels[field] = ok
    return channels


def _mod_rows(mod):
    try:
        binding = mod['modified_Q4_pressure_C3_bindings']
        rows = mod['exact_graph_nodes']
        first = binding['charts']['first_switch']['physical']
        second = binding['charts']['second_switch']['physical']
        result = {
            'typed_physical_ports': (_prefix(first.get('U'), [148085, 148086])
                and _prefix(first.get('m'), [148097, 148098, 148099])
                and _prefix(first.get('Q'), [148180, 148184])
                and _prefix(second.get('U'), [154247, 154248])
                and _prefix(second.get('m'), [154259, 154260, 154261])
                and _prefix(second.get('Q'), [154342, 154346])),
            'physical_U_Pstar': all(_row_is(rows, node, operation='product', arguments=[source, 33572])
                for node, source in ((148085, 106932), (148086, 106933),
                                     (154247, 109340), (154248, 109341))),
            'finite_N_m_memory': (all(_row_is(rows, node, operation='sum', arguments=[base, history])
                for node, base, history in zip((107095, 107096, 107097),
                                                (87591, 87592, 87593),
                                                (106970, 106971, 106972)))
                and all(_row_is(rows, node, operation='sum', arguments=[base, history])
                for node, base, history in zip((109504, 109505, 109506),
                                                (90109, 90110, 90111),
                                                (109379, 109380, 109381)))),
            'physical_m_Pstar_once': (all(_row_is(rows, node, operation='product', arguments=[source, 33572])
                for node, source in zip((148097, 148098, 148099), (107095, 107096, 107097)))
                and all(_row_is(rows, node, operation='product', arguments=[source, 33572])
                for node, source in zip((154259, 154260, 154261), (109504, 109505, 109506)))),
            'Q_numerators': all(_row_is(rows, node, operation='sum', arguments=args)
                for node, args in ((148175, [148127, 148170]), (148183, [148176, 148182]),
                                   (154337, [154289, 154332]), (154345, [154338, 154344]))),
            'Q_positive_L': all(_row_is(rows, node, operation='positive_quotient', numerator=numerator,
                denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
                for node, numerator in ((148180, 148175), (148184, 148183),
                                        (154342, 154337), (154346, 154345))),
            'common_L': _row_is(rows, 25272, operation='sum', arguments=[1, 25271]),
        }
        result['all'] = all(result.values())
        return result
    except (KeyError, IndexError, TypeError, AttributeError):
        return {name: False for name in ('typed_physical_ports', 'physical_U_Pstar',
            'finite_N_m_memory', 'physical_m_Pstar_once', 'Q_numerators', 'Q_positive_L',
            'common_L', 'all')}


def _endpoint_rows(background, loop, mod):
    try:
        cb = background['prescribed_R100_Rm_C6_bindings']
        lb = loop['prescribed_loop_C5_bindings']
        mb = mod['modified_Q4_pressure_C3_bindings']
        rows = background['exact_graph_nodes']
        loop_rows = loop['exact_graph_nodes']
        first, second = cb['charts']['first_switch'], cb['charts']['second_switch']
        r1 = cb['endpoints']['R1']
        first_l, second_l = lb['charts']['first_switch'], lb['charts']['second_switch']
        first_m, second_m = mb['charts']['first_switch'], mb['charts']['second_switch']
        out = {
            'chart_types': (first.get('phase') == 73190 and second.get('phase') == 75949
                and first.get('domain') == second.get('domain') == [0, 1]
                and first.get('radius') == 85905 and second.get('radius') == 88600
                and first_l.get('coordinate') == 73190 and second_l.get('coordinate') == 75949
                and first_l.get('radius') == 85905 and second_l.get('radius') == 88600),
            'endpoint_actual_rows': (_prefix(first['actual']['phi'], [86947, 86950])
                and _prefix(first['actual']['U'], [87154, 87155])
                and _prefix(first['actual']['phi_y'], [87196, 87197])
                and _prefix(first['actual']['U_y'], [87126, 87127])
                and second['actual']['U_y'] == [0, 0, 0, 0, 0, 0, 0]
                and second['radial']['V_y'] == [0, 0, 0, 0, 0, 0, 0]),
            'second_actual_transport_rows': (_prefix(second['actual']['phi'], [89644, 89647])
                and _prefix(second['actual']['phi_y'], [89687, 89690])),
            'R1_alias_rows': (_prefix(r1['actual']['U'], [88437, 88438])
                and _prefix(second['incoming_histories']['M'], [88465, 88466, 88467])
                and _prefix(r1['actual']['phi'], [88430, 88431])
                and _prefix(r1['actual']['phi_y'], [88444, 88445])
                and _prefix(r1['actual']['U_y'], [88451, 88452])
                and _prefix(r1['fields']['E'], [88500, 88501])
                and _prefix(r1['fields']['V'], [88507, 88508])
                and _prefix(r1['generic']['p2'], [88592, 88593])
                and _prefix(second['actual']['U'], [88437, 88438])),
            'bridge_field_rows': (_prefix(first['fields']['E'], [87550, 87553])
                and _prefix(first['fields']['V'], [87584, 87585])
                and _prefix(first['radial']['E_y'], [87790, 87793])
                and _prefix(first['radial']['V_y'], [87824, 87825])
                and _prefix(first['generic']['p2'], [88423, 88424])
                and _prefix(second['fields']['E'], [90068, 90071])
                and _prefix(second['fields']['V'], [90102, 90103])
                and _prefix(second['radial']['E_y'], [90308, 90311])
                and _prefix(second['generic']['p2'], [90934, 90935])),
            'c5_root_operand_rows': (_prefix(first_l['roots']['a'], [104562, 104566])
                and _prefix(first_l['roots']['b'], [104597, 104601])
                and _prefix(first_l['roots']['p1'], [104650, 104654])
                and _prefix(second_l['roots']['a'], [107184, 107188])
                and _prefix(second_l['roots']['b'], [0, 0])
                and _prefix(second_l['roots']['p1'], [107237, 107241])),
            'loop_root_rows': (_prefix(first_l['roots']['E'], [87550, 87553])
                and _prefix(first_l['roots']['V'], [87584, 87585])
                and _prefix(second_l['roots']['E'], [90068, 90071])
                and _prefix(second_l['roots']['V'], [90102, 90103])
                and _prefix(first_l['roots']['p2'], [88423, 88424])
                and _prefix(second_l['roots']['p2'], [90934, 90935])),
            'p2_cross_bindings': (first_l['roots']['p2'] == first['generic']['p2']
                and second_l['roots']['p2'] == second['generic']['p2']
                and _row_is(loop_rows, 88423, operation='product', arguments=[88388, 85905])
                and _row_is(loop_rows, 88424, operation='product', arguments=[88392, 85905])
                and _row_is(loop_rows, 90934, operation='product', arguments=[90899, 88600])
                and _row_is(loop_rows, 90935, operation='product', arguments=[90903, 88600])),
            'physical_ports': (_prefix(first_m['physical']['U'], [148085, 148086])
                and _prefix(second_m['physical']['U'], [154247, 154248])
                and _prefix(first_m['physical']['m'], [148097, 148098, 148099])
                and _prefix(second_m['physical']['m'], [154259, 154260, 154261])
                and _prefix(first_m['physical']['Q'], [148180, 148184])
                and _prefix(second_m['physical']['Q'], [154342, 154346])),
            'radius_comparison': (first.get('radius') == 85905 and second.get('radius') == 88600
                and first.get('known_comparison_direction', {}).get('D_over_R', [])[:2] == [86464, 86468]
                and second.get('known_comparison_direction', {}).get('D_over_R', [])[:2] == [89159, 89163]),
            'comparison_direction_rows': (_prefix(
                first.get('known_comparison_direction', {}).get('drive_hydro', []), [86650, 86654])
                and _prefix(first.get('known_comparison_direction', {}).get('drive_pressure', []), [68391, 68395])
                and _prefix(first.get('known_comparison_direction', {}).get('drive_swirl', []), [86853, 86857])
                and _prefix(second.get('known_comparison_direction', {}).get('drive_hydro', []), [89345, 89349])
                and _prefix(second.get('known_comparison_direction', {}).get('drive_pressure', []), [68391, 68395])
                and _prefix(second.get('known_comparison_direction', {}).get('drive_swirl', []), [89548, 89552])),
            'fixed_substitutions': all(_row_is(rows, number,
                operation='function_substitution', expression=expression,
                variable=73190, value=1, Z_independent_substitution=True)
                for number, expression in ((88430, 86947), (88431, 86950),
                (88437, 87154), (88438, 87155), (88444, 87196), (88445, 87197),
                (88451, 87126), (88452, 87127), (88500, 87550), (88501, 87553),
                (88507, 87584), (88508, 87585), (88592, 88423), (88593, 88424))),
        }
        out['all'] = all(out.values())
        return out
    except (KeyError, IndexError, TypeError, AttributeError):
        return {name: False for name in ('chart_types', 'endpoint_actual_rows', 'R1_alias_rows',
            'second_actual_transport_rows', 'bridge_field_rows', 'c5_root_operand_rows',
            'loop_root_rows', 'p2_cross_bindings', 'physical_ports', 'radius_comparison',
            'comparison_direction_rows', 'fixed_substitutions', 'all')}


def _archive_basics(background, loop, mod, raw_hashes):
    return {
        'background_shape': _archive_shape(background, 'prescribed_R100_Rm_C6_bindings', 97295),
        'loop_shape': _archive_shape(loop, 'prescribed_loop_C5_bindings', 121834),
        'mod_shape': _archive_shape(mod, 'modified_Q4_pressure_C3_bindings', 183700),
        'raw_sha256': isinstance(raw_hashes, _Mapping)
            and raw_hashes.get('background') == RECORD_HASHES[BACKGROUND]
            and raw_hashes.get('loop') == RECORD_HASHES[LOOP]
            and raw_hashes.get('mod') == RECORD_HASHES[MOD],
        'Pstar': (background.get('prescribed_R100_Rm_C6_bindings', {}).get('Pstar') == 33572
            and _row_is(mod.get('exact_graph_nodes'), 33572,
                        operation='analytic_unary', name='exp', argument=5)),
        'common_L': _row_is(mod.get('exact_graph_nodes'), 25272, operation='sum', arguments=[1, 25271]),
        'source_family': all(value.get('source_family') == FAMILY
            for value in (background, loop, mod) if isinstance(value, _Mapping)),
    }


def _fresh_ru(saved, cutoff, background, loop, mod, c6_text, bridge_text,
              loop_text, mod_text, raw_hashes, parent_pins):
    saved_preds = saved.get('ru_predicates', {})
    saved_positive = bool(saved.get('passed')) and all(
        type(value) is bool and value is True
        for key, value in saved_preds.items() if key not in V2_FALSE_RU_KEYS)
    bridge = _bridge_source_guard(bridge_text)
    mod_creator = _creator_formula_guard(mod_text, MOD_CREATOR, RECORD_HASHES[MOD_CREATOR], (
        "U=a.scale(fields['V'],Pstar)", "m=a.scale(fields['m'],Pstar)",
        'L=b.sub(one,b.mul(de,b.mul(z,z)))', 'Q=div(numerator(U,m))',
        'Qy=div(numerator(Uy,my))'))
    finite = _finite_rows(loop['exact_graph_nodes'])
    finite_all = all(finite.values())
    mods = _mod_rows(mod)
    endpoint = _endpoint_rows(background, loop, mod)
    basics = _archive_basics(background, loop, mod, raw_hashes)
    creator_sources = {
        'c6_sha256': _sha_text(c6_text) == RECORD_HASHES[C6_CREATOR],
        'loop_sha256': _sha_text(loop_text) == RECORD_HASHES[LOOP_CREATOR],
        'mod_sha256': mod_creator['source_sha256'],
        'bridge_sha256': bridge['source_sha256'],
    }
    predicates = {
        'saved_prefix_positive': saved_positive and saved.get('predicates', {}).get('diagnostic_chain') is True,
        'saved_cutoff_positive': saved_positive and cutoff.get('passed') is True,
        'saved_R_U_positive_predicates': saved_positive,
        'saved_failure_flags_preserved': saved.get('saved_failure_flags') ==
            {key: False for key in sorted(V2_FALSE_RU_KEYS)},
        'archive_basics': all(basics.values()),
        'creator_sources': all(creator_sources.values()),
        'bridge.body_sha256.method.CurrentOriginalPrescribedBridgeC6Functions.__init__': bridge['body_sha256'],
        'bridge.body_sha256_not_stale': bridge['body_is_not_stale'],
        'bridge.histories_and_generic_rows': bridge['bridge_function_needles'],
        'MOD.creator_actual_physical_formula': mod_creator['passed'],
        'LOOP.finite_N_h_k_e_p_m_rows_Z2': finite_all,
        'MOD.finite_N_m_memory_rows': mods['finite_N_m_memory'],
        'MOD.physical_U_Pstar_rows': mods['physical_U_Pstar'],
        'MOD.physical_m_Pstar_once_rows': mods['physical_m_Pstar_once'],
        'MOD.Q_numerator_rows': mods['Q_numerators'],
        'MOD.Q_positive_L_rows': mods['Q_positive_L'],
        'MOD.common_L_recipe': mods['common_L'],
        'source_join.endpoint_semantic_law': endpoint['endpoint_actual_rows'],
        'source_join.second_actual_transport_rows': endpoint['second_actual_transport_rows'],
        'source_join.R1_alias_rows': endpoint['R1_alias_rows'],
        'source_join.bridge_field_rows': endpoint['bridge_field_rows'],
        'source_join.c5_root_operand_rows': endpoint['c5_root_operand_rows'],
        'source_join.loop_root_rows': endpoint['loop_root_rows'],
        'source_join.p2_cross_bindings': endpoint['p2_cross_bindings'],
        'source_join.physical_ports': endpoint['physical_ports'],
        'source_join.common_radius_comparison': endpoint['radius_comparison'],
        'source_join.comparison_direction_rows': endpoint['comparison_direction_rows'],
        'source_join.fixed_substitutions': endpoint['fixed_substitutions'],
        'source_semantic.sigma_cutoff_law': cutoff.get('passed') is True,
        'R_U.corrected_memory_lookup_key': (mods.get('finite_N_m_memory') is True
            and saved_preds.get('MOD.finite_N_m_memory_rows') is True),
    }
    predicates['R_U.full_unexecuted_suffix_guards'] = all(predicates.values())
    passed = all(type(value) is bool and value is True for value in predicates.values())
    return {
        'schema': 'native_J3_R_U_source_transport_v3',
        'status': 'conditional_source_theorem_proved' if passed else 'source_guard_failure',
        'source_only': True, 'accepted_identity_claim': False,
        'J3_R_C1_proved': passed, 'J3_U_C1_proved': passed,
        'ordinary_Z_orders': [0, 1], 'ordinary_Z_not_temporal_n': True,
        'owned_Z_cell': ['0', '1/2'], 'native_six_signed_channels': ['H', 'M', 'K', 'A', 'B', 'C'],
        'generic_five_channels': ['m', 'h', 'k', 'e', 'p'],
        'endpoint_projection': dict(ENDPOINT_PROJECTION),
        'source_semantic_endpoint_law': {
            'phi_U_phi_y_U_y': {
                'first_switch': {'phi': [86947, 86950], 'U': [87154, 87155],
                                 'phi_y': [87196, 87197], 'U_y': [87126, 87127]},
                'R1': {'phi': [88430, 88431], 'U': [88437, 88438],
                       'phi_y': [88444, 88445], 'U_y': [88451, 88452]},
                'second_switch': {'phi': [89644, 89647], 'phi_y': [89687, 89690],
                                  'U_y': [0, 0, 0, 0, 0, 0, 0]},
            },
            'bridge_E_V_E_y_V_y': {'first': {'E': [87550, 87553], 'V': [87584, 87585],
                                              'E_y': [87790, 87793], 'V_y': [87824, 87825]},
                                   'second': {'E': [90068, 90071], 'V': [90102, 90103],
                                                'E_y': [90308, 90311], 'V_y': [0, 0]}},
            'C5_p1_p2': {
                'first': {'a': [104562, 104566], 'b': [104597, 104601],
                          'p1': [104650, 104654], 'p2': [88423, 88424]},
                'second': {'a': [107184, 107188], 'b': [0, 0],
                           'p1': [107237, 107241], 'p2': [90934, 90935]},
            },
            'p2_cross_bindings': {
                'first': {'c6_generic_p2': [88423, 88424], 'c5_root_p2': [88423, 88424],
                          'denominator_rows': [[88423, [88388, 85905]],
                                               [88424, [88392, 85905]]]},
                'second': {'c6_generic_p2': [90934, 90935], 'c5_root_p2': [90934, 90935],
                           'denominator_rows': [[90934, [90899, 88600]],
                                                [90935, [90903, 88600]]]},
            },
            'common_comparison_radius': {'first_radius': 85905, 'second_radius': 88600,
                                         'first_comparison': [86464, 86468],
                                         'second_comparison': [89159, 89163],
                                         'first_directions': {
                                             'D_over_R': [86464, 86468],
                                             'drive_hydro': [86650, 86654],
                                             'drive_pressure': [68391, 68395],
                                             'drive_swirl': [86853, 86857]},
                                         'second_directions': {
                                             'D_over_R': [89159, 89163],
                                             'drive_hydro': [89345, 89349],
                                             'drive_pressure': [68391, 68395],
                                             'drive_swirl': [89548, 89552]}},
            'symbolic_only': True,
        },
        'exact_U_row_pairs': dict(RU_ENDPOINT_PAIRS),
        'source_hashes': dict(SOURCE_HASHES), 'numeric_values_computed': False,
        'source_constructor_calls': 0, 'source_producer_function_calls': 0,
        'old_proof_functions_called': 0, 'generic_DAG_walks': 0,
        'fragment_was_unexecuted_before_frozen_entry': True,
        'predicates': predicates, 'saved_failure_flags': dict(saved.get('saved_failure_flags', {})),
        'saved_R_U_predicates': dict(saved_preds),
        'source_derivation': {
            'law': 'Use saved positive source predicates, correct the MOD finite-N producer key, and bind exact C6 endpoint, comparison, bridge, C5 root, and p2 denominator rows through authenticated source formulas.',
            'steps': [
                'Authenticate the saved Boolean prefix and exact cutoff source laws.',
                'Recheck five finite-N channels through ordinary Z orders 0, 1, and 2, including the zero-rate p branch.',
                'Recheck the six modified-m memory sums before the Pstar product rows and Q numerator/positive-L rows.',
                'Apply sigma(1)=1 and sigma(0)=0 to the C6 endpoint aliases, then bind second phi/phi_y, full comparison directions, bridge fields, C5 operands/p1, and C6/C5 p2 denominator rows.',
            ],
            'dependencies': {key: predicates.get(key) for key in (
                'saved_R_U_positive_predicates',
                'bridge.body_sha256.method.CurrentOriginalPrescribedBridgeC6Functions.__init__',
                'MOD.finite_N_m_memory_rows', 'LOOP.finite_N_h_k_e_p_m_rows_Z2',
                'source_join.endpoint_semantic_law', 'source_join.second_actual_transport_rows',
                'source_join.R1_alias_rows', 'source_join.bridge_field_rows',
                'source_join.c5_root_operand_rows', 'source_join.loop_root_rows',
                'source_join.p2_cross_bindings', 'source_join.physical_ports',
                'source_join.common_radius_comparison', 'source_join.comparison_direction_rows',
                'source_join.fixed_substitutions')},
            'corrected_lookup': 'mods[finite_N_m_memory] is bound to saved MOD.finite_N_m_memory_rows and feeds R_U.full_unexecuted_suffix_guards',
        },
        'nonclaims': {'numeric_containment': False, 'physical_m_Q_C1': False,
                      'accepted_identity_claim': False},
    }


def _fresh_mq(saved, ru, loop, mod, background, mod_text):
    mods = _mod_rows(mod)
    finite = all(_finite_rows(loop['exact_graph_nodes']).values())
    saved_preds = saved.get('ru_predicates', {})
    predicates = {
        'R_U.contract': ru.get('J3_R_C1_proved') is True and ru.get('J3_U_C1_proved') is True,
        'saved_C6_native_six_rows': saved_preds.get('C6.native_six_history_rows_Z2') is True,
        'saved_C6_generic_five_rows': saved_preds.get('C6.generic_h_k_e_p_rows') is True,
        'native_six_channels_distinct': ru.get('native_six_signed_channels') == ['H', 'M', 'K', 'A', 'B', 'C'],
        'generic_five_channels_distinct': ru.get('generic_five_channels') == ['m', 'h', 'k', 'e', 'p'],
        'saved_source_composition': all(value is True for key, value in saved_preds.items()
                                        if key not in V2_FALSE_RU_KEYS),
        'finite_N_signed_memory_preserved': finite and mods['finite_N_m_memory'],
        'physical_m_pairs': mods['physical_m_Pstar_once'],
        'physical_U_pairs': mods['physical_U_Pstar'],
        'physical_Q_pairs': mods['Q_positive_L'],
        'Q_numerator_rows': mods['Q_numerators'],
        'common_L_recipe': mods['common_L'],
        'actual_MOD_formula': _creator_formula_guard(mod_text, MOD_CREATOR,
            RECORD_HASHES[MOD_CREATOR], ("U=a.scale(fields['V'],Pstar)",
            "m=a.scale(fields['m'],Pstar)", 'Q=div(numerator(U,m))'))['passed'],
        'ordinary_Z_scope': ru.get('ordinary_Z_orders') == [0, 1]
            and ru.get('owned_Z_cell') == ['0', '1/2'],
    }
    predicates['m_Q.full_unexecuted_suffix_guards'] = all(predicates.values())
    passed = all(type(value) is bool and value is True for value in predicates.values())
    return {
        'schema': 'native_J3_m_Q_source_transport_v3',
        'status': 'source_m_Q_C1_proved' if passed else 'source_guard_failure',
        'source_only': True, 'accepted_identity_claim': False,
        'J3_Q_C1_proved': passed, 'finite_N_m_Z2_proved': passed,
        'physical_m_Z2_proved': passed, 'ordinary_Z_orders': [0, 1],
        'owned_Z_cell': ['0', '1/2'], 'native_six_signed_channels': ['H', 'M', 'K', 'A', 'B', 'C'],
        'generic_five_channels': ['m', 'h', 'k', 'e', 'p'],
        'source_pairs': dict(MQ_SOURCE_PAIRS), 'predicates': predicates,
        'numeric_values_computed': False, 'source_constructor_calls': 0,
        'source_producer_function_calls': 0, 'old_proof_functions_called': 0,
        'generic_DAG_walks': 0, 'fragment_was_unexecuted_before_frozen_entry': True,
        'source_derivation': {
            'law': 'Finite-N memory is preserved before one physical Pstar scaling; Q is the same positive-L quotient of the exact numerator rows.',
            'dependencies': {key: predicates.get(key) for key in (
                'R_U.contract', 'finite_N_signed_memory_preserved', 'physical_m_pairs',
                'physical_U_pairs', 'physical_Q_pairs', 'Q_numerator_rows',
                'common_L_recipe', 'actual_MOD_formula')},
            'rows': {
                'finite_N_m_memory': [[107095, 87591, 106970], [107096, 87592, 106971],
                                      [107097, 87593, 106972], [109504, 90109, 109379],
                                      [109505, 90110, 109380], [109506, 90111, 109381]],
                'physical_m': [[148097, 107095], [148098, 107096], [148099, 107097],
                               [154259, 109504], [154260, 109505], [154261, 109506]],
                'Q_numerators': [[148175, [148127, 148170]], [148183, [148176, 148182]],
                                 [154337, [154289, 154332]], [154345, [154338, 154344]]],
            },
        },
        'nonclaims': {'Q_y_C1_identity': False, 'weighted_J3_jump_zero': False,
                      'numeric_containment': False, 'accepted_identity_claim': False},
    }


def _weighted_stage(parent, weighted, ru, mq):
    old = parent.get('component_result', {}) if isinstance(parent, _Mapping) else {}
    rows = weighted.get('definition_graph', {}).get('nodes') if isinstance(weighted, _Mapping) else None
    weighted_result = weighted.get('component_result', {}) if isinstance(weighted, _Mapping) else {}
    def row_is(number, **expected):
        item = _row(rows, number)
        return item is not None and all(item.get(key) == value for key, value in expected.items())
    predicates = {
        'weighted_graph_shape': isinstance(rows, list) and len(rows) == 545
            and rows == weighted_result.get('definition_graph', {}).get('nodes')
            and weighted_result.get('actual_native_active_S1_weighted_IBP_C1_identity_proved') is True,
        'J3_jump_row': row_is(523, operation='actual_internal_primitive_jump', arguments=[377, 367],
            jump_index=3, formula='J_i=B_next_lower-B_prev_upper', is_actual_function=True,
            retained=True, zero_assumption=False, source_identity_proved=False),
        'B_endpoint_rows': (row_is(367, operation='actual_B_endpoint_function', chart='first_switch',
                endpoint='upper', coordinate=75, domain=76, radius=77, endpoint_bound='domain[1]',
                formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
            and row_is(377, operation='actual_B_endpoint_function', chart='second_switch',
                endpoint='lower', coordinate=91, domain=92, radius=93, endpoint_bound='domain[0]',
                formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
            and _row(rows, 367).get('physical_ports', {}).get('U') == 80
            and _row(rows, 367).get('physical_ports', {}).get('Q') == 84
            and _row(rows, 377).get('physical_ports', {}).get('U') == 96
            and _row(rows, 377).get('physical_ports', {}).get('Q') == 100),
        'source_ports': all(row_is(number, source_nodes=source_nodes) for number, source_nodes in (
            (75, [73190]), (76, [0, 1]), (77, [85905]),
            (80, [148085, 148086, 148087, 148088, 148089, 148090]),
            (84, [148180, 148184, 148190, 148196, 148202]),
            (91, [75949]), (92, [0, 1]), (93, [88600]),
            (96, [154247, 154248, 154249, 154250, 154251, 154252]),
            (100, [154342, 154346, 154352, 154358, 154364]))),
        'source_port_types': all(row_is(number, operation=operation, source_record=source,
                chart=chart, quantity=quantity) for number, operation, source, chart, quantity in (
            (75, 'actual_native_geometry_coordinate', LOOP, 'first_switch', 'coordinate'),
            (77, 'actual_native_geometry_radius', LOOP, 'first_switch', 'radius'),
            (91, 'actual_native_geometry_coordinate', LOOP, 'second_switch', 'coordinate'),
            (93, 'actual_native_geometry_radius', LOOP, 'second_switch', 'radius'),
            (80, 'actual_native_physical_U', MOD, 'first_switch', 'U'),
            (84, 'actual_native_physical_Q', MOD, 'first_switch', 'Q'),
            (96, 'actual_native_physical_U', MOD, 'second_switch', 'U'),
            (100, 'actual_native_physical_Q', MOD, 'second_switch', 'Q'))),
        'common_L_and_measure': (all(row_is(number, source_nodes=[9], source_record=LOOP)
            for number in (78, 79, 94, 95))
            and _row(rows, 367).get('measure') == _row(rows, 377).get('measure') == 'R^2*Jacobian*dnative'
            and weighted_result.get('source_interface', {}).get('actual_L', {}).get('source_record') == MOD
            and weighted_result.get('source_interface', {}).get('actual_L', {}).get('positive_quotient_denominator') == 25272),
        'uniform_B_creator_recipe': old.get('component_proofs', {}).get('J2_B', {}).get('source_ast_hashes', {}).get(
            'active_S1_creator_sha256') == RECORD_HASHES[WEIGHTED_CREATOR]
            and old.get('native_weighted_B_definition', {}).get('formula') ==
            'R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)',
        'RU_and_MQ': ru.get('J3_R_C1_proved') is True and ru.get('J3_U_C1_proved') is True
            and mq.get('J3_Q_C1_proved') is True,
    }
    passed = all(type(value) is bool and value is True for value in predicates.values())
    return {
        'schema': 'native_J3_B_C1_source_composition_v3',
        'status': 'source_B_C1_proved' if passed else 'source_guard_failure',
        'source_only': True, 'accepted_identity_claim': False,
        'J3_B_C1_proved': passed, 'J3_native_weighted_jump_zero_proved': passed,
        'jump_node': 523, 'jump_index': 3, 'jump_arguments': [377, 367],
        'predicates': predicates, 'source_ports': {'left_B': 367, 'right_B': 377,
            'first_coordinate': 73190, 'first_value': 1, 'second_coordinate': 75949,
            'second_value': 0, 'common_L': 25272},
        'source_derivation': {
            'law': 'The weighted graph node 523 is B377 minus B367 with identical source C1 ports on the common original domain.',
            'dependencies': {'R_U': ru.get('J3_R_C1_proved') is True,
                             'M_Q': mq.get('J3_Q_C1_proved') is True,
                             'weighted_node_523': predicates.get('J3_jump_row'),
                             'endpoint_rows': predicates.get('B_endpoint_rows'),
                             'source_ports': predicates.get('source_ports')},
            'no_numeric_or_domain_claim': True,
        },
        'numeric_values_computed': False, 'source_constructor_calls': 0,
        'source_producer_function_calls': 0, 'old_proof_functions_called': 0,
        'generic_DAG_walks': 0,
    }


def _emit(callback, events, stage, payload):
    event = {'stage': stage, 'payload': payload}
    events.append(event)
    if callback is not None:
        callback(stage, payload)


def derive(records, expected_input_digest, diagnostic_callback=None):
    """Authenticate the V2 failure chain and derive the fresh J3 source leaf."""
    require(isinstance(records, dict) and set(records) == set(RECORD_NAMES),
            'Changed record inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    checks, events, decoded = {}, [], {}
    record_checks = {name: _record_ok(name, records[name]) for name in RECORD_NAMES}
    auth_payload = {
        'schema': 'native_J3_v3_saved_failure_artifacts_v1',
        'source_only': True, 'record_checks': dict(record_checks),
        'failure_artifact_inventory': dict(FAILURE_ARTIFACT_INVENTORY),
        'raw_sha256_before_decode': True,
        'saved_V2_code_and_runner': all(record_checks[name] for name in (V2_SOURCE_CODE, V2_RUNNER)),
        'saved_V2_source_entry_result_diagnostics': all(record_checks[name]
            for name in (V2_SOURCE_PACKET, V2_ENTERED, V2_RESULT, V2_DIAGNOSTICS)),
    }
    _emit(diagnostic_callback, events, 'saved_failure_artifacts', auth_payload)
    require(all(record_checks.values()), 'Changed authenticated source record')
    for name in RECORD_NAMES:
        text = records[name]['source_text']
        decoded[name] = text if name.endswith('.py') or name == V2_DIAGNOSTICS else _json_object(text, name)
    saved = _saved_chain(records)
    _emit(diagnostic_callback, events, 'saved_V2_certificate', saved)
    require(saved.get('passed') is True, 'Saved V2 failure chain changed')
    parent = decoded[PARENT]
    accepted = _accepted_parent_ok(parent)
    failed_report = decoded[FAILED_REPORT]
    failed_parent = _failed_parent_ok(failed_report, saved['result'])
    failure_evidence = _failure_evidence_ok(decoded[FAILURE_EVIDENCE], saved['result'])
    source_checks = {
        'accepted_parent.' + key: value for key, value in accepted.items()
    }
    source_checks.update({'failed_parent.' + key: value for key, value in failed_parent.items()})
    source_checks.update({'failure_evidence.' + key: value for key, value in failure_evidence.items()})
    source_pins = parent.get('source_pins', {})
    source_checks['inherited_old_records'] = all(
        isinstance(source_pins, _Mapping) and source_pins.get(name) == RECORD_HASHES[name]
        for name in (LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, C6_CREATOR, BRIDGE_CREATOR,
                     BACKGROUND, WEIGHTED, WEIGHTED_CREATOR, GAMMA_CREATOR))
    source_checks['same_source_family'] = all(
        isinstance(value, _Mapping) and value.get('source_family') == FAMILY
        for value in (decoded[LOOP], decoded[MOD], decoded[BACKGROUND], decoded[WEIGHTED]))
    source_checks['archive_creator_pins'] = (
        decoded[LOOP].get('input_hashes', {}).get(LOOP_CREATOR) == RECORD_HASHES[LOOP_CREATOR]
        and decoded[MOD].get('input_hashes', {}).get(MOD_CREATOR) == RECORD_HASHES[MOD_CREATOR]
        and decoded[BACKGROUND].get('input_hashes', {}).get(C6_CREATOR) == RECORD_HASHES[C6_CREATOR]
        and decoded[BACKGROUND].get('input_hashes', {}).get(BRIDGE_CREATOR) == RECORD_HASHES[BRIDGE_CREATOR])
    source_checks['failed_parent_inventory'] = failed_parent['inventory']
    _emit(diagnostic_callback, events, 'accepted_sources', {
        'schema': 'native_J3_v3_accepted_source_inputs_v1',
        'source_only': True, 'predicates': dict(source_checks),
        'accepted_J2_B_source_pins': 1850,
        'failed_V2_source_pins': INHERITED_SOURCE_COUNT,
        'future_source_union': NEW_SOURCE_UNION_COUNT,
        'failed_artifact_union': FAILED_ARTIFACT_UNION_COUNT,
        'checkpoint_count': INHERITED_CHECKPOINT_COUNT,
    })
    require(all(value is True for value in source_checks.values()),
            'Accepted source or published V2 failure report changed')
    prefix_payload = dict(saved.get('prefix', {}))
    _emit(diagnostic_callback, events, 'boolean_prefix', prefix_payload)
    require(saved['predicates'].get('prefix_positive') is not False
            and saved['predicates'].get('diagnostic_chain') is True,
            'Saved Boolean prefix changed')
    cutoff = _sigma_source_guard(decoded[GAMMA_CREATOR], decoded[C6_CREATOR],
                                 decoded[BACKGROUND], source_pins)
    _emit(diagnostic_callback, events, 'cutoff_source_laws', cutoff)
    require(cutoff.get('passed') is True, 'Exact cutoff endpoint source laws changed')
    raw_hashes = {'loop': records[LOOP]['raw_sha256'], 'mod': records[MOD]['raw_sha256'],
                  'background': records[BACKGROUND]['raw_sha256']}
    ru = _fresh_ru(saved, cutoff, decoded[BACKGROUND], decoded[LOOP], decoded[MOD],
                   decoded[C6_CREATOR], decoded[BRIDGE_CREATOR], decoded[LOOP_CREATOR],
                   decoded[MOD_CREATOR], raw_hashes, source_pins)
    _emit(diagnostic_callback, events, 'R_U', ru)
    require(ru.get('J3_R_C1_proved') is True and ru.get('J3_U_C1_proved') is True,
            'Fresh R_U source certificate failed')
    mq = _fresh_mq(saved, ru, decoded[LOOP], decoded[MOD], decoded[BACKGROUND], decoded[MOD_CREATOR])
    _emit(diagnostic_callback, events, 'M_Q', mq)
    require(mq.get('J3_Q_C1_proved') is True and mq.get('finite_N_m_Z2_proved') is True,
            'Fresh m_Q source certificate failed')
    weighted = _weighted_stage(parent, decoded[WEIGHTED], ru, mq)
    _emit(diagnostic_callback, events, 'B_C1', weighted)
    require(weighted.get('J3_B_C1_proved') is True, 'Fresh B523 source composition failed')
    old = parent.get('component_result', {})
    kept = parent.get('retained_construction_data')
    stage_flags = dict.fromkeys(STAGE_FLAGS, True)
    false_flags = dict.fromkeys(FALSE_FLAGS, False)
    graph = [
        {'id': 'saved_V2_failure_chain', 'rule': 'authenticate the saved source, entry, result, and diagnostic bytes'},
        {'id': 'R_U_transport', 'rule': 'combine saved positive rows with the repaired MOD memory key and fresh bridge fingerprint'},
        {'id': 'm_Q_transport', 'rule': 'preserve finite-N memory before physical m and positive-L Q rows', 'depends_on': ['R_U_transport']},
        {'id': 'B523_composition', 'rule': 'bind actual weighted node 523 to B377 minus B367', 'depends_on': ['R_U_transport', 'm_Q_transport']},
        {'id': 'native_IN_reduction', 'rule': 'substitute only the new source J3 zero into the inherited source identity', 'depends_on': ['B523_composition']},
    ]
    source_checks.update({
        'R_U': ru.get('J3_R_C1_proved') is True,
        'm_Q': mq.get('J3_Q_C1_proved') is True,
        'B523': weighted.get('J3_B_C1_proved') is True,
        'all_stage_sources': all(value is True for value in source_checks.values()),
    })
    result = dict(schema=SCHEMA, component='endpoint_transport', passed=True,
        source_only=True, **stage_flags, **false_flags,
        stage_flags=dict(stage_flags), nonclaims=dict(false_flags),
        J3_B_identity_proved=True, J3_R_identity_proved=True, J3_U_identity_proved=True,
        J3_Q_identity_proved=True, ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        owned_Z_cell=['0', '1/2'], input_digest=expected_input_digest,
        native_six_signed_channels=['H', 'M', 'K', 'A', 'B', 'C'],
        generic_five_channels=['m', 'h', 'k', 'e', 'p'], source_family=FAMILY,
        source_function_identity_only=True,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_checks=source_checks, proof_graph=graph, proof_graph_sha256=digest(graph),
        failure_artifact_inventory=dict(FAILURE_ARTIFACT_INVENTORY),
        component_proofs={'saved_V2_certificate': saved, 'cutoff_source_laws': cutoff,
                          'R_U': ru, 'physical_m_Q': mq, 'B523': weighted},
        diagnostics=list(events), diagnostic_events=list(events), diagnostic_payloads=list(events),
        retained_construction_data=kept, source_binding_scope=old.get('source_binding_scope'),
        retained_integral_term_inventory=parent.get('retained_integral_term_inventory'),
        old_errors_frozen=parent.get('old_errors_frozen'), quantity_distinctions=old.get('quantity_distinctions'),
        preserved_term_counts={'native_six_history_channels': ['H', 'M', 'K', 'A', 'B', 'C'],
            'six_signed_inner_lower_ledger_channels_retained': kept.get('six_signed_channels_retained') is True,
            'integral_K_terms': 17, 'outer_B_terms': 2, 'reserved_contribution': 1,
            'all_old_error_data_retained': True},
        inherited_R_identity_jumps=old.get('inherited_R_identity_jumps'),
        B_identity_scope=dict(weighted_jump_node=523, jump_index=3, left_B_node=367, right_B_node=377,
            jump_arguments=[377, 367], ordinary_Z_orders=[0, 1], endpoint_substitutions={'73190': 1, '75949': 0},
            common_L_node=25272, newly_proved_B_identity_nodes=[523], inherited_B_identity_nodes=[520, 521, 522],
            cumulative_B_identity_nodes=[520, 521, 522, 523], unresolved_B_identity_nodes=list(range(524, 536)),
            remaining_B_identity_count=12, all_native_B_jumps_zero_proved=False),
        unresolved_B_identity_nodes=list(range(524, 536)), remaining_B_identity_nodes=list(range(524, 536)),
        remaining_B_identity_count=12, native_weighted_B_definition=old.get('native_weighted_B_definition'),
        B_C1_derivative_rule=old.get('component_proofs', {}).get('J2_B', {}).get('C1_derivative_rule'),
        native_IN_inherited_identity=old.get('native_IN_reduced_source_identity',
            'I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)'),
        native_IN_reduced_source_identity='I_N=B_lastupper-B_firstlower-sum(J_4..J_15)+sum17(integral K_i dy)',
        native_IN_transition='INsumJ3 ->sumJ4 sourceonly',
        native_IN_reduction_is_source_only=True, archived_all_16_J_terms_retained=True,
        source_jump_zero_does_not_delete_quantitative_errors=True, unresolved_gates=[],
        remaining_construction_gates=[
            'the other twelve native B jumps and same-source quantitative packet transport',
            'independent terminal functional moments and retained quantitative errors',
            'numerical containment and domain admission',
            'heat stress flat remainder time recursion oscillatory correction independent Cartesian validation',
        ],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_source_consumer_calls=0,
            accepted_ancestor_operator_calls=0, old_proof_functions_called=0, generic_DAG_walks=0,
            numerical_evaluations=0, quadratures=0, recurrence_calls=0, runtime_model_metadata=None),
        runtime_model_metadata=None)
    return result
