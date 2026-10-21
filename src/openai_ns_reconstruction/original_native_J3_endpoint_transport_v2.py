"""Unexecuted native J3 successor source kernel.

This successor consumes the accepted J2 B report, the two immutable records
that preserve the failed J3 attempt, and the same C6/bridge/loop/modified
archives used by the bounded source seam.  The failed attempt is projected to
Boolean prefix facts after authenticating its five raw artifacts.  Downstream
guards consume that certificate directly; no old runtime-shaped state object
is reconstructed and no frozen producer or proof function is imported.
"""

import ast as _ast
import base64 as _base64
import hashlib as _hashlib
import json as _json
from collections.abc import Mapping as _Mapping


SCHEMA = 'original_native_J3_endpoint_transport_v2'
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
FAILED_REPORT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J3_endpoint_transport.json'
FAILURE_EVIDENCE = 'experiments/root_st073/original_native_J3_endpoint_transport_failure_evidence.json'
GAMMA_CREATOR = 'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py'
RECORD_NAMES = [PARENT, LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, C6_CREATOR,
                BRIDGE_CREATOR, BACKGROUND, WEIGHTED, WEIGHTED_CREATOR,
                FAILED_REPORT, FAILURE_EVIDENCE, GAMMA_CREATOR]
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
    FAILED_REPORT: 'cc276ae3ab77c8cd6730e4f44b3ea60e95ea4af4b1e50d2153915781b26b16d0',
    FAILURE_EVIDENCE: 'f0b3869589559344ca2fa2513f81104e50654d25e6560bf9474c682fdc82319d',
    GAMMA_CREATOR: '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01',
}
STAGE_FLAGS = ['actual_J3_source_ports_bound', 'original_five_moments_Z2_endpoint_identity_proved',
               'finite_N_m_Z2_endpoint_identity_proved', 'J3_R_C1_endpoint_identity_proved',
               'J3_U_C1_endpoint_identity_proved', 'J3_Q_C1_endpoint_identity_proved',
               'J3_B_C1_endpoint_identity_proved', 'J3_native_weighted_jump_zero_proved',
               'native_IN_J3_source_reduction_proved', 'all_old_terms_retained',
               'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['all_native_B_jumps_zero_proved', 'J4_through_J15_B_identities_proved',
               'quantitative_packet_error_reduction_proved', 'numeric_values_computed',
               'J3_Qy_identity_proved', 'source_C5_installation_claimed',
               'numeric_containment_claimed', 'new_loop_domain_admission_claimed',
               'admissible_cone_admission_claimed', 'physical_layer_completed',
               'old_error_terms_removed', 'global_five_moment_closure',
               'time_order_recursion_claimed', 'full_goal_complete',
               'uninterpreted_graph_parameter_identity_proved', 'old_operators_replayed']
FAMILY = {
    'actual_five_defect_family_sha256': '3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894',
    'implicit_source_sha256': '5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae',
    'datum_enclosure_sha256': 'dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010',
}

# The five raw artifacts are the exact frozen source, entry, and result bytes
# captured by the failure evidence.  They are decoded only after this map is
# checked against the evidence inventory.
J3V2__FROZEN_ARTIFACT_SHA256 = {
    'kernel': '65bad0c995cd6e508cb76a7986269dc5312fd0f3761ccaad14515a8863d24956',
    'runner': 'ebbe29d91389b4a531e548ee6f18cc4cd705406a801bbeadf5c89087acf4aebf',
    'source': '5c7a83cafb2cca70ad6625520f1ab54884a5e566e1d215e61c409ba4618f6267',
    'entered': '4a7fa9f818bf54b80b65dd3155be7cf1c386a9d591f1b587bad87897d4a425f3',
    'result': '393cbc84869592cfb0d412225e81ad85054532e410799e463d98e106cf466c0d',
}
J3V2__FROZEN_INPUT_DIGEST = 'dd8a0b64b35fdf59211f0d7cacac16a08c2c4e8efa7a00dbbb6e014d4e61ca22'
J3V2__FROZEN_PARENT_SHA256 = RECORD_HASHES[PARENT]
J3V2__FROZEN_RECORD_HASHES = {name: RECORD_HASHES[name] for name in RECORD_NAMES[:10]}
J3V2__FROZEN_DECODED_RECORD_HASHES = {
    PARENT: '37f20f105dcb58691b32deb7f5fca59f1a25ff7641086614aa61885ddc090655',
    LOOP: 'f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2',
    MOD: '0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75',
    LOOP_CREATOR: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    MOD_CREATOR: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
    C6_CREATOR: '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32',
    BRIDGE_CREATOR: '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605',
    BACKGROUND: '18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a',
    WEIGHTED: 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28',
    WEIGHTED_CREATOR: 'cdf33563c34e211a00128e3901f49875e50bdea3b4cd9224cc00f81e3288ed8d',
}
J3V2__FROZEN_STAGE_FLAGS = [
    'actual_J3_source_ports_bound', 'original_five_moments_Z2_endpoint_identity_proved',
    'finite_N_m_Z2_endpoint_identity_proved', 'J3_R_C1_endpoint_identity_proved',
    'J3_U_C1_endpoint_identity_proved', 'J3_Q_C1_endpoint_identity_proved',
    'J3_B_C1_endpoint_identity_proved', 'J3_native_weighted_jump_zero_proved',
    'native_IN_J3_source_reduction_proved', 'all_old_terms_retained',
    'accepted_ancestors_not_reexecuted']
J3V2__FROZEN_FALSE_FLAGS = [
    'all_native_B_jumps_zero_proved', 'J4_through_J15_B_identities_proved',
    'quantitative_packet_error_reduction_proved', 'numeric_values_computed',
    'J3_Qy_identity_proved', 'source_C5_installation_claimed',
    'numeric_containment_claimed', 'new_loop_domain_admission_claimed',
    'admissible_cone_admission_claimed', 'physical_layer_completed',
    'old_error_terms_removed', 'global_five_moment_closure',
    'time_order_recursion_claimed', 'full_goal_complete',
    'uninterpreted_graph_parameter_identity_proved', 'old_operators_replayed']
J3V2__FROZEN_AST_BODY_SHA256 = {
    'kernel.require': '9f1cfc7e95bb0eb793def0a0d06f35844705272b4df5821841a291ac020db0bf',
    'kernel.prove_J3_endpoint_source_algebra': '8876c6d7d43b062eb87f6d49ec728c83cb67228de0bec6f7e3d6a39ea3c169d9',
    'kernel.derive': '1d0aeddf784783bb843eefc4690aae2d824a5ab1744142db4f228cd84d445635',
    'kernel.J3MQ__original_state_contract_ok': '9c4034ac94c540975df759e276733b1ca604508a9846f53c7f0ae397e3689c16',
    'runner.require': '9f1cfc7e95bb0eb793def0a0d06f35844705272b4df5821841a291ac020db0bf',
    'runner.read_record': '5f461e72d9a7cd8c817745cecabbf2d7af42a35a2b60d9f0f4c724d7c5e39baf',
    'runner.setup': '62e0ccd0980182b82c843357ca4cf5e9e3e9dea5b5962e1e1020a1d770c7bb85',
    'runner.validate': '2b5059da78f93b0207713909f7e11bb0b234e081074cf1d07355d5d3d5b8553d',
    'runner.run': '7c172ede9078023feb671000454c00e9018722a2058914f1a10d7032d904b111',
}
J3V2__FROZEN_MQ_SCHEMA = 'native_J3_m_Q_source_fragment_v2_guarded'
J3V2__FROZEN_MQ_STATUS = 'original_state_source_guard_ready'
J3V2__FROZEN_MQ_GUARD_NAMES = (
    'decoded_background_and_loop_shapes', 'shared_authenticated_source_family',
    'authenticated_C6_constructor_AST_and_calls', 'authenticated_bridge_histories_AST_and_calls',
    'typed_C6_actual_moment_and_loop_chart_bindings', 'C6_coordinate_radius_Pstar_and_L_source_rows',
    'C6_first_switch_upper_actual_U_endpoint_alias', 'C6_first_switch_upper_moment_endpoint_alias',
    'C6_first_switch_original_M_history_through_Z2', 'C6_second_switch_M_history_through_Z2',
    'C6_original_m_fields_use_inverse_Pstar_once', 'C6_second_switch_actual_phi_transport_rows',
    'C6_first_switch_upper_all_native_moment_endpoint_aliases',
    'C6_native_six_history_transport_M_H_K_A_B_C_through_Z2',
    'C6_generic_h_k_e_p_rows_after_native_normalization',
    'C6_common_W_FF_F0_Pstar_radius_source_rows', 'LOOP_shared_C6_prefix_fixed_source_rows')
J3V2__FROZEN_MQ_SOURCE_HASHES = {
    'parent_report_sha256': J3V2__FROZEN_PARENT_SHA256,
    'c6_creator_sha256': RECORD_HASHES[C6_CREATOR],
    'bridge_creator_sha256': RECORD_HASHES[BRIDGE_CREATOR],
    'c6_raw_archive_sha256': RECORD_HASHES[BACKGROUND],
    'loop_raw_archive_sha256': RECORD_HASHES[LOOP],
}
J3V2__ENDPOINT_PROJECTION = {
    'first_switch_upper': 'coordinate node 73190 -> source value 1; radius node 85905',
    'second_switch_lower': 'coordinate node 75949 -> source value 0; radius node 88600',
    'radius_rule': '100*exp(9*1) = 100*exp(9*(1+0))',
    'phase_rule': 'same source fractional_part phase after equal-radius substitution',
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return _hashlib.sha256(_json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def _J3V2__sha(raw):
    return _hashlib.sha256(raw).hexdigest()


def _J3V2__json(raw, label):
    try:
        value = _json.loads(raw.decode('utf-8'))
    except (UnicodeDecodeError, _json.JSONDecodeError) as exc:
        raise ValueError(label + ': invalid UTF-8 JSON') from exc
    if not isinstance(value, dict):
        raise ValueError(label + ': top-level object required')
    return value


def _J3V2__raw(value, label):
    if isinstance(value, bytes):
        raw = value
    elif isinstance(value, bytearray):
        raw = bytes(value)
    elif isinstance(value, _Mapping) and value.get('encoding') == 'base64':
        if set(value) != {'encoding', 'data'} or type(value.get('data')) is not str:
            raise ValueError(label + ': malformed base64 envelope')
        # Bound the encoded representation before allocating a decoded buffer.
        if len(value['data']) > 2_666_668:
            raise ValueError(label + ': bounded encoded artifact size exceeded')
        try:
            raw = _base64.b64decode(value['data'].encode('ascii'), validate=True)
        except (ValueError, UnicodeEncodeError) as exc:
            raise ValueError(label + ': invalid base64') from exc
    else:
        raise TypeError(label + ': require bytes or base64 envelope')
    if len(raw) > 2_000_000:
        raise ValueError(label + ': bounded artifact size exceeded')
    return raw


def _J3V2__row(rows, number):
    if not isinstance(rows, list) or not isinstance(number, int):
        return None
    if number < 0 or number >= len(rows) or not isinstance(rows[number], dict):
        return None
    return rows[number]


def _J3V2__row_is(rows, number, **expected):
    row = _J3V2__row(rows, number)
    return row is not None and all(row.get(key) == value for key, value in expected.items())


def _J3V2__array(value, prefix):
    return isinstance(value, list) and value[:len(prefix)] == list(prefix)


def _J3V2__compact(text):
    return ''.join(str(text).split())


def _J3V2__top_functions(tree):
    return {node.name: node for node in tree.body
            if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))}


def _J3V2__body_sha(node):
    dumped = _ast.dump(node, annotate_fields=True, include_attributes=False)
    return _hashlib.sha256(dumped.encode('utf-8')).hexdigest()


def _J3V2__call_name(call):
    if not isinstance(call, _ast.Call):
        return None
    node, parts = call.func, []
    while isinstance(node, _ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, _ast.Name):
        return None
    parts.append(node.id)
    return '.'.join(reversed(parts))


def _J3V2__direct_call(statement):
    if isinstance(statement, _ast.Expr) and isinstance(statement.value, _ast.Call):
        return statement.value
    if isinstance(statement, _ast.Assign) and isinstance(statement.value, _ast.Call):
        return statement.value
    return None


def _J3V2__first_string(call):
    if isinstance(call, _ast.Call) and call.args:
        value = call.args[0]
        if isinstance(value, _ast.Constant) and type(value.value) is str:
            return value.value
    return None


def _J3V2__assignment_literals(tree):
    result = {}
    for statement in tree.body:
        if isinstance(statement, _ast.Assign):
            for target in statement.targets:
                if isinstance(target, _ast.Name):
                    result[target.id] = statement.value
    return result


def _J3V2__resolve_literal(node, assignments, seen=(), depth=0):
    if depth > 8:
        raise ValueError('literal depth exceeded')
    if isinstance(node, _ast.Constant):
        return node.value
    if isinstance(node, _ast.Name):
        if node.id in seen or node.id not in assignments:
            raise ValueError('unbound literal name')
        return _J3V2__resolve_literal(assignments[node.id], assignments,
                                      seen + (node.id,), depth + 1)
    if isinstance(node, _ast.List):
        return [_J3V2__resolve_literal(item, assignments, seen, depth + 1)
                for item in node.elts]
    if isinstance(node, _ast.Tuple):
        return tuple(_J3V2__resolve_literal(item, assignments, seen, depth + 1)
                     for item in node.elts)
    if isinstance(node, _ast.Dict):
        if any(key is None for key in node.keys):
            raise ValueError('dictionary unpack is not a literal')
        return {_J3V2__resolve_literal(key, assignments, seen, depth + 1):
                _J3V2__resolve_literal(value, assignments, seen, depth + 1)
                for key, value in zip(node.keys, node.values)}
    raise ValueError('nonliteral AST node')


def _J3V2__failure_ok(result):
    trace = result.get('traceback') if isinstance(result.get('traceback'), str) else ''
    markers = ('run_original_native_J3_endpoint_transport.py", line 186, in run',
               'original_native_J3_endpoint_transport.py", line 226, in derive',
               'original_native_J3_endpoint_transport.py", line 170, in check',
               'original_native_J3_endpoint_transport.py", line 59, in require',
               'AssertionError: Source prerequisite: new_R_U_C1')
    return bool(result.get('passed') is False
                and result.get('failed_attempt_preserved') is True
                and result.get('error_type') == 'AssertionError'
                and result.get('error') == 'Source prerequisite: new_R_U_C1'
                and result.get('runner_caught_exception') is True
                and result.get('runner_full_goal_complete') is False
                and all(marker in trace for marker in markers))


def _J3V2__source_packet_ok(source):
    return bool(source.get('schema') == 'native_J3_endpoint_transport_inputs_v1'
                and source.get('component') == 'endpoint_transport'
                and source.get('input_digest') == J3V2__FROZEN_INPUT_DIGEST
                and source.get('record_hashes') == J3V2__FROZEN_RECORD_HASHES
                and source.get('decoded_record_hashes') == J3V2__FROZEN_DECODED_RECORD_HASHES
                and source.get('implementation_hashes') == {
                    'scripts/run_original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['runner'],
                    'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['kernel']}
                and source.get('schema_expected') == 'original_native_J3_endpoint_transport_v1'
                and source.get('stage_flags') == J3V2__FROZEN_STAGE_FLAGS
                and source.get('false_flags') == J3V2__FROZEN_FALSE_FLAGS
                and source.get('ancestor_operators_replayed') is False
                and source.get('kernel_imported_during_setup') is False
                and source.get('numeric_values_computed') is False
                and source.get('full_goal_complete') is False)


def _J3V2__entry_chain_ok(source, entered, result, source_sha, entered_sha):
    return bool(entered.get('component') == 'endpoint_transport'
                and entered.get('source_sha256') == source_sha == J3V2__FROZEN_ARTIFACT_SHA256['source']
                and entered.get('input_digest') == source.get('input_digest') == J3V2__FROZEN_INPUT_DIGEST
                and entered.get('implementation_hashes') == source.get('implementation_hashes')
                and entered.get('kernel_sha256') == J3V2__FROZEN_ARTIFACT_SHA256['kernel']
                and result.get('source_checkpoint_sha256') == source_sha
                and result.get('entry_checkpoint_sha256') == entered_sha
                and result.get('source_checkpoint_sha256') == J3V2__FROZEN_ARTIFACT_SHA256['source']
                and result.get('entry_checkpoint_sha256') == J3V2__FROZEN_ARTIFACT_SHA256['entered'])


def _J3V2__derive_sequence_ok(tree):
    functions = _J3V2__top_functions(tree)
    derive = functions.get('derive')
    require_node = functions.get('require')
    algebra = functions.get('prove_J3_endpoint_source_algebra')
    contract = functions.get('J3MQ__original_state_contract_ok')
    if not all(node is not None for node in (derive, require_node, algebra, contract)):
        return False
    if any(_J3V2__body_sha(node) != J3V2__FROZEN_AST_BODY_SHA256[key]
           for node, key in ((require_node, 'kernel.require'),
                             (algebra, 'kernel.prove_J3_endpoint_source_algebra'),
                             (derive, 'kernel.derive'),
                             (contract, 'kernel.J3MQ__original_state_contract_ok'))):
        return False
    direct = [(i, _J3V2__direct_call(statement), statement)
              for i, statement in enumerate(derive.body)]
    calls = [(i, call) for i, call, _ in direct]
    algebra_calls = [i for i, call in calls if _J3V2__call_name(call) == 'prove_J3_endpoint_source_algebra']
    original_calls = [i for i, call in calls if _J3V2__call_name(call) == 'prove_original_switch_transport']
    ru_calls = [(i, call) for i, call in calls if _J3V2__call_name(call) == 'prove_J3_R_U_transport']
    checks = [(i, call, statement) for i, call, statement in direct
              if _J3V2__call_name(call) == 'check']
    fresh = [(i, call, statement) for i, call, statement in checks
             if _J3V2__first_string(call) == 'fresh_original_switch_certificate_contract']
    failed = [(i, call, statement) for i, call, statement in checks
              if _J3V2__first_string(call) == 'new_R_U_C1']
    if not (len(algebra_calls) == len(original_calls) == len(ru_calls) == 1
            and len(fresh) == len(failed) == 1):
        return False
    if not (algebra_calls[0] < original_calls[0] < fresh[0][0] < ru_calls[0][0] < failed[0][0]):
        return False
    if fresh[0][2].lineno != 213 or failed[0][2].lineno != 226:
        return False
    if len(fresh[0][1].args) < 2 or _J3V2__call_name(fresh[0][1].args[1]) != 'J3MQ__original_state_contract_ok':
        return False
    if len(failed[0][1].args) != 2 or _J3V2__first_string(failed[0][1]) != 'new_R_U_C1':
        return False
    if any(isinstance(statement, _ast.Return) for statement in derive.body[:failed[0][0]]):
        return False
    nested = next((node for node in derive.body
                   if isinstance(node, _ast.FunctionDef) and node.name == 'check'), None)
    return bool(nested is not None and nested.body
                and _J3V2__call_name(_J3V2__direct_call(nested.body[0])) == 'require')


def _J3V2__runner_sequence_ok(tree):
    functions = _J3V2__top_functions(tree)
    required = ('require', 'read_record', 'setup', 'validate', 'run')
    if any(functions.get(name) is None for name in required):
        return False
    if any(_J3V2__body_sha(functions[name]) != J3V2__FROZEN_AST_BODY_SHA256['runner.' + name]
           for name in required):
        return False
    run = functions['run']
    direct = [(i, _J3V2__direct_call(statement), statement)
              for i, statement in enumerate(run.body)]
    try_nodes = [statement for _, _, statement in direct if isinstance(statement, _ast.Try)]
    writes = [i for i, call, _ in direct if _J3V2__call_name(call) == 'write_new']
    if len(try_nodes) != 1 or len(writes) != 2:
        return False
    try_node = try_nodes[0]
    names = [_J3V2__call_name(_J3V2__direct_call(statement)) for statement in try_node.body]
    expected = ['importlib.util.spec_from_file_location', 'importlib.util.module_from_spec',
                'spec.loader.exec_module', 'module.derive']
    if names != expected or len(try_node.handlers) != 1:
        return False
    handler = try_node.handlers[0]
    return bool(isinstance(handler.type, _ast.Name) and handler.type.id == 'Exception'
                and writes[0] < next(i for i, _, statement in direct if statement is try_node) < writes[1])


def _J3V2__mq_static_ok(tree, text):
    functions = _J3V2__top_functions(tree)
    contract = functions.get('J3MQ__original_state_contract_ok')
    if contract is None or _J3V2__body_sha(contract) != J3V2__FROZEN_AST_BODY_SHA256['kernel.J3MQ__original_state_contract_ok']:
        return False
    assignments = _J3V2__assignment_literals(tree)
    try:
        schema = _J3V2__resolve_literal(assignments['J3MQ_SCHEMA'], assignments)
        guards = _J3V2__resolve_literal(assignments['J3MQ__ORIGINAL_GUARD_NAMES'], assignments)
        hashes = _J3V2__resolve_literal(assignments['J3MQ__ORIGINAL_SOURCE_HASHES'], assignments)
    except (KeyError, ValueError, TypeError):
        return False
    needles = (
        "original_state.get('schema') == J3MQ_SCHEMA",
        "original_state.get('status') == 'original_state_source_guard_ready'",
        "original_state.get('source_only') is True",
        "original_state.get('accepted_identity_claim') is False",
        "original_state.get('original_state_Z2_proved') is True",
        "original_state.get('original_five_moments_Z2_proved') is True",
        "original_state.get('ordinary_Z_orders') == [0, 1, 2]",
        "original_state.get('owned_cell') == {'lower': 0, 'upper': '1/2'}",
        "set(guards) == set(J3MQ__ORIGINAL_GUARD_NAMES)",
        "all((guards.get(name) is True for name in J3MQ__ORIGINAL_GUARD_NAMES))")
    return bool(schema == J3V2__FROZEN_MQ_SCHEMA
                and guards == J3V2__FROZEN_MQ_GUARD_NAMES
                and hashes == J3V2__FROZEN_MQ_SOURCE_HASHES
                and all(_J3V2__compact(needle) in _J3V2__compact(text) for needle in needles)
                and _J3V2__compact('J3MQ__original_state_contract_ok(original)') in _J3V2__compact(text)
                and J3V2__FROZEN_MQ_STATUS.replace(' ', '') in _J3V2__compact(text))


def _J3V2__static_prefix_ok(kernel_text, runner_text):
    try:
        kernel_tree = _ast.parse(kernel_text)
        runner_tree = _ast.parse(runner_text)
        return {
            'ast.kernel_runner_body_pins': _J3V2__derive_sequence_ok(kernel_tree)
                and _J3V2__runner_sequence_ok(runner_tree),
            'ast.fail_fast_require_check_derive_sequence': _J3V2__derive_sequence_ok(kernel_tree),
            'ast.runner_entry_try_result_sequence': _J3V2__runner_sequence_ok(runner_tree),
            'old_J3MQ.contract_static_shape': _J3V2__mq_static_ok(kernel_tree, kernel_text),
            'source_algebra.static_R_phase_rules': all(
                _J3V2__compact(needle) in _J3V2__compact(kernel_text) for needle in (
                    "check('endpoint_algebra:actual_fixed_radius_phase_spine'",
                    "check('endpoint_algebra:typed_endpoint_domains'",
                    "check('endpoint_algebra:Z_independent_geometry_parameters'",
                    "check('endpoint_algebra:radius_congruence_after_source_substitution'",
                    "check('endpoint_algebra:global_phase_congruence_after_source_substitution'",
                    'left_bound == 1 and right_bound == 0', 'left_radius == right_radius',
                    'left_phase == right_phase')),
        }
    except (SyntaxError, TypeError, ValueError, KeyError):
        return {name: False for name in ('ast.kernel_runner_body_pins',
            'ast.fail_fast_require_check_derive_sequence', 'ast.runner_entry_try_result_sequence',
            'old_J3MQ.contract_static_shape', 'source_algebra.static_R_phase_rules')}


def J3V2__project_boolean_prefix(kernel_raw, runner_raw, source_raw, entered_raw,
                                 result_raw, parent_pins, check=None):
    """Authenticate the frozen five raw bytes and return Boolean prefix facts.

    This is a static evidence projector.  It parses only the authenticated
    frozen source and runner text and never imports, executes, or reconstructs
    any frozen runtime object.
    """
    raw = {'kernel': _J3V2__raw(kernel_raw, 'kernel'),
           'runner': _J3V2__raw(runner_raw, 'runner'),
           'source': _J3V2__raw(source_raw, 'source'),
           'entered': _J3V2__raw(entered_raw, 'entered'),
           'result': _J3V2__raw(result_raw, 'result')}
    actual = {name: _J3V2__sha(value) for name, value in raw.items()}
    predicates = {
        'pins.parent_shape': (isinstance(parent_pins, _Mapping)
            and set(parent_pins) == {'artifact_raw_sha256', 'implementation_hashes',
                                     'input_digest', 'record_hashes',
                                     'decoded_record_hashes', 'ast_body_sha256'})
            and parent_pins.get('artifact_raw_sha256') == J3V2__FROZEN_ARTIFACT_SHA256
            and parent_pins.get('implementation_hashes') == {
                'scripts/run_original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['runner'],
                'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['kernel']}
            and parent_pins.get('input_digest') == J3V2__FROZEN_INPUT_DIGEST
            and parent_pins.get('record_hashes') == J3V2__FROZEN_RECORD_HASHES
            and parent_pins.get('decoded_record_hashes') == J3V2__FROZEN_DECODED_RECORD_HASHES
            and parent_pins.get('ast_body_sha256') == J3V2__FROZEN_AST_BODY_SHA256,
        'artifacts.raw_sha256_exact': actual == J3V2__FROZEN_ARTIFACT_SHA256,
    }
    if predicates['artifacts.raw_sha256_exact'] is True:
        try:
            source = _J3V2__json(raw['source'], 'source')
            entered = _J3V2__json(raw['entered'], 'entered')
            result = _J3V2__json(raw['result'], 'result')
            predicates['artifacts.source_packet_shape'] = _J3V2__source_packet_ok(source)
            predicates['artifacts.entry_result_chain'] = _J3V2__entry_chain_ok(
                source, entered, result, actual['source'], actual['entered'])
            predicates['frozen_failure.new_R_U_C1'] = _J3V2__failure_ok(result)
            static = _J3V2__static_prefix_ok(raw['kernel'].decode('utf-8'), raw['runner'].decode('utf-8'))
            predicates.update(static)
        except (UnicodeDecodeError, _json.JSONDecodeError, ValueError, TypeError, KeyError):
            predicates.update({
                'artifacts.source_packet_shape': False, 'artifacts.entry_result_chain': False,
                'frozen_failure.new_R_U_C1': False, 'ast.kernel_runner_body_pins': False,
                'ast.fail_fast_require_check_derive_sequence': False,
                'ast.runner_entry_try_result_sequence': False,
                'old_J3MQ.contract_static_shape': False,
                'source_algebra.static_R_phase_rules': False})
    else:
        predicates.update({
            'artifacts.source_packet_shape': False, 'artifacts.entry_result_chain': False,
            'frozen_failure.new_R_U_C1': False, 'ast.kernel_runner_body_pins': False,
            'ast.fail_fast_require_check_derive_sequence': False,
            'ast.runner_entry_try_result_sequence': False,
            'old_J3MQ.contract_static_shape': False,
            'source_algebra.static_R_phase_rules': False})
    predicates['prefix.source_algebra_R_phase'] = bool(
        predicates.get('source_algebra.static_R_phase_rules') is True
        and predicates.get('ast.fail_fast_require_check_derive_sequence') is True
        and predicates.get('frozen_failure.new_R_U_C1') is True)
    predicates['prefix.original_native6_generic5_Z2'] = bool(
        predicates.get('old_J3MQ.contract_static_shape') is True
        and predicates.get('ast.fail_fast_require_check_derive_sequence') is True
        and predicates.get('frozen_failure.new_R_U_C1') is True)
    predicates['prefix.fresh_certificate_before_RU'] = bool(
        predicates.get('old_J3MQ.contract_static_shape') is True
        and predicates.get('ast.fail_fast_require_check_derive_sequence') is True
        and predicates.get('frozen_failure.new_R_U_C1') is True)
    predicates = {name: value is True for name, value in predicates.items()}
    missing = next((name for name, value in predicates.items() if value is not True), None)
    accepted = missing is None
    if accepted and callable(check):
        try:
            accepted = check(dict(predicates)) is True
        except Exception:
            accepted = False
            missing = 'projector_check_raised'
        if not accepted and missing is None:
            missing = 'projector_check_rejected'
    return {
        'schema': 'native_J3_failed_boolean_prefix_projector_v1',
        'status': 'boolean_prefix_projected' if accepted else 'missing_gate',
        'missing_gate': missing,
        'source_algebra_R_phase_prefix_proved': accepted and predicates.get('prefix.source_algebra_R_phase') is True,
        'original_native6_generic5_Z2_prefix_proved': accepted and predicates.get('prefix.original_native6_generic5_Z2') is True,
        'fresh_certificate_before_RU': accepted and predicates.get('prefix.fresh_certificate_before_RU') is True,
        'new_R_U_C1_checked_and_failed': predicates.get('frozen_failure.new_R_U_C1') is True,
        'ordinary_Z_metadata': ({'ordinary_Z_orders': [0, 1],
                                 'original_Z2_premise_orders': [0, 1, 2],
                                 'owned_Z_cell': {'lower': '0', 'upper': '1/2'},
                                 'ordinary_Z_not_temporal_n': True,
                                 'proved': True} if accepted else
                                {'ordinary_Z_orders': None,
                                 'original_Z2_premise_orders': None,
                                 'owned_Z_cell': None,
                                 'ordinary_Z_not_temporal_n': None,
                                 'proved': False}),
        'projection': (dict(J3V2__ENDPOINT_PROJECTION) if accepted else
                       {'proved': False}),
        'source_hashes': {'artifacts': actual, 'implementation': dict(J3V2__FROZEN_ARTIFACT_SHA256),
                          'records': dict(J3V2__FROZEN_RECORD_HASHES),
                          'decoded_records': dict(J3V2__FROZEN_DECODED_RECORD_HASHES),
                          'ast_body': dict(J3V2__FROZEN_AST_BODY_SHA256)},
        'predicates': predicates,
        'provenance': {'static_only': True, 'runtime_intermediate_objects_reconstructed': False,
                       'source_constructor_calls': 0, 'source_producer_function_calls': 0,
                       'old_proof_functions_called': 0, 'generic_DAG_walks': 0,
                       'numeric_evaluations': 0, 'fragment_was_unexecuted_before_frozen_entry': True,
                       'failure_path': 'runner.run -> module.derive -> check(new_R_U_C1) -> require'},
        'runtime_model_metadata': None,
    }


def _J3V2__inventory_raw(evidence, path, expected):
    inventory = evidence.get('artifact_inventory') if isinstance(evidence, dict) else None
    if not isinstance(inventory, list):
        return None
    item = next((entry for entry in inventory if isinstance(entry, dict)
                 and entry.get('path') == path), None)
    if item is None or item.get('sha256') != expected or item.get('encoding') != 'base64':
        return None
    encoded = item.get('base64_bytes')
    if not isinstance(encoded, str) or len(encoded) > 2_666_668:
        return None
    try:
        raw = _base64.b64decode(encoded.encode('ascii'), validate=True)
    except (KeyError, ValueError, UnicodeEncodeError):
        return None
    return raw if len(raw) == item.get('bytes') and _J3V2__sha(raw) == expected else None


# These body digests are copied from the bounded source guard.  The bridge
# ``generic_fields`` needle below intentionally uses the repaired algebra API
# call; the creator body digest remains the frozen value.
J3V2__CREATOR_AST_SHA256 = {
    'c6': {
        ('function', 'first_switch_fields'): '0bd697b785e7241b2ea8520d74267342e4993f23dffc90fe49617d13430d68e2',
        ('function', 'second_switch_fields'): 'c6c08d15f7b5058f9ff48358ea904eed294d21046ca016deb7d1c440abd4d438',
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', '__init__'): 'd7e974bcc4e2156f8c05b3a67192aefe1df25cb0857df1e7a47e8bd1f6eb2b32',
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'add_chart'): '18f1dbdd7128ac2e72098ab5d7249659923a5d7166640057b2587acfcf2f3812',
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'end'): '3a9ca663af71f7dcf34fa945662ea1535de41a0d8bc2a53a2aacf2954d5b6040',
    },
    'loop': {
        ('function', 'loop_functions'): 'fec28c8c3ff3bf20a5cc64bc777b0b0ed304387b95f71cb6a8f9ec30fd9e0559',
        ('function', 'generic_roots'): '3ba21f81b7aae33751e6534e15e97311c1b851f103b522259a80d1170fad7847',
        ('function', 'finite_N_functions'): '219500c1ac0665e869303a421303a9d2a4c05bfa3121a9fe1f06f022649ed4d2',
        ('function', 'transport'): 'fd5f45c38d246f2555c045ca2d386f5b49bbfaa7699981156efbf99c1869e9b3',
        ('method', 'CurrentOriginalPrescribedLoopC5Functions', '__init__'): '4e6a13cf1025d85e3bf0d7c5abfe9b90d89e9cf45692ee6f5b6161b4ed8546c2',
    },
    'mod': {
        ('function', 'physical_pressure_fields'): '619f51980979fc35ad965386a0038fc2785bf15c99ef2a1b61d749a2d6833a14',
        ('method', 'CurrentOriginalModifiedQ4PressureC3Functions', '__init__'): '8040fa64b5da2de7fa1e2e29cd2583bb803b5b0735d36357050200ad87d645d4',
    },
    'bridge': {
        ('function', 'histories'): '9dc010e43b2a9b1c498380b4452f226f467df7f5529f4c026ed5b8f7a6aaaf5d',
        ('function', 'generic_fields'): '41f3b37e36843799b913ec6af1da1d59ef8dcd668c34a0b3584621fe7e791cdf',
        ('method', 'CurrentOriginalPrescribedBridgeC6Functions', '__init__'): 'df1b4b43d63a07bc8e44fa945662ea1535de41a0d8bc2a53a2aacf2954d5b6040',
    },
}
# Correct bridge constructor digest is kept as the authenticated frozen value.
J3V2__CREATOR_NEEDLES = {
    'c6': {
        ('function', 'first_switch_fields'): ('Uy=a.scale(a.mul(ratio,drive),g.neg(chi_axial))', 'chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))', "U=a.add(initial['U'],a.integral(", 'source_measure_applied_once=True', 'return dict(phi=phi,U=U,phi_y=phiy,U_y=Uy'),
        ('function', 'second_switch_fields'): ('sig=sigma(g,phase)', 'original_second_switch_full_current_radius_direction_and_power=True', "return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)"),
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', '__init__'): ('known,barphi=self.comparison_at_radius(R)', "initial={key:bridge_owner.R100['actual'][key] for key in ('phi','U')}", "self.add_chart('first_switch',s,R,hb,actual,bridge_owner.R100['moments'],known=known)", "self.endpoints['R1']=self.end('first_switch',g.one)", "incoming=self.endpoints['R1']", "actual=second_switch_fields(a,known,incoming['actual'],t,hb,R)", "self.add_chart('second_switch',t,R,hb,actual,incoming['moments'],known=known)"),
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'add_chart'): ("J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,", 'fields,radial,generic=previous.generic_fields('),
        ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'end'): ('Z_independent_substitution=True', "actual={key:at(chart['actual'][key])", "fields={k:at(jet) for k,jet in chart['fields'].items()", "generic={k:Algebra(self.g,5).at(jet,chart['phase'],value) for k,jet in chart['generic'].items()}"),
    },
    'bridge': {
        ('function', 'histories'): ('def histories(a,phi,U,initial,phase,jacobian,label):', 'bodies=densities(a,phi,U);out={}', 'for name,rate in RATES.items():', 'actual_six_history_rate=rate', 'original_incoming_memory_retained=True', 'source_measure_applied_once=True', 'prescribed_bridge_chart=label'),
        # Exact stale-needle repair required by the new V2 guard.
        ('function', 'generic_fields'): ("fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",),
    },
    'loop': {
        ('function', 'loop_functions'): ("shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)", 'g.quotient(g.one,eta,\'same positive canonical eta\')', "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))", "inverse=g.node('original_monotone_phase_inverse'", 'return dict(q=q,inverse=psi,A=A,B_over_Pstar=B', 'angle_endpoints_independent_of_Z=True'),
        ('function', 'generic_roots'): ("E,V=(a.trim(field[key]) for key in ('E','V'))", "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)", 'aa,bb=divide(C),divide(B)', "p1=divide(a.scale(a.add(generic['theta_linear'],a.scale(generic['theta_quadratic'],Pstar)),chart['radius']))", "p2=generic['p2']"),
        ('function', 'finite_N_functions'): ("s=a.scale(loop['B_over_Pstar'],invN)", 'return dict(E_N=a.mul(E,ex),V_N=a.add(V,s),'),
        ('function', 'transport'): ('actual_finite_N_signed_correction=True', 'physical_Jacobian_applied_once=True', 'predecessor_memory_never_reset=True'),
        ('method', 'CurrentOriginalPrescribedLoopC5Functions', '__init__'): ("self.eta=self.unique('current_original_source_parameter','eta')", "self.dstar=self.unique('current_original_source_parameter','d_star')", 'phase=g.unary(\'fractional_part\',g.mul(self.N,offset))', 'roots=generic_roots(g,chart,self.Pstar)', 'loop=loop_functions(g,roots,phase,self.eta,self.dstar,name,coordinate,self.sc)', 'finite=finite_N_functions(g,roots,loop,self.N)', "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)", 'transport(g,finite[\'five_signed_density_increments\'],incoming,', 'fields={key:a.trim(jet) for key,jet in chart[\'fields\'].items()}', 'fields.update(E=finite[\'E_N\'],V=finite[\'V_N\'])', 'self.endpoints[name]={key:a.at(jet,coordinate,upper) for key,jet in history.items()}'),
    },
    'mod': {
        ('function', 'physical_pressure_fields'): ("U=a.scale(fields['V'],Pstar)", "m=a.scale(fields['m'],Pstar)", 'L=b.sub(one,b.mul(de,b.mul(z,z)))', 'same original L=1-delta*Z^2', 'Q=div(numerator(U,m));Qy=div(numerator(Uy,my))', 'return dict(U=U,U_y=Uy,m=m,m_y=my,Q=Q,Q_y=Qy)'),
        ('method', 'CurrentOriginalModifiedQ4PressureC3Functions', '__init__'): ("chart=loop_owner.charts[name];bg=chart['source_chart'];coordinate=chart['coordinate'];jac=chart['physical_measure']", 'physical=physical_pressure_fields(g,chart[\'modified_fields\'],radial,self.background.delta,self.background.Pstar)', 'self.charts[name]=dict(source_chart=chart'),
    },
}


def _J3V2__definition(tree, spec):
    if spec[0] == 'function':
        return next((node for node in tree.body
                     if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))
                     and node.name == spec[1]), None)
    if spec[0] == 'method':
        cls = next((node for node in tree.body if isinstance(node, _ast.ClassDef)
                    and node.name == spec[1]), None)
        if cls is None:
            return None
        return next((node for node in cls.body
                     if isinstance(node, (_ast.FunctionDef, _ast.AsyncFunctionDef))
                     and node.name == spec[2]), None)
    return None


def _J3V2__source_body_sha(node):
    module = _ast.Module(body=list(node.body), type_ignores=[])
    dumped = _ast.dump(module, annotate_fields=True, include_attributes=False)
    return _hashlib.sha256(dumped.encode('utf-8')).hexdigest()


def _J3V2__creator_ast_guard(text, family):
    """Return source and repaired-needle predicates without executing text."""
    predicates = {}
    expected_sha = RECORD_HASHES[{'c6': C6_CREATOR, 'bridge': BRIDGE_CREATOR,
                                 'loop': LOOP_CREATOR, 'mod': MOD_CREATOR}[family]]
    predicates[family + '.source_sha256'] = (isinstance(text, str)
        and _hashlib.sha256(text.encode('utf-8')).hexdigest() == expected_sha)
    try:
        tree = _ast.parse(text, filename=family)
    except (SyntaxError, TypeError, ValueError):
        predicates[family + '.ast_parse'] = False
        return predicates
    predicates[family + '.ast_parse'] = True
    for spec, expected_body_sha in J3V2__CREATOR_AST_SHA256[family].items():
        node = _J3V2__definition(tree, spec)
        name = '.'.join(spec)
        segment = _ast.get_source_segment(text, node) if node is not None else ''
        predicates[family + '.body_present.' + name] = node is not None
        predicates[family + '.body_sha256.' + name] = node is not None and _J3V2__source_body_sha(node) == expected_body_sha
        predicates[family + '.body_needles.' + name] = node is not None and all(
            needle in segment for needle in J3V2__CREATOR_NEEDLES.get(family, {}).get(spec, ()))
    return predicates


def _J3V2__archive_shape(archive, wrapper, node_count):
    return (isinstance(archive, _Mapping)
            and isinstance(archive.get('exact_graph_nodes'), list)
            and len(archive['exact_graph_nodes']) == node_count
            and isinstance(archive.get(wrapper), _Mapping))


def _J3V2__archive_guards(background, loop, mod, archive_raw_sha256,
                           prefix_certificate=None):
    """Authenticate bounded source rows needed by the R/U and m/Q guards."""
    p = {}
    p['archive.background_shape'] = _J3V2__archive_shape(
        background, 'prescribed_R100_Rm_C6_bindings', 97295)
    p['archive.loop_shape'] = _J3V2__archive_shape(
        loop, 'prescribed_loop_C5_bindings', 121834)
    p['archive.mod_shape'] = _J3V2__archive_shape(
        mod, 'modified_Q4_pressure_C3_bindings', 183700)
    p['archive.raw_sha256'] = (isinstance(archive_raw_sha256, _Mapping)
        and archive_raw_sha256.get('background') == RECORD_HASHES[BACKGROUND]
        and archive_raw_sha256.get('loop') == RECORD_HASHES[LOOP]
        and archive_raw_sha256.get('mod') == RECORD_HASHES[MOD])
    if not all(p.values()):
        return p
    try:
        bg = background['prescribed_R100_Rm_C6_bindings']
        lp = loop['prescribed_loop_C5_bindings']
        mp = mod['modified_Q4_pressure_C3_bindings']
        bg_rows, lp_rows, mp_rows = (background['exact_graph_nodes'],
                                      loop['exact_graph_nodes'], mod['exact_graph_nodes'])
        families = [background.get('source_family'), loop.get('source_family'), mod.get('source_family')]
        p['source_family.shared'] = all(isinstance(family, dict) and family == FAMILY for family in families)
        p['source_family.shape'] = all(isinstance(family, dict) and len(family) == 3
                                       and all(isinstance(k, str) and isinstance(v, str) and len(v) == 64
                                               for k, v in family.items()) for family in families)
        first, second = bg['charts']['first_switch'], bg['charts']['second_switch']
        first_l, second_l = lp['charts']['first_switch'], lp['charts']['second_switch']
        first_m, second_m = mp['charts']['first_switch'], mp['charts']['second_switch']
        p['bindings.C6_chart_types'] = (
            first.get('phase') == 73190 and second.get('phase') == 75949
            and first.get('domain') == second.get('domain') == [0, 1]
            and first.get('radius') == 85905 and second.get('radius') == 88600
            and first.get('physical_measure') == second.get('physical_measure') == 9
            and first_l.get('coordinate') == 73190 and second_l.get('coordinate') == 75949
            and first_l.get('domain') == second_l.get('domain') == [0, 1]
            and first_l.get('radius') == 85905 and second_l.get('radius') == 88600)
        p['bindings.C6_actual_and_M_rows'] = (
            _J3V2__array(first.get('actual', {}).get('U'), [87154, 87155])
            and _J3V2__array(second.get('actual', {}).get('U'), [88437, 88438])
            and _J3V2__array(bg.get('endpoints', {}).get('R1', {}).get('actual', {}).get('U'), [88437, 88438])
            and _J3V2__array(first.get('moments', {}).get('M'), [87392, 87393, 87394])
            and _J3V2__array(second.get('moments', {}).get('M'), [89910, 89911, 89912])
            and _J3V2__array(bg.get('endpoints', {}).get('R1', {}).get('moments', {}).get('M'), [88465, 88466, 88467]))
        p['bindings.native_six_channels'] = all(
            _J3V2__array(first.get('moments', {}).get(name), prefix)
            and _J3V2__array(second.get('moments', {}).get(name), suffix)
            for name, prefix, suffix in (
                ('H', [87352, 87353, 87354], [89870, 89871, 89872]),
                ('M', [87392, 87393, 87394], [89910, 89911, 89912]),
                ('K', [87427, 87428, 87429], [89945, 89946, 89947]),
                ('A', [87462, 87463, 87464], [89980, 89981, 89982]),
                ('B', [87497, 87498, 87499], [90015, 90016, 90017]),
                ('C', [87525, 87526, 87527], [90043, 90044, 90045])))
        p['bindings.generic_five_channels'] = all(
            _J3V2__array(first.get('fields', {}).get(name), left)
            and _J3V2__array(second.get('fields', {}).get(name), right)
            for name, left, right in (
                ('m', [87591, 87592, 87593], [90109, 90110, 90111]),
                ('h', [87632, 87633, 87634], [90150, 90151, 90152]),
                ('k', [87673, 87674, 87675], [90191, 90192, 90193]),
                ('e', [87728, 87729, 87730], [90246, 90247, 90248]),
                ('p', [87735, 87738, 87742], [90253, 90256, 90260])))
        p['bindings.C6_fields_to_loop_roots'] = all(
            _J3V2__array(first_l.get('roots', {}).get(name), values)
            and _J3V2__array(second_l.get('roots', {}).get(name), values2)
            for name, values, values2 in (
                ('V', [87584, 87585], [90102, 90103]),
                ('p2', [88423, 88424], [90934, 90935])))
        p['bindings.MOD_physical_ports'] = (
            _J3V2__array(first_m.get('physical', {}).get('U'), [148085, 148086])
            and _J3V2__array(second_m.get('physical', {}).get('U'), [154247, 154248])
            and _J3V2__array(first_m.get('physical', {}).get('m'), [148097, 148098, 148099])
            and _J3V2__array(second_m.get('physical', {}).get('m'), [154259, 154260, 154261])
            and _J3V2__array(first_m.get('physical', {}).get('Q'), [148180, 148184])
            and _J3V2__array(second_m.get('physical', {}).get('Q'), [154342, 154346]))
        p['rows.fixed_geometry'] = all(_J3V2__row_is(lp_rows, number, **expected)
            for number, expected in (
                (73190, {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'}),
                (75949, {'operation': 'bound_variable', 'name': 'true_prescribed_second_switch_t'}),
                (85905, {'operation': 'product', 'arguments': [31709, 85904]}),
                (88600, {'operation': 'product', 'arguments': [31709, 88599]}),
                (33572, {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5}),
                (69245, {'operation': 'positive_quotient', 'numerator': 1, 'denominator': 33572,
                         'source_positive_certificate': 'same strictly positive Pstar'})))
        p['rows.common_L_and_Pstar'] = (
            _J3V2__row_is(mp_rows, 25272, operation='sum', arguments=[1, 25271])
            and _J3V2__row_is(mp_rows, 33572, operation='analytic_unary', name='exp', argument=5)
            and _J3V2__row_is(mp_rows, 69245, operation='positive_quotient', numerator=1,
                              denominator=33572, source_positive_certificate='same strictly positive Pstar'))
        p['rows.U_source_pairs'] = (
            _J3V2__row_is(mp_rows, 148085, operation='product', arguments=[106932, 33572])
            and _J3V2__row_is(mp_rows, 148086, operation='product', arguments=[106933, 33572])
            and _J3V2__row_is(mp_rows, 154247, operation='product', arguments=[109340, 33572])
            and _J3V2__row_is(mp_rows, 154248, operation='product', arguments=[109341, 33572]))
        p['rows.Q_source_pairs'] = (
            _J3V2__row_is(mp_rows, 148180, operation='positive_quotient', numerator=148175,
                          denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
            and _J3V2__row_is(mp_rows, 148184, operation='positive_quotient', numerator=148183,
                              denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
            and _J3V2__row_is(mp_rows, 154342, operation='positive_quotient', numerator=154337,
                              denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
            and _J3V2__row_is(mp_rows, 154346, operation='positive_quotient', numerator=154345,
                              denominator=25272, source_positive_certificate='same original L=1-delta*Z^2'))
    except (KeyError, IndexError, TypeError, AttributeError):
        for name in ('source_family.shared', 'source_family.shape', 'bindings.C6_chart_types',
                     'bindings.C6_actual_and_M_rows', 'bindings.native_six_channels',
                     'bindings.generic_five_channels', 'bindings.C6_fields_to_loop_roots',
                     'bindings.MOD_physical_ports', 'rows.fixed_geometry',
                     'rows.common_L_and_Pstar', 'rows.U_source_pairs', 'rows.Q_source_pairs'):
            p[name] = False
    if prefix_certificate is not None:
        try:
            p.update(_J3V2__full_suffix_guards(
                background, loop, mod, prefix_certificate))
        except (KeyError, IndexError, TypeError, AttributeError, ValueError):
            p['suffix.guard_exception'] = False
    return p


# These are copied source-row recipes from the unexecuted suffix of the
# frozen endpoint attempt.  They are data-only guards: no creator, operator,
# or proof function is called by the V2 kernel.
J3V2__C6_HISTORY_SPECS = (
    ('H', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73043, 73044, 73045), 87320, (87338, 87339, 87340),
     (87345, 87346, 87347), (87352, 87353, 87354), 2),
    ('M', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73050, 73051, 73052), 87360, (87378, 87379, 87380),
     (87385, 87386, 87387), (87392, 87393, 87394), 1),
    ('K', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73057, 73058, 73059), 87320, (87413, 87414, 87415),
     (87420, 87421, 87422), (87427, 87428, 87429), 2),
    ('A', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73064, 73065, 73066), 87360, (87448, 87449, 87450),
     (87455, 87456, 87457), (87462, 87463, 87464), 1),
    ('B', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73071, 73072, 73073), 87320, (87483, 87484, 87485),
     (87490, 87491, 87492), (87497, 87498, 87499), 2),
    ('C', 'true_first_switch', 'true_first_switch_history_t', 73190,
     (73078, 73079, 73080), 87360, (87511, 87512, 87513),
     (87518, 87519, 87520), (87525, 87526, 87527), 1),
    ('H', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88458, 88459, 88460), 89838, (89856, 89857, 89858),
     (89863, 89864, 89865), (89870, 89871, 89872), 2),
    ('M', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88465, 88466, 88467), 89878, (89896, 89897, 89898),
     (89903, 89904, 89905), (89910, 89911, 89912), 1),
    ('K', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88472, 88473, 88474), 89838, (89931, 89932, 89933),
     (89938, 89939, 89940), (89945, 89946, 89947), 2),
    ('A', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88479, 88480, 88481), 89878, (89966, 89967, 89968),
     (89973, 89974, 89975), (89980, 89981, 89982), 1),
    ('B', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88486, 88487, 88488), 89838, (90001, 90002, 90003),
     (90008, 90009, 90010), (90015, 90016, 90017), 2),
    ('C', 'true_second_switch', 'true_second_switch_history_t', 75949,
     (88493, 88494, 88495), 89878, (90029, 90030, 90031),
     (90036, 90037, 90038), (90043, 90044, 90045), 1),
)

J3V2__CORRECTION_SPECS = {
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

J3V2__SIGMA_SOURCE = (
    'lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py')
J3V2__SIGMA_SOURCE_SHA256 = (
    '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01')
J3V2__FLAT_SIGMA_SOURCE = (
    'lei_ren_part1_paper_compliant_flat_pulse_derivatives.py')
J3V2__FLAT_SIGMA_SOURCE_SHA256 = (
    '32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7')
J3V2__SIGMA_CUTOFF_DEFINITION = (
    '0 for x<=0; 1 for x>=1; '
    'exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside')


def _J3V2__integral_ok(row, variable, upper, order, *, rate=None,
                       history_name=None, bridge_chart=None, finite=False):
    if not isinstance(row, dict):
        return False
    expected = {
        'operation': 'definite_integral', 'variable': variable, 'lower': 0,
        'upper': upper, 'exact_function_integral': True,
        'numerical_value_not_installed': True,
        'endpoint_independent_of_Z': True,
        'ordinary_slow_Z_derivative_order': order,
        'true_core_inlet_complete_history': True,
        'numerical_cap_not_defining_integrand': True,
    }
    if any(row.get(key) != value for key, value in expected.items()):
        return False
    if finite:
        return all(row.get(key) == value for key, value in {
            'actual_finite_N_signed_correction': True,
            'physical_Jacobian_applied_once': True,
            'zero_correction_inlet_only_at_original_rminus': True,
            'predecessor_memory_never_reset': True,
        }.items())
    return all(row.get(key) == value for key, value in {
        'actual_six_history_rate': rate,
        'actual_history_name': history_name,
        'original_incoming_memory_retained': True,
        'source_measure_applied_once': True,
        'prescribed_bridge_chart': bridge_chart,
    }.items())


def _J3V2__c6_typed_bindings(background, loop):
    """Check typed C6/bridge fields and the first three ordinary orders."""
    try:
        cb = background['prescribed_R100_Rm_C6_bindings']
        lb = loop['prescribed_loop_C5_bindings']
        first, second = cb['charts']['first_switch'], cb['charts']['second_switch']
        r1 = cb['endpoints']['R1']
        lf, ls = lb['charts']['first_switch'], lb['charts']['second_switch']
        p0 = [32693, 32694, 32695, 32696, 66509, 66510, 66511]
        native = {
            'H': ([87352, 87353, 87354], [89870, 89871, 89872]),
            'M': ([87392, 87393, 87394], [89910, 89911, 89912]),
            'K': ([87427, 87428, 87429], [89945, 89946, 89947]),
            'A': ([87462, 87463, 87464], [89980, 89981, 89982]),
            'B': ([87497, 87498, 87499], [90015, 90016, 90017]),
            'C': ([87525, 87526, 87527], [90043, 90044, 90045]),
        }
        incoming = {
            'H': ([73043, 73044, 73045], [88458, 88459, 88460]),
            'M': ([73050, 73051, 73052], [88465, 88466, 88467]),
            'K': ([73057, 73058, 73059], [88472, 88473, 88474]),
            'A': ([73064, 73065, 73066], [88479, 88480, 88481]),
            'B': ([73071, 73072, 73073], [88486, 88487, 88488]),
            'C': ([73078, 73079, 73080], [88493, 88494, 88495]),
        }
        generic = {
            'm': ([87591, 87592, 87593], [90109, 90110, 90111]),
            'h': ([87632, 87633, 87634], [90150, 90151, 90152]),
            'k': ([87673, 87674, 87675], [90191, 90192, 90193]),
            'e': ([87728, 87729, 87730], [90246, 90247, 90248]),
            'p': ([87735, 87738, 87742], [90253, 90256, 90260]),
        }
        native_ok = all(
            _J3V2__array(first['moments'][name], left)
            and _J3V2__array(second['moments'][name], right)
            and _J3V2__array(first['incoming_histories'][name], incoming[name][0])
            and _J3V2__array(second['incoming_histories'][name], incoming[name][1])
            and _J3V2__array(r1['moments'][name], incoming[name][1])
            and r1['moments'][name][:3] == second['incoming_histories'][name][:3]
            for name, (left, right) in native.items())
        generic_ok = all(
            _J3V2__array(first['fields'][name], left)
            and _J3V2__array(second['fields'][name], right)
            for name, (left, right) in generic.items())
        shape_ok = (
            cb.get('Pstar') == 33572 and _J3V2__array(cb.get('F0'), [35657, 35658, 35661])
            and _J3V2__array(cb.get('P0'), p0)
            and first.get('phase') == 73190 and second.get('phase') == 75949
            and first.get('domain') == second.get('domain') == [0, 1]
            and first.get('radius') == 85905 and second.get('radius') == 88600
            and first.get('physical_measure') == second.get('physical_measure') == 9
            and _J3V2__array(first['actual']['U'], [87154, 87155, 87156])
            and _J3V2__array(first['moments']['M'], [87392, 87393, 87394])
            and _J3V2__array(second['actual']['U'], [88437, 88438, 88439])
            and _J3V2__array(r1['actual']['U'], [88437, 88438, 88439])
            and _J3V2__array(r1['moments']['M'], [88465, 88466, 88467])
            and lf.get('coordinate') == 73190 and ls.get('coordinate') == 75949
            and lf.get('radius') == 85905 and ls.get('radius') == 88600
            and lf.get('domain') == ls.get('domain') == [0, 1]
            and lf.get('physical_measure') == ls.get('physical_measure') == 9
            and lf.get('same_background_pressure_datum') is True
            and ls.get('same_background_pressure_datum') is True
            and lf.get('total_Z_phase_derivative_exact_zero') is True
            and ls.get('total_Z_phase_derivative_exact_zero') is True
            and lf.get('P0') == p0 and ls.get('P0') == p0)
        endpoint = {
            field: [107131, 107132, 107133] if field == 'm' else
                   [107137, 107138, 107139] if field == 'h' else
                   [107143, 107144, 107145] if field == 'k' else
                   [107149, 107150, 107151] if field == 'e' else
                   [107155, 107156, 107157]
            for field in ('m', 'h', 'k', 'e', 'p')}
        endpoints_ok = all(
            _J3V2__array(lb['endpoints']['first_switch'][field], rows)
            and _J3V2__array(lb['endpoints']['second_switch']['incoming'][field], rows)
            for field, rows in endpoint.items())
        return bool(shape_ok and native_ok and generic_ok and endpoints_ok)
    except (KeyError, IndexError, TypeError, AttributeError):
        return False


def _J3V2__c6_native_history_rows(rows):
    for name, chart, variable, upper, incoming, decay, products, integrals, final, rate in J3V2__C6_HISTORY_SPECS:
        if not all(_J3V2__row_is(rows, product, operation='product',
                                  arguments=[source, decay])
                   for product, source in zip(products, incoming)):
            return False
        if not all(_J3V2__integral_ok(_J3V2__row(rows, node), variable,
                                     upper, order, rate=rate,
                                     history_name=name, bridge_chart=chart)
                   for node, order in zip(integrals, (0, 1, 2))):
            return False
        if not all(_J3V2__row_is(rows, node, operation='sum',
                                  arguments=[product, integral])
                   for node, product, integral in zip(final, products, integrals)):
            return False
    return True


def _J3V2__c6_source_rows(background):
    try:
        cb = background['prescribed_R100_Rm_C6_bindings']
        rows = background['exact_graph_nodes']
        first, second = cb['charts']['first_switch'], cb['charts']['second_switch']
        r1 = cb['endpoints']['R1']
        fixed = {
            73190: {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'},
            75949: {'operation': 'bound_variable', 'name': 'true_prescribed_second_switch_t'},
            85905: {'operation': 'product', 'arguments': [31709, 85904]},
            88600: {'operation': 'product', 'arguments': [31709, 88599]},
            33572: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5},
            69245: {'operation': 'positive_quotient', 'numerator': 1,
                    'denominator': 33572,
                    'source_positive_certificate': 'same strictly positive Pstar'},
            25272: {'operation': 'sum', 'arguments': [1, 25271]},
        }
        fixed_ok = all(_J3V2__row_is(rows, number, **expected)
                       for number, expected in fixed.items())
        alias_ok = all(_J3V2__row_is(
            rows, endpoint, operation='function_substitution',
            expression=expression, variable=73190, value=1,
            Z_independent_substitution=True)
            for endpoint, expression in (
                (88430, 86947), (88431, 86950), (88432, 86954),
                (88437, 87154), (88438, 87155), (88439, 87156),
                (88444, 87196), (88445, 87197), (88446, 87198)))
        native_alias_ok = all(_J3V2__row_is(
            rows, endpoint, operation='function_substitution',
            expression=expression, variable=73190, value=1,
            Z_independent_substitution=True)
            for endpoint, expression in (
                (88458, 87352), (88459, 87353), (88460, 87354),
                (88465, 87392), (88466, 87393), (88467, 87394),
                (88472, 87427), (88473, 87428), (88474, 87429),
                (88479, 87462), (88480, 87463), (88481, 87464),
                (88486, 87497), (88487, 87498), (88488, 87499),
                (88493, 87525), (88494, 87526), (88495, 87527)))
        first_integrals = all(_J3V2__integral_ok(
            _J3V2__row(rows, node), 'true_first_switch_history_t', 73190,
            order, rate=1, history_name='M', bridge_chart='true_first_switch')
            for node, order in zip((87385, 87386, 87387), (0, 1, 2)))
        first_m = all(_J3V2__row_is(rows, node, operation='sum',
                                     arguments=[initial, integral])
                      for node, initial, integral in zip(
                          (87392, 87393, 87394), (87378, 87379, 87380),
                          (87385, 87386, 87387)))
        second_logratio = all(_J3V2__row_is(
            rows, node, operation='definite_integral',
            variable='true_second_switch_angular_t', lower=0, upper=75949,
            exact_function_integral=True, numerical_value_not_installed=True,
            endpoint_independent_of_Z=True,
            ordinary_slow_Z_derivative_order=order,
            true_core_inlet_complete_history=True,
            numerical_cap_not_defining_integrand=True,
            original_second_switch_full_current_radius_direction_and_power=True)
            for node, order in zip((89610, 89611, 89612), (0, 1, 2)))
        second_integrals = all(_J3V2__integral_ok(
            _J3V2__row(rows, node), 'true_second_switch_history_t', 75949,
            order, rate=1, history_name='M', bridge_chart='true_second_switch')
            for node, order in zip((89903, 89904, 89905), (0, 1, 2)))
        second_m = all(_J3V2__row_is(rows, node, operation='sum',
                                     arguments=[decay, integral])
                      for node, decay, integral in zip(
                          (89910, 89911, 89912), (89896, 89897, 89898),
                          (89903, 89904, 89905)))
        incoming_m = all(_J3V2__row_is(rows, node, operation='product',
                                       arguments=[incoming, 89878])
                         for node, incoming in zip(
                             (89896, 89897, 89898), (88465, 88466, 88467)))
        second_phi = all(_J3V2__row_is(rows, node, operation=operation,
                                       arguments=arguments)
                         for node, operation, arguments in (
                             (89644, 'product', [88430, 89617]),
                             (89645, 'product', [88430, 89618]),
                             (89646, 'product', [88431, 89617]),
                             (89647, 'sum', [89645, 89646]),
                             (89648, 'product', [88430, 89621]),
                             (89649, 'product', [48, 88431, 89618]),
                             (89650, 'product', [88432, 89617]),
                             (89651, 'sum', [89648, 89649, 89650])))
        m_scale = all(_J3V2__row_is(rows, node, operation='product',
                                    arguments=[moment, 69245])
                      for node, moment in zip(
                          (87591, 87592, 87593), (87392, 87393, 87394)))
        m_scale = m_scale and all(_J3V2__row_is(
            rows, node, operation='product', arguments=[moment, 69245])
            for node, moment in zip(
                (90109, 90110, 90111), (89910, 89911, 89912)))
        generic = (
            all(_J3V2__row_is(rows, node, operation='product',
                              arguments=[source, 28]) for node, source in (
                                  (87632, 87598), (87633, 87601),
                                  (87634, 87605), (90150, 90116),
                                  (90151, 90119), (90152, 90123)))
            and all(_J3V2__row_is(rows, node, operation='product',
                                  arguments=[source, 69387]) for node, source in (
                                      (87673, 87639), (87674, 87642),
                                      (87675, 87646), (90191, 90157),
                                      (90192, 90160), (90193, 90164)))
            and all(_J3V2__row_is(rows, node, operation='negative', argument=argument)
                    for node, argument in ((87721, 87687), (87722, 87690),
                                           (87723, 87694), (90239, 90205),
                                           (90240, 90208), (90241, 90212)))
            and all(_J3V2__row_is(rows, node, operation='sum',
                                  arguments=[source, negative]) for node, source, negative in (
                                      (87728, 87680, 87721),
                                      (87729, 87681, 87722),
                                      (87730, 87682, 87723),
                                      (90246, 90198, 90239),
                                      (90247, 90199, 90240),
                                      (90248, 90200, 90241)))
            and all(_J3V2__row_is(rows, node, operation=operation,
                                  arguments=arguments) for node, operation, arguments in (
                                      (87735, 'product', [87543, 87525]),
                                      (87738, 'sum', [87736, 87737]),
                                      (87742, 'sum', [87739, 87740, 87741]),
                                      (90253, 'product', [90061, 90043]),
                                      (90256, 'sum', [90254, 90255]),
                                      (90260, 'sum', [90257, 90258, 90259]))))
        # Keep the complete common F0/Pstar/radius recipe explicit.  The
        # compact row map below carries the same first/second switch proof.
        scales = all(_J3V2__row_is(rows, node, operation=operation,
                                   **kwargs) for node, operation, kwargs in (
            (87532, 'product', {'arguments': [48, 85905]}),
            (87533, 'analytic_unary', {'name': 'sqrt', 'argument': 87532}),
            (87534, 'product', {'arguments': [87533, 69245]}),
            (87535, 'product', {'arguments': [35657, 87534]}),
            (87536, 'product', {'arguments': [35658, 87534]}),
            (87537, 'product', {'arguments': [35661, 87534]}),
            (87542, 'product', {'arguments': [85905, 69245, 69245]}),
            (87543, 'product', {'arguments': [35693, 87542]}),
            (87544, 'product', {'arguments': [35696, 87542]}),
            (87545, 'product', {'arguments': [35700, 87542]}),
            (90050, 'product', {'arguments': [48, 88600]}),
            (90051, 'analytic_unary', {'name': 'sqrt', 'argument': 90050}),
            (90052, 'product', {'arguments': [90051, 69245]}),
            (90053, 'product', {'arguments': [35657, 90052]}),
            (90054, 'product', {'arguments': [35658, 90052]}),
            (90055, 'product', {'arguments': [35661, 90052]}),
            (90060, 'product', {'arguments': [88600, 69245, 69245]}),
            (90061, 'product', {'arguments': [35693, 90060]}),
            (90062, 'product', {'arguments': [35696, 90060]}),
            (90063, 'product', {'arguments': [35700, 90060]})))
        return {
            'C6.coordinate_radius_Pstar_L_rows': fixed_ok,
            'C6.actual_endpoint_alias_rows': alias_ok,
            'C6.native_endpoint_alias_rows': native_alias_ok,
            'C6.first_M_history_rows': first_integrals and first_m,
            'C6.second_M_history_rows': second_logratio and second_integrals and second_m and incoming_m,
            'C6.second_phi_transport_rows': second_phi,
            'C6.inverse_Pstar_m_rows': m_scale,
            'C6.native_six_history_rows_Z2': _J3V2__c6_native_history_rows(rows),
            'C6.generic_h_k_e_p_rows': generic,
            'C6.common_scale_rows': scales,
        }
    except (KeyError, IndexError, TypeError, AttributeError):
        return {name: False for name in (
            'C6.coordinate_radius_Pstar_L_rows', 'C6.actual_endpoint_alias_rows',
            'C6.native_endpoint_alias_rows', 'C6.first_M_history_rows',
            'C6.second_M_history_rows', 'C6.second_phi_transport_rows',
            'C6.inverse_Pstar_m_rows', 'C6.native_six_history_rows_Z2',
            'C6.generic_h_k_e_p_rows', 'C6.common_scale_rows')}

J3V2__RU_ENDPOINT_PAIRS = {
    'c6_R1_and_second_lower': [[88437, 88438]],
    'c5_first_V_N': [[106932, 106933]],
    'c5_second_V_N': [[109340, 109341]],
    'c3_first_physical_U': [[148085, 148086]],
    'c3_second_physical_U': [[154247, 154248]],
    'common_Pstar': 33572,
    'common_L': 25272,
    'common_N': 18,
}
J3V2__MQ_SOURCE_PAIRS = {
    'physical_m_through_Z2': [[148097, 154259], [148098, 154260], [148099, 154261]],
    'physical_U_premise': [[148085, 154247], [148086, 154248]],
    'physical_Q': [[148180, 154342], [148184, 154346]],
    'physical_Q_numerators': [[148175, 154337], [148183, 154345]],
}
J3V2__SOURCE_HASHES = {
    'parent_report_sha256': RECORD_HASHES[PARENT],
    'c6_creator_sha256': RECORD_HASHES[C6_CREATOR],
    'bridge_creator_sha256': RECORD_HASHES[BRIDGE_CREATOR],
    'loop_creator_sha256': RECORD_HASHES[LOOP_CREATOR],
    'mod_creator_sha256': RECORD_HASHES[MOD_CREATOR],
    'c6_raw_archive_sha256': RECORD_HASHES[BACKGROUND],
    'loop_raw_archive_sha256': RECORD_HASHES[LOOP],
    'mod_raw_archive_sha256': RECORD_HASHES[MOD],
}


def _J3V2__finite_N_rows(loop_rows, endpoint_aliases):
    """Check every finite-N h/k/e/p/m row through ordinary Z order 2."""
    fields_ok = True
    first_history_rows, second_history_rows = {}, {}
    for field, spec in J3V2__CORRECTION_SPECS.items():
        first_history_rows[field] = [spec['first_history'] + i for i in range(3)]
        second_history_rows[field] = [spec['second_history'] + i for i in range(3)]
        for order in range(3):
            fi = spec['first_integral'] + order
            si = spec['second_integral'] + order
            frow, srow = _J3V2__row(loop_rows, fi), _J3V2__row(loop_rows, si)
            fields_ok = (fields_ok
                and _J3V2__integral_ok(
                    frow, 'true_loop_C5_transport_first_switch', 73190,
                    order, finite=True)
                and _J3V2__integral_ok(
                    srow, 'true_loop_C5_transport_second_switch', 75949,
                    order, finite=True)
                and isinstance(frow, dict) and frow.get('original_rate') == spec['rate']
                and isinstance(frow, dict) and frow.get('original_chart') == 'first_switch'
                and isinstance(srow, dict) and srow.get('original_rate') == spec['rate']
                and isinstance(srow, dict) and srow.get('original_chart') == 'second_switch')
            first_head = (spec['first_decay'] + order
                          if spec['rate'] != '0' else spec['first_incoming'][order])
            second_head = (spec['second_decay'] + order
                           if spec['rate'] != '0' else spec['second_incoming'][order])
            fields_ok = (fields_ok
                and _J3V2__row_is(loop_rows, spec['first_history'] + order,
                                  operation='sum',
                                  arguments=[first_head, spec['first_integral'] + order])
                and _J3V2__row_is(loop_rows, spec['second_history'] + order,
                                  operation='sum',
                                  arguments=[second_head, spec['second_integral'] + order]))
            if spec['rate'] != '0':
                fields_ok = (fields_ok
                    and _J3V2__row_is(loop_rows, spec['first_decay'] + order,
                                      operation='product',
                                      arguments=[spec['first_incoming'][order],
                                                  spec['first_memory']])
                    and _J3V2__row_is(loop_rows, spec['second_decay'] + order,
                                      operation='product',
                                      arguments=[spec['second_incoming'][order],
                                                  spec['second_memory']]))
            else:
                fields_ok = (fields_ok
                    and spec['first_decay'] == spec['first_history']
                    and spec['second_decay'] == spec['second_history'])
        for order, node in enumerate(endpoint_aliases[field]):
            fields_ok = (fields_ok and _J3V2__row_is(
                loop_rows, node, operation='function_substitution',
                expression=spec['first_history'] + order,
                variable=73190, value=1, Z_independent_substitution=True))
    return fields_ok, first_history_rows, second_history_rows


def _J3V2__loop_source_rows(loop, background):
    """Check common N, canonical parameters, phase/inverse, and C5 roots."""
    try:
        lp = loop['prescribed_loop_C5_bindings']
        rows = loop['exact_graph_nodes']
        first, second = lp['charts']['first_switch'], lp['charts']['second_switch']
        endpoint_aliases = {
            'm': [107131, 107132, 107133],
            'h': [107137, 107138, 107139],
            'k': [107143, 107144, 107145],
            'e': [107149, 107150, 107151],
            'p': [107155, 107156, 107157],
        }
        typed = (first.get('coordinate') == 73190 and first.get('radius') == 85905
            and first.get('domain') == [0, 1] and first.get('physical_measure') == 9
            and second.get('coordinate') == 75949 and second.get('radius') == 88600
            and second.get('domain') == [0, 1] and second.get('physical_measure') == 9
            and first.get('same_background_pressure_datum') is True
            and second.get('same_background_pressure_datum') is True
            and first.get('total_Z_phase_derivative_exact_zero') is True
            and second.get('total_Z_phase_derivative_exact_zero') is True
            and all(_J3V2__array(lp['endpoints']['first_switch'][field], values)
                    and J3V2__CORRECTION_SPECS[field]['second_incoming'] == values
                    for field, values in endpoint_aliases.items()))
        roots = {
            'first': {'E': [87550, 87553], 'V': [87584, 87585],
                      'a': [104562, 104566], 'b': [104597, 104601],
                      'p1': [104650, 104654], 'p2': [88423, 88424]},
            'second': {'E': [90068, 90071], 'V': [90102, 90103],
                       'a': [107184, 107188], 'b': [0, 0],
                       'p1': [107237, 107241], 'p2': [90934, 90935]},
        }
        root_ok = all(_J3V2__array(chart.get('roots', {}).get(name), values)
                      for chart, expected in ((first, roots['first']),
                                              (second, roots['second']))
                      for name, values in expected.items())
        fixed = {
            18: {'operation': 'exact_positive_integer_expression',
                 'type': 'DoubleDyadicInteger', 'definition': 'N=2^(2^J)',
                 'exact_positive_integer': True, 'Z_independent': True,
                 'one_common_N_for_all_charts_and_control_denominators': True},
            99411: {'operation': 'positive_quotient', 'numerator': 1,
                    'denominator': 18,
                    'source_positive_certificate': 'same exact positive integer N on every source chart'},
            104535: {'operation': 'sum', 'arguments': [104534, 97298]},
            104536: {'operation': 'product', 'arguments': [18, 104535]},
            104537: {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 104536},
            107163: {'operation': 'sum', 'arguments': [107162, 97298]},
            107164: {'operation': 'product', 'arguments': [18, 107163]},
            107165: {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 107164},
            87584: {'operation': 'product', 'arguments': [87154, 69245]},
            87585: {'operation': 'product', 'arguments': [87155, 69245]},
            90102: {'operation': 'product', 'arguments': [88437, 69245]},
            90103: {'operation': 'product', 'arguments': [88438, 69245]},
            106691: {'operation': 'product', 'arguments': [106437, 99411]},
            106692: {'operation': 'product', 'arguments': [106438, 99411]},
            106932: {'operation': 'sum', 'arguments': [87584, 106691]},
            106933: {'operation': 'sum', 'arguments': [87585, 106692]},
            108857: {'operation': 'original_lazy_flat_branch', 'active_body': 108851,
                     'Delta': 107272, 'eta': 21, 'flat_value': 0,
                     'flat_predicate': 107273, 'original_collar_coordinate': None,
                     'original_collar_upper': None,
                     'active_body_not_evaluated_on_flat_branch': True,
                     'ordinary_slow_Z_derivative_order': 0,
                     'same_original_flat_extension_through_order5': True},
            108858: {'operation': 'original_lazy_flat_branch', 'active_body': 108852,
                     'Delta': 107272, 'eta': 21, 'flat_value': 0,
                     'flat_predicate': 107273, 'original_collar_coordinate': None,
                     'original_collar_upper': None,
                     'active_body_not_evaluated_on_flat_branch': True,
                     'ordinary_slow_Z_derivative_order': 1,
                     'same_original_flat_extension_through_order5': True},
            109099: {'operation': 'product', 'arguments': [108857, 99411]},
            109100: {'operation': 'product', 'arguments': [108858, 99411]},
            109340: {'operation': 'sum', 'arguments': [90102, 109099]},
            109341: {'operation': 'sum', 'arguments': [90103, 109100]},
            106437: {'operation': 'original_lazy_flat_branch', 'active_body': 106431,
                     'Delta': 104793, 'eta': 21, 'flat_value': 0,
                     'flat_predicate': 104794, 'original_collar_coordinate': None,
                     'original_collar_upper': None,
                     'active_body_not_evaluated_on_flat_branch': True,
                     'ordinary_slow_Z_derivative_order': 0,
                     'same_original_flat_extension_through_order5': True},
            106438: {'operation': 'original_lazy_flat_branch', 'active_body': 106432,
                     'Delta': 104793, 'eta': 21, 'flat_value': 0,
                     'flat_predicate': 104794, 'original_collar_coordinate': None,
                     'original_collar_upper': None,
                     'active_body_not_evaluated_on_flat_branch': True,
                     'ordinary_slow_Z_derivative_order': 1,
                     'same_original_flat_extension_through_order5': True},
        }
        fixed_ok = all(_J3V2__row_is(rows, number, **expected)
                       for number, expected in fixed.items())
        p2_ok = all(_J3V2__row_is(rows, number, **expected) for number, expected in {
            88388: {'operation': 'positive_quotient', 'numerator': 88382,
                    'denominator': 87550,
                    'source_positive_certificate': 'same actual positive core inlet swirl E'},
            88392: {'operation': 'positive_quotient', 'numerator': 88391,
                    'denominator': 87550,
                    'source_positive_certificate': 'same actual positive core inlet swirl E'},
            88423: {'operation': 'product', 'arguments': [88388, 85905]},
            88424: {'operation': 'product', 'arguments': [88392, 85905]},
            90899: {'operation': 'positive_quotient', 'numerator': 90893,
                    'denominator': 90068,
                    'source_positive_certificate': 'same actual positive core inlet swirl E'},
            90903: {'operation': 'positive_quotient', 'numerator': 90902,
                    'denominator': 90068,
                    'source_positive_certificate': 'same actual positive core inlet swirl E'},
            90934: {'operation': 'product', 'arguments': [90899, 88600]},
            90935: {'operation': 'product', 'arguments': [90903, 88600]},
        }.items())
        canonical = lambda number, name, source_attribute: (
            isinstance(_J3V2__row(rows, number), dict)
            and _J3V2__row(rows, number).get('operation') == 'current_original_source_parameter'
            and _J3V2__row(rows, number).get('name') == name
            and _J3V2__row(rows, number).get('source_attribute') == source_attribute
            and _J3V2__row(rows, number).get('quantity_is_exact_original_parameter_not_range_endpoint') is True
            and _J3V2__row(rows, number).get('Z_independent') is True
            and _J3V2__row(rows, number).get('source_family') == loop.get('source_family')
            and isinstance(_J3V2__row(rows, number).get('canonical_parameter_proof'), dict)
            and _J3V2__row(rows, number)['canonical_parameter_proof'].get('source_family') == loop.get('source_family')
            and _J3V2__row(rows, number)['canonical_parameter_proof'].get('source_provider_types') == [
                'WholeZO3RcFiniteN', 'WholeZSwitchFiniteN',
                'WholeZLongRmFiniteN', 'WholeZRmPatchFiniteN']
            and _J3V2__row(rows, number)['canonical_parameter_proof'].get('exact_original_parameter_tuples_equal') is True
            and _J3V2__row(rows, number)['canonical_parameter_proof'].get('parameter_ranges_not_source_field_values') is True)
        canonical_ok = canonical(21, 'eta', 'eta_log') and canonical(22, 'd_star', 'dstar_log')
        inverse_ok = all(_J3V2__row_is(rows, number,
            operation='original_monotone_phase_inverse', phase_function=phase,
            angle_variable=angle, target_phase=target, angle_lower=0,
            angle_upper=117, phase_endpoints=[0, 1],
            derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction',
            flat_inverse='2*pi*phase', Delta=delta, eta=21,
            branch_inverse_is_original_not_support_cap=True,
            true_prescribed_source_C5_implicit_inverse=True)
            for number, phase, angle, target, delta in (
                (106017, 105891, 'true_loop_C5_angle_first_switch', 104537, 104793),
                (108455, 108329, 'true_loop_C5_angle_second_switch', 107165, 107272)))
        phase_ok = (lp.get('N') == 18 and lp.get('phase_origin') == 24954
                    and fixed_ok and inverse_ok)
        finite_ok, first_hist, second_hist = _J3V2__finite_N_rows(
            rows, endpoint_aliases)
        return {
            'LOOP.typed_switch_bindings': typed,
            'LOOP.common_N_and_phase_origin': phase_ok,
            'LOOP.root_operands_and_p2': root_ok and p2_ok,
            'LOOP.canonical_eta_dstar': canonical_ok,
            'LOOP.phase_inverse_and_lazy_flat': inverse_ok and fixed_ok,
            'LOOP.finite_N_h_k_e_p_m_rows_Z2': finite_ok,
            'LOOP.finite_N_history_rows': finite_ok,
        }
    except (KeyError, IndexError, TypeError, AttributeError):
        return {
            'LOOP.typed_switch_bindings': False,
            'LOOP.common_N_and_phase_origin': False,
            'LOOP.root_operands_and_p2': False,
            'LOOP.canonical_eta_dstar': False,
            'LOOP.phase_inverse_and_lazy_flat': False,
            'LOOP.finite_N_h_k_e_p_m_rows_Z2': False,
            'LOOP.finite_N_history_rows': False,
        }


def _J3V2__mod_source_rows(mod):
    try:
        mb = mod['modified_Q4_pressure_C3_bindings']
        rows = mod['exact_graph_nodes']
        first, second = mb['charts']['first_switch']['physical'], mb['charts']['second_switch']['physical']
        bindings = (
            _J3V2__array(first.get('U'), [148085, 148086])
            and _J3V2__array(first.get('m'), [148097, 148098, 148099])
            and _J3V2__array(first.get('Q'), [148180, 148184])
            and _J3V2__array(second.get('U'), [154247, 154248])
            and _J3V2__array(second.get('m'), [154259, 154260, 154261])
            and _J3V2__array(second.get('Q'), [154342, 154346]))
        physical_u = all(_J3V2__row_is(rows, node, operation='product',
                                        arguments=[source, 33572])
                         for node, source in ((148085, 106932),
                                              (148086, 106933),
                                              (154247, 109340),
                                              (154248, 109341)))
        modified_m = (all(_J3V2__row_is(rows, node, operation='sum',
                                         arguments=[base, history])
                          for node, base, history in zip(
                              (107095, 107096, 107097),
                              (87591, 87592, 87593),
                              (106970, 106971, 106972)))
                       and all(_J3V2__row_is(rows, node, operation='sum',
                                             arguments=[base, history])
                               for node, base, history in zip(
                                   (109504, 109505, 109506),
                                   (90109, 90110, 90111),
                                   (109379, 109380, 109381))))
        physical_m = (all(_J3V2__row_is(rows, node, operation='product',
                                        arguments=[source, 33572])
                          for node, source in zip(
                              (148097, 148098, 148099),
                              (107095, 107096, 107097)))
                      and all(_J3V2__row_is(rows, node, operation='product',
                                            arguments=[source, 33572])
                              for node, source in zip(
                                  (154259, 154260, 154261),
                                  (109504, 109505, 109506))))
        numerator = all(_J3V2__row_is(rows, node, operation='sum',
                                      arguments=args) for node, args in (
            (148175, [148127, 148170]), (148183, [148176, 148182]),
            (154337, [154289, 154332]), (154345, [154338, 154344])))
        qroots = all(_J3V2__row_is(
            rows, node, operation='positive_quotient', numerator=numerator_node,
            denominator=25272,
            source_positive_certificate='same original L=1-delta*Z^2')
            for node, numerator_node in (
                (148180, 148175), (148184, 148183),
                (154342, 154337), (154346, 154345)))
        common_l = _J3V2__row_is(rows, 25272, operation='sum', arguments=[1, 25271])
        return {
            'MOD.typed_physical_m_U_Q_bindings': bindings,
            'MOD.physical_U_Pstar_rows': physical_u,
            'MOD.finite_N_m_memory_rows': modified_m,
            'MOD.physical_m_Pstar_once_rows': physical_m,
            'MOD.Q_numerator_rows': numerator,
            'MOD.Q_positive_L_rows': qroots,
            'MOD.common_L_recipe': common_l,
        }
    except (KeyError, IndexError, TypeError, AttributeError):
        return {name: False for name in (
            'MOD.typed_physical_m_U_Q_bindings', 'MOD.physical_U_Pstar_rows',
            'MOD.finite_N_m_memory_rows', 'MOD.physical_m_Pstar_once_rows',
            'MOD.Q_numerator_rows', 'MOD.Q_positive_L_rows',
            'MOD.common_L_recipe')}


def _J3V2__full_suffix_guards(background, loop, mod, prefix_certificate):
    """Return the complete unexecuted RU/MQ suffix predicate set."""
    p = {
        'premise.boolean_prefix_certificate': _J3V2__certificate_ok(prefix_certificate),
        'premise.fresh_RU_boundary': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('fresh_certificate_before_RU') is True
            and prefix_certificate.get('new_R_U_C1_checked_and_failed') is True),
        'scope.ordinary_Z_orders': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('ordinary_Z_metadata', {}).get('ordinary_Z_orders') == [0, 1]),
        'scope.original_Z2_orders': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('ordinary_Z_metadata', {}).get('original_Z2_premise_orders') == [0, 1, 2]),
        'scope.owned_cell': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('ordinary_Z_metadata', {}).get('owned_Z_cell') == {'lower': '0', 'upper': '1/2'}),
    }
    c6 = _J3V2__c6_source_rows(background)
    loop_checks = _J3V2__loop_source_rows(loop, background)
    mod_checks = _J3V2__mod_source_rows(mod)
    p.update(c6)
    p.update(loop_checks)
    p.update({key: value for key, value in mod_checks.items()
              if key != 'MOD.finite_N_memory_rows'})
    p['MOD.finite_N_memory_rows'] = bool(
        loop_checks.get('LOOP.finite_N_h_k_e_p_m_rows_Z2')
        and mod_checks.get('MOD.finite_N_memory_rows'))
    try:
        cb = background['prescribed_R100_Rm_C6_bindings']
        lp = loop['prescribed_loop_C5_bindings']
        mp = mod['modified_Q4_pressure_C3_bindings']
        first_c6, second_c6 = cb['charts']['first_switch'], cb['charts']['second_switch']
        r1 = cb['endpoints']['R1']
        first_loop, second_loop = lp['charts']['first_switch'], lp['charts']['second_switch']
        first_mod, second_mod = mp['charts']['first_switch'], mp['charts']['second_switch']
        p['source_join.endpoint_semantic_law'] = (
            _J3V2__array(first_c6['actual']['phi'], [86947, 86950])
            and _J3V2__array(first_c6['actual']['U'], [87154, 87155])
            and _J3V2__array(first_c6['actual']['phi_y'], [87196, 87197])
            and _J3V2__array(first_c6['actual']['U_y'], [87126, 87127])
            and _J3V2__array(first_c6['fields']['E'], [87550, 87553])
            and _J3V2__array(first_c6['fields']['V'], [87584, 87585])
            and _J3V2__array(first_c6['radial']['E_y'], [87790, 87793])
            and _J3V2__array(first_c6['radial']['V_y'], [87824, 87825])
            and _J3V2__array(first_c6['generic']['p2'], [88423, 88424])
            and first_c6['actual']['U_y'][:2] == [87126, 87127]
            and second_c6['actual']['U_y'] == [0, 0, 0, 0, 0, 0, 0]
            and second_c6['radial']['V_y'] == [0, 0, 0, 0, 0, 0, 0]
            and r1['actual']['U'][:2] == [88437, 88438]
            and second_c6['actual']['U'][:2] == r1['actual']['U'][:2]
            and r1['actual']['phi'][:2] == [88430, 88431]
            and r1['actual']['phi_y'][:2] == [88444, 88445]
            and r1['actual']['U_y'][:2] == [88451, 88452]
            and r1['fields']['E'][:2] == [88500, 88501]
            and r1['fields']['V'][:2] == [88507, 88508]
            and r1['generic']['p2'][:2] == [88592, 88593]
            and second_c6['actual']['phi'][:2] == [89644, 89647]
            and second_c6['actual']['phi_y'][:2] == [89687, 89690]
            and second_c6['fields']['E'][:2] == [90068, 90071]
            and second_c6['fields']['V'][:2] == [90102, 90103]
            and second_c6['radial']['E_y'][:2] == [90308, 90311]
            and all(_J3V2__row_is(background['exact_graph_nodes'], number,
                operation='function_substitution', expression=expression,
                variable=73190, value=1, Z_independent_substitution=True)
                for number, expression in (
                    (88430, 86947), (88431, 86950),
                    (88437, 87154), (88438, 87155),
                    (88444, 87196), (88445, 87197),
                    (88451, 87126), (88452, 87127),
                    (88500, 87550), (88501, 87553),
                    (88507, 87584), (88508, 87585),
                    (88592, 88423), (88593, 88424)))
            and _J3V2__array(first_loop['roots']['E'], [87550, 87553])
            and _J3V2__array(first_loop['roots']['V'], [87584, 87585])
            and _J3V2__array(second_loop['roots']['E'], [90068, 90071])
            and _J3V2__array(second_loop['roots']['V'], [90102, 90103])
            and _J3V2__array(first_mod['physical']['U'], [148085, 148086])
            and _J3V2__array(second_mod['physical']['U'], [154247, 154248]))
        p['source_join.common_comparison_radius'] = (
            first_c6.get('radius') == 85905 and second_c6.get('radius') == 88600
            and first_c6.get('known_comparison_direction', {}).get('D_over_R', [])[:2]
                == [86464, 86468]
            and second_c6.get('known_comparison_direction', {}).get('D_over_R', [])[:2]
                == [89159, 89163]
            and first_c6.get('known_comparison_direction', {}).get('drive_hydro', [])[:2]
                == [86650, 86654]
            and second_c6.get('known_comparison_direction', {}).get('drive_hydro', [])[:2]
                == [89345, 89349]
            and first_c6.get('known_comparison_direction', {}).get('drive_pressure', [])[:2]
                == [68391, 68395]
            and second_c6.get('known_comparison_direction', {}).get('drive_pressure', [])[:2]
                == [68391, 68395]
            and first_c6.get('known_comparison_direction', {}).get('drive_swirl', [])[:2]
                == [86853, 86857]
            and second_c6.get('known_comparison_direction', {}).get('drive_swirl', [])[:2]
                == [89548, 89552])
        p['source_join.c6_fields_loop_roots_p2'] = (
            first_loop['roots']['p2'] == first_c6['generic']['p2']
            and second_loop['roots']['p2'] == second_c6['generic']['p2'])
        p['source_join.p2_actual_denominator_inputs'] = (
            first_c6['fields']['E'][:2] == [87550, 87553]
            and second_c6['fields']['E'][:2] == [90068, 90071]
            and first_c6['generic']['p2'][:2] == [88423, 88424]
            and second_c6['generic']['p2'][:2] == [90934, 90935]
            and first_loop['roots']['p2'][:2] == first_c6['generic']['p2'][:2]
            and second_loop['roots']['p2'][:2] == second_c6['generic']['p2'][:2]
            and loop_checks.get('LOOP.root_operands_and_p2') is True)
        p['source_join.c5_roots_actual_operand_rows'] = (
            first_c6['radial']['E_y'][:2] == [87790, 87793]
            and first_c6['radial']['V_y'][:2] == [87824, 87825]
            and second_c6['radial']['E_y'][:2] == [90308, 90311]
            and second_c6['radial']['V_y'][:2] == [0, 0]
            and first_loop['roots']['a'][:2] == [104562, 104566]
            and first_loop['roots']['b'][:2] == [104597, 104601]
            and first_loop['roots']['p1'][:2] == [104650, 104654]
            and second_loop['roots']['a'][:2] == [107184, 107188]
            and second_loop['roots']['b'][:2] == [0, 0]
            and second_loop['roots']['p1'][:2] == [107237, 107241])
        p['source_join.physical_m_pairs'] = (
            first_mod['physical']['m'][:3] == [148097, 148098, 148099]
            and second_mod['physical']['m'][:3] == [154259, 154260, 154261])
        p['source_join.physical_Q_pairs'] = (
            first_mod['physical']['Q'][:2] == [148180, 148184]
            and second_mod['physical']['Q'][:2] == [154342, 154346])
        # The C6 delegate/end constructor and the exact flat sigma law are
        # authenticated by creator body needles plus inherited source pins.
        p['source_semantic.sigma_cutoff_law'] = (
            J3V2__SIGMA_CUTOFF_DEFINITION ==
            '0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'
            and J3V2__SIGMA_SOURCE_SHA256 ==
            '38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01')
        p['source_semantic.sigma_archive_branch_recipe'] = all(
            _J3V2__row_is(background['exact_graph_nodes'], number,
                operation='exact_original_source_flat_sigma', argument=argument,
                definition=J3V2__SIGMA_CUTOFF_DEFINITION,
                exact_function_not_its_interval_cap=True)
            for number, argument in ((74428, 73190), (76947, 75949)))
    except (KeyError, IndexError, TypeError, AttributeError):
        for name in ('source_join.endpoint_semantic_law',
                     'source_join.c6_fields_loop_roots_p2',
                     'source_join.physical_m_pairs',
                     'source_join.physical_Q_pairs',
                     'source_join.c5_roots_actual_operand_rows',
                     'source_join.p2_actual_denominator_inputs',
                     'source_join.common_comparison_radius',
                     'source_semantic.sigma_cutoff_law',
                     'source_semantic.sigma_archive_branch_recipe',
                     'source_semantic.c6_delegate_and_end'):
            p[name] = False
    return p


def _J3V2__sigma_source_laws(gamma_text, c6_text, background, inherited_pins):
    """Read two fixed source rows and prove only the exact cutoff endpoint law."""
    predicates = {
        'gamma_source_raw_sha': _hashlib.sha256(gamma_text.encode()).hexdigest() == RECORD_HASHES[GAMMA_CREATOR],
        'c6_source_raw_sha': _hashlib.sha256(c6_text.encode()).hexdigest() == RECORD_HASHES[C6_CREATOR],
        'inherited_gamma_pin': inherited_pins.get(GAMMA_CREATOR) == RECORD_HASHES[GAMMA_CREATOR],
        'inherited_flat_binding_pin': inherited_pins.get('lei_ren_part1_paper_compliant_outer_initial.py') == 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657',
    }
    gamma_fn = next((n for n in _ast.parse(gamma_text).body if isinstance(n, _ast.FunctionDef) and n.name == 'sigma'), None)
    delegate = next((n for n in _ast.parse(c6_text).body if isinstance(n, _ast.FunctionDef) and n.name == 'sigma'), None)
    def body_hash(node):
        return _hashlib.sha256(_ast.dump(_ast.Module(body=node.body, type_ignores=[]), include_attributes=False).encode()).hexdigest() if node is not None else None
    predicates['gamma_exact_definition_body'] = body_hash(gamma_fn) == '90710db1c1ee4bdec06cb4515c082198c1f3fc006d7f1ec58daebd0cddc2476a'
    predicates['c6_same_Gamma_delegate_body'] = body_hash(delegate) == '1f72db601deda38dc3a82db0b52df5d6ee8c7ecad4d6bb983b718c2906217ea4'
    rows = background['exact_graph_nodes']
    binding = dict(AST_sha256='4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c', function='stable_sigma', module='lei_ren_part1_paper_compliant_outer_initial.py', module_sha256='a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657')
    predicates['fixed_first_second_cutoff_rows'] = all(
        _J3V2__row_is(rows, number, operation='exact_original_source_flat_sigma', argument=argument, definition=J3V2__SIGMA_CUTOFF_DEFINITION, source_binding=binding, exact_function_not_its_interval_cap=True)
        for number, argument in ((74428, 73190), (76947, 75949)))
    passed = all(value is True for value in predicates.values())
    return dict(schema='native_J3_exact_cutoff_endpoint_source_laws_v1', passed=passed,
        source_only=True, predicates=predicates, source_sha256=RECORD_HASHES[GAMMA_CREATOR],
        cutoff_definition=J3V2__SIGMA_CUTOFF_DEFINITION,
        fixed_cutoff_rows=[74428, 76947], endpoint_source_values={'at_0': 0, 'at_1': 1} if passed else None,
        ordinary_Z_orders=[0, 1], endpoint_inputs_Z_independent=True,
        derivation=['Use the exact piecewise definition, taking x=0 in x<=0 and x=1 in x>=1.', 'Both substituted arguments are source constants, so their ordinary Z derivatives are zero.'],
        source_constructor_calls=0, old_proof_functions_called=0, numeric_values_computed=False,
        generic_DAG_walks=0, source_body_hashes={'Gamma.sigma': body_hash(gamma_fn), 'C6.sigma': body_hash(delegate)})


def _J3V2__sigma_law_contract(value):
    return (isinstance(value, dict) and value.get('schema') == 'native_J3_exact_cutoff_endpoint_source_laws_v1'
        and value.get('passed') is True and value.get('source_only') is True
        and value.get('source_sha256') == RECORD_HASHES[GAMMA_CREATOR]
        and value.get('endpoint_source_values') == {'at_0': 0, 'at_1': 1}
        and value.get('fixed_cutoff_rows') == [74428, 76947]
        and value.get('cutoff_definition') == J3V2__SIGMA_CUTOFF_DEFINITION
        and isinstance(value.get('predicates'), dict) and bool(value['predicates'])
        and all(flag is True for flag in value['predicates'].values()))


def _J3V2__certificate_ok(certificate):
    """Validate only the projected Boolean prefix, never a runtime state."""
    if not isinstance(certificate, _Mapping):
        return False
    metadata = certificate.get('ordinary_Z_metadata')
    provenance = certificate.get('provenance')
    hashes = certificate.get('source_hashes')
    predicates = certificate.get('predicates')
    return bool(
        certificate.get('schema') == 'native_J3_failed_boolean_prefix_projector_v1'
        and certificate.get('status') == 'boolean_prefix_projected'
        and certificate.get('source_algebra_R_phase_prefix_proved') is True
        and certificate.get('original_native6_generic5_Z2_prefix_proved') is True
        and certificate.get('fresh_certificate_before_RU') is True
        and certificate.get('new_R_U_C1_checked_and_failed') is True
        and isinstance(metadata, _Mapping)
        and metadata.get('ordinary_Z_orders') == [0, 1]
        and metadata.get('original_Z2_premise_orders') == [0, 1, 2]
        and metadata.get('owned_Z_cell') == {'lower': '0', 'upper': '1/2'}
        and metadata.get('ordinary_Z_not_temporal_n') is True
        and isinstance(hashes, _Mapping)
        and hashes.get('implementation') == J3V2__FROZEN_ARTIFACT_SHA256
        and hashes.get('records') == J3V2__FROZEN_RECORD_HASHES
        and hashes.get('decoded_records') == J3V2__FROZEN_DECODED_RECORD_HASHES
        and hashes.get('ast_body') == J3V2__FROZEN_AST_BODY_SHA256
        and isinstance(provenance, _Mapping)
        and provenance.get('static_only') is True
        and provenance.get('runtime_intermediate_objects_reconstructed') is False
        and provenance.get('source_constructor_calls') == 0
        and provenance.get('source_producer_function_calls') == 0
        and provenance.get('old_proof_functions_called') == 0
        and provenance.get('generic_DAG_walks') == 0
        and provenance.get('numeric_evaluations') == 0
        and provenance.get('fragment_was_unexecuted_before_frozen_entry') is True
        and isinstance(predicates, dict) and bool(predicates)
        and all(type(value) is bool and value is True for value in predicates.values()))


def _J3V2__ru_result_contract_ok(ru_result):
    if not isinstance(ru_result, _Mapping):
        return False
    predicates = ru_result.get('predicates')
    return bool(
        ru_result.get('schema') == 'native_J3_R_U_source_transport_guarded_v4'
        and ru_result.get('status') == 'conditional_source_theorem_proved'
        and ru_result.get('source_only') is True
        and ru_result.get('accepted_identity_claim') is False
        and ru_result.get('J3_R_C1_proved') is True
        and ru_result.get('J3_U_C1_proved') is True
        and ru_result.get('ordinary_Z_orders') == [0, 1]
        and ru_result.get('owned_Z_cell') == ['0', '1/2']
        and ru_result.get('exact_U_row_pairs') == J3V2__RU_ENDPOINT_PAIRS
        and ru_result.get('source_hashes') == J3V2__SOURCE_HASHES
        and ru_result.get('numeric_values_computed') is False
        and ru_result.get('source_constructor_calls') == 0
        and ru_result.get('old_proof_functions_called') == 0
        and ru_result.get('source_producer_function_calls') == 0
        and ru_result.get('generic_DAG_walks') == 0
        and ru_result.get('fragment_was_unexecuted_before_frozen_entry') is True
        and isinstance(predicates, dict) and bool(predicates)
        and all(type(value) is bool and value is True for value in predicates.values()))


def _J3V2__prefix_predicates(prefix_certificate):
    metadata = (prefix_certificate.get('ordinary_Z_metadata')
                if isinstance(prefix_certificate, _Mapping) else None)
    if not isinstance(metadata, _Mapping):
        metadata = {}
    return {
        'prefix.boolean_certificate': _J3V2__certificate_ok(prefix_certificate),
        'prefix.source_algebra_R_phase': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('source_algebra_R_phase_prefix_proved') is True),
        'prefix.original_native6_generic5_Z2': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('original_native6_generic5_Z2_prefix_proved') is True),
        'prefix.fresh_certificate_before_RU': (isinstance(prefix_certificate, _Mapping)
            and prefix_certificate.get('fresh_certificate_before_RU') is True),
    }


def J3V2__prove_R_U_transport(loop, mod, background, loop_creator_text,
                               mod_creator_text, c6_creator_text,
                               bridge_creator_text, prefix_certificate,
                               archive_raw_sha256, source_semantic_pins=None):
    """Collect bounded source predicates for R/U without a runtime state."""
    predicates = _J3V2__prefix_predicates(prefix_certificate)
    metadata = (prefix_certificate.get('ordinary_Z_metadata')
                if isinstance(prefix_certificate, _Mapping) else None)
    if not isinstance(metadata, _Mapping):
        metadata = {}
    try:
        predicates.update(_J3V2__creator_ast_guard(c6_creator_text, 'c6'))
        predicates.update(_J3V2__creator_ast_guard(bridge_creator_text, 'bridge'))
        predicates.update(_J3V2__creator_ast_guard(loop_creator_text, 'loop'))
        predicates.update(_J3V2__creator_ast_guard(mod_creator_text, 'mod'))
    except (TypeError, ValueError, KeyError):
        predicates['creator.guard_exception'] = False
    predicates['source_semantic.c6_delegate_and_end'] = all(
        predicates.get(name) is True for name in (
            'c6.body_needles.function.first_switch_fields',
            'c6.body_needles.function.second_switch_fields',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.__init__',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.end'))
    predicates['source_semantic.inherited_sigma_pins'] = (
        isinstance(source_semantic_pins, _Mapping)
        and source_semantic_pins.get('gamma_source_sha256') == J3V2__SIGMA_SOURCE_SHA256
        and source_semantic_pins.get('flat_sigma_source_sha256') == J3V2__FLAT_SIGMA_SOURCE_SHA256
        and _J3V2__sigma_law_contract(source_semantic_pins.get('cutoff_source_laws')))
    try:
        predicates.update(_J3V2__archive_guards(
            background, loop, mod, archive_raw_sha256,
            prefix_certificate=prefix_certificate))
    except (TypeError, ValueError, KeyError, IndexError, AttributeError):
        predicates['archive.guard_exception'] = False
    predicates['source_join.bridge_field_formula_rows'] = all(
        predicates.get(name) is True for name in (
            'bridge.body_needles.function.generic_fields',
            'bridge.body_needles.function.histories',
            'loop.body_needles.function.generic_roots',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.add_chart'))
    predicates['source_join.endpoint_symbolic_derivation'] = all(
        predicates.get(name) is True for name in (
            'source_semantic.c6_delegate_and_end',
            'source_semantic.inherited_sigma_pins',
            'source_semantic.sigma_cutoff_law',
            'source_semantic.sigma_archive_branch_recipe',
            'source_join.endpoint_semantic_law',
            'source_join.common_comparison_radius',
            'source_join.bridge_field_formula_rows',
            'C6.native_six_history_rows_Z2'))
    predicates.update({
        'R_U.endpoint_projection': (J3V2__ENDPOINT_PROJECTION.get('first_switch_upper')
            == 'coordinate node 73190 -> source value 1; radius node 85905'
            and J3V2__ENDPOINT_PROJECTION.get('second_switch_lower')
            == 'coordinate node 75949 -> source value 0; radius node 88600'),
        'R_U.ordinary_Z_orders': metadata.get('ordinary_Z_orders') == [0, 1],
        'R_U.owned_Z_cell': metadata.get('owned_Z_cell') == {'lower': '0', 'upper': '1/2'},
        'R_U.common_Pstar_and_L': (J3V2__RU_ENDPOINT_PAIRS['common_Pstar'] == 33572
            and J3V2__RU_ENDPOINT_PAIRS['common_L'] == 25272
            and predicates.get('rows.common_L_and_Pstar') is True),
        'R_U.common_positive_N_row': (J3V2__RU_ENDPOINT_PAIRS['common_N'] == 18
            and predicates.get('LOOP.common_N_and_phase_origin') is True),
        'R_U.distinct_native_six_generic_five_channels': (
            predicates.get('bindings.native_six_channels') is True
            and predicates.get('bindings.generic_five_channels') is True),
        'R_U.signed_inner_lower_ledger_channels_retained': (
            predicates.get('bindings.native_six_channels') is True
            and predicates.get('C6.native_six_history_rows_Z2') is True),
        'R_U.fresh_source_composition_rows': all(
            predicates.get(name) is True for name in (
                'bindings.C6_actual_and_M_rows', 'bindings.C6_fields_to_loop_roots',
                'bindings.MOD_physical_ports', 'rows.U_source_pairs',
                'rows.Q_source_pairs', 'source_join.endpoint_semantic_law',
                'source_join.c6_fields_loop_roots_p2',
                'source_join.c5_roots_actual_operand_rows',
                'source_join.p2_actual_denominator_inputs',
                'source_join.common_comparison_radius',
                'LOOP.root_operands_and_p2', 'LOOP.canonical_eta_dstar',
                'LOOP.phase_inverse_and_lazy_flat')),
        'R_U.full_unexecuted_suffix_guards': all(
            predicates.get(name) is True for name in (
                'C6.coordinate_radius_Pstar_L_rows',
                'C6.actual_endpoint_alias_rows', 'C6.native_endpoint_alias_rows',
                'C6.first_M_history_rows', 'C6.second_M_history_rows',
                'C6.second_phi_transport_rows', 'C6.inverse_Pstar_m_rows',
                'C6.native_six_history_rows_Z2', 'C6.generic_h_k_e_p_rows',
                'C6.common_scale_rows', 'LOOP.typed_switch_bindings',
                'LOOP.common_N_and_phase_origin', 'LOOP.root_operands_and_p2',
                'LOOP.canonical_eta_dstar', 'LOOP.phase_inverse_and_lazy_flat',
                'LOOP.finite_N_h_k_e_p_m_rows_Z2',
                'MOD.typed_physical_m_U_Q_bindings',
                'MOD.physical_U_Pstar_rows', 'MOD.finite_N_memory_rows',
                'MOD.physical_m_Pstar_once_rows', 'MOD.Q_numerator_rows',
                'MOD.Q_positive_L_rows', 'MOD.common_L_recipe',
                'source_join.c5_roots_actual_operand_rows',
                'source_join.p2_actual_denominator_inputs',
                'source_join.common_comparison_radius',
                'source_join.bridge_field_formula_rows',
                'source_join.endpoint_symbolic_derivation')),
    })
    passed = bool(predicates) and all(type(value) is bool and value is True
                                      for value in predicates.values())
    return {
        'schema': 'native_J3_R_U_source_transport_guarded_v4',
        'status': 'conditional_source_theorem_proved' if passed else 'source_guard_failure',
        'source_only': True,
        'accepted_identity_claim': False,
        'J3_R_C1_proved': passed,
        'J3_U_C1_proved': passed,
        'ordinary_Z_orders': [0, 1],
        'ordinary_Z_not_temporal_n': True,
        'owned_Z_cell': ['0', '1/2'],
        'native_six_signed_channels': ['H', 'M', 'K', 'A', 'B', 'C'],
        'generic_five_channels': ['m', 'h', 'k', 'e', 'p'],
        'endpoint_projection': dict(J3V2__ENDPOINT_PROJECTION),
        'exact_U_row_pairs': dict(J3V2__RU_ENDPOINT_PAIRS),
        'conditional_premises': ['boolean prefix certificate',
                                 'authenticated C6/bridge/C5/C3 source creators',
                                 'authenticated archive rows and common Pstar/L/N'],
        'numeric_values_computed': False,
        'source_constructor_calls': 0,
        'old_proof_functions_called': 0,
        'source_producer_function_calls': 0,
        'generic_DAG_walks': 0,
        'fragment_was_unexecuted_before_frozen_entry': True,
        'source_hashes': dict(J3V2__SOURCE_HASHES),
        'proof_guards': dict(predicates),
        'predicates': dict(predicates),
        'source_derivation': {
            'law': 'Gamma sigma(1)=1 and sigma(0)=0; C6 end carries first endpoint into second; bridge maps phi,U,phi_y,U_y to E,V,E_y,V_y; C5 roots derive p1,p2 at common radius',
            'proof_steps': [
                'exact Gamma flat-sigma cutoff pin and C6 delegate/end body needles',
                'first upper U_y and second lower U_y source endpoint rows',
                'R1 actual U equals second-switch incoming U through Z orders 0 and 1',
                'bridge E,V,E_y,V_y and C5 p1,p2 typed root rows',
                'common radius, comparison, Pstar, and L source rows',
            ],
            'guard_dependencies': {
                key: predicates.get(key) for key in (
                    'source_semantic.inherited_sigma_pins',
                    'source_semantic.c6_delegate_and_end',
                    'source_semantic.sigma_cutoff_law',
                    'source_semantic.sigma_archive_branch_recipe',
                    'source_join.endpoint_semantic_law',
                    'source_join.common_comparison_radius',
                    'source_join.endpoint_symbolic_derivation',
                    'C6.native_six_history_rows_Z2',
                    'LOOP.root_operands_and_p2',
                )},
            'source_hashes': dict(J3V2__SOURCE_HASHES),
            'gamma_source_sha256': J3V2__SIGMA_SOURCE_SHA256,
            'flat_sigma_source_sha256': J3V2__FLAT_SIGMA_SOURCE_SHA256,
        },
        'prefix_certificate_schema': (prefix_certificate.get('schema')
                                      if isinstance(prefix_certificate, _Mapping) else None),
        'nonclaims': {'numeric_containment': False, 'physical_m_Q_C1': False,
                      'accepted_identity_claim': False},
    }


def J3V2__prove_physical_Q_transport(loop, mod, background, loop_creator_text,
                                     mod_creator_text, c6_creator_text,
                                     bridge_creator_text, prefix_certificate,
                                     ru_result, archive_raw_sha256,
                                     source_semantic_pins=None):
    """Collect physical m/Q source predicates from the fresh Boolean premise."""
    predicates = _J3V2__prefix_predicates(prefix_certificate)
    predicates['R_U.contract'] = _J3V2__ru_result_contract_ok(ru_result)
    try:
        predicates.update(_J3V2__creator_ast_guard(c6_creator_text, 'c6'))
        predicates.update(_J3V2__creator_ast_guard(bridge_creator_text, 'bridge'))
        predicates.update(_J3V2__creator_ast_guard(loop_creator_text, 'loop'))
        predicates.update(_J3V2__creator_ast_guard(mod_creator_text, 'mod'))
    except (TypeError, ValueError, KeyError):
        predicates['creator.guard_exception'] = False
    predicates['source_semantic.c6_delegate_and_end'] = all(
        predicates.get(name) is True for name in (
            'c6.body_needles.function.first_switch_fields',
            'c6.body_needles.function.second_switch_fields',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.__init__',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.end'))
    predicates['source_semantic.inherited_sigma_pins'] = (
        isinstance(source_semantic_pins, _Mapping)
        and source_semantic_pins.get('gamma_source_sha256') == J3V2__SIGMA_SOURCE_SHA256
        and source_semantic_pins.get('flat_sigma_source_sha256') == J3V2__FLAT_SIGMA_SOURCE_SHA256
        and _J3V2__sigma_law_contract(source_semantic_pins.get('cutoff_source_laws')))
    try:
        predicates.update(_J3V2__archive_guards(
            background, loop, mod, archive_raw_sha256,
            prefix_certificate=prefix_certificate))
    except (TypeError, ValueError, KeyError, IndexError, AttributeError):
        predicates['archive.guard_exception'] = False
    predicates['source_join.bridge_field_formula_rows'] = all(
        predicates.get(name) is True for name in (
            'bridge.body_needles.function.generic_fields',
            'bridge.body_needles.function.histories',
            'loop.body_needles.function.generic_roots',
            'c6.body_needles.method.CurrentOriginalPrescribedR100RmC6Functions.add_chart'))
    predicates['source_join.endpoint_symbolic_derivation'] = all(
        predicates.get(name) is True for name in (
            'source_semantic.c6_delegate_and_end',
            'source_semantic.inherited_sigma_pins',
            'source_semantic.sigma_cutoff_law',
            'source_semantic.sigma_archive_branch_recipe',
            'source_join.endpoint_semantic_law',
            'source_join.common_comparison_radius',
            'source_join.bridge_field_formula_rows',
            'C6.native_six_history_rows_Z2'))
    predicates.update({
        'm_Q.physical_m_pairs': (predicates.get('MOD.physical_m_Pstar_once_rows') is True
            and predicates.get('source_join.physical_m_pairs') is True),
        'm_Q.physical_U_pairs': (predicates.get('MOD.physical_U_Pstar_rows') is True
            and predicates.get('bindings.MOD_physical_ports') is True),
        'm_Q.physical_Q_pairs': (predicates.get('MOD.Q_positive_L_rows') is True
            and predicates.get('source_join.physical_Q_pairs') is True),
        'm_Q.common_L_and_Pstar': (J3V2__RU_ENDPOINT_PAIRS['common_L'] == 25272
            and J3V2__RU_ENDPOINT_PAIRS['common_Pstar'] == 33572
            and predicates.get('rows.common_L_and_Pstar') is True),
        'm_Q.original_L_recipe': predicates.get('rows.Q_source_pairs') is True,
        'm_Q.finite_N_signed_memory_preserved': (
            predicates.get('LOOP.finite_N_h_k_e_p_m_rows_Z2') is True
            and predicates.get('MOD.finite_N_memory_rows') is True
            and predicates.get('C6.native_six_history_rows_Z2') is True),
        'm_Q.native_six_channels_distinct': (
            predicates.get('bindings.native_six_channels') is True
            and predicates.get('C6.native_six_history_rows_Z2') is True),
        'm_Q.generic_five_channels_distinct': (
            predicates.get('bindings.generic_five_channels') is True
            and predicates.get('C6.generic_h_k_e_p_rows') is True),
        'm_Q.fresh_source_composition_rows': all(
            predicates.get(name) is True for name in (
                'bindings.MOD_physical_ports', 'rows.U_source_pairs',
                'rows.Q_source_pairs', 'rows.common_L_and_Pstar',
                'MOD.typed_physical_m_U_Q_bindings',
                'MOD.finite_N_memory_rows', 'MOD.physical_m_Pstar_once_rows',
                'MOD.Q_numerator_rows', 'MOD.Q_positive_L_rows',
                'MOD.common_L_recipe')),
        'm_Q.full_unexecuted_suffix_guards': all(
            predicates.get(name) is True for name in (
                'LOOP.finite_N_h_k_e_p_m_rows_Z2',
                'LOOP.finite_N_history_rows', 'MOD.finite_N_memory_rows',
                'MOD.physical_m_Pstar_once_rows', 'MOD.Q_numerator_rows',
                'MOD.Q_positive_L_rows', 'MOD.common_L_recipe',
                'source_join.physical_m_pairs', 'source_join.physical_Q_pairs',
                'source_join.endpoint_symbolic_derivation')),
    })
    physical_m = bool(predicates.get('R_U.contract')
                      and predicates.get('m_Q.physical_m_pairs'))
    finite_m = bool(physical_m and predicates.get('m_Q.finite_N_signed_memory_preserved')
                    and predicates.get('loop.ast_parse', False))
    passed = bool(predicates) and all(type(value) is bool and value is True
                                      for value in predicates.values())
    return {
        'schema': 'native_J3_m_Q_source_fragment_v2_guarded',
        'status': 'J3_Q_source_guard_ready' if passed else 'source_guard_failure',
        'source_only': True,
        'accepted_identity_claim': False,
        'finite_N_m_Z2_proved': finite_m and passed,
        'physical_m_Z2_proved': physical_m and passed,
        'J3_Q_C1_proved': passed,
        'ordinary_Z_orders': [0, 1],
        'm_ordinary_Z_orders': [0, 1, 2],
        'ordinary_Z_not_temporal_n': True,
        'owned_Z_cell': ['0', '1/2'],
        'native_six_signed_channels': ['H', 'M', 'K', 'A', 'B', 'C'],
        'generic_five_channels': ['m', 'h', 'k', 'e', 'p'],
        'Q_C1_conditional_on': ['boolean_prefix_certificate', 'J3_R_C1_proved',
                                'J3_U_C1_proved'],
        'source_pairs': dict(J3V2__MQ_SOURCE_PAIRS),
        'common_source_factors': {'Pstar': 33572, 'inverse_Pstar': 69245,
                                  'L': 25272,
                                  'L_recipe': {'operation': 'sum', 'arguments': [1, 25271]}},
        'proof_guards': dict(predicates),
        'predicates': dict(predicates),
        'source_derivation': {
            'law': 'finite-N signed memory rows feed physical m once through Pstar; Q numerators feed positive-L quotient at the same C5 endpoints',
            'proof_steps': [
                'all finite-N h,k,e,p,m correction integrals carry exact Z orders 0,1,2 and predecessor memory',
                'modified m rows add each finite-N history before the single physical Pstar scale',
                'physical U rows preserve the R/U source pairs at both switches',
                'Q numerator rows and the common L=1-delta*Z^2 quotient are exact source operations',
            ],
            'guard_dependencies': {
                key: predicates.get(key) for key in (
                    'R_U.contract', 'source_join.endpoint_symbolic_derivation',
                    'LOOP.finite_N_h_k_e_p_m_rows_Z2',
                    'MOD.finite_N_memory_rows',
                    'MOD.physical_m_Pstar_once_rows',
                    'MOD.Q_numerator_rows', 'MOD.Q_positive_L_rows',
                    'm_Q.full_unexecuted_suffix_guards',
                )},
            'source_hashes': dict(J3V2__SOURCE_HASHES),
        },
        'source_hashes': dict(J3V2__SOURCE_HASHES),
        'numeric_values_computed': False,
        'source_constructor_calls': 0,
        'source_producer_function_calls': 0,
        'old_proof_functions_called': 0,
        'generic_DAG_walks': 0,
        'fragment_was_unexecuted_before_frozen_entry': True,
        'nonclaims': {'Q_y_C1_identity': False, 'weighted_J3_jump_zero': False,
                      'numeric_containment': False, 'accepted_identity_claim': False},
    }


def _J3V2__failed_wrapper_guards(report, evidence, evidence_source_text):
    """Return Boolean authentication facts for the saved failed attempt."""
    if not isinstance(report, _Mapping) or not isinstance(evidence, _Mapping):
        return {'wrapper.input_shape': False}
    p = {}
    p['report.schema'] = (isinstance(report, _Mapping)
                          and report.get('schema') ==
                          'original_native_J3_endpoint_transport_failure_stage_v1')
    p['report.scope'] = (report.get('component') == 'endpoint_transport'
                         and report.get('source_only') is True
                         and report.get('accepted_identity') is False
                         and report.get('accepted_identity_flag') is False
                         and report.get('component_result_is_projection') is False
                         and report.get('component_result_status') ==
                         'full_failed_result_saved_once')
    p['report.failure_stage'] = report.get('failure_stage') == 'new_R_U_C1'
    p['report.parent'] = (report.get('parent_report') == PARENT
                          and report.get('parent_report_path') ==
                          'experiments/root_st073/' + PARENT
                          and report.get('parent_report_sha256') == RECORD_HASHES[PARENT])
    p['report.attempt_counters'] = (report.get('source_pin_count') == 1857
        and report.get('checkpoint_count') == 3
        and report.get('new_current_checkpoint_count') == 3
        and report.get('total_checkpoint_count') == 136
        and report.get('preserved_failed_artifact_count') == 5
        and report.get('preserved_failed_artifact_count_union') == 20
        and report.get('source_consumer_attempted') is True
        and report.get('source_consumer_attempt_count') == 1
        and report.get('accepted_source_consumer_calls') == 0
        and report.get('accepted_ancestor_operator_calls') == 0
        and report.get('old_proof_functions_called') == 0
        and report.get('old_operators_replayed') is False
        and report.get('old_algebra_replayed') is False)
    p['report.failed_result'] = (isinstance(report.get('component_result'), dict)
                                 and _J3V2__failure_ok(report['component_result']))
    pointer = {'path': FAILURE_EVIDENCE, 'sha256': RECORD_HASHES[FAILURE_EVIDENCE],
               'bytes': 566962, 'encoding': 'json', 'raw_bytes_preserved': True,
               'artifact_inventory_selector': '$.artifact_inventory',
               'result_selector': '$.full_result'}
    inventory_pointer = {'selector': '$.failure_evidence.artifact_inventory',
                         'path': FAILURE_EVIDENCE,
                         'sha256': RECORD_HASHES[FAILURE_EVIDENCE],
                         'raw_bytes_lossless': True, 'artifact_count': 5}
    p['report.evidence_pointer'] = report.get('failure_evidence') == pointer
    p['report.inventory_pointer'] = report.get('failure_artifact_inventory') == inventory_pointer
    p['report.false_flags'] = report.get('false_flags') == dict.fromkeys(FALSE_FLAGS, False)
    p['report.input_summary'] = False
    summary = report.get('component_input_summary')
    if isinstance(summary, dict):
        expected_new = {
            'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J2_B_source_closure.json': RECORD_HASHES[PARENT],
            'scripts/run_original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['runner'],
            'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['kernel'],
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_source.json': J3V2__FROZEN_ARTIFACT_SHA256['source'],
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_entered.json': J3V2__FROZEN_ARTIFACT_SHA256['entered'],
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_result.json': J3V2__FROZEN_ARTIFACT_SHA256['result'],
            FAILURE_EVIDENCE: RECORD_HASHES[FAILURE_EVIDENCE],
        }
        p['report.input_summary'] = (summary.get('component') == 'endpoint_transport'
            and summary.get('inherited_source_count') == 1850
            and summary.get('input_digest') == J3V2__FROZEN_INPUT_DIGEST
            and summary.get('record_hashes') == J3V2__FROZEN_RECORD_HASHES
            and summary.get('new_input_hashes') == expected_new)
    p['report.implementation_hashes'] = report.get('implementation_hashes') == {
        'scripts/run_original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['runner'],
        'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['kernel']}
    config = report.get('kernel_literal_configuration')
    p['report.literal_configuration'] = (isinstance(config, dict)
        and config.get('SCHEMA') == 'original_native_J3_endpoint_transport_v1'
        and config.get('RECORD_NAMES') == list(J3V2__FROZEN_RECORD_HASHES)
        and config.get('RECORD_HASHES') == J3V2__FROZEN_RECORD_HASHES
        and config.get('STAGE_FLAGS') == J3V2__FROZEN_STAGE_FLAGS
        and config.get('FALSE_FLAGS') == J3V2__FROZEN_FALSE_FLAGS)
    if isinstance(evidence_source_text, str):
        evidence_bytes = evidence_source_text.encode('utf-8')
        p['report.evidence_bytes'] = len(evidence_bytes) == 566962
    else:
        p['report.evidence_bytes'] = False
    p['evidence.schema'] = (isinstance(evidence, _Mapping)
                            and evidence.get('schema') ==
                            'original_native_J3_endpoint_transport_failure_evidence_v1')
    p['evidence.scope'] = (evidence.get('component') == 'endpoint_transport'
        and evidence.get('metadata_only') is True
        and evidence.get('failed_attempt_preserved') is True
        and evidence.get('raw_bytes_preserved') is True
        and evidence.get('runtime_objects_recovered') is False
        and evidence.get('numeric_values_recovered') is False
        and evidence.get('accepted_identity') is False
        and evidence.get('source_prefix_objects_recovered') is False
        and evidence.get('do_not_reenter') is True
        and evidence.get('error_stage') == 'new_R_U_C1'
        and evidence.get('attempt_count') == 1)
    p['evidence.full_result'] = (isinstance(evidence.get('full_result'), dict)
                                 and _J3V2__failure_ok(evidence['full_result']))
    p['evidence.freeze'] = evidence.get('freeze_metadata') == {
        'accepted_ancestor_operator_calls': 0,
        'accepted_source_consumer_calls': 0,
        'do_not_reenter': True,
        'full_result_saved': True,
        'numeric_evaluations': 0,
        'old_operators_replayed': False,
        'old_proof_functions_called': 0,
        'parent_commit': 'a45b034ddd9f1982e69d6ac6836838667493b3da',
        'parent_report_sha256': RECORD_HASHES[PARENT],
        'runner_caught_exception': True,
        'source_consumer_attempts': 1,
        'source_execution_occurred': True,
    }
    expected_inventory = {
        'scripts/run_original_native_J3_endpoint_transport.py':
            (J3V2__FROZEN_ARTIFACT_SHA256['runner'], 13060),
        'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py':
            (J3V2__FROZEN_ARTIFACT_SHA256['kernel'], 130268),
        'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_source.json':
            (J3V2__FROZEN_ARTIFACT_SHA256['source'], 275391),
        'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_entered.json':
            (J3V2__FROZEN_ARTIFACT_SHA256['entered'], 662),
        'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_result.json':
            (J3V2__FROZEN_ARTIFACT_SHA256['result'], 1762),
    }
    items = evidence.get('artifact_inventory') if isinstance(evidence, _Mapping) else None
    by_path = {item.get('path'): item for item in items
               if isinstance(item, dict)} if isinstance(items, list) else {}
    p['evidence.inventory_shape'] = (isinstance(items, list)
        and len(items) == len(expected_inventory)
        and set(by_path) == set(expected_inventory)
        and all(item.get('sha256') == expected[0]
                and item.get('bytes') == expected[1]
                and item.get('encoding') == 'base64'
                and isinstance(item.get('base64_bytes'), str)
                and item.get('base64_bytes')
                for path, expected in expected_inventory.items()
                for item in [by_path.get(path)]))
    return p


def _J3V2__emit_diagnostic(callback, events, stage, payload):
    event = {'stage': stage, 'payload': payload}
    events.append(event)
    if callback is not None:
        require(callable(callback), 'diagnostic_callback must be callable')
        callback(stage, payload)


def derive(records, expected_input_digest, diagnostic_callback=None):
    """Authenticate the twelve records and derive the source-only J3 leaf."""
    require(isinstance(records, dict) and set(records) == set(RECORD_NAMES),
            'Changed record inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    decoded, checks, events = {}, {}, []

    def check(name, condition):
        require(bool(condition), 'Source prerequisite: ' + name)
        checks[name] = True
        return True

    for name in RECORD_NAMES:
        item = records[name]
        path = name if '/' in name else 'experiments/root_st073/' + name
        check('authenticated:' + name,
              isinstance(item, dict)
              and item.get('raw_sha256') == RECORD_HASHES[name]
              and item.get('raw_authenticated_before_decode') is True
              and item.get('source_path') == path
              and item.get('is_gzip') is name.endswith('.gz')
              and isinstance(item.get('source_text'), str)
              and _hashlib.sha256(item['source_text'].encode('utf-8')).hexdigest() ==
                  item.get('decoded_sha256'))
        decoded[name] = (item['source_text'] if name.endswith('.py')
                         else _json.loads(item['source_text']))

    parent = decoded[PARENT]
    loop, mod = decoded[LOOP], decoded[MOD]
    background, weighted = decoded[BACKGROUND], decoded[WEIGHTED]
    report, evidence = decoded[FAILED_REPORT], decoded[FAILURE_EVIDENCE]
    old = parent.get('component_result', {}) if isinstance(parent, dict) else {}
    check('accepted_parent_scope',
          parent.get('source_only') is True
          and old.get('passed') is True
          and parent.get('full_goal_complete') is False
          and len(parent.get('source_pins', {})) == 1850
          and old.get('J2_B_identity_proved') is True
          and old.get('ordinary_Z_orders') == [0, 1]
          and old.get('B_identity_scope', {}).get('cumulative_B_identity_nodes') == [520, 521, 522]
          and old.get('B_identity_scope', {}).get('unresolved_B_identity_nodes') == list(range(523, 536))
          and old.get('native_IN_reduced_source_identity') ==
              'I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)')
    check('inherited_input_pins', all(
        parent['source_pins'].get(name) == RECORD_HASHES[name]
        for name in RECORD_NAMES[:10] if name != PARENT))
    check('inherited_sigma_gamma_pins',
          parent['source_pins'].get(J3V2__SIGMA_SOURCE) == J3V2__SIGMA_SOURCE_SHA256
          and parent['source_pins'].get(J3V2__FLAT_SIGMA_SOURCE) ==
              J3V2__FLAT_SIGMA_SOURCE_SHA256)
    check('same_source_family', all(
        isinstance(value, dict) and value == FAMILY
        for value in (old.get('source_family'), loop.get('source_family'),
                      mod.get('source_family'), weighted.get('source_family'),
                      background.get('source_family'))))
    check('archive_creator_input_pins',
          loop.get('input_hashes', {}).get(LOOP_CREATOR) == RECORD_HASHES[LOOP_CREATOR]
          and mod.get('input_hashes', {}).get(MOD_CREATOR) == RECORD_HASHES[MOD_CREATOR]
          and background.get('input_hashes', {}).get(C6_CREATOR) == RECORD_HASHES[C6_CREATOR]
          and background.get('input_hashes', {}).get(BRIDGE_CREATOR) == RECORD_HASHES[BRIDGE_CREATOR])
    check('retained_ledger', all(
        parent.get('retained_construction_data', {}).get(key) == value
        for key, value in {'cumulative_terms': 212, 'raw_signed_terms': 213,
                           'reduced_terms': 183, 'six_signed_channels_retained': True,
                           'all_old_error_data_retained': True,
                           'reserved_contribution_registered_once': True,
                           'inherited_errors_deleted': False}.items()))
    wrapper_predicates = _J3V2__failed_wrapper_guards(
        report, evidence, records[FAILURE_EVIDENCE].get('source_text'))
    check('authenticated_failed_report_and_evidence',
          bool(wrapper_predicates) and all(wrapper_predicates.values()))
    parent_pins = {
        'artifact_raw_sha256': dict(J3V2__FROZEN_ARTIFACT_SHA256),
        'implementation_hashes': {
            'scripts/run_original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['runner'],
            'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py': J3V2__FROZEN_ARTIFACT_SHA256['kernel']},
        'input_digest': J3V2__FROZEN_INPUT_DIGEST,
        'record_hashes': dict(J3V2__FROZEN_RECORD_HASHES),
        'decoded_record_hashes': dict(J3V2__FROZEN_DECODED_RECORD_HASHES),
        'ast_body_sha256': dict(J3V2__FROZEN_AST_BODY_SHA256),
    }
    raw_artifacts = {
        'kernel': _J3V2__inventory_raw(evidence,
            'src/openai_ns_reconstruction/original_native_J3_endpoint_transport.py',
            J3V2__FROZEN_ARTIFACT_SHA256['kernel']),
        'runner': _J3V2__inventory_raw(evidence,
            'scripts/run_original_native_J3_endpoint_transport.py',
            J3V2__FROZEN_ARTIFACT_SHA256['runner']),
        'source': _J3V2__inventory_raw(evidence,
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_source.json',
            J3V2__FROZEN_ARTIFACT_SHA256['source']),
        'entered': _J3V2__inventory_raw(evidence,
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_entered.json',
            J3V2__FROZEN_ARTIFACT_SHA256['entered']),
        'result': _J3V2__inventory_raw(evidence,
            'tmp/original_native_J3_endpoint_transport/new_endpoint_transport_result.json',
            J3V2__FROZEN_ARTIFACT_SHA256['result']),
    }
    prefix = J3V2__project_boolean_prefix(raw_artifacts['kernel'], raw_artifacts['runner'],
        raw_artifacts['source'], raw_artifacts['entered'], raw_artifacts['result'],
        parent_pins)
    _J3V2__emit_diagnostic(diagnostic_callback, events, 'boolean_prefix', prefix)
    check('failed_prefix_boolean_certificate', _J3V2__certificate_ok(prefix))
    cutoff_laws = _J3V2__sigma_source_laws(decoded[GAMMA_CREATOR], decoded[C6_CREATOR], background, parent['source_pins'])
    _J3V2__emit_diagnostic(diagnostic_callback, events, 'cutoff_source_laws', cutoff_laws)
    check('exact_cutoff_endpoint_laws', _J3V2__sigma_law_contract(cutoff_laws))
    ru = J3V2__prove_R_U_transport(
        loop, mod, background, decoded[LOOP_CREATOR], decoded[MOD_CREATOR],
        decoded[C6_CREATOR], decoded[BRIDGE_CREATOR], prefix,
        {'loop': records[LOOP]['raw_sha256'], 'mod': records[MOD]['raw_sha256'],
         'background': records[BACKGROUND]['raw_sha256']},
        {'gamma_source_sha256': parent['source_pins'][J3V2__SIGMA_SOURCE],
         'flat_sigma_source_sha256': parent['source_pins'][J3V2__FLAT_SIGMA_SOURCE],
         'cutoff_source_laws': cutoff_laws})
    _J3V2__emit_diagnostic(diagnostic_callback, events, 'R_U', ru)
    check('new_R_U_C1', ru.get('J3_R_C1_proved') is True
          and ru.get('J3_U_C1_proved') is True
          and _J3V2__ru_result_contract_ok(ru))
    mq = J3V2__prove_physical_Q_transport(
        loop, mod, background, decoded[LOOP_CREATOR], decoded[MOD_CREATOR],
        decoded[C6_CREATOR], decoded[BRIDGE_CREATOR], prefix, ru,
        {'loop': records[LOOP]['raw_sha256'], 'mod': records[MOD]['raw_sha256'],
         'background': records[BACKGROUND]['raw_sha256']},
        {'gamma_source_sha256': parent['source_pins'][J3V2__SIGMA_SOURCE],
         'flat_sigma_source_sha256': parent['source_pins'][J3V2__FLAT_SIGMA_SOURCE],
         'cutoff_source_laws': cutoff_laws})
    _J3V2__emit_diagnostic(diagnostic_callback, events, 'M_Q', mq)
    check('new_physical_m_Z2_and_Q_C1', mq.get('physical_m_Z2_proved') is True
          and mq.get('finite_N_m_Z2_proved') is True
          and mq.get('J3_Q_C1_proved') is True)

    rows = weighted.get('definition_graph', {}).get('nodes')
    weighted_result = weighted.get('component_result', {})
    check('weighted_source_graph_scope', isinstance(rows, list) and len(rows) == 545
          and rows == weighted_result.get('definition_graph', {}).get('nodes')
          and weighted_result.get('actual_native_active_S1_weighted_IBP_C1_identity_proved') is True)

    def row_is(number, **expected):
        return all(_J3V2__row(rows, number).get(key) == value
                   for key, value in expected.items()) if _J3V2__row(rows, number) else False

    check('actual_J3_weighted_jump', row_is(523, operation='actual_internal_primitive_jump',
          arguments=[377, 367], jump_index=3, formula='J_i=B_next_lower-B_prev_upper',
          is_actual_function=True, retained=True, zero_assumption=False,
          source_identity_proved=False))
    check('actual_J3_B_endpoint_recipe',
          row_is(367, operation='actual_B_endpoint_function', chart='first_switch',
                 endpoint='upper', coordinate=75, domain=76, radius=77,
                 endpoint_bound='domain[1]',
                 formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
          and row_is(377, operation='actual_B_endpoint_function', chart='second_switch',
                     endpoint='lower', coordinate=91, domain=92, radius=93,
                     endpoint_bound='domain[0]',
                     formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
          and _J3V2__row(rows, 367).get('physical_ports', {}).get('U') == 80
          and _J3V2__row(rows, 367).get('physical_ports', {}).get('Q') == 84
          and _J3V2__row(rows, 377).get('physical_ports', {}).get('U') == 96
          and _J3V2__row(rows, 377).get('physical_ports', {}).get('Q') == 100)
    check('actual_J3_source_ports', all(row_is(i, source_nodes=values) for i, values in (
        (75, [73190]), (76, [0, 1]), (77, [85905]),
        (80, [148085, 148086, 148087, 148088, 148089, 148090]),
        (84, [148180, 148184, 148190, 148196, 148202]),
        (91, [75949]), (92, [0, 1]), (93, [88600]),
        (96, [154247, 154248, 154249, 154250, 154251, 154252]),
        (100, [154342, 154346, 154352, 154358, 154364]))))
    check('actual_J3_port_archive_types', all(
        row_is(i, operation=operation, source_record=source, chart=chart,
               quantity=quantity)
        for i, operation, source, chart, quantity in (
            (75, 'actual_native_geometry_coordinate', LOOP, 'first_switch', 'coordinate'),
            (77, 'actual_native_geometry_radius', LOOP, 'first_switch', 'radius'),
            (91, 'actual_native_geometry_coordinate', LOOP, 'second_switch', 'coordinate'),
            (93, 'actual_native_geometry_radius', LOOP, 'second_switch', 'radius'),
            (80, 'actual_native_physical_U', MOD, 'first_switch', 'U'),
            (84, 'actual_native_physical_Q', MOD, 'first_switch', 'Q'),
            (96, 'actual_native_physical_U', MOD, 'second_switch', 'U'),
            (100, 'actual_native_physical_Q', MOD, 'second_switch', 'Q'))))
    actual_l = weighted_result.get('source_interface', {}).get('actual_L', {})
    check('actual_J3_physical_measure_and_common_L',
          all(row_is(i, source_nodes=[9], source_record=LOOP) for i in (78, 79, 94, 95))
          and _J3V2__row(rows, 367).get('measure') ==
              _J3V2__row(rows, 377).get('measure') == 'R^2*Jacobian*dnative'
          and all(actual_l.get(key) == value for key, value in {
              'source_record': MOD, 'positive_quotient_denominator': 25272,
              'formula': '1-delta*Z^2', 'Z_independent_delta_port': True}.items()))
    check('accepted_uniform_B_creator_recipe',
          old.get('component_proofs', {}).get('J2_B', {}).get('source_ast_hashes', {}).get(
              'active_S1_creator_sha256') == RECORD_HASHES[WEIGHTED_CREATOR]
          and old.get('native_weighted_B_definition', {}).get('formula') ==
              'R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)')
    graph = [
        {'id': 'failed_prefix_projection', 'rule': 'authenticate frozen failed attempt and project Boolean prefix'},
        {'id': 'R_U_transport', 'rule': 'fresh source root and loop composition at common R1 endpoint',
         'depends_on': ['failed_prefix_projection']},
        {'id': 'physical_m_Q_transport', 'rule': 'finite-N memory then physical m and positive-L quotient',
         'depends_on': ['failed_prefix_projection', 'R_U_transport']},
        {'id': 'B_C1_composition', 'rule': 'identical R U Q C1 ports in the uniform actual B recipe',
         'depends_on': ['R_U_transport', 'physical_m_Q_transport']},
        {'id': 'signed_J3', 'rule': 'actual right-lower minus left-upper [377,367]',
         'depends_on': ['B_C1_composition']},
        {'id': 'native_IN_reduction', 'rule': 'substitute proved zero J3 into accepted native source identity',
         'depends_on': ['signed_J3']}]
    kept = parent['retained_construction_data']
    stage_flags = dict.fromkeys(STAGE_FLAGS, True)
    false_flags = dict.fromkeys(FALSE_FLAGS, False)
    result = dict(schema=SCHEMA, component='endpoint_transport', passed=True,
        source_only=True, **stage_flags, **false_flags,
        stage_flags=dict(stage_flags), nonclaims=dict(false_flags),
        J3_B_identity_proved=True, J3_R_identity_proved=True,
        J3_U_identity_proved=True, J3_Q_identity_proved=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        owned_Z_cell=['0', '1/2'], input_digest=expected_input_digest,
        native_six_signed_channels=['H', 'M', 'K', 'A', 'B', 'C'],
        generic_five_channels=['m', 'h', 'k', 'e', 'p'],
        source_family=FAMILY, source_function_identity_only=True,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_checks=checks, proof_graph=graph, proof_graph_sha256=digest(graph),
        component_proofs={'failed_prefix': prefix, 'cutoff_source_laws': cutoff_laws, 'R_U': ru, 'physical_m_Q': mq},
        diagnostics=list(events), diagnostic_events=list(events),
        diagnostic_payloads=list(events),
        retained_construction_data=kept,
        retained_integral_term_inventory=parent.get('retained_integral_term_inventory'),
        old_errors_frozen=parent.get('old_errors_frozen'),
        quantity_distinctions=old.get('quantity_distinctions'),
        inherited_R_identity_jumps=old.get('inherited_R_identity_jumps'),
        B_identity_scope=dict(weighted_jump_node=523, jump_index=3,
            left_B_node=367, right_B_node=377, jump_arguments=[377, 367],
            ordinary_Z_orders=[0, 1], endpoint_substitutions={'73190': 1, '75949': 0},
            common_L_node=25272, newly_proved_B_identity_nodes=[523],
            inherited_B_identity_nodes=[520, 521, 522],
            cumulative_B_identity_nodes=[520, 521, 522, 523],
            unresolved_B_identity_nodes=list(range(524, 536)),
            remaining_B_identity_count=12, all_native_B_jumps_zero_proved=False),
        unresolved_B_identity_nodes=list(range(524, 536)), remaining_B_identity_count=12,
        native_weighted_B_definition=old.get('native_weighted_B_definition'),
        B_C1_derivative_rule=old.get('component_proofs', {}).get('J2_B', {}).get('C1_derivative_rule'),
        native_IN_inherited_identity=old.get('native_IN_inherited_identity',
            'I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)'),
        native_IN_reduced_source_identity='I_N=B_lastupper-B_firstlower-sum(J_4..J_15)+sum17(integral K_i dy)',
        native_IN_reduction_is_source_only=True, archived_all_16_J_terms_retained=True,
        source_jump_zero_does_not_delete_quantitative_errors=True,
        unresolved_gates=[], remaining_construction_gates=[
            'the other twelve native B jumps and same-source quantitative packet transport',
            'independent terminal functional moments and retained quantitative errors',
            'numerical containment and domain admission',
            'heat stress flat remainder time recursion oscillatory correction independent Cartesian validation'],
        source_inspection_audit=dict(new_source_consumer_calls=1,
            accepted_source_consumer_calls=0,
            accepted_ancestor_operator_calls=0, old_proof_functions_called=0,
            generic_DAG_walks=0, numerical_evaluations=0, quadratures=0,
            recurrence_calls=0, runtime_model_metadata=None), runtime_model_metadata=None)
    return result
