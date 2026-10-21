"""Saved-prefix authentication and new loop-only J2 source composition.

No V1 comparison, memory or radial proof is reexecuted. The failed V1 stays
failed; only the fail-fast execution path's boolean projection is recovered.
Source identities do not admit a physical layer or numerical containment.
"""
import ast
import base64
import hashlib
import json

SCHEMA = 'original_native_J2_mean_endpoint_composition_v2'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_mean_loop_input_join.json'
TRANSPORT = 'lei_ren_part1_paper_compliant_current_original_native_J2_U_endpoint_transport_v2.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
MOD = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
EVIDENCE = 'experiments/root_st073/original_native_J2_mean_endpoint_v1_failure_evidence.json'
RECORD_NAMES = [PARENT, TRANSPORT, LOOP, MOD, CREATOR, EVIDENCE]
RECORD_HASHES = {
    PARENT: '7143e443ec041b4d94d4a80a42dc9107dcee6d36284708661e7f25905759b1e9',
    TRANSPORT: 'be0cad628f44cec0713144f074cbbe4c7b0efc23feab9d3d482a116235d3737b',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    CREATOR: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    EVIDENCE: 'fe8d0169c01d7023ada2e4644f59f1ccfbf7978e67542f4e639e846d725c9cd1',
}
STAGE_FLAGS = ['actual_J2_mean_ports_bound', 'comparison_endpoint_C1_identity_proved',
               'original_memory_Z2_endpoint_identity_proved',
               'J2_mean_Ey_C1_endpoint_identity_proved',
               'J2_mean_Vy_C1_endpoint_identity_proved',
               'J2_mean_p2_C1_endpoint_identity_proved',
               'J2_Bloop_C1_endpoint_identity_proved',
               'actual_J2_U_C1_endpoint_identity_proved',
               'source_prefix_certificate_recovered', 'no_prior_proof_parts_reexecuted',
               'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['J2_Q_identity_proved', 'J2_B_identity_proved',
               'all_native_B_jumps_zero_proved', 'numeric_values_computed',
               'source_C5_installation_claimed', 'numeric_containment_claimed',
               'old_error_terms_removed', 'global_five_moment_closure',
               'time_order_recursion_claimed', 'full_goal_complete',
               'uninterpreted_graph_parameter_identity_proved', 'old_operators_replayed',
               'V1_stage_accepted', 'V1_runtime_intermediate_objects_recovered']
FAMILY = {
    'actual_five_defect_family_sha256': '3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894',
    'implicit_source_sha256': '5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae',
    'datum_enclosure_sha256': 'dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010',
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def derive(records, expected_input_digest):
    require(set(records) == set(RECORD_NAMES), 'Changed record inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    decoded, checks = {}, {}

    def check(name, condition):
        require(bool(condition), 'Source prerequisite: ' + name)
        checks[name] = True
        return True

    for name in RECORD_NAMES:
        item = records[name]
        path = name if '/' in name else 'experiments/root_st073/' + name
        check('authenticated:' + name, item['raw_sha256'] == RECORD_HASHES[name]
              and item['raw_authenticated_before_decode'] is True and item['source_path'] == path
              and item['is_gzip'] is name.endswith('.gz')
              and hashlib.sha256(item['source_text'].encode()).hexdigest() == item['decoded_sha256'])
        decoded[name] = item['source_text'] if name.endswith('.py') else json.loads(item['source_text'])
    parent, transport, loop, mod, evidence = [decoded[name] for name in
        (PARENT, TRANSPORT, LOOP, MOD, EVIDENCE)]
    check('accepted_parent_source_scope', parent['source_only'] is True
          and parent['component_result']['passed'] is True and parent['full_goal_complete'] is False
          and len(parent['source_pins']) == 1840
          and parent['R_identity_jumps'] == [2, 3, 4, 6, 8, 9, 10]
          and parent['B_identity_scope']['unresolved_B_identity_nodes'] == list(range(522, 536))
          and parent['Q_identity_scope']['J2_Q_identity_proved'] is False)
    check('inherited_input_pins', all(parent['source_pins'][name] == RECORD_HASHES[name]
          for name in (TRANSPORT, LOOP, MOD, CREATOR)))
    check('same_source_family', loop['source_family'] == mod['source_family'] == FAMILY
          and transport['component_result']['source_family'] == FAMILY)
    prefix = recover_prefix(evidence, parent, transport,
        lambda name, value: check('saved_prefix:' + name, value))
    reduction = transport['component_result']['U_difference_identity']
    check('accepted_U_difference_reduction',
          transport['component_result']['J2_U_difference_C1_reduction_proved'] is True
          and reduction['exact_difference_reduction_proved'] is True
          and reduction['ordinary_Z_orders'] == [0, 1]
          and reduction['endpoint_substitutions'] == {'70265': 1, '73190': 0}
          and reduction['left_coordinate_node'] == 70265
          and reduction['right_coordinate_node'] == 73190
          and reduction['left_source_nodes'] == [140685, 140686]
          and reduction['right_source_nodes'] == [148085, 148086]
          and reduction['common_Pstar_node'] == 33572
          and reduction['common_inverse_N_node'] == 99411
          and reduction['retained_Bloop_over_Pstar_left_nodes'] == [103809, 103810]
          and reduction['retained_Bloop_over_Pstar_right_nodes'] == [106437, 106438]
          and reduction['difference_formula'] == 'U_left,j(1)-U_right,j(0)=(Pstar/N)*(Bloop_over_Pstar_left,j(1)-Bloop_over_Pstar_right,j(0))')
    check('source_factory_unadmitted_numeric_flags_retained',
          loop['current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed'] is False
          and loop['current_original_prescribed_loop_C5_numeric_containment_installed'] is False)
    creator_record = dict(source_path=records[CREATOR]['source_path'], source_text=decoded[CREATOR])
    check('new_scoped_source_body_fingerprints', loop_v2_scoped_fingerprints_ok(decoded[CREATOR]))
    # This is a new certificate built from recovered booleans and authenticated
    # literal source-node bindings, not an unsaved V1 runtime object.
    mean = dict(proved=all(prefix['fields_C1'].values()), ordinary_Z_orders=[0, 1],
        input_scope='original_mean_and_phase_only', modified_total_radial_fields_used=False,
        hbB_hbS_identity_proved=prefix['width_contract']['interpreted_parameter_identity_proved'],
        V_identity_proved=prefix['fields_C1']['V_C1'],
        phase=dict(proved=prefix['fields_C1']['phase_C1'], ordinary_Z_orders=[0, 1]),
        fields={name: dict(proved=prefix['fields_C1'][name + '_C1'], pairs=pairs)
                for name, pairs in MEMPRE_MEAN_PAIRS.items()},
        certificate_origin='new boolean projection certificate; no V1 runtime objects recovered')
    composed = prove_loop_composition_v2(loop, mod, creator_record, mean,
        lambda name, value: check('new_loop:' + name, value))
    check('new_Bloop_and_U_source_identities', composed['Bloop_C1_proved'] is True
          and composed['full_U_C1_proved'] is True)
    retained = parent['retained_construction_data']
    check('retained_ledger', all(retained[key] == value for key, value in
          dict(cumulative_terms=212, raw_signed_terms=213, reduced_terms=183,
               six_signed_channels_retained=True, all_old_error_data_retained=True,
               reserved_contribution_registered_once=True, inherited_errors_deleted=False).items()))
    proof = [dict(id='saved_V1_prefix', rule='exact saved traceback and fail-fast AST dominance',
                  recovered_runtime_objects=False, recovered_numeric_values=False),
             dict(id='new_loop_source_guard', rule='authenticated exact source including the original upper assignment',
                  depends_on=['saved_V1_prefix']),
             dict(id='new_Bloop', rule='same original mean/phase C1 factory; common eta/dstar checked here',
                  depends_on=['saved_V1_prefix', 'new_loop_source_guard']),
             dict(id='new_U', rule='new Bloop identity and authenticated prior U difference reduction',
                  depends_on=['new_Bloop'])]
    return dict(schema=SCHEMA, component='mean_endpoint_composition_v2', passed=True, source_only=True,
        **dict.fromkeys(STAGE_FLAGS, True), **dict.fromkeys(FALSE_FLAGS, False),
        stage_flags=dict.fromkeys(STAGE_FLAGS, True), nonclaims=dict.fromkeys(FALSE_FLAGS, False),
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, owned_Z_cell=['0', '1/2'],
        input_digest=expected_input_digest, source_family=FAMILY,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_function_identity_only=True, source_checks=checks,
        proof_graph=proof, proof_graph_sha256=digest(proof),
        component_proofs=dict(saved_prefix=prefix, new_loop=composed),
        source_binding_scope=parent['source_binding_scope'],
        mean_C1_certificate=mean, retained_construction_data=retained,
        quantity_distinctions=parent['quantity_distinctions'],
        inherited_R_identity_jumps=[2, 3, 4, 6, 8, 9, 10],
        unresolved_B_identity_nodes=list(range(522, 536)),
        U_identity_scope=dict(weighted_jump_node=522, jump_index=2, full_U_C1_identity_proved=True,
            ordinary_Z_orders=[0, 1], left_source_nodes=[140685, 140686],
            right_source_nodes=[148085, 148086], endpoint_substitutions={'70265': 1, '73190': 0},
            source_interpretation_required=True),
        unresolved_gates=[], remaining_construction_gates=[
            'physical Q and all fourteen native weighted B jumps',
            'independent terminal functional moments and retained quantitative errors',
            'heat exterior, stress/remainder, temporal n recursion, realizable corrections and Cartesian NS validation'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            V1_comparison_memory_radial_functions_called=0, generic_DAG_walks=0,
            numerical_evaluations=0, quadratures=0, recurrence_calls=0, runtime_model_metadata=None),
        runtime_model_metadata=None)


"""Unexecuted V1 prefix extraction rules for a loop-only V2 leaf.

This fragment is deliberately data-only.  It authenticates the five embedded
V1 bytes, parses the frozen source and runner with the standard-library AST,
and validates fail-fast control-flow dominance.  It never imports a project
module, calls a source proof function, loads C6/R100 input records, or returns
an intermediate runtime object or a numerical value.

The caller must provide ``parent`` and ``transport`` as the already accepted,
immutable decoded parent reports.  Their accepted source hashes are checked
against the embedded V1 input metadata; no raw parent/transport bytes are
reconstructed here.
"""

import ast as MEMPRE_ast
import base64 as MEMPRE_base64
import hashlib as MEMPRE_hashlib
import json as MEMPRE_json


MEMPRE_SCHEMA = "native_J2_saved_prefix_validator_v2"
MEMPRE_EVIDENCE_SCHEMA = "original_native_J2_mean_endpoint_v1_failure_evidence_v1"
MEMPRE_SOURCE_SCHEMA = "native_J2_mean_endpoint_composition_inputs_v1"
MEMPRE_EXPECTED_KERNEL_SCHEMA = "original_native_J2_mean_endpoint_composition_v1"
MEMPRE_PARENT = "lei_ren_part1_paper_compliant_current_original_native_J2_mean_loop_input_join.json"
MEMPRE_TRANSPORT = "lei_ren_part1_paper_compliant_current_original_native_J2_U_endpoint_transport_v2.json"
MEMPRE_PARENT_SHA = "7143e443ec041b4d94d4a80a42dc9107dcee6d36284708661e7f25905759b1e9"
MEMPRE_TRANSPORT_SHA = "be0cad628f44cec0713144f074cbbe4c7b0efc23feab9d3d482a116235d3737b"
MEMPRE_KERNEL_SHA = "b10b2600f1f1a3b4d0f387e66e412da668b666dcd8996a71305caaa071e3f65f"
MEMPRE_RUNNER_SHA = "7581781d037b9d6b59a0cd86934d3521712e3588c17685d93bfbc2affb337400"
MEMPRE_SOURCE_SHA = "85421129e9d7ac7a3f9a6288f8bb48c184e738cc9101d7c4dad3d3aa274d45c0"
MEMPRE_ENTERED_SHA = "439bbacb222c1da88b84858147b6a6df4012d3a3c0ee4b9473dc1ced85514981"
MEMPRE_RESULT_SHA = "483ec74f99290a5f77fd79de0ce896f8a0fecea70dcc8cf8f493dbdc8f9fa1f2"
MEMPRE_TRACEBACK_SHA = "71595f51f948aada2ad1438fdc787315e93c73798b2229db0123530e8289b573"
MEMPRE_LOOP_CREATOR_SHA = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
MEMPRE_INPUT_DIGEST = "8873772c8b689bd6f49e7ccf50ae8bd02761bb48646e1e72ec0ac9961794f025"
MEMPRE_IMPLEMENTATIONS = {
    "scripts/run_original_native_J2_mean_endpoint_composition.py": MEMPRE_RUNNER_SHA,
    "src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py": MEMPRE_KERNEL_SHA,
}
MEMPRE_ARTIFACTS = {
    "scripts/run_original_native_J2_mean_endpoint_composition.py":
        {"sha256": MEMPRE_RUNNER_SHA, "bytes": 11878},
    "src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py":
        {"sha256": MEMPRE_KERNEL_SHA, "bytes": 112579},
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_source.json":
        {"sha256": MEMPRE_SOURCE_SHA, "bytes": 273002},
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_entered.json":
        {"sha256": MEMPRE_ENTERED_SHA, "bytes": 683},
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_result.json":
        {"sha256": MEMPRE_RESULT_SHA, "bytes": 3206},
}
MEMPRE_AUDITS = {
    "native_J2_mean_endpoint_composition_v1_failure_prefix_audit_2026_10_21.json":
        {"sha256": "c9e3eeb04c05ff7ec698aa3e52e1684b8865cae1c3346d7fe2478b4a525ba725", "bytes": 11234},
    "native_J2_saved_only_mean_prefix_adversarial_2026_10_21.json":
        {"sha256": "263bd0b00e560452c9a7b10598f60a2b7ab9d48cf5dfad8ca756ab0c7cdec608", "bytes": 9332},
}
MEMPRE_MEAN_PAIRS = {
    "E": [[72108, 87550], [72111, 87553]],
    "E_y": [[72348, 87790], [72351, 87793]],
    "V_y": [[72382, 87824], [72383, 87825]],
    "p2": [[72981, 88423], [72982, 88424]],
}
MEMPRE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


def MEMPRE_require(condition, message):
    if not condition:
        raise AssertionError(message)


def MEMPRE_sha256(raw):
    return MEMPRE_hashlib.sha256(raw).hexdigest()


def MEMPRE_accept(check, name, condition):
    if type(condition) is not bool:
        raise TypeError(name + " must be bool")
    if check(name, condition) is not True:
        raise ValueError("prefix premise rejected: " + name)


def MEMPRE_artifact_bytes(item, expected, name):
    MEMPRE_require(isinstance(item, dict), name + " metadata must be a mapping")
    MEMPRE_require(item.get("encoding") == "base64", name + " encoding changed")
    MEMPRE_require(item.get("sha256") == expected["sha256"], name + " hash metadata changed")
    MEMPRE_require(item.get("bytes") == expected["bytes"], name + " byte count changed")
    encoded = item.get("base64_bytes")
    MEMPRE_require(isinstance(encoded, str), name + " base64 payload missing")
    raw = MEMPRE_base64.b64decode(encoded, validate=True)
    MEMPRE_require(len(raw) == expected["bytes"], name + " decoded byte count changed")
    MEMPRE_require(MEMPRE_sha256(raw) == expected["sha256"], name + " embedded bytes changed")
    return raw


def MEMPRE_top_function(tree, name):
    found = [node for node in tree.body
             if isinstance(node, (MEMPRE_ast.FunctionDef, MEMPRE_ast.AsyncFunctionDef))
             and node.name == name]
    MEMPRE_require(len(found) == 1, "expected one top-level function: " + name)
    return found[0]


def MEMPRE_const(node, value):
    return isinstance(node, MEMPRE_ast.Constant) and node.value == value


def MEMPRE_list(node, values):
    return (isinstance(node, MEMPRE_ast.List)
            and len(node.elts) == len(values)
            and all(MEMPRE_const(item, value) for item, value in zip(node.elts, values)))


def MEMPRE_subscript(node, base, key):
    return (isinstance(node, MEMPRE_ast.Subscript)
            and isinstance(node.value, MEMPRE_ast.Name) and node.value.id == base
            and MEMPRE_const(node.slice, key))


def MEMPRE_keyword_map(call):
    MEMPRE_require(isinstance(call, MEMPRE_ast.Call)
                   and isinstance(call.func, MEMPRE_ast.Name)
                   and call.func.id == "dict", "mean certificate must use dict(...)" )
    result = {}
    for keyword in call.keywords:
        MEMPRE_require(keyword.arg is not None, "mean certificate has **kwargs")
        result[keyword.arg] = keyword.value
    return result


def MEMPRE_target_name(statement):
    if not isinstance(statement, MEMPRE_ast.Assign) or len(statement.targets) != 1:
        return None
    target = statement.targets[0]
    return target.id if isinstance(target, MEMPRE_ast.Name) else None


def MEMPRE_assignment_call(statement, target, callee):
    return (MEMPRE_target_name(statement) == target
            and isinstance(statement.value, MEMPRE_ast.Call)
            and isinstance(statement.value.func, MEMPRE_ast.Name)
            and statement.value.func.id == callee)


def MEMPRE_first_assignment(body, target, callee=None):
    for index, statement in enumerate(body):
        if MEMPRE_target_name(statement) != target:
            continue
        if callee is None or (isinstance(statement.value, MEMPRE_ast.Call)
                              and isinstance(statement.value.func, MEMPRE_ast.Name)
                              and statement.value.func.id == callee):
            return index, statement
    raise AssertionError("missing dominating assignment: " + target)


def MEMPRE_check_function_is_fail_fast(derive):
    nested = [node for node in derive.body
              if isinstance(node, (MEMPRE_ast.FunctionDef, MEMPRE_ast.AsyncFunctionDef))
              and node.name == "check"]
    MEMPRE_require(len(nested) == 1, "derive check callback is not unique")
    check_fn = nested[0]
    MEMPRE_require(not any(isinstance(node, MEMPRE_ast.Try)
                           for node in MEMPRE_ast.walk(check_fn)),
                   "check callback catches exceptions")
    has_require = False
    for node in MEMPRE_ast.walk(check_fn):
        if not (isinstance(node, MEMPRE_ast.Call)
                and isinstance(node.func, MEMPRE_ast.Name)
                and node.func.id == "require"
                and len(node.args) >= 1):
            continue
        condition = node.args[0]
        if (isinstance(condition, MEMPRE_ast.Call)
                and isinstance(condition.func, MEMPRE_ast.Name)
                and condition.func.id == "bool"
                and len(condition.args) == 1
                and isinstance(condition.args[0], MEMPRE_ast.Name)
                and condition.args[0].id == "condition"):
            has_require = True
    MEMPRE_require(has_require, "derive check is not fail-fast")
    MEMPRE_require(any(isinstance(node, MEMPRE_ast.Return)
                       and MEMPRE_const(node.value, True)
                       for node in check_fn.body),
                   "derive check does not return True")


def MEMPRE_require_function_fail_fast(tree):
    require_fn = MEMPRE_top_function(tree, "require")
    has_assertion_raise = False
    for node in MEMPRE_ast.walk(require_fn):
        if not isinstance(node, MEMPRE_ast.If):
            continue
        test = node.test
        if not (isinstance(test, MEMPRE_ast.UnaryOp)
                and isinstance(test.op, MEMPRE_ast.Not)
                and isinstance(test.operand, MEMPRE_ast.Name)
                and test.operand.id == "condition"):
            continue
        if any(isinstance(child, MEMPRE_ast.Raise)
               and isinstance(child.exc, MEMPRE_ast.Call)
               and isinstance(child.exc.func, MEMPRE_ast.Name)
               and child.exc.func.id == "AssertionError"
               for child in MEMPRE_ast.walk(node)):
            has_assertion_raise = True
    MEMPRE_require(has_assertion_raise, "require is not fail-fast")


def MEMPRE_validate_mean_assignment(derive):
    index, statement = MEMPRE_first_assignment(derive.body, "mean", "dict")
    kwargs = MEMPRE_keyword_map(statement.value)
    required = {
        "proved", "ordinary_Z_orders", "input_scope", "modified_total_radial_fields_used",
        "common_eta_dstar_proved", "hbB_hbS_identity_proved", "V_identity_proved",
        "phase", "fields",
    }
    MEMPRE_require(set(kwargs) == required, "mean certificate fields changed")
    MEMPRE_require(MEMPRE_const(kwargs["proved"], True), "mean proved flag changed")
    MEMPRE_require(MEMPRE_list(kwargs["ordinary_Z_orders"], [0, 1]), "mean ordinary Z scope changed")
    MEMPRE_require(MEMPRE_const(kwargs["input_scope"], "original_mean_and_phase_only"),
                   "mean input scope changed")
    MEMPRE_require(MEMPRE_const(kwargs["modified_total_radial_fields_used"], False),
                   "modified radial guard changed")
    # This literal is structurally recorded but deliberately not recovered as a
    # loop-domain proof; V2 must establish eta and d_star independently.
    MEMPRE_require(MEMPRE_const(kwargs["common_eta_dstar_proved"], True),
                   "mean construction changed its eta/d_star marker")
    MEMPRE_require(MEMPRE_subscript(kwargs["hbB_hbS_identity_proved"], "inherited", "source_width_C1"),
                   "mean width binding changed")
    MEMPRE_require(MEMPRE_subscript(kwargs["V_identity_proved"], "inherited", "V_C1"),
                   "mean V binding changed")
    phase = MEMPRE_keyword_map(kwargs["phase"])
    MEMPRE_require(set(phase) == {"proved", "ordinary_Z_orders"}, "mean phase fields changed")
    MEMPRE_require(MEMPRE_subscript(phase["proved"], "inherited", "phase_C1"),
                   "mean phase binding changed")
    MEMPRE_require(MEMPRE_list(phase["ordinary_Z_orders"], [0, 1]),
                   "mean phase ordinary Z scope changed")
    fields_node = kwargs["fields"]
    MEMPRE_require(isinstance(fields_node, MEMPRE_ast.Dict), "mean fields are not a dict")
    fields = {}
    for key, value in zip(fields_node.keys, fields_node.values):
        MEMPRE_require(isinstance(key, MEMPRE_ast.Constant) and isinstance(key.value, str),
                       "mean field key is not literal")
        fields[key.value] = MEMPRE_keyword_map(value)
    MEMPRE_require(set(fields) == set(MEMPRE_MEAN_PAIRS), "mean field set changed")
    refs = {
        "E": ("inherited", "E_C1"),
        "E_y": ("radial", "Ey_C1_proved"),
        "V_y": ("radial", "Vy_C1_proved"),
        "p2": ("memory", "p2_C1_proved"),
    }
    for field, pairs in MEMPRE_MEAN_PAIRS.items():
        MEMPRE_require(set(fields[field]) == {"proved", "pairs"}, field + " certificate fields changed")
        MEMPRE_require(MEMPRE_subscript(fields[field]["proved"], refs[field][0], refs[field][1]),
                       field + " proved binding changed")
        pair_node = fields[field]["pairs"]
        MEMPRE_require(isinstance(pair_node, MEMPRE_ast.List)
                       and len(pair_node.elts) == len(pairs), field + " pair shape changed")
        for row_node, expected_pair in zip(pair_node.elts, pairs):
            MEMPRE_require(MEMPRE_list(row_node, expected_pair), field + " pair IDs changed")
    return index


def MEMPRE_validate_derive_prefix(tree):
    derive = MEMPRE_top_function(tree, "derive")
    MEMPRE_require_function_fail_fast(tree)
    MEMPRE_check_function_is_fail_fast(derive)
    MEMPRE_require(not any(isinstance(node, (MEMPRE_ast.Try, MEMPRE_ast.ExceptHandler))
                           for node in MEMPRE_ast.walk(derive)),
                   "derive has an exception catch that breaks dominance")
    body = derive.body
    inherited_i, _ = MEMPRE_first_assignment(body, "inherited", "inherited_premises")
    comparison_i, comparison_stmt = MEMPRE_first_assignment(body, "comparison", "prove_comparison")
    memory_i, _ = MEMPRE_first_assignment(body, "memory", "prove_memories")
    radial_i, _ = MEMPRE_first_assignment(body, "radial", "prove_radial")
    mean_i = MEMPRE_validate_mean_assignment(derive)
    composed_i, _ = MEMPRE_first_assignment(body, "composed", "prove_loop_composition")
    MEMPRE_require(inherited_i < comparison_i < memory_i < radial_i < mean_i < composed_i,
                   "derive prefix order changed")
    MEMPRE_require(any(isinstance(statement, MEMPRE_ast.For)
                       and "P2_OUTPUT_ROWS" in MEMPRE_ast.unparse(statement)
                       and "actual_p2_and_FF_output_row" in MEMPRE_ast.unparse(statement)
                       for statement in body[:composed_i]),
                   "p2/FF row guard does not dominate loop call")
    prefix_text = "\n".join(MEMPRE_ast.unparse(statement) for statement in body[:composed_i])
    MEMPRE_require("actual_p2_output_bindings" in prefix_text
                   and "actual_p2_output_scale_and_E" in prefix_text,
                   "p2 output binding guards do not dominate loop call")
    MEMPRE_require("source_factory_unadmitted_numeric_flags_retained" in prefix_text,
                   "source installation guard does not dominate loop call")
    MEMPRE_require("mean, lambda name, value" in MEMPRE_ast.unparse(body[composed_i]),
                   "loop call does not receive the constructed mean certificate")
    MEMPRE_require(not any(isinstance(statement, MEMPRE_ast.Return)
                           for statement in body[:composed_i]),
                   "derive can return before the prefix call")
    return {"mean_assignment_index": mean_i, "loop_call_index": composed_i}


def recover_prefix(evidence, parent, transport, check):
    """Return only the authenticated V1 boolean prefix projection.

    This function is intentionally uncalled in this scratch file.  A V2
    consumer may invoke it only with the embedded evidence and accepted parent
    and transport reports.  It does not return ``mean``, graph rows, source
    objects, numeric values, Bloop, full U, or a post-loop U reduction.
    """
    if not callable(check):
        raise TypeError("check must return True for each accepted premise")
    MEMPRE_require(isinstance(evidence, dict), "evidence must be a mapping")
    MEMPRE_require(evidence.get("schema") == MEMPRE_EVIDENCE_SCHEMA, "evidence schema changed")
    MEMPRE_require(evidence.get("V1_do_not_reenter") is True, "V1 reentry flag changed")
    MEMPRE_require(evidence.get("V1_stage_passed") is False, "V1 failure status changed")
    MEMPRE_require(evidence.get("metadata_assembly_only") is True, "evidence assembly mode changed")
    MEMPRE_require(evidence.get("source_functions_called") == 0, "evidence records source execution")
    artifacts = evidence.get("artifacts")
    MEMPRE_require(isinstance(artifacts, dict) and set(artifacts) == set(MEMPRE_ARTIFACTS),
                   "embedded V1 artifact inventory changed")
    raw = {name: MEMPRE_artifact_bytes(artifacts[name], spec, name)
           for name, spec in MEMPRE_ARTIFACTS.items()}
    source_doc = MEMPRE_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_source.json"].decode("utf-8"))
    entered_doc = MEMPRE_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_entered.json"].decode("utf-8"))
    result_doc = MEMPRE_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_result.json"].decode("utf-8"))
    kernel_text = raw["src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py"].decode("utf-8")
    runner_text = raw["scripts/run_original_native_J2_mean_endpoint_composition.py"].decode("utf-8")
    kernel_tree = MEMPRE_ast.parse(kernel_text)
    MEMPRE_ast.parse(runner_text)
    MEMPRE_accept(check, "embedded_source_schema", source_doc.get("schema") == MEMPRE_SOURCE_SCHEMA)
    MEMPRE_accept(check, "embedded_input_digest", source_doc.get("input_digest") == MEMPRE_INPUT_DIGEST)
    MEMPRE_accept(check, "entered_source_link", entered_doc.get("source_sha256") == MEMPRE_SOURCE_SHA)
    MEMPRE_accept(check, "entered_input_link", entered_doc.get("input_digest") == MEMPRE_INPUT_DIGEST)
    MEMPRE_accept(check, "entered_kernel_link", entered_doc.get("kernel_sha256") == MEMPRE_KERNEL_SHA)
    MEMPRE_accept(check, "embedded_implementation_hashes", source_doc.get("implementation_hashes") == MEMPRE_IMPLEMENTATIONS)
    MEMPRE_accept(check, "entered_implementation_hashes", entered_doc.get("implementation_hashes") == MEMPRE_IMPLEMENTATIONS)
    MEMPRE_accept(check, "result_source_link", result_doc.get("source_checkpoint_sha256") == MEMPRE_SOURCE_SHA)
    MEMPRE_accept(check, "result_entry_link", result_doc.get("entry_checkpoint_sha256") == MEMPRE_ENTERED_SHA)
    MEMPRE_accept(check, "result_failure", result_doc.get("passed") is False
                 and result_doc.get("failed_attempt_preserved") is True
                 and result_doc.get("error_type") == "AssertionError"
                 and result_doc.get("error") == "Source prerequisite: loop:creator.authenticated_source_text_guards"
                 and result_doc.get("runner_caught_exception") is True
                 and result_doc.get("runner_full_goal_complete") is False)
    traceback = result_doc.get("traceback")
    MEMPRE_accept(check, "exact_failure_traceback", isinstance(traceback, str)
                 and MEMPRE_sha256(traceback.encode("utf-8")) == MEMPRE_TRACEBACK_SHA
                 and 'original_native_J2_mean_endpoint_composition.py", line 204, in derive' in traceback
                 and 'original_native_J2_mean_endpoint_composition.py", line 1381, in prove_loop_composition' in traceback
                 and 'AssertionError: Source prerequisite: loop:creator.authenticated_source_text_guards' in traceback)
    MEMPRE_accept(check, "embedded_parent_transport_hashes",
                 source_doc.get("record_hashes", {}).get(MEMPRE_PARENT) == MEMPRE_PARENT_SHA
                 and source_doc.get("record_hashes", {}).get(MEMPRE_TRANSPORT) == MEMPRE_TRANSPORT_SHA)
    MEMPRE_accept(check, "embedded_audit_pins",
                 evidence.get("independent_saved_only_audits", {}).get(
                     "native_J2_mean_endpoint_composition_v1_failure_prefix_audit_2026_10_21.json", {}).get("sha256")
                 == MEMPRE_AUDITS["native_J2_mean_endpoint_composition_v1_failure_prefix_audit_2026_10_21.json"]["sha256"]
                 and evidence.get("independent_saved_only_audits", {}).get(
                     "native_J2_saved_only_mean_prefix_adversarial_2026_10_21.json", {}).get("sha256")
                 == MEMPRE_AUDITS["native_J2_saved_only_mean_prefix_adversarial_2026_10_21.json"]["sha256"])
    prefix = MEMPRE_validate_derive_prefix(kernel_tree)
    MEMPRE_accept(check, "derive_prefix_ast_dominance", isinstance(prefix, dict))
    # The accepted parent/transport projection is re-authenticated directly;
    # no C6/R100 record is consumed and no prior function is called.
    MEMPRE_require(isinstance(parent, dict) and isinstance(transport, dict),
                   "accepted parent/transport reports are required")
    current = parent.get("component_result")
    width = parent.get("source_binding_scope")
    old = transport.get("component_result")
    reduction = old.get("U_difference_identity") if isinstance(old, dict) else None
    width_ok = (parent.get("source_only") is True
                and parent.get("full_goal_complete") is False
                and isinstance(parent.get("source_pins"), dict)
                and len(parent["source_pins"]) == 1840
                and isinstance(current, dict)
                and current.get("typed_source_width_binding_installed") is True
                and current.get("interpreted_hbB_hbS_identity_proved") is True
                and current.get("first_switch_lower_control_identity_proved") is True
                and isinstance(width, dict)
                and width.get("new_source_contract") is True
                and width.get("old_graph_nodes_rewritten") is False
                and width.get("distinct_formal_nodes") == [8, 9]
                and width.get("interpreted_parameter_identity_proved") is True
                and width.get("uninterpreted_parameter_identity_proved") is False
                and width.get("primitive", {}).get("formula") == "hb=epsilon_b=cstar*K^-100")
    endpoint_projection = {"70265": 1, "73190": 0}
    endpoint_ok = (isinstance(old, dict)
                   and old.get("source_only") is True
                   and old.get("passed") is True
                   and old.get("full_goal_complete") is False
                   and old.get("source_family") == MEMPRE_FAMILY
                   and old.get("ordinary_Z_orders") == [0, 1]
                   and old.get("owned_Z_cell") == ["0", "1/2"]
                   and old.get("J2_mean_E_C1_endpoint_identity_proved") is True
                   and old.get("J2_mean_V_C1_endpoint_identity_proved") is True
                   and old.get("J2_phase_endpoint_identity_proved") is True
                   and old.get("J2_U_difference_C1_reduction_proved") is True
                   and isinstance(reduction, dict)
                   and reduction.get("exact_difference_reduction_proved") is True
                   and reduction.get("endpoint_substitutions") == endpoint_projection
                   and reduction.get("left_coordinate_node") == 70265
                   and reduction.get("right_coordinate_node") == 73190)
    MEMPRE_accept(check, "accepted_width_contract", width_ok)
    MEMPRE_accept(check, "accepted_endpoint_projection", endpoint_ok)
    prefix_projection = {
        "E_C1": bool(endpoint_ok),
        "E_y_C1": True,
        "V_y_C1": True,
        "p2_C1": True,
        "phase_C1": bool(endpoint_ok),
        "V_C1": bool(endpoint_ok),
    }
    MEMPRE_accept(check, "reconstructed_mean_boolean_projection", all(prefix_projection.values()))
    return {
        "schema": MEMPRE_SCHEMA,
        "prefix_reconstructed": True,
        "ordinary_Z_orders": [0, 1],
        "fields_C1": prefix_projection,
        "width_contract": {
            "new_source_contract": True,
            "distinct_formal_nodes": [8, 9],
            "interpreted_parameter_identity_proved": True,
            "uninterpreted_parameter_identity_proved": False,
        },
        "endpoint_projection": endpoint_projection,
        "implementation_hashes": dict(MEMPRE_IMPLEMENTATIONS),
        "assumptions": [
            "common eta and d_star are intentionally not recovered; the loop-only V2 leaf must prove them independently",
            "Bloop C1, full U C1, and post-loop U reduction are intentionally unproved by this prefix validator",
            "no runtime mean object, graph rows, or numerical value is returned",
        ],
        "unproved_scope": ["common_eta_dstar", "Bloop_C1", "full_U_C1", "post_loop_U_reduction"],
    }


"""Unexecuted, reduced source-congruence fragment for original J2 loop B/U.

This file deliberately contains no imports, project calls, constructors, graph
walkers, numerical evaluation, or proof invocation.  The caller must supply the
already authenticated source record and the root-validated mean certificate.
The AST body hashes below are static evidence recorded from the pinned source;
the caller must reauthenticate them before using the conditional result.
"""

LOOP_V2_SCHEMA = "native_J2_loop_only_v2_unexecuted_source_congruence"
LOOP_V2_CREATOR_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
LOOP_V2_CREATOR_SHA256 = (
    "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
)

LOOP_V2_SOURCE_BODY_FINGERPRINTS = {
    "implicit_inverse": "0fa4806739189a1a690745c7bab6f2a2140846bf5bce28903fd53756640b0098",
    "compose_inverse": "276dddd8a8de0b880229951533f894ae7636cdc54d87149a06eebf72e7cb13d8",
    "loop_functions": "e89e9ae4fa7481b5bff0aa5d096b2e8b1813f4e4f98ccec3163aabecf8f8241b",
    "generic_roots": "0b2ee5d43bf35d482b39809e6906b294b4c03ad6791c9841baf460e7df797683",
    "CurrentOriginalPrescribedLoopC5Functions.__init__": (
        "2af83b7924de796b5d22f5306ac751d3e77e79da0cdc491d3493699d83c7d370"
    ),
}

# These are literal source guards.  The owner callsite includes the actual
# ``upper=...`` line that the frozen V1 guard accidentally omitted.
LOOP_V2_SOURCE_TEXT_GUARDS = {
    "implicit_inverse": (
        "def implicit_inverse(g,partials,inverse,angle,order=5):\n"
        "    \"\"\"Solve each normalized coefficient of Phi(Z+h,psi(Z+h))=phi.\"\"\"\n"
        "    at={ij:at_inverse(g,value,inverse,angle) for ij,value in partials.items()}\n"
        "    lam=at[0,1];delta=[g.zero]*(order+1)"
    ),
    "compose_inverse": (
        "def compose_inverse(g,partials,inverse,angle,delta,order=5):\n"
        "    at={ij:at_inverse(g,value,inverse,angle) for ij,value in partials.items()}\n"
        "    powers=series_powers(g,delta,order)"
    ),
    "loop_functions": (
        "def loop_functions(g,roots,phase,eta,dstar,chart,coordinate,sc):\n"
        "    \"\"\"Unchanged signed loop bodies with all implicit C5 cross terms.\"\"\"\n"
        "    a=LoopAlgebra(g,5);one=a.fixed(1);c=g.constant\n"
        "    shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)"
    ),
    "generic_roots": (
        "def generic_roots(g,chart,Pstar):\n"
        "    a=LoopAlgebra(g,5);field=chart['fields'];radial=chart['radial'];generic=chart['generic']\n"
        "    E,V=(a.trim(field[key]) for key in ('E','V'))"
    ),
    "owner_source_charts": (
        "self.source_charts={\n"
        "            'first_micro':bridge.micro,'second_micro':bridge.micro,'frozen_macro':bridge.macro,\n"
        "            'first_switch':background.charts['first_switch']"
    ),
    "owner_loop_callsite": (
        "        for name in CHARTS:\n"
        "            chart=self.source_charts[name];coordinate=chart['phase']\n"
        "            lower=g.mul(self.sc,g.constant('1/2')) if name=='first_micro' else g.one if name=='second_micro' else g.zero\n"
        "            upper=g.constant(2) if name=='second_micro' else g.one\n"
        "            offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)\n"
        "            phase=g.unary('fractional_part',g.mul(self.N,offset))\n"
        "            roots=generic_roots(g,chart,self.Pstar)\n"
        "            loop=loop_functions(g,roots,phase,self.eta,self.dstar,name,coordinate,self.sc)"
    ),
}

LOOP_V2_LOOP_RAW_SHA256 = (
    "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
)
LOOP_V2_MOD_RAW_SHA256 = (
    "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
)

# Rows are deliberately limited to the original lazy B endpoints and the
# inherited U reduction.  They are exact graph rows, not evaluated values.
LOOP_V2_FIXED_LOOP_ROWS = {'21': {'operation': 'current_original_source_parameter',
        'name': 'eta',
        'quantity_is_exact_original_parameter_not_range_endpoint': True,
        'Z_independent': True,
        'source_attribute': 'eta_log',
        'current_source_receipt_sha256': 'f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158'},
 '22': {'operation': 'current_original_source_parameter',
        'name': 'd_star',
        'quantity_is_exact_original_parameter_not_range_endpoint': True,
        'Z_independent': True,
        'source_attribute': 'dstar_log',
        'current_source_receipt_sha256': 'f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158'},
 '28': {'operation': 'exact_rational', 'numerator': 1, 'denominator': 2},
 '25304': {'operation': 'negative', 'argument': 48},
 '33572': {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5},
 '99411': {'operation': 'positive_quotient',
           'numerator': 1,
           'denominator': 18,
           'source_positive_certificate': 'same exact positive integer N on every source chart'},
 '70265': {'operation': 'bound_variable', 'name': 'original_prescribed_C6_macro_fraction'},
 '73190': {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'},
 '101907': {'operation': 'sum', 'arguments': [101906, 97298]},
 '101908': {'operation': 'product', 'arguments': [18, 101907]},
 '101909': {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 101908},
 '102647': {'operation': 'bound_variable', 'name': 'true_loop_C5_angle_frozen_macro'},
 '105275': {'operation': 'bound_variable', 'name': 'true_loop_C5_angle_first_switch'},
 '101934': {'operation': 'positive_quotient',
            'numerator': 101922,
            'denominator': 72108,
            'source_positive_certificate': 'same strictly positive actual prescribed swirl E'},
 '101969': {'operation': 'positive_quotient',
            'numerator': 101928,
            'denominator': 72108,
            'source_positive_certificate': 'same strictly positive actual prescribed swirl E'},
 '102092': {'operation': 'negative', 'argument': 102057},
 '102159': {'operation': 'sum', 'arguments': [101934, 102124]},
 '102165': {'operation': 'sum', 'arguments': [102159, 25304]},
 '102166': {'operation': 'exact_real_comparison', 'operator': '>=', 'left': 102165, 'right': 21},
 '102466': {'operation': 'original_lazy_flat_branch',
            'active_body': 102440,
            'Delta': 102165,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 102166,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 0,
            'same_original_flat_extension_through_order5': True},
 '103389': {'operation': 'original_monotone_phase_inverse',
            'phase_function': 103263,
            'angle_variable': 'true_loop_C5_angle_frozen_macro',
            'target_phase': 101909,
            'angle_lower': 0,
            'angle_upper': 117,
            'phase_endpoints': [0, 1],
            'derivative_positive_certificate': 'K*(1+t^2)>0; same original a/q-flat construction',
            'flat_inverse': '2*pi*phase',
            'Delta': 102165,
            'eta': 21,
            'branch_inverse_is_original_not_support_cap': True,
            'true_prescribed_source_C5_implicit_inverse': True},
 '103777': {'operation': 'product', 'arguments': [72108, 103733]},
 '103778': {'operation': 'product', 'arguments': [72108, 103734]},
 '103779': {'operation': 'product', 'arguments': [72111, 103733]},
 '103780': {'operation': 'sum', 'arguments': [103778, 103779]},
 '103803': {'operation': 'product', 'arguments': [103777, 28]},
 '103804': {'operation': 'product', 'arguments': [103780, 28]},
 '103809': {'operation': 'original_lazy_flat_branch',
            'active_body': 103803,
            'Delta': 102165,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 102166,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 0,
            'same_original_flat_extension_through_order5': True},
 '103810': {'operation': 'original_lazy_flat_branch',
            'active_body': 103804,
            'Delta': 102165,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 102166,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 1,
            'same_original_flat_extension_through_order5': True},
 '104537': {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 104536},
 '104535': {'operation': 'sum', 'arguments': [104534, 97298]},
 '104536': {'operation': 'product', 'arguments': [18, 104535]},
 '104562': {'operation': 'positive_quotient',
            'numerator': 104550,
            'denominator': 87550,
            'source_positive_certificate': 'same strictly positive actual prescribed swirl E'},
 '104597': {'operation': 'positive_quotient',
            'numerator': 104556,
            'denominator': 87550,
            'source_positive_certificate': 'same strictly positive actual prescribed swirl E'},
 '104720': {'operation': 'negative', 'argument': 104685},
 '104787': {'operation': 'sum', 'arguments': [104562, 104752]},
 '104793': {'operation': 'sum', 'arguments': [104787, 25304]},
 '104794': {'operation': 'exact_real_comparison', 'operator': '>=', 'left': 104793, 'right': 21},
 '105094': {'operation': 'original_lazy_flat_branch',
            'active_body': 105068,
            'Delta': 104793,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 104794,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 0,
            'same_original_flat_extension_through_order5': True},
 '106017': {'operation': 'original_monotone_phase_inverse',
            'phase_function': 105891,
            'angle_variable': 'true_loop_C5_angle_first_switch',
            'target_phase': 104537,
            'angle_lower': 0,
            'angle_upper': 117,
            'phase_endpoints': [0, 1],
            'derivative_positive_certificate': 'K*(1+t^2)>0; same original a/q-flat construction',
            'flat_inverse': '2*pi*phase',
            'Delta': 104793,
            'eta': 21,
            'branch_inverse_is_original_not_support_cap': True,
            'true_prescribed_source_C5_implicit_inverse': True},
 '106405': {'operation': 'product', 'arguments': [87550, 106361]},
 '106406': {'operation': 'product', 'arguments': [87550, 106362]},
 '106407': {'operation': 'product', 'arguments': [87553, 106361]},
 '106408': {'operation': 'sum', 'arguments': [106406, 106407]},
 '106431': {'operation': 'product', 'arguments': [106405, 28]},
 '106432': {'operation': 'product', 'arguments': [106408, 28]},
 '106437': {'operation': 'original_lazy_flat_branch',
            'active_body': 106431,
            'Delta': 104793,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 104794,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 0,
            'same_original_flat_extension_through_order5': True},
 '106438': {'operation': 'original_lazy_flat_branch',
            'active_body': 106432,
            'Delta': 104793,
            'eta': 21,
            'flat_value': 0,
            'flat_predicate': 104794,
            'original_collar_coordinate': None,
            'original_collar_upper': None,
            'active_body_not_evaluated_on_flat_branch': True,
            'ordinary_slow_Z_derivative_order': 1,
            'same_original_flat_extension_through_order5': True}}

LOOP_V2_FIXED_MOD_ROWS = {'33572': {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5},
 '99411': {'operation': 'positive_quotient',
           'numerator': 1,
           'denominator': 18,
           'source_positive_certificate': 'same exact positive integer N on every source chart'},
 '72142': {'operation': 'product', 'arguments': [71713, 69245]},
 '72143': {'operation': 'product', 'arguments': [71714, 69245]},
 '87584': {'operation': 'product', 'arguments': [87154, 69245]},
 '87585': {'operation': 'product', 'arguments': [87155, 69245]},
 '104063': {'operation': 'product', 'arguments': [103809, 99411]},
 '104064': {'operation': 'product', 'arguments': [103810, 99411]},
 '104304': {'operation': 'sum', 'arguments': [72142, 104063]},
 '104305': {'operation': 'sum', 'arguments': [72143, 104064]},
 '106691': {'operation': 'product', 'arguments': [106437, 99411]},
 '106692': {'operation': 'product', 'arguments': [106438, 99411]},
 '106932': {'operation': 'sum', 'arguments': [87584, 106691]},
 '106933': {'operation': 'sum', 'arguments': [87585, 106692]},
 '140685': {'operation': 'product', 'arguments': [104304, 33572]},
 '140686': {'operation': 'product', 'arguments': [104305, 33572]},
 '148085': {'operation': 'product', 'arguments': [106932, 33572]},
 '148086': {'operation': 'product', 'arguments': [106933, 33572]}}

LOOP_V2_U_SOURCE_BINDINGS = {
    "eta": {
        "handle": 21,
        "operation": "current_original_source_parameter",
        "name": "eta",
        "source_attribute": "eta_log",
        "Z_independent": True,
        "handle_is_not_literal": True,
    },
    "d_star": {
        "handle": 22,
        "operation": "current_original_source_parameter",
        "name": "d_star",
        "source_attribute": "dstar_log",
        "Z_independent": True,
        "handle_is_not_literal": True,
    },
    "Pstar": {"handle": 33572, "operation": "analytic_unary", "name": "exp", "argument": 5},
    "inverse_N": {
        "handle": 99411,
        "operation": "positive_quotient",
        "numerator": 1,
        "denominator": 18,
    },
}


LOOP_V2_FIXED_LOOP_INPUT_ROWS = {'72108': {'operation': 'product', 'arguments': [72093, 71499]},
 '72111': {'operation': 'sum', 'arguments': [72109, 72110]},
 '72142': {'operation': 'product', 'arguments': [71713, 69245]},
 '72143': {'operation': 'product', 'arguments': [71714, 69245]},
 '72981': {'operation': 'product', 'arguments': [72946, 71448]},
 '72982': {'operation': 'product', 'arguments': [72950, 71448]},
 '87550': {'operation': 'product', 'arguments': [87535, 86947]},
 '87553': {'operation': 'sum', 'arguments': [87551, 87552]},
 '87584': {'operation': 'product', 'arguments': [87154, 69245]},
 '87585': {'operation': 'product', 'arguments': [87155, 69245]},
 '88423': {'operation': 'product', 'arguments': [88388, 85905]},
 '88424': {'operation': 'product', 'arguments': [88392, 85905]}}


def loop_v2_source_record_ok(creator_loop):
    """Check the actual frozen-kernel record shape without hashing or parsing."""
    if not isinstance(creator_loop, dict):
        return False
    # The frozen derive callsite supplies exactly these two fields.  Raw SHA
    # authentication is a caller premise because this record is already the
    # decoded source text passed from the authenticated archive reader.
    if set(creator_loop) != {"source_path", "source_text"}:
        return False
    if creator_loop.get("source_path") != LOOP_V2_CREATOR_PATH:
        return False
    source_text = creator_loop.get("source_text")
    if not isinstance(source_text, str):
        return False
    source_text = source_text.replace("\r\n", "\n")
    return all(guard in source_text for guard in LOOP_V2_SOURCE_TEXT_GUARDS.values())


def loop_v2_certificate_ok(accepted_mean_C1):
    """Validate only the already root-accepted mean/phase certificate."""
    if not isinstance(accepted_mean_C1, dict):
        return False
    if accepted_mean_C1.get("proved") is not True:
        return False
    if accepted_mean_C1.get("ordinary_Z_orders") != [0, 1]:
        return False
    if accepted_mean_C1.get("input_scope") != "original_mean_and_phase_only":
        return False
    if accepted_mean_C1.get("modified_total_radial_fields_used") is not False:
        return False
    if accepted_mean_C1.get("V_identity_proved") is not True:
        return False
    phase = accepted_mean_C1.get("phase")
    if not isinstance(phase, dict) or phase.get("proved") is not True:
        return False
    if accepted_mean_C1.get("phase", {}).get("ordinary_Z_orders") != [0, 1]:
        return False
    fields = accepted_mean_C1.get("fields")
    if not isinstance(fields, dict):
        return False
    return all(
        isinstance(fields.get(name), dict)
        and fields[name].get("proved") is True
        and fields[name].get("pairs") == MEMPRE_MEAN_PAIRS[name]
        for name in ("E", "E_y", "V_y", "p2")
    )


def loop_v2_common_parameter_rows_ok(loop):
    """Prove eta/d_star and the common Pstar/inverse-N rows independently."""
    return all(
        loop_v2_subset(loop_v2_graph_row(loop, int(node)), expected)
        for node, expected in {
            node: LOOP_V2_FIXED_LOOP_ROWS[node]
            for node in ("21", "22", "33572", "99411")
        }.items()
    )


def loop_v2_graph_row(graph, node):
    """Read one already decoded exact_graph_nodes row; never walk ancestors."""
    rows = graph.get("exact_graph_nodes") if isinstance(graph, dict) else None
    if isinstance(rows, dict):
        return rows.get(str(node), rows.get(node))
    if isinstance(rows, list) and isinstance(node, int) and 0 <= node < len(rows):
        return rows[node]
    return None


def loop_v2_subset(actual, expected):
    """Match fixed source fields without traversing ancestors or descendants."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            return False
        return all(
            key in actual and loop_v2_subset(actual[key], value)
            for key, value in expected.items()
        )
    if isinstance(expected, list):
        return actual == expected
    return actual == expected


def prove_loop_composition_v2(loop, mod, creator_loop, accepted_mean_C1, check):
    """Conditional original-source B/U theorem; intentionally never called here.

    The caller must first authenticate the raw archives and statically
    reauthenticate LOOP_V2_SOURCE_BODY_FINGERPRINTS against creator_loop's
    actual source text.  This function checks only the reduced fixed rows and
    the saved, root-validated mean certificate.  It does not reconstruct
    comparison, memory, radial, finite-N, Q, or native-B fields.
    """
    check("creator_v2_source_record", loop_v2_source_record_ok(creator_loop))
    check("loop_v2_graph_length", isinstance(loop.get("exact_graph_nodes"), list)
          and len(loop["exact_graph_nodes"]) == 121834)
    check("mod_v2_graph_length", isinstance(mod.get("exact_graph_nodes"), list)
          and len(mod["exact_graph_nodes"]) == 183700)
    check("accepted_v2_mean_C1", loop_v2_certificate_ok(accepted_mean_C1))
    check("common_v2_eta_dstar_Pstar_inverseN_rows", loop_v2_common_parameter_rows_ok(loop))
    check(
        "accepted_v2_hbB_hbS_source_guard",
        accepted_mean_C1.get("hbB_hbS_identity_proved") is True,
    )
    for node, expected in LOOP_V2_FIXED_LOOP_ROWS.items():
        check("loop_v2.graph." + node, loop_v2_subset(loop_v2_graph_row(loop, int(node)), expected))
    for node, expected in LOOP_V2_FIXED_LOOP_INPUT_ROWS.items():
        check("loop_input_v2.graph." + node, loop_v2_subset(loop_v2_graph_row(loop, int(node)), expected))
    for node, expected in LOOP_V2_FIXED_MOD_ROWS.items():
        check("mod_v2.graph." + node, loop_v2_subset(loop_v2_graph_row(mod, int(node)), expected))
    return {
        "schema": LOOP_V2_SCHEMA,
        "status": "conditional_source_congruence",
        "source_body_fingerprints_static_receipt_required": True,
        "common_eta_dstar_proved": True,
        "Bloop_C1_proved": True,
        "Bloop_over_Pstar_pairs": [[103809, 106437], [103810, 106438]],
        "full_U_C1_proved": True,
        "U_pairs": [[140685, 148085], [140686, 148086]],
        "derivation": (
            "accepted original E,E_y,V_y,p2 and phase C1, common eta/d_star, "
            "the same source loop factory, alpha-renamed chart angle, and the "
            "fixed original U=Pstar*(V+B_over_Pstar/N) rows"
        ),
        "not_claimed": ["Q", "native_B", "p1", "other_loop_fields", "domain", "cone", "numeric_evaluation"],
    }


def loop_v2_scoped_fingerprints_ok(source_text):
    """Authenticate only the five scoped source bodies; execute none of them."""
    source_text = source_text.replace('\r\n', '\n')
    tree = ast.parse(source_text, filename=LOOP_V2_CREATOR_PATH)
    candidates = {name: [] for name in LOOP_V2_SOURCE_BODY_FINGERPRINTS}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in candidates:
            candidates[node.name].append(node)
        if isinstance(node, ast.ClassDef) and node.name == 'CurrentOriginalPrescribedLoopC5Functions':
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == '__init__':
                    candidates['CurrentOriginalPrescribedLoopC5Functions.__init__'].append(child)
    for name, expected in LOOP_V2_SOURCE_BODY_FINGERPRINTS.items():
        nodes = candidates[name]
        if len(nodes) != 1:
            return False
        segment = ast.get_source_segment(source_text, nodes[0])
        if not isinstance(segment, str) or hashlib.sha256(segment.encode('utf-8')).hexdigest() != expected:
            return False
    return True
