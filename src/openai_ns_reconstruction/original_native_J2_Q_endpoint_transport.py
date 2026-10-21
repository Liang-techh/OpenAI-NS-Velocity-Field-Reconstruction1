"""Source-only transport of physical Q at the macro/first-switch seam.

This new leaf inherits accepted source identities without replay. It checks
the finite-N m memory handoff through ordinary Z order 2 and then the actual
physical Q quotient through order 1. No numerical or physical admission.
"""
import ast
import base64
import hashlib
import json

SCHEMA = 'original_native_J2_Q_endpoint_transport_v1'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_mean_endpoint_composition_v2.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
MOD = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
MOD_CREATOR = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
EVIDENCE = 'experiments/root_st073/original_native_J2_mean_endpoint_v1_failure_evidence.json'
RECORD_NAMES = [PARENT, LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, EVIDENCE]
RECORD_HASHES = {
    PARENT: 'c0d3021331fa967f4504cf915287e95af401cb8a38bb544ed0ef8ed3a3cbcd8b',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP_CREATOR: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
    MOD_CREATOR: '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7',
    EVIDENCE: 'fe8d0169c01d7023ada2e4644f59f1ccfbf7978e67542f4e639e846d725c9cd1',
}
STAGE_FLAGS = ['actual_J2_Q_ports_bound', 'original_m_Z2_premise_authenticated',
               'finite_N_m_correction_Z2_endpoint_identity_proved',
               'physical_m_Z2_endpoint_identity_proved',
               'J2_Q_C1_endpoint_identity_proved',
               'inherited_R_U_identities_retained', 'all_old_terms_retained',
               'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['J2_B_identity_proved', 'all_native_B_jumps_zero_proved',
               'J2_Qy_identity_proved', 'numeric_values_computed',
               'source_C5_installation_claimed', 'numeric_containment_claimed',
               'new_loop_domain_admission_claimed', 'admissible_cone_admission_claimed',
               'old_error_terms_removed', 'global_five_moment_closure',
               'time_order_recursion_claimed', 'full_goal_complete',
               'uninterpreted_graph_parameter_identity_proved', 'old_operators_replayed']
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
    parent, loop, mod, evidence = [decoded[name] for name in (PARENT, LOOP, MOD, EVIDENCE)]
    old = parent['component_result']
    check('accepted_parent_scope', parent['source_only'] is True
          and old['passed'] is True and parent['full_goal_complete'] is False
          and len(parent['source_pins']) == 1844
          and parent['inherited_R_identity_jumps'] == [2, 3, 4, 6, 8, 9, 10]
          and parent['B_identity_scope']['unresolved_B_identity_nodes'] == list(range(522, 536))
          and parent['Q_identity_scope']['J2_Q_identity_proved'] is False
          and old['actual_J2_U_C1_endpoint_identity_proved'] is True
          and old['original_memory_Z2_endpoint_identity_proved'] is True
          and old['ordinary_Z_orders'] == [0, 1] and old['owned_Z_cell'] == ['0', '1/2']
          and old['source_function_identity_only'] is True
          and old['numeric_values_computed'] is False
          and old['U_identity_scope']['left_source_nodes'] == [140685, 140686]
          and old['U_identity_scope']['right_source_nodes'] == [148085, 148086]
          and old['U_identity_scope']['endpoint_substitutions'] == {'70265': 1, '73190': 0})
    check('inherited_input_pins', all(parent['source_pins'][name] == RECORD_HASHES[name]
          for name in (LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, EVIDENCE)))
    check('same_source_family', loop['source_family'] == mod['source_family'] == old['source_family'] == FAMILY)
    check('unadmitted_source_factory_flags_retained',
          loop['current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed'] is False
          and loop['current_original_prescribed_loop_C5_numeric_containment_installed'] is False)
    # New bounded proof fragments are integrated only after static review.
    original_m = original_m_Z2_certificate(parent, evidence,
        lambda name, value: check('original_m:' + name, value))
    transported = prove_physical_Q_transport(loop, mod, decoded[LOOP_CREATOR],
        decoded[MOD_CREATOR], original_m, old,
        lambda name, value: check('Q_transport:' + name, value))
    check('fresh_m_and_Q_source_identities', transported['finite_N_m_Z2_proved'] is True
          and transported['physical_m_Z2_proved'] is True and transported['Q_C1_proved'] is True)
    retained = parent['retained_construction_data']
    check('retained_ledger', all(retained[key] == value for key, value in
          dict(cumulative_terms=212, raw_signed_terms=213, reduced_terms=183,
               six_signed_channels_retained=True, all_old_error_data_retained=True,
               reserved_contribution_registered_once=True, inherited_errors_deleted=False).items()))
    proof = [dict(id='original_m_Z2', rule='authenticate accepted saved-source ordinary-Z2 m binding'),
             dict(id='finite_N_m_Z2', rule='same predecessor incoming memory, decay one, equal-bound integrals zero through Z2',
                  depends_on=['original_m_Z2']),
             dict(id='physical_m_Z2', rule='original m plus transported finite-N correction; physical Pstar applied once',
                  depends_on=['finite_N_m_Z2']),
             dict(id='physical_Q_C1', rule='actual signed numerator and differentiated common-L quotient',
                  depends_on=['physical_m_Z2'])]
    return dict(schema=SCHEMA, component='Q_endpoint_transport', passed=True, source_only=True,
        **dict.fromkeys(STAGE_FLAGS, True), **dict.fromkeys(FALSE_FLAGS, False),
        stage_flags=dict.fromkeys(STAGE_FLAGS, True), nonclaims=dict.fromkeys(FALSE_FLAGS, False),
        J2_Q_identity_proved=True, ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        m_ordinary_Z_orders=[0, 1, 2], owned_Z_cell=['0', '1/2'],
        input_digest=expected_input_digest, source_family=FAMILY,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_function_identity_only=True, source_checks=checks,
        proof_graph=proof, proof_graph_sha256=digest(proof),
        component_proofs=dict(original_m=original_m, Q_transport=transported),
        source_binding_scope=parent['source_binding_scope'],
        retained_construction_data=retained, quantity_distinctions=parent['quantity_distinctions'],
        inherited_R_identity_jumps=[2, 3, 4, 6, 8, 9, 10],
        unresolved_B_identity_nodes=list(range(522, 536)),
        Q_identity_scope=dict(weighted_jump_node=522, jump_index=2, full_Q_C1_identity_proved=True,
            ordinary_Z_orders=[0, 1], left_source_nodes=[140780, 140784],
            right_source_nodes=[148180, 148184], endpoint_substitutions={'70265': 1, '73190': 0},
            common_L_denominator_node=25272, source_interpretation_required=True),
        inherited_U_identity_scope=parent['U_identity_scope'],
        unresolved_gates=[], remaining_construction_gates=[
            'native weighted J2 B source identity and the other thirteen native B jumps',
            'independent terminal functional moments and retained quantitative errors',
            'heat exterior, stress/remainder, temporal n recursion, realizable corrections and Cartesian NS validation'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_proof_functions_called=0, generic_DAG_walks=0,
            numerical_evaluations=0, quadratures=0, recurrence_calls=0, runtime_model_metadata=None),
        runtime_model_metadata=None)


"""Unexecuted original J2 m ordinary-Z endpoint projection.

This is a data-only adapter for an accepted parent certificate.  It does not
rerun the V1 proof, recover an intermediate object, read an archive, execute a
source function, or select a numerical source value.  It authenticates the
embedded V1 bytes and checks the frozen V1 AST shape before returning the
already-reviewed typed m-row projection.
"""
import ast as MZ_ast
import base64 as MZ_base64
import hashlib as MZ_hashlib
import json as MZ_json


MZ_SCHEMA = "native_J2_original_m_Z2_endpoint_certificate_v1"
MZ_EVIDENCE_SCHEMA = "original_native_J2_mean_endpoint_v1_failure_evidence_v1"
MZ_SOURCE_SCHEMA = "native_J2_mean_endpoint_composition_inputs_v1"
MZ_EVIDENCE_SHA256 = "fe8d0169c01d7023ada2e4644f59f1ccfbf7978e67542f4e639e846d725c9cd1"
MZ_FROZEN_KERNEL_SHA256 = "b10b2600f1f1a3b4d0f387e66e412da668b666dcd8996a71305caaa071e3f65f"
MZ_FROZEN_RUNNER_SHA256 = "7581781d037b9d6b59a0cd86934d3521712e3588c17685d93bfbc2affb337400"
MZ_SOURCE_CHECKPOINT_SHA256 = "85421129e9d7ac7a3f9a6288f8bb48c184e738cc9101d7c4dad3d3aa274d45c0"
MZ_ENTERED_CHECKPOINT_SHA256 = "439bbacb222c1da88b84858147b6a6df4012d3a3c0ee4b9473dc1ced85514981"
MZ_RESULT_CHECKPOINT_SHA256 = "483ec74f99290a5f77fd79de0ce896f8a0fecea70dcc8cf8f493dbdc8f9fa1f2"
MZ_TRACEBACK_SHA256 = "71595f51f948aada2ad1438fdc787315e93c73798b2229db0123530e8289b573"
MZ_INPUT_DIGEST = "8873772c8b689bd6f49e7ccf50ae8bd02761bb48646e1e72ec0ac9961794f025"
MZ_PARENT_RECORD = "lei_ren_part1_paper_compliant_current_original_native_J2_mean_loop_input_join.json"
MZ_TRANSPORT_RECORD = "lei_ren_part1_paper_compliant_current_original_native_J2_U_endpoint_transport_v2.json"
MZ_PARENT_RECORD_SHA256 = "7143e443ec041b4d94d4a80a42dc9107dcee6d36284708661e7f25905759b1e9"
MZ_TRANSPORT_RECORD_SHA256 = "be0cad628f44cec0713144f074cbbe4c7b0efc23feab9d3d482a116235d3737b"

MZ_IMPLEMENTATIONS = {
    "scripts/run_original_native_J2_mean_endpoint_composition.py": MZ_FROZEN_RUNNER_SHA256,
    "src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py": MZ_FROZEN_KERNEL_SHA256,
}
MZ_ARTIFACTS = {
    "scripts/run_original_native_J2_mean_endpoint_composition.py": {
        "sha256": MZ_FROZEN_RUNNER_SHA256,
        "bytes": 11878,
    },
    "src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py": {
        "sha256": MZ_FROZEN_KERNEL_SHA256,
        "bytes": 112579,
    },
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_source.json": {
        "sha256": MZ_SOURCE_CHECKPOINT_SHA256,
        "bytes": 273002,
    },
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_entered.json": {
        "sha256": MZ_ENTERED_CHECKPOINT_SHA256,
        "bytes": 683,
    },
    "tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_result.json": {
        "sha256": MZ_RESULT_CHECKPOINT_SHA256,
        "bytes": 3206,
    },
}
MZ_PREFIX_CHECKS = (
    "saved_prefix:embedded_source_schema",
    "saved_prefix:embedded_input_digest",
    "saved_prefix:entered_source_link",
    "saved_prefix:entered_input_link",
    "saved_prefix:entered_kernel_link",
    "saved_prefix:embedded_implementation_hashes",
    "saved_prefix:entered_implementation_hashes",
    "saved_prefix:result_source_link",
    "saved_prefix:result_entry_link",
    "saved_prefix:result_failure",
    "saved_prefix:exact_failure_traceback",
    "saved_prefix:embedded_parent_transport_hashes",
    "saved_prefix:embedded_audit_pins",
    "saved_prefix:derive_prefix_ast_dominance",
    "saved_prefix:accepted_width_contract",
    "saved_prefix:accepted_endpoint_projection",
    "saved_prefix:reconstructed_mean_boolean_projection",
)
MZ_LEFT_M_ROWS = [72149, 72150, 72151]
MZ_RIGHT_M_ROWS = [87591, 87592, 87593]
MZ_PROJECTION_ROWS = [73099, 73100, 73101]
MZ_ORDINARY_Z = [0, 1, 2]
MZ_PSTAR_NODE = 33572
MZ_PSTAR_INVERSE_NODE = 69245


def MZ_require(condition, message):
    if not condition:
        raise AssertionError(message)


def MZ_accept(check, name, condition):
    MZ_require(type(condition) is bool, name + " must be bool")
    MZ_require(check(name, condition) is True, "projection premise rejected: " + name)


def MZ_sha256(raw):
    return MZ_hashlib.sha256(raw).hexdigest()


def MZ_decode_artifact(item, expected, name):
    MZ_require(isinstance(item, dict), name + " artifact metadata changed")
    MZ_require(item.get("encoding") == "base64", name + " encoding changed")
    MZ_require(item.get("sha256") == expected["sha256"], name + " metadata hash changed")
    MZ_require(item.get("bytes") == expected["bytes"], name + " metadata byte count changed")
    encoded = item.get("base64_bytes")
    MZ_require(isinstance(encoded, str), name + " base64 payload missing")
    raw = MZ_base64.b64decode(encoded, validate=True)
    MZ_require(len(raw) == expected["bytes"], name + " decoded byte count changed")
    MZ_require(MZ_sha256(raw) == expected["sha256"], name + " embedded raw bytes changed")
    return raw


def MZ_const(node, value):
    return isinstance(node, MZ_ast.Constant) and node.value == value


def MZ_literal_sequence(node, values):
    return (
        isinstance(node, (MZ_ast.List, MZ_ast.Tuple))
        and len(node.elts) == len(values)
        and all(MZ_const(item, value) for item, value in zip(node.elts, values))
    )


def MZ_top_function(tree, name):
    found = [
        node for node in tree.body
        if isinstance(node, (MZ_ast.FunctionDef, MZ_ast.AsyncFunctionDef))
        and node.name == name
    ]
    MZ_require(len(found) == 1, "frozen V1 function scope changed: " + name)
    return found[0]


def MZ_direct_assignment(function, target):
    found = []
    for node in function.body:
        if not isinstance(node, MZ_ast.Assign) or len(node.targets) != 1:
            continue
        if isinstance(node.targets[0], MZ_ast.Name) and node.targets[0].id == target:
            found.append(node)
    MZ_require(len(found) == 1, "frozen V1 assignment scope changed: " + target)
    return found[0]


def MZ_return_mapping(function):
    found = [node for node in function.body if isinstance(node, MZ_ast.Return)]
    MZ_require(len(found) == 1 and isinstance(found[0].value, MZ_ast.Dict),
               "frozen V1 prove_memories return scope changed")
    keys = {}
    for key, value in zip(found[0].value.keys, found[0].value.values):
        MZ_require(isinstance(key, MZ_ast.Constant) and isinstance(key.value, str),
                   "frozen V1 return key is not a literal")
        MZ_require(key.value not in keys, "duplicate frozen V1 return key")
        keys[key.value] = value
    return keys


def MZ_validate_memory_ast(kernel_text):
    tree = MZ_ast.parse(kernel_text, filename="original_native_J2_mean_endpoint_composition.py")
    memory = MZ_top_function(tree, "prove_memories")

    z_orders = MZ_direct_assignment(memory, "z_orders")
    MZ_require(
        isinstance(z_orders.value, MZ_ast.Call)
        and isinstance(z_orders.value.func, MZ_ast.Name)
        and z_orders.value.func.id == "range"
        and len(z_orders.value.args) == 1
        and MZ_const(z_orders.value.args[0], 3),
        "original memory Z scope is not range(3)",
    )

    endpoint_fields = MZ_direct_assignment(memory, "endpoint_fields")
    MZ_require(MZ_literal_sequence(endpoint_fields.value, ["m", "e", "p", "pressure"]),
               "original memory m endpoint field binding changed")

    memory_assignment = MZ_direct_assignment(memory, "memory_Z2_proved")
    MZ_require(
        isinstance(memory_assignment.value, MZ_ast.Name)
        and memory_assignment.value.id == "internal_source_checks",
        "original memory Z2 return binding changed",
    )
    p2_assignment = MZ_direct_assignment(memory, "p2_C1_proved")
    MZ_require(
        isinstance(p2_assignment.value, MZ_ast.Name)
        and p2_assignment.value.id == "internal_source_checks",
        "original p2 return binding changed",
    )

    returned = MZ_return_mapping(memory)
    MZ_require(MZ_literal_sequence(returned["ordinary_Z_orders"], MZ_ORDINARY_Z),
               "original memory ordinary Z return scope changed")
    MZ_require(isinstance(returned["actual_field_rows"], MZ_ast.Dict),
               "original memory actual field row return changed")
    returned_text = MZ_ast.unparse(returned["actual_field_rows"])
    MZ_require("c6_macro_fields[field]" in returned_text
               and "c6_r100_fields[field]" in returned_text
               and "r_fields[field]" in returned_text,
               "original m endpoint return bindings changed")
    MZ_require(MZ_const(returned["source_only"], True),
               "original memory source-only scope changed")
    MZ_require(MZ_const(returned["numeric_values_computed"], False),
               "original memory numeric-value scope changed")
    return True


def original_m_Z2_certificate(parent, evidence, check):
    """Project the accepted original m endpoint certificate without reentry."""
    MZ_require(callable(check), "check must record boolean premises")
    MZ_require(isinstance(parent, dict), "accepted parent report is required")
    MZ_require(isinstance(evidence, dict), "frozen V1 evidence is required")

    parent_component = parent.get("component_result", parent)
    MZ_require(isinstance(parent_component, dict), "accepted parent component is missing")
    MZ_accept(check, "parent_original_memory_Z2_endpoint_identity_proved",
              parent_component.get("original_memory_Z2_endpoint_identity_proved") is True)
    MZ_accept(check, "parent_source_prefix_certificate_recovered",
              parent_component.get("source_prefix_certificate_recovered") is True)
    MZ_accept(check, "parent_source_only",
              parent_component.get("source_only") is True
              and parent_component.get("full_goal_complete") is False)
    saved_checks = parent_component.get("source_checks")
    MZ_require(isinstance(saved_checks, dict), "accepted saved-prefix checks are missing")
    MZ_accept(check, "saved_prefix_check_inventory",
              set(saved_checks) >= set(MZ_PREFIX_CHECKS)
              and len([name for name in MZ_PREFIX_CHECKS if saved_checks.get(name) is True]) == 17)
    MZ_accept(check, "saved_prefix_checks_all_true",
              all(saved_checks.get(name) is True for name in MZ_PREFIX_CHECKS))

    MZ_require(evidence.get("schema") == MZ_EVIDENCE_SCHEMA, "V1 evidence schema changed")
    MZ_require(evidence.get("V1_do_not_reenter") is True, "V1 reentry flag changed")
    MZ_require(evidence.get("V1_stage_passed") is False, "V1 failure flag changed")
    MZ_require(evidence.get("metadata_assembly_only") is True, "V1 evidence mode changed")
    MZ_require(evidence.get("source_functions_called") == 0, "V1 evidence records source execution")
    artifacts = evidence.get("artifacts")
    MZ_require(isinstance(artifacts, dict) and set(artifacts) == set(MZ_ARTIFACTS),
               "V1 artifact inventory changed")
    raw = {
        name: MZ_decode_artifact(artifacts[name], expected, name)
        for name, expected in MZ_ARTIFACTS.items()
    }
    source_doc = MZ_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_source.json"].decode("utf-8"))
    entered_doc = MZ_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_entered.json"].decode("utf-8"))
    result_doc = MZ_json.loads(raw["tmp/original_native_J2_mean_endpoint_composition/new_mean_endpoint_composition_result.json"].decode("utf-8"))
    kernel_text = raw["src/openai_ns_reconstruction/original_native_J2_mean_endpoint_composition.py"].decode("utf-8")

    MZ_accept(check, "embedded_source_schema", source_doc.get("schema") == MZ_SOURCE_SCHEMA)
    MZ_accept(check, "embedded_source_input_digest", source_doc.get("input_digest") == MZ_INPUT_DIGEST)
    MZ_accept(check, "embedded_source_implementation_hashes",
              source_doc.get("implementation_hashes") == MZ_IMPLEMENTATIONS)
    MZ_accept(check, "embedded_record_hashes",
              source_doc.get("record_hashes", {}).get(MZ_PARENT_RECORD) == MZ_PARENT_RECORD_SHA256
              and source_doc.get("record_hashes", {}).get(MZ_TRANSPORT_RECORD) == MZ_TRANSPORT_RECORD_SHA256)
    MZ_accept(check, "entered_checkpoint_linkage",
              entered_doc.get("source_sha256") == MZ_SOURCE_CHECKPOINT_SHA256
              and entered_doc.get("input_digest") == MZ_INPUT_DIGEST
              and entered_doc.get("kernel_sha256") == MZ_FROZEN_KERNEL_SHA256
              and entered_doc.get("implementation_hashes") == MZ_IMPLEMENTATIONS)
    traceback = result_doc.get("traceback")
    MZ_accept(check, "result_failure_linkage",
              result_doc.get("source_checkpoint_sha256") == MZ_SOURCE_CHECKPOINT_SHA256
              and result_doc.get("entry_checkpoint_sha256") == MZ_ENTERED_CHECKPOINT_SHA256
              and result_doc.get("passed") is False
              and result_doc.get("failed_attempt_preserved") is True
              and result_doc.get("error_type") == "AssertionError"
              and result_doc.get("error") == "Source prerequisite: loop:creator.authenticated_source_text_guards"
              and result_doc.get("runner_caught_exception") is True
              and result_doc.get("runner_full_goal_complete") is False)
    MZ_accept(check, "result_traceback_linkage",
              isinstance(traceback, str)
              and MZ_sha256(traceback.encode("utf-8")) == MZ_TRACEBACK_SHA256
              and 'original_native_J2_mean_endpoint_composition.py", line 204, in derive' in traceback
              and 'original_native_J2_mean_endpoint_composition.py", line 1381, in prove_loop_composition' in traceback)
    MZ_accept(check, "frozen_memory_ast_scope", MZ_validate_memory_ast(kernel_text))

    return {
        "source_only": True,
        "memory_Z2_proved": True,
        "ordinary_Z_orders": list(MZ_ORDINARY_Z),
        "left_m_rows": list(MZ_LEFT_M_ROWS),
        "right_m_rows": list(MZ_RIGHT_M_ROWS),
        "projection_rows": list(MZ_PROJECTION_ROWS),
        "common_Pstar_node": MZ_PSTAR_NODE,
        "Pstar_inverse_node": MZ_PSTAR_INVERSE_NODE,
        "typed_row_registry_requires_archive_authentication": True,
        "recovered_runtime_objects": False,
        "numeric_values_computed": False,
    }


"""Unexecuted source leaf for the physical J2 Q seam.

The caller supplies decoded, already authenticated archives and the accepted
original-m Z2 certificate.  This module only contains source-text and fixed
row guards; it is deliberately not invoked here.  The endpoint argument is
the frozen-macro upper point 70265=1 and the first-switch lower point
73190=0.  The m transport is checked through ordinary Z order 2 because the
physical Q formula consumes m, m_Z, and m_ZZ in its C1 rows.
"""

import ast as _ast
import hashlib as _hashlib


Q_FRAGMENT_SCHEMA = "native_J2_Q_numerator_transport_fragment_v1"
LOOP_CREATOR_SHA256 = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
MOD_CREATOR_SHA256 = "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"
MOD_RAW_GZIP_SHA256 = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
LOOP_RAW_GZIP_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
MOD_DECODED_SHA256 = "0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75"
LOOP_DECODED_SHA256 = "f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2"


def prove_physical_Q_transport(loop, mod, loop_creator_text, mod_creator_text,
                               original_m_certificate, accepted_parent_component,
                               check):
    """Build the unexecuted source proof result for physical Q.

    ``check`` is the caller's recording assertion.  No source producer,
    ancestor, graph walker, evaluator, or quadrature is called by this leaf.
    It reads only the supplied decoded records and parses the two supplied
    creator texts.
    """
    checks = {}

    def guard(name, condition):
        value = bool(condition)
        check(name, value)
        checks[name] = value
        return value

    def node(rows, number):
        if not isinstance(rows, list) or number < 0 or number >= len(rows):
            return None
        value = rows[number]
        return value if isinstance(value, dict) else None

    def row_is(rows, number, **expected):
        actual = node(rows, number)
        return actual is not None and all(actual.get(key) == value
                                         for key, value in expected.items())

    def rows_are(rows, numbers):
        return isinstance(numbers, list) and all(
            isinstance(number, int) and node(rows, number) is not None
            for number in numbers)

    def compact(text):
        return ''.join(str(text).split())

    loop_rows = loop.get('exact_graph_nodes', []) if isinstance(loop, dict) else []
    mod_rows = mod.get('exact_graph_nodes', []) if isinstance(mod, dict) else []
    loop_bind = loop.get('prescribed_loop_C5_bindings', {}) if isinstance(loop, dict) else {}
    mod_bind = mod.get('modified_Q4_pressure_C3_bindings', {}) if isinstance(mod, dict) else {}

    guard('loop_graph_length', len(loop_rows) == 121834)
    guard('mod_graph_length', len(mod_rows) == 183700)

    parent_ok = (
        isinstance(accepted_parent_component, dict)
        and accepted_parent_component.get('passed') is True
        and accepted_parent_component.get('source_only') is True
        and accepted_parent_component.get('full_goal_complete') is False
        and accepted_parent_component.get('numeric_values_computed') is False
        and accepted_parent_component.get('J2_Q_identity_proved') is False
        and accepted_parent_component.get('original_memory_Z2_endpoint_identity_proved') is True
        and accepted_parent_component.get('actual_J2_U_C1_endpoint_identity_proved') is True
        and 2 in accepted_parent_component.get('inherited_R_identity_jumps', [])
        and accepted_parent_component.get('ordinary_Z_orders') == [0, 1]
    )
    guard('accepted_parent_R_U_and_open_Q_scope', parent_ok)

    family_ok = (
        isinstance(loop, dict) and isinstance(mod, dict)
        and loop.get('source_family') == mod.get('source_family')
        and loop.get('source_family') == accepted_parent_component.get('source_family')
    )
    guard('shared_actual_source_family', family_ok)

    cert = original_m_certificate if isinstance(original_m_certificate, dict) else {}
    cert_memory = cert.get('memory_Z2_proved') is True or cert.get('original_m_Z2_proved') is True
    cert_rows = (
        cert.get('left_m_rows', []) == [72149, 72150, 72151]
        and cert.get('right_m_rows', []) == [87591, 87592, 87593]
    )
    cert_scope = (
        cert.get('source_only') is True
        and cert_memory
        and cert.get('numeric_values_computed') is False
        and cert.get('recovered_runtime_objects', False) is False
        and cert.get('ordinary_Z_orders') == [0, 1, 2]
        and cert.get('common_Pstar_node') == 33572
        and cert.get('Pstar_inverse_node') == 69245
    )
    guard('accepted_original_m_Z2_certificate', cert_scope and cert_rows)

    try:
        loop_tree = _ast.parse(loop_creator_text)
        mod_tree = _ast.parse(mod_creator_text)
        loop_parse = True
        mod_parse = True
    except (SyntaxError, TypeError, ValueError):
        loop_tree = mod_tree = None
        loop_parse = mod_parse = False
    loop_text = compact(loop_creator_text)
    mod_text = compact(mod_creator_text)
    guard('authenticated_loop_creator_text',
          loop_parse and _hashlib.sha256(loop_creator_text.encode()).hexdigest() == LOOP_CREATOR_SHA256)
    guard('authenticated_modified_Q_creator_text',
          mod_parse and _hashlib.sha256(mod_creator_text.encode()).hexdigest() == MOD_CREATOR_SHA256)

    loop_defs = [] if loop_tree is None else [n for n in _ast.walk(loop_tree)
        if isinstance(n, (_ast.FunctionDef, _ast.AsyncFunctionDef)) and n.name == 'transport']
    mod_defs = [] if mod_tree is None else [n for n in _ast.walk(mod_tree)
        if isinstance(n, (_ast.FunctionDef, _ast.AsyncFunctionDef)) and n.name == 'physical_pressure_fields']
    transport_ast = (
        len(loop_defs) == 1
        and 'width=g.mul(jacobian,g.sub(coordinate,lower))' in loop_text
        and "decay=g.unary('exp',g.mul(g.constant(-rate),width))ifrateelseg.one" in loop_text
        and 'contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True,' in loop_text
        and 'physical_Jacobian_applied_once=True' in loop_text
        and 'predecessor_memory_never_reset=True' in loop_text
        and 'out[name]=a.add(a.scale(incoming[name],decay),contribution)' in loop_text
        and "lower=g.mul(self.sc,g.constant('1/2'))ifname=='first_micro'elseg.oneifname=='second_micro'elseg.zero" in loop_text
        and 'incoming=self.endpoints[name]' in loop_text
        and 'fields[key]=a.add(a.trim(chart[\'fields\'][key]),history[key])' in loop_text
        and 'self.endpoints[name]={key:a.at(jet,coordinate,upper)forkey,jetinhistory.items()}' in loop_text
    )
    guard('loop_transport_and_predecessor_ast', transport_ast)

    q_ast = (
        len(mod_defs) == 1
        and 'a=Algebra(g,5);b=Algebra(g,4)' in mod_text
        and "U=a.scale(fields['V'],Pstar)" in mod_text
        and "m=a.scale(fields['m'],Pstar)" in mod_text
        and 'd=b.sub(one,b.mul(z,z));L=b.sub(one,b.mul(de,b.mul(z,z)))' in mod_text
        and "sameoriginalL=1-delta*Z^2" in mod_text
        and 'numerator=lambdau,mean:b.sub(b.scale(b.mul(z,b.trim(u)),2),b.add(b.mul(b.sub(one,de),b.mul(z,b.trim(mean))),b.mul(d,b.shift(mean))))' in mod_text
        and 'Q=div(numerator(U,m));Qy=div(numerator(Uy,my))' in mod_text
        and "chart=loop_owner.charts[name];bg=chart['source_chart'];coordinate=chart['coordinate'];jac=chart['physical_measure']" in mod_text
        and "physical=physical_pressure_fields(g,chart['modified_fields'],radial,self.background.delta,self.background.Pstar)" in mod_text
        and "returndict(U=U,U_y=Uy,m=m,m_y=my,Q=Q,Q_y=Qy)" in mod_text
    )
    guard('physical_Q_creator_ast', q_ast)

    try:
        lcharts = loop_bind['charts']
        macro = lcharts['frozen_macro']
        switch = lcharts['first_switch']
        mendpoints = loop_bind['endpoints']
        mchart = mod_bind['charts']['frozen_macro']
        schart = mod_bind['charts']['first_switch']
        loop_binding_shape = (
            macro['coordinate'] == 70265 and macro['domain'] == [0, 1] and macro['radius'] == 71448
            and switch['coordinate'] == 73190 and switch['domain'] == [0, 1] and switch['radius'] == 85905
            and macro['incoming']['m'][:3] == [101875, 101876, 101877]
            and macro['history']['m'][:3] == [104342, 104343, 104344]
            and macro['modified_fields']['m'][:3] == [104467, 104468, 104469]
            and mendpoints['frozen_macro']['m'][:3] == [104503, 104504, 104505]
            and switch['incoming']['m'][:3] == [104503, 104504, 104505]
            and switch['history']['m'][:3] == [106970, 106971, 106972]
            and switch['contributions']['m'][:3] == [106958, 106959, 106960]
            and switch['memory']['m'] == 106940
            and switch['modified_fields']['m'][:3] == [107095, 107096, 107097]
        )
        mod_binding_shape = (
            mchart['physical']['U'][:2] == [140685, 140686]
            and schart['physical']['U'][:2] == [148085, 148086]
            and mchart['physical']['m'][:3] == [140697, 140698, 140699]
            and schart['physical']['m'][:3] == [148097, 148098, 148099]
            and mchart['physical']['Q'][:2] == [140780, 140784]
            and schart['physical']['Q'][:2] == [148180, 148184]
        )
    except (KeyError, TypeError, IndexError):
        loop_binding_shape = mod_binding_shape = False
    guard('typed_loop_endpoint_bindings', loop_binding_shape)
    guard('typed_physical_Q_endpoint_bindings', mod_binding_shape)

    loop_endpoint_rows = {
        104467: {'operation': 'sum', 'arguments': [72149, 104342]},
        104468: {'operation': 'sum', 'arguments': [72150, 104343]},
        104469: {'operation': 'sum', 'arguments': [72151, 104344]},
        104503: {'operation': 'function_substitution', 'expression': 104342,
                 'variable': 70265, 'value': 1, 'Z_independent_substitution': True},
        104504: {'operation': 'function_substitution', 'expression': 104343,
                 'variable': 70265, 'value': 1, 'Z_independent_substitution': True},
        104505: {'operation': 'function_substitution', 'expression': 104344,
                 'variable': 70265, 'value': 1, 'Z_independent_substitution': True},
        106939: {'operation': 'product', 'arguments': [220, 85903]},
        106940: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 106939},
        85903: {'operation': 'product', 'arguments': [9, 73190]},
        106964: {'operation': 'product', 'arguments': [104503, 106940]},
        106965: {'operation': 'product', 'arguments': [104504, 106940]},
        106966: {'operation': 'product', 'arguments': [104505, 106940]},
        106970: {'operation': 'sum', 'arguments': [106964, 106958]},
        106971: {'operation': 'sum', 'arguments': [106965, 106959]},
        106972: {'operation': 'sum', 'arguments': [106966, 106960]},
        107095: {'operation': 'sum', 'arguments': [87591, 106970]},
        107096: {'operation': 'sum', 'arguments': [87592, 106971]},
        107097: {'operation': 'sum', 'arguments': [87593, 106972]},
    }
    loop_rows_ok = all(row_is(rows, number, **expected)
                       for rows in (loop_rows, mod_rows)
                       for number, expected in loop_endpoint_rows.items())
    guard('loop_incoming_alias_decay_and_m_source_rows', loop_rows_ok)

    integral_rows_ok = all(
        row_is(loop_rows, number, operation='definite_integral',
               variable='true_loop_C5_transport_first_switch', lower=0, upper=73190,
               exact_function_integral=True, numerical_value_not_installed=True,
               endpoint_independent_of_Z=True, ordinary_slow_Z_derivative_order=order,
               true_core_inlet_complete_history=True, numerical_cap_not_defining_integrand=True,
               actual_finite_N_signed_correction=True, original_rate='1',
               original_chart='first_switch', physical_Jacobian_applied_once=True,
               zero_correction_inlet_only_at_original_rminus=True,
               predecessor_memory_never_reset=True)
        for number, order in ((106958, 0), (106959, 1), (106960, 2))
    )
    guard('first_switch_equal_bound_integrals_through_Z2', integral_rows_ok)

    # The accepted certificate names the original source-chart m rows, but
    # the row operations are authenticated again here from both decoded
    # archives.  This prevents a literal registry entry from becoming its
    # own source proof.
    original_m_source_rows_ok = all(
        row_is(rows, 72149, operation='product', arguments=[71950, 69245])
        and row_is(rows, 72150, operation='product', arguments=[71951, 69245])
        and row_is(rows, 72151, operation='product', arguments=[71952, 69245])
        and row_is(rows, 87591, operation='product', arguments=[87392, 69245])
        and row_is(rows, 87592, operation='product', arguments=[87393, 69245])
        and row_is(rows, 87593, operation='product', arguments=[87394, 69245])
        and row_is(rows, 69245, operation='positive_quotient', numerator=1,
                   denominator=33572, source_positive_certificate='same strictly positive Pstar')
        and row_is(rows, 33572, operation='analytic_unary', name='exp', argument=5)
        for rows in (loop_rows, mod_rows)
    )
    original_m_source_rows_ok = original_m_source_rows_ok and all(
        row_is(rows, 71950, operation='sum', arguments=[71936, 71943])
        and row_is(rows, 71951, operation='sum', arguments=[71937, 71944])
        and row_is(rows, 71952, operation='sum', arguments=[71938, 71945])
        and row_is(rows, 87392, operation='sum', arguments=[87378, 87385])
        and row_is(rows, 87393, operation='sum', arguments=[87379, 87386])
        and row_is(rows, 87394, operation='sum', arguments=[87380, 87387])
        for rows in (loop_rows, mod_rows)
    )
    guard('original_m_source_chart_rows_authenticated', original_m_source_rows_ok)

    pstar_rows_ok = (
        row_is(mod_rows, 140697, operation='product', arguments=[104467, 33572])
        and row_is(mod_rows, 140698, operation='product', arguments=[104468, 33572])
        and row_is(mod_rows, 140699, operation='product', arguments=[104469, 33572])
        and row_is(mod_rows, 148097, operation='product', arguments=[107095, 33572])
        and row_is(mod_rows, 148098, operation='product', arguments=[107096, 33572])
        and row_is(mod_rows, 148099, operation='product', arguments=[107097, 33572])
        and row_is(mod_rows, 33572, operation='analytic_unary', name='exp', argument=5)
        and row_is(mod_rows, 69245, operation='positive_quotient', numerator=1,
                   denominator=33572, source_positive_certificate='same strictly positive Pstar')
    )
    guard('physical_m_Pstar_applied_once', pstar_rows_ok)

    # Only the order 0 and order 1 Q numerator rows are needed for ordinary
    # Z C1.  Each pair is checked independently; the distinct node numbers
    # remain distinct source records.
    q_row_pairs = [
        (140685, 148085, {'operation': 'product', 'arguments': [104304, 33572]},
                      {'operation': 'product', 'arguments': [106932, 33572]}),
        (140686, 148086, {'operation': 'product', 'arguments': [104305, 33572]},
                      {'operation': 'product', 'arguments': [106933, 33572]}),
        (140715, 148115, {'operation': 'product', 'arguments': [20, 140685]},
                      {'operation': 'product', 'arguments': [20, 148085]}),
        (140716, 148116, {'operation': 'product', 'arguments': [20, 140686]},
                      {'operation': 'product', 'arguments': [20, 148086]}),
        (140717, 148117, {'operation': 'sum', 'arguments': [140716, 140685]},
                      {'operation': 'sum', 'arguments': [148116, 148085]}),
        (140727, 148127, {'operation': 'product', 'arguments': [140715, 48]},
                      {'operation': 'product', 'arguments': [148115, 48]}),
        (140728, 148128, {'operation': 'product', 'arguments': [140717, 48]},
                      {'operation': 'product', 'arguments': [148117, 48]}),
        (140732, 148132, {'operation': 'product', 'arguments': [20, 140697]},
                      {'operation': 'product', 'arguments': [20, 148097]}),
        (140733, 148133, {'operation': 'product', 'arguments': [20, 140698]},
                      {'operation': 'product', 'arguments': [20, 148098]}),
        (140734, 148134, {'operation': 'sum', 'arguments': [140733, 140697]},
                      {'operation': 'sum', 'arguments': [148133, 148097]}),
        (140735, 148135, {'operation': 'product', 'arguments': [20, 140699]},
                      {'operation': 'product', 'arguments': [20, 148099]}),
        (140736, 148136, {'operation': 'product', 'arguments': [48, 140698]},
                      {'operation': 'product', 'arguments': [48, 148098]}),
        (140737, 148137, {'operation': 'sum', 'arguments': [140735, 140736]},
                      {'operation': 'sum', 'arguments': [148135, 148136]}),
        (140744, 148144, {'operation': 'product', 'arguments': [25222, 140732]},
                      {'operation': 'product', 'arguments': [25222, 148132]}),
        (140745, 148145, {'operation': 'product', 'arguments': [25222, 140734]},
                      {'operation': 'product', 'arguments': [25222, 148134]}),
        (140749, 148149, {'operation': 'product', 'arguments': [25251, 140698]},
                      {'operation': 'product', 'arguments': [25251, 148098]}),
        (140750, 148150, {'operation': 'product', 'arguments': [25251, 140699]},
                      {'operation': 'product', 'arguments': [25251, 148099]}),
        (140751, 148151, {'operation': 'product', 'arguments': [25247, 140698]},
                      {'operation': 'product', 'arguments': [25247, 148098]}),
        (140752, 148152, {'operation': 'sum', 'arguments': [140750, 140751]},
                      {'operation': 'sum', 'arguments': [148150, 148151]}),
        (140765, 148165, {'operation': 'sum', 'arguments': [140744, 140749]},
                      {'operation': 'sum', 'arguments': [148144, 148149]}),
        (140766, 148166, {'operation': 'sum', 'arguments': [140745, 140752]},
                      {'operation': 'sum', 'arguments': [148145, 148152]}),
        (140770, 148170, {'operation': 'negative', 'argument': 140765},
                      {'operation': 'negative', 'argument': 148165}),
        (140771, 148171, {'operation': 'negative', 'argument': 140766},
                      {'operation': 'negative', 'argument': 148166}),
        (140775, 148175, {'operation': 'sum', 'arguments': [140727, 140770]},
                      {'operation': 'sum', 'arguments': [148127, 148170]}),
        (140776, 148176, {'operation': 'sum', 'arguments': [140728, 140771]},
                      {'operation': 'sum', 'arguments': [148128, 148171]}),
        (140781, 148181, {'operation': 'product', 'arguments': [25275, 140780]},
                      {'operation': 'product', 'arguments': [25275, 148180]}),
        (140782, 148182, {'operation': 'negative', 'argument': 140781},
                      {'operation': 'negative', 'argument': 148181}),
        (140783, 148183, {'operation': 'sum', 'arguments': [140776, 140782]},
                      {'operation': 'sum', 'arguments': [148176, 148182]}),
    ]
    q_pairs_ok = all(row_is(mod_rows, left, **left_expected)
                     and row_is(mod_rows, right, **right_expected)
                     for left, right, left_expected, right_expected in q_row_pairs)
    guard('physical_Q_numerator_source_pair_rows', q_pairs_ok)

    factor_rows = {
        20: {'operation': 'bound_variable', 'name': 'Z'},
        48: {'operation': 'exact_rational', 'numerator': 2, 'denominator': 1},
        25220: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 25219},
        25221: {'operation': 'negative', 'argument': 25220},
        25222: {'operation': 'sum', 'arguments': [1, 25221]},
        25246: {'operation': 'sum', 'arguments': [20, 20]},
        25247: {'operation': 'negative', 'argument': 25246},
        25249: {'operation': 'product', 'arguments': [20, 20]},
        25250: {'operation': 'negative', 'argument': 25249},
        25251: {'operation': 'sum', 'arguments': [1, 25250]},
        25270: {'operation': 'product', 'arguments': [25220, 25249]},
        25271: {'operation': 'negative', 'argument': 25270},
        25272: {'operation': 'sum', 'arguments': [1, 25271]},
        25274: {'operation': 'product', 'arguments': [25220, 25246]},
        25275: {'operation': 'negative', 'argument': 25274},
    }
    guard('common_Z_delta_d_L_and_derivative_source_rows',
          all(row_is(mod_rows, number, **expected)
              for number, expected in factor_rows.items()))

    q_roots_ok = (
        row_is(mod_rows, 140780, operation='positive_quotient', numerator=140775,
               denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
        and row_is(mod_rows, 140784, operation='positive_quotient', numerator=140783,
                   denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
        and row_is(mod_rows, 148180, operation='positive_quotient', numerator=148175,
                   denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
        and row_is(mod_rows, 148184, operation='positive_quotient', numerator=148183,
                   denominator=25272, source_positive_certificate='same original L=1-delta*Z^2')
        and row_is(mod_rows, 25272, operation='sum', arguments=[1, 25271])
    )
    guard('common_positive_L_source_quotient', q_roots_ok)

    finite_N_m_Z2_proved = all((
        checks.get('accepted_parent_R_U_and_open_Q_scope', False),
        checks.get('loop_transport_and_predecessor_ast', False),
        checks.get('typed_loop_endpoint_bindings', False),
        checks.get('loop_incoming_alias_decay_and_m_source_rows', False),
        checks.get('first_switch_equal_bound_integrals_through_Z2', False),
    ))
    physical_m_Z2_proved = all((
        finite_N_m_Z2_proved,
        checks.get('accepted_original_m_Z2_certificate', False),
        checks.get('original_m_source_chart_rows_authenticated', False),
        checks.get('physical_m_Pstar_applied_once', False),
    ))
    Q_C1_proved = all((
        physical_m_Z2_proved,
        checks.get('physical_Q_creator_ast', False),
        checks.get('typed_physical_Q_endpoint_bindings', False),
        checks.get('physical_Q_numerator_source_pair_rows', False),
        checks.get('common_positive_L_source_quotient', False),
    ))

    return {
        'schema': Q_FRAGMENT_SCHEMA,
        'status': 'source_Q_C1_leaf_ready' if Q_C1_proved else 'source_guard_failure',
        'source_only': True,
        'numeric_values_computed': False,
        'finite_N_m_Z2_proved': finite_N_m_Z2_proved,
        'physical_m_Z2_proved': physical_m_Z2_proved,
        'Q_C1_proved': Q_C1_proved,
        'J2_Q_identity_proved': Q_C1_proved,
        'ordinary_Z_orders': [0, 1],
        'm_ordinary_Z_orders': [0, 1, 2],
        'Q_C1_conditional_on_original_m_Z2_and_accepted_R_U': True,
        'endpoint_projections': {
            'frozen_macro_upper': {'coordinate_node': 70265, 'coordinate_value': 1,
                                  'radius_node': 71448},
            'first_switch_lower': {'coordinate_node': 73190, 'coordinate_value': 0,
                                   'radius_node': 85905},
        },
        'source_pairs': {
            'physical_Q': [[140780, 148180], [140784, 148184]],
            'physical_Q_numerators': [[140775, 148175], [140783, 148183]],
            'physical_U': [[140685, 148085], [140686, 148086]],
            'physical_m_through_Z2': [[140697, 148097], [140698, 148098], [140699, 148099]],
            'loop_modified_m_through_Z2': [[104467, 107095], [104468, 107096], [104469, 107097]],
        },
        'common_source_factors': {
            'Pstar': 33572,
            'inverse_Pstar': 69245,
            'L': 25272,
            'L_recipe': {'operation': 'sum', 'arguments': [1, 25271],
                         'source_positive_certificate': 'same original L=1-delta*Z^2'},
        },
        'derivation': [
            'first_switch lower coordinate 73190=0 uses incoming m rows 104503..104505 from frozen-macro endpoint substitutions 104503..104505',
            'transport rows 106958..106960 are exact definite integrals with lower 0 and upper 73190 for ordinary Z orders 0,1,2; equal-bound endpoint projection makes each correction zero',
            'decay row 106940 is the source exp(-rate*Jacobian*(coordinate-lower)); its first-switch lower endpoint is one by the exact transport AST and lower projection',
            'modified m rows 107095..107097 are baseline first-switch normalized m rows 87591..87593 plus the transported history',
            'physical m rows 148097..148099 apply the common Pstar row 33572 once, matching the macro physical m rows 140697..140699',
            'the actual physical Q creator uses the same Pstar-scaled V and m, signed numerator, and common positive L quotient; rows 140780/140784 and 148180/148184 consume the paired numerator rows through ordinary Z order 1',
        ],
        'source_ast_hashes': {
            'loop_creator_sha256': LOOP_CREATOR_SHA256,
            'modified_Q_creator_sha256': MOD_CREATOR_SHA256,
        },
        'proof_guards': checks,
        'discharged_caller_premises': [
            'the caller must pass the authenticated original m Z2 certificate with its typed source rows and common Pstar handles',
            'the caller must pass the accepted source-only R and U endpoint identities; this leaf does not infer them from mean E/V',
            'the first-switch lower endpoint projection is a source endpoint rule, not a numeric evaluation',
        ],
        'nonclaims': {
            'Q_y_C1_identity_proved': False,
            'J2_B_identity_proved': False,
            'all_fourteen_B_identities_open': True,
            'numeric_containment_claimed': False,
            'new_loop_domain_admission_claimed': False,
        },
        'source_producer_function_calls': 0,
        'ancestor_replays': 0,
        'runtime_model_metadata': None,
    }
