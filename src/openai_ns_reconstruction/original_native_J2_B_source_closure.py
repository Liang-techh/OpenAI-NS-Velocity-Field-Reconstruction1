"""Native weighted J2 B source identity from accepted R/U/Q C1 joins.

Only this new source leaf is executed. Numerical packet errors remain intact.
The source identity does not admit a physical layer, cone or numerical domain.
"""
import ast
import hashlib
import json

SCHEMA = 'original_native_J2_B_source_closure_v1'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_Q_endpoint_transport.json'
RADIUS = 'lei_ren_part1_paper_compliant_current_original_native_J2_radius_identity.json'
WEIGHTED = 'lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json'
CREATOR = 'src/openai_ns_reconstruction/original_native_active_S1_weighted_ibp.py'
MOD = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
PRIOR_B = 'lei_ren_part1_paper_compliant_current_original_native_micro_macro_B_jump_source_v2.json'
RECORD_NAMES = [PARENT, RADIUS, WEIGHTED, CREATOR, MOD, LOOP, PRIOR_B]
RECORD_HASHES = {
    PARENT: '9fa26ecc482b422a4be1c323c50efe9339efdeb2a438a21f45cc3117ef6a2450',
    RADIUS: 'ccb7afca508bc5b31706356fd9058290c587a527ec534ec701824aaba8da431e',
    WEIGHTED: 'aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28',
    CREATOR: 'cdf33563c34e211a00128e3901f49875e50bdea3b4cd9224cc00f81e3288ed8d',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    PRIOR_B: '2505e298aacbac9e15444e93dd7b054deef201c8d02193938e5391ec7b528048',
}
STAGE_FLAGS = ['accepted_R_U_Q_C1_premises_authenticated', 'actual_J2_weighted_B_ports_bound',
               'J2_B_C1_endpoint_identity_proved', 'J2_native_weighted_jump_zero_proved',
               'native_IN_J2_source_reduction_proved',
               'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['all_native_B_jumps_zero_proved', 'J3_through_J15_B_identities_proved',
               'quantitative_packet_error_reduction_proved', 'numeric_values_computed',
               'J2_Qy_identity_proved', 'source_C5_installation_claimed',
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
    parent, radius, weighted, mod, loop = [decoded[name] for name in
        (PARENT, RADIUS, WEIGHTED, MOD, LOOP)]
    old = parent['component_result']
    check('accepted_parent_scope', parent['source_only'] is True and old['passed'] is True
          and parent['full_goal_complete'] is False and len(parent['source_pins']) == 1847
          and old['J2_Q_identity_proved'] is True and old['J2_B_identity_proved'] is False
          and old['ordinary_Z_orders'] == [0, 1] and old['owned_Z_cell'] == ['0', '1/2']
          and old['source_function_identity_only'] is True and old['numeric_values_computed'] is False
          and parent['inherited_R_identity_jumps'] == [2, 3, 4, 6, 8, 9, 10]
          and parent['B_identity_scope']['unresolved_B_identity_nodes'] == list(range(522, 536)))
    check('inherited_input_pins', all(parent['source_pins'][name] == RECORD_HASHES[name]
          for name in (RADIUS, WEIGHTED, CREATOR, MOD, LOOP, PRIOR_B)))
    check('same_source_family', old['source_family'] == radius['component_result']['source_family']
          == weighted['source_family'] == mod['source_family'] == loop['source_family'] == FAMILY)
    prior = decoded[PRIOR_B]['component_result']
    check('accepted_inherited_J0_J1_source_certificates', prior['passed'] is True
          and prior['source_only'] is True and prior['source_family'] == FAMILY
          and prior['ordinary_Z_orders'] == [0, 1]
          and prior['actual_native_two_B_jump_source_identities_bound'] is True
          and prior['source_interface']['proved_weighted_jump_nodes'] == [520, 521]
          and prior['source_interface']['old_archived_jump_rows_unchanged'] is True
          and prior['source_interface']['old_packet_errors_deleted'] is False)
    check('inherited_native_IN_identity', decoded[PRIOR_B]['IN_identity'] ==
          'I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)')
    # The fresh bounded source-composition fragment is integrated before freeze.
    proof_result = prove_J2_B_source(parent, radius, weighted, decoded[CREATOR], mod, loop,
        lambda name, value: check('J2_B:' + name, value))
    check('fresh_B_and_signed_jump', proof_result['J2_B_C1_proved'] is True
          and proof_result['J2_jump_zero_proved'] is True)
    retained = parent['retained_construction_data']
    check('retained_ledger', all(retained[key] == value for key, value in
          dict(cumulative_terms=212, raw_signed_terms=213, reduced_terms=183,
               six_signed_channels_retained=True, all_old_error_data_retained=True,
               reserved_contribution_registered_once=True, inherited_errors_deleted=False).items()))
    graph = [dict(id='R_U_Q_premises', rule='authenticate accepted exact endpoint C1 certificates'),
             dict(id='source_binding', rule='actual weighted B creator, endpoint ports and common L',
                  depends_on=['R_U_Q_premises']),
             dict(id='common_C1_B', rule='same rational C1 composition on the common declared domain',
                  depends_on=['source_binding']),
             dict(id='signed_jump', rule='right B minus identical left B, ordinary Z orders 0 and 1',
                  depends_on=['common_C1_B']),
             dict(id='native_IN_reduction', rule='substitute proved zero J2 into the accepted native source identity',
                  depends_on=['signed_jump'])]
    return dict(schema=SCHEMA, component='B_source_closure', passed=True, source_only=True,
        **dict.fromkeys(STAGE_FLAGS, True), **dict.fromkeys(FALSE_FLAGS, False),
        stage_flags=dict.fromkeys(STAGE_FLAGS, True), nonclaims=dict.fromkeys(FALSE_FLAGS, False),
        J2_B_identity_proved=True, J2_Q_identity_proved=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, owned_Z_cell=['0', '1/2'],
        input_digest=expected_input_digest, source_family=FAMILY,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_function_identity_only=True, source_checks=checks,
        proof_graph=graph, proof_graph_sha256=digest(graph), component_proofs=dict(J2_B=proof_result),
        source_binding_scope=parent['source_binding_scope'],
        retained_construction_data=retained, quantity_distinctions=parent['quantity_distinctions'],
        inherited_R_identity_jumps=[2, 3, 4, 6, 8, 9, 10],
        B_identity_scope=dict(weighted_jump_node=522, jump_index=2, left_B_node=356, right_B_node=366,
            jump_arguments=[366, 356], ordinary_Z_orders=[0, 1],
            endpoint_substitutions={'70265': 1, '73190': 0}, common_L_node=25272,
            newly_proved_B_identity_nodes=[522], inherited_B_identity_nodes=[520, 521],
            cumulative_B_identity_nodes=[520, 521, 522],
            unresolved_B_identity_nodes=list(range(523, 536)), remaining_B_identity_count=13,
            all_fourteen_B_identities_open=False, all_native_B_jumps_zero_proved=False),
        unresolved_B_identity_nodes=list(range(523, 536)),
        remaining_B_identity_count=13, all_fourteen_B_identities_open=False,
        inherited_Q_identity_scope=parent['Q_identity_scope'],
        inherited_U_identity_scope=parent['inherited_U_identity_scope'],
        native_weighted_B_definition=dict(quantity_tag='native_weighted_IBP_boundary_primitive',
            formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)',
            not_loop_B_over_Pstar=True, not_cone_shear_2Vy=True),
        source_jump_zero_does_not_delete_quantitative_errors=True,
        native_IN_inherited_identity='I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)',
        native_IN_reduced_source_identity='I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)',
        native_IN_reduction_is_source_only=True, archived_all_16_J_terms_retained=True,
        unresolved_gates=[], remaining_construction_gates=[
            'the other thirteen native B jumps and same-source quantitative packet transport',
            'independent terminal functional moments and retained quantitative errors',
            'heat exterior, stress/remainder, temporal n recursion, realizable corrections and Cartesian NS validation'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_proof_functions_called=0, generic_DAG_walks=0, numerical_evaluations=0,
            quadratures=0, recurrence_calls=0, runtime_model_metadata=None),
        runtime_model_metadata=None)


"""Unexecuted source leaf for the native weighted J2 B seam.

This fragment consumes decoded records that the parent closure has already
authenticated. It only parses the supplied active-S1 creator text and checks
fixed graph rows and typed endpoint bindings. It never calls a source
producer, ancestor proof, graph walker, evaluator, operator, constructor,
quadrature, recurrence, or numerical routine. The function is intentionally
left uncalled for the root freeze/entry step.
"""

import ast as _J2B_ast
import hashlib as _J2B_hashlib


J2B_SCHEMA = "native_J2_B_source_fragment_v1"
J2B_PARENT_SHA256 = "9fa26ecc482b422a4be1c323c50efe9339efdeb2a438a21f45cc3117ef6a2450"
J2B_RADIUS_SHA256 = "ccb7afca508bc5b31706356fd9058290c587a527ec534ec701824aaba8da431e"
J2B_WEIGHTED_SHA256 = "aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28"
J2B_CREATOR_SHA256 = "cdf33563c34e211a00128e3901f49875e50bdea3b4cd9224cc00f81e3288ed8d"
J2B_MOD_RAW_GZIP_SHA256 = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
J2B_LOOP_RAW_GZIP_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"
J2B_B_CREATOR_LOWER_CALL_SHA256 = "f0fb3ed9db5a246f7a89929ed7ac1c98c845e5cca8dee2670f17a396a2dd153e"
J2B_B_CREATOR_UPPER_CALL_SHA256 = "6f5637e8ba84bd88694fc71d33e01282fab264ff3535923997bf253d55416760"
J2B_B_CREATOR_JUMP_CALL_SHA256 = "ff37a9d2c214a531d9329a67a79a407a150d159deccec573d037511dd9034a5f"
J2B_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256":
        "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256":
        "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256":
        "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

J2B_RADIUS_PIN = "lei_ren_part1_paper_compliant_current_original_native_J2_radius_identity.json"
J2B_WEIGHTED_PIN = "lei_ren_part1_paper_compliant_current_original_native_active_S1_weighted_ibp.json"
J2B_CREATOR_PIN = "src/openai_ns_reconstruction/original_native_active_S1_weighted_ibp.py"
J2B_MOD_PIN = "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
J2B_LOOP_PIN = "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"


def prove_J2_B_source(parent, radius_report, weighted_ibp_report,
                      creator_text, mod, loop, check):
    """Return an unexecuted source-only C1 certificate for weighted J2 B.

    The parent supplies accepted R/U/Q C1 certificates. This leaf binds those
    premises to the actual weighted B ports, the active-S1 B creator AST, the
    two decoded physical archives, and the retained jump row. It performs no
    mathematical evaluation; equalities are source-composition implications
    on the already declared common domain.
    """
    checks = {}

    def guard(name, condition):
        value = bool(condition)
        check(name, value)
        checks[name] = value
        return value

    def node(rows, number):
        if not isinstance(rows, list) or not isinstance(number, int):
            return None
        if number < 0 or number >= len(rows):
            return None
        value = rows[number]
        return value if isinstance(value, dict) else None

    def row_is(rows, number, **expected):
        actual = node(rows, number)
        return actual is not None and all(
            actual.get(key) == value for key, value in expected.items())

    parent_component = (
        parent.get("component_result", {}) if isinstance(parent, dict) else {})
    radius_component = (
        radius_report.get("component_result", {})
        if isinstance(radius_report, dict) else {})
    weighted_component = (
        weighted_ibp_report.get("component_result", {})
        if isinstance(weighted_ibp_report, dict) else {})
    weighted_graph = (
        weighted_component.get("definition_graph", {})
        if isinstance(weighted_component, dict) else {})
    weighted_rows = weighted_graph.get("nodes", [])
    mod_rows = mod.get("exact_graph_nodes", []) if isinstance(mod, dict) else []
    loop_rows = loop.get("exact_graph_nodes", []) if isinstance(loop, dict) else []
    mod_bindings = (
        mod.get("modified_Q4_pressure_C3_bindings", {})
        if isinstance(mod, dict) else {})
    loop_bindings = (
        loop.get("prescribed_loop_C5_bindings", {})
        if isinstance(loop, dict) else {})

    guard("weighted_graph_length", len(weighted_rows) == 545 and weighted_graph.get("node_count") == 545)
    guard("mod_graph_length", len(mod_rows) == 183700)
    guard("loop_graph_length", len(loop_rows) == 121834)
    guard("accepted_weighted_C1_definition_scope",
          weighted_ibp_report.get("source_only") is True
          and weighted_component.get("passed") is True
          and weighted_component.get("ordinary_Z_orders") == [0, 1]
          and weighted_component.get("actual_native_active_S1_weighted_IBP_C1_identity_proved") is True
          and weighted_ibp_report.get("definition_graph") == weighted_graph)
    actual_L = weighted_component.get("source_interface", {}).get("actual_L", {})
    guard("weighted_actual_L_binding", actual_L.get("source_record") == J2B_MOD_PIN
          and actual_L.get("positive_quotient_denominator") == 25272
          and actual_L.get("formula") == "1-delta*Z^2"
          and actual_L.get("Z_independent_delta_port") is True)

    pins = parent.get("source_pins", {}) if isinstance(parent, dict) else {}
    pins_ok = (
        isinstance(pins, dict)
        and len(pins) == 1847
        and pins.get(J2B_RADIUS_PIN) == J2B_RADIUS_SHA256
        and pins.get(J2B_WEIGHTED_PIN) == J2B_WEIGHTED_SHA256
        and pins.get(J2B_CREATOR_PIN) == J2B_CREATOR_SHA256
        and pins.get(J2B_MOD_PIN) == J2B_MOD_RAW_GZIP_SHA256
        and pins.get(J2B_LOOP_PIN) == J2B_LOOP_RAW_GZIP_SHA256
    )
    guard("authenticated_parent_source_pins", pins_ok)

    accepted_parent = (
        isinstance(parent, dict)
        and parent.get("source_only") is True
        and parent.get("full_goal_complete") is False
        and parent_component.get("passed") is True
        and parent_component.get("source_only") is True
        and parent_component.get("numeric_values_computed") is False
        and parent_component.get("J2_Q_identity_proved") is True
        and parent_component.get("J2_B_identity_proved") is False
        and parent_component.get("ordinary_Z_orders") == [0, 1]
        and parent_component.get("owned_Z_cell") == ["0", "1/2"]
        and parent_component.get("source_function_identity_only") is True
        and parent.get("inherited_R_identity_jumps") == [2, 3, 4, 6, 8, 9, 10]
        and parent.get("B_identity_scope", {}).get("unresolved_B_identity_nodes")
            == list(range(522, 536))
        and parent_component.get("all_native_B_jumps_zero_proved") is False
    )
    guard("accepted_parent_R_U_Q_scope", accepted_parent)

    family_ok = (
        isinstance(parent_component, dict)
        and parent_component.get("source_family") == J2B_SOURCE_FAMILY
        and radius_component.get("source_family") == J2B_SOURCE_FAMILY
        and weighted_ibp_report.get("source_family") == J2B_SOURCE_FAMILY
        and isinstance(mod, dict) and mod.get("source_family") == J2B_SOURCE_FAMILY
        and isinstance(loop, dict) and loop.get("source_family") == J2B_SOURCE_FAMILY
    )
    guard("shared_actual_source_family", family_ok)

    radius_identity = radius_component.get("radius_identity", {})
    radius_ok = (
        radius_component.get("passed") is True
        and radius_component.get("source_only") is True
        and radius_component.get("numeric_values_computed") is False
        and radius_component.get("J2_radius_C1_identity_proved") is True
        and radius_component.get("J2_B_identity_proved") is False
        and radius_component.get("ordinary_Z_orders") == [0, 1]
        and radius_component.get("owned_Z_cell") == ["0", "1/2"]
        and radius_identity.get("field") == "R"
        and radius_identity.get("weighted_jump_node") == 522
        and radius_identity.get("jump_index") == 2
        and radius_identity.get("source_function_identity_proved") is True
        and radius_identity.get("source_parameters_aliased") is False
        and radius_identity.get("left_source_node") == 71448
        and radius_identity.get("right_source_node") == 85905
        and radius_identity.get("left_coordinate_node") == 70265
        and radius_identity.get("right_coordinate_node") == 73190
        and radius_identity.get("ordinary_Z_orders") == [0, 1]
        and radius_identity.get("exact_endpoint_rational")
            == {"numerator": 100, "denominator": 1}
        and radius_identity.get("exact_ordinary_Z_derivative")
            == {"numerator": 0, "denominator": 1}
        and radius_identity.get("endpoint_substitutions")
            == {"70265": 1, "73190": 0}
        and radius_identity.get("numeric_value") is None
    )
    guard("accepted_radius_C1_scope", radius_ok)

    u_scope = (
        parent.get("inherited_U_identity_scope", {})
        if isinstance(parent, dict) else {})
    q_scope = (
        parent.get("Q_identity_scope", {})
        if isinstance(parent, dict) else {})
    accepted_u = (
        u_scope.get("weighted_jump_node") == 522
        and u_scope.get("jump_index") == 2
        and u_scope.get("full_U_C1_identity_proved") is True
        and u_scope.get("ordinary_Z_orders") == [0, 1]
        and u_scope.get("left_source_nodes") == [140685, 140686]
        and u_scope.get("right_source_nodes") == [148085, 148086]
        and u_scope.get("endpoint_substitutions")
            == {"70265": 1, "73190": 0}
    )
    accepted_q = (
        q_scope.get("weighted_jump_node") == 522
        and q_scope.get("jump_index") == 2
        and q_scope.get("full_Q_C1_identity_proved") is True
        and q_scope.get("ordinary_Z_orders") == [0, 1]
        and q_scope.get("left_source_nodes") == [140780, 140784]
        and q_scope.get("right_source_nodes") == [148180, 148184]
        and q_scope.get("common_L_denominator_node") == 25272
        and q_scope.get("endpoint_substitutions")
            == {"70265": 1, "73190": 0}
    )
    guard("accepted_U_C1_certificate", accepted_u)
    guard("accepted_Q_C1_certificate", accepted_q)

    creator_ok = False
    endpoint_body_hashes = {}
    jump_body_hash = None
    try:
        creator_tree = _J2B_ast.parse(creator_text)
        creator_ok = (
            _J2B_hashlib.sha256(creator_text.encode()).hexdigest()
            == J2B_CREATOR_SHA256)
        endpoint_calls = []
        jump_calls = []
        for item in _J2B_ast.walk(creator_tree):
            if not isinstance(item, _J2B_ast.Call):
                continue
            if not isinstance(item.func, _J2B_ast.Attribute):
                continue
            if item.func.attr != "node" or not item.args:
                continue
            first = item.args[0]
            if isinstance(first, _J2B_ast.Constant):
                if first.value == "actual_B_endpoint_function":
                    endpoint_calls.append(item)
                elif first.value == "actual_internal_primitive_jump":
                    jump_calls.append(item)

        def kw(call, name):
            for item in call.keywords:
                if item.arg == name:
                    return item.value
            return None

        def const(call, name):
            value = kw(call, name)
            return value.value if isinstance(value, _J2B_ast.Constant) else None

        def name(call, name):
            value = kw(call, name)
            return value.id if isinstance(value, _J2B_ast.Name) else None

        def expr(call, name):
            value = kw(call, name)
            return _J2B_ast.get_source_segment(creator_text, value)

        endpoint_semantics = []
        for item in endpoint_calls:
            endpoint = const(item, "endpoint")
            endpoint_bound = const(item, "endpoint_bound")
            endpoint_semantics.append({
                "endpoint": endpoint,
                "endpoint_bound": endpoint_bound,
                "formula": const(item, "formula"),
                "measure": const(item, "measure"),
                "integration_variable": const(item, "integration_variable"),
                "limits_Z_independent": const(item, "limits_Z_independent"),
                "source_function_identity_not_claimed":
                    const(item, "source_function_identity_not_claimed"),
                "physical_ports": expr(item, "physical_ports"),
                "radius": expr(item, "radius"),
                "Jacobian": expr(item, "Jacobian"),
                "physical_measure": expr(item, "physical_measure"),
                "bounds": expr(item, "bounds"),
                "domain": expr(item, "domain"),
            })
            source = _J2B_ast.get_source_segment(creator_text, item)
            if source is not None:
                endpoint_body_hashes[endpoint] = _J2B_hashlib.sha256(
                    (source + "\n").encode()).hexdigest()

        jump_body = jump_calls[0] if len(jump_calls) == 1 else None
        if jump_body is not None:
            source = _J2B_ast.get_source_segment(creator_text, jump_body)
            if source is not None:
                jump_body_hash = _J2B_hashlib.sha256(
                    (source + "\n").encode()).hexdigest()

        expected_common = {
            "formula": "R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)",
            "measure": "R^2*Jacobian*dnative",
            "integration_variable": "native_coordinate",
            "limits_Z_independent": True,
            "source_function_identity_not_claimed": True,
            "physical_ports": "physical",
            "radius": "radius",
            "Jacobian": "geometry.get('Jacobian')",
            "physical_measure": "measure",
            "bounds": "geometry.get('domain')",
            "domain": "geometry.get('domain')",
        }
        creator_endpoint_semantics = (
            len(endpoint_calls) == 2
            and sorted(item.get("endpoint") for item in endpoint_semantics)
                == ["lower", "upper"]
            and all(
                all(item.get(key) == value
                    for key, value in expected_common.items())
                for item in endpoint_semantics)
            and endpoint_body_hashes.get("lower")
                == J2B_B_CREATOR_LOWER_CALL_SHA256
            and endpoint_body_hashes.get("upper")
                == J2B_B_CREATOR_UPPER_CALL_SHA256
        )
        creator_jump_semantics = (
            len(jump_calls) == 1
            and jump_body_hash == J2B_B_CREATOR_JUMP_CALL_SHA256
            and const(jump_body, "formula") == "J_i=B_next_lower-B_prev_upper"
            and name(jump_body, "jump_index") == "index"
            and const(jump_body, "is_actual_function") is True
            and const(jump_body, "retained") is True
            and const(jump_body, "zero_assumption") is False
            and const(jump_body, "source_identity_proved") is False
        )
        creator_ok = creator_ok and creator_endpoint_semantics and creator_jump_semantics
    except (SyntaxError, TypeError, ValueError):
        creator_ok = False
    guard("authenticated_active_S1_B_creator_AST", creator_ok)

    weighted_ports = (
        row_is(weighted_rows, 59, operation="actual_native_geometry_coordinate",
               chart="frozen_macro", quantity="coordinate", source_nodes=[70265],
               source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 60, operation="actual_native_geometry_domain",
                   chart="frozen_macro", quantity="domain", source_nodes=[0, 1],
                   source_record=J2B_LOOP_PIN, ordinary_Z_orders=[0, 1])
        and row_is(weighted_rows, 61, operation="actual_native_geometry_radius",
                   chart="frozen_macro", quantity="radius", source_nodes=[71448],
                   source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 62, operation="actual_native_geometry_physical_measure",
                   chart="frozen_macro", quantity="physical_measure",
                   source_nodes=[66987], source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 63, operation="actual_native_geometry_Jacobian",
                   chart="frozen_macro", quantity="Jacobian", source_nodes=[66987],
                   source_name="physical_measure", source_path="geometry.physical_measure",
                   source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 64, operation="actual_native_physical_U",
                   chart="frozen_macro", quantity="U",
                   source_nodes=[140685, 140686, 140687, 140688, 140689, 140690],
                   source_record=J2B_MOD_PIN)
        and row_is(weighted_rows, 68, operation="actual_native_physical_Q",
                   chart="frozen_macro", quantity="Q",
                   source_nodes=[140780, 140784, 140790, 140796, 140802],
                   source_record=J2B_MOD_PIN)
        and row_is(weighted_rows, 75, operation="actual_native_geometry_coordinate",
                   chart="first_switch", quantity="coordinate", source_nodes=[73190],
                   source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 76, operation="actual_native_geometry_domain",
                   chart="first_switch", quantity="domain", source_nodes=[0, 1],
                   source_record=J2B_LOOP_PIN, ordinary_Z_orders=[0, 1])
        and row_is(weighted_rows, 77, operation="actual_native_geometry_radius",
                   chart="first_switch", quantity="radius", source_nodes=[85905],
                   source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 78, operation="actual_native_geometry_physical_measure",
                   chart="first_switch", quantity="physical_measure",
                   source_nodes=[9], source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 79, operation="actual_native_geometry_Jacobian",
                   chart="first_switch", quantity="Jacobian", source_nodes=[9],
                   source_name="physical_measure", source_path="geometry.physical_measure",
                   source_record=J2B_LOOP_PIN)
        and row_is(weighted_rows, 80, operation="actual_native_physical_U",
                   chart="first_switch", quantity="U",
                   source_nodes=[148085, 148086, 148087, 148088, 148089, 148090],
                   source_record=J2B_MOD_PIN)
        and row_is(weighted_rows, 84, operation="actual_native_physical_Q",
                   chart="first_switch", quantity="Q",
                   source_nodes=[148180, 148184, 148190, 148196, 148202],
                   source_record=J2B_MOD_PIN)
    )
    guard("authenticated_weighted_B_port_graph", weighted_ports)

    weighted_B = (
        row_is(weighted_rows, 356, operation="actual_B_endpoint_function",
               chart="frozen_macro", endpoint="upper",
               endpoint_bound="domain[1]", coordinate=59, domain=60,
               radius=61, Jacobian=63, physical_measure=62,
               bounds=60, integration_variable="native_coordinate",
               limits_Z_independent=True, measure="R^2*Jacobian*dnative",
               physical_ports={
                   "U": 64, "U_y": 65, "m": 66, "m_y": 67,
                   "Q": 68, "Q_y": 69, "S": 70,
                   "local_integrals.S1": 71,
                   "cumulative_integrals.S1": 72,
                   "incoming_integrals.S1": 73, "R_Q_y": 74},
               formula="R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)",
               source_function_identity_not_claimed=True)
        and row_is(weighted_rows, 366, operation="actual_B_endpoint_function",
                   chart="first_switch", endpoint="lower",
                   endpoint_bound="domain[0]", coordinate=75, domain=76,
                   radius=77, Jacobian=79, physical_measure=78,
                   bounds=76, integration_variable="native_coordinate",
                   limits_Z_independent=True, measure="R^2*Jacobian*dnative",
                   physical_ports={
                       "U": 80, "U_y": 81, "m": 82, "m_y": 83,
                       "Q": 84, "Q_y": 85, "S": 86,
                       "local_integrals.S1": 87,
                       "cumulative_integrals.S1": 88,
                       "incoming_integrals.S1": 89, "R_Q_y": 90},
                   formula="R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)",
                   source_function_identity_not_claimed=True)
    )
    guard("authenticated_weighted_B_endpoint_rows", weighted_B)

    weighted_jump = row_is(
        weighted_rows, 522, operation="actual_internal_primitive_jump",
        arguments=[366, 356], jump_index=2,
        formula="J_i=B_next_lower-B_prev_upper",
        is_actual_function=True, retained=True, zero_assumption=False,
        source_identity_proved=False)
    guard("authenticated_weighted_J2_jump_row", weighted_jump)

    try:
        macro_loop = loop_bindings["charts"]["frozen_macro"]
        switch_loop = loop_bindings["charts"]["first_switch"]
        macro_mod = mod_bindings["charts"]["frozen_macro"]
        switch_mod = mod_bindings["charts"]["first_switch"]
        typed_bindings = (
            macro_loop.get("coordinate") == 70265
            and macro_loop.get("domain") == [0, 1]
            and macro_loop.get("radius") == 71448
            and switch_loop.get("coordinate") == 73190
            and switch_loop.get("domain") == [0, 1]
            and switch_loop.get("radius") == 85905
            and macro_mod.get("physical", {}).get("U")
                == [140685, 140686, 140687, 140688, 140689, 140690]
            and switch_mod.get("physical", {}).get("U")
                == [148085, 148086, 148087, 148088, 148089, 148090]
            and macro_mod.get("physical", {}).get("Q")
                == [140780, 140784, 140790, 140796, 140802]
            and switch_mod.get("physical", {}).get("Q")
                == [148180, 148184, 148190, 148196, 148202]
        )
    except (KeyError, TypeError, IndexError):
        typed_bindings = False
    guard("typed_macro_switch_B_endpoint_bindings", typed_bindings)

    geometry_rows = (
        row_is(loop_rows, 0, operation="exact_rational", numerator=0, denominator=1)
        and row_is(loop_rows, 1, operation="exact_rational", numerator=1, denominator=1)
        and row_is(loop_rows, 70265, operation="bound_variable",
                   name="original_prescribed_C6_macro_fraction")
        and row_is(loop_rows, 73190, operation="bound_variable",
                   name="true_prescribed_first_switch_s")
        and row_is(loop_rows, 71448, operation="product",
                   arguments=[65604, 71447])
        and row_is(loop_rows, 85905, operation="product",
                   arguments=[31709, 85904])
        and row_is(loop_rows, 66987, operation="sum",
                   arguments=[66984, 66986])
        and row_is(loop_rows, 9, operation="current_original_source_parameter",
                   name="hbS", Z_independent=True)
    )
    guard("authenticated_radius_coordinate_source_rows", geometry_rows)

    u_rows = (
        row_is(mod_rows, 140685, operation="product",
               arguments=[104304, 33572])
        and row_is(mod_rows, 140686, operation="product",
                   arguments=[104305, 33572])
        and row_is(mod_rows, 148085, operation="product",
                   arguments=[106932, 33572])
        and row_is(mod_rows, 148086, operation="product",
                   arguments=[106933, 33572])
    )
    q_rows = (
        row_is(mod_rows, 140780, operation="positive_quotient",
               numerator=140775, denominator=25272,
               source_positive_certificate="same original L=1-delta*Z^2")
        and row_is(mod_rows, 140784, operation="positive_quotient",
                   numerator=140783, denominator=25272,
                   source_positive_certificate="same original L=1-delta*Z^2")
        and row_is(mod_rows, 148180, operation="positive_quotient",
                   numerator=148175, denominator=25272,
                   source_positive_certificate="same original L=1-delta*Z^2")
        and row_is(mod_rows, 148184, operation="positive_quotient",
                   numerator=148183, denominator=25272,
                   source_positive_certificate="same original L=1-delta*Z^2")
    )
    common_L = row_is(mod_rows, 25272, operation="sum", arguments=[1, 25271])
    guard("authenticated_physical_U_endpoint_rows", u_rows)
    guard("authenticated_physical_Q_endpoint_rows", q_rows)
    guard("authenticated_common_L_source_row", common_L)

    endpoint_proven = all((
        checks.get("accepted_parent_R_U_Q_scope", False),
        checks.get("accepted_radius_C1_scope", False),
        checks.get("accepted_U_C1_certificate", False),
        checks.get("accepted_Q_C1_certificate", False),
        checks.get("authenticated_active_S1_B_creator_AST", False),
        checks.get("authenticated_weighted_B_port_graph", False),
        checks.get("authenticated_weighted_B_endpoint_rows", False),
        checks.get("typed_macro_switch_B_endpoint_bindings", False),
        checks.get("authenticated_radius_coordinate_source_rows", False),
        checks.get("authenticated_physical_U_endpoint_rows", False),
        checks.get("authenticated_physical_Q_endpoint_rows", False),
        checks.get("authenticated_common_L_source_row", False),
    ))
    jump_proven = all((
        endpoint_proven,
        checks.get("authenticated_weighted_J2_jump_row", False),
    ))

    return {
        "schema": J2B_SCHEMA,
        "status": "source_J2_B_C1_leaf_ready" if jump_proven else "source_guard_failure",
        "source_only": True,
        "numeric_values_computed": False,
        "J2_B_C1_proved": endpoint_proven,
        "J2_jump_zero_proved": jump_proven,
        "J2_B_identity_proved": jump_proven,
        "ordinary_Z_orders": [0, 1],
        "ordinary_Z_not_temporal_n": True,
        "owned_Z_cell": ["0", "1/2"],
        "endpoint_projections": {
            "frozen_macro_upper": {
                "coordinate_node": 70265, "coordinate_value": 1,
                "radius_node": 71448, "weighted_B_node": 356},
            "first_switch_lower": {
                "coordinate_node": 73190, "coordinate_value": 0,
                "radius_node": 85905, "weighted_B_node": 366},
        },
        "source_pairs": {
            "radius": [[71448, 85905]],
            "physical_U": [[140685, 148085], [140686, 148086]],
            "physical_Q": [[140780, 148180], [140784, 148184]],
            "weighted_B": [[356, 366]],
            "weighted_jump": [366, 356],
        },
        "typed_binding_nodes": {
            "left": {
                "coordinate": 70265, "domain": [0, 1], "radius": 71448,
                "weighted_ports": {"coordinate": 59, "domain": 60,
                                  "radius": 61, "Jacobian": 63,
                                  "physical_measure": 62, "U": 64, "Q": 68},
                "B": 356,
            },
            "right": {
                "coordinate": 73190, "domain": [0, 1], "radius": 85905,
                "weighted_ports": {"coordinate": 75, "domain": 76,
                                  "radius": 77, "Jacobian": 79,
                                  "physical_measure": 78, "U": 80, "Q": 84},
                "B": 366,
            },
            "common_L": 25272,
        },
        "C1_derivative_rule": "B_Z=2*R*R_Z*H+R^2*(Q_Z/L-Q*L_Z/L^2+Q*Q_Z-4*Z*U^2/L^2-4*Z^2*U*U_Z/L^2+4*Z^2*U^2*L_Z/L^3), H=Q/L+Q^2/2-2*Z^2*U^2/L^2",
        "archived_source_status_flags_rewritten": False,
        "quantitative_packet_error_reduction_proved": False,
        "derivation": [
            "The accepted radius C1 certificate binds R nodes 71448 and 85905 at coordinates 70265=1 and 73190=0, for ordinary Z orders 0 and 1, with the source endpoint rational and derivative metadata retained symbolically.",
            "The accepted U and Q C1 certificates bind the exact physical port pairs used by weighted B, while the weighted graph independently binds those ports to the two actual B endpoint rows.",
            "The authenticated active-S1 creator supplies one literal R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2) endpoint recipe at each endpoint, with Z-independent limits and the retained physical endpoint port dictionary.",
            "The decoded source rows bind both U pairs to Pstar products, both Q pairs to the same positive L quotient, and L to row 25272; no numeric value is evaluated.",
            "Weighted jump 522 is the actual right-lower minus left-upper call [366,356]. Therefore the source-composed endpoint C1 identity discharges only J2 and leaves nodes 523 through 535 open.",
        ],
        "proof_guards": checks,
        "source_ast_hashes": {
            "active_S1_creator_sha256": J2B_CREATOR_SHA256,
            "lower_endpoint_call_sha256": J2B_B_CREATOR_LOWER_CALL_SHA256,
            "upper_endpoint_call_sha256": J2B_B_CREATOR_UPPER_CALL_SHA256,
            "jump_call_sha256": J2B_B_CREATOR_JUMP_CALL_SHA256,
        },
        "source_archive_raw_sha256": {
            "modified_Q4_pressure_C3_functions.json.gz": J2B_MOD_RAW_GZIP_SHA256,
            "prescribed_loop_C5_functions.json.gz": J2B_LOOP_RAW_GZIP_SHA256,
        },
        "unresolved": {
            "other_B_jump_nodes": list(range(523, 536)),
            "remaining_B_identity_count": 13,
            "all_fourteen_B_identities_open": False,
            "quantitative_packet_error_transport": True,
            "global_five_moment_closure": True,
            "numeric_values_computed": False,
            "Q_y_C1_identity": False,
            "new_loop_domain_admission": False,
            "admissible_cone_admission": False,
        },
        "source_producer_function_calls": 0,
        "ancestor_replays": 0,
        "graph_walks": 0,
        "normalization_calls": 0,
        "numeric_evaluations": 0,
        "runtime_model_metadata": None,
    }
