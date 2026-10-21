"""Source transport across the actual first-switch / second-switch seam.

This fresh consumer keeps original and finite-N histories separate and proves
only source-function endpoint identities. Quantitative packet errors remain.
"""
import ast
import hashlib
import json
from fractions import Fraction

SCHEMA = 'original_native_J3_endpoint_transport_v1'
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
RECORD_NAMES = [PARENT, LOOP, MOD, LOOP_CREATOR, MOD_CREATOR, C6_CREATOR,
                BRIDGE_CREATOR, BACKGROUND, WEIGHTED, WEIGHTED_CREATOR]
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


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def prove_J3_endpoint_source_algebra(loop, background, check):
    """Apply exact 0/1 substitution rules to a fixed, authenticated spine.

    This is a new bounded source theorem. It never walks or normalizes the
    inherited DAG and never evaluates exp/log/fractional_part or any field.
    """
    nodes = loop['exact_graph_nodes']
    charts = loop['prescribed_loop_C5_bindings']['charts']
    c6 = background['prescribed_R100_Rm_C6_bindings']['charts']
    expected = {
        0: dict(operation='exact_rational', numerator=0, denominator=1),
        1: dict(operation='exact_rational', numerator=1, denominator=1),
        31709: dict(operation='exact_rational', numerator=100, denominator=1),
        73190: dict(operation='bound_variable', name='true_prescribed_first_switch_s'),
        75949: dict(operation='bound_variable', name='true_prescribed_second_switch_t'),
        75950: dict(operation='sum', arguments=[1, 75949]),
        85903: dict(operation='product', arguments=[9, 73190]),
        85904: dict(operation='analytic_unary', name='exp', argument=85903),
        85905: dict(operation='product', arguments=[31709, 85904]),
        88598: dict(operation='product', arguments=[9, 75950]),
        88599: dict(operation='analytic_unary', name='exp', argument=88598),
        88600: dict(operation='product', arguments=[31709, 88599]),
        104533: dict(operation='positive_quotient', numerator=85905, denominator=65604,
                     source_positive_certificate='same positive actual R/Ra'),
        104534: dict(operation='analytic_unary', name='log', argument=104533),
        104535: dict(operation='sum', arguments=[104534, 97298]),
        104536: dict(operation='product', arguments=[18, 104535]),
        104537: dict(operation='analytic_unary', name='fractional_part', argument=104536),
        107161: dict(operation='positive_quotient', numerator=88600, denominator=65604,
                     source_positive_certificate='same positive actual R/Ra'),
        107162: dict(operation='analytic_unary', name='log', argument=107161),
        107163: dict(operation='sum', arguments=[107162, 97298]),
        107164: dict(operation='product', arguments=[18, 107163]),
        107165: dict(operation='analytic_unary', name='fractional_part', argument=107164),
        24954: dict(operation='product', arguments=[8, 10, 28]),
        97298: dict(operation='negative', argument=24954),
        28: dict(operation='exact_rational', numerator=1, denominator=2),
        13: dict(operation='exact_rational', numerator=4, denominator=1),
        14: dict(operation='product', arguments=[13, 5]),
        643: dict(operation='exact_rational', numerator=1000, denominator=1),
        35562: dict(operation='sum', arguments=[14, 643]),
        35564: dict(operation='negative', argument=35562),
        35565: dict(operation='analytic_unary', name='exp', argument=35564),
        65604: dict(operation='product', arguments=[13, 35565]),
    }
    check('endpoint_algebra:actual_fixed_radius_phase_spine',
          all(nodes[number] == row for number, row in expected.items()))
    check('endpoint_algebra:typed_endpoint_domains',
          charts['first_switch']['coordinate'] == c6['first_switch']['phase'] == 73190
          and charts['second_switch']['coordinate'] == c6['second_switch']['phase'] == 75949
          and charts['first_switch']['domain'] == c6['first_switch']['domain'] == [0, 1]
          and charts['second_switch']['domain'] == c6['second_switch']['domain'] == [0, 1]
          and charts['first_switch']['radius'] == c6['first_switch']['radius'] == 85905
          and charts['second_switch']['radius'] == c6['second_switch']['radius'] == 88600
          and charts['first_switch']['global_offset'] == 104535
          and charts['second_switch']['global_offset'] == 107163
          and charts['first_switch']['actual_phase'] == 104537
          and charts['second_switch']['actual_phase'] == 107165)
    check('endpoint_algebra:Z_independent_geometry_parameters',
          all(nodes[number].get('Z_independent') is True for number in (8, 9, 10, 18))
          and nodes[2]['operation'] == 'exact_rational'
          and nodes[3] == dict(operation='analytic_unary', name='exp', argument=2)
          and nodes[4] == dict(operation='exact_rational', numerator=11, denominator=1)
          and nodes[5] == dict(operation='sum', arguments=[3, 4]))
    left_bound = charts['first_switch']['domain'][1]
    right_bound = charts['second_switch']['domain'][0]
    left_coefficient = Fraction(nodes[left_bound]['numerator'], nodes[left_bound]['denominator'])
    right_coefficient = Fraction(nodes[1]['numerator'], nodes[1]['denominator']) + Fraction(
        nodes[right_bound]['numerator'], nodes[right_bound]['denominator'])

    def radius_form(coefficient):
        return ['product', ['rational', 100, 1], ['exp', ['product', ['source', 9],
            ['rational', coefficient.numerator, coefficient.denominator]]]]

    left_radius, right_radius = radius_form(left_coefficient), radius_form(right_coefficient)
    check('endpoint_algebra:radius_congruence_after_source_substitution',
          left_bound == 1 and right_bound == 0 and left_coefficient == right_coefficient == 1
          and left_radius == right_radius)

    def phase_form(radius):
        return ['fractional_part', ['product', ['source', 18], ['sum',
            ['log', ['positive_quotient', radius, ['source', 65604]]],
            ['negative', ['source', 24954]]]]]

    left_phase, right_phase = phase_form(left_radius), phase_form(right_radius)
    check('endpoint_algebra:global_phase_congruence_after_source_substitution', left_phase == right_phase)
    return dict(source_only=True, J3_radius_source_identity_proved=True,
        J3_phase_source_identity_proved=True, ordinary_Z_orders=[0, 1],
        endpoint_substitutions={'73190': 1, '75949': 0},
        radius_left_normal_form=left_radius, radius_right_normal_form=right_radius,
        phase_left_normal_form=left_phase, phase_right_normal_form=right_phase,
        radius_and_phase_Z_independent_proved=True, exact_endpoint_rational_arithmetic_used=True,
        hbB_and_hbS_not_aliased=True, physical_or_source_function_numeric_evaluations=0,
        rule='fixed source substitution; rational 1*hbS=(1+0)*hbS; constructor congruence',
        generic_DAG_walks=0, inherited_normalizer_calls=0)


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
    parent, loop, mod, weighted = [decoded[name] for name in (PARENT, LOOP, MOD, WEIGHTED)]
    old = parent['component_result']
    check('accepted_parent_scope', parent['source_only'] is True and old['passed'] is True
          and parent['full_goal_complete'] is False and len(parent['source_pins']) == 1850
          and old['J2_B_identity_proved'] is True and old['ordinary_Z_orders'] == [0, 1]
          and old['B_identity_scope']['cumulative_B_identity_nodes'] == [520, 521, 522]
          and old['unresolved_B_identity_nodes'] == list(range(523, 536))
          and old['native_IN_reduced_source_identity'] ==
          'I_N=B_lastupper-B_firstlower-sum(J_3..J_15)+sum17(integral K_i dy)')
    check('inherited_input_pins', all(parent['source_pins'][name] == RECORD_HASHES[name]
          for name in RECORD_NAMES if name != PARENT))
    check('same_source_family', old['source_family'] == loop['source_family']
          == mod['source_family'] == weighted['source_family'] == decoded[BACKGROUND]['source_family'] == FAMILY)
    check('archive_creator_input_pins',
          loop['input_hashes'][LOOP_CREATOR] == RECORD_HASHES[LOOP_CREATOR]
          and mod['input_hashes'][MOD_CREATOR] == RECORD_HASHES[MOD_CREATOR]
          and decoded[BACKGROUND]['input_hashes'][C6_CREATOR] == RECORD_HASHES[C6_CREATOR]
          and decoded[BACKGROUND]['input_hashes'][BRIDGE_CREATOR] == RECORD_HASHES[BRIDGE_CREATOR])
    check('common_Z_independent_Pstar_source', all(
          archive['exact_graph_nodes'][2]['operation'] == 'exact_rational'
          and archive['exact_graph_nodes'][3] == dict(operation='analytic_unary', name='exp', argument=2)
          and archive['exact_graph_nodes'][4] == dict(operation='exact_rational', numerator=11, denominator=1)
          and archive['exact_graph_nodes'][5] == dict(operation='sum', arguments=[3, 4])
          and archive['exact_graph_nodes'][33572] == dict(operation='analytic_unary', name='exp', argument=5)
          for archive in (decoded[BACKGROUND], loop, mod)))
    algebra = prove_J3_endpoint_source_algebra(loop, decoded[BACKGROUND], check)
    original = prove_original_switch_transport(decoded[BACKGROUND], loop,
        decoded[C6_CREATOR], decoded[BRIDGE_CREATOR],
        lambda name, value: check('J3_original:' + name, value))
    check('original_switch_state_Z2', original['original_state_Z2_proved'] is True
          and original['original_five_moments_Z2_proved'] is True)
    check('fresh_original_switch_certificate_contract', J3MQ__original_state_contract_ok(original))

    def ru_check(predicates):
        require(isinstance(predicates, dict) and bool(predicates), 'Missing fresh R/U predicates')
        for name, value in predicates.items():
            check('J3_R_U:' + name, type(value) is bool and value is True)
        return True

    ru = prove_J3_R_U_transport(loop, mod, decoded[LOOP_CREATOR], decoded[C6_CREATOR],
        original, ru_check, background=decoded[BACKGROUND], mod_creator_text=decoded[MOD_CREATOR],
        bridge_creator_text=decoded[BRIDGE_CREATOR],
        archive_raw_sha256={'loop': records[LOOP]['raw_sha256'], 'mod': records[MOD]['raw_sha256'],
                            'background': records[BACKGROUND]['raw_sha256']})
    check('new_R_U_C1', ru['J3_R_C1_proved'] is True and ru['J3_U_C1_proved'] is True)
    mq = prove_J3_physical_Q_transport(decoded[BACKGROUND], loop, mod, decoded[LOOP_CREATOR],
        decoded[MOD_CREATOR], original, ru, lambda name, value: check('J3_m_Q:' + name, value))
    check('new_physical_m_Z2_and_Q_C1', mq['physical_m_Z2_proved'] is True
          and mq['finite_N_m_Z2_proved'] is True and mq['J3_Q_C1_proved'] is True)
    rows = weighted['definition_graph']['nodes']
    check('weighted_source_graph_scope', len(rows) == 545
          and rows == weighted['component_result']['definition_graph']['nodes']
          and weighted['component_result']['actual_native_active_S1_weighted_IBP_C1_identity_proved'] is True)

    def row_is(number, **expected):
        return all(rows[number].get(key) == value for key, value in expected.items())

    check('actual_J3_weighted_jump', row_is(523, operation='actual_internal_primitive_jump',
          arguments=[377, 367], jump_index=3, formula='J_i=B_next_lower-B_prev_upper',
          is_actual_function=True, retained=True, zero_assumption=False, source_identity_proved=False))
    check('actual_J3_B_endpoint_recipe',
          row_is(367, operation='actual_B_endpoint_function', chart='first_switch', endpoint='upper',
                 coordinate=75, domain=76, radius=77, endpoint_bound='domain[1]',
                 formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
          and row_is(377, operation='actual_B_endpoint_function', chart='second_switch', endpoint='lower',
                     coordinate=91, domain=92, radius=93, endpoint_bound='domain[0]',
                     formula='R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)', limits_Z_independent=True)
          and rows[367]['physical_ports']['U'] == 80 and rows[367]['physical_ports']['Q'] == 84
          and rows[377]['physical_ports']['U'] == 96 and rows[377]['physical_ports']['Q'] == 100)
    check('actual_J3_source_ports', all(row_is(i, source_nodes=values) for i, values in
          [(75, [73190]), (76, [0, 1]), (77, [85905]), (80, [148085, 148086, 148087, 148088, 148089, 148090]),
           (84, [148180, 148184, 148190, 148196, 148202]), (91, [75949]), (92, [0, 1]),
           (93, [88600]), (96, [154247, 154248, 154249, 154250, 154251, 154252]),
           (100, [154342, 154346, 154352, 154358, 154364])]))
    check('actual_J3_port_archive_types', all(
          row_is(i, operation=operation, source_record=source, chart=chart, quantity=quantity)
          for i, operation, source, chart, quantity in [
              (75, 'actual_native_geometry_coordinate', LOOP, 'first_switch', 'coordinate'),
              (77, 'actual_native_geometry_radius', LOOP, 'first_switch', 'radius'),
              (91, 'actual_native_geometry_coordinate', LOOP, 'second_switch', 'coordinate'),
              (93, 'actual_native_geometry_radius', LOOP, 'second_switch', 'radius'),
              (80, 'actual_native_physical_U', MOD, 'first_switch', 'U'),
              (84, 'actual_native_physical_Q', MOD, 'first_switch', 'Q'),
              (96, 'actual_native_physical_U', MOD, 'second_switch', 'U'),
              (100, 'actual_native_physical_Q', MOD, 'second_switch', 'Q')]))
    check('actual_J3_physical_measure_and_common_L',
          all(row_is(i, source_nodes=[9], source_record=LOOP) for i in (78, 79, 94, 95))
          and rows[367]['measure'] == rows[377]['measure'] == 'R^2*Jacobian*dnative'
          and all(weighted['component_result']['source_interface']['actual_L'].get(key) == value
                  for key, value in dict(source_record=MOD, positive_quotient_denominator=25272,
                      formula='1-delta*Z^2', Z_independent_delta_port=True).items()))
    check('accepted_uniform_B_creator_recipe',
          old['component_proofs']['J2_B']['source_ast_hashes']['active_S1_creator_sha256']
          == RECORD_HASHES[WEIGHTED_CREATOR]
          and old['native_weighted_B_definition']['formula'] ==
          'R^2*(Q/L+Q^2/2-2*Z^2*U^2/L^2)')
    kept = parent['retained_construction_data']
    check('retained_ledger', all(kept[key] == value for key, value in
          dict(cumulative_terms=212, raw_signed_terms=213, reduced_terms=183,
               six_signed_channels_retained=True, all_old_error_data_retained=True,
               reserved_contribution_registered_once=True, inherited_errors_deleted=False).items()))
    graph = [dict(id='original_switch_transport', rule='actual endpoint handoff and equal-bound original histories'),
             dict(id='R_U_transport', rule='source root and loop composition at the common R1 endpoint',
                  depends_on=['original_switch_transport']),
             dict(id='physical_m_Q_transport', rule='finite-N memory then physical m and positive-L quotient',
                  depends_on=['original_switch_transport', 'R_U_transport']),
             dict(id='B_C1_composition', rule='identical R U Q C1 ports in the uniform actual B recipe',
                  depends_on=['R_U_transport', 'physical_m_Q_transport']),
             dict(id='signed_J3', rule='actual right-lower minus left-upper [377,367]',
                  depends_on=['B_C1_composition']),
             dict(id='native_IN_reduction', rule='substitute proved zero J3 into accepted native source identity',
                  depends_on=['signed_J3'])]
    return dict(schema=SCHEMA, component='endpoint_transport', passed=True, source_only=True,
        **dict.fromkeys(STAGE_FLAGS, True), **dict.fromkeys(FALSE_FLAGS, False),
        stage_flags=dict.fromkeys(STAGE_FLAGS, True), nonclaims=dict.fromkeys(FALSE_FLAGS, False),
        J3_B_identity_proved=True, J3_R_identity_proved=True, J3_U_identity_proved=True, J3_Q_identity_proved=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, owned_Z_cell=['0', '1/2'],
        input_digest=expected_input_digest, source_family=FAMILY, source_function_identity_only=True,
        identity_domain_scope='common domain of original source definitions; no new loop-domain, cone or numerical admission',
        source_checks=checks, proof_graph=graph, proof_graph_sha256=digest(graph),
        component_proofs=dict(original_switch=original, R_U=ru, physical_m_Q=mq),
        endpoint_source_algebra_proof=algebra,
        source_binding_scope=old['source_binding_scope'], retained_construction_data=kept,
        quantity_distinctions=old['quantity_distinctions'], inherited_R_identity_jumps=old['inherited_R_identity_jumps'],
        B_identity_scope=dict(weighted_jump_node=523, jump_index=3, left_B_node=367, right_B_node=377,
            jump_arguments=[377, 367], ordinary_Z_orders=[0, 1],
            endpoint_substitutions={'73190': 1, '75949': 0}, common_L_node=25272,
            newly_proved_B_identity_nodes=[523], inherited_B_identity_nodes=[520, 521, 522],
            cumulative_B_identity_nodes=[520, 521, 522, 523],
            unresolved_B_identity_nodes=list(range(524, 536)), remaining_B_identity_count=12,
            all_native_B_jumps_zero_proved=False),
        unresolved_B_identity_nodes=list(range(524, 536)), remaining_B_identity_count=12,
        native_weighted_B_definition=old['native_weighted_B_definition'],
        B_C1_derivative_rule=old['component_proofs']['J2_B']['C1_derivative_rule'],
        native_IN_inherited_identity=old['native_IN_reduced_source_identity'],
        native_IN_reduced_source_identity='I_N=B_lastupper-B_firstlower-sum(J_4..J_15)+sum17(integral K_i dy)',
        native_IN_reduction_is_source_only=True, archived_all_16_J_terms_retained=True,
        source_jump_zero_does_not_delete_quantitative_errors=True,
        unresolved_gates=[], remaining_construction_gates=[
            'the other twelve native B jumps and same-source quantitative packet transport',
            'independent terminal functional moments and retained quantitative errors',
            'numerical containment and domain admission',
            'heat stress flat remainder time recursion oscillatory correction independent Cartesian validation'],
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_proof_functions_called=0, generic_DAG_walks=0, numerical_evaluations=0,
            quadratures=0, recurrence_calls=0, runtime_model_metadata=None), runtime_model_metadata=None)


# Integrated fresh source fragment: native_J3_m_Q_source_fragment_v3_2026_10_21.py
"""Guarded, source-only J3 m/Q transport leaf.

The two public functions inspect the caller-supplied decoded archive records,
typed chart bindings, fixed graph rows, and authenticated creator AST bodies.
They never execute a producer, proof, operator, checker, provider, or
normalizer.  Preparation of this file did not call either entry point.

The original C6 archive is a required ``background`` argument because the
LOOP snapshot intentionally omits the source-chart bindings.  The functions
return conditional source guards only; root acceptance remains separate.
"""
import ast as J3MQ__ast
import hashlib as J3MQ__hashlib
J3MQ_SCHEMA = 'native_J3_m_Q_source_fragment_v2_guarded'
J3MQ_STATUS = 'unexecuted_guarded_source_leaf'
J3MQ_RUNTIME_MODEL_METADATA = None
J3MQ_ACCEPTED_IDENTITY_CLAIM = False
J3MQ_PARENT_REPORT = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_J2_B_source_closure.json'
J3MQ_PARENT_REPORT_SHA256 = '37f20f105dcb58691b32deb7f5fca59f1a25ff7641086614aa61885ddc090655'
J3MQ_C6_CREATOR_SHA256 = '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32'
J3MQ_BRIDGE_CREATOR_SHA256 = '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605'
J3MQ_LOOP_CREATOR_SHA256 = '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d'
J3MQ_MOD_CREATOR_SHA256 = '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7'
J3MQ_C6_RAW_ARCHIVE_SHA256 = '5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a'
J3MQ_LOOP_RAW_ARCHIVE_SHA256 = '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa'
J3MQ_MOD_RAW_ARCHIVE_SHA256 = 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a'
J3MQ_C6_AST_SHA256 = {'first_switch_fields': 'ac140c0eaee3e5e3130d8bdd11fb5d27a3a25430233cd41d107dcd11ca23e6ac', 'second_switch_fields': 'af0fcc282e2569690f21f253e8b6a163390e63358db858c393fc0a447eb2115e', 'add_chart': '50194431065c07b29545a5f53c997f731f0ee79ee7e23e32641d789fae8ab2a9', 'end': 'e2fa8acbb0bd0d73ddd946fdb39e35d2542d1815c03528c381fda508fc10a460'}
J3MQ_BRIDGE_AST_SHA256 = {'histories': 'f868b63668f662a5c7a07f44a09f47c7b1724b19fb474518a0f62fea072469c7', 'generic_fields': '03ad6648602c0aafbd5a304b29fd7b8e0e25e16c6ebe61039ee6e52ee0a46dbc'}
J3MQ_LOOP_AST_SHA256 = {'finite_N_functions': '6221e9b2fb9ba6173f0693f0a51dcc868a1ff76a5a45f489fde115e9f74f270d', 'transport': '2b7815e4252db2157c9f6142bb0e3d80c125a2865f68cbefc5b051f4ac8e0141'}
J3MQ_MOD_AST_SHA256 = {'physical_pressure_fields': '06b5af9f5c843afaf2068b9af4d485fadff4351f4753d4838275db85e61857cc'}
J3MQ_FIRST_COORDINATE = 73190
J3MQ_SECOND_COORDINATE = 75949
J3MQ_FIRST_RADIUS = 85905
J3MQ_SECOND_RADIUS = 88600
J3MQ_PSTAR = 33572
J3MQ_PSTAR_INVERSE = 69245
J3MQ_COMMON_L = 25272
J3MQ__ORIGINAL_ENDPOINT_PROJECTIONS = {'first_switch_upper': {'coordinate_node': J3MQ_FIRST_COORDINATE, 'coordinate_value': 1, 'radius_node': J3MQ_FIRST_RADIUS}, 'second_switch_lower': {'coordinate_node': J3MQ_SECOND_COORDINATE, 'coordinate_value': 0, 'radius_node': J3MQ_SECOND_RADIUS}}
J3MQ__ORIGINAL_SOURCE_HASHES = {'parent_report_sha256': J3MQ_PARENT_REPORT_SHA256, 'c6_creator_sha256': J3MQ_C6_CREATOR_SHA256, 'bridge_creator_sha256': J3MQ_BRIDGE_CREATOR_SHA256, 'c6_raw_archive_sha256': J3MQ_C6_RAW_ARCHIVE_SHA256, 'loop_raw_archive_sha256': J3MQ_LOOP_RAW_ARCHIVE_SHA256}
J3MQ__NATIVE6_CERTIFICATE = {'channels': ['M', 'H', 'K', 'A', 'B', 'C'], 'ordinary_Z_orders': [0, 1, 2], 'first_switch_upper_endpoint_aliases': [[87392, 88465], [87352, 88458], [87427, 88472], [87462, 88479], [87497, 88486], [87525, 88493]], 'second_switch_lower_incoming_aliases': [[88465, 88465], [88458, 88458], [88472, 88472], [88479, 88479], [88486, 88486], [88493, 88493]]}
J3MQ__GENERIC5_CERTIFICATE = {'fields': ['m', 'h', 'k', 'e', 'p'], 'ordinary_Z_orders': [0, 1, 2], 'first_field_rows': {'m': [87591, 87592, 87593], 'h': [87632, 87633, 87634], 'k': [87673, 87674, 87675], 'e': [87728, 87729, 87730], 'p': [87735, 87738, 87742]}, 'second_field_rows': {'m': [90109, 90110, 90111], 'h': [90150, 90151, 90152], 'k': [90191, 90192, 90193], 'e': [90246, 90247, 90248], 'p': [90253, 90256, 90260]}}
J3MQ__ORIGINAL_GUARD_NAMES = ('decoded_background_and_loop_shapes', 'shared_authenticated_source_family', 'authenticated_C6_constructor_AST_and_calls', 'authenticated_bridge_histories_AST_and_calls', 'typed_C6_actual_moment_and_loop_chart_bindings', 'C6_coordinate_radius_Pstar_and_L_source_rows', 'C6_first_switch_upper_actual_U_endpoint_alias', 'C6_first_switch_upper_moment_endpoint_alias', 'C6_first_switch_original_M_history_through_Z2', 'C6_second_switch_M_history_through_Z2', 'C6_original_m_fields_use_inverse_Pstar_once', 'C6_second_switch_actual_phi_transport_rows', 'C6_first_switch_upper_all_native_moment_endpoint_aliases', 'C6_native_six_history_transport_M_H_K_A_B_C_through_Z2', 'C6_generic_h_k_e_p_rows_after_native_normalization', 'C6_common_W_FF_F0_Pstar_radius_source_rows', 'LOOP_shared_C6_prefix_fixed_source_rows')
J3MQ__RU_SOURCE_HASHES = {'creator_text': {'loop_creator_text': J3MQ_LOOP_CREATOR_SHA256, 'mod_creator_text': J3MQ_MOD_CREATOR_SHA256, 'c6_creator_text': J3MQ_C6_CREATOR_SHA256, 'bridge_creator_text': J3MQ_BRIDGE_CREATOR_SHA256}, 'raw_archives': {'loop': J3MQ_LOOP_RAW_ARCHIVE_SHA256, 'mod': J3MQ_MOD_RAW_ARCHIVE_SHA256, 'background': J3MQ_C6_RAW_ARCHIVE_SHA256, 'parent_bridge': '5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92'}, 'ast_body': {'c6.first_switch_fields': '0bd697b785e7241b2ea8520d74267342e4993f23dffc90fe49617d13430d68e2', 'c6.second_switch_fields': 'c6c08d15f7b5058f9ff48358ea904eed294d21046ca016deb7d1c440abd4d438', 'c6.CurrentOriginalPrescribedR100RmC6Functions.__init__': 'd7e974bcc4e2156f8c05b3a67192aefe1df25cb0857df1e7a47e8bd1f6eb2b32', 'c6.CurrentOriginalPrescribedR100RmC6Functions.add_chart': '18f1dbdd7128ac2e72098ab5d7249659923a5d7166640057b2587acfcf2f3812', 'c6.CurrentOriginalPrescribedR100RmC6Functions.end': '3a9ca663af71f7dcf34fa945662ea1535de41a0d8bc2a53a2aacf2954d5b6040', 'loop.loop_functions': 'fec28c8c3ff3bf20a5cc64bc777b0b0ed304387b95f71cb6a8f9ec30fd9e0559', 'loop.generic_roots': '3ba21f81b7aae33751e6534e15e97311c1b851f103b522259a80d1170fad7847', 'loop.finite_N_functions': '219500c1ac0665e869303a421303a9d2a4c05bfa3121a9fe1f06f022649ed4d2', 'loop.transport': 'fd5f45c38d246f2555c045ca2d386f5b49bbfaa7699981156efbf99c1869e9b3', 'loop.CurrentOriginalPrescribedLoopC5Functions.__init__': '4e6a13cf1025d85e3bf0d7c5abfe9b90d89e9cf45692ee6f5b6161b4ed8546c2', 'mod.physical_pressure_fields': '619f51980979fc35ad965386a0038fc2785bf15c99ef2a1b61d749a2d6833a14', 'mod.CurrentOriginalModifiedQ4PressureC3Functions.__init__': '8040fa64b5da2de7fa1e2e29cd2583bb803b5b0735d36357050200ad87d645d4', 'bridge.histories': '9dc010e43b2a9b1c498380b4452f226f467df7f5529f4c026ed5b8f7a6aaaf5d', 'bridge.generic_fields': '41f3b37e36843799b913ec6af1da1d59ef8dcd668c34a0b3584621fe7e791cdf', 'bridge.CurrentOriginalPrescribedBridgeC6Functions.__init__': 'df1b4b43d63a07bc8e44aa8c5db8bed8c09d73ac31626b3befbd69e750f93e5b'}}
J3MQ__RU_ENDPOINT_PAIRS = {'c6_R1_and_second_lower': [[88437, 88438]], 'c5_first_V_N': [[106932, 106933]], 'c5_second_V_N': [[109340, 109341]], 'c3_first_physical_U': [[148085, 148086]], 'c3_second_physical_U': [[154247, 154248]], 'common_Pstar': J3MQ_PSTAR, 'common_L': J3MQ_COMMON_L, 'common_N': 18}

def J3MQ__compact(text):
    return ''.join(str(text).split())

def J3MQ__function_nodes(text, names):
    """Parse creator text and return exactly named AST function nodes."""
    if not isinstance(text, str):
        return (False, {}, '')
    try:
        tree = J3MQ__ast.parse(text)
    except (SyntaxError, TypeError, ValueError):
        return (False, {}, J3MQ__compact(text))
    found = {name: [] for name in names}
    for node in J3MQ__ast.walk(tree):
        if isinstance(node, (J3MQ__ast.FunctionDef, J3MQ__ast.AsyncFunctionDef)):
            if node.name in found:
                found[node.name].append(node)
    return (True, found, J3MQ__compact(text))

def J3MQ__ast_pin(text, expected_sha, function_hashes, snippets):
    names = tuple(function_hashes)
    parsed, found, compact = J3MQ__function_nodes(text, names)
    if not parsed:
        return False
    if J3MQ__hashlib.sha256(text.encode()).hexdigest() != expected_sha:
        return False
    for name, expected in function_hashes.items():
        nodes = found.get(name, ())
        if len(nodes) != 1:
            return False
        dump = J3MQ__ast.dump(nodes[0], annotate_fields=True, include_attributes=False)
        if J3MQ__hashlib.sha256(dump.encode()).hexdigest() != expected:
            return False
    return all((needle in compact for needle in snippets))

def J3MQ__rows(archive):
    return archive.get('exact_graph_nodes', ()) if isinstance(archive, dict) else ()

def J3MQ__row(rows, number):
    if not isinstance(rows, list) or not isinstance(number, int):
        return None
    if number < 0 or number >= len(rows):
        return None
    value = rows[number]
    return value if isinstance(value, dict) else None

def J3MQ__row_is(rows, number, **expected):
    actual = J3MQ__row(rows, number)
    return actual is not None and all((actual.get(key) == value for key, value in expected.items()))

def J3MQ__array(value, prefix):
    return isinstance(value, list) and value[:len(prefix)] == list(prefix)

def J3MQ__emit(check, name, value):
    value = bool(value)
    if callable(check):
        check(name, value)
    elif isinstance(check, dict):
        check[name] = value
    elif hasattr(check, 'record') and callable(check.record):
        check.record(name, value)
    else:
        raise TypeError('check must be a callback, dict, or recorder')
    return value

def J3MQ__guard(checks, check, name, value):
    value = bool(value)
    J3MQ__emit(check, name, value)
    checks[name] = value
    return value

def J3MQ__safe(condition):
    try:
        return bool(condition())
    except (KeyError, IndexError, TypeError, AttributeError, ValueError):
        return False

def J3MQ__integral_ok(row, variable, upper, order, *, rate=None, history_name=None, bridge_chart=None, finite=False):
    if not isinstance(row, dict):
        return False
    expected = {'operation': 'definite_integral', 'variable': variable, 'lower': 0, 'upper': upper, 'exact_function_integral': True, 'numerical_value_not_installed': True, 'endpoint_independent_of_Z': True, 'ordinary_slow_Z_derivative_order': order, 'true_core_inlet_complete_history': True, 'numerical_cap_not_defining_integrand': True}
    if any((row.get(key) != value for key, value in expected.items())):
        return False
    if finite:
        return all((row.get(key) == value for key, value in {'actual_finite_N_signed_correction': True, 'physical_Jacobian_applied_once': True, 'zero_correction_inlet_only_at_original_rminus': True, 'predecessor_memory_never_reset': True}.items()))
    return all((row.get(key) == value for key, value in {'actual_six_history_rate': rate, 'actual_history_name': history_name, 'original_incoming_memory_retained': True, 'source_measure_applied_once': True, 'prescribed_bridge_chart': bridge_chart}.items()))

def J3MQ__c6_source_ast_ok(c6_text, bridge_text):
    c6_needles = ("initial={key:bridge_owner.R100['actual'][key]forkeyin('phi','U')}", "self.endpoints['R1']=self.end('first_switch',g.one)", "incoming=self.endpoints['R1']", "actual=second_switch_fields(a,known,incoming['actual'],t,hb,R)", "self.add_chart('second_switch',t,R,hb,actual,incoming['moments'],known=known)", 'logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase,original_second_switch_full_current_radius_direction_and_power=True)', "returndict(phi=phi,U=initial['U']", "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)", "actual={key:at(chart['actual'][key])")
    bridge_needles = ('defhistories(a,phi,U,initial,phase,jacobian,label):', 'forname,rateinRATES.items():', 'out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,', 'actual_six_history_rate=rate', 'original_incoming_memory_retained=True', 'source_measure_applied_once=True', 'prescribed_bridge_chart=label', 'defgeneric_fields(a,actual,J,F0,P0,Pstar,radius,delta):', "W=a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))", 'FF=a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))', "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)")
    c6_ok = J3MQ__ast_pin(c6_text, J3MQ_C6_CREATOR_SHA256, J3MQ_C6_AST_SHA256, c6_needles)
    bridge_ok = J3MQ__ast_pin(bridge_text, J3MQ_BRIDGE_CREATOR_SHA256, J3MQ_BRIDGE_AST_SHA256, bridge_needles)
    return (c6_ok, bridge_ok)

def J3MQ__loop_source_ast_ok(loop_text):
    needles = ('deffinite_N_functions(g,roots,loop,N):', "s=a.scale(loop['B_over_Pstar'],invN)", 'five_signed_density_increments=inc', 'deftransport(g,densities,incoming,coordinate,lower,jacobian,label):', 'width=g.mul(jacobian,g.sub(coordinate,lower))', 'contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True,', 'physical_Jacobian_applied_once=True', 'predecessor_memory_never_reset=True', 'out[name]=a.add(a.scale(incoming[name],decay),contribution)', 'incoming=self.endpoints[name]', "fields[key]=a.add(a.trim(chart['fields'][key]),history[key])", 'self.endpoints[name]={key:a.at(jet,coordinate,upper)forkey,jetinhistory.items()}')
    return J3MQ__ast_pin(loop_text, J3MQ_LOOP_CREATOR_SHA256, J3MQ_LOOP_AST_SHA256, needles)

def J3MQ__mod_source_ast_ok(mod_text):
    needles = ('defphysical_pressure_fields(g,fields,radial,delta,Pstar):', "U=a.scale(fields['V'],Pstar)", "m=a.scale(fields['m'],Pstar)", 'L=b.sub(one,b.mul(de,b.mul(z,z)))', 'sameoriginalL=1-delta*Z^2', 'Q=div(numerator(U,m));Qy=div(numerator(Uy,my))', 'returndict(U=U,Uy=Uy,m=m,my=my,Q=Q,Qy=Qy)')
    return J3MQ__ast_pin(mod_text, J3MQ_MOD_CREATOR_SHA256, J3MQ_MOD_AST_SHA256, needles)

def J3MQ__archive_shapes(background, loop, mod=None):
    if not isinstance(background, dict) or not isinstance(loop, dict):
        return False
    if not isinstance(background.get('exact_graph_nodes'), list):
        return False
    if not isinstance(loop.get('exact_graph_nodes'), list):
        return False
    if mod is not None and (not isinstance(mod, dict)):
        return False
    return True

def J3MQ__shared_source_family(*archives):
    """Require the authenticated three-digest family object on every archive."""
    if not archives or not all((isinstance(archive, dict) for archive in archives)):
        return False
    families = [archive.get('source_family') for archive in archives]
    family = families[0]
    return isinstance(family, dict) and len(family) == 3 and all((isinstance(key, str) and isinstance(value, str) for key, value in family.items())) and all((other == family for other in families[1:]))

def J3MQ__c6_typed_bindings(background, loop):
    try:
        cb = background['prescribed_R100_Rm_C6_bindings']
        lb = loop['prescribed_loop_C5_bindings']
        first = cb['charts']['first_switch']
        second = cb['charts']['second_switch']
        r1 = cb['endpoints']['R1']
        first_fields = first['fields']
        second_fields = second['fields']
        lf = lb['charts']['first_switch']
        ls = lb['charts']['second_switch']
        native_prefixes = {'H': ([87352, 87353, 87354], [89870, 89871, 89872]), 'M': ([87392, 87393, 87394], [89910, 89911, 89912]), 'K': ([87427, 87428, 87429], [89945, 89946, 89947]), 'A': ([87462, 87463, 87464], [89980, 89981, 89982]), 'B': ([87497, 87498, 87499], [90015, 90016, 90017]), 'C': ([87525, 87526, 87527], [90043, 90044, 90045])}
        incoming_prefixes = {'H': ([73043, 73044, 73045], [88458, 88459, 88460]), 'M': ([73050, 73051, 73052], [88465, 88466, 88467]), 'K': ([73057, 73058, 73059], [88472, 88473, 88474]), 'A': ([73064, 73065, 73066], [88479, 88480, 88481]), 'B': ([73071, 73072, 73073], [88486, 88487, 88488]), 'C': ([73078, 73079, 73080], [88493, 88494, 88495])}
        generic_prefixes = {'m': ([87591, 87592, 87593], [90109, 90110, 90111]), 'h': ([87632, 87633, 87634], [90150, 90151, 90152]), 'k': ([87673, 87674, 87675], [90191, 90192, 90193]), 'e': ([87728, 87729, 87730], [90246, 90247, 90248]), 'p': ([87735, 87738, 87742], [90253, 90256, 90260])}
        native_ok = all((J3MQ__array(first['moments'][name], prefixes[0]) and J3MQ__array(second['moments'][name], prefixes[1]) and J3MQ__array(first['incoming_histories'][name], incoming_prefixes[name][0]) and J3MQ__array(second['incoming_histories'][name], incoming_prefixes[name][1]) and J3MQ__array(r1['moments'][name], incoming_prefixes[name][1]) and (r1['moments'][name][:3] == second['incoming_histories'][name][:3]) for name, prefixes in native_prefixes.items()))
        generic_ok = all((J3MQ__array(first_fields[name], prefixes[0]) and J3MQ__array(second_fields[name], prefixes[1]) for name, prefixes in generic_prefixes.items()))
        p0 = [32693, 32694, 32695, 32696, 66509, 66510, 66511]
        if not (cb.get('Pstar') == J3MQ_PSTAR and J3MQ__array(cb.get('F0'), [35657, 35658, 35661]) and J3MQ__array(cb.get('P0'), p0) and (first.get('phase') == J3MQ_FIRST_COORDINATE) and (first.get('domain') == [0, 1]) and (first.get('radius') == J3MQ_FIRST_RADIUS) and (first.get('physical_measure') == 9) and (second.get('phase') == J3MQ_SECOND_COORDINATE) and (second.get('domain') == [0, 1]) and (second.get('radius') == J3MQ_SECOND_RADIUS) and (second.get('physical_measure') == 9) and J3MQ__array(r1['actual']['U'], [88437, 88438, 88439]) and J3MQ__array(r1['moments']['M'], [88465, 88466, 88467]) and J3MQ__array(first['actual']['U'], [87154, 87155, 87156]) and J3MQ__array(first['moments']['M'], [87392, 87393, 87394]) and J3MQ__array(second['actual']['U'], [88437, 88438, 88439]) and J3MQ__array(second['incoming_histories']['M'], [88465, 88466, 88467]) and (lf.get('coordinate') == J3MQ_FIRST_COORDINATE) and (lf.get('radius') == J3MQ_FIRST_RADIUS) and (lf.get('domain') == [0, 1]) and (ls.get('coordinate') == J3MQ_SECOND_COORDINATE) and (ls.get('radius') == J3MQ_SECOND_RADIUS) and (ls.get('domain') == [0, 1]) and (first_fields.get('P0') == p0) and (second_fields.get('P0') == p0) and native_ok and generic_ok):
            return False
    except (KeyError, IndexError, TypeError):
        return False
    return True

def J3MQ__c6_rows(background):
    rows = J3MQ__rows(background)
    basic = all((J3MQ__row_is(rows, 73190, operation='bound_variable', name='true_prescribed_first_switch_s'), J3MQ__row_is(rows, 75949, operation='bound_variable', name='true_prescribed_second_switch_t'), J3MQ__row_is(rows, 85905, operation='product', arguments=[31709, 85904]), J3MQ__row_is(rows, 88600, operation='product', arguments=[31709, 88599]), J3MQ__row_is(rows, 33572, operation='analytic_unary', name='exp', argument=5), J3MQ__row_is(rows, 69245, operation='positive_quotient', numerator=1, denominator=33572, source_positive_certificate='same strictly positive Pstar'), J3MQ__row_is(rows, 25272, operation='sum', arguments=[1, 25271])))
    r1_alias = all((J3MQ__row_is(rows, endpoint, operation='function_substitution', expression=expression, variable=J3MQ_FIRST_COORDINATE, value=1, Z_independent_substitution=True) for endpoint, expression in zip((88437, 88438, 88439), (87154, 87155, 87156))))
    r1_m_alias = all((J3MQ__row_is(rows, endpoint, operation='function_substitution', expression=expression, variable=J3MQ_FIRST_COORDINATE, value=1, Z_independent_substitution=True) for endpoint, expression in zip((88465, 88466, 88467), (87392, 87393, 87394))))
    r1_native_alias = all((J3MQ__row_is(rows, endpoint, operation='function_substitution', expression=expression, variable=J3MQ_FIRST_COORDINATE, value=1, Z_independent_substitution=True) for endpoint, expression in ((88458, 87352), (88459, 87353), (88460, 87354), (88465, 87392), (88466, 87393), (88467, 87394), (88472, 87427), (88473, 87428), (88474, 87429), (88479, 87462), (88480, 87463), (88481, 87464), (88486, 87497), (88487, 87498), (88488, 87499), (88493, 87525), (88494, 87526), (88495, 87527))))
    first_integrals = all((J3MQ__integral_ok(J3MQ__row(rows, node), 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, order, rate=1, history_name='M', bridge_chart='true_first_switch') for node, order in zip((87385, 87386, 87387), (0, 1, 2))))
    first_m = all((J3MQ__row_is(rows, node, operation='sum', arguments=[initial, integral]) for node, initial, integral in zip((87392, 87393, 87394), (87378, 87379, 87380), (87385, 87386, 87387))))
    second_logratio = all((J3MQ__row_is(rows, node, operation='definite_integral', variable='true_second_switch_angular_t', lower=0, upper=J3MQ_SECOND_COORDINATE, exact_function_integral=True, numerical_value_not_installed=True, endpoint_independent_of_Z=True, ordinary_slow_Z_derivative_order=order, true_core_inlet_complete_history=True, numerical_cap_not_defining_integrand=True, original_second_switch_full_current_radius_direction_and_power=True) for node, order in zip((89610, 89611, 89612), (0, 1, 2))))
    second_integrals = all((J3MQ__integral_ok(J3MQ__row(rows, node), 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, order, rate=1, history_name='M', bridge_chart='true_second_switch') for node, order in zip((89903, 89904, 89905), (0, 1, 2))))
    second_m = all((J3MQ__row_is(rows, node, operation='sum', arguments=[decay, integral]) for node, decay, integral in zip((89910, 89911, 89912), (89896, 89897, 89898), (89903, 89904, 89905))))
    incoming_m = all((J3MQ__row_is(rows, node, operation='product', arguments=[incoming, 89878]) for node, incoming in zip((89896, 89897, 89898), (88465, 88466, 88467))))
    second_phi = J3MQ__row_is(rows, 89644, operation='product', arguments=[88430, 89617]) and J3MQ__row_is(rows, 89645, operation='product', arguments=[88430, 89618]) and J3MQ__row_is(rows, 89646, operation='product', arguments=[88431, 89617]) and J3MQ__row_is(rows, 89647, operation='sum', arguments=[89645, 89646]) and J3MQ__row_is(rows, 89648, operation='product', arguments=[88430, 89621]) and J3MQ__row_is(rows, 89649, operation='product', arguments=[48, 88431, 89618]) and J3MQ__row_is(rows, 89650, operation='product', arguments=[88432, 89617]) and J3MQ__row_is(rows, 89651, operation='sum', arguments=[89648, 89649, 89650])
    m_scale = all((J3MQ__row_is(rows, node, operation='product', arguments=[moment, J3MQ_PSTAR_INVERSE]) for node, moment in zip((87591, 87592, 87593), (87392, 87393, 87394)))) and all((J3MQ__row_is(rows, node, operation='product', arguments=[moment, J3MQ_PSTAR_INVERSE]) for node, moment in zip((90109, 90110, 90111), (89910, 89911, 89912))))
    return (basic, r1_alias, r1_m_alias, first_integrals and first_m, second_logratio and second_integrals and second_m and incoming_m, m_scale, second_phi, r1_native_alias)

def J3MQ__c6_native_history_rows(background):
    """Check the six native bridge histories at both switch endpoints, Z 0..2."""
    rows = J3MQ__rows(background)
    specs = (('H', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73043, 73044, 73045), 87320, (87338, 87339, 87340), (87345, 87346, 87347), (87352, 87353, 87354), 2), ('M', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73050, 73051, 73052), 87360, (87378, 87379, 87380), (87385, 87386, 87387), (87392, 87393, 87394), 1), ('K', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73057, 73058, 73059), 87320, (87413, 87414, 87415), (87420, 87421, 87422), (87427, 87428, 87429), 2), ('A', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73064, 73065, 73066), 87360, (87448, 87449, 87450), (87455, 87456, 87457), (87462, 87463, 87464), 1), ('B', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73071, 73072, 73073), 87320, (87483, 87484, 87485), (87490, 87491, 87492), (87497, 87498, 87499), 2), ('C', 'true_first_switch', 'true_first_switch_history_t', J3MQ_FIRST_COORDINATE, (73078, 73079, 73080), 87360, (87511, 87512, 87513), (87518, 87519, 87520), (87525, 87526, 87527), 1), ('H', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88458, 88459, 88460), 89838, (89856, 89857, 89858), (89863, 89864, 89865), (89870, 89871, 89872), 2), ('M', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88465, 88466, 88467), 89878, (89896, 89897, 89898), (89903, 89904, 89905), (89910, 89911, 89912), 1), ('K', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88472, 88473, 88474), 89838, (89931, 89932, 89933), (89938, 89939, 89940), (89945, 89946, 89947), 2), ('A', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88479, 88480, 88481), 89878, (89966, 89967, 89968), (89973, 89974, 89975), (89980, 89981, 89982), 1), ('B', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88486, 88487, 88488), 89838, (90001, 90002, 90003), (90008, 90009, 90010), (90015, 90016, 90017), 2), ('C', 'true_second_switch', 'true_second_switch_history_t', J3MQ_SECOND_COORDINATE, (88493, 88494, 88495), 89878, (90029, 90030, 90031), (90036, 90037, 90038), (90043, 90044, 90045), 1))
    for name, chart, variable, upper, incoming, decay, products, integrals, final, rate in specs:
        if not all((J3MQ__row_is(rows, product, operation='product', arguments=[source, decay]) for product, source in zip(products, incoming))):
            return False
        if not all((J3MQ__integral_ok(J3MQ__row(rows, node), variable, upper, order, rate=rate, history_name=name, bridge_chart=chart) for node, order in zip(integrals, (0, 1, 2)))):
            return False
        if not all((J3MQ__row_is(rows, node, operation='sum', arguments=[product, integral]) for node, product, integral in zip(final, products, integrals))):
            return False
    return True

def J3MQ__c6_generic_field_rows(background):
    """Check bounded generic h/k/e/p source rows after native normalization."""
    rows = J3MQ__rows(background)
    h_ok = all((J3MQ__row_is(rows, node, operation='product', arguments=[source, 28]) for node, source in ((87632, 87598), (87633, 87601), (87634, 87605), (90150, 90116), (90151, 90119), (90152, 90123))))
    k_ok = all((J3MQ__row_is(rows, node, operation='product', arguments=[source, 69387]) for node, source in ((87673, 87639), (87674, 87642), (87675, 87646), (90191, 90157), (90192, 90160), (90193, 90164))))
    e_ok = all((J3MQ__row_is(rows, node, operation='negative', argument=argument) for node, argument in ((87721, 87687), (87722, 87690), (87723, 87694), (90239, 90205), (90240, 90208), (90241, 90212)))) and all((J3MQ__row_is(rows, node, operation='sum', arguments=[source, negative]) for node, source, negative in ((87728, 87680, 87721), (87729, 87681, 87722), (87730, 87682, 87723), (90246, 90198, 90239), (90247, 90199, 90240), (90248, 90200, 90241))))
    p_ok = all((J3MQ__row_is(rows, node, operation=operation, arguments=arguments) for node, operation, arguments in ((87735, 'product', [87543, 87525]), (87738, 'sum', [87736, 87737]), (87742, 'sum', [87739, 87740, 87741]), (90253, 'product', [90061, 90043]), (90256, 'sum', [90254, 90255]), (90260, 'sum', [90257, 90258, 90259]))))
    return h_ok and k_ok and e_ok and p_ok

def J3MQ__c6_common_scale_rows(background):
    """Check the shared F0/Pstar/radius W and FF source rows at both seams."""
    rows = J3MQ__rows(background)
    return all((J3MQ__row_is(rows, 87532, operation='product', arguments=[48, J3MQ_FIRST_RADIUS]), J3MQ__row_is(rows, 87533, operation='analytic_unary', name='sqrt', argument=87532), J3MQ__row_is(rows, 87534, operation='product', arguments=[87533, J3MQ_PSTAR_INVERSE]), J3MQ__row_is(rows, 87535, operation='product', arguments=[35657, 87534]), J3MQ__row_is(rows, 87536, operation='product', arguments=[35658, 87534]), J3MQ__row_is(rows, 87537, operation='product', arguments=[35661, 87534]), J3MQ__row_is(rows, 87542, operation='product', arguments=[J3MQ_FIRST_RADIUS, J3MQ_PSTAR_INVERSE, J3MQ_PSTAR_INVERSE]), J3MQ__row_is(rows, 87543, operation='product', arguments=[35693, 87542]), J3MQ__row_is(rows, 87544, operation='product', arguments=[35696, 87542]), J3MQ__row_is(rows, 87545, operation='product', arguments=[35700, 87542]), J3MQ__row_is(rows, 90050, operation='product', arguments=[48, J3MQ_SECOND_RADIUS]), J3MQ__row_is(rows, 90051, operation='analytic_unary', name='sqrt', argument=90050), J3MQ__row_is(rows, 90052, operation='product', arguments=[90051, J3MQ_PSTAR_INVERSE]), J3MQ__row_is(rows, 90053, operation='product', arguments=[35657, 90052]), J3MQ__row_is(rows, 90054, operation='product', arguments=[35658, 90052]), J3MQ__row_is(rows, 90055, operation='product', arguments=[35661, 90052]), J3MQ__row_is(rows, 90060, operation='product', arguments=[J3MQ_SECOND_RADIUS, J3MQ_PSTAR_INVERSE, J3MQ_PSTAR_INVERSE]), J3MQ__row_is(rows, 90061, operation='product', arguments=[35693, 90060]), J3MQ__row_is(rows, 90062, operation='product', arguments=[35696, 90060]), J3MQ__row_is(rows, 90063, operation='product', arguments=[35700, 90060])))

def J3MQ__loop_typed_bindings(loop):
    try:
        lb = loop['prescribed_loop_C5_bindings']
        first = lb['charts']['first_switch']
        second = lb['charts']['second_switch']
        first_end = lb['endpoints']['first_switch']
        fields = ('m', 'h', 'k', 'e', 'p')
        expected_endpoint = {'m': [107131, 107132, 107133], 'h': [107137, 107138, 107139], 'k': [107143, 107144, 107145], 'e': [107149, 107150, 107151], 'p': [107155, 107156, 107157]}
        chart_shape = first.get('coordinate') == J3MQ_FIRST_COORDINATE and first.get('radius') == J3MQ_FIRST_RADIUS and (first.get('domain') == [0, 1]) and (first.get('physical_measure') == 9) and (second.get('coordinate') == J3MQ_SECOND_COORDINATE) and (second.get('radius') == J3MQ_SECOND_RADIUS) and (second.get('domain') == [0, 1]) and (second.get('physical_measure') == 9) and (first.get('same_background_pressure_datum') is True) and (second.get('same_background_pressure_datum') is True) and (first.get('total_Z_phase_derivative_exact_zero') is True) and (second.get('total_Z_phase_derivative_exact_zero') is True)
        aliases = all((J3MQ__array(first_end[field], expected_endpoint[field]) and J3MQ__array(second['incoming'][field], expected_endpoint[field]) for field in fields))
        return (chart_shape and aliases, expected_endpoint)
    except (KeyError, IndexError, TypeError):
        return (False, {})
J3MQ__CORRECTION_SPECS = {'m': {'rate': '1', 'first_integral': 106958, 'first_history': 106970, 'first_decay': 106964, 'first_memory': 106940, 'second_integral': 109367, 'second_history': 109379, 'second_decay': 109373, 'second_memory': 109349, 'first_modified': [107095, 107096, 107097], 'first_base': [87591, 87592, 87593], 'second_modified': [109504, 109505, 109506], 'second_base': [90109, 90110, 90111], 'first_incoming': [104503, 104504, 104505], 'second_incoming': [107131, 107132, 107133]}, 'h': {'rate': '3/2', 'first_integral': 106993, 'first_history': 107005, 'first_decay': 106999, 'first_memory': 106977, 'second_integral': 109402, 'second_history': 109414, 'second_decay': 109408, 'second_memory': 109386, 'first_modified': [107101, 107102, 107103], 'first_base': [87601, 87602, 87603], 'second_modified': [109510, 109511, 109512], 'second_base': [90119, 90120, 90121], 'first_incoming': [104509, 104510, 104511], 'second_incoming': [107137, 107138, 107139]}, 'k': {'rate': '3/2', 'first_integral': 107023, 'first_history': 107035, 'first_decay': 107029, 'first_memory': 106977, 'second_integral': 109432, 'second_history': 109444, 'second_decay': 109438, 'second_memory': 109386, 'first_modified': [107107, 107108, 107109], 'first_base': [87611, 87612, 87613], 'second_modified': [109516, 109517, 109518], 'second_base': [90129, 90130, 90131], 'first_incoming': [104515, 104516, 104517], 'second_incoming': [107143, 107144, 107145]}, 'e': {'rate': '1', 'first_integral': 107053, 'first_history': 107065, 'first_decay': 107059, 'first_memory': 106940, 'second_integral': 109462, 'second_history': 109474, 'second_decay': 109468, 'second_memory': 109349, 'first_modified': [107113, 107114, 107115], 'first_base': [87621, 87622, 87623], 'second_modified': [109522, 109523, 109524], 'second_base': [90139, 90140, 90141], 'first_incoming': [104521, 104522, 104523], 'second_incoming': [107149, 107150, 107151]}, 'p': {'rate': '0', 'first_integral': 107083, 'first_history': 107089, 'first_decay': 107089, 'first_memory': 1, 'second_integral': 109492, 'second_history': 109498, 'second_decay': 109498, 'second_memory': 1, 'first_modified': [107119, 107120, 107121], 'first_base': [87631, 87632, 87633], 'second_modified': [109528, 109529, 109530], 'second_base': [90149, 90150, 90151], 'first_incoming': [104527, 104528, 104529], 'second_incoming': [107155, 107156, 107157]}}

def J3MQ__finite_N_rows(loop, endpoint_aliases):
    rows = J3MQ__rows(loop)
    fields_ok = True
    first_history_rows = {}
    second_history_rows = {}
    for field, spec in J3MQ__CORRECTION_SPECS.items():
        first_history_rows[field] = [spec['first_history'] + i for i in range(3)]
        second_history_rows[field] = [spec['second_history'] + i for i in range(3)]
        for order in range(3):
            fi = spec['first_integral'] + order
            si = spec['second_integral'] + order
            fields_ok = fields_ok and J3MQ__integral_ok(J3MQ__row(rows, fi), 'true_loop_C5_transport_first_switch', J3MQ_FIRST_COORDINATE, order, finite=True) and (J3MQ__row(rows, fi).get('original_rate') == spec['rate']) and (J3MQ__row(rows, fi).get('original_chart') == 'first_switch')
            fields_ok = fields_ok and J3MQ__integral_ok(J3MQ__row(rows, si), 'true_loop_C5_transport_second_switch', J3MQ_SECOND_COORDINATE, order, finite=True) and (J3MQ__row(rows, si).get('original_rate') == spec['rate']) and (J3MQ__row(rows, si).get('original_chart') == 'second_switch')
            first_head = spec['first_decay'] + order if spec['rate'] != '0' else spec['first_incoming'][order]
            second_head = spec['second_decay'] + order if spec['rate'] != '0' else spec['second_incoming'][order]
            fields_ok = fields_ok and J3MQ__row_is(rows, spec['first_history'] + order, operation='sum', arguments=[first_head, spec['first_integral'] + order])
            fields_ok = fields_ok and J3MQ__row_is(rows, spec['second_history'] + order, operation='sum', arguments=[second_head, spec['second_integral'] + order])
            if spec['rate'] != '0':
                fields_ok = fields_ok and J3MQ__row_is(rows, spec['first_decay'] + order, operation='product', arguments=[spec['first_incoming'][order], spec['first_memory']])
                fields_ok = fields_ok and J3MQ__row_is(rows, spec['second_decay'] + order, operation='product', arguments=[spec['second_incoming'][order], spec['second_memory']])
            else:
                fields_ok = fields_ok and spec['first_decay'] == spec['first_history']
                fields_ok = fields_ok and spec['second_decay'] == spec['second_history']
        endpoint = endpoint_aliases[field]
        for order, node in enumerate(endpoint):
            fields_ok = fields_ok and J3MQ__row_is(rows, node, operation='function_substitution', expression=spec['first_history'] + order, variable=J3MQ_FIRST_COORDINATE, value=1, Z_independent_substitution=True)
    return (fields_ok, first_history_rows, second_history_rows)

def J3MQ__physical_m_rows(mod):
    rows = J3MQ__rows(mod)
    direct = all((J3MQ__row_is(rows, node, operation='product', arguments=[source, J3MQ_PSTAR]) for node, source in zip((148097, 148098, 148099), (107095, 107096, 107097)))) and all((J3MQ__row_is(rows, node, operation='product', arguments=[source, J3MQ_PSTAR]) for node, source in zip((154259, 154260, 154261), (109504, 109505, 109506))))
    modified = all((J3MQ__row_is(rows, node, operation='sum', arguments=[base, history]) for node, base, history in zip((107095, 107096, 107097), (87591, 87592, 87593), (106970, 106971, 106972)))) and all((J3MQ__row_is(rows, node, operation='sum', arguments=[base, history]) for node, base, history in zip((109504, 109505, 109506), (90109, 90110, 90111), (109379, 109380, 109381))))
    return direct and modified

def J3MQ__physical_u_rows(mod):
    rows = J3MQ__rows(mod)
    return all((J3MQ__row_is(rows, node, operation='product', arguments=[source, J3MQ_PSTAR]) for node, source in zip((148085, 148086), (106932, 106933)))) and all((J3MQ__row_is(rows, node, operation='product', arguments=[source, J3MQ_PSTAR]) for node, source in zip((154247, 154248), (109340, 109341))))

def J3MQ__q_rows(mod):
    rows = J3MQ__rows(mod)
    roots = all((J3MQ__row_is(rows, node, operation='positive_quotient', numerator=numerator, denominator=J3MQ_COMMON_L, source_positive_certificate='same original L=1-delta*Z^2') for node, numerator in ((148180, 148175), (148184, 148183), (154342, 154337), (154346, 154345))))
    numerators = all((J3MQ__row_is(rows, node, operation='sum', arguments=args) for node, args in ((148175, [148127, 148170]), (148183, [148176, 148182]), (154337, [154289, 154332]), (154345, [154338, 154344]))))
    common = J3MQ__row_is(rows, J3MQ_COMMON_L, operation='sum', arguments=[1, 25271])
    return roots and numerators and common

def J3MQ__physical_bindings(mod):
    try:
        mb = mod['modified_Q4_pressure_C3_bindings']
        first = mb['charts']['first_switch']['physical']
        second = mb['charts']['second_switch']['physical']
        return J3MQ__array(first['U'], [148085, 148086]) and J3MQ__array(first['m'], [148097, 148098, 148099]) and J3MQ__array(first['Q'], [148180, 148184]) and J3MQ__array(second['U'], [154247, 154248]) and J3MQ__array(second['m'], [154259, 154260, 154261]) and J3MQ__array(second['Q'], [154342, 154346])
    except (KeyError, IndexError, TypeError):
        return False

def prove_original_switch_transport(background, loop, c6_creator_text, bridge_creator_text, check):
    """Guard original C6 actual/moment/m source transport through Z order 2.

    ``background`` is the decoded C6 archive and ``loop`` is the decoded LOOP
    archive.  Both are inspected directly; no producer or proof function is
    called.  The result is a conditional source theorem record for the later
    finite-N physical m/Q leaf.
    """
    checks = {}
    shapes = J3MQ__archive_shapes(background, loop)
    J3MQ__guard(checks, check, 'decoded_background_and_loop_shapes', shapes)
    source_family = shapes and J3MQ__shared_source_family(background, loop)
    J3MQ__guard(checks, check, 'shared_authenticated_source_family', source_family)
    c6_ast, bridge_ast = J3MQ__c6_source_ast_ok(c6_creator_text, bridge_creator_text)
    J3MQ__guard(checks, check, 'authenticated_C6_constructor_AST_and_calls', c6_ast)
    J3MQ__guard(checks, check, 'authenticated_bridge_histories_AST_and_calls', bridge_ast)
    typed = J3MQ__safe(lambda: J3MQ__c6_typed_bindings(background, loop))
    J3MQ__guard(checks, check, 'typed_C6_actual_moment_and_loop_chart_bindings', typed)
    if shapes:
        try:
            c6_basic, c6_alias, c6_m_alias, c6_first_m, c6_second_m, c6_m_scale, c6_second_phi, c6_native_alias = J3MQ__c6_rows(background)
            c6_native = J3MQ__c6_native_history_rows(background)
            c6_generic = J3MQ__c6_generic_field_rows(background)
            c6_scales = J3MQ__c6_common_scale_rows(background)
        except (KeyError, IndexError, TypeError, AttributeError, ValueError):
            c6_basic, c6_alias, c6_m_alias, c6_first_m, c6_second_m, c6_m_scale, c6_second_phi, c6_native_alias = (False, False, False, False, False, False, False, False)
            c6_native, c6_generic, c6_scales = (False, False, False)
    else:
        c6_basic, c6_alias, c6_m_alias, c6_first_m, c6_second_m, c6_m_scale, c6_second_phi, c6_native_alias = (False, False, False, False, False, False, False, False)
        c6_native, c6_generic, c6_scales = (False, False, False)
    J3MQ__guard(checks, check, 'C6_coordinate_radius_Pstar_and_L_source_rows', c6_basic)
    J3MQ__guard(checks, check, 'C6_first_switch_upper_actual_U_endpoint_alias', c6_alias)
    J3MQ__guard(checks, check, 'C6_first_switch_upper_moment_endpoint_alias', c6_m_alias)
    J3MQ__guard(checks, check, 'C6_first_switch_original_M_history_through_Z2', c6_first_m)
    J3MQ__guard(checks, check, 'C6_second_switch_M_history_through_Z2', c6_second_m)
    J3MQ__guard(checks, check, 'C6_original_m_fields_use_inverse_Pstar_once', c6_m_scale)
    J3MQ__guard(checks, check, 'C6_second_switch_actual_phi_transport_rows', c6_second_phi)
    J3MQ__guard(checks, check, 'C6_first_switch_upper_all_native_moment_endpoint_aliases', c6_native_alias)
    J3MQ__guard(checks, check, 'C6_native_six_history_transport_M_H_K_A_B_C_through_Z2', c6_native)
    J3MQ__guard(checks, check, 'C6_generic_h_k_e_p_rows_after_native_normalization', c6_generic)
    J3MQ__guard(checks, check, 'C6_common_W_FF_F0_Pstar_radius_source_rows', c6_scales)
    loop_rows = J3MQ__rows(loop)
    loop_prefix = all((J3MQ__row_is(loop_rows, node, **expected) for node, expected in ((73190, {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'}), (75949, {'operation': 'bound_variable', 'name': 'true_prescribed_second_switch_t'}), (85905, {'operation': 'product', 'arguments': [31709, 85904]}), (88600, {'operation': 'product', 'arguments': [31709, 88599]}), (33572, {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5}), (69245, {'operation': 'positive_quotient', 'numerator': 1, 'denominator': 33572, 'source_positive_certificate': 'same strictly positive Pstar'}))))
    J3MQ__guard(checks, check, 'LOOP_shared_C6_prefix_fixed_source_rows', loop_prefix)
    original_state_Z2_proved = all((checks.get(name, False) for name in ('decoded_background_and_loop_shapes', 'shared_authenticated_source_family', 'authenticated_C6_constructor_AST_and_calls', 'authenticated_bridge_histories_AST_and_calls', 'typed_C6_actual_moment_and_loop_chart_bindings', 'C6_coordinate_radius_Pstar_and_L_source_rows', 'C6_first_switch_upper_actual_U_endpoint_alias', 'C6_first_switch_upper_moment_endpoint_alias', 'C6_first_switch_original_M_history_through_Z2', 'C6_second_switch_M_history_through_Z2', 'C6_original_m_fields_use_inverse_Pstar_once', 'C6_second_switch_actual_phi_transport_rows', 'C6_first_switch_upper_all_native_moment_endpoint_aliases', 'C6_native_six_history_transport_M_H_K_A_B_C_through_Z2', 'C6_generic_h_k_e_p_rows_after_native_normalization', 'C6_common_W_FF_F0_Pstar_radius_source_rows', 'LOOP_shared_C6_prefix_fixed_source_rows')))
    return {'schema': J3MQ_SCHEMA, 'status': 'original_state_source_guard_ready' if original_state_Z2_proved else 'source_guard_failure', 'source_only': True, 'numeric_values_computed': False, 'accepted_identity_claim': False, 'original_state_Z2_proved': original_state_Z2_proved, 'original_five_moments_Z2_proved': bool(checks.get('typed_C6_actual_moment_and_loop_chart_bindings', False) and checks.get('C6_first_switch_upper_all_native_moment_endpoint_aliases', False) and checks.get('C6_native_six_history_transport_M_H_K_A_B_C_through_Z2', False) and checks.get('C6_generic_h_k_e_p_rows_after_native_normalization', False) and checks.get('C6_original_m_fields_use_inverse_Pstar_once', False) and checks.get('C6_common_W_FF_F0_Pstar_radius_source_rows', False)), 'ordinary_Z_orders': [0, 1, 2], 'owned_cell': {'lower': 0, 'upper': '1/2'}, 'endpoint_projections': J3MQ__ORIGINAL_ENDPOINT_PROJECTIONS, 'native6_certificate': J3MQ__NATIVE6_CERTIFICATE, 'generic5_certificate': J3MQ__GENERIC5_CERTIFICATE, 'source_pairs': {'physical_m': [[87591, 90109], [87592, 90110], [87593, 90111]], 'C6_actual_U': [[88437, 88437], [88438, 88438], [88439, 88439]], 'C6_moments_M': [[88465, 88465], [88466, 88466], [88467, 88467]], 'C6_native_moments_H': [[87352, 89870], [87353, 89871], [87354, 89872]], 'C6_native_moments_K': [[87427, 89945], [87428, 89946], [87429, 89947]], 'C6_native_moments_A': [[87462, 89980], [87463, 89981], [87464, 89982]], 'C6_native_moments_B': [[87497, 90015], [87498, 90016], [87499, 90017]], 'C6_native_moments_C': [[87525, 90043], [87526, 90044], [87527, 90045]]}, 'common_source_factors': {'Pstar': J3MQ_PSTAR, 'inverse_Pstar': J3MQ_PSTAR_INVERSE, 'L': J3MQ_COMMON_L}, 'proof_guards': checks, 'source_hashes': J3MQ__ORIGINAL_SOURCE_HASHES, 'remaining_premises': ['finite-N correction transport and MOD physical Q guards are supplied by prove_J3_physical_Q_transport', 'RU U C1 remains an external worker premise', 'root must perform final acceptance; this function claims no accepted identity'], 'nonclaims': {'physical_m_Z2_proved': False, 'J3_Q_C1_proved': False, 'Q_y_C1_identity': False, 'numeric_containment': False}}

def J3MQ__original_state_contract_ok(original_state):
    if not isinstance(original_state, dict):
        return False
    guards = original_state.get('proof_guards')
    return original_state.get('schema') == J3MQ_SCHEMA and original_state.get('status') == 'original_state_source_guard_ready' and (original_state.get('source_only') is True) and (original_state.get('accepted_identity_claim') is False) and (original_state.get('original_state_Z2_proved') is True) and (original_state.get('original_five_moments_Z2_proved') is True) and (original_state.get('ordinary_Z_orders') == [0, 1, 2]) and (original_state.get('owned_cell') == {'lower': 0, 'upper': '1/2'}) and (original_state.get('endpoint_projections') == J3MQ__ORIGINAL_ENDPOINT_PROJECTIONS) and (original_state.get('source_hashes') == J3MQ__ORIGINAL_SOURCE_HASHES) and (original_state.get('native6_certificate') == J3MQ__NATIVE6_CERTIFICATE) and (original_state.get('generic5_certificate') == J3MQ__GENERIC5_CERTIFICATE) and isinstance(guards, dict) and (set(guards) == set(J3MQ__ORIGINAL_GUARD_NAMES)) and all((guards.get(name) is True for name in J3MQ__ORIGINAL_GUARD_NAMES))

def J3MQ__ru_result_contract_ok(ru_result):
    if not isinstance(ru_result, dict):
        return False
    predicates = ru_result.get('predicates')
    return ru_result.get('schema') == 'J3_R_U_source_transport_guarded_v3' and ru_result.get('status') == 'conditional_source_theorem_proved' and (ru_result.get('source_only') is True) and (ru_result.get('accepted_identity_claim', False) is False) and (ru_result.get('J3_R_C1_proved') is True) and (ru_result.get('J3_U_C1_proved') is True) and (ru_result.get('ordinary_Z_orders') == [0, 1]) and (ru_result.get('owned_Z_cell') == [0, 0.5]) and (ru_result.get('exact_U_row_pairs') == J3MQ__RU_ENDPOINT_PAIRS) and (ru_result.get('source_hashes') == J3MQ__RU_SOURCE_HASHES) and (ru_result.get('numeric_values_computed') is False) and (ru_result.get('source_constructor_calls') == 0) and (ru_result.get('old_proof_functions_called') == 0) and (ru_result.get('source_producer_function_calls') == 0) and (ru_result.get('generic_DAG_walks') == 0) and isinstance(predicates, dict) and bool(predicates) and all((value is True for value in predicates.values()))

def prove_J3_physical_Q_transport(background, loop, mod, loop_creator_text, mod_creator_text, original_state, ru_result, check):
    """Guard finite-N physical m through Z2 and Q through ordinary Z1.

    This function requires the independently guarded original state and RU/U
    premise.  It checks the actual LOOP correction rows, one Pstar physical
    scale, actual MOD Q numerator/quotient rows, and the common positive L
    source row.  No S1 alias or registry self-comparison is used as a proof.
    """
    checks = {}
    shapes = J3MQ__archive_shapes(background, loop, mod)
    J3MQ__guard(checks, check, 'decoded_background_loop_and_mod_shapes', shapes)
    source_family = shapes and J3MQ__shared_source_family(background, loop, mod)
    J3MQ__guard(checks, check, 'shared_authenticated_source_family', source_family)
    original_ok = J3MQ__original_state_contract_ok(original_state)
    J3MQ__guard(checks, check, 'accepted_original_state_Z2_premise', original_ok)
    ru_ok = J3MQ__ru_result_contract_ok(ru_result)
    J3MQ__guard(checks, check, 'accepted_external_RU_U_C1_premise', ru_ok)
    loop_ast = J3MQ__loop_source_ast_ok(loop_creator_text)
    mod_ast = J3MQ__mod_source_ast_ok(mod_creator_text)
    J3MQ__guard(checks, check, 'authenticated_loop_finite_N_transport_AST_and_calls', loop_ast)
    J3MQ__guard(checks, check, 'authenticated_MOD_physical_pressure_AST_and_calls', mod_ast)
    typed_loop, endpoint_aliases = J3MQ__loop_typed_bindings(loop) if shapes else (False, {})
    J3MQ__guard(checks, check, 'typed_loop_first_second_switch_bindings', typed_loop)
    finite_rows, first_history_rows, second_history_rows = J3MQ__finite_N_rows(loop, endpoint_aliases) if shapes else (False, {}, {})
    J3MQ__guard(checks, check, 'finite_N_m_h_k_e_p_transport_rows_through_Z2', finite_rows)
    physical_bindings = J3MQ__physical_bindings(mod) if shapes else False
    J3MQ__guard(checks, check, 'typed_MOD_physical_m_Q_endpoint_bindings', physical_bindings)
    physical_u = J3MQ__physical_u_rows(mod) if shapes else False
    J3MQ__guard(checks, check, 'physical_U_Pstar_source_rows', physical_u)
    physical_m = J3MQ__physical_m_rows(mod) if shapes else False
    J3MQ__guard(checks, check, 'physical_m_Pstar_applied_once_through_Z2', physical_m)
    q_rows = J3MQ__q_rows(mod) if shapes else False
    J3MQ__guard(checks, check, 'physical_Q_numerator_and_common_L_source_rows', q_rows)
    finite_N_m_Z2_proved = all((checks.get(name, False) for name in ('decoded_background_loop_and_mod_shapes', 'shared_authenticated_source_family', 'accepted_original_state_Z2_premise', 'authenticated_loop_finite_N_transport_AST_and_calls', 'typed_loop_first_second_switch_bindings', 'finite_N_m_h_k_e_p_transport_rows_through_Z2')))
    physical_m_Z2_proved = all((checks.get(name, False) for name in ('decoded_background_loop_and_mod_shapes', 'shared_authenticated_source_family', 'accepted_original_state_Z2_premise', 'authenticated_loop_finite_N_transport_AST_and_calls', 'typed_loop_first_second_switch_bindings', 'finite_N_m_h_k_e_p_transport_rows_through_Z2', 'typed_MOD_physical_m_Q_endpoint_bindings', 'physical_m_Pstar_applied_once_through_Z2')))
    J3_Q_C1_proved = all((physical_m_Z2_proved, checks.get('accepted_external_RU_U_C1_premise', False), checks.get('authenticated_MOD_physical_pressure_AST_and_calls', False), checks.get('physical_U_Pstar_source_rows', False), checks.get('physical_Q_numerator_and_common_L_source_rows', False)))
    return {'schema': J3MQ_SCHEMA, 'status': 'J3_Q_source_guard_ready' if J3_Q_C1_proved else 'source_guard_failure', 'source_only': True, 'numeric_values_computed': False, 'accepted_identity_claim': False, 'finite_N_m_Z2_proved': finite_N_m_Z2_proved, 'physical_m_Z2_proved': physical_m_Z2_proved, 'J3_Q_C1_proved': J3_Q_C1_proved, 'ordinary_Z_orders': [0, 1], 'm_ordinary_Z_orders': [0, 1, 2], 'owned_cell': {'lower': 0, 'upper': '1/2'}, 'Q_C1_conditional_on': ['original_state_Z2_proved', 'J3_U_C1_proved'], 'endpoint_projections': J3MQ__ORIGINAL_ENDPOINT_PROJECTIONS, 'source_pairs': {'physical_m_through_Z2': [[148097, 154259], [148098, 154260], [148099, 154261]], 'physical_U_premise': [[148085, 154247], [148086, 154248]], 'physical_Q': [[148180, 154342], [148184, 154346]], 'physical_Q_numerators': [[148175, 154337], [148183, 154345]], 'loop_modified_m': [[107095, 109504], [107096, 109505], [107097, 109506]]}, 'common_source_factors': {'Pstar': J3MQ_PSTAR, 'inverse_Pstar': J3MQ_PSTAR_INVERSE, 'L': J3MQ_COMMON_L, 'L_recipe': {'operation': 'sum', 'arguments': [1, 25271]}}, 'history_rows': {'first_switch': first_history_rows, 'second_switch': second_history_rows}, 'proof_guards': checks, 'source_hashes': {'c6_raw_archive_sha256': J3MQ_C6_RAW_ARCHIVE_SHA256, 'loop_raw_archive_sha256': J3MQ_LOOP_RAW_ARCHIVE_SHA256, 'mod_raw_archive_sha256': J3MQ_MOD_RAW_ARCHIVE_SHA256, 'loop_creator_sha256': J3MQ_LOOP_CREATOR_SHA256, 'mod_creator_sha256': J3MQ_MOD_CREATOR_SHA256}, 'nonclaims': {'Q_y_C1_identity': False, 'weighted_J3_jump_zero': False, 'numeric_containment': False, 'full_C5_physical_install': False, 'accepted_identity_claim': False}}


# Integrated fresh source fragment: native_J3_R_U_source_fragment_v4_2026_10_21.py
"""Guarded, unexecuted source proof fragment for the J3 R/U seam.

This is the fourth version of the J3 R/U fragment.  The earlier static
outline and guarded v2/v3 are intentionally left untouched.  This file contains
executable guards, but it was not imported or called while preparing the
artifact.

The guarded entry point accepts decoded LOOP, MOD, and BACKGROUND C6
archives, plus the creator source texts.  It checks exact wrapper names,
archive lengths, selected fixed graph rows, chart binding rows, source AST
body hashes, the six native C6 histories through Z2, and the MQ
``original_state_Z2_proved`` plus ``original_five_moments_Z2_proved`` premise.
It invokes the parent-supplied ``check`` callback only after those actual
predicates pass.  The two proof flags are emitted as true only after that
callback passes.

The guards are source-only.  They do not evaluate a symbolic graph, execute a
source constructor, call an operator/provider/normalizer, perform arithmetic,
run quadrature, or inspect the forbidden continuation or density archive.
"""
import ast as J3RU_ast
import hashlib as J3RU_hashlib
from collections.abc import Mapping as J3RU_Mapping
J3RU_FRAGMENT_ID = 'native_J3_R_U_source_fragment_v4_2026_10_21'
J3RU_SCHEMA = 'J3_R_U_source_transport_guarded_v3'
J3RU_STATUS_UNEXECUTED = 'unexecuted_guarded_source_fragment'
J3RU_RUNTIME_MODEL_METADATA = None
J3RU_EXPECTED_SOURCE_SHA256 = {'loop_creator_text': '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d', 'mod_creator_text': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'c6_creator_text': '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32', 'bridge_creator_text': '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605'}
J3RU_EXPECTED_RAW_ARCHIVE_SHA256 = {'loop': '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa', 'mod': 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a', 'background': '5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a', 'parent_bridge': '5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92'}
J3RU_EXPECTED_AST_BODY_SHA256 = {'c6_creator_text': {('function', 'first_switch_fields'): '0bd697b785e7241b2ea8520d74267342e4993f23dffc90fe49617d13430d68e2', ('function', 'second_switch_fields'): 'c6c08d15f7b5058f9ff48358ea904eed294d21046ca016deb7d1c440abd4d438', ('method', 'CurrentOriginalPrescribedR100RmC6Functions', '__init__'): 'd7e974bcc4e2156f8c05b3a67192aefe1df25cb0857df1e7a47e8bd1f6eb2b32', ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'add_chart'): '18f1dbdd7128ac2e72098ab5d7249659923a5d7166640057b2587acfcf2f3812', ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'end'): '3a9ca663af71f7dcf34fa945662ea1535de41a0d8bc2a53a2aacf2954d5b6040'}, 'loop_creator_text': {('function', 'loop_functions'): 'fec28c8c3ff3bf20a5cc64bc777b0b0ed304387b95f71cb6a8f9ec30fd9e0559', ('function', 'generic_roots'): '3ba21f81b7aae33751e6534e15e97311c1b851f103b522259a80d1170fad7847', ('function', 'finite_N_functions'): '219500c1ac0665e869303a421303a9d2a4c05bfa3121a9fe1f06f022649ed4d2', ('function', 'transport'): 'fd5f45c38d246f2555c045ca2d386f5b49bbfaa7699981156efbf99c1869e9b3', ('method', 'CurrentOriginalPrescribedLoopC5Functions', '__init__'): '4e6a13cf1025d85e3bf0d7c5abfe9b90d89e9cf45692ee6f5b6161b4ed8546c2'}, 'mod_creator_text': {('function', 'physical_pressure_fields'): '619f51980979fc35ad965386a0038fc2785bf15c99ef2a1b61d749a2d6833a14', ('method', 'CurrentOriginalModifiedQ4PressureC3Functions', '__init__'): '8040fa64b5da2de7fa1e2e29cd2583bb803b5b0735d36357050200ad87d645d4'}, 'bridge_creator_text': {('function', 'histories'): '9dc010e43b2a9b1c498380b4452f226f467df7f5529f4c026ed5b8f7a6aaaf5d', ('function', 'generic_fields'): '41f3b37e36843799b913ec6af1da1d59ef8dcd668c34a0b3584621fe7e791cdf', ('method', 'CurrentOriginalPrescribedBridgeC6Functions', '__init__'): 'df1b4b43d63a07bc8e44aa8c5db8bed8c09d73ac31626b3befbd69e750f93e5b'}}
J3RU_EXPECTED_C6_CHARTS = {'first_switch': {'phase': 73190, 'domain': [0, 1], 'radius': 85905, 'physical_measure': 9, 'actual': {'phi': [86947, 86950, 86954, 86959, 86965, 86972, 86980], 'U': [87154, 87155, 87156, 87157, 87158, 87159, 87160], 'phi_y': [87196, 87197, 87198, 87199, 87200, 87201, 87202], 'U_y': [87126, 87127, 87128, 87129, 87130, 87131, 87132]}, 'fields': {'E': [87550, 87553, 87557, 87562, 87568, 87575, 87583], 'V': [87584, 87585, 87586, 87587, 87588, 87589, 87590]}, 'radial': {'E_y': [87790, 87793, 87797, 87802, 87808, 87815, 87823], 'V_y': [87824, 87825, 87826, 87827, 87828, 87829, 87830]}, 'generic': {'p2': [88423, 88424, 88425, 88426, 88427, 88428]}, 'moments': {'M': [87392, 87393, 87394, 87395, 87396, 87397, 87398]}}, 'second_switch': {'phase': 75949, 'domain': [0, 1], 'radius': 88600, 'physical_measure': 9, 'actual': {'phi': [89644, 89647, 89651, 89656, 89662, 89669, 89677], 'U': [88437, 88438, 88439, 88440, 88441, 88442, 88443], 'phi_y': [89687, 89690, 89694, 89699, 89705, 89712, 89720], 'U_y': [0, 0, 0, 0, 0, 0, 0]}, 'fields': {'E': [90068, 90071, 90075, 90080, 90086, 90093, 90101], 'V': [90102, 90103, 90104, 90105, 90106, 90107, 90108]}, 'radial': {'E_y': [90308, 90311, 90315, 90320, 90326, 90333, 90341], 'V_y': [0, 0, 0, 0, 0, 0, 0]}, 'generic': {'p2': [90934, 90935, 90936, 90937, 90938, 90939]}, 'moments': {'M': [89910, 89911, 89912, 89913, 89914, 89915, 89916]}}}
J3RU_EXPECTED_C6_MOMENT_ROWS = {'first_switch': {'H': [87352, 87353, 87354, 87355, 87356, 87357, 87358], 'M': [87392, 87393, 87394, 87395, 87396, 87397, 87398], 'K': [87427, 87428, 87429, 87430, 87431, 87432, 87433], 'A': [87462, 87463, 87464, 87465, 87466, 87467, 87468], 'B': [87497, 87498, 87499, 87500, 87501, 87502, 87503], 'C': [87525, 87526, 87527, 87528, 87529, 87530, 87531]}, 'second_switch': {'H': [89870, 89871, 89872, 89873, 89874, 89875, 89876], 'M': [89910, 89911, 89912, 89913, 89914, 89915, 89916], 'K': [89945, 89946, 89947, 89948, 89949, 89950, 89951], 'A': [89980, 89981, 89982, 89983, 89984, 89985, 89986], 'B': [90015, 90016, 90017, 90018, 90019, 90020, 90021], 'C': [90043, 90044, 90045, 90046, 90047, 90048, 90049]}}
J3RU_EXPECTED_C6_R1_MOMENT_ROWS = {'H': [88458, 88459, 88460, 88461, 88462, 88463, 88464], 'M': [88465, 88466, 88467, 88468, 88469, 88470, 88471], 'K': [88472, 88473, 88474, 88475, 88476, 88477, 88478], 'A': [88479, 88480, 88481, 88482, 88483, 88484, 88485], 'B': [88486, 88487, 88488, 88489, 88490, 88491, 88492], 'C': [88493, 88494, 88495, 88496, 88497, 88498, 88499]}
J3RU_EXPECTED_C6_KNOWN_COMPARISON = {'first_switch': {'D_over_R': [86464, 86468, 86474, 86480, 86486, 86492, 86498], 'drive_hydro': [86650, 86654, 86660, 86666, 86672, 86678, 86684], 'drive_pressure': [68391, 68395, 68401, 68407, 68413, 68419, 68425], 'drive_swirl': [86853, 86857, 86863, 86870, 86878, 86887, 86897]}, 'second_switch': {'D_over_R': [89159, 89163, 89169, 89175, 89181, 89187, 89193], 'drive_hydro': [89345, 89349, 89355, 89361, 89367, 89373, 89379], 'drive_pressure': [68391, 68395, 68401, 68407, 68413, 68419, 68425], 'drive_swirl': [89548, 89552, 89558, 89565, 89573, 89582, 89592]}}
J3RU_EXPECTED_C6_R1_FIELDS = {'E': [88500, 88501, 88502, 88503, 88504, 88505, 88506], 'V': [88507, 88508, 88509, 88510, 88511, 88512, 88513]}
J3RU_EXPECTED_C6_R1_GENERIC = {'p2': [88592, 88593, 88594, 88595, 88596, 88597]}
J3RU_EXPECTED_C6_HISTORY_Z2_ROWS = (('H', 'true_first_switch', (87345, 87346, 87347), 2, 'true_first_switch_history_t', 73190), ('M', 'true_first_switch', (87385, 87386, 87387), 1, 'true_first_switch_history_t', 73190), ('K', 'true_first_switch', (87420, 87421, 87422), 2, 'true_first_switch_history_t', 73190), ('A', 'true_first_switch', (87455, 87456, 87457), 1, 'true_first_switch_history_t', 73190), ('B', 'true_first_switch', (87490, 87491, 87492), 2, 'true_first_switch_history_t', 73190), ('C', 'true_first_switch', (87518, 87519, 87520), 1, 'true_first_switch_history_t', 73190), ('H', 'true_second_switch', (89863, 89864, 89865), 2, 'true_second_switch_history_t', 75949), ('M', 'true_second_switch', (89903, 89904, 89905), 1, 'true_second_switch_history_t', 75949), ('K', 'true_second_switch', (89938, 89939, 89940), 2, 'true_second_switch_history_t', 75949), ('A', 'true_second_switch', (89973, 89974, 89975), 1, 'true_second_switch_history_t', 75949), ('B', 'true_second_switch', (90008, 90009, 90010), 2, 'true_second_switch_history_t', 75949), ('C', 'true_second_switch', (90036, 90037, 90038), 1, 'true_second_switch_history_t', 75949))
J3RU_EXPECTED_LOOP_CHARTS = {'first_switch': {'coordinate': 73190, 'domain': [0, 1], 'radius': 85905, 'physical_measure': 9, 'global_offset': 104535, 'actual_phase': 104537, 'roots': {'E': [87550, 87553, 87557, 87562, 87568, 87575], 'V': [87584, 87585, 87586, 87587, 87588, 87589], 'a': [104562, 104566, 104572, 104579, 104587, 104596], 'b': [104597, 104601, 104607, 104614, 104622, 104631], 'p1': [104650, 104654, 104660, 104667, 104675, 104684], 'p2': [88423, 88424, 88425, 88426, 88427, 88428]}, 'loop': {'B_over_Pstar': [106437, 106438, 106439, 106440, 106441, 106442], 'A': [106399, 106400, 106401, 106402, 106403, 106404], 'q': [105094, 105095, 105096, 105097, 105098, 105099], 'inverse': [106017, 106040], 'flat_predicate': 104794, 'phase_partial_coefficients': [{'Z': 0, 'angle': 0, 'node': 105891}, {'Z': 1, 'angle': 0, 'node': 105894}]}, 'finite_N': {'V_N': [106932, 106933, 106934, 106935, 106936, 106937], 'delta_V': [106691, 106692, 106693, 106694, 106695, 106696]}}, 'second_switch': {'coordinate': 75949, 'domain': [0, 1], 'radius': 88600, 'physical_measure': 9, 'global_offset': 107163, 'actual_phase': 107165, 'roots': {'E': [90068, 90071, 90075, 90080, 90086, 90093], 'V': [90102, 90103, 90104, 90105, 90106, 90107], 'a': [107184, 107188, 107194, 107201, 107209, 107218], 'b': [0, 0, 0, 0, 0, 0], 'p1': [107237, 107241, 107247, 107254, 107262, 107271], 'p2': [90934, 90935, 90936, 90937, 90938, 90939]}, 'loop': {'B_over_Pstar': [108857, 108858, 108859, 108860, 108861, 108862], 'A': [108819, 108820, 108821, 108822, 108823, 108824], 'q': [107573, 107574, 107575, 107576, 107577, 107578], 'inverse': [108455, 108478], 'flat_predicate': 107273, 'phase_partial_coefficients': [{'Z': 0, 'angle': 0, 'node': 108329}, {'Z': 1, 'angle': 0, 'node': 108332}]}, 'finite_N': {'V_N': [109340, 109341, 109342, 109343, 109344, 109345], 'delta_V': [109099, 109100, 109101, 109102, 109103, 109104]}}}
J3RU_EXPECTED_MOD_PHYSICAL_U = {'first_switch': [148085, 148086, 148087, 148088, 148089, 148090], 'second_switch': [154247, 154248, 154249, 154250, 154251, 154252]}
J3RU__MISSING = object()

def J3RU__path(value, *keys):
    current = value
    for key in keys:
        if not isinstance(current, J3RU_Mapping) or key not in current:
            return J3RU__MISSING
        current = current[key]
    return current

def J3RU__row_has(nodes, index, expected):
    if not isinstance(nodes, list) or not isinstance(index, int):
        return False
    if index < 0 or index >= len(nodes) or (not isinstance(nodes[index], J3RU_Mapping)):
        return False
    row = nodes[index]
    return all((key in row and row[key] == value for key, value in expected.items()))

def J3RU__equal_path(value, keys, expected):
    return J3RU__path(value, *keys) == expected

def J3RU__archive_shape(archive, wrapper, node_count):
    return isinstance(archive, J3RU_Mapping) and isinstance(archive.get('exact_graph_nodes'), list) and (len(archive['exact_graph_nodes']) == node_count) and isinstance(archive.get(wrapper), J3RU_Mapping)

def J3RU__source_body_hash(node):
    body = J3RU_ast.Module(body=list(node.body), type_ignores=[])
    dumped = J3RU_ast.dump(body, annotate_fields=True, include_attributes=False)
    return J3RU_hashlib.sha256(dumped.encode('utf-8')).hexdigest()

def J3RU__definition(tree, spec):
    if spec[0] == 'function':
        for node in tree.body:
            if isinstance(node, (J3RU_ast.FunctionDef, J3RU_ast.AsyncFunctionDef)) and node.name == spec[1]:
                return node
        return None
    if spec[0] == 'method':
        for node in tree.body:
            if isinstance(node, J3RU_ast.ClassDef) and node.name == spec[1]:
                for child in node.body:
                    if isinstance(child, (J3RU_ast.FunctionDef, J3RU_ast.AsyncFunctionDef)) and child.name == spec[2]:
                        return child
        return None
    return None

def J3RU__creator_ast_guards(texts):
    predicates = {}
    details = {}
    needles = {'c6_creator_text': {('function', 'first_switch_fields'): ('Uy=a.scale(a.mul(ratio,drive),g.neg(chi_axial))', 'chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))', "U=a.add(initial['U'],a.integral(", 'source_measure_applied_once=True', 'return dict(phi=phi,U=U,phi_y=phiy,U_y=Uy'), ('function', 'second_switch_fields'): ('sig=sigma(g,phase)', 'original_second_switch_full_current_radius_direction_and_power=True', "return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)"), ('method', 'CurrentOriginalPrescribedR100RmC6Functions', '__init__'): ('known,barphi=self.comparison_at_radius(R)', "initial={key:bridge_owner.R100['actual'][key] for key in ('phi','U')}", "self.add_chart('first_switch',s,R,hb,actual,bridge_owner.R100['moments'],known=known)", "self.endpoints['R1']=self.end('first_switch',g.one)", "incoming=self.endpoints['R1']", "actual=second_switch_fields(a,known,incoming['actual'],t,hb,R)", "self.add_chart('second_switch',t,R,hb,actual,incoming['moments'],known=known)"), ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'add_chart'): ("J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,", 'fields,radial,generic=previous.generic_fields('), ('method', 'CurrentOriginalPrescribedR100RmC6Functions', 'end'): ('Z_independent_substitution=True', "actual={key:at(chart['actual'][key])", "fields={k:at(jet) for k,jet in chart['fields'].items()}", "generic={k:Algebra(self.g,5).at(jet,chart['phase'],value) for k,jet in chart['generic'].items()}")}, 'bridge_creator_text': {('function', 'histories'): ('def histories(a,phi,U,initial,phase,jacobian,label):', 'bodies=densities(a,phi,U);out={}', 'for name,rate in RATES.items():', 'actual_six_history_rate=rate', 'original_incoming_memory_retained=True', 'source_measure_applied_once=True', 'prescribed_bridge_chart=label'), ('function', 'generic_fields'): ("fields=dict(E=W*actual['phi'],V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse),", "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))")}, 'loop_creator_text': {('function', 'loop_functions'): ("shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)", "g.quotient(g.one,eta,'same positive canonical eta')", "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))", "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))", "inverse=g.node('original_monotone_phase_inverse'", 'return dict(q=q,inverse=psi,A=A,B_over_Pstar=B', 'angle_endpoints_independent_of_Z=True'), ('function', 'generic_roots'): ("E,V=(a.trim(field[key]) for key in ('E','V'))", "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)", 'aa,bb=divide(C),divide(B)', "p1=divide(a.scale(a.add(generic['theta_linear'],a.scale(generic['theta_quadratic'],Pstar)),chart['radius']))", "p2=generic['p2']"), ('function', 'finite_N_functions'): ("s=a.scale(loop['B_over_Pstar'],invN)", 'return dict(E_N=a.mul(E,ex),V_N=a.add(V,s),'), ('method', 'CurrentOriginalPrescribedLoopC5Functions', '__init__'): ("self.eta=self.unique('current_original_source_parameter','eta')", "self.dstar=self.unique('current_original_source_parameter','d_star')", "phase=g.unary('fractional_part',g.mul(self.N,offset))", 'roots=generic_roots(g,chart,self.Pstar)', 'loop=loop_functions(g,roots,phase,self.eta,self.dstar,name,coordinate,self.sc)', 'finite=finite_N_functions(g,roots,loop,self.N)', "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)", "transport(g,finite['five_signed_density_increments'],incoming,", "fields={key:a.trim(jet) for key,jet in chart['fields'].items()}", "fields.update(E=finite['E_N'],V=finite['V_N'])", 'self.endpoints[name]={key:a.at(jet,coordinate,upper) for key,jet in history.items()}'), ('function', 'transport'): ('actual_finite_N_signed_correction=True', 'physical_Jacobian_applied_once=True', 'predecessor_memory_never_reset=True')}, 'mod_creator_text': {('function', 'physical_pressure_fields'): ("U=a.scale(fields['V'],Pstar)", "Uy=a.scale(radial['V_y'],Pstar)"), ('method', 'CurrentOriginalModifiedQ4PressureC3Functions', '__init__'): ("chart=loop_owner.charts[name];bg=chart['source_chart'];coordinate=chart['coordinate'];jac=chart['physical_measure']", "physical=physical_pressure_fields(g,chart['modified_fields'],radial,self.background.delta,self.background.Pstar)", 'self.charts[name]=dict(source_chart=chart')}}
    for key, text in texts.items():
        expected_sha = J3RU_EXPECTED_SOURCE_SHA256[key]
        source_ok = isinstance(text, str) and J3RU_hashlib.sha256(text.encode('utf-8')).hexdigest() == expected_sha
        predicates[key + '.source_sha256'] = source_ok
        if not source_ok:
            details[key] = {'source_sha256': 'mismatch'}
            continue
        try:
            tree = J3RU_ast.parse(text, filename=key)
        except (SyntaxError, TypeError, ValueError) as exc:
            predicates[key + '.ast_parse'] = False
            details[key] = {'ast_parse': type(exc).__name__}
            continue
        predicates[key + '.ast_parse'] = True
        for spec, expected_hash in J3RU_EXPECTED_AST_BODY_SHA256[key].items():
            node = J3RU__definition(tree, spec)
            spec_name = '.'.join(spec)
            present = node is not None
            predicates[key + '.body_present.' + spec_name] = present
            if not present:
                details.setdefault(key, {}).setdefault('missing_definitions', []).append(spec_name)
                continue
            digest = J3RU__source_body_hash(node)
            predicates[key + '.body_sha256.' + spec_name] = digest == expected_hash
            segment = J3RU_ast.get_source_segment(text, node) or ''
            expected_needles = needles.get(key, {}).get(spec, ())
            predicates[key + '.body_needles.' + spec_name] = all((needle in segment for needle in expected_needles))
            if digest != expected_hash or not all((needle in segment for needle in expected_needles)):
                details.setdefault(key, {}).setdefault('body_failures', []).append(spec_name)
    return (predicates, details)

def J3RU__archive_guards(loop, mod, background, archive_raw_sha256, original_state=None):
    predicates = {}
    details = {}
    loop_nodes = loop.get('exact_graph_nodes') if isinstance(loop, J3RU_Mapping) else None
    mod_nodes = mod.get('exact_graph_nodes') if isinstance(mod, J3RU_Mapping) else None
    c6_nodes = background.get('exact_graph_nodes') if isinstance(background, J3RU_Mapping) else None
    predicates['archive.loop_shape'] = J3RU__archive_shape(loop, 'prescribed_loop_C5_bindings', 121834)
    predicates['archive.mod_shape'] = J3RU__archive_shape(mod, 'modified_Q4_pressure_C3_bindings', 183700)
    predicates['archive.background_shape'] = J3RU__archive_shape(background, 'prescribed_R100_Rm_C6_bindings', 97295)
    predicates['archive.raw_sha256'] = isinstance(archive_raw_sha256, J3RU_Mapping) and all((isinstance(archive_raw_sha256.get(name), str) and len(archive_raw_sha256.get(name, '')) == 64 for name in ('loop', 'mod', 'background'))) and (archive_raw_sha256.get('loop') == J3RU_EXPECTED_RAW_ARCHIVE_SHA256['loop']) and (archive_raw_sha256.get('mod') == J3RU_EXPECTED_RAW_ARCHIVE_SHA256['mod']) and (archive_raw_sha256.get('background') == J3RU_EXPECTED_RAW_ARCHIVE_SHA256['background'])
    if not all(predicates.values()):
        return (predicates, details)
    c6 = background['prescribed_R100_Rm_C6_bindings']
    lp = loop['prescribed_loop_C5_bindings']
    mp = mod['modified_Q4_pressure_C3_bindings']
    source_family = background.get('source_family')
    predicates['source_family.mapping_shape'] = isinstance(source_family, J3RU_Mapping) and all((isinstance(source_family.get(name), str) and len(source_family.get(name, '')) == 64 for name in ('actual_five_defect_family_sha256', 'implicit_source_sha256', 'datum_enclosure_sha256')))
    predicates['source_family.join'] = predicates['source_family.mapping_shape'] and source_family == loop.get('source_family') == mod.get('source_family')
    predicates['bindings.chart_maps'] = isinstance(c6.get('charts'), J3RU_Mapping) and isinstance(c6.get('endpoints'), J3RU_Mapping) and isinstance(lp.get('charts'), J3RU_Mapping) and isinstance(mp.get('charts'), J3RU_Mapping)
    if not predicates['bindings.chart_maps']:
        return (predicates, details)
    c6_rows = {9: {'operation': 'current_original_source_parameter', 'name': 'hbS', 'Z_independent': True}, 31709: {'operation': 'exact_rational', 'numerator': 100, 'denominator': 1}, 73190: {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'}, 75949: {'operation': 'bound_variable', 'name': 'true_prescribed_second_switch_t'}, 75950: {'operation': 'sum', 'arguments': [1, 75949]}, 85903: {'operation': 'product', 'arguments': [9, 73190]}, 85904: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 85903}, 85905: {'operation': 'product', 'arguments': [31709, 85904]}, 88598: {'operation': 'product', 'arguments': [9, 75950]}, 88599: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 88598}, 88600: {'operation': 'product', 'arguments': [31709, 88599]}, 88430: {'operation': 'function_substitution', 'expression': 86947, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88431: {'operation': 'function_substitution', 'expression': 86950, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88432: {'operation': 'function_substitution', 'expression': 86954, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88437: {'operation': 'function_substitution', 'expression': 87154, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88438: {'operation': 'function_substitution', 'expression': 87155, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88444: {'operation': 'function_substitution', 'expression': 87196, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88445: {'operation': 'function_substitution', 'expression': 87197, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 88446: {'operation': 'function_substitution', 'expression': 87198, 'variable': 73190, 'value': 1, 'Z_independent_substitution': True}, 89644: {'operation': 'product', 'arguments': [88430, 89617]}, 89645: {'operation': 'product', 'arguments': [88430, 89618]}, 89646: {'operation': 'product', 'arguments': [88431, 89617]}, 89647: {'operation': 'sum', 'arguments': [89645, 89646]}, 89648: {'operation': 'product', 'arguments': [88430, 89621]}, 89649: {'operation': 'product', 'arguments': [48, 88431, 89618]}, 89650: {'operation': 'product', 'arguments': [88432, 89617]}, 89651: {'operation': 'sum', 'arguments': [89648, 89649, 89650]}, 89687: {'operation': 'product', 'arguments': [89644, 89686]}, 89690: {'operation': 'sum', 'arguments': [89688, 89689]}, 89694: {'operation': 'sum', 'arguments': [89691, 89692, 89693]}}
    predicates['graph.c6_fixed_rows'] = all((J3RU__row_has(c6_nodes, i, expected) for i, expected in c6_rows.items()))
    predicates['bindings.c6_Pstar'] = c6.get('Pstar') == 33572
    for chart_name, expected in J3RU_EXPECTED_C6_CHARTS.items():
        chart = c6.get('charts', {}).get(chart_name)
        predicates['bindings.c6_chart.' + chart_name] = isinstance(chart, J3RU_Mapping) and chart.get('phase') == expected['phase'] and (chart.get('domain') == expected['domain']) and (chart.get('radius') == expected['radius']) and (chart.get('physical_measure') == expected['physical_measure']) and all((chart.get('actual', {}).get(key) == value for key, value in expected['actual'].items())) and (chart.get('fields', {}).get('E') == expected['fields']['E']) and (chart.get('fields', {}).get('V') == expected['fields']['V']) and (chart.get('radial', {}).get('E_y') == expected['radial']['E_y']) and (chart.get('radial', {}).get('V_y') == expected['radial']['V_y']) and (chart.get('generic', {}).get('p2') == expected['generic']['p2']) and all((chart.get('moments', {}).get(channel) == rows for channel, rows in J3RU_EXPECTED_C6_MOMENT_ROWS[chart_name].items())) and (chart.get('known_comparison_direction') == J3RU_EXPECTED_C6_KNOWN_COMPARISON[chart_name])
    predicates['bindings.c6_U_y_cutoff_flat_law'] = c6.get('charts', {}).get('first_switch', {}).get('actual', {}).get('U_y') == [87126, 87127, 87128, 87129, 87130, 87131, 87132] and c6.get('charts', {}).get('second_switch', {}).get('actual', {}).get('U_y') == [0, 0, 0, 0, 0, 0, 0] and (c6.get('charts', {}).get('second_switch', {}).get('radial', {}).get('V_y') == [0, 0, 0, 0, 0, 0, 0])
    r1 = c6.get('endpoints', {}).get('R1')
    predicates['bindings.c6_R1_endpoint'] = isinstance(r1, J3RU_Mapping) and r1.get('actual', {}).get('phi') == [88430, 88431, 88432, 88433, 88434, 88435, 88436] and (r1.get('actual', {}).get('U') == J3RU_EXPECTED_C6_CHARTS['second_switch']['actual']['U']) and (r1.get('actual', {}).get('phi_y') == [88444, 88445, 88446, 88447, 88448, 88449, 88450]) and (r1.get('actual', {}).get('U_y') == [88451, 88452, 88453, 88454, 88455, 88456, 88457]) and (r1.get('fields', {}).get('E') == J3RU_EXPECTED_C6_R1_FIELDS['E']) and (r1.get('fields', {}).get('V') == J3RU_EXPECTED_C6_R1_FIELDS['V']) and (r1.get('generic', {}).get('p2') == J3RU_EXPECTED_C6_R1_GENERIC['p2']) and all((r1.get('moments', {}).get(channel) == rows for channel, rows in J3RU_EXPECTED_C6_R1_MOMENT_ROWS.items())) and (c6.get('charts', {}).get('second_switch', {}).get('incoming_histories') == r1.get('moments'))
    predicates['bindings.c6_known_comparison_and_R1_source'] = all((c6.get('charts', {}).get(name, {}).get('known_comparison_direction') == J3RU_EXPECTED_C6_KNOWN_COMPARISON[name] for name in ('first_switch', 'second_switch'))) and predicates['bindings.c6_R1_endpoint']
    predicates['seam.c6_first_to_second_U'] = isinstance(r1, J3RU_Mapping) and r1.get('actual', {}).get('U') == c6.get('charts', {}).get('second_switch', {}).get('actual', {}).get('U') and (r1.get('actual', {}).get('U', [])[:2] == [88437, 88438])
    loop_chart_checks = []
    for chart_name, expected in J3RU_EXPECTED_LOOP_CHARTS.items():
        chart = lp.get('charts', {}).get(chart_name)
        if not isinstance(chart, J3RU_Mapping):
            loop_chart_checks.append(False)
            continue
        loop_chart_checks.append(chart.get('coordinate') == expected['coordinate'] and chart.get('domain') == expected['domain'] and (chart.get('radius') == expected['radius']) and (chart.get('physical_measure') == expected['physical_measure']) and (chart.get('global_offset') == expected['global_offset']) and (chart.get('actual_phase') == expected['actual_phase']) and all((chart.get('roots', {}).get(key) == value for key, value in expected['roots'].items())) and (chart.get('loop', {}).get('A') == expected['loop']['A']) and (chart.get('loop', {}).get('q') == expected['loop']['q']) and (chart.get('loop', {}).get('B_over_Pstar') == expected['loop']['B_over_Pstar']) and (chart.get('loop', {}).get('inverse', [])[:2] == expected['loop']['inverse']) and (chart.get('loop', {}).get('flat_predicate') == expected['loop']['flat_predicate']) and (chart.get('loop', {}).get('phase_partial_coefficients', [])[:2] == expected['loop']['phase_partial_coefficients']) and (chart.get('finite_N', {}).get('V_N') == expected['finite_N']['V_N']) and (chart.get('finite_N', {}).get('delta_V') == expected['finite_N']['delta_V']) and (chart.get('same_background_pressure_datum') is True) and (chart.get('total_Z_phase_derivative_exact_zero') is True))
    predicates['bindings.loop_switch_charts'] = all(loop_chart_checks)
    predicates['bindings.loop_common_N'] = lp.get('N') == 18
    predicates['bindings.loop_phase_origin'] = lp.get('phase_origin') == 24954
    loop_rows = {18: {'operation': 'exact_positive_integer_expression', 'type': 'DoubleDyadicInteger', 'definition': 'N=2^(2^J)', 'exact_positive_integer': True, 'Z_independent': True, 'one_common_N_for_all_charts_and_control_denominators': True}, 21: {'operation': 'current_original_source_parameter', 'name': 'eta', 'definition': 'same canonical original source eta_log; all current owners share identical source tuples', 'quantity_is_exact_original_parameter_not_range_endpoint': True, 'Z_independent': True, 'source_attribute': 'eta_log'}, 22: {'operation': 'current_original_source_parameter', 'name': 'd_star', 'definition': 'same canonical original source dstar_log; all current owners share identical source tuples', 'quantity_is_exact_original_parameter_not_range_endpoint': True, 'Z_independent': True, 'source_attribute': 'dstar_log'}, 33572: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5}, 69245: {'operation': 'positive_quotient', 'numerator': 1, 'denominator': 33572, 'source_positive_certificate': 'same strictly positive Pstar'}, 99411: {'operation': 'positive_quotient', 'numerator': 1, 'denominator': 18, 'source_positive_certificate': 'same exact positive integer N on every source chart'}, 104535: {'operation': 'sum', 'arguments': [104534, 97298]}, 104536: {'operation': 'product', 'arguments': [18, 104535]}, 104537: {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 104536}, 107163: {'operation': 'sum', 'arguments': [107162, 97298]}, 107164: {'operation': 'product', 'arguments': [18, 107163]}, 107165: {'operation': 'analytic_unary', 'name': 'fractional_part', 'argument': 107164}, 87584: {'operation': 'product', 'arguments': [87154, 69245]}, 87585: {'operation': 'product', 'arguments': [87155, 69245]}, 90102: {'operation': 'product', 'arguments': [88437, 69245]}, 90103: {'operation': 'product', 'arguments': [88438, 69245]}, 106691: {'operation': 'product', 'arguments': [106437, 99411]}, 106692: {'operation': 'product', 'arguments': [106438, 99411]}, 106932: {'operation': 'sum', 'arguments': [87584, 106691]}, 106933: {'operation': 'sum', 'arguments': [87585, 106692]}, 106437: {'operation': 'original_lazy_flat_branch', 'active_body': 106431, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 106438: {'operation': 'original_lazy_flat_branch', 'active_body': 106432, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 108857: {'operation': 'original_lazy_flat_branch', 'active_body': 108851, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 108858: {'operation': 'original_lazy_flat_branch', 'active_body': 108852, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 109099: {'operation': 'product', 'arguments': [108857, 99411]}, 109100: {'operation': 'product', 'arguments': [108858, 99411]}, 109340: {'operation': 'sum', 'arguments': [90102, 109099]}, 109341: {'operation': 'sum', 'arguments': [90103, 109100]}}
    loop_rows.update({104538: {'operation': 'product', 'arguments': [87790, 48]}, 104539: {'operation': 'product', 'arguments': [87793, 48]}, 104544: {'operation': 'negative', 'argument': 104538}, 104545: {'operation': 'negative', 'argument': 104539}, 104550: {'operation': 'sum', 'arguments': [87550, 104544]}, 104551: {'operation': 'sum', 'arguments': [87553, 104545]}, 104556: {'operation': 'product', 'arguments': [87824, 48]}, 104557: {'operation': 'product', 'arguments': [87825, 48]}, 104565: {'operation': 'sum', 'arguments': [104551, 104564]}, 104600: {'operation': 'sum', 'arguments': [104557, 104599]}, 104653: {'operation': 'sum', 'arguments': [104645, 104652]}, 104562: {'operation': 'positive_quotient', 'numerator': 104550, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 104566: {'operation': 'positive_quotient', 'numerator': 104565, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 104597: {'operation': 'positive_quotient', 'numerator': 104556, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 104601: {'operation': 'positive_quotient', 'numerator': 104600, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 104638: {'operation': 'sum', 'arguments': [87997, 104632]}, 104639: {'operation': 'sum', 'arguments': [88001, 104633]}, 104644: {'operation': 'product', 'arguments': [104638, 85905]}, 104645: {'operation': 'product', 'arguments': [104639, 85905]}, 104650: {'operation': 'positive_quotient', 'numerator': 104644, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 104654: {'operation': 'positive_quotient', 'numerator': 104653, 'denominator': 87550, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 105094: {'operation': 'original_lazy_flat_branch', 'active_body': 105068, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 105095: {'operation': 'original_lazy_flat_branch', 'active_body': 105071, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 106399: {'operation': 'original_lazy_flat_branch', 'active_body': 106393, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 106400: {'operation': 'original_lazy_flat_branch', 'active_body': 106394, 'Delta': 104793, 'eta': 21, 'flat_value': 0, 'flat_predicate': 104794, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 107166: {'operation': 'product', 'arguments': [90308, 48]}, 107167: {'operation': 'product', 'arguments': [90311, 48]}, 107187: {'operation': 'sum', 'arguments': [107179, 107186]}, 107172: {'operation': 'negative', 'argument': 107166}, 107173: {'operation': 'negative', 'argument': 107167}, 107178: {'operation': 'sum', 'arguments': [90068, 107172]}, 107179: {'operation': 'sum', 'arguments': [90071, 107173]}, 107184: {'operation': 'positive_quotient', 'numerator': 107178, 'denominator': 90068, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 107188: {'operation': 'positive_quotient', 'numerator': 107187, 'denominator': 90068, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 107225: {'operation': 'sum', 'arguments': [90508, 107219]}, 107226: {'operation': 'sum', 'arguments': [90512, 107220]}, 107240: {'operation': 'sum', 'arguments': [107232, 107239]}, 107231: {'operation': 'product', 'arguments': [107225, 88600]}, 107232: {'operation': 'product', 'arguments': [107226, 88600]}, 107237: {'operation': 'positive_quotient', 'numerator': 107231, 'denominator': 90068, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 107241: {'operation': 'positive_quotient', 'numerator': 107240, 'denominator': 90068, 'source_positive_certificate': 'same strictly positive actual prescribed swirl E'}, 107573: {'operation': 'original_lazy_flat_branch', 'active_body': 107547, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 107574: {'operation': 'original_lazy_flat_branch', 'active_body': 107550, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 108819: {'operation': 'original_lazy_flat_branch', 'active_body': 108813, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 108820: {'operation': 'original_lazy_flat_branch', 'active_body': 108814, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 108857: {'operation': 'original_lazy_flat_branch', 'active_body': 108851, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 0, 'same_original_flat_extension_through_order5': True}, 108858: {'operation': 'original_lazy_flat_branch', 'active_body': 108852, 'Delta': 107272, 'eta': 21, 'flat_value': 0, 'flat_predicate': 107273, 'original_collar_coordinate': None, 'original_collar_upper': None, 'active_body_not_evaluated_on_flat_branch': True, 'ordinary_slow_Z_derivative_order': 1, 'same_original_flat_extension_through_order5': True}, 88388: {'operation': 'positive_quotient', 'numerator': 88382, 'denominator': 87550, 'source_positive_certificate': 'same actual positive core inlet swirl E'}, 88392: {'operation': 'positive_quotient', 'numerator': 88391, 'denominator': 87550, 'source_positive_certificate': 'same actual positive core inlet swirl E'}, 88423: {'operation': 'product', 'arguments': [88388, 85905]}, 88424: {'operation': 'product', 'arguments': [88392, 85905]}, 90899: {'operation': 'positive_quotient', 'numerator': 90893, 'denominator': 90068, 'source_positive_certificate': 'same actual positive core inlet swirl E'}, 90903: {'operation': 'positive_quotient', 'numerator': 90902, 'denominator': 90068, 'source_positive_certificate': 'same actual positive core inlet swirl E'}, 90934: {'operation': 'product', 'arguments': [90899, 88600]}, 90935: {'operation': 'product', 'arguments': [90903, 88600]}})
    loop_rows[106017] = {'operation': 'original_monotone_phase_inverse', 'phase_function': 105891, 'angle_variable': 'true_loop_C5_angle_first_switch', 'target_phase': 104537, 'angle_lower': 0, 'angle_upper': 117, 'phase_endpoints': [0, 1], 'derivative_positive_certificate': 'K*(1+t^2)>0; same original a/q-flat construction', 'flat_inverse': '2*pi*phase', 'Delta': 104793, 'eta': 21, 'branch_inverse_is_original_not_support_cap': True, 'true_prescribed_source_C5_implicit_inverse': True}
    loop_rows[108455] = {'operation': 'original_monotone_phase_inverse', 'phase_function': 108329, 'angle_variable': 'true_loop_C5_angle_second_switch', 'target_phase': 107165, 'angle_lower': 0, 'angle_upper': 117, 'phase_endpoints': [0, 1], 'derivative_positive_certificate': 'K*(1+t^2)>0; same original a/q-flat construction', 'flat_inverse': '2*pi*phase', 'Delta': 107272, 'eta': 21, 'branch_inverse_is_original_not_support_cap': True, 'true_prescribed_source_C5_implicit_inverse': True}
    predicates['graph.loop_fixed_rows'] = all((J3RU__row_has(loop_nodes, i, expected) for i, expected in loop_rows.items()))

    def _canonical_parameter_row(index, name, source_attribute):
        row = loop_nodes[index] if isinstance(loop_nodes, list) and index < len(loop_nodes) else None
        proof = row.get('canonical_parameter_proof') if isinstance(row, J3RU_Mapping) else None
        return isinstance(row, J3RU_Mapping) and row.get('operation') == 'current_original_source_parameter' and (row.get('name') == name) and (row.get('source_attribute') == source_attribute) and (row.get('quantity_is_exact_original_parameter_not_range_endpoint') is True) and (row.get('Z_independent') is True) and (row.get('source_family') == source_family) and isinstance(proof, J3RU_Mapping) and (proof.get('source_family') == source_family) and (proof.get('source_provider_types') == ['WholeZO3RcFiniteN', 'WholeZSwitchFiniteN', 'WholeZLongRmFiniteN', 'WholeZRmPatchFiniteN']) and (proof.get('exact_original_parameter_tuples_equal') is True) and (proof.get('parameter_ranges_not_source_field_values') is True)
    predicates['source_join.canonical_eta_dstar_rows'] = _canonical_parameter_row(21, 'eta', 'eta_log') and _canonical_parameter_row(22, 'd_star', 'dstar_log')
    native_history_rows_ok = all((J3RU__row_has(c6_nodes, row_index, {'operation': 'definite_integral', 'lower': 0, 'upper': upper, 'variable': variable, 'exact_function_integral': True, 'numerical_value_not_installed': True, 'endpoint_independent_of_Z': True, 'ordinary_slow_Z_derivative_order': order, 'true_core_inlet_complete_history': True, 'numerical_cap_not_defining_integrand': True, 'actual_six_history_rate': rate, 'actual_history_name': history_name, 'original_incoming_memory_retained': True, 'source_measure_applied_once': True, 'prescribed_bridge_chart': chart_label}) for history_name, chart_label, row_indices, rate, variable, upper in J3RU_EXPECTED_C6_HISTORY_Z2_ROWS for order, row_index in enumerate(row_indices)))
    predicates['graph.c6_native_history_Z2'] = native_history_rows_ok
    predicates['source_join.native6_history_Z2_to_MQ'] = native_history_rows_ok and isinstance(original_state, J3RU_Mapping) and (original_state.get('original_state_Z2_proved') is True) and (original_state.get('original_five_moments_Z2_proved') is True) and (original_state.get('ordinary_Z_orders') == [0, 1, 2]) and (original_state.get('owned_cell') == {'lower': 0, 'upper': '1/2'})
    first_loop = lp.get('charts', {}).get('first_switch', {})
    second_loop = lp.get('charts', {}).get('second_switch', {})
    first_c6 = c6.get('charts', {}).get('first_switch', {})
    second_c6 = c6.get('charts', {}).get('second_switch', {})
    predicates['source_join.c5_roots_actual_operand_rows'] = first_c6.get('radial', {}).get('E_y', [])[:2] == [87790, 87793] and first_c6.get('radial', {}).get('V_y', [])[:2] == [87824, 87825] and (second_c6.get('radial', {}).get('E_y', [])[:2] == [90308, 90311]) and (second_c6.get('radial', {}).get('V_y', [])[:2] == [0, 0]) and (first_loop.get('roots', {}).get('a', [])[:2] == [104562, 104566]) and (first_loop.get('roots', {}).get('b', [])[:2] == [104597, 104601]) and (first_loop.get('roots', {}).get('p1', [])[:2] == [104650, 104654]) and (first_loop.get('roots', {}).get('p2', [])[:2] == [88423, 88424]) and (second_loop.get('roots', {}).get('a', [])[:2] == [107184, 107188]) and (second_loop.get('roots', {}).get('b', [])[:2] == [0, 0]) and (second_loop.get('roots', {}).get('p1', [])[:2] == [107237, 107241]) and (second_loop.get('roots', {}).get('p2', [])[:2] == [90934, 90935])
    predicates['source_join.p2_actual_denominator_inputs'] = J3RU__row_has(loop_nodes, 88388, {'operation': 'positive_quotient', 'numerator': 88382, 'denominator': 87550, 'source_positive_certificate': 'same actual positive core inlet swirl E'}) and J3RU__row_has(loop_nodes, 88392, {'operation': 'positive_quotient', 'numerator': 88391, 'denominator': 87550, 'source_positive_certificate': 'same actual positive core inlet swirl E'}) and J3RU__row_has(loop_nodes, 88423, {'operation': 'product', 'arguments': [88388, 85905]}) and J3RU__row_has(loop_nodes, 88424, {'operation': 'product', 'arguments': [88392, 85905]}) and J3RU__row_has(loop_nodes, 90899, {'operation': 'positive_quotient', 'numerator': 90893, 'denominator': 90068, 'source_positive_certificate': 'same actual positive core inlet swirl E'}) and J3RU__row_has(loop_nodes, 90903, {'operation': 'positive_quotient', 'numerator': 90902, 'denominator': 90068, 'source_positive_certificate': 'same actual positive core inlet swirl E'}) and J3RU__row_has(loop_nodes, 90934, {'operation': 'product', 'arguments': [90899, 88600]}) and J3RU__row_has(loop_nodes, 90935, {'operation': 'product', 'arguments': [90903, 88600]}) and (first_c6.get('fields', {}).get('E', [])[:2] == [87550, 87553]) and (second_c6.get('fields', {}).get('E', [])[:2] == [90068, 90071]) and (first_c6.get('generic', {}).get('p2', [])[:2] == [88423, 88424]) and (second_c6.get('generic', {}).get('p2', [])[:2] == [90934, 90935]) and (first_loop.get('roots', {}).get('p2', [])[:2] == first_c6.get('generic', {}).get('p2', [])[:2]) and (second_loop.get('roots', {}).get('p2', [])[:2] == second_c6.get('generic', {}).get('p2', [])[:2])
    mod_rows = {25272: {'operation': 'sum', 'arguments': [1, 25271]}, 33572: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 5}, 106932: {'operation': 'sum', 'arguments': [87584, 106691]}, 106933: {'operation': 'sum', 'arguments': [87585, 106692]}, 109340: {'operation': 'sum', 'arguments': [90102, 109099]}, 109341: {'operation': 'sum', 'arguments': [90103, 109100]}, 148085: {'operation': 'product', 'arguments': [106932, 33572]}, 148086: {'operation': 'product', 'arguments': [106933, 33572]}, 154247: {'operation': 'product', 'arguments': [109340, 33572]}, 154248: {'operation': 'product', 'arguments': [109341, 33572]}}
    predicates['graph.mod_fixed_rows'] = all((J3RU__row_has(mod_nodes, i, expected) for i, expected in mod_rows.items()))
    mod_charts = mp.get('charts', {})
    predicates['bindings.mod_physical_U'] = all((isinstance(mod_charts.get(name), J3RU_Mapping) and mod_charts[name].get('physical', {}).get('U') == rows for name, rows in J3RU_EXPECTED_MOD_PHYSICAL_U.items()))
    predicates['source_join.c6_fields_to_loop_roots'] = all((lp.get('charts', {}).get(name, {}).get('roots', {}).get('E') == c6.get('charts', {}).get(name, {}).get('fields', {}).get('E', [])[:6] and lp.get('charts', {}).get(name, {}).get('roots', {}).get('V') == c6.get('charts', {}).get(name, {}).get('fields', {}).get('V', [])[:6] and (lp.get('charts', {}).get(name, {}).get('roots', {}).get('p2') == c6.get('charts', {}).get(name, {}).get('generic', {}).get('p2')) for name in ('first_switch', 'second_switch')))
    predicates['source_join.exact_U_row_pairs'] = c6.get('charts', {}).get('first_switch', {}).get('actual', {}).get('U', [])[:2] == [87154, 87155] and c6.get('endpoints', {}).get('R1', {}).get('actual', {}).get('U', [])[:2] == [88437, 88438] and (c6.get('charts', {}).get('second_switch', {}).get('actual', {}).get('U', [])[:2] == [88437, 88438]) and (mod_charts.get('first_switch', {}).get('physical', {}).get('U', [])[:2] == [148085, 148086]) and (mod_charts.get('second_switch', {}).get('physical', {}).get('U', [])[:2] == [154247, 154248])
    predicates['scope.ordinary_Z_orders'] = isinstance(original_state, J3RU_Mapping) and original_state.get('ordinary_Z_orders') == [0, 1, 2] and all((order in original_state.get('ordinary_Z_orders', []) for order in (0, 1)))
    predicates['scope.owned_Z_cell'] = isinstance(original_state, J3RU_Mapping) and original_state.get('owned_cell') == {'lower': 0, 'upper': '1/2'}
    return (predicates, details)

def J3RU__mq_guard(original_state):
    if not isinstance(original_state, J3RU_Mapping):
        return (False, {'reason': 'original_state must be a mapping'})
    proved = original_state.get('source_only') is True and original_state.get('accepted_identity_claim') is False and (original_state.get('original_state_Z2_proved') is True) and (original_state.get('original_five_moments_Z2_proved') is True) and (original_state.get('ordinary_Z_orders') == [0, 1, 2]) and (original_state.get('owned_cell') == {'lower': 0, 'upper': '1/2'})
    details = {'original_state_Z2_proved': original_state.get('original_state_Z2_proved') is True, 'original_five_moments_Z2_proved': original_state.get('original_five_moments_Z2_proved') is True, 'source_only': original_state.get('source_only') is True, 'accepted_identity_claim': original_state.get('accepted_identity_claim'), 'ordinary_Z_orders': original_state.get('ordinary_Z_orders'), 'owned_cell': original_state.get('owned_cell'), 'moments_are_external_premise_only': True}
    return (proved, details)

def J3RU__failure(predicates, details, gate):
    return {'schema': J3RU_SCHEMA, 'fragment_id': J3RU_FRAGMENT_ID, 'status': 'missing_gate', 'runtime_model_metadata': J3RU_RUNTIME_MODEL_METADATA, 'J3_R_C1_proved': False, 'J3_U_C1_proved': False, 'missing_gate': gate, 'predicates': predicates, 'details': details, 'source_only': True, 'Q_C1_proved': False, 'B_C1_proved': False, 'weighted_J3_jump_zero_proved': False, 'numeric_values_computed': False, 'source_constructor_calls': 0, 'source_producer_function_calls': 0, 'old_proof_functions_called': 0, 'fragment_was_unexecuted_before_frozen_entry': True, 'generic_DAG_walks': 0}

def J3RU__check_passed(outcome):
    if outcome is True:
        return True
    return isinstance(outcome, J3RU_Mapping) and outcome.get('passed') is True

def prove_J3_R_U_transport(loop, mod, loop_creator_text, c6_creator_text, original_state, check, *, background=None, mod_creator_text=None, bridge_creator_text=None, weighted_creator_text=None, archive_raw_sha256=None):
    """Authenticate and conditionally prove the source-only J3 R/U seam.

    ``background`` is the decoded C6 archive.  ``archive_raw_sha256`` must
    carry the pre-decode raw archive pins for ``loop``, ``mod``, and
    ``background``.  ``mod_creator_text`` and ``bridge_creator_text`` are
    required source pins; ``weighted_creator_text`` is accepted for the parent
    wrapper's inventory but is not used for this R/U leaf.  The callback contract is ``check(predicates)`` and it must
    return ``True`` or ``{'passed': True}``.  No caller has invoked this
    function during preparation of the fragment.
    """
    texts = {'loop_creator_text': loop_creator_text, 'mod_creator_text': mod_creator_text, 'c6_creator_text': c6_creator_text, 'bridge_creator_text': bridge_creator_text}
    predicates, details = J3RU__creator_ast_guards(texts)
    mq_ok, mq_details = J3RU__mq_guard(original_state)
    predicates['premise.original_state_Z2_proved'] = mq_ok
    details['original_state'] = mq_details
    archive_predicates, archive_details = J3RU__archive_guards(loop, mod, background, archive_raw_sha256, original_state=original_state)
    predicates.update(archive_predicates)
    details['archives'] = archive_details
    if not all(predicates.values()):
        first_gate = next((name for name, value in predicates.items() if value is not True))
        return J3RU__failure(predicates, details, first_gate)
    if not callable(check):
        return J3RU__failure(predicates, details, 'check_callable_required')
    try:
        check_result = check(predicates)
    except Exception as exc:
        details['check_exception'] = type(exc).__name__
        return J3RU__failure(predicates, details, 'check_raised_exception')
    predicates['parent_check_passed'] = J3RU__check_passed(check_result)
    details['parent_check_result_type'] = type(check_result).__name__
    if not predicates['parent_check_passed']:
        return J3RU__failure(predicates, details, 'parent_check_rejected_actual_predicates')
    return {'schema': J3RU_SCHEMA, 'fragment_id': J3RU_FRAGMENT_ID, 'status': 'conditional_source_theorem_proved', 'runtime_model_metadata': J3RU_RUNTIME_MODEL_METADATA, 'J3_R_C1_proved': True, 'J3_U_C1_proved': True, 'source_only': True, 'accepted_identity_claim': False, 'ordinary_Z_orders': [order for order in original_state['ordinary_Z_orders'] if order in (0, 1)], 'ordinary_Z2_premise_orders': list(original_state['ordinary_Z_orders']), 'owned_Z_cell': [original_state['owned_cell']['lower'], 0.5], 'owned_Z_cell_source_typed': dict(original_state['owned_cell']), 'source_steps': ['C6 first_switch endpoint at coordinate 73190 = 1 has Z-independent U substitutions 88437/88438.', 'C6 constructor passes end(first_switch, 1) into second_switch_fields at coordinate 75949 = 0; second chart preserves initial U and has U_y = 0.', 'Bridge source maps actual U to V and actual U_y to V_y using common Pstar row 33572.', 'C5 source roots authenticate E,V,a,b,p1,p2; finite-N V_N uses the same phase and exact positive N row 18 with the lazy flat implicit inverse recipe.', 'C3 source maps V_N back to physical U by common Pstar, yielding the exact U row pairs 148085/148086 and 154247/154248.'], 'exact_U_row_pairs': {'c6_R1_and_second_lower': [[88437, 88438]], 'c5_first_V_N': [[106932, 106933]], 'c5_second_V_N': [[109340, 109341]], 'c3_first_physical_U': [[148085, 148086]], 'c3_second_physical_U': [[154247, 154248]], 'common_Pstar': 33572, 'common_L': 25272, 'common_N': 18}, 'conditional_premises': ['original_state_Z2_proved', 'original_five_moments_Z2_proved', 'source C1 guards for E,V,E_y,V_y,p1,p2 from actual C6/bridge/C5 bindings', 'common positive Pstar and common positive N graph rows', 'parent check accepted the actual predicates'], 'Q_C1_proved': False, 'B_C1_proved': False, 'weighted_J3_jump_zero_proved': False, 'numeric_values_computed': False, 'source_constructor_calls': 0, 'source_producer_function_calls': 0, 'old_proof_functions_called': 0, 'fragment_was_unexecuted_before_frozen_entry': True, 'generic_DAG_walks': 0, 'predicates': predicates, 'details': details, 'source_hashes': {'creator_text': dict(J3RU_EXPECTED_SOURCE_SHA256), 'raw_archives': dict(J3RU_EXPECTED_RAW_ARCHIVE_SHA256), 'ast_body': {family: {'.'.join(spec): digest for spec, digest in mapping.items()} for family, mapping in J3RU_EXPECTED_AST_BODY_SHA256.items()}}}
