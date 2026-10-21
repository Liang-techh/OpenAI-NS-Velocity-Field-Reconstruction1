"""Fixed-source real exp/log proof of J2 radius equality; no DAG evaluator.

This certificate checks the literal archived radius definitions, then applies
elementary real identities with explicit domain guards. It does not compare U,
physical Q, or B, and does not execute any accepted source consumer.
"""
import base64
import hashlib
import json

SCHEMA = 'original_native_J2_radius_identity_v1'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_remaining_B_seam_sources.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
GAP = 'tmp/original_native_remaining_B_seam_sources/native_J2_saved_radius_DAG_gap_2026_10_21.json'
CREATOR = 'lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py'
WIDTH = 'lei_ren_part1_paper_compliant_actual_bridge_integrals.json'
RECORD_NAMES = [PARENT, LOOP, GAP, CREATOR, WIDTH]
RECORD_HASHES = {
    PARENT: 'f765b451fcb962cbf59824db2d0f7e994b3de72c30c3ee4bf46ab9f01c0cc1e6',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    GAP: '2dcd59604323857b8d0851c1ac2c301cc491ad3fd1657a048e49e2ecca2329d1',
    CREATOR: '74511ee00a655184021b0516478cffb7c8fdca830b5dc01ab44dfde2e9ad648b',
    WIDTH: '7e25fa846ef77e408dea54acf331df081e96b45a702ff3c9bc87136221144455',
}
STAGE_FLAGS = ['literal_radius_sources_bound', 'positive_real_log_domain_proved',
               'same_source_hbB_cancellation_proved', 'J2_radius_C1_identity_proved',
               'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['numeric_values_computed', 'old_error_terms_removed',
               'global_five_moment_closure', 'time_order_recursion_claimed',
               'full_goal_complete', 'J2_B_identity_proved',
               'all_native_B_jumps_zero_proved', 'old_operators_replayed']
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
    """Enter once on five authenticated records; return a source-only proof."""
    require(set(records) == set(RECORD_NAMES), 'Changed input inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    decoded, checks = {}, {}

    def check(name, condition):
        require(bool(condition), 'Source prerequisite: ' + name)
        checks[name] = True

    for name in RECORD_NAMES:
        packet = records[name]
        check('authenticated:' + name, packet['raw_sha256'] == RECORD_HASHES[name]
              and packet['raw_authenticated_before_decode'] is True
              and hashlib.sha256(packet['source_text'].encode('utf-8')).hexdigest()
              == packet['decoded_sha256'])
        decoded[name] = packet['source_text'] if name == CREATOR else json.loads(packet['source_text'])
    parent, loop, gap = (decoded[name] for name in [PARENT, LOOP, GAP])
    check('accepted_parent', parent['component_result']['passed'] is True
          and parent['full_goal_complete'] is False
          and parent['component_result_is_projection'] is True)
    check('parent_inventory', len(parent['source_pins']) == 1830
          and parent['source_pins'][LOOP] == RECORD_HASHES[LOOP])
    check('finite_width_contract_pins', all(parent['source_pins'][name] == RECORD_HASHES[name]
          for name in [CREATOR, WIDTH]))
    creator, width = decoded[CREATOR], decoded[WIDTH]
    check('finite_width_creator_guard', all(fragment in creator for fragment in [
        'class MacroFlow:',
        'self.logs=tuple(c.mpf(v) for v in (logh,logP2,logF02,logRa,0))',
        'if any(not mp.isfinite(x) for row in self.logs+(self.Y,self.weight) for x in ep(row)):',
        "raise ValueError('Finite fixed source log records required')",
        'self.h=self.factor((1,0,0,0,0))',
        'return prior.ScaledEnclosure(prior.FormalScale(self.logs,powers,offset),1,self.ledger)',
        "original['source_log_hb_enclosure']"]))
    check('exact_width_not_cap', width['exact_source_width'] == 'hb=epsilon_b=cstar*K^-100'
          and all(width[name] == FAMILY[name] for name in FAMILY))
    log_h = width['source_log_hb_enclosure']
    for side in ['lower', 'upper']:
        # A finite mpf tuple is (-1)^s*m*2^e. Inspect its finite integer
        # representation only: never materialize 2^e, h, or an endpoint value.
        row = log_h[side + '_exact_mpf_tuple']
        check('finite_real_log_h:' + side, type(row) is list and len(row) == 4
              and all(type(item) is int for item in row)
              and row[0] in [0, 1] and row[1] > 0
              and row[3] == row[1].bit_length())
    check('all_B_still_open', parent['all_fourteen_B_identities_open'] is True
          and parent['newly_proved_weighted_jump_nodes'] == []
          and parent['unresolved_weighted_jump_nodes'] == list(range(522, 536)))
    check('inherited_J0_J1', parent['source_interface']['proved_weighted_jump_nodes'] == [520, 521]
          and parent['source_interface']['old_17_K_integrals_and_2_outer_B_ports_retained'] is True)
    check('retained_ledger', parent['retained_construction_data'] == {
        'cumulative_terms': 212, 'raw_signed_terms': 213, 'reduced_terms': 183,
        'six_signed_channels_retained': True, 'all_old_error_data_retained': True,
        'reserved_contribution_registered_once': True, 'inherited_errors_deleted': False,
        'P1_differs_from_m15': True, 'physical_Q_differs_from_loop_q': True,
        'accepted_physical_baseline_commit': '8f36f19e3c524b424a9286c46cc8e0b855f75307'})
    gap_raw = base64.b64decode(parent['saved_gap_evidence']['raw_base64'], validate=True)
    check('gap_exact_parent_bytes', hashlib.sha256(gap_raw).hexdigest() == RECORD_HASHES[GAP]
          and gap_raw == records[GAP]['source_text'].encode('utf-8'))
    check('gap_scope', gap['schema'] == 'native_J2_saved_radius_DAG_gap_v1'
          and gap['seam']['weighted_jump_node'] == 522
          and gap['seam']['ordinary_Z_order'] == 0 and gap['seam']['field'] == 'R')
    nodes = loop['exact_graph_nodes']
    check('loop_source_family', loop['source_family'] == FAMILY and len(nodes) == 121834)
    # The fixed pattern is the whole relevant definition. No recursive walk or
    # canonicalizer is used, and no exp/log is numerically evaluated.
    patterns = {
        0: {'operation': 'exact_rational', 'numerator': 0, 'denominator': 1},
        1: {'operation': 'exact_rational', 'numerator': 1, 'denominator': 1},
        2: {'operation': 'exact_rational', 'numerator': 40, 'denominator': 1},
        3: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 2},
        4: {'operation': 'exact_rational', 'numerator': 11, 'denominator': 1},
        5: {'operation': 'sum', 'arguments': [3, 4]},
        13: {'operation': 'exact_rational', 'numerator': 4, 'denominator': 1},
        14: {'operation': 'product', 'arguments': [13, 5]},
        48: {'operation': 'exact_rational', 'numerator': 2, 'denominator': 1},
        643: {'operation': 'exact_rational', 'numerator': 1000, 'denominator': 1},
        31709: {'operation': 'exact_rational', 'numerator': 100, 'denominator': 1},
        35562: {'operation': 'sum', 'arguments': [14, 643]},
        35564: {'operation': 'negative', 'argument': 35562},
        35565: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 35564},
        65604: {'operation': 'product', 'arguments': [13, 35565]},
        66981: {'operation': 'analytic_unary', 'name': 'log', 'argument': 31709},
        66982: {'operation': 'analytic_unary', 'name': 'log', 'argument': 65604},
        66983: {'operation': 'negative', 'argument': 66982},
        66984: {'operation': 'sum', 'arguments': [66981, 66983]},
        66985: {'operation': 'product', 'arguments': [48, 8]},
        66986: {'operation': 'negative', 'argument': 66985},
        66987: {'operation': 'sum', 'arguments': [66984, 66986]},
        70265: {'operation': 'bound_variable', 'name': 'original_prescribed_C6_macro_fraction'},
        70266: {'operation': 'product', 'arguments': [66987, 70265]},
        71446: {'operation': 'sum', 'arguments': [66985, 70266]},
        71447: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 71446},
        71448: {'operation': 'product', 'arguments': [65604, 71447]},
        73190: {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'},
        85903: {'operation': 'product', 'arguments': [9, 73190]},
        85904: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 85903},
        85905: {'operation': 'product', 'arguments': [31709, 85904]},
    }
    for node, pattern in patterns.items():
        check('literal_node:' + str(node), nodes[node] == pattern)
    for node, name, definition in (
        (8, 'hbB', 'original bridge_owner.flow.h; positive source ledger'),
        (9, 'hbS', 'original switch_owner.flow.h; not merged with hbB by name')):
        parameter = nodes[node]
        check('source_parameter:' + name,
              parameter['operation'] == 'current_original_source_parameter'
              and parameter['name'] == name and parameter['definition'] == definition
              and parameter['source_family'] == FAMILY
              and parameter['quantity_is_exact_original_parameter_not_range_endpoint'] is True
              and parameter['Z_independent'] is True
              and type(parameter['fixed_leading_input_N0']) is int
              and parameter['fixed_leading_input_N0'] > 0
              and parameter['current_source_receipt_sha256']
              == 'f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158')
    check('parameters_not_aliased', nodes[8] != nodes[9]
          and nodes[8]['fixed_leading_input_N0'] == nodes[9]['fixed_leading_input_N0'])
    charts = loop['prescribed_loop_C5_bindings']['charts']
    seam = parent['component_result']['seam_comparisons'][0]
    check('actual_J2_seam', seam['jump_index'] == 2 and seam['weighted_jump_node'] == 522
          and seam['left_chart'] == 'frozen_macro' and seam['right_chart'] == 'first_switch'
          and seam['function_identity_proved'] is False)
    for side, chart, endpoint, coord, radius, port, endpoint_value in (
        ('left', 'frozen_macro', 'upper', 70265, 71448, 61, 1),
        ('right', 'first_switch', 'lower', 73190, 85905, 77, 0)):
        binding = charts[chart]
        source = seam[side + '_endpoint_source']
        saved = seam['endpoint_rows']['R_Z0'][side]
        gap_side = gap['seam'][side]
        check('endpoint_binding:' + side, binding['coordinate'] == coord
              and binding['domain'] == [0, 1] and binding['radius'] == radius
              and source['coordinate']['source_nodes'] == [coord]
              and source['domain']['source_nodes'] == [0, 1]
              and source['radius']['source_nodes'] == [radius]
              and source['radius']['weighted_node'] == port
              and source['radius']['source_record'] == LOOP
              and saved['source_node'] == radius and saved['coordinate_nodes'] == [coord]
              and saved['endpoint'] == endpoint and saved['port_weighted_node'] == port
              and saved['canonical_supported'] is True
              and saved['Z_independent_source'] is True
              and gap_side['source_node'] == radius
              and gap_side['canonical_digest'] == saved['canonical_digest'])
        check('endpoint_constant:' + side, nodes[binding['domain'][endpoint_value]]
              == {'operation': 'exact_rational', 'numerator': endpoint_value, 'denominator': 1})
    # Each proof step states the real-algebra rule and earlier premises. The
    # source checks above establish the complete fixed syntax for these rules.
    proof = []

    def step(name, rule, premises, conclusion, source_nodes):
        require(all(item in [row['id'] for row in proof] for item in premises), 'Forward proof premise')
        proof.append(dict(id=name, rule=rule, premises=premises,
                          conclusion=conclusion, source_nodes=source_nodes))

    step('parameters', 'original MacroFlow finite-log guard and exact h=exp(log h) source contract', [],
         'hbB and hbS bind finite positive real original flow.h; their source IDs remain distinct', [8, 9])
    step('domain', 'finite real constants; real exp is finite and positive', [],
         'D=1000+4*(exp(40)+11) is real; A=4*exp(-D)>0; 100>0', [35562, 65604, 31709])
    step('logs', 'real log domain for positive arguments', ['domain'],
         'log(A) and log(100) are defined; A is nonzero', [66981, 66982])
    step('left_endpoint', 'substitute native macro coordinate 70265=1', ['logs', 'parameters'],
         'R_left_upper(Z)=A*exp(2*hbB+((log(100)-log(A))-2*hbB))', [70265, 71448])
    step('same_source_cancel', 'associativity and additive inverse of identical node 66985',
         ['left_endpoint'], 'R_left_upper(Z)=A*exp(log(100)-log(A))', [66985, 66986, 71446])
    step('exp_log', 'exp(a-b)=exp(a)/exp(b); exp(log(x))=x for x>0',
         ['logs', 'same_source_cancel'], 'R_left_upper(Z)=A*(100/A)', [71447, 66981, 66982])
    step('factor_cancel', 'field cancellation at A nonzero', ['exp_log', 'domain'],
         'R_left_upper(Z)=100', [65604, 71448])
    step('right_endpoint', 'substitute native switch coordinate 73190=0; exp(0)=1', ['parameters'],
         'R_right_lower(Z)=100*exp(hbS*0)=100', [73190, 85903, 85904, 85905])
    step('radius_equal', 'transitivity of exact real function equality',
         ['factor_cancel', 'right_endpoint'], 'R_left_upper(Z)=R_right_lower(Z)=100', [71448, 85905])
    step('ordinary_Z', 'differentiate constant function; endpoints and parameters Z independent',
         ['radius_equal'], 'd_Z R_left_upper=d_Z R_right_lower=0 on the inherited owned Z cell',
         [8, 9, 0, 1, 71448, 85905])
    identity = dict(weighted_jump_node=522, jump_index=2, field='R',
        ordinary_Z_orders=[0, 1], exact_endpoint_rational={'numerator': 100, 'denominator': 1},
        exact_ordinary_Z_derivative={'numerator': 0, 'denominator': 1},
        source_function_identity_proved=True, numeric_value=None,
        left_source_node=71448, right_source_node=85905,
        left_coordinate_node=70265, right_coordinate_node=73190,
        endpoint_substitutions={'70265': 1, '73190': 0},
        hbB_source_node=8, hbS_source_node=9, source_parameters_aliased=False,
        saved_parent_R_row=seam['endpoint_rows']['R_Z0'])
    return dict(schema=SCHEMA, component='J2_radius_identity', passed=True,
        source_only=True, ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        owned_Z_cell=parent['owned_Z_cell'], input_digest=expected_input_digest,
        source_family=FAMILY, source_checks=checks, proof_graph=proof,
        proof_graph_digest=digest(proof), radius_identity=identity,
        finite_real_parameter_contract=dict(creator=CREATOR, creator_sha256=RECORD_HASHES[CREATOR],
            exact_width_record=WIDTH, exact_width_record_sha256=RECORD_HASHES[WIDTH],
            exact_source_width=width['exact_source_width'], finite_log_h_witness=log_h,
            finite_enclosure_used_only_for_domain_not_as_function_value=True,
            exact_source_parameter_values_selected=False, source_width_cap_used=False,
            parameter_values_aliased=False),
        fixed_raw_source_rows={str(i): nodes[i] for i in list(patterns) + [8, 9]},
        inherited_context=dict(parent_report=PARENT, parent_report_sha256=RECORD_HASHES[PARENT],
            inherited_R_identity_jump_indices=[3, 4, 6, 8, 9, 10],
            accepted_R_identity_jump_indices=[2, 3, 4, 6, 8, 9, 10],
            native_B_jump_zero_nodes=[520, 521],
            unresolved_B_weighted_jump_nodes=list(range(522, 536)),
            native_IN_formula=parent['source_interface']['native_IN_formula']),
        retained_construction_data=parent['retained_construction_data'],
        stage_flags={flag: True for flag in STAGE_FLAGS},
        nonclaims={flag: False for flag in FALSE_FLAGS},
        source_pins=dict(RECORD_HASHES), source_pin_count=len(RECORD_HASHES),
        source_inspection_audit=dict(new_source_consumer_calls=1,
            accepted_ancestor_operator_calls=0, general_DAG_walker_calls=0,
            numeric_exp_log_evaluations=0, endpoint_canonicalizer_calls=0),
        unresolved_gates=[], remaining_construction_gates=[
            'J2 actual U and physical Q endpoint C1 source identities',
            'all fourteen remaining B source identities',
            'native/old density identities and terminal functional closure',
            'admissible stress, temporal recursion, corrections, independent Cartesian NS'],
        runtime_model_metadata=None, **{flag: True for flag in STAGE_FLAGS},
        **{flag: False for flag in FALSE_FLAGS})
