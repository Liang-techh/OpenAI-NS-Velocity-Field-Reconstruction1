"""Fixed original-source interpretation for the J2 mean-input seam.

The new typed contract extends the denotation of two old, distinct formal
parameters. It does not rewrite their graph IDs or prove equality for arbitrary
uninterpreted assignments. No enclosure endpoint defines an exact width here.
Only the frozen, entered-once source consumer may call ``derive``.
"""
import ast
import hashlib
import json
from dataclasses import dataclass

SCHEMA = 'original_native_J2_mean_loop_input_join_v1'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_U_endpoint_transport_v2.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
C6 = 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.json.gz'
R100 = 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz'
WIDTH = 'lei_ren_part1_paper_compliant_actual_bridge_integrals.json'
MACRO = 'lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py'
MOMENTS = 'lei_ren_part1_paper_compliant_current_original_bridge_macro_moments.py'
ENDPOINT = 'lei_ren_part1_paper_compliant_current_original_R100_endpoint.py'
SWITCH = 'lei_ren_part1_paper_compliant_current_original_first_switch_functions.py'
C6_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py'
R100_CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py'
CONNECTION = 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_control_connection.py'
SIGMA = 'lei_ren_part1_paper_compliant_outer_initial.py'
RECORD_NAMES = [PARENT, LOOP, C6, R100, WIDTH, MACRO, MOMENTS, ENDPOINT,
                SWITCH, C6_CREATOR, R100_CREATOR, CONNECTION, SIGMA]
RECORD_HASHES = {
    PARENT: 'be0cad628f44cec0713144f074cbbe4c7b0efc23feab9d3d482a116235d3737b',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    C6: '5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92',
    R100: '5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a',
    WIDTH: '7e25fa846ef77e408dea54acf331df081e96b45a702ff3c9bc87136221144455',
    MACRO: '74511ee00a655184021b0516478cffb7c8fdca830b5dc01ab44dfde2e9ad648b',
    MOMENTS: 'ab0f25f7a9c9436f8fbaba9952d538a4c861695911d4e6596eebc219ac124b85',
    ENDPOINT: 'f3453e14fac85530698bc7d3a394cc30dd22a5bbc8954d98723794a9dd6b0fef',
    SWITCH: 'a98753313fb1689f784c7a56cf09dc06df83e0b52cb738b5d4531c19f5402496',
    C6_CREATOR: '8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605',
    R100_CREATOR: '70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32',
    CONNECTION: 'aed09a89972ef9de0234e72691c4b5d33dfaa06d82d07b8df0781ef0f7a3807d',
    SIGMA: 'a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657',
}
STAGE_FLAGS = ['actual_J2_mean_input_ports_bound', 'typed_source_width_binding_installed',
               'interpreted_hbB_hbS_identity_proved',
               'first_switch_lower_control_identity_proved',
               'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['J2_mean_Ey_C1_endpoint_identity_proved',
               'J2_mean_Vy_C1_endpoint_identity_proved',
               'J2_mean_p2_C1_endpoint_identity_proved',
               'J2_Bloop_C1_endpoint_identity_proved',
               'actual_J2_U_C1_endpoint_identity_proved', 'J2_Q_identity_proved',
               'J2_B_identity_proved', 'uninterpreted_graph_parameter_identity_proved',
               'numeric_values_computed', 'old_error_terms_removed',
               'global_five_moment_closure', 'time_order_recursion_claimed',
               'full_goal_complete', 'all_native_B_jumps_zero_proved', 'old_operators_replayed']
FAMILY = {
    'actual_five_defect_family_sha256': '3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894',
    'implicit_source_sha256': '5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae',
    'datum_enclosure_sha256': 'dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010',
}
PARAMETER_DEFINITIONS = {
    8: ('hbB', 'original bridge_owner.flow.h; positive source ledger'),
    9: ('hbS', 'original switch_owner.flow.h; not merged with hbB by name'),
}
SOURCE_RECEIPT = 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json'
SOURCE_RECEIPT_SHA = 'f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158'
WIDTH_FORMULA = 'hb=epsilon_b=cstar*K^-100'
MEAN_PAIRS = {'Ey': [[72348, 87790], [72351, 87793]],
              'Vy': [[72382, 87824], [72383, 87825]],
              'p2': [[72981, 88423], [72982, 88424]]}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def source_ast(text):
    return ast.parse(text)


def has_statement(tree, scope, statement):
    """Exact syntax guard; never executes the authenticated producer."""
    for name in scope:
        matches = [node for node in tree.body
                   if isinstance(node, (ast.ClassDef, ast.FunctionDef)) and node.name == name]
        if len(matches) != 1:
            return False
        tree = matches[0]
    expected = ast.dump(ast.parse(statement).body[0], include_attributes=False)
    return sum(ast.dump(node, include_attributes=False) == expected
               for node in ast.walk(tree) if isinstance(node, ast.stmt)) == 1


@dataclass(frozen=True)
class ExactWidthDenotation:
    """An exact source recipe, deliberately without a numeric conversion API."""
    source_record_sha256: str
    source_field: str
    formula: str
    source_family_digest: str
    fixed_N0_digest: str

    @property
    def primitive_id(self):
        return digest(self.__dict__)


@dataclass(frozen=True)
class BoundWidthParameter:
    node_id: int
    name: str
    owner_path: tuple
    denotation: ExactWidthDenotation


class OriginalWidthRegistry:
    """Fixed adapter for actual original owners, not a caller-selected witness.

    Installation is a new explicit source interpretation. An old graph with
    arbitrary independent hbB/hbS assignments is outside this contract.
    Authentication and source-route guards must complete before construction.
    """
    def __init__(self, denotation, parameter_rows):
        self._rows = {i: dict(parameter_rows[i]) for i in (8, 9)}
        self._bindings = {
            8: BoundWidthParameter(8, 'hbB',
                ('OriginalBridgeMacroFunctions.owner', 'MacroFlow', 'flow.h'), denotation),
            9: BoundWidthParameter(9, 'hbS',
                ('OriginalFirstSwitchFunctions.owner', 'OriginalR100Endpoint.owner',
                 'OriginalBridgeMacroMoments.owner', 'OriginalBridgeMacroFunctions.owner',
                 'MacroFlow', 'flow.h'), denotation),
        }

    def resolve(self, node_id, authenticated_row):
        require(node_id in self._bindings, 'Width registry has no other parameter bindings')
        require(authenticated_row == self._rows[node_id], 'Width source row changed')
        return self._bindings[node_id]


def derive(records, expected_input_digest):
    require(set(records) == set(RECORD_NAMES), 'Changed record inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    checks, decoded, proof = {}, {}, []

    def check(name, condition):
        require(bool(condition), 'Source prerequisite: ' + name)
        checks[name] = True

    for name in RECORD_NAMES:
        packet = records[name]
        check('authenticated:' + name, packet['raw_sha256'] == RECORD_HASHES[name]
              and packet['raw_authenticated_before_decode'] is True
              and packet['source_path'] == 'experiments/root_st073/' + name
              and packet['is_gzip'] is name.endswith('.gz')
              and hashlib.sha256(packet['source_text'].encode('utf-8')).hexdigest()
              == packet['decoded_sha256'])
        decoded[name] = (packet['source_text'] if name.endswith('.py')
                         else json.loads(packet['source_text']))
    parent, loop, c6, r100, width = [decoded[n] for n in (PARENT, LOOP, C6, R100, WIDTH)]
    check('accepted_parent', parent['component_result']['passed'] is True
          and parent['source_only'] is True and parent['component_result']['source_only'] is True
          and parent['full_goal_complete'] is False and len(parent['source_pins']) == 1837
          and parent['R_identity_jumps'] == [2, 3, 4, 6, 8, 9, 10]
          and parent['B_identity_scope'] == dict(all_fourteen_B_identities_open=True,
              newly_proved_B_identity_nodes=[], unresolved_B_identity_nodes=list(range(522, 536)))
          and parent['Q_identity_scope']['J2_Q_identity_proved'] is False
          and all(parent['component_result'][flag] is False for flag in
              ['numeric_values_computed', 'old_error_terms_removed', 'global_five_moment_closure',
               'time_order_recursion_claimed', 'J2_B_identity_proved', 'J2_Q_identity_proved',
               'all_native_B_jumps_zero_proved', 'old_operators_replayed'])
          and parent['component_result']['all_old_terms_retained'] is True
          and parent['component_result']['source_inspection_audit']['accepted_ancestor_operator_calls'] == 0)
    check('inherited_source_pins', all(parent['source_pins'][n] == h
          for n, h in RECORD_HASHES.items() if n != PARENT))
    check('retained_ledger', all(parent['retained_construction_data'][k] == v for k, v in {
        'cumulative_terms': 212, 'raw_signed_terms': 213, 'reduced_terms': 183,
        'six_signed_channels_retained': True, 'all_old_error_data_retained': True,
        'reserved_contribution_registered_once': True}.items()))
    check('inherited_U_scope_open', parent['component_result']['J2_mean_E_C1_endpoint_identity_proved']
          and parent['component_result']['J2_mean_V_C1_endpoint_identity_proved']
          and parent['component_result']['J2_phase_endpoint_identity_proved']
          and parent['component_result']['actual_J2_U_C1_endpoint_identity_proved'] is False)
    check('source_family', loop['source_family'] == FAMILY and c6['source_family'] == FAMILY
          and r100['source_family'] == FAMILY
          and {k: width[k] for k in FAMILY} == FAMILY)
    ln, cn, rn = [d['exact_graph_nodes'] for d in (loop, c6, r100)]
    parameters = {i: ln[i] for i in (8, 9)}
    for i, (name, definition) in PARAMETER_DEFINITIONS.items():
        row = parameters[i]
        check('separate_literal_width_parameter:' + name,
              set(row) == {'operation', 'name', 'definition', 'source_family',
              'fixed_leading_input_N0', 'quantity_is_exact_original_parameter_not_range_endpoint',
              'Z_independent', 'current_source_receipt', 'current_source_receipt_sha256'}
              and row['operation'] == 'current_original_source_parameter'
              and row['name'] == name and row['definition'] == definition
              and row['source_family'] == FAMILY and row['Z_independent'] is True
              and row['quantity_is_exact_original_parameter_not_range_endpoint'] is True
              and row['current_source_receipt'] == SOURCE_RECEIPT
              and row['current_source_receipt_sha256'] == SOURCE_RECEIPT_SHA
              and type(row['fixed_leading_input_N0']) is int and row['fixed_leading_input_N0'] > 0
              and row == cn[i] == rn[i])
    check('common_fixed_N0', parameters[8]['fixed_leading_input_N0']
          == parameters[9]['fixed_leading_input_N0'])
    check('exact_width_recipe', width['exact_source_width'] == WIDTH_FORMULA
          and 'source_log_hb_enclosure' in width)

    syntax = {
        MACRO: ["original = self.records['actual_bridge_integrals']",
                "logh = previous.read_interval(c, original['source_log_hb_enclosure'])",
                "flow = MacroFlow(c, logh, logP2, logF02, logRa, c.ln(100) - logRa, weight)",
                "self.logs = tuple(c.mpf(v) for v in (logh, logP2, logF02, logRa, 0))",
                "self.h = self.factor((1, 0, 0, 0, 0))",
                "return prior.ScaledEnclosure(prior.FormalScale(self.logs, powers, offset), 1, self.ledger)"],
        MOMENTS: ["flow, proof = self.fields.owner(label)",
                  "owner = ActualMacroMoments(flow, inlet)", "self.flow = flow"],
        ENDPOINT: ["moment_owner, proof = self.upstream.owner(label)",
                   "flow = moment_owner.flow", "self.owners[label] = (flow, proof, value)"],
        SWITCH: ["flow, proof, value = self.upstream.owner(label)",
                 "self.owners[label] = FirstSwitchFunctions(flow, value['normalized_actual_fields'], value['normalized_actual_own_six_moments'], self.cells)",
                 "self.flow = flow"],
        C6_CREATOR: ["self.hb = hb = inlet_owner.previous.lower.parameters['hb']"],
        R100_CREATOR: ["self.bridge_hb = bridge_owner.hb",
                       "self.hb = hb = self.native_parameter('hbS', 'original switch_owner.flow.h; not merged with hbB by name')",
                       "chi_axial = g.mul(hb, g.sub(g.one, sigma(g, phase)))"],
    }
    scopes = {
        MACRO: [('OriginalBridgeMacroFunctions', 'owner')] * 3
               + [('MacroFlow', '__init__')] * 2 + [('MacroFlow', 'factor')],
        MOMENTS: [('OriginalBridgeMacroMoments', 'owner')] * 2
                 + [('ActualMacroMoments', '__init__')],
        ENDPOINT: [('OriginalR100Endpoint', 'owner')] * 3,
        SWITCH: [('OriginalFirstSwitchFunctions', 'owner')] * 2
                + [('FirstSwitchFunctions', '__init__')],
        C6_CREATOR: [('CurrentOriginalPrescribedBridgeC6Functions', '__init__')],
        R100_CREATOR: [('CurrentOriginalPrescribedR100RmC6Functions', '__init__')] * 2
                      + [('first_switch_fields',)],
    }
    for name, statements in syntax.items():
        tree = source_ast(decoded[name])
        for i, statement in enumerate(statements):
            check('source_route_syntax:' + name + ':' + str(i),
                  has_statement(tree, scopes[name][i], statement))
    check('original_bridge_parameter_definition_declared',
          PARAMETER_DEFINITIONS[8][1] in decoded[CONNECTION])
    check('separate_switch_nodes_preserved',
          r100['prescribed_R100_Rm_C6_bindings']['definition'][
              'switch_hbS_and_preceding_bridge_hbB_source_nodes_not_aliased'] is True)
    proof.append(dict(rule='authenticated_original_owner_route',
        conclusion='Both actual owners use MacroFlow.h from the global width recipe, independent of label',
        enclosure_values_compared=False, old_graph_alias_inferred=False))
    denotation = ExactWidthDenotation(RECORD_HASHES[WIDTH], 'exact_source_width',
        WIDTH_FORMULA, digest(FAMILY), digest(parameters[8]['fixed_leading_input_N0']))
    registry = OriginalWidthRegistry(denotation, parameters)
    left, right = [registry.resolve(i, parameters[i]) for i in (8, 9)]
    check('typed_registry_separate_nodes', left.node_id == 8 and right.node_id == 9
          and left.name != right.name and left.owner_path != right.owner_path)
    check('interpreted_common_width_primitive', left.denotation == right.denotation
          and left.denotation.primitive_id == right.denotation.primitive_id)
    proof.append(dict(rule='install_fixed_actual_source_interpretation',
        premise='authenticated original owner routes and exact source recipe',
        conclusion='interpret(hbB)=interpret(hbS)=cstar*K^-100',
        primitive_id=denotation.primitive_id, distinct_formal_nodes=[8, 9],
        new_source_contract=True, equality_of_uninterpreted_graph_parameters=False))

    fixed_rows = {
        0: dict(operation='exact_rational', numerator=0, denominator=1),
        1: dict(operation='exact_rational', numerator=1, denominator=1),
        73190: dict(operation='bound_variable', name='true_prescribed_first_switch_s'),
        74429: dict(operation='negative', argument=74428),
        74430: dict(operation='sum', arguments=[1, 74429]),
        87090: dict(operation='product', arguments=[9, 74430]),
    }
    for i, row in fixed_rows.items():
        check('first_switch_control_row:' + str(i), ln[i] == row and rn[i] == row)
    sig = ln[74428]
    check('literal_sigma_definition', sig == rn[74428]
          and sig['operation'] == 'exact_original_source_flat_sigma'
          and sig['argument'] == 73190 and sig['exact_function_not_its_interval_cap'] is True
          and sig['definition'] == '0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'
          and sig['source_binding'] == dict(module=SIGMA, module_sha256=RECORD_HASHES[SIGMA],
              function='stable_sigma', AST_sha256='4eca34dd2d6bd23b5bff76496682bece54924d3145ac99e95ed9066d9572ab3c'))
    sigma_ast = next(n for n in source_ast(decoded[SIGMA]).body
                     if isinstance(n, ast.FunctionDef) and n.name == 'stable_sigma')
    check('sigma_source_AST_binding', hashlib.sha256(ast.dump(ast.Module(
          body=[sigma_ast], type_ignores=[])).encode('utf-8')).hexdigest()
          == sig['source_binding']['AST_sha256'])
    proof.append(dict(rule='exact_sigma_lower_branch_then_typed_substitution',
        endpoint=dict(variable=73190, value=0, Z_independent=True),
        sigma_node=74428, sigma_value=0, control_node=87090,
        control_expression='interpret(hbS)*(1-sigma(0))=interpret(hbB)',
        ordinary_Z_orders=[0, 1], width_Z1=0, control_Z1=0,
        no_drive_or_radial_input_identity_inferred=True))

    macro_source = c6['prescribed_bridge_C6_bindings']['macro']
    switch_source = r100['prescribed_R100_Rm_C6_bindings']['charts']['first_switch']
    mean_roots = loop['prescribed_loop_C5_bindings']['charts']
    check('mean_Ey_ports', macro_source['radial']['E_y'][:2] == [72348, 72351]
          and switch_source['radial']['E_y'][:2] == [87790, 87793])
    check('mean_Vy_ports', macro_source['radial']['V_y'][:2] == [72382, 72383]
          and switch_source['radial']['V_y'][:2] == [87824, 87825])
    check('mean_p2_ports', macro_source['generic']['p2'][:2] == [72981, 72982]
          and switch_source['generic']['p2'][:2] == [88423, 88424]
          and mean_roots['frozen_macro']['roots']['p2'][:2] == [72981, 72982]
          and mean_roots['first_switch']['roots']['p2'][:2] == [88423, 88424])
    binding_scope = dict(new_source_contract=True, old_graph_nodes_rewritten=False,
        distinct_formal_nodes=[8, 9], interpreted_parameter_identity_proved=True,
        uninterpreted_parameter_identity_proved=False,
        primitive=dict(**denotation.__dict__, primitive_id=denotation.primitive_id),
        bindings=[dict(node_id=v.node_id, name=v.name, owner_path=list(v.owner_path),
                       primitive_id=v.denotation.primitive_id) for v in (left, right)],
        numeric_witness_selected=False, enclosure_endpoint_comparison=False,
        positive_exact_width_source_recipe=True, finite_positive_value_not_materialized=True)
    remaining = ['prove original comparison D_over_R and signed drive body congruence at J2',
        'join original mean Ey,Vy,p2 C1; then same-loop composition and full U',
        'physical Q and all fourteen native weighted B identities remain open',
        'terminal functional five-moment closure and all dependent physical layers remain incomplete']
    return dict(schema=SCHEMA, component='mean_loop_input_join', source_only=True, passed=True,
        **dict.fromkeys(STAGE_FLAGS, True), **dict.fromkeys(FALSE_FLAGS, False),
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        input_digest=expected_input_digest, runtime_model_metadata=None,
        owned_Z_cell=parent['component_result']['owned_Z_cell'],
        stage_flags={flag: True for flag in STAGE_FLAGS},
        nonclaims={flag: False for flag in FALSE_FLAGS},
        source_family=FAMILY, source_binding_scope=binding_scope,
        six_input_scope=dict(pairs=MEAN_PAIRS, all_six_joined=False,
            finite_N_modified_total_fields_used=False, body_congruence_still_required=True),
        source_checks=checks, proof_graph=proof, proof_graph_sha256=digest(proof),
        unresolved_gates=[], remaining_construction_gates=remaining,
        inherited_U_transport_scope=parent['U_transport_scope'],
        inherited_R_identity_jumps=[2, 3, 4, 6, 8, 9, 10],
        unresolved_B_identity_nodes=list(range(522, 536)),
        retained_construction_data=parent['retained_construction_data'],
        source_inspection_audit=dict(new_source_consumer_calls=1,
            accepted_ancestor_operator_calls=0, raw_gzip_authenticated_before_decode=True,
            generic_DAG_walks=0, numerical_evaluations=0, quadratures=0,
            recurrence_calls=0, runtime_model_metadata=None))
