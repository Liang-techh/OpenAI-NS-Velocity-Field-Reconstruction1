"""Exact J2 endpoint transport and remaining physical U difference.

Only a new frozen source consumer may call derive. No ancestor operators,
quadratures, numerical evaluators or broad endpoint DAG walks are used.
"""
import hashlib
import json

SCHEMA = 'original_native_J2_U_endpoint_transport_v2'
PARENT = 'lei_ren_part1_paper_compliant_current_original_native_J2_radius_identity.json'
LOOP = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
MOD = 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
CREATOR = 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py'
RECORD_NAMES = [PARENT, LOOP, MOD, CREATOR]
RECORD_HASHES = {
    PARENT: 'ccb7afca508bc5b31706356fd9058290c587a527ec534ec701824aaba8da431e',
    LOOP: '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa',
    MOD: 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a',
    CREATOR: '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d',
}
STAGE_FLAGS = ['actual_J2_U_ports_bound', 'J2_mean_V_C1_endpoint_identity_proved',
               'J2_mean_E_C1_endpoint_identity_proved', 'J2_phase_endpoint_identity_proved',
               'J2_U_difference_C1_reduction_proved',
               'all_old_terms_retained', 'accepted_ancestors_not_reexecuted']
FALSE_FLAGS = ['numeric_values_computed', 'old_error_terms_removed',
               'global_five_moment_closure', 'time_order_recursion_claimed',
               'full_goal_complete', 'J2_B_identity_proved', 'J2_Q_identity_proved',
               'all_native_B_jumps_zero_proved', 'old_operators_replayed',
               'actual_J2_U_C1_endpoint_identity_proved',
               'J2_finite_N_delta_V_C1_endpoint_identity_proved']
FAMILY = {
    'actual_five_defect_family_sha256': '3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894',
    'implicit_source_sha256': '5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae',
    'datum_enclosure_sha256': 'dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010',
}

BASE_ROWS = {70152: {'operation': 'function_substitution',
         'expression': 68910,
         'variable': 66988,
         'value': 48,
         'Z_independent_substitution': True},
 70153: {'operation': 'function_substitution',
         'expression': 68911,
         'variable': 66988,
         'value': 48,
         'Z_independent_substitution': True},
 70265: {'operation': 'bound_variable', 'name': 'original_prescribed_C6_macro_fraction'},
 71699: {'operation': 'definite_integral',
         'integrand': 71692,
         'variable': 'macro_actual_U_t',
         'lower': 0,
         'upper': 70265,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'numerical_cap_not_defining_integrand': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 0,
         'true_core_inlet_complete_history': True,
         'actual_prescribed_axial_equation': True,
         'all_hydro_pressure_swirl_sectors_retained': True,
         'source_measure_applied_once': True,
         'prescribed_bridge_chart': 'macro'},
 71700: {'operation': 'definite_integral',
         'integrand': 71693,
         'variable': 'macro_actual_U_t',
         'lower': 0,
         'upper': 70265,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'numerical_cap_not_defining_integrand': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 1,
         'true_core_inlet_complete_history': True,
         'actual_prescribed_axial_equation': True,
         'all_hydro_pressure_swirl_sectors_retained': True,
         'source_measure_applied_once': True,
         'prescribed_bridge_chart': 'macro'},
 71706: {'operation': 'negative', 'argument': 71699},
 71707: {'operation': 'negative', 'argument': 71700},
 71713: {'operation': 'sum', 'arguments': [70152, 71706]},
 71714: {'operation': 'sum', 'arguments': [70153, 71707]},
 72994: {'operation': 'function_substitution',
         'expression': 71713,
         'variable': 70265,
         'value': 1,
         'Z_independent_substitution': True},
 72995: {'operation': 'function_substitution',
         'expression': 71714,
         'variable': 70265,
         'value': 1,
         'Z_independent_substitution': True},
 73190: {'operation': 'bound_variable', 'name': 'true_prescribed_first_switch_s'},
 87133: {'operation': 'function_substitution',
         'expression': 87126,
         'variable': 73190,
         'value': 74235,
         'Z_independent_substitution': True},
 87134: {'operation': 'function_substitution',
         'expression': 87127,
         'variable': 73190,
         'value': 74235,
         'Z_independent_substitution': True},
 87140: {'operation': 'product', 'arguments': [87133, 9]},
 87141: {'operation': 'product', 'arguments': [87134, 9]},
 87147: {'operation': 'definite_integral',
         'integrand': 87140,
         'variable': 'true_first_switch_axial_t',
         'lower': 0,
         'upper': 73190,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'numerical_cap_not_defining_integrand': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 0,
         'true_core_inlet_complete_history': True,
         'original_first_switch_axial_control': True,
         'source_measure_applied_once': True},
 87148: {'operation': 'definite_integral',
         'integrand': 87141,
         'variable': 'true_first_switch_axial_t',
         'lower': 0,
         'upper': 73190,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'numerical_cap_not_defining_integrand': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 1,
         'true_core_inlet_complete_history': True,
         'original_first_switch_axial_control': True,
         'source_measure_applied_once': True},
 87154: {'operation': 'sum', 'arguments': [72994, 87147]},
 87155: {'operation': 'sum', 'arguments': [72995, 87148]},
 72142: {'operation': 'product', 'arguments': [71713, 69245]},
 72143: {'operation': 'product', 'arguments': [71714, 69245]},
 87584: {'operation': 'product', 'arguments': [87154, 69245]},
 87585: {'operation': 'product', 'arguments': [87155, 69245]}}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode('utf-8')).hexdigest()


def derive(records, expected_input_digest):
    require(set(records) == set(RECORD_NAMES), 'Changed record inventory')
    require(digest(records) == expected_input_digest, 'Changed input digest')
    checks, decoded = {}, {}

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
    parent, loop, mod = (decoded[name] for name in [PARENT, LOOP, MOD])
    check('accepted_parent', parent['component_result']['passed'] is True
          and parent['full_goal_complete'] is False and parent['R_identity_jumps'] == [2,3,4,6,8,9,10])
    check('inherited_actual_radius_identity', parent['component_result']['J2_radius_C1_identity_proved'] is True
          and parent['radius_scope']['left']['source_node'] == 71448
          and parent['radius_scope']['right']['source_node'] == 85905
          and parent['radius_scope']['left']['radius_value'] == 100
          and parent['radius_scope']['right']['radius_value'] == 100
          and parent['component_result']['positive_real_log_domain_proved'] is True)
    check('inherited_inventory', len(parent['source_pins']) == 1834
          and all(parent['source_pins'][name] == RECORD_HASHES[name] for name in [LOOP, MOD, CREATOR]))
    check('same_source_family', parent['component_result']['source_family'] == FAMILY
          and loop['source_family'] == FAMILY and mod['source_family'] == FAMILY)
    check('all_B_still_open', parent['B_identity_scope']['unresolved_B_identity_nodes'] == list(range(522,536))
          and parent['B_identity_scope']['newly_proved_B_identity_nodes'] == [])
    check('retained_ledger', parent['retained_construction_data']['cumulative_terms'] == 212
          and parent['retained_construction_data']['raw_signed_terms'] == 213
          and parent['retained_construction_data']['reduced_terms'] == 183
          and parent['retained_construction_data']['six_signed_channels_retained'] is True
          and parent['retained_construction_data']['all_old_error_data_retained'] is True
          and parent['retained_construction_data']['reserved_contribution_registered_once'] is True)
    lnodes, mnodes = loop['exact_graph_nodes'], mod['exact_graph_nodes']
    for node, row in BASE_ROWS.items():
        check('literal_mean_row:' + str(node), lnodes[node] == row and mnodes[node] == row)
    common = {
        0: dict(operation='exact_rational',numerator=0,denominator=1),
        1: dict(operation='exact_rational',numerator=1,denominator=1),
        2: dict(operation='exact_rational',numerator=40,denominator=1),
        3: dict(operation='analytic_unary',name='exp',argument=2),
        4: dict(operation='exact_rational',numerator=11,denominator=1),
        5: dict(operation='sum',arguments=[3,4]),
        33572: dict(operation='analytic_unary',name='exp',argument=5),
        69245: dict(operation='positive_quotient',numerator=1,denominator=33572,
                    source_positive_certificate='same strictly positive Pstar'),
        99411: dict(operation='positive_quotient',numerator=1,denominator=18,
                    source_positive_certificate='same exact positive integer N on every source chart'),
        104304: dict(operation='sum',arguments=[72142,104063]),
        104305: dict(operation='sum',arguments=[72143,104064]),
        106932: dict(operation='sum',arguments=[87584,106691]),
        106933: dict(operation='sum',arguments=[87585,106692]),
        104063: dict(operation='product',arguments=[103809,99411]),
        104064: dict(operation='product',arguments=[103810,99411]),
        106691: dict(operation='product',arguments=[106437,99411]),
        106692: dict(operation='product',arguments=[106438,99411]),
    }
    for node, row in common.items():
        check('literal_common_row:' + str(node), lnodes[node] == row and mnodes[node] == row)
    N = lnodes[18]
    check('exact_finite_positive_N', N == mnodes[18]
          and N['operation'] == 'exact_positive_integer_expression'
          and N['type'] == 'DoubleDyadicInteger' and type(N['J']) is int and N['J'] >= 0
          and N['definition'] == 'N=2^(2^J)' and N['exact_positive_integer'] is True
          and N['Z_independent'] is True and N['outer_integer_and_inner_power_not_materialized'] is True
          and N['one_common_N_for_all_charts_and_control_denominators'] is True)
    creator = decoded[CREATOR]
    check('producer_serializes_unadmitted_flags_false',
          '**dict.fromkeys(GATES+OPEN,False)' in creator
          and 'prescribed_loop_C5_bindings=owner.snapshot' in creator
          and loop['current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed'] is False
          and loop['current_original_prescribed_loop_C5_numeric_containment_installed'] is False)
    definition = loop['prescribed_loop_C5_bindings']['definition']
    for flag in ['true_prescribed_background_source_not_legacy_recipe_derivative_padding',
                 'all_ten_signed_root_q_inverse_A_B_and_phase_primitive_rows_through_Z5',
                 'inverse_rows_from_exact_implicit_normalized_series_not_interval_endpoints',
                 'one_original_common_N_eta_dstar_and_global_log_R_rminus_phase',
                 'all_five_signed_finite_N_density_cross_terms_and_expm1_retained',
                 'new_functions_preserve_original_fixed_leading_source_family',
                 'numeric_C5_containment_and_mixed_y_bounds_not_certified']:
        check('actual_source_definition:' + flag, definition[flag] is True)
    physical = mod['modified_Q4_pressure_C3_bindings']['charts']
    geometry = loop['prescribed_loop_C5_bindings']['charts']
    physical_rows = {
        140685: dict(operation='product',arguments=[104304,33572]),
        140686: dict(operation='product',arguments=[104305,33572]),
        148085: dict(operation='product',arguments=[106932,33572]),
        148086: dict(operation='product',arguments=[106933,33572]),
    }
    for node, row in physical_rows.items():
        check('physical_U_row:' + str(node), mnodes[node] == row)
    for chart, coord, radius, U, V, deltaV in (
        ('frozen_macro',70265,71448,[140685,140686],[72142,72143],[104063,104064]),
        ('first_switch',73190,85905,[148085,148086],[87584,87585],[106691,106692])):
        check('actual_chart_binding:' + chart, geometry[chart]['coordinate'] == coord
              and geometry[chart]['domain'] == [0,1] and geometry[chart]['radius'] == radius
              and physical[chart]['physical']['U'][:2] == U
              and geometry[chart]['roots']['V'][:2] == V
              and geometry[chart]['finite_N']['delta_V'][:2] == deltaV)
    proof = []

    def step(name, rule, premises, conclusion, nodes):
        require(all(p in [item['id'] for item in proof] for p in premises), 'Forward proof premise')
        proof.append(dict(id=name,rule=rule,premises=premises,conclusion=conclusion,source_nodes=nodes))

    step('source_domain','exact positive finite exp source and fixed finite integer N',[],
         'Pstar=exp(exp(40)+11)>0 and N=2^(2^J)>0 are Z independent; no value is evaluated',
         [33572,69245,18,99411])
    step('inherited_radius','authenticate and reuse the accepted parent source theorem; no reexecution',[],
         'R_left_upper=R_right_lower=100, their ordinary Z derivatives vanish, and Ra>0',
         [71448,85905,65604])
    step('empty_switch_integrals','standard exact integral over identical bounds; native s=0',[],
         'I_87147(Z,0)=I_87148(Z,0)=0; source measure applied once; no nonempty integral evaluated',
         [73190,87147,87148])
    step('same_inlet','nested same-source endpoint substitution at macro x=1',
         ['empty_switch_integrals'],
         '87154(Z,0)=72994(Z)=71713(Z,1); 87155(Z,0)=72995(Z)=71714(Z,1)',
         [72994,72995,71713,71714,87154,87155])
    step('mean_V','multiplication by common inverse Pstar; explicit ordinary Z0 and Z1 rows',
         ['source_domain','same_inlet'],
         'Vbase_left_upper=Vbase_right_lower through ordinary Z order1',
         [72142,72143,87584,87585])
    phase_rows = {
        101905: dict(operation='positive_quotient',numerator=71448,denominator=65604,
                     source_positive_certificate='same positive actual R/Ra'),
        104533: dict(operation='positive_quotient',numerator=85905,denominator=65604,
                     source_positive_certificate='same positive actual R/Ra'),
        101906: dict(operation='analytic_unary',name='log',argument=101905),
        104534: dict(operation='analytic_unary',name='log',argument=104533),
        101907: dict(operation='sum',arguments=[101906,97298]),
        104535: dict(operation='sum',arguments=[104534,97298]),
        101908: dict(operation='product',arguments=[18,101907]),
        104536: dict(operation='product',arguments=[18,104535]),
        101909: dict(operation='analytic_unary',name='fractional_part',argument=101908),
        104537: dict(operation='analytic_unary',name='fractional_part',argument=104536),
        97298: dict(operation='negative',argument=24954),
    }
    for node, row in phase_rows.items():
        check('literal_phase_row:' + str(node), lnodes[node] == row and mnodes[node] == row)
    for chart, offset, phase in [('frozen_macro',101907,101909),('first_switch',104535,104537)]:
        check('phase_binding:' + chart, geometry[chart]['global_offset'] == offset
              and geometry[chart]['actual_phase'] == phase
              and geometry[chart]['total_Z_phase_derivative_exact_zero'] is True)
    step('same_phase','reuse accepted radius equality and positive Ra; congruence of log, sum, product and fractional_part',
         ['source_domain','inherited_radius'],
         'Both endpoint phases are fractional_part(N*(log(100/Ra)+97298)); phase_Z=0. No phase value is selected.',
         list(phase_rows))
    mean_E = _mean_E_certificate(lnodes, mnodes, geometry, check, step)
    # Preserve the unresolved source branch; common mean transport cannot
    # establish equality of the full finite-N physical velocity.
    oscillatory = _oscillatory_certificate(lnodes, mnodes, geometry, check, step)
    step('U_difference','subtract common mean V, then factor the same Pstar and inverse N',
         ['mean_V',oscillatory['final_proof_node']],
         'For j=0,1: U_left,j(1)-U_right,j(0)=(Pstar/N)*(Bloop_over_Pstar_left,j(1)-Bloop_over_Pstar_right,j(0))',
         [140685,140686,148085,148086,103809,103810,106437,106438])
    return dict(schema=SCHEMA,component='J2_U_endpoint_transport_v2',passed=True,
        source_only=True,ordinary_Z_orders=[0,1],ordinary_Z_not_temporal_n=True,
        owned_Z_cell=parent['owned_Z_cell'],input_digest=expected_input_digest,source_family=FAMILY,
        source_checks=checks,proof_graph=proof,proof_graph_digest=digest(proof),
        source_definition_contract=dict(producer=CREATOR,producer_sha256=RECORD_HASHES[CREATOR],
            definition_metadata=definition,top_level_C5_installed=False,
            numeric_containment_installed=False,serialized_false_flags_preserved=True,
            source_definitions_distinct_from_numeric_or_physical_admission=True),
        U_difference_identity=dict(weighted_jump_node=522,jump_index=2,field='U',ordinary_Z_orders=[0,1],
            exact_difference_reduction_proved=True,full_U_endpoint_equality_proved=False,
            source_difference_is_symbolic=True,left_source_nodes=[140685,140686],
            right_source_nodes=[148085,148086],left_coordinate_node=70265,right_coordinate_node=73190,
            endpoint_substitutions={'70265':1,'73190':0},numeric_value=None,
            difference_formula='U_left,j(1)-U_right,j(0)=(Pstar/N)*(Bloop_over_Pstar_left,j(1)-Bloop_over_Pstar_right,j(0))',
            common_Pstar_node=33572,common_inverse_N_node=99411,
            retained_Bloop_over_Pstar_left_nodes=[103809,103810],
            retained_Bloop_over_Pstar_right_nodes=[106437,106438],
            Bloop_over_Pstar_differs_from_native_weighted_primitive_B=True),
        mean_E_certificate=mean_E, oscillatory_certificate=oscillatory,
        fixed_raw_source_rows={str(i):lnodes[i] for i in set(BASE_ROWS)|set(common)|set(phase_rows)|{18}},
        fixed_physical_U_rows={str(i):mnodes[i] for i in physical_rows},
        inherited_context=dict(parent_report=PARENT,parent_report_sha256=RECORD_HASHES[PARENT],
            accepted_R_identity_jump_indices=[2,3,4,6,8,9,10],native_B_jump_zero_nodes=[520,521],
            unresolved_B_weighted_jump_nodes=list(range(522,536)),
            native_IN_formula='I_N=B_lastupper-B_firstlower-sum(J_2..J_15)+sum17(integral K_i dy)'),
        retained_construction_data=parent['retained_construction_data'],
        stage_flags={flag:True for flag in STAGE_FLAGS},nonclaims={flag:False for flag in FALSE_FLAGS},
        source_pins=dict(RECORD_HASHES),source_pin_count=len(RECORD_HASHES),
        source_inspection_audit=dict(new_source_consumer_calls=1,accepted_ancestor_operator_calls=0,
            general_DAG_walker_calls=0,numeric_exp_log_evaluations=0,endpoint_canonicalizer_calls=0),
        unresolved_gates=[],remaining_construction_gates=['J2 mean radial E_y, V_y and p2 endpoint input identities',
            'J2 finite-N oscillatory source join and full U endpoint equality', 'J2 physical Q endpoint source identity',
            'all fourteen remaining B source identities','native/old density and terminal functional closure',
            'stress, temporal recursion, corrections, independent Cartesian NS'],
        runtime_model_metadata=None,**{flag:True for flag in STAGE_FLAGS},**{flag:False for flag in FALSE_FLAGS})


MEAN_E_ROWS = {48: {'operation': 'exact_rational', 'numerator': 2, 'denominator': 1},
 72090: {'operation': 'product', 'arguments': [48, 71448]},
 72091: {'operation': 'analytic_unary', 'name': 'sqrt', 'argument': 72090},
 87532: {'operation': 'product', 'arguments': [48, 85905]},
 87533: {'operation': 'analytic_unary', 'name': 'sqrt', 'argument': 87532},
 72092: {'operation': 'product', 'arguments': [72091, 69245]},
 87534: {'operation': 'product', 'arguments': [87533, 69245]},
 72093: {'operation': 'product', 'arguments': [35657, 72092]},
 72094: {'operation': 'product', 'arguments': [35658, 72092]},
 87535: {'operation': 'product', 'arguments': [35657, 87534]},
 87536: {'operation': 'product', 'arguments': [35658, 87534]},
 72108: {'operation': 'product', 'arguments': [72093, 71499]},
 72109: {'operation': 'product', 'arguments': [72093, 71502]},
 72110: {'operation': 'product', 'arguments': [72094, 71499]},
 72111: {'operation': 'sum', 'arguments': [72109, 72110]},
 87550: {'operation': 'product', 'arguments': [87535, 86947]},
 87551: {'operation': 'product', 'arguments': [87535, 86950]},
 87552: {'operation': 'product', 'arguments': [87536, 86947]},
 87553: {'operation': 'sum', 'arguments': [87551, 87552]},
 72987: {'operation': 'function_substitution',
         'expression': 71499,
         'variable': 70265,
         'value': 1,
         'Z_independent_substitution': True},
 72988: {'operation': 'function_substitution',
         'expression': 71502,
         'variable': 70265,
         'value': 1,
         'Z_independent_substitution': True},
 86913: {'operation': 'definite_integral',
         'integrand': 86906,
         'variable': 'true_first_switch_angular_t',
         'lower': 0,
         'upper': 73190,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 0,
         'true_core_inlet_complete_history': True,
         'numerical_cap_not_defining_integrand': True,
         'original_first_switch_full_current_radius_direction': True},
 86914: {'operation': 'definite_integral',
         'integrand': 86907,
         'variable': 'true_first_switch_angular_t',
         'lower': 0,
         'upper': 73190,
         'exact_function_integral': True,
         'numerical_value_not_installed': True,
         'endpoint_independent_of_Z': True,
         'ordinary_slow_Z_derivative_order': 1,
         'true_core_inlet_complete_history': True,
         'numerical_cap_not_defining_integrand': True,
         'original_first_switch_full_current_radius_direction': True},
 86920: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 86913},
 86921: {'operation': 'product', 'arguments': [86914, 86920]},
 86947: {'operation': 'product', 'arguments': [72987, 86920]},
 86948: {'operation': 'product', 'arguments': [72987, 86921]},
 86949: {'operation': 'product', 'arguments': [72988, 86920]},
 86950: {'operation': 'sum', 'arguments': [86948, 86949]},
 35657: {'operation': 'analytic_unary', 'name': 'exp', 'argument': 35656},
 35658: {'operation': 'product', 'arguments': [35657, 35653]}}

OSCILLATORY_ROWS = {102165: {'operation': 'sum', 'arguments': [102159, 25304]},
 102166: {'operation': 'exact_real_comparison', 'operator': '>=', 'left': 102165, 'right': 21},
 104793: {'operation': 'sum', 'arguments': [104787, 25304]},
 104794: {'operation': 'exact_real_comparison', 'operator': '>=', 'left': 104793, 'right': 21},
 103809: {'operation': 'original_lazy_flat_branch',
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
 103810: {'operation': 'original_lazy_flat_branch',
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
 106437: {'operation': 'original_lazy_flat_branch',
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
 106438: {'operation': 'original_lazy_flat_branch',
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


def _mean_E_certificate(lnodes, mnodes, geometry, check, step):
    for node, row in MEAN_E_ROWS.items():
        check('literal_mean_E_row:' + str(node), lnodes[node] == row and mnodes[node] == row)
    check('mean_E_chart_roots', geometry['frozen_macro']['roots']['E'][:2] == [72108,72111]
          and geometry['first_switch']['roots']['E'][:2] == [87550,87553])
    step('empty_angular_integrals','exact equal-bound integrals at native first_switch s=0',[],
         '86913=86914=0, exp(86913)=1 and 86914*exp(86913)=0 at s=0',
         [73190,86913,86914,86920,86921])
    step('same_F_inlet','exact macro endpoint aliases; multiply and add the empty-integral factors',
         ['empty_angular_integrals'],
         '86947(Z,0)=72987(Z)=71499(Z,1); 86950(Z,0)=72988(Z)=71502(Z,1)',
         [72987,72988,86947,86950])
    step('mean_E','reuse inherited R=100 and ordinary Z derivative0; positive sqrt(2R) and common source factors',
         ['source_domain','inherited_radius','same_F_inlet'],
         'Emean_left,j(1)=Emean_right,j(0), j=0,1; no nonempty macro angular integral evaluated',
         [72090,72091,87532,87533,72108,72111,87550,87553])
    return dict(source_function_identity_proved=True,ordinary_Z_orders=[0,1],
                left_source_nodes=[72108,72111],right_source_nodes=[87550,87553],
                nonempty_macro_integral_preserved=True,
                fixed_raw_source_rows={str(i):lnodes[i] for i in MEAN_E_ROWS})


def _oscillatory_certificate(lnodes, mnodes, geometry, check, step):
    for node, row in OSCILLATORY_ROWS.items():
        check('literal_retained_loop_row:' + str(node), lnodes[node] == row and mnodes[node] == row)
    for chart, B, predicate in [('frozen_macro',[103809,103810],102166),
                                 ('first_switch',[106437,106438],104794)]:
        check('retained_loop_binding:' + chart, geometry[chart]['loop']['B_over_Pstar'][:2] == B
              and geometry[chart]['loop']['flat_predicate'] == predicate
              and geometry[chart]['loop']['original_collar_upper'] is None)
    eta = lnodes[21]
    check('eta_is_source_parameter_not_literal_21', eta == mnodes[21]
          and eta['operation'] == 'current_original_source_parameter' and eta['name'] == 'eta'
          and eta['source_family'] == FAMILY and eta['Z_independent'] is True
          and eta['quantity_is_exact_original_parameter_not_range_endpoint'] is True)
    step('finite_N_terms_retained','literal finite-N source product; no unproved lazy predicate evaluated',[],
         'deltaV_left,j=Bloop_over_Pstar_left,j/N and deltaV_right,j=Bloop_over_Pstar_right,j/N, j=0,1; Bloop endpoint difference remains unproved',
         [104063,104064,106691,106692,103809,103810,106437,106438,99411])
    return dict(final_proof_node='finite_N_terms_retained',
        endpoint_identity_proved=False,flat_predicate_truth='unproved',
        common_loop_input_join='unproved',numeric_values_selected=False,
        lazy_active_bodies_evaluated=False,source_parameters_approximated=False,
        eta_source_parameter_node=21,eta_is_literal_21=False,
        retained_left_Bloop_over_Pstar_nodes=[103809,103810],
        retained_right_Bloop_over_Pstar_nodes=[106437,106438],
        next_required_inputs=['mean E_y C1','mean V_y C1','generic p2 C1'],
        source_Bloop_is_distinct_from_native_weighted_primitive_B=True,
        fixed_raw_source_rows={str(i):lnodes[i] for i in OSCILLATORY_ROWS})
