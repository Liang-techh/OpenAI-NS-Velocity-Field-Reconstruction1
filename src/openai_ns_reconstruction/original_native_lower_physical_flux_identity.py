"""Source-bound lower physical RQy congruence, with explicit axial scope.

Read archived source definitions and saved premises as data. Only the new
universal C1 numerator identity is normalized. No archived helper is called.
"""
from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_compact_2Rc_flux_identity.json'
SCALE = PREFIX + 'native_physical_scale_source_identity.json'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
TRUE = PREFIX + 'modified_Q4_pressure_C3_functions.json.gz'
TC = PREFIX + 'modified_Q4_pressure_C3_functions_check.json'
LOOP = PREFIX + 'prescribed_loop_C5_functions.json.gz'
LC = PREFIX + 'prescribed_loop_C5_functions_check.json'
CORR = PREFIX + 'corrected_reserved_m15_source.json'
SEAM = PREFIX + 'n1_first_micro_seam.json'
EC = PREFIX + 'n1_preRc_endpoint_C1_check.json'
PRES = PREFIX + 'n1_preRp_moment_functions.py'
TS = PREFIX + 'modified_Q4_pressure_C3_functions.py'
LS = PREFIX + 'prescribed_loop_C5_functions.py'
PART = PREFIX + 'n1_partial_mean_Q.py'
GEN = PREFIX + 'Rm_generic_inputs.py'
SWITCH = PREFIX + 'switch_finite_N.py'
ES = PREFIX + 'n1_preRc_endpoint_C1.py'
SS = 'src/openai_ns_reconstruction/original_n1_first_micro_seam.py'
FS = 'src/openai_ns_reconstruction/original_n1_first_micro_prefix.py'
CS = 'src/openai_ns_reconstruction/original_native_collar_transport.py'
SCHEMA = 'original_native_lower_physical_flux_identity_v1'
RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_native_compact_2Rc_flux_identity.json': 'f5a7ccd09b895d89e360d6b94f9a71b091e4c6002246d86bbd2ad288ec15b717', 'lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json': 'b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz': 'b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json': '189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz': '8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json': '2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7', 'lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json': 'fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5', 'lei_ren_part1_paper_compliant_current_original_n1_first_micro_seam.json': '5d03e68d02b322a7d1e8ebee9529f4aeb57199afcf939dfa9fc7bb566a1edefa', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1_check.json': '02dec83e6e5b4c52d8aadc7e5f8dc359b2322dc62e5aa3ce895f60e8f855ecb2', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py': '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py': '69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d', 'lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q.py': 'd826c150d021fa318201ce9b393b8efb92b459546e4ce8377c2e852f2d65cbc6', 'lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs.py': '62df6a6513b6e9c521cb6cf7cff25710dec17bf04f7b23f4bc0f15dbc1a57a88', 'lei_ren_part1_paper_compliant_current_original_switch_finite_N.py': '15a7844a1171518cd73912301407b7b7a137ecfb2f63219605283870c26e3b8d', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.py': '828afe7bd79f91d8ccbeb989711147306b41abe400344b4b9fa3f55aaa14eb37', 'src/openai_ns_reconstruction/original_n1_first_micro_seam.py': '7decf0492712c53742cd49f37a75238294760cffa1332568841e4d61b56fe236', 'src/openai_ns_reconstruction/original_n1_first_micro_prefix.py': 'f48b46035b1d9a3129bff8fc5125c9da8f1fdaf6667a71bd78a00cd55a603892', 'src/openai_ns_reconstruction/original_native_collar_transport.py': '2545bf6864a89f7ddb6a6ce040df2fd5b26b6a747ae1a176597b705727fefdff'}
DECODED_RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz': '0a35efead3cb3ea587d1c44a57d70ec78b11b4a5a5305483c8d5b7a041d91f75', 'lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz': 'f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2'}
STAGE_FLAGS = ('actual_native_initial_physical_RQy_bound', 'corrected_original_lower_physical_RQy_bound_once', 'lower_physical_RQy_source_C1_identity_proved_on_owned_Z_cell', 'signed_terminal_Rc_flux_correlation_preserved')
FALSE_FLAGS = ('native_compact_m15_source_identity_proved', 'native_compact_V_m_source_identity_proved', 'native_active_m15_source_identity_proved', 'terminal_difference_zero_proved', 'old_error_terms_removed', 'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1', 'global_five_moment_closure', 'legacy_feedback_reference_alias_proved', 'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved', 'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed', 'full_goal_complete', 'old212_adapter_source_identity_proved', 'independent_terminal_function_match_complete', 'native_vs_registered_preRp_m15_identity_proved')

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True, allow_nan=False).encode()).hexdigest()

def syntax(tree, function, *needles):
    bodies = [ast.unparse(row) for row in ast.walk(tree) if isinstance(row, ast.FunctionDef) and row.name == function]
    return any((all((needle in body for needle in needles)) for body in bodies))

class FreshGraph:
    """Small proof graph. Source ports are atomic, never expanded."""

    def __init__(self):
        self.nodes = []

    def node(self, operation, arguments=(), **metadata):
        ref = len(self.nodes)
        self.nodes.append(dict(operation=operation, arguments=list(arguments), **metadata))
        return ref

    def constant(self, value):
        q = Fraction(value)
        return self.node('exact_rational', numerator=q.numerator, denominator=q.denominator)

    def add(self, *args):
        return self.node('sum', args)

    def mul(self, *args):
        return self.node('product', args)

    def sub(self, left, right):
        return self.add(left, self.mul(self.constant(-1), right))

    def polynomial(self, ref):
        row = self.nodes[ref]
        if row['operation'] == 'exact_rational':
            return {(): Fraction(row['numerator'], row['denominator'])}
        if row['operation'] == 'sum':
            out = {}
            for arg in row['arguments']:
                for key, value in self.polynomial(arg).items():
                    out[key] = out.get(key, Fraction(0)) + value
            return {key: value for key, value in out.items() if value}
        if row['operation'] == 'product':
            out = {(): Fraction(1)}
            for arg in row['arguments']:
                product = {}
                for a, av in out.items():
                    for b, bv in self.polynomial(arg).items():
                        key = tuple(sorted(a + b))
                        product[key] = product.get(key, Fraction(0)) + av * bv
                out = product
            return {key: value for key, value in out.items() if value}
        return {(ref,): Fraction(1)}


def derive(source_records, *, expected_input_digest):
    checks, data, trees = {}, {}, {}

    def check(key, condition):
        checks[key] = bool(condition)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    for name, expected in RECORD_HASHES.items():
        row = source_records[name]
        check('authenticated:' + name, row['raw_sha256'] == expected
            and row['raw_authenticated_before_decode'] is True
            and row['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, expected))
        if name.endswith('.py'):
            trees[name] = ast.parse(row['source_text'])
        else:
            data[name] = json.loads(row['source_text'])
    parent = data[PARENT]['component_result']
    scale = data[SCALE]['component_result']
    pre, true, loop = (data[k] for k in (PRE, TRUE, LOOP))
    pn, tn = pre['exact_graph_nodes'], true['exact_graph_nodes']
    f = true['modified_Q4_pressure_C3_bindings']['charts']['first_micro']
    init = true['modified_Q4_pressure_C3_bindings']['initial']
    lc = loop['prescribed_loop_C5_bindings']['charts']['first_micro']
    ld = loop['prescribed_loop_C5_bindings']['definition']
    corr = data[CORR]
    inputs, construction = (corr[k] for k in ('corrected_terminal_source_inputs', 'corrected_reserved_m15_source'))
    proof = inputs['physical_endpoint_source_proof']
    seam_report = data[SEAM]
    seam = seam_report['actual_first_micro_seam']
    family = parent['source_family']
    check('accepted_parent_upper_flux_cancellation_and_retained_Rc_term', parent['passed'] is True
        and parent['common_K_2Rc_flux_difference_C1_zero_proved'] is True
        and parent['signed_terminal_with_Rc_flux_correlation_defined'] is True)
    check('same_exact_source_family', all(row['source_family'] == family
        for row in (pre, true, loop, seam_report)) and proof['source_family'] == family)
    check('accepted_actual_physical_source_and_native_derivative', data[TC]['all_passed'] is True
        and data[TC]['checks']['actual_source']['exact_existing_seven_background_seed_rows_and_unique_core_source'] is True
        and data[TC]['checks']['actual_source']['Pstar_once_and_physical_R_Rsquared_measures_once'] is True
        and data[TC]['checks']['total_radial']['nonunit_Jacobian_and_N_cancellation_in_Aphi_Bphi_checked'] is True)
    check('actual_common_prefix_and_native_initial_flux_ports', tn[:len(pn)] == pn
        and init['R_Q_y'][:2] == [129567, 129568]
        and f['physical']['Q_y'][:2] == [129400, 129404]
        and f['physical']['U_y'][:2] == [129223, 129224]
        and f['physical']['m'][:3] == [129229, 129230, 129231])
    check('actual_native_RQy_product_and_Z_independent_pullback',
        tn[129559] == dict(operation='product', arguments=[129400, 68644])
        and tn[129560] == dict(operation='product', arguments=[129404, 68644])
        and all(tn[node] == dict(operation='function_substitution', expression=expr,
            variable=66988, value=97295, Z_independent_substitution=True)
            for node, expr in ((129567, 129559), (129568, 129560))))
    check('actual_native_total_physical_y_and_mass_ODE',
        syntax(trees[TS], 'total_radial_functions',
            "a.trim(bg['radial']['V_y'])", "a.scale(slow['B_over_Pstar'], invN)",
            "chart['loop']['B_over_Pstar_phi']", 'physical_y_equals_native_partial_divided_by_Jacobian=True')
        and syntax(trees[TS], 'physical_pressure_fields', "U = a.scale(fields['V'], Pstar)",
            "Uy = a.scale(radial['V_y'], Pstar)", 'my = a.sub(U, m)',
            'Q = div(numerator(U, m))', 'Qy = div(numerator(Uy, my))'))
    check('original_loop_flat_branch_is_radial_not_a_Z_sample', data[LC]['all_passed'] is True
        and data[LC]['checks']['lazy_flat']['lazy_original_collar_and_shear_flat_rows'] == 144
        and ld['original_first_micro_start_sc_half_and_collar_three_sc_quarters'] is True
        and ld['correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged'] is True
        and ld['new_functions_preserve_original_fixed_leading_source_family'] is True
        and lc['coordinate'] == 66988 and lc['domain'] == [97295, 1]
        and syntax(trees[LS], 'loop_functions', "collar = g.mul(c('3/4'), sc)",
            "operator='<=', left=coordinate.node, right=collar.node")
        and syntax(trees[LS], '__init__', "incoming = {key: a.fixed(0) for key in RATES}",
            "fields[key] = a.add(a.trim(chart['fields'][key]), history[key])"))
    check('owned_axial_cell_and_exact_finite_N_initial_predicates',
        proof['endpoint'] == 'first_micro/sc_half singleton'
        and proof['incoming_history']['exact_Z_cell'] == ['0', '.5']
        and all(proof[key] is True for key in ('source_owned_collar_q_C0_Z_exact_zero',
            'finite_N_initial_history_zero', 'finite_N_predecessor_history_retained',
            'finite_N_local_driver_exact_zero', 'actual_V_y_row_0_used',
            'actual_V_y_row_1_retained', 'actual_V_y_rows_0_through_3_factorial_converted_once',
            'common_Pstar_applied_once'))
        and proof['incoming_history']['all_four_incoming_N_scaled_mean_rows_exact_zero'] is True
        and all(row['exact_zero'] is True for row in proof['incoming_history']['incoming_rows'])
        and proof['flat_driver']['source_owned_initial_condition']['exact_zero_correction_follows_from_defining_Duhamel_lower_bound'] is True
        and proof['flat_driver']['source_owned_initial_condition']['actual_background_incoming_histories_remain_nonzero'] is True
        and proof['flat_driver']['source_owned_initial_condition']['q_A_B_and_Z_rows_flat_on_same_source_two_sided_left_collar'] is True)
    check('accepted_seam_source_map_not_its_numerical_packet', seam_report['checks']['passed'] is True
        and corr['input_hashes'][SEAM] == RECORD_HASHES[SEAM]
        and seam['passed'] is True and seam['exact_algebra']['passed'] is True
        and seam['new_cross_map']['smooth_native_identity_on_closed_first_window'] is True
        and seam['new_cross_map']['finite_N_correction_lower_jet_zero'] is True
        and seam['new_cross_map']['derivative_and_own_mean_source_not_replaced_by_comparison'] is True
        and seam_report['actual_native_collar_transport']['native_defining_leading_function_equals_canonical_collar'] is True
        and seam_report['actual_native_collar_transport']['actual_own_mean_and_pressure_equal_canonical_actual_histories'] is True
        and seam['old_six_numeric_terms_and_all_error_caps_unchanged'] is True)
    check('same_physical_normalization_from_actual_source_recipe',
        scale['passed'] is True and scale['old_preRp_Pstar_source_identity_proved'] is True
        and scale['source_aliases']['Pstar']['provider_node'] == 25203
        and scale['source_aliases']['Pstar']['underlying_formula_node'] == 33572
        and syntax(trees[GEN], 'recover_inputs', "invS = f.factor((0, -0.5, 0, 0, 0))",
            "V = scale(raw['velocity']['axial'][0], invS)",
            "Vy = scale(raw['velocity']['axial'][1], invS)")
        and syntax(trees[SWITCH], 'recover_source', "axial=[V, source['V_y']]"))
    check('old_lower_physical_source_quantities_are_defining_rows',
        all(pn[node]['operation'] == 'current_original_physical_source_quantity'
            and pn[node]['quantity'] == 'V_y' and pn[node]['coordinate'] == 29
            and pn[node]['ordinary_slow_Z_derivative_order'] == k
            and pn[node]['source_family'] == family
            and pn[node]['defining_quantity_not_a_range_value'] is True
            and pn[node]['finite_N_lower_initial_history_and_local_driver_exact_zero'] is True
            for k, node in enumerate((30363, 30364)))
        and syntax(trees[PRES], 'lower_velocity_y_source',
            "values = packet['actual_background_source']['V_y']",
            'math.factorial(k)', "endpoint.c0._zero_first_incoming_proof(fields, tuple(ends))")
        and syntax(trees[PRES], 'active_functions', 'qy = algebra.sub(qp, q)')
        and syntax(trees[PRES], 'q_plus_q_y', 'Uz = source.C1Function(U.Z, U.ZZ)'))
    patch, mapping = inputs['saved_patch'], construction['mapping']
    check('saved_physical_scale_patch_consumed_once_and_relocated',
        patch == seam['corrected_physical_endpoint_graph']
        and inputs['saved_patch_digest'] == seam['corrected_physical_endpoint_graph_digest']
        and patch['exact_replacements'] == {'30365': 30363, '30366': 30364}
        and patch['corrected_pairs']['lower_R_Qy'] == [31228, 31240]
        and mapping['31228'] == 53450 and mapping['31240'] == 53462
        and patch['old_graph_unchanged'] is True
        and patch['corrected_Uy'] == 'actual_background_source.V_y; already physical, no additional Pstar'
        and inputs['old_graphs_unchanged'] is True
        and construction['passed'] is True and construction['source_only'] is True)
    delta_rows = construction['delta_nodes']
    base = construction['graph_patch']['base_graph_node_count']
    check('actual_final_lower_flux_rows_in_retained_namespace', base == 53444
        and delta_rows[53450-base] == dict(operation='product', arguments=[30361, 53449])
        and delta_rows[53462-base] == dict(operation='product', arguments=[30361, 53461])
        and pre['defining_bindings']['active']['finite_endpoints'][0]['R_q_y'] == [30442, 30443])
    check('same_physical_radius_and_ordinary_endpoints',
        proof['source_quantity_paths']['R'] == 'original_generic_source.original_physical_radius'
        and seam['new_cross_map']['physical_radius'] == 'R=epsilon*rho'
        and seam_report['actual_native_collar_transport']['domain']['ordinary_Z_independent_endpoints'] is True
        and patch['Pstar_and_radius_Z_independent'] is True
        and pn[30361] == dict(operation='function_substitution', expression=25202,
            variable=25199, value=29, Z_independent_substitution=True))

    # A fresh universal source lemma. Ordinary jets are derivatives, not
    # Taylor coefficients; U_ZZ and m_ZZ are required for the Z1 row.
    g = FreshGraph()
    one, zero = g.constant(1), g.constant(0)
    z = g.node('fresh_bound_ordinary_Z', domain=['0', '1/2'])
    delta = g.node('common_original_delta', Z_independent=True)
    atom = lambda name: g.node('fresh_common_source_field_row', label=name,
        bound_or_midpoint=False, leading_history_not_reset=True)
    U = tuple(atom('U_Z' + str(k)) for k in range(3))
    M = tuple(atom('physical_own_m_Z' + str(k)) for k in range(3))
    Uy = tuple(atom('physical_Uy_Z' + str(k)) for k in range(2))
    Z, D, O = (z, one), (delta, zero), (one, zero)
    add = lambda a, b: tuple(g.add(a[k], b[k]) for k in range(2))
    sub = lambda a, b: tuple(g.sub(a[k], b[k]) for k in range(2))
    mul = lambda a, b: (g.mul(a[0], b[0]), g.add(g.mul(a[1], b[0]), g.mul(a[0], b[1])))
    scalej = lambda value, a: tuple(g.mul(g.constant(value), v) for v in a)
    d, L = sub(O, mul(Z, Z)), sub(O, mul(D, mul(Z, Z)))
    numerator = lambda V, m, mz: sub(sub(scalej(2, mul(Z, V)), mul(mul(sub(O, D), Z), m)), mul(d, mz))
    Nq = numerator(U[:2], M[:2], M[1:3])
    My, Myz = sub(U[:2], M[:2]), sub(U[1:3], M[1:3])
    Ny = numerator(Uy, My, Myz)
    Nplus = sub(add(mul(mul(add(O, D), Z), U[:2]), scalej(2, mul(Z, Uy))), mul(d, U[1:3]))
    witnesses = [g.sub(Nplus[k], g.add(Nq[k], Ny[k])) for k in range(2)]
    check('new_universal_full_vs_reduced_Qy_numerator_C1',
        all(not g.polynomial(ref) for ref in witnesses))
    R = atom('same_physical_Rminus_Z_independent')
    Qy = (g.node('fresh_positive_quotient', (Ny[0], L[0]), source_positive_L=True),
          g.node('fresh_positive_quotient',
            (g.sub(g.mul(Ny[1], L[0]), g.mul(Ny[0], L[1])), g.mul(L[0], L[0])),
            source_positive_L_squared=True))
    G = tuple(g.mul(R, row) for row in Qy)
    parent_graph = parent['definition_graph']
    check('accepted_parent_graph_immutable', digest(parent_graph['nodes']) == parent_graph['nodes_digest'])
    outputs = []
    for k in range(2):
        left = g.node('actual_native_initial_physical_flux_port', source_archive=TRUE,
            source_node=init['R_Q_y'][k], ordinary_Z_order=k, bound_or_midpoint=False)
        right = g.node('actual_corrected_original_lower_physical_flux_port', source_report=CORR,
            source_node=(53450, 53462)[k], saved_patch_node=(31228, 31240)[k],
            ordinary_Z_order=k, patch_applied_again=False, bound_or_midpoint=False)
        cert = g.node('qualified_lower_physical_C1_source_congruence', (left, right, G[k]),
            admitted=False, ordinary_Z_order=k, source_owned_Z_cell=['0', '1/2'],
            exact_source_predicates_not_sample_agreement=True,
            same_mean_ODE='m_y=U-m', numerical_packets_not_promoted=True)
        prior = parent_graph['outputs'][k]
        terminal = g.node('accepted_parent_terminal_function_port', parent_report=PARENT,
            parent_graph_digest=parent_graph['nodes_digest'],
            parent_node=prior['terminal_with_Rc_flux_correlation'], ordinary_Z_order=k)
        outputs.append(dict(ordinary_Z_order=k, native_initial_RQy=left,
            corrected_original_lower_RQy=right, common_lower_function=G[k],
            exact_source_congruence_certificate=cert,
            unchanged_terminal_with_Rc_flux_correlation=terminal,
            lower_flux_difference_exact='0 on source-owned [0,1/2]',
            active_Rc_compact_inlet_flux_difference_retained=True,
            cumulative_m15_equality_claimed=False))
    passed = all(checks.values())
    for node in g.nodes:
        if 'admitted' in node:
            node['admitted'] = passed
    result = dict(schema=SCHEMA, passed=passed, source_only=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, source_family=family,
        source_prerequisite_audit=checks,
        definition_graph=dict(namespace=SCHEMA, nodes=g.nodes, node_count=len(g.nodes),
            nodes_digest=digest(g.nodes), outputs=outputs, proof_graph_not_numeric_evaluator=True,
            archived_graphs_rebuilt=False),
        source_interface=dict(native_initial_RQy=[129567, 129568],
            original_uncorrected_lower_RQy=[30442, 30443], corrected_patch_RQy=[31228, 31240],
            corrected_final_RQy=[53450, 53462], physical_native_fields=f['physical'],
            common_function='G=Rminus*(((1+delta)*Z*U+2*Z*Uy-(1-Z^2)*UZ)/(1-delta*Z^2)-Q)',
            ordinary_Z1='GZ=Rminus*(Nplus_Z/L+2*delta*Z*Nplus/L^2-QZ)',
            native_formula='Qy=numerator(Uy,U-m)/L with the shifted ordinary mZ rows',
            source_owned_Z_cell=['0', '1/2'], source_cell_proof=proof['incoming_history']['exact_Z_cell'],
            old_quantitative_packet_selected_Z_interval=seam['selected_Z_interval'],
            old_quantitative_packet_whole_cell_bound=False,
            source_scope_not_obtained_from_numeric_packet=True,
            universal_numerator_C1_zero_witnesses=witnesses,
            Pstar_consumed_once=True, lower_y_is_physical=True, mass_ODE='m_y=U-m',
            retained_own_mean_pressure_and_nonzero_incoming_histories=True,
            required_alias_keys=[], required_identity_keys=[],
            direct_graph_node_equality=False,
            native_Rc_compact_join_report='native_Rc_2Rc_pressure_join_lift.json',
            old_active_upper_vs_compact_lower_join_proved=False),
        retained_construction_data=parent['retained_construction_data'],
        unresolved_gates=[key for key, value in checks.items() if not value],
        remaining_construction_gates=['actual corrected active upper Rc flux versus original compact inlet flux source identity',
            'active density and pressure cocycle with retained inlet',
            'actual compact density integral identity', 'old212 same-source packet error theorem',
            'global five-moment closure'],
        legacy_numeric_chain_reconfigured=False, new_unentered_expression_simulations=0,
        runtime_model_metadata=None,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_operators_replayed=False, numeric_values_computed=False,
            archived_definition_graphs_rebuilt=False,
            new_universal_C1_polynomial_only=True))
    result.update({flag: passed for flag in STAGE_FLAGS})
    result.update({flag: False for flag in FALSE_FLAGS})
    return result
