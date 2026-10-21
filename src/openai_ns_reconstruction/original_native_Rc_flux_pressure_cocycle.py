"""Actual Rc flux and pressure cocycle source binding.

Archived source graphs stay opaque. Only a new affine source contraction is
normalized. This proof does not select values from old numerical enclosures.
"""
from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_lower_physical_flux_identity.json'
UPPER = PREFIX + 'native_compact_2Rc_flux_identity.json'
POLAR = PREFIX + 'native_compact_pressure_polarization.json'
JOIN = PREFIX + 'native_Rc_2Rc_pressure_join_lift.json'
SCALE = PREFIX + 'native_physical_scale_source_identity.json'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
NATIVE = PREFIX + 'true_outer_C6_C5_functions.json.gz'
NC = PREFIX + 'true_outer_C6_C5_functions_check.json'
NS = PREFIX + 'true_outer_C6_C5_functions.py'
BAND = PREFIX + 'C3_repair_band.json.gz'
BC = PREFIX + 'C3_repair_band_check.json'
BS = PREFIX + 'C3_repair_band.py'
B2S = PREFIX + 'C2_repair_band.py'
BRIDGE = PREFIX + 'outer_leading_C3_bridge.json.gz'
BRC = PREFIX + 'outer_leading_C3_bridge_check.json'
BRS = PREFIX + 'outer_leading_C3_bridge.py'
CORRECTION = PREFIX + 'outer_to_repair_C3_bridge.json.gz'
CC = PREFIX + 'outer_to_repair_C3_bridge_check.json'
CCS = PREFIX + 'outer_to_repair_C3_bridge.py'
CORR = PREFIX + 'corrected_reserved_m15_source.json'
PRES = PREFIX + 'n1_preRp_moment_functions.py'
END = PREFIX + 'n1_preRc_endpoint_C1.json.gz'
EC = PREFIX + 'n1_preRc_endpoint_C1_check.json'
ES = PREFIX + 'n1_preRc_endpoint_C1.py'
PHYS = PREFIX + 'n1_all_chart_physical_pressure_fields.py'
PART = PREFIX + 'n1_partial_mean_Q.py'
C0S = PREFIX + 'n1_preRc_pressure_endpoint.py'
QUIETCHECK = PREFIX + 'whole_Z_all_N_outer_Rc_functions_check.json'
SCHEMA = 'original_native_Rc_flux_pressure_cocycle_v1'
RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_native_lower_physical_flux_identity.json': 'f1f98828579a76fc2e71f812b7da0f356620f73e5c707bb91edb1bf789792035', 'lei_ren_part1_paper_compliant_current_original_native_compact_2Rc_flux_identity.json': 'f5a7ccd09b895d89e360d6b94f9a71b091e4c6002246d86bbd2ad288ec15b717', 'lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json': 'be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610', 'lei_ren_part1_paper_compliant_current_original_native_Rc_2Rc_pressure_join_lift.json': '1f4b6aedf6683299117fb3bc2cddb90241a216f2b83296dddba04f60b53df39f', 'lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json': 'b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz': '75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.json': 'ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py': '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': 'cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band_check.json': '00149597ab6206348a62dc39037ff53b7c5bcfec5d7eec6ae2591444f0f91a23', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.py': '62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4', 'lei_ren_part1_paper_compliant_current_original_C2_repair_band.py': 'e0de87544765904c065eb8edbe8fc3c71567329a2b7e49c6cd1205d428f88226', 'lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge.json.gz': '86bcc0c06ba4958081875d9f1260f078b745cd7c2d1071f27dfdc586c3a6ef76', 'lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge_check.json': 'dfdfa501b5d6ffdafdb141bb858bf1783597215a9ffafd3cdd470b5fe68e9eb6', 'lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge.py': '6f23b6b3e5d3bf451f03e23d33a53e28525f2f29a1dbeafe7add1725b94b03cd', 'lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge.json.gz': '7d1984552d6fce11ef9fbd86eb76997311bb4096eb5e4fd43047b5940740b36e', 'lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge_check.json': 'a475fe2e8d927645b799a06e72629ac0e0179da24ff00e330ad8a86418f1dea1', 'lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge.py': '6b7922662e42956e790f8fa970070ae7a570662151ed4ee64dfa73e10df3d353', 'lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json': 'fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py': '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.json.gz': 'e691db8090b4400537228db14d387a46e3163e73c94880a319a04f576cf8ec78', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1_check.json': '02dec83e6e5b4c52d8aadc7e5f8dc359b2322dc62e5aa3ce895f60e8f855ecb2', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.py': '828afe7bd79f91d8ccbeb989711147306b41abe400344b4b9fa3f55aaa14eb37', 'lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields.py': 'df9d282e3e1b5c7c74d5101c6e37046b6e28138d40b08ca05024dbe790c59f28', 'lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q.py': 'd826c150d021fa318201ce9b393b8efb92b459546e4ce8377c2e852f2d65cbc6', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_pressure_endpoint.py': 'aa36f0d9b59071f2d32c75602a9c32475a1835a997a70b1f2a61858a4ccf8a20', 'lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json': 'f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158'}
DECODED_RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz': '53b1bbf15c8c2c5a8400d7e22ef940c6a502397df8cf26bb23e51ff3869e9139', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': '845bade42d86d80e0ef82964e449c9167b18fc26c10f3849ae88768204f8698a', 'lei_ren_part1_paper_compliant_current_original_outer_leading_C3_bridge.json.gz': '971bccc00d5d8574d71bd3d9ef55abd2e19489a92990b9c6170fa922bc06f744', 'lei_ren_part1_paper_compliant_current_original_outer_to_repair_C3_bridge.json.gz': '7ecffa99a071aaa58e93c77e610fd0a0029934a94f0ac52dc1f66f042054a3df', 'lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1.json.gz': '5b21d170f1cfcc681f0895d745e60d38bbfe428dae161111d752925b4f9877cb'}
STAGE_FLAGS = ('actual_original_Rc_active_compact_flux_C1_identity_proved', 'actual_native_Rc_cumulative_pressure_cocycle_bound', 'actual_corrected_active_pressure_cocycle_bound', 'internal_Rc_and_lower_fluxes_cancelled_by_source_identities', 'signed_terminal_reduced_to_actual_density_integral_difference')
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

    def check(name, condition):
        checks[name] = bool(condition)

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
    upper = data[UPPER]['component_result']
    polar = data[POLAR]['component_result']
    joins = data[JOIN]['component_result']
    pre, native, band, bridge = (data[name] for name in (PRE, NATIVE, BAND, BRIDGE))
    pn, nn = pre['exact_graph_nodes'], native['exact_graph_nodes']
    bn = band['exact_graph_nodes']
    nf = native['true_outer_C6_C5_bindings']
    ob = band['actual_C3_repair_band_functions']
    correction = data[CORRECTION]['actual_correction_interface_function_identity']
    cc, bc, brc = (data[name] for name in (CC, BC, BRC))
    corr = data[CORR]
    patch = corr['corrected_terminal_source_inputs']['saved_patch']
    construction = corr['corrected_reserved_m15_source']
    family = parent['source_family']
    check('accepted_parent_lower_physical_C1_identity_on_owned_cell', parent['passed'] is True
        and parent['lower_physical_RQy_source_C1_identity_proved_on_owned_Z_cell'] is True
        and parent['source_interface']['source_owned_Z_cell'] == ['0', '1/2'])
    check('accepted_upper_and_actual_compact_density_representation', upper['passed'] is True
        and upper['common_K_2Rc_flux_difference_C1_zero_proved'] is True
        and polar['passed'] is True
        and polar['exact_compact_m15_defect_representation_proved'] is True
        and polar['common_physical_R_squared_log_measure_proved'] is True)
    check('actual_common_source_family_and_graph_prefix',
        all(row['source_family'] == family for row in (pre, native, band, bridge))
        and pn[:len(bn)] == bn and nn[:len(pn)] == pn)
    check('actual_native_pressure_geometry_cocycle_source', data[NC]['all_passed'] is True
        and syntax(trees[NS], 'pressure_geometry',
            "cumulative['S1']", "initial['R_Q_y']", "boundary['R_Q_y']",
            "m15=a.sub(a.scale(cumulative['S1'], '1/2')")
        and syntax(trees[PRES], 'pressure_functions',
            "boundary = lambda name: a.sub(upper[name], lower[name])",
            "m15=a.sub(half(S1), boundary('R_q_y'))"))
    check('actual_native_Rc_endpoint_pressure_integral_and_flux_ports',
        nf['endpoints']['O3_power']['integrals']['S1'][:2] == [258837, 258838]
        and nf['endpoints']['O3_power']['pressure']['m15'][:2] == [258845, 258846]
        and nf['endpoints']['O3_power']['boundary']['R_Q_y'][:2] == [258853, 258854]
        and nf['initial']['R_Q_y'][:2] == [129567, 129568])
    # Check the saved constructor's literal arithmetic, not its value.
    check('native_cumulative_S1_half_and_initial_boundary_literal_rows',
        nn[258614]['operation'] == 'sum' and nn[258615]['operation'] == 'sum'
        and nn[258837] == dict(operation='function_substitution', expression=258614,
            variable=250215, value=48, Z_independent_substitution=True)
        and nn[258838] == dict(operation='function_substitution', expression=258615,
            variable=250215, value=48, Z_independent_substitution=True))
    check('corrected_active_source_pressure_ports_and_raw_reduced_S1',
        pre['signed_pieces']['active']['m15'] == [30707, 30708]
        and pn[30699] == dict(operation='product', arguments=[28, 30687])
        and pn[30700] == dict(operation='product', arguments=[28, 30688])
        and patch['corrected_pairs']['active_m15'] == [31244, 31248]
        and construction['mapping']['31244'] == 53466
        and construction['mapping']['31248'] == 53470)
    check('lower_patch_does_not_change_actual_Rc_upper_flux_or_density',
        patch['old_graph_unchanged'] is True
        and patch['exact_replacements'] == {'30365': 30363, '30366': 30364}
        and not any(str(node) in patch['changed_original_to_appended']
            for node in (30479, 30480, 30687, 30688)))
    endpoint = data[END]['samples']["('0', '.5')"]['upper']['proof']
    check('actual_active_power_exit_flat_physical_jet_source_theorem',
        data[EC]['all_passed'] is True and endpoint['source_family'] == family
        and all(endpoint[key] is True for key in ('physical_U_Uy_UZ_UZZ_UyZ_exact_zero',
            'q_y_equals_minus_q', 'q_y_Z_equals_minus_q_Z', 'finite_N_own_mean_q_cap_retained'))
        and endpoint['physical_mean_ode'] == 'm_y=U-m=-m'
        and syntax(trees[C0S], '_upper_caps', "_flat_driver_proof(fields, ends, 'O3_power')")
        and syntax(trees[PRES], 'active_functions', 'qy = algebra.neg(q)'))
    quiet_source = pn[30140]['original_leading_source_recipe']['quiet_source_proof']
    check('actual_power_loop_primitive_flat_not_physical_Q_zero',
        data[QUIETCHECK]['all_passed'] is True
        and quiet_source['checked_receipt_sha256'] == RECORD_HASHES[QUIETCHECK]
        and quiet_source['source_identity'] == '2*mu>=eta implies q=q_Z=0 throughout current power offset [0,2]'
        and quiet_source['original_primitive_A_B_zero_but_incoming_histories_retained'] is True)
    profile = ob['profiles_at_band_x']
    check('actual_compact_inlet_flat_profiles_and_complete_mass_ODE',
        bc['all_passed'] is True
        and bc['accepted_C2_aliases']['existing_second_functions_and_independent_P0_exactly_retained'] is True
        and profile['compact_support_all_endpoint_jets_zero'] is True
        and ob['actual_reserved_band'] == [1, 2]
        and syntax(trees[B2S], 'profiles_C2', "lazy_outside_support=True",
            "compact_support_all_endpoint_jets_zero=True")
        and syntax(trees[BS], 'build', "m=alg.add(alg.scale(complete['m'], -1), V)"))
    rc_interface = brc['independent_complete_Rc_interface']
    check('accepted_actual_complete_Rc_to_C3_compact_inlet_source_theorem',
        brc['all_passed'] is True
        and brc['current_outer_complete_history_to_repair_inlet_function_identity_installed'] is True
        and brc['current_outer_leading_endpoint_band_seed_function_identity_installed'] is True
        and rc_interface['independent_actual_canonical_Rc_seed_function_rows'] == 20
        and rc_interface['independent_complete_Rc_to_actual_band_inlet_rows'] == 20
        and rc_interface['same_original_N_divided_once'] is True
        and rc_interface['zero_length_local_integrals_with_nonzero_incoming_kept'] is True
        and rc_interface['source_geometry_Rw_exp2_equals_Rc_and_actual_band_x1'] is True
        and bridge['resolved_function_graph_contract']['same_original_source_family_and_common_V_m_k_Pstar_normalization'] is True)
    mass_rows = [row for row in correction['rows'] if row['history'] == 'm']
    check('accepted_exact_actual_nonzero_correction_inlet_not_a_bound',
        cc['all_passed'] is True
        and correction['outer_N_scaled_equals_N_times_actual_repair_inlet'] is True
        and correction['graph_node_numbering_and_range_overlap_not_used_as_identity'] is True
        and len(mass_rows) == 4
        and all(row['same_defining_function'] is True for row in mass_rows)
        and [row['outer_N_scaled_node'] for row in mass_rows[:3]] == [4920, 4921, 12042]
        and [row['actual_repair_inlet_node'] for row in mass_rows[:3]] == [4972, 4973, 13578])
    common_recipe_keys = ('module', 'module_sha256', 'function', 'AST_sha256',
                          'defining_provider', 'native_chart')
    check('active_leading_mass_same_actual_provider_path_and_ordinary_rows',
        all(pn[a]['operation'] == 'current_original_physical_source_quantity'
            and pn[a]['quantity'] == 'm' and pn[a]['Z_order'] == k
            and pn[a]['Taylor_coefficient_factorial'] == (1, 1, 2)[k]
            and pn[a]['defining_quantity_not_a_range_value'] is True
            and pn[a]['raw_quantity_path'] == 'original_generic_source.common_own_five_histories_axial5["m"][' + str(k) + ']'
            and pn[b]['quantity'] == 'history_m' and pn[b]['coordinate'] == 48
            and pn[b].get('Taylor_coefficient_factorial', 1) == (1, 1, 2)[k]
            and all(pn[a]['original_leading_source_recipe'][key] == pn[b]['recipe'][key]
                for key in common_recipe_keys)
            for k, (a, b) in enumerate(((30140, 13087), (30150, 13088), (30160, 13608))))
        and syntax(trees[PHYS], '_field_functions',
            "delta_m = owner.parent.functions[chart]['mean']",
            'm = alg.mul(scale, alg.add(baseline_m, delta_m))',
            'q = partial.recover_Q_increment(g, U, m, owner.delta_node)')
        and syntax(trees[PART], 'partial_mean_functions', 'incoming_mean_not_reset=True',
            'mean = alg.mul(eps, c3.C3Function(*scaled))'))
    # Actual power-chart mass memory is exp(-(y-y_left))*incoming/N.
    # The free-query and complete terminal forms use the same saved incoming
    # ports. No old Duhamel integral, signature reducer or normalizer runs.
    check('actual_active_power_mass_memory_terminal_and_original_N_ports',
        pn[19]['numerator'] == pn[25225]['numerator'] == 1
        and pn[19]['denominator'] == pn[25225]['denominator'] == 18
        and pn[30127]['body'] == 4891 and pn[30127]['native_variable'] == 'native_O3_power'
        and pn[30127]['query_coordinate'] == 30126
        and pn[30142] == dict(operation='sum', arguments=[30127, 4894])
        and pn[30143] == dict(operation='negative', argument=30142)
        and pn[30144] == dict(operation='analytic_unary', name='exp', argument=30143)
        and pn[4918] == dict(operation='negative', argument=4895)
        and pn[4919] == dict(operation='analytic_unary', name='exp', argument=4918)
        and pn[4895] == dict(operation='sum', arguments=[4893, 4894])
        and all(pn[a] == dict(operation='product', arguments=[30144, incoming])
            and pn[b] == dict(operation='product', arguments=[4919, incoming])
            for a, b, incoming in ((30145, 4920, 4820), (30152, 4921, 4821), (30162, 12042, 12010))))
    check('actual_active_complete_physical_mass_rows_scale_and_sum_once',
        all(pn[scaled_delta] == dict(operation='product', arguments=[25225, memory])
            and pn[complete] == dict(operation='sum', arguments=[leading, scaled_delta])
            and pn[physical] == dict(operation='product', arguments=[25203, complete])
            for leading, memory, scaled_delta, complete, physical in (
                (30141, 30145, 30146, 30147, 30148),
                (30151, 30152, 30153, 30154, 30155),
                (30161, 30162, 30163, 30164, 30165))))
    check('actual_power_coordinate_and_absolute_Rc_alias_recipes',
        pn[4584] == dict(operation='exact_rational', numerator=1012, denominator=1)
        and pn[4892] == dict(operation='exact_rational', numerator=1014, denominator=1)
        and pn[48] == dict(operation='exact_rational', numerator=2, denominator=1)
        and pn[4891]['arguments'] == [4584, 4890, 2084, 2091, 649, 3, 645, 1507]
        and pn[4893]['arguments'] == [4892, 2084, 2091, 649, 3, 645, 1507]
        and pn[30128] == dict(operation='sum', arguments=[24955, 30127])
        and pn[24956] == dict(operation='sum', arguments=[4893, 24955])
        and pn[30828] == dict(operation='analytic_unary', name='exp', argument=24956))
    check('actual_active_and_compact_Rc_flux_definition_ports',
        pre['defining_bindings']['active']['finite_endpoints'][1]['R_q_y'] == [30479, 30480]
        and pre['defining_bindings']['compact']['finite_endpoints'][0]['R_q_y'] == [30903, 30904]
        and pn[30479] == dict(operation='product', arguments=[30448, 30450])
        and pn[30480] == dict(operation='product', arguments=[30448, 30451])
        and pn[30450] == dict(operation='product', arguments=[220, 30446])
        and pn[30451] == dict(operation='product', arguments=[220, 30447])
        and syntax(trees[PRES], 'compact_functions',
            "f['complete_histories']['m']", "f['complete_history_y']['m']",
            'recover_Q_increment(g, Uy, my, delta)', 'for edge in (g.one, g.constant(2))'))
    check('accepted_native_Rc_physical_flux_join_including_ordinary_Z1',
        joins['passed'] is True
        and all(any(row['seam'] == 'Rc' and row['moment'] == 'm15'
            and row['ordinary_Z_order'] == k and row['actual_left_handle'] == (258853, 258854)[k]
            and row['actual_right_handle'] == (263155, 263156)[k] and row['admitted'] is True
            for row in joins['actual_jump_zero_rows']) for k in range(2)))

    g = FreshGraph()
    one, half = g.constant(1), g.constant('1/2')
    coordinate_witness = g.sub(g.add(g.constant(pn[4584]['numerator']),
        g.constant(pn[48]['numerator'])), g.constant(pn[4892]['numerator']))
    check('new_actual_power_query_right_endpoint_affine_identity', not g.polynomial(coordinate_witness))
    # Reuse the accepted universal full/reduced-Qy numerator theorem as data.
    # At the actual flat inlet it specializes to Qy=-Q for arbitrary own m.
    flat_theorem = g.node('qualified_actual_old_Rc_flux_source_congruence',
        native_chart='O3_power', original_active_RQy=[30479, 30480],
        original_compact_RQy=[30903, 30904], ordinary_Z_orders=[0, 1],
        source_cell=['0', '1/2'], common_Pstar_node=25203,
        source_complete_m_rows=ob['complete_histories']['m'][:3],
        accepted_quotient_C1_theorem_report=PARENT,
        accepted_quotient_C1_theorem_graph=parent['definition_graph']['nodes_digest'],
        accepted_complete_history_bridge=BRIDGE,
        actual_power_endpoint_coordinate_witness=coordinate_witness,
        formula='U=Uy=UZ=0, my=-m imply Qy=-Q; same complete m and physical R imply same RQy',
        nonzero_memory_not_reset=True, m_equals_Zc_assumed=False,
        direct_graph_handle_equality=False, admitted=False)
    parent_graph, compact_graph = parent['definition_graph'], upper['definition_graph']
    check('accepted_parent_and_compact_graph_hashes',
        digest(parent_graph['nodes']) == parent_graph['nodes_digest']
        and digest(compact_graph['nodes']) == compact_graph['nodes_digest'])
    outputs = []
    for k in range(2):
        def port(label, archive, node):
            return g.node('actual_source_function_port', label=label,
                source_record=archive, source_node=node, ordinary_Z_order=k,
                bound_or_midpoint=False)

        IN = port('native_Rc_cumulative_S1', NATIVE, (258837, 258838)[k])
        IO = port('original_active_reduced_S1', PRE, (30687, 30688)[k])
        IC = g.node('accepted_compact_density_difference_integral_port',
            source_record=UPPER, source_graph_digest=compact_graph['nodes_digest'],
            source_node=compact_graph['outputs'][k]['weighted_integral'],
            ordinary_Z_order=k, measure='Rc^2*x dx', bound_or_midpoint=False)
        GN = port('native_Rc_RQy', NATIVE, (258853, 258854)[k])
        GO = port('original_active_Rc_RQy', PRE, (30479, 30480)[k])
        GC = port('original_compact_inlet_RQy', PRE, (30903, 30904)[k])
        gN = port('native_initial_RQy', NATIVE, (129567, 129568)[k])
        gO = port('corrected_original_lower_RQy', CORR, (53450, 53462)[k])
        mn, mo = (port(label, archive, node) for label, archive, node in (
            ('native_Rc_m15', NATIVE, (258845, 258846)[k]),
            ('corrected_active_m15', CORR, (53466, 53470)[k])))
        n_formula = g.add(g.sub(g.mul(half, IN), GN), gN)
        o_formula = g.add(g.sub(g.mul(half, IO), GO), gO)
        n_cert = g.node('qualified_actual_native_pressure_cocycle', (mn, n_formula),
            extra_initial_m15=False, ordinary_Z_order=k, admitted=False)
        o_cert = g.node('qualified_actual_corrected_active_pressure_cocycle', (mo, o_formula),
            actual_saved_patch_consumed_once=True, ordinary_Z_order=k, admitted=False)
        full = g.add(g.sub(n_formula, o_formula), g.mul(half, IC), g.sub(GN, GC))
        density = g.mul(half, g.add(g.sub(IN, IO), IC))
        seam_defect = g.sub(GO, GC)
        lower_defect = g.sub(gN, gO)
        witness = g.sub(g.sub(full, density), g.add(seam_defect, lower_defect))
        check('new_actual_signed_pressure_flux_contraction_Z' + str(k), not g.polynomial(witness))
        terminal = g.node('accepted_parent_terminal_source_function_port',
            source_record=PARENT, source_graph_digest=parent_graph['nodes_digest'],
            source_node=parent_graph['outputs'][k]['unchanged_terminal_with_Rc_flux_correlation'],
            ordinary_Z_order=k)
        certificate = g.node('qualified_actual_terminal_density_difference_identity',
            (terminal, density, n_cert, o_cert, flat_theorem),
            accepted_lower_flux_identity_report=PARENT, ordinary_Z_order=k,
            source_owned_Z_cell=['0', '1/2'], inherited_numeric_errors_deleted=False,
            density_difference_zero_claimed=False, admitted=False)
        outputs.append(dict(ordinary_Z_order=k, native_cumulative_S1=IN,
            old_active_reduced_S1=IO, compact_density_difference_integral=IC,
            terminal_density_difference=density, density_source_identity_certificate=certificate,
            affine_zero_witness=witness, cancelled_original_Rc_flux_difference=seam_defect,
            cancelled_lower_flux_difference=lower_defect,
            exact_identity='F-T=1/2*(I_N_Rc-I_O_active+Rc^2*integral_1^2 x*DeltaS dx)',
            density_difference_zero_claimed=False))
    passed = all(checks.values())
    for row in g.nodes:
        if 'admitted' in row:
            row['admitted'] = passed
    result = dict(schema=SCHEMA, passed=passed, source_only=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, source_family=family,
        source_prerequisite_audit=checks,
        definition_graph=dict(namespace=SCHEMA, nodes=g.nodes, node_count=len(g.nodes),
            nodes_digest=digest(g.nodes), outputs=outputs, proof_graph_not_numeric_evaluator=True,
            archived_graphs_rebuilt=False),
        source_interface=dict(original_active_Rc_RQy=[30479, 30480],
            original_compact_inlet_RQy=[30903, 30904], native_Rc_RQy=[258853, 258854],
            native_compact_inlet_RQy=[263155, 263156], native_S1_Rc=[258837, 258838],
            original_active_S1=[30687, 30688], native_Rc_m15=[258845, 258846],
            corrected_active_m15=[53466, 53470], source_owned_Z_cell=['0', '1/2'],
            pressure_cocycle='m15=S1/2-upper_RQy+lower_RQy', initial_m15_added=False,
            mass_ODE='m_y=U-m', original_Rc_mass_arbitrary_C3_not_linear_Z=True,
            physical_measure='S1 integrates R^2 S dy; compact difference uses Rc^2*x dx',
            new_actual_coordinate_zero_witness=coordinate_witness,
            original_Rc_source_identity_certificate=flat_theorem,
            complete_mass_leading_source_nodes=[13087, 13088, 13608],
            complete_mass_correction_inlet_nodes=[4972, 4973, 13578],
            ordinary_mass_rows_needed_for_flux_Z1=[0, 1, 2],
            old_quantitative_packet_whole_cell_bound=False,
            retained_nonzero_complete_histories=True, old_signed_terms_or_errors_rewritten=False,
            required_alias_keys=[], required_identity_keys=[]),
        retained_construction_data=parent['retained_construction_data'],
        unresolved_gates=[key for key, value in checks.items() if not value],
        remaining_construction_gates=['actual native versus original active density source identity',
            'actual compact incoming and control density function identity',
            'same-source old212 numerical packet error theorem', 'global five-moment closure'],
        legacy_numeric_chain_reconfigured=False, new_unentered_expression_simulations=0,
        runtime_model_metadata=None,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_operators_replayed=False, numeric_values_computed=False,
            archived_definition_graphs_rebuilt=False, archived_quotient_lemma_reexecuted=False,
            new_actual_endpoint_affine_and_terminal_contraction_only=True))
    result.update({key: passed for key in STAGE_FLAGS})
    result.update({key: False for key in FALSE_FLAGS})
    return result
