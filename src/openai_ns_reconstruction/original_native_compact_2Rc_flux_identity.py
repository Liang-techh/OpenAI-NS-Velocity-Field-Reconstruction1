"""Exact actual compact 2Rc flux cancellation; archived functions stay opaque.

Only new elementary polynomial witnesses and a fresh affine contraction run.
The source theorem concerns exact functions, not inherited numerical errors.
"""
from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_compact_pressure_polarization.json'
SCALE = PREFIX + 'native_physical_scale_source_identity.json'
NATIVE = PREFIX + 'true_compact_C5_functions.json.gz'
OLD = PREFIX + 'C3_repair_band.json.gz'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
POWER = PREFIX + 'C3_power_to_Rp.json.gz'
FRAME = PREFIX + 'C3_Rp_native_frame.json.gz'
NC = PREFIX + 'true_compact_C5_functions_check.json'
QC = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions_check.json'
PC = PREFIX + 'C3_power_to_Rp_check.json'
FC = PREFIX + 'C3_Rp_native_frame_check.json'
JOIN = PREFIX + 'native_Rc_2Rc_pressure_join_lift.json'
NS = PREFIX + 'true_compact_C5_functions.py'
QS = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.py'
PS = PREFIX + 'C3_power_to_Rp.py'
PRES = PREFIX + 'n1_preRp_moment_functions.py'
PHYS = PREFIX + 'modified_Q4_pressure_C3_functions.py'
OS = PREFIX + 'C3_repair_band.py'
SCHEMA = 'original_native_compact_2Rc_flux_identity_v1'
RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json': 'be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610', 'lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json': 'b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz': '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': 'cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz': 'bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4', 'lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.json.gz': 'cb2a69ca971965a6e487bd0fa818204216747094f37819a2cc3e42a3af79911a', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.json.gz': '859b59d422104d184092bcf60f8afa95267335ce87e8b41f2c271c18130e9a28', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json': '6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.json': '3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0', 'lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp_check.json': '49eb13920a5a57acfecb237b9e606ce80afb58726be8423273940fbfc1bbf362', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame_check.json': 'd34d0e4f90b131e68bbde16caad33029ff89b75174741ad02d03d0a6e135d9e8', 'lei_ren_part1_paper_compliant_current_original_native_Rc_2Rc_pressure_join_lift.json': '1f4b6aedf6683299117fb3bc2cddb90241a216f2b83296dddba04f60b53df39f', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py': 'f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.py': '0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068', 'lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.py': 'f1630fb589f9e5dc74d2f4c5ac3bff4f3e396cc91be141c7dfe6fce5bd95ec94', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py': '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.py': '62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4'}
DECODED_RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz': 'b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': '845bade42d86d80e0ef82964e449c9167b18fc26c10f3849ae88768204f8698a', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz': '3bc0f9984a0aeb6f447d177c40c25116b8a59c3b90ee09b8661f029cbe032a40', 'lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.json.gz': 'cb768abbc4c71e1122250c3c37e7d43e003f1e1457f617802d3b225cf58e0b01', 'lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.json.gz': '1d563553e3c7e6f9115a024a52e61f3b308d822104aaf557f3d235c54f26738d'}
STAGE_FLAGS = ('actual_native_compact_2Rc_flux_equals_native_K_C1_proved', 'actual_original_compact_2Rc_flux_equals_original_K_C1_proved', 'common_K_2Rc_flux_difference_C1_zero_proved', 'compact_pressure_upper_endpoint_removed_by_source_identity', 'signed_terminal_with_Rc_flux_correlation_defined')
FALSE_FLAGS = ('native_compact_m15_source_identity_proved', 'native_compact_V_m_source_identity_proved', 'native_active_m15_source_identity_proved', 'terminal_difference_zero_proved', 'old_error_terms_removed', 'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1', 'global_five_moment_closure', 'legacy_feedback_reference_alias_proved', 'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved', 'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed', 'full_goal_complete', 'old212_adapter_source_identity_proved', 'independent_terminal_function_match_complete', 'native_vs_registered_preRp_m15_identity_proved')

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def syntax(tree, function, *needles):
    bodies = [ast.unparse(row) for row in ast.walk(tree)
        if isinstance(row, ast.FunctionDef) and row.name == function]
    return any(all(needle in body for needle in needles) for body in bodies)


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
                        key = tuple(sorted(a+b))
                        product[key] = product.get(key, Fraction(0)) + av*bv
                out = product
            return {key: value for key, value in out.items() if value}
        return {(ref,): Fraction(1)}


def derive(source_records, *, expected_input_digest):
    checks, data, trees = {}, {}, {}

    def check(key, condition):
        checks[key] = bool(condition)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    for name, pin in RECORD_HASHES.items():
        row = source_records[name]
        check('authenticated:' + name, row['raw_sha256'] == pin
            and row['raw_authenticated_before_decode'] is True
            and row['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, pin))
        if name.endswith('.py'):
            trees[name] = ast.parse(row['source_text'])
        else:
            data[name] = json.loads(row['source_text'])
    parent = data[PARENT]['component_result']
    scale = data[SCALE]['component_result']
    native, old, pre, quiet, power, frame = (data[k] for k in (NATIVE, OLD, PRE, QUIET, POWER, FRAME))
    n = native['true_compact_C5_bindings']
    o = old['actual_C3_repair_band_functions']
    p = pre['defining_bindings']['compact']
    q = quiet['true_quiet_Rp_C5_C3_bindings']
    oq = pre['defining_bindings']['quiet']
    pw = power['actual_C3_repaired_exit_power_and_pulse_functions']
    fr = frame['actual_Rp_native_frame_functions']
    ng, og, pg, qg = (v['exact_graph_nodes'] for v in (native, old, pre, quiet))
    family = parent['source_family']
    check('same_source_family', all(v['source_family'] == family
        for v in (native, old, pre, quiet, power, frame)))
    check('accepted_parent_actual_compact_representation', parent['passed'] is True
        and parent['exact_compact_density_polarization_C1_proved'] is True
        and parent['exact_compact_m15_defect_representation_proved'] is True)
    check('accepted_common_physical_Rc_Pstar_K_source_aliases', scale['passed'] is True
        and scale['native_original_absolute_Rc_Rp_identity_proved'] is True
        and scale['old_preRp_Pstar_source_identity_proved'] is True
        and scale['native_original_quiet_K_identity_proved'] is True
        and scale['source_aliases']['quietK']['native_node'] == q['K'] == 265180
        and scale['source_aliases']['quietK']['old_kappa_node'] == oq['kappa'] == 31204
        and scale['source_aliases']['Pstar']['provider_node'] == 25203
        and scale['source_aliases']['Pstar']['underlying_formula_node'] == n['Pstar'] == 33572)
    check('archived_graph_prefixes_exact', pg[:len(og)] == og
        and ng[:len(pg)] == pg and qg[:len(ng)] == ng)
    check('actual_old_compact_complete_mass_and_physical_scale',
        p['common_Pstar_applied_once'] is True
        and p['actual_C3_original_band_variable'] == o['original_band_variable'] == 5556
        and syntax(trees[PRES], 'compact_functions', "f['complete_histories']['m']",
            "f['complete_history_y']['m']", "profile[0]['V']", "profile[1]['V']",
            'recover_Q_increment(g, Uy, my, delta)'))
    check('actual_physical_Qy_and_mass_ODE_source',
        syntax(trees[PHYS], 'physical_pressure_fields', 'my = a.sub(U, m)',
            'Qy = div(numerator(Uy, my))')
        and syntax(trees[OS], 'build', "m=alg.add(alg.scale(complete['m'], -1), V)"))
    nc, qc, pc, fc = (data[k] for k in (NC, QC, PC, FC))
    check('native_actual_flat_compact_profile_theorem', nc['all_passed'] is True
        and nc['checks']['actual_profiles']['passed'] is True
        and nc['checks']['actual_profiles']['actual_original_raw_beta_end_flat_rows'] == 18
        and nc['checks']['actual_source']['original_compact_log_transport_measure_and_relative_only_terminal_theorem'] is True)
    theorem = qg[q['linear_mass_theorem']]
    check('native_actual_complete_linear_mass_and_Z_independent_K',
        qc['all_passed'] is True and qc['checks']['native_mass_identity']['passed'] is True
        and qc['checks']['actual_source_geometry']['original_constant_K_source_Z_independent'] is True
        and qc['checks']['actual_source_geometry']['true_compact_exit_and_all_original_nonzero_histories_preserved'] is True
        and theorem['operation'] == 'proved_original_quiet_linear_mass_source_identity'
        and theorem['actual_compact_complete_mass_rows'] == n['terminal']['modified_fields']['m']
        and theorem['compact_relative_zero'] == n['terminal_certificates']['m']
        and theorem['physical_scale'] == n['Pstar']
        and theorem['constant_K'] == q['K']
        and theorem['current_quiet_identity'] == 'physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q'
        and theorem['source_constant_not_a_saved_range_endpoint'] is True
        and qg[q['K']] == dict(operation='product', arguments=[q['R2'], q['kappa2']]))
    joins = data[JOIN]['component_result']
    native_flux = n['terminal']['boundary']['R_Q_y'][:2]
    quiet_flux = q['inlet']['boundary']['R_Q_y'][:2]
    check('actual_native_2Rc_RQy_join_not_only_pressure_cumulative',
        joins['passed'] is True and joins['native_2Rc_pressure_join_C3_proved'] is True
        and native_flux == [262914, 262915] and quiet_flux == [265671, 265672]
        and all(any(row['seam'] == '2Rc' and row['moment'] == 'm15'
            and row['ordinary_Z_order'] == k and row['actual_left_handle'] == native_flux[k]
            and row['actual_right_handle'] == quiet_flux[k] and row['admitted'] is True
            for row in joins['actual_jump_zero_rows']) for k in range(2)))
    check('native_exact_2Rc_radius_substitution', q['R2'] == n['terminal']['radius']
        and ng[q['R2']] == dict(operation='function_substitution',
            expression=n['background']['radius'], variable=n['background']['phase'],
            value=48, Z_independent_substitution=True)
        and ng[48] == dict(operation='exact_rational', numerator=2, denominator=1))
    old_exit = pc['independent_actual_compact_repaired_exit']
    check('old_actual_flat_complete_compact_exit_theorem', pc['all_passed'] is True
        and old_exit['actual_compact_support_nodes'] == 9
        and set(old_exit['actual_beta_arguments_at_2Rc']) == {'8', '20', '32'}
        and old_exit['actual_profile_y0_y1_y2_Z0_Z1_Z2_Z3_flat_exit_rows'] == 24
        and old_exit['same_source_compact_beta_all_higher_radial_derivatives_flat'] is True
        and old_exit['complete_absolute_histories_and_independent_P0_not_zeroed'] is True)
    check('old_complete_mass_power_flow_and_inlet_semigroup',
        syntax(trees[PS], 'power_flow', "m=alg.scale(own['m'], decay(1))", 'V=zero', 'm=zero')
        and syntax(trees[PS], 'build', "for key, row in bf['complete_histories'].items()",
            'actual_flow = power_flow(g, actualE2, complete2, mu, s)',
            'Rp = g.mul(R2, g.unary(\'exp\', length))')
        and pc['independent_quiet_mixed_calculus']['all_five_actual_ODE_and_inlet_identities'] is True
        and pc['independent_quiet_mixed_calculus']['all_five_complete_state_semigroup_identities'] is True
        and pw['compact_repair_at_Rc_to_2Rc_consumed'] is True)
    check('old_complete_Rp_frame_matches_actual_power_exit', fc['all_passed'] is True
        and fc['independent_actual_complete_chain_and_Rp_frame']['genuine_compact_band_and_complete_exit_consumed'] is True
        and fc['independent_actual_complete_chain_and_Rp_frame']['actual_current_quiet_Rp_native_frame_C3_rows'] == 32
        and fc['independent_actual_native_radius']['exact_same_current_quiet_Rc_and_Tw_defining_expressions'] is True
        and oq['constant_Z_identity']['same_original_endpoint_function_nodes']['quiet'] == pw
        and oq['constant_Z_identity']['same_original_endpoint_function_nodes']['native'] == fr)
    check('old_frame_linear_mass_constant_kappa_and_exact_import',
        oq['constant_Z_identity']['exact_Z_independent_kappa'] is True
        and oq['constant_Z_identity']['kappa_Z_derivative'] == '0'
        and oq['constant_Z_identity']['exact_q_common_at_Rp_is_minus_kappa'] is True
        and fr['actual_identified_Rp_histories_C3']['m'] == [26262, 26263, 26264, 26264]
        and oq['original_native_frame_import_mapping']['26263'] == oq['kappa']
        and pg[oq['K']] == dict(operation='product', arguments=[oq['upper'], 25203, oq['kappa']])
        and syntax(trees[PRES], 'quiet_functions', "f['actual_identified_Rp_histories_C3']['m'][1]",
            'K = g.mul(Rp, scale, kappa)'))
    old_flux = p['finite_endpoints'][1]['R_q_y']
    check('old_actual_2Rc_RQy_source_ports', old_flux == [30939, 30940]
        and all(pg[node]['operation'] == 'product' and 30911 in pg[node]['arguments']
            for node in old_flux)
        and pg[30911] == dict(operation='product', arguments=[30828, 48]))

    # Elementary proof on new formal symbols only: Qy at a flat linear-mass exit.
    g = FreshGraph()
    one, zero = g.constant(1), g.constant(0)
    z = g.node('fresh_formal_Z')
    delta = g.node('common_original_delta', native_node=n['delta'], Z_independent=True)
    z2 = g.mul(z, z)
    L = g.sub(one, g.mul(delta, z2))
    numerator_factor = g.add(g.mul(g.sub(one, delta), z2), g.sub(one, z2))
    witness0 = g.sub(numerator_factor, L)
    witness1 = g.sub(g.add(g.mul(g.constant(2), g.sub(one, delta), z),
        g.mul(g.constant(-2), z)), g.mul(g.constant(-2), delta, z))
    check('new_flat_exit_Qy_c_numerator_identity_Z0_Z1',
        not g.polynomial(witness0) and not g.polynomial(witness1))
    # exp(length) is positive and has an exact multiplicative inverse. Do not
    # evaluate astronomical radii. Old m2=Z*kappa*exp(length), R2=Rp/exp(length).
    check('old_mass_back_transport_has_exact_positive_inverse',
        pw['coordinate_identity'] == 's=log(R/(2Rc)); 0<=s<=Tw-2-log2; Rp=Rc*exp(Tw-2)'
        and pc['independent_native_frame_and_geometry']['exact_Rc_2Rc_Rp_and_quiet_length_identity'] is True)
    prereq = all(checks.values())
    endpoint_proofs = []
    for role, archive, ports, K in (('native', NATIVE, native_flux, q['K']),
                                  ('original', PRE, old_flux, oq['K'])):
        endpoint_proofs.append(dict(role=role, archive=archive, RQy_ports=ports,
            K_archive=QUIET if role == 'native' else PRE, K_source_node=K,
            exact_C1_flux=[dict(value='K', ordinary_Z_order=0),
                           dict(value='0', ordinary_Z_order=1)],
            ordinary_Z_orders=[0, 1], admitted=prereq,
            source_route=('native compact RQy = actual quiet inlet RQy; actual complete mass theorem'
                if role == 'native' else
                'actual flat C3 compact exit; complete mass power-flow inlet and semigroup; actual linear Rp frame; one Pstar'),
            lemma='U=Uy=0; m=Z*c; c_Z=0; my=-m gives Qy=c and R2*Qy=K',
            Z0_axis_extension='linear source coefficient; no m/Z numerical division',
            complete_histories_not_reset=True))
    parent_graph = parent['definition_graph']
    check('accepted_parent_graph_digest', digest(parent_graph['nodes']) == parent_graph['nodes_digest'])
    rc = g.node('common_positive_Rc', source_parent=SCALE,
        native_node=n['Rc'], old_node=30828, Z_independent=True)
    outputs = []
    for k, saved in enumerate(parent_graph['outputs']):
        check('parent_C1_output_order:' + str(k), saved['ordinary_Z_order'] == k)

        def port(ref, label):
            return g.node('accepted_parent_function_port', parent_report=PARENT,
                parent_graph_digest=parent_graph['nodes_digest'], parent_node=ref,
                label=label, ordinary_Z_order=k, bound_or_midpoint=False)

        endpoints = [(i, row) for i, row in enumerate(parent_graph['nodes'])
            if row['operation'] == 'new_exact_source_endpoint_pullback'
            and row['ordinary_Z_order'] == k]
        check('actual_parent_two_endpoint_pullbacks:' + str(k),
            len(endpoints) == 2 and {row['value'] for _, row in endpoints} == {'1', '2'})
        lower = port(next(i for i, row in endpoints if row['value'] == '1'), 'DeltaQy_at_Rc')
        upper = port(next(i for i, row in endpoints if row['value'] == '2'), 'DeltaQy_at_2Rc')
        integral = port(saved['weighted_integral'], 'Rc_squared_integral_x_DeltaS')
        active = port(saved['active_difference'], 'native_Rc_minus_corrected_active')
        original_defect = port(saved['signed_compact_m15_defect'], 'accepted_compact_defect')
        gamma2 = g.mul(g.constant(2), rc, upper)
        gammaRc = g.mul(rc, lower)
        # The endpoint source theorem, not a small interval, supplies exact zero.
        upper_certificate = g.node('qualified_same_source_2Rc_flux_difference_zero',
            (gamma2, zero), admitted=prereq, ordinary_Z_order=k,
            actual_native_endpoint=native_flux[k], actual_original_endpoint=old_flux[k],
            exact_common_K_source_alias_report=SCALE,
            premises=['both actual flat linear-mass exits', 'accepted common original K'],
            inherited_numeric_errors_removed=False)
        retained = g.add(g.mul(g.constant('1/2'), integral), gammaRc)
        full = g.add(g.sub(g.mul(g.constant('1/2'), integral), gamma2), gammaRc)
        contraction = g.sub(g.sub(full, retained), g.mul(g.constant(-1), gamma2))
        check('new_affine_upper_endpoint_contraction:' + str(k), not g.polynomial(contraction))
        cert = g.node('qualified_compact_upper_endpoint_contraction',
            (original_defect, retained, upper_certificate), admitted=prereq,
            ordinary_Z_order=k, only_actual_upper_endpoint_removed=True)
        terminal = g.add(active, retained)
        outputs.append(dict(ordinary_Z_order=k, accepted_parent_compact_defect=original_defect,
            upper_2Rc_flux_difference=gamma2, upper_source_zero_certificate=upper_certificate,
            retained_Rc_flux_difference=gammaRc, weighted_integral=integral,
            compact_defect_after_upper_source_cancellation=retained,
            compact_contraction_certificate=cert, active_difference=active,
            terminal_with_Rc_flux_correlation=terminal,
            exact_identity='F-T=(native_Rc-corrected_active)+Rc^2/2*int_1^2 x*DeltaS dx+Rc*DeltaQy(1)',
            active_plus_Rc_flux_not_cancelled=True, compact_difference_zero_claimed=False))
    passed = all(checks.values())
    for row in g.nodes:
        if 'admitted' in row:
            row['admitted'] = passed
    for proof in endpoint_proofs:
        proof['admitted'] = passed
    result = dict(schema=SCHEMA, passed=passed, source_only=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True, source_family=family,
        source_prerequisite_audit=checks,
        definition_graph=dict(namespace=SCHEMA, nodes=g.nodes, node_count=len(g.nodes),
            nodes_digest=digest(g.nodes), outputs=outputs, proof_graph_not_numeric_evaluator=True,
            archived_graphs_rebuilt=False),
        source_interface=dict(endpoint_proofs=endpoint_proofs,
            new_universal_Qy_c_numerator_zero_witnesses=[witness0, witness1],
            native_actual_RQy=native_flux, original_actual_RQy=old_flux,
            common_K_alias_parent=SCALE, common_K_alias_parent_sha256=RECORD_HASHES[SCALE],
            native_K=q['K'], original_K=oq['K'], original_kappa=oq['kappa'],
            ordinary_Z1_K_derivative_zero=True,
            old_mass_transport='m_Rp=m_2Rc*exp(-length); m_Rp=Z*kappa; c2=Pstar*kappa*exp(length); K=2Rc*c2=Rp*Pstar*kappa',
            original_mass_frame_rows=fr['actual_identified_Rp_histories_C3']['m'],
            same_physical_measure='R^2 dy=Rc^2*x dx',
            required_alias_keys=[], required_identity_keys=[],
            parent_compact_graph_digest=parent_graph['nodes_digest'],
            field_interior_or_inlet_equality_assumed=False),
        retained_construction_data=parent['retained_construction_data'],
        unresolved_gates=[key for key, value in checks.items() if not value],
        remaining_construction_gates=['actual compact density integral and Rc active flux correlation',
            'native first-micro initial RQy and corrected original lower endpoint identity',
            'old212 same-source packet error theorem', 'global five-moment closure'],
        legacy_numeric_chain_reconfigured=False, new_unentered_expression_simulations=0,
        runtime_model_metadata=None,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_operators_replayed=False, numeric_values_computed=False,
            archived_definition_graphs_rebuilt=False,
            new_elementary_polynomial_and_affine_witnesses_only=True))
    result.update({key: passed for key in STAGE_FLAGS})
    result.update({key: False for key in FALSE_FLAGS})
    return result
