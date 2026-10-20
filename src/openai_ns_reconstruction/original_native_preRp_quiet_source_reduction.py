"""Resolve both exact quiet primitives in the native pre-Rp defect.

Only new affine function expressions are calculated. Archived source graphs,
their pressure constructors and the original K identity are read as data.
Their Z1 difference is zero; the physical Z0 difference is retained.
"""
from __future__ import annotations

from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_Rp_m15_source_comparison.json'
SIGNED = PREFIX + 'native_Rc_Rp_signed_pressure_source.json'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
QUIET_CHECK = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions_check.json'
QUIET_SOURCE = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.py'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
PRE_SOURCE = PREFIX + 'n1_preRp_moment_functions.py'
CORRECTED = PREFIX + 'corrected_reserved_m15_source.json'
PRIMITIVE_SOURCE = PREFIX + 'n1_quiet_pressure_primitives.py'
PRIMITIVE_CHECK = PREFIX + 'n1_quiet_pressure_primitives_check.json'
SCHEMA = 'original_native_preRp_quiet_source_reduction_v1'
RECORD_HASHES = {
    PARENT: 'c38540a019963d2023ba9fc7d4e9c292fd7af4d23fe103018d3f5009687a6586',
    SIGNED: 'a8d76d71b7d527c65c4eb29cc4e8340fc6592723246a87df65e496d3bce43a59',
    QUIET: 'bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4',
    QUIET_CHECK: '3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0',
    QUIET_SOURCE: '0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068',
    PRE: '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce',
    PRE_SOURCE: '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e',
    CORRECTED: 'fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5',
    PRIMITIVE_SOURCE: '8e024c3c1047c60e08bc249c99211d407c2865cb1d7b6de377ab5408e42db732',
    PRIMITIVE_CHECK: '8316c8bf06f3f524ef6c57aa7313f4b5546399f5bd59f2bc209a806cf71bdf0d',
}
DECODED_RECORD_HASHES = {
    QUIET: '3bc0f9984a0aeb6f447d177c40c25116b8a59c3b90ee09b8661f029cbe032a40',
    PRE: '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64',
}
STAGE_FLAGS = (
    'native_quiet_increment_closed_form_proved',
    'corrected_preRp_piece_decomposition_proved',
    'native_preRp_quiet_defect_decomposition_proved',
    'ordinary_Z1_quiet_defect_zero_proved',
)
FALSE_FLAGS = (
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
)
NATIVE_GRAPH = '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56'
CORRECTED_GRAPH = 'bb8c47eb39b4bfce588a7ffe4311d07acf23f7098d8f5d49ee07859ea2525834'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


class FunctionDifferenceGraph:
    """Exact arithmetic only on this new graph's formal function ports."""

    def __init__(self):
        self.nodes, self.forms, self.ports = [], {}, {}

    def append(self, operation, arguments, form=None, **data):
        ref = SCHEMA + ':s' + str(len(self.nodes))
        self.nodes.append(dict(id=ref, operation=operation, arguments=arguments, **data))
        self.forms[ref] = {ref: Fraction(1)} if form is None else form
        return ref

    def port(self, namespace, handle, order):
        key = (namespace, handle, order)
        if key not in self.ports:
            self.ports[key] = self.append('physical_source_function_port', [],
                namespace=namespace, source_node=handle, ordinary_Z_order=order,
                source_graph_digest=NATIVE_GRAPH if namespace == 'native_quiet_master' else CORRECTED_GRAPH,
                bound_or_midpoint=False, physical_m15=True)
        return self.ports[key]

    def linear(self, terms):
        form = {}
        for coefficient, ref in terms:
            for atom, value in self.forms[ref].items():
                form[atom] = form.get(atom, Fraction(0)) + Fraction(coefficient) * value
        return self.append('signed_function_sum', [ref for _, ref in terms],
            {key:value for key,value in form.items() if value},
            coefficients=[str(Fraction(c)) for c,_ in terms])

    def normal(self, ref):
        return {key:str(value) for key,value in sorted(self.forms[ref].items()) if value}

    def parameter(self, namespace, handle, name):
        return self.append('Z_independent_source_parameter',[],namespace=namespace,
            source_node=handle,parameter=name,ordinary_Z1_exact_zero=True,
            source_graph_digest=NATIVE_GRAPH if namespace == 'native_quiet_master' else CORRECTED_GRAPH,
            bound_or_midpoint=False)

    def quiet_primitive(self, K, length, order, owner):
        if order == 1:
            return self.append('exact_zero_quiet_Z1',[],{},owner=owner,
                K_and_length_Z_independent=True)
        return self.append('closed_quiet_pressure_primitive',[K,length],owner=owner,
            coefficient='-1/4',K_power=2,length_power=1,physical_radius_factors_retained=True)


def derive(source_records, *, expected_input_digest):
    checks = {}

    def check(key, value):
        checks[key] = bool(value)

    check('input_digest', digest(source_records) == expected_input_digest)
    for name, expected in RECORD_HASHES.items():
        record = source_records[name]
        check('authenticated:' + name, record['raw_authenticated_before_decode'] is True
            and record['raw_sha256'] == expected
            and record['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, expected)
            and hashlib.sha256(record['source_text'].encode()).hexdigest() == record['decoded_sha256'])
    data = {name:json.loads(value['source_text']) for name,value in source_records.items() if not name.endswith('.py')}
    parent = data[PARENT]['component_result']
    signed = data[SIGNED]['component_result']
    archive, receipt, pre = data[QUIET], data[QUIET_CHECK], data[PRE]
    nodes, old_nodes = archive['exact_graph_nodes'], pre['exact_graph_nodes']
    q = archive['true_quiet_Rp_C5_C3_bindings']
    interface = signed['function_source_interface']
    c = interface['compact']['source']
    corrected = data[CORRECTED]['corrected_reserved_m15_source']
    patch, roots = corrected['graph_patch'], corrected['physical_C1_new_roots']
    old = pre['signed_pieces']
    check('accepted_actual_comparison', parent['passed'] is True
        and parent['native_m15_signed_comparison_defined'] is True
        and parent['registered_reserved_m15_source_bound_once'] is True
        and parent['source_witness']['current_vs_corrected_preRp_source_identity_unproved'] is True)
    check('actual_native_quiet_owner', q == interface['quiet']['source']
        and len(nodes) == 265856 and digest(nodes) == NATIVE_GRAPH
        and parent['source_interface']['native_Rp_pressure'] == q['Rp_endpoint']['pressure'])
    check('actual_old_preRp_prefix', len(old_nodes) == 31222 and nodes[:31222] == old_nodes
        and digest(old_nodes) == corrected['base_prefix']['sha256'])
    check('source_family', archive['source_family'] == pre['source_family'] == parent['source_family'])
    check('incoming_pressure_and_integrals_retained', q['incoming_pressure'] == c['terminal']['pressure']
        and q['incoming_integrals'] == c['terminal']['integrals'] and q['initial'] == c['initial']
        and q['definition']['original_full_composer_lower_preRp_pulse_entries_not_changed'] is True)
    mass = nodes[q['linear_mass_theorem']]
    check('actual_linear_mass_source', mass['operation'] == 'proved_original_quiet_linear_mass_source_identity'
        and mass['constant_K'] == q['K'] and mass['physical_scale'] == q['Pstar']
        and mass['actual_compact_complete_mass_rows'] == q['incoming']['modified_fields']['m']
        and mass['current_quiet_identity'] == 'physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q'
        and mass['source_constant_not_a_saved_range_endpoint'] is True)
    check('saved_native_K_recipe_identity', receipt['all_passed'] is True
        and receipt['checks']['native_mass_identity']['passed'] is True
        and receipt['checks']['native_mass_identity']['actual_compact_derived_K_equals_original_Rp_Pstar_kappa'] is True
        and receipt['checks']['native_mass_identity']['common_Rp_Rc_ratio_and_finite_mu_Tw_recipe_bound'] is True
        and receipt['checks']['native_mass_identity']['no_numeric_kernel_value_or_range_endpoint_used'] is True)
    pressure_receipt = receipt['checks']['actual_physical_pressure']
    check('saved_generic_quiet_pressure_calculus', pressure_receipt['passed'] is True
        and pressure_receipt['actual_general_pressure_equals_closed_original_quiet_primitive_C3_rows'] == 8
        and pressure_receipt['finite_nonzero_initial_pressure_and_original_rminus_boundary_retained'] is True
        and pressure_receipt['physical_R_and_Rweighted_signs_checked'] is True)
    primitive = data[PRIMITIVE_CHECK]
    check('original_quiet_primitive_checked', primitive['all_passed'] is True
        and q['bindings']['original_quiet_primitive']['module_sha256'] == RECORD_HASHES[PRIMITIVE_SOURCE]
        and mass['original_quiet_identity_receipt'] == PRIMITIVE_CHECK)
    for name, required in ((QUIET_SOURCE, ('CurrentOriginalTrueQuietRpC5PressureC3Functions',)),
                           (PRE_SOURCE, ('quiet_functions','pressure_functions'))):
        tree = ast.parse(source_records[name]['source_text'])
        names = {node.name for node in ast.walk(tree) if isinstance(node,(ast.ClassDef,ast.FunctionDef))}
        check('source_owner:' + name, all(item in names for item in required))
    old_quiet = pre['defining_bindings']['quiet']
    old_q = old['quiet']['m15']
    check('old_actual_quiet_formula', old_q == [31215,0]
        and old_nodes[old_q[0]] == dict(operation='product',arguments=[31214,old_quiet['K'],old_quiet['K'],31213])
        and old_nodes[31214] == dict(operation='exact_rational',numerator=-1,denominator=4)
        and old_nodes[0] == dict(operation='exact_rational',numerator=0,denominator=1))
    check('native_actual_quiet_formula', q['quiet_closed']['m15'][1:] == [0,0,0]
        and nodes[q['quiet_closed']['m15'][0]] == dict(operation='product',arguments=[31214,265182,q['s']])
        and nodes[265182] == dict(operation='product',arguments=[q['K'],q['K']])
        and q['definition']['original_R2_and_Tw_minus2_minus_log2_geometry'] is True)
    check('corrected_patch_namespace', corrected['passed'] is True and corrected['source_only'] is True
        and patch['full_canonical_nodes_sha256'] == CORRECTED_GRAPH
        and patch['base_graph_node_count'] == 53444 and patch['full_node_count'] == 53488
        and roots['corrected_preRp_m15'] == [53472,53473]
        and parent['source_interface']['registered_corrected_physical_roots'] == roots)
    rows = {53444+i:row for i,row in enumerate(patch['appended_nodes'])}
    check('corrected_active_compact_quiet_decomposition',
        rows[53471] == dict(operation='sum',arguments=[53466,old['compact']['m15'][0]])
        and rows[53472] == dict(operation='sum',arguments=[53471,old_q[0]])
        and rows[53473] == dict(operation='sum',arguments=[53470,old['compact']['m15'][1]]))
    check('original_registration_and_errors_unchanged', parent['legacy_numeric_chain_reconfigured'] is False
        and parent['old_error_terms_removed'] is False and parent['extra_reserved_contribution_added'] is False
        and [parent['retained_construction_data'][k] for k in ('cumulative_terms','raw_signed_terms','reduced_terms')] == [212,213,183])
    qualified = all(checks.values())
    graph, outputs = FunctionDifferenceGraph(), []
    native_K = graph.parameter('native_quiet_master',q['K'],'physical native K')
    native_length = graph.parameter('native_quiet_master',q['length'],'native log(Rp/(2Rc))')
    old_K = graph.parameter('corrected_reserved_physical',old_quiet['K'],'physical old K')
    old_length = graph.parameter('corrected_reserved_physical',31213,'old log(Rp/(2Rc))')
    for order in (0,1):
        n = lambda handle: graph.port('native_quiet_master',handle,order)
        o = lambda handle: graph.port('corrected_reserved_physical',handle,order)
        rp, two_rc, rc = [n(endpoint['pressure']['m15'][order]) for endpoint in (q['Rp_endpoint'],c['terminal'],interface['native_Rc'])]
        corrected_pre, corrected_active, old_compact = o(roots['corrected_preRp_m15'][order]), o((53466,53470)[order]), o(old['compact']['m15'][order])
        old_quiet_port = o(old_q[order])
        native_quiet_formula = graph.quiet_primitive(native_K,native_length,order,'native')
        old_quiet_formula = graph.quiet_primitive(old_K,old_length,order,'old')
        native_quiet_increment = graph.linear([(1,rp),(-1,two_rc)])
        native_compact_increment = graph.linear([(1,two_rc),(-1,rc)])
        original_difference = graph.linear([(1,rp),(-1,corrected_pre)])
        corrected_expansion = graph.linear([(1,corrected_active),(1,old_compact),(1,old_quiet_formula)])
        expanded_difference = graph.linear([(1,two_rc),(1,native_quiet_formula),(-1,corrected_expansion)])
        reduced = graph.linear([(1,two_rc),(-1,corrected_active),(-1,old_compact)])
        quiet_defect = graph.linear([(1,native_quiet_formula),(-1,old_quiet_formula)])
        full_residual = graph.linear([(1,reduced),(1,quiet_defect)])
        active_difference = graph.linear([(1,rc),(-1,corrected_active)])
        compact_difference = graph.linear([(1,native_compact_increment),(-1,old_compact)])
        split = graph.linear([(1,active_difference),(1,compact_difference),(1,quiet_defect)])
        algebra = graph.normal(expanded_difference) == graph.normal(full_residual) == graph.normal(split)
        check('new_quiet_decomposition_Z' + str(order),algebra)
        if order == 1:
            check('quiet_derivative_exact_zero',graph.normal(quiet_defect) == {})
        quiet_proof = graph.append('qualified_exact_quiet_closed_forms',
            [native_quiet_increment,old_quiet_port,native_quiet_formula,old_quiet_formula],
            admitted=qualified,ordinary_Z_order=order,formula='-K^2*log(Rp/(2Rc))/4',
            actual_native_K=q['K'],actual_old_K=pre['defining_bindings']['quiet']['K'],
            native_length=q['length'],old_logratio=31213,common_physical_K_identity_proved=False,
            absolute_radius_alias_assumed=False,K_Z_exact_zero=True,increment_Z1_exact_zero=True,
            physical_quiet_Z0_difference_retained=True,source_domain_and_linear_mass_identity_required=True)
        reduction = graph.append('qualified_preRp_defect_quiet_reduction',
            [original_difference,corrected_expansion,quiet_proof,full_residual,split],
            admitted=qualified and algebra,ordinary_Z_order=order,remaining_difference_zero_claimed=False)
        outputs.append(dict(ordinary_Z_order=order,original_difference=original_difference,
            corrected_piece_expansion=corrected_expansion,native_quiet_increment=native_quiet_increment,
            old_quiet_increment=old_quiet_port,quiet_identity=quiet_proof,expanded_difference=expanded_difference,
            quiet_geometric_defect=quiet_defect,active_compact_residual=reduced,
            quiet_reduced_difference=full_residual,active_difference=active_difference,
            native_compact_increment=native_compact_increment,compact_difference=compact_difference,
            active_compact_split=split,canonical_difference=graph.normal(full_residual),
            canonical_quiet_defect=graph.normal(quiet_defect),reduction=reduction,
            residual_zero_proved=False,compact_difference_zero_proved=False,active_difference_zero_proved=False))
    passed = all(checks.values())
    result = dict(schema=SCHEMA,passed=passed,source_only=True,ordinary_Z_orders=[0,1],
        ordinary_Z_not_temporal_n=True,input_digest=expected_input_digest,source_family=parent['source_family'],
        source_prerequisite_audit=checks,definition_graph=dict(namespace=SCHEMA,nodes=graph.nodes,
            node_count=len(graph.nodes),nodes_digest=digest(graph.nodes),outputs=outputs,proof_graph_not_numeric_evaluator=True),
        quiet_source_identity=dict(formula='each quiet increment_i=-K_i^2*log(Rp_i/(2Rc_i))/4',
            ordinary_Z1='0',K_is_original_source_constant=True,source_fixture_values_not_original_field_values=True,
            actual_source_mass_theorem=q['linear_mass_theorem'],native_K=q['K'],old_K=pre['defining_bindings']['quiet']['K'],
            actual_source_regularities_and_common_domain_required=True,
            common_physical_K_identity_proved=False,absolute_radius_alias_assumed=False,
            physical_Z0_defect='-(K_native^2*length_native-K_old^2*length_old)/4'),
        source_interface=dict(native_Rc=interface['native_Rc']['pressure']['m15'][:2],
            native_2Rc=c['terminal']['pressure']['m15'][:2],native_Rp=q['Rp_endpoint']['pressure']['m15'][:2],
            corrected_preRp=roots['corrected_preRp_m15'],corrected_active=[53466,53470],
            old_compact=old['compact']['m15'],old_quiet=old_q,initial_RQ_y=q['initial']['R_Q_y'][:2],
            native_graph_digest=NATIVE_GRAPH,corrected_graph_digest=CORRECTED_GRAPH,
            physical_scales_already_applied_once=True,m15_is_not_P1=True),
        remaining_function='native_2Rc-(corrected_active+old_compact)+quiet_Z0_defect; ordinary Z1 quiet defect=0',
        retained_construction_data=parent['retained_construction_data'],legacy_numeric_chain_reconfigured=False,
        old_error_terms_removed=False,extra_reserved_contribution_added=False,compact_source_identity_proved=False,
        corrected_active_source_identity_proved=False,physical_quiet_Z0_identity_proved=False,
        source_inspection_audit=dict(new_source_consumer_calls=1,accepted_ancestor_operator_calls=0,
            old_definition_graphs_rebuilt=False,old_operators_replayed=False,numeric_values_computed=False,
            quadrature_or_norm_calls=0,source_value_reducer_calls=0,temporal_recursion_calls=0),
        unresolved_gates=[key for key,value in checks.items() if not value],runtime_model_metadata=None,
        open_obligations=['physical quiet K and absolute radius source correspondence for Z0',
            'native compact increment versus old compact physical pressure source',
            'native Rc cumulative pressure versus corrected active source and first-micro baseline',
            'old212 packet identity and independent terminal match','remaining full physical construction layers'])
    result.update({key:passed for key in STAGE_FLAGS})
    result.update({key:False for key in FALSE_FLAGS})
    return result
