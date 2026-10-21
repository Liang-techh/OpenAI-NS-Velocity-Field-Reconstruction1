"""Current-frame pulse source with current physical radial measures.

An append-only expression substitution reads the frozen complete pulse DAG.
It binds both frame jets and physical radius ports, without evaluating any
archived function, constructor, pressure normalizer or improper integral.
"""
from __future__ import annotations

from copy import deepcopy
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_preRp_quiet_source_reduction.json'
SIGNED = PREFIX + 'native_Rc_Rp_signed_pressure_source.json'
PULSE = PREFIX + 'Rp_native_C5_pulse_source.json.gz'
PULSE_CHECK = PREFIX + 'Rp_native_C5_pulse_source_check.json'
PULSE_SOURCE = PREFIX + 'Rp_native_C5_pulse_source.py'
HIGH = PREFIX + 'n1_pulse_pressure_C3_functions.json.gz'
HIGH_CHECK = PREFIX + 'n1_pulse_pressure_C3_functions_check.json'
OWN = PREFIX + 'n1_pulse_own_fields.json.gz'
OWN_SOURCE = PREFIX + 'n1_pulse_own_fields.py'
PRESSURE_SOURCE = PREFIX + 'n1_pulse_pressure_functions.py'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
SCHEMA = 'original_native_physical_pulse_source_v1'
RECORD_HASHES = {
    PARENT: 'c3fd24955c5297c457de7718710525cf855f42e5e5e672b71171baa5235b0c38',
    SIGNED: 'a8d76d71b7d527c65c4eb29cc4e8340fc6592723246a87df65e496d3bce43a59',
    PULSE: '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
    PULSE_CHECK: '320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15',
    PULSE_SOURCE: '369e6404d885dd4548e4f6451f2d5a5ed4160c0abb79a4aca2849520ede03541',
    HIGH: 'bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d',
    HIGH_CHECK: 'cbac4685fc9031a2ed5c5a7ef624b438e1a4ce3f32b85b569512aecb1d38feb8',
    OWN: '2c0f83a9ac36f4ea7970b7938bac31d36af7c19df03e9c1584a6a11aee14713b',
    OWN_SOURCE: '1400066d2da5373a95f6cd30edad19359256211d0358b454203a92fac53c5ced',
    PRESSURE_SOURCE: '8cc7442487dbbd3054e29211ea799a04d27c43fac3f060f74f5fadd929737594',
}
DECODED_RECORD_HASHES = {
    PULSE: '56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1',
    HIGH: 'a4a76f5e109db0b351c1bbc149e52a5f366e966e10a9356672342c792b0ea11e',
    OWN: '8cfee4eb2ad845c98df7cb8e7dafac6a32df473c56f02eef56751902293715a3',
}
STAGE_FLAGS = (
    'actual_native_pulse_geometry_source_bound',
    'actual_native_physical_pulse_C3_functions_defined',
    'ordinary_Z_derivatives_preserved_under_geometry_binding',
    'independent_complete_physical_pulse_m15_target_defined',
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
BASE_GRAPH = '5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e'
FRAME_KEYS = ('u','m1','m2','X','energy','Mp','P0','pressure')
OPERANDS = {
    'exact_rational': (), 'bound_variable': (), 'mathematical_pi': (),
    'mathematical_positive_infinity': (), 'current_original_source_parameter': (),
    'original_source_parameter': (), 'exact_positive_integer_expression': (),
    'shared_positive_integer': (), 'current_original_physical_source_scale': (),
    'sum': ('arguments',), 'product': ('arguments',), 'negative': ('argument',),
    'analytic_unary': ('argument',), 'positive_quotient': ('numerator','denominator'),
    'exact_real_comparison': ('left','right'), 'logical_or': ('arguments',),
    'compact_raw_beta': ('argument',), 'exact_original_source_flat_sigma': ('argument',),
    'exact_original_source_heat_phi': ('argument',),
    'current_original_leading_function_recipe': ('coordinate','Z_variable'),
    'current_original_physical_source_quantity': ('coordinate','Z_variable'),
    'current_original_physical_radius': ('log_radius','free_query'),
    'substitute_original_native_coordinate': ('body','query_coordinate'),
    'substitute_free_query_at_native_coordinate': ('body','free_query'),
    'substitute_original_inverse_angle': ('body','inverse_angle'),
    'function_substitution': ('expression','variable','value'),
    'definite_integral': ('integrand','lower','upper'),
    'original_lazy_flat_branch': ('active_body','Delta','eta','flat_predicate',
        'original_collar_coordinate','original_collar_upper'),
    'original_monotone_phase_inverse': ('phase_function','target_phase','angle_lower',
        'angle_upper','Delta','eta'),
}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),
        ensure_ascii=True,allow_nan=False).encode()).hexdigest()


def references(row):
    op = row['operation']
    if op not in OPERANDS:
        raise ValueError('Unsupported physical pulse expression: ' + op)
    result = []
    for key in OPERANDS[op]:
        value = row.get(key)
        if value is not None:
            result.extend(value if isinstance(value,list) else [value])
    if any(type(value) is not int for value in result):
        raise ValueError('Noninteger expression edge')
    return result


class CurrentPhysicalPulseExpressions:
    """Copy only changed expressions; no old graph algebra or simplification."""

    def __init__(self, base, protected):
        self.base, self.delta, self.protected = base, [], set(protected)
        self.mapping = {ref:ref for ref in protected}
        self.changed_rows, self.visited_rows = [], set()

    def append(self, row):
        ref = len(self.base) + len(self.delta)
        if any(child < 0 or child >= ref for child in references(row)):
            raise ValueError('New physical pulse graph has forward edge')
        self.delta.append(deepcopy(row))
        return ref

    def bind(self, old, current):
        if old in (0,1) or old in self.protected and current != old:
            raise ValueError('Cannot replace a global constant or current source port')
        if old in self.mapping and self.mapping[old] != current:
            raise ValueError('Ambiguous physical source binding')
        self.mapping[old] = current

    def rewrite(self, root):
        stack, visiting = [(root,False)], set()
        while stack:
            ref, ready = stack.pop()
            if ref in self.mapping:
                continue
            if not 0 <= ref < len(self.base):
                raise ValueError('Original pulse source handle outside frozen base')
            row = self.base[ref]
            children = references(row)
            if not ready:
                if ref in visiting:
                    raise ValueError('Cyclic original pulse expression')
                visiting.add(ref)
                stack.append((ref,True))
                stack.extend((child,False) for child in reversed(children) if child not in self.mapping)
                continue
            new = deepcopy(row)
            for key in OPERANDS[row['operation']]:
                value = row.get(key)
                if value is not None:
                    new[key] = [self.mapping[n] for n in value] if isinstance(value,list) else self.mapping[value]
            self.mapping[ref] = ref if new == row else self.append(new)
            self.visited_rows.add(ref)
            if new != row:
                self.changed_rows.append(ref)
            visiting.remove(ref)
        return self.mapping[root]

    def tree(self, value):
        if isinstance(value,dict):
            return {key:self.tree(item) for key,item in value.items()}
        if isinstance(value,list):
            return [self.tree(item) for item in value]
        if type(value) is not int:
            raise ValueError('Only explicit function-root trees can be substituted')
        return self.rewrite(value)

    def reachable(self, roots):
        seen, stack = set(), list(roots)
        while stack:
            ref = stack.pop()
            if ref in seen:
                continue
            seen.add(ref)
            if ref in self.protected:
                continue
            row = self.base[ref] if ref < len(self.base) else self.delta[ref-len(self.base)]
            stack.extend(references(row))
        return seen


def flatten(value):
    if isinstance(value,dict):
        return [n for child in value.values() for n in flatten(child)]
    if isinstance(value,list):
        return [n for child in value for n in flatten(child)]
    if type(value) is not int:
        raise ValueError('A scalar expression handle is required')
    return [value]


def derive(source_records, *, expected_input_digest):
    checks = {}

    def check(name, value):
        checks[name] = bool(value)

    check('input_digest',digest(source_records) == expected_input_digest)
    for name, expected in RECORD_HASHES.items():
        record = source_records[name]
        check('authenticated:' + name,record['raw_authenticated_before_decode'] is True
            and record['raw_sha256'] == expected
            and record['decoded_sha256'] == DECODED_RECORD_HASHES.get(name,expected)
            and hashlib.sha256(record['source_text'].encode()).hexdigest() == record['decoded_sha256'])
    data = {name:json.loads(record['source_text']) for name,record in source_records.items() if not name.endswith('.py')}
    parent, signed = data[PARENT], data[SIGNED]['component_result']
    pulse, high, own = data[PULSE], data[HIGH], data[OWN]
    pb, h, fields = pulse['Rp_native_C5_pulse_source_bindings'], high['C5_pulse_C3_pressure_bindings'], own['own_main_gap_end_fields']
    base, old, own_nodes = pulse['exact_graph_nodes'], high['exact_graph_nodes'], own['exact_graph_nodes']
    q = pb['inputs'][QUIET]
    check('accepted_parent_and_open_geometry_gate',parent['component_result']['passed'] is True
        and parent['component_result']['ordinary_Z1_quiet_defect_zero_proved'] is True
        and parent['parent_complete_pulse_physical_geometry_binding_proved'] is False)
    check('authenticated_master_and_function_prefixes',len(base) == 280182 and digest(base) == BASE_GRAPH
        and len(old) == 62795 and base[:len(old)] == old
        and len(own_nodes) == 34295 and old[:len(own_nodes)] == own_nodes)
    check('actual_current_frame_and_geometry',q == signed['function_source_interface']['quiet']['source']
        and pb['current_frame'] == q['Rp_frame']
        and pb['inputs'][HIGH] == h
        and pb['geometry']['Rp'] == q['Rp'] == 263349
        and [pb['geometry'][key] for key in ('Pstar','mu','N','delta')] == [33572,17,18,25220])
    check('whole_function_native_frame_qualification',data[PULSE_CHECK]['all_passed'] is True
        and data[PULSE_CHECK]['current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions'] is True
        and set(pb['current_frame']) == set(FRAME_KEYS)
        and all(len(rows) == 6 for rows in pb['current_frame'].values()))
    check('complete_C3_source_qualification',data[HIGH_CHECK]['all_passed'] is True
        and h['definition']['actual_signed_n1_pulse_pressure_C3_integral_functions'] is True
        and h['definition']['physical_Pstar_radial_measures_and_finite_endpoints_retained'] is True
        and all(len(rows) == 4 for rows in h['pressure'].values())
        and set(h['integrals']) == set(h['leading_fields']) == {'main','gap','end'})
    check('common_pulse_Pstar_is_an_exact_port',base[33572] == dict(operation='analytic_unary',name='exp',argument=5)
        and own_nodes[33573] == dict(operation='product',arguments=[h['native']['u'][0],33572]))
    check('old_radius_ports_are_actual_expressions',old[31162] == dict(operation='analytic_unary',name='exp',argument=31161)
        and old[34349] == dict(operation='product',arguments=[31162,31162]))
    check('source_family',parent['component_result']['source_family'] == pulse['source_family'] == high['source_family'] == own['source_family'])
    for name, required in ((OWN_SOURCE,('pulse_fields','radial_Q')),
                           (PRESSURE_SOURCE,('pulse_pressure_functions','actual_outer_endpoints'))):
        tree = ast.parse(source_records[name]['source_text'])
        names = {node.name for node in ast.walk(tree) if isinstance(node,ast.FunctionDef)}
        check('actual_definition:' + name,all(item in names for item in required))
    protected = set(flatten(pb['current_frame']) + list(pb['geometry'].values()) + pb['current_P0'])
    expr = CurrentPhysicalPulseExpressions(base,protected)
    logRp = expr.append(dict(operation='analytic_unary',name='log',argument=q['Rp'],
        source_positive_certificate='actual native Rp>0',ordinary_Z1_exact_zero=True))
    geometry_bindings = {31161:logRp,31162:q['Rp']}
    for chart, field in fields.items():
        coordinate = pb['pulse_coordinates'][chart]
        check('actual_chart_coordinate:' + chart,field['query'] == coordinate['query'] and field['y'] == coordinate['y'])
        geometry_bindings[field['radius']] = coordinate['radius']
    for old_ref, new_ref in geometry_bindings.items():
        expr.bind(old_ref,new_ref)
    frame_bindings = {}
    for key, rows in h['native'].items():
        for order, old_ref in enumerate(rows):
            if old_ref not in (0,1):
                current = pb['current_frame'][key][order]
                expr.bind(old_ref,current)
                frame_bindings[old_ref] = current
    groups = {name:expr.tree(h[name]) for name in ('selected_pulse','leading_fields','integrals','endpoints','global_endpoints','pressure','pieces')}
    radii = {chart:expr.rewrite(field['radius']) for chart,field in fields.items()}
    check('all_physical_radius_ports_are_current',all(radii[chart] == pb['pulse_coordinates'][chart]['radius'] for chart in radii))
    roots = flatten(groups)
    closure = expr.reachable(roots)
    forbidden = {31161,31162,24956,30828,33805,33982,34272,34349}
    check('old_physical_geometry_absent_outside_current_ports',not bool(closure & forbidden))
    check('actual_current_radius_reaches_m15',q['Rp'] in expr.reachable(groups['pressure']['m15']))
    check('current_input_ports_unchanged',all(expr.mapping[ref] == ref for ref in protected)
        and pb['current_frame'] == q['Rp_frame'])
    check('geometry_binding_Z_independent',base[q['N']]['Z_independent'] is True
        and q['definition']['original_R2_and_Tw_minus2_minus_log2_geometry'] is True
        and all(field['domain']['coordinate'] == ('s=log(R/Rv)' if chart=='end' else 'xi=mu*log(R/Rp)') for chart,field in fields.items()))
    Rv = expr.append(dict(operation='function_substitution',expression=radii['end'],
        variable=fields['end']['query'],value=0,Z_independent_substitution=True,
        definition='actual Rv=actual Rp*exp(13/mu), end coordinate s=0'))
    negative_target = [expr.append(dict(operation='negative',argument=ref)) for ref in groups['pressure']['m15']]
    pulse_correction = []
    for current, previous in zip(groups['pressure']['m15'],pb['pressure']['m15']):
        negative = expr.append(dict(operation='negative',argument=previous))
        pulse_correction.append(expr.append(dict(operation='sum',arguments=[current,negative])))
    passed = all(checks.values())
    joined_digest = digest(base + expr.delta)
    comparison = [dict(ordinary_Z_order=order,
        preRp_defect_source=dict(report=PARENT,report_sha256=RECORD_HASHES[PARENT],
            function_ref=parent['component_result']['definition_graph']['outputs'][order]['original_difference']),
        physical_pulse_geometry_correction=pulse_correction[order],
        source_identity='Delta_physical=Delta_preRp+(pulse_physical-pulse_frame_only)',
        complete_terminal_difference_claimed_zero=False) for order in (0,1)]
    result = dict(schema=SCHEMA,passed=passed,source_only=True,ordinary_Z_orders=[0,1,2,3],
        ordinary_Z_not_temporal_n=True,input_digest=expected_input_digest,source_family=pulse['source_family'],
        source_prerequisite_audit=checks,
        definition_graph=dict(base_node_count=len(base),base_nodes_digest=BASE_GRAPH,appended_nodes=expr.delta,
            appended_node_count=len(expr.delta),full_node_count=len(base)+len(expr.delta),full_nodes_digest=joined_digest,
            old_prefix_exactly_preserved=True,source_graph_not_numeric_evaluator=True),
        physical_pulse_bindings=dict(current_frame=pb['current_frame'],current_P0=pb['current_P0'],
            geometry=pb['geometry'],physical_radii=radii,coordinates=pb['pulse_coordinates'],
            current_logRp=logRp,current_Rv=Rv,geometry_substitution={str(k):v for k,v in geometry_bindings.items()},
            frame_substitution={str(k):v for k,v in frame_bindings.items()},**groups,
            negative_complete_m15_target=negative_target,pulse_geometry_correction_m15=pulse_correction),
        substitution_receipt=dict(affected_original_rows=len(expr.changed_rows),visited_original_rows=len(expr.visited_rows),
            changed_row_mapping={str(ref):expr.mapping[ref] for ref in expr.changed_rows},
            protected_current_source_ports=sorted(protected),forbidden_old_geometry_ports=sorted(forbidden),
            old_geometry_outside_current_ports_present=False,source_coefficients_or_field_values_evaluated=False,
            ordinary_Z_preserved_by_exact_C5_inlet_and_Z_independent_geometry=True),
        physical_scale_receipt=dict(Pstar=33572,old_pulse_Pstar=33572,common_exact_expression_port=True,
            extra_Pstar_factor_applied=False,old_preRp_physical_scale_projection_25203_not_aliased=True,
            radial_measures='R/mu and R^2/mu on main/gap; R and R^2 on end',
            all_three_charts_and_finite_endpoints_retained=True,m15_paper_sign_retained=True,target_negative_applied_once=True,
            positive_actual_Rp_and_original_common_source_domain_required=True),
        terminal_comparison=comparison,retained_construction_data=parent['retained_construction_data'],
        old_physical_pulse_geometry_alias_proved=False,new_physical_pulse_equals_old_pulse_proved=False,
        legacy_numeric_chain_reconfigured=False,extra_reserved_contribution_added=False,old_error_terms_removed=False,
        source_inspection_audit=dict(new_source_consumer_calls=1,accepted_ancestor_operator_calls=0,
            old_definition_graphs_rebuilt=False,old_operators_replayed=False,numeric_values_computed=False,
            quadrature_or_norm_calls=0,source_value_reducer_calls=0,temporal_recursion_calls=0,
            new_changed_expression_nodes=len(expr.delta)),
        unresolved_gates=[name for name,value in checks.items() if not value],runtime_model_metadata=None,
        open_obligations=['physical pulse geometry correction versus registered reserved source',
            'native versus corrected preRp active/compact/quiet Z0 function correspondence',
            'old212 packet identity and independent terminal matching',
            'original-scale bounds, global physical seams, heat, stress, temporal recursion and corrected Cartesian NS'])
    result.update({name:passed for name in STAGE_FLAGS})
    result.update({name:False for name in FALSE_FLAGS})
    return result
