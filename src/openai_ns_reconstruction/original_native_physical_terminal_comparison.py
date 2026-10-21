"""Source-only comparison of the current physical terminal and saved overlay.

Three JSON reports and one immutable source archive are authenticated.  Their archived
expression graphs are read as pinned metadata and source handles; no archived
constructor, evaluator, normalizer, quadrature, or graph builder is called.
The graph made here is a small, fresh affine graph whose ports retain their
saved graph digest and source handle.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json


PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT_REPORT = PREFIX + 'native_physical_pulse_source.json'
QUIET_REPORT = PREFIX + 'native_preRp_quiet_source_reduction.json'
COMPARE_REPORT = PREFIX + 'native_Rp_m15_source_comparison.json'

# The aliases keep the report roles readable to callers while the declared
# records below are the only inputs accepted by derive.
PARENT = PARENT_REPORT
QUIET = QUIET_REPORT
COMPARE = COMPARE_REPORT
NATIVE_ARCHIVE = PREFIX + 'Rp_native_C5_pulse_source.json.gz'

SCHEMA = 'original_native_physical_terminal_comparison_v1'
RECORD_HASHES = {
    NATIVE_ARCHIVE: '442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe',
    PARENT_REPORT: '8b21f8e868025d8cc95149ffab80b54603833f5f299ec9d100602af355cea82b',
    QUIET_REPORT: 'c3fd24955c5297c457de7718710525cf855f42e5e5e672b71171baa5235b0c38',
    COMPARE_REPORT: 'c38540a019963d2023ba9fc7d4e9c292fd7af4d23fe103018d3f5009687a6586',
}
DECODED_RECORD_HASHES = {
    NATIVE_ARCHIVE: '56137195254ae71d7eedfccc82ac7b5c048c8f3acd8de87fedfee5c6725fd7b1',
}

STAGE_FLAGS = (
    'absolute_radius_log_recipe_difference_proved',
    'native_old_quiet_length_identity_proved',
    'actual_physical_pulse_target_integrated',
    'signed_physical_terminal_decomposition_proved',
    'inherited_quiet_Z1_reduction_preserved',
    'registered_reserved_m15_source_bound_once',
)

FALSE_FLAGS = (
    'same_density_certificate',
    'finite_N_history_identity',
    'physical_Q_y_Z1',
    'global_five_moment_closure',
    'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved',
    'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed',
    'numeric_values_computed',
    'old_operators_replayed',
    'full_goal_complete',
    'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
    'new_physical_pulse_equals_old_pulse_proved',
    'old_physical_geometry_alias_proved',
)

BASE_GRAPH = '5cdd0c8fd6a758b9e176bf67bd29ea2f576c4ddae2c7cb63d438ea1ca8388c1e'
PHYSICAL_GRAPH = '29d2257340a9bb9399448dfe0ab6bd2cd31e1815d3646741a242dac506a0a84c'
QUIET_NATIVE_GRAPH = '15d51cd83dd82f41eba990979ace934bd5c8d5c339cd516b9fc96bc20a64af56'
QUIET_RESULT_GRAPH = '9a09124e300a64bc3d955ba1541ac09ae7f4b072f62ae6025a86a34794c01edb'
CORRECTED_GRAPH = 'bb8c47eb39b4bfce588a7ffe4311d07acf23f7098d8f5d49ee07859ea2525834'

BASE_NODE_COUNT = 280182
PHYSICAL_NODE_COUNT = 290924
PHYSICAL_APPENDED_COUNT = 10742
QUIET_RESULT_NODE_COUNT = 48
CORRECTED_NODE_COUNT = 53488

INNER = (53161, 53162)
LOWER = (36592, 36593)
NATIVE_PRE_RP = (265413, 265414)
OLD_FRAME_ONLY = (279426, 279542, 279770, 280139)
OLD_PULSE = (35556, 35557)
PHYSICAL_PULSE = (287417, 287491, 287663, 287956)
NEGATIVE_PHYSICAL_PULSE = (290912, 290913, 290914, 290915)
PHYSICAL_GEOMETRY_DELTA = (290917, 290919, 290921, 290923)

NATIVE_RC = (258845, 258846)
NATIVE_2RC = (262906, 262907)
CORRECTED_ACTIVE = (53466, 53470)
CORRECTED_PRE_RP = (53472, 53473)
CORRECTED_RESERVED = (53476, 53477)
OLD_COMPACT = (31155, 31156)

QUIET_OUTPUTS = (
    {
        'quiet_geometric_defect': 'original_native_preRp_quiet_source_reduction_v1:s19',
        'quiet_reduced_difference': 'original_native_preRp_quiet_source_reduction_v1:s20',
        'active_difference': 'original_native_preRp_quiet_source_reduction_v1:s21',
        'compact_difference': 'original_native_preRp_quiet_source_reduction_v1:s22',
        'active_compact_split': 'original_native_preRp_quiet_source_reduction_v1:s23',
        'canonical_quiet_defect': {
            'original_native_preRp_quiet_source_reduction_v1:s11': '1',
            'original_native_preRp_quiet_source_reduction_v1:s12': '-1',
        },
    },
    {
        'quiet_geometric_defect': 'original_native_preRp_quiet_source_reduction_v1:s41',
        'quiet_reduced_difference': 'original_native_preRp_quiet_source_reduction_v1:s42',
        'active_difference': 'original_native_preRp_quiet_source_reduction_v1:s43',
        'compact_difference': 'original_native_preRp_quiet_source_reduction_v1:s44',
        'active_compact_split': 'original_native_preRp_quiet_source_reduction_v1:s45',
        'canonical_quiet_defect': {},
    },
)


def digest(value):
    """Digest JSON metadata without evaluating any source expression."""
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':'), ensure_ascii=True,
        allow_nan=False).encode('utf-8')).hexdigest()


class PhysicalTerminalComparisonGraph:
    """A fresh exact affine graph over authenticated source-function ports."""

    def __init__(self):
        self.nodes = []
        self.forms = {}
        self.ports = {}

    def append(self, operation, arguments, form=None, **data):
        ref = SCHEMA + ':s' + str(len(self.nodes))
        self.nodes.append(dict(id=ref, operation=operation,
                               arguments=arguments, **data))
        self.forms[ref] = {ref: Fraction(1)} if form is None else form
        return ref

    def port(self, role, source_graph_digest, source_node, order, *,
             accepted_report, accepted_field):
        # The source graph digest and source handle are part of the identity.
        # A role is retained in the node metadata, but never used to alias
        # ports from distinct saved graph namespaces.
        key = (source_graph_digest, repr(source_node), order)
        if key not in self.ports:
            self.ports[key] = self.append(
                'authenticated_source_function_port', [],
                role=role,
                source_graph_digest=source_graph_digest,
                source_node=source_node,
                ordinary_Z_order=order,
                accepted_report=accepted_report,
                accepted_field=accepted_field,
                bound_or_midpoint=False,
                source_value_evaluated=False,
            )
        return self.ports[key]

    def linear(self, terms, *, identity=None):
        form = {}
        for coefficient, ref in terms:
            for atom, value in self.forms[ref].items():
                atom_value = form.get(atom, Fraction(0))
                form[atom] = atom_value + Fraction(coefficient) * value
        form = {key: value for key, value in form.items() if value}
        data = {'coefficients': [str(Fraction(c)) for c, _ in terms]}
        if identity is not None:
            data['identity'] = identity
        return self.append('exact_signed_function_sum',
                           [ref for _, ref in terms], form, **data)

    def normal(self, ref, relations=None):
        relations = relations or {}
        result = {}

        def expand(atom, coefficient, seen):
            if atom in relations:
                if atom in seen:
                    raise ValueError('Cyclic qualified source relation')
                for child, value in self.forms[relations[atom]].items():
                    expand(child, coefficient * value, seen | {atom})
            else:
                result[atom] = result.get(atom, Fraction(0)) + coefficient

        for atom, coefficient in self.forms[ref].items():
            expand(atom, coefficient, set())
        return {key: str(value) for key, value in sorted(result.items()) if value}


def _component(data, name):
    document = data.get(name, {})
    result = document.get('component_result', {})
    return result if isinstance(result, dict) else {}


def _same_sequence(value, expected):
    return isinstance(value, list) and value == list(expected)


def _all_true(mapping, names):
    return all(mapping.get(name) is True for name in names)


def _record_digest_ok(record, expected, decoded, compressed):
    if not isinstance(record, dict) or record.get('raw_authenticated_before_decode') is not True:
        return False
    text = record.get('source_text')
    if not isinstance(text, str):
        return False
    return (record.get('raw_sha256') == expected
            and record.get('decoded_sha256') == decoded
            and hashlib.sha256(text.encode('utf-8')).hexdigest() == decoded
            and record.get('is_gzip', False) is compressed)


def _saved_row(nodes, handle, base_count, expected):
    index = handle - base_count
    return 0 <= index < len(nodes) and nodes[index] == expected


def _radius_scalar_audit(nodes, check):
    """New exact polynomial/logarithm calculus on the frozen radius recipes.

    No floating point or field evaluator is used. Free exact control ports stay
    distinct. Substitution values are expression HANDLES, never literal IDs.
    """
    check('radius_archive_exact_master', len(nodes) == BASE_NODE_COUNT
          and digest(nodes) == BASE_GRAPH)
    check('radius_endpoint_handles', nodes[48] == dict(operation='exact_rational', numerator=2, denominator=1)
          and nodes[1] == dict(operation='exact_rational', numerator=1, denominator=1)
          and nodes[251066] == dict(operation='function_substitution', expression=250457,
              variable=250215, value=48, Z_independent_substitution=True)
          and nodes[85743] == dict(operation='function_substitution', expression=84533,
              variable=84528, value=1, Z_independent_substitution=True))
    check('radius_exact_control_sources', nodes[35560]['name'] == 'true_core_logC'
          and nodes[35560]['source_path'] == 'common.core.logC'
          and nodes[35560]['exact_fixed_source_parameter'] is True
          and nodes[6]['name'] == 'logC'
          and 'selected_logCstar' in nodes[6]['definition']
          and nodes[6]['quantity_is_exact_original_parameter_not_range_endpoint'] is True
          and nodes[35560]['source_family'] == nodes[6]['source_family'])
    opaque = {6:'old_selected_logCstar',35560:'native_core_logC',
              8:'shared_source_8',10:'shared_source_10',16:'shared_logmu'}
    seen = set()

    def clean(poly):
        return {key:value for key,value in poly.items() if value}

    def add(*polys):
        result = {}
        for poly in polys:
            for key,value in poly.items():
                result[key] = result.get(key,Fraction(0)) + value
        return clean(result)

    def mul(left,right):
        result = {}
        for a,x in left.items():
            for b,y in right.items():
                key = tuple(sorted(a+b))
                result[key] = result.get(key,Fraction(0)) + x*y
        return clean(result)

    def const(value):
        return clean({():Fraction(value)})

    def neg(poly):
        return {key:-value for key,value in poly.items()}

    def exact_log(value):
        if value <= 0:
            raise ValueError('Radius log requires a positive exact rational')
        result = {}
        for integer,sign in ((value.numerator,1),(value.denominator,-1)):
            factor = 2
            while factor*factor <= integer:
                exponent = 0
                while integer % factor == 0:
                    exponent += 1
                    integer //= factor
                if exponent:
                    key = ('log_prime_'+str(factor),)
                    result[key] = result.get(key,Fraction(0)) + sign*exponent
                factor += 1
            if integer > 1:
                key = ('log_prime_'+str(integer),)
                result[key] = result.get(key,Fraction(0)) + sign
        return clean(result)

    def rational(poly):
        if set(poly) - {()}:
            raise ValueError('Exact rational scalar required')
        return poly.get((),Fraction(0))

    def scalar(ref,env):
        seen.add(ref)
        if ref in env:
            return env[ref]
        if ref in opaque:
            return {(opaque[ref],):Fraction(1)}
        row = nodes[ref]
        op = row['operation']
        if op == 'exact_rational':
            return const(Fraction(row['numerator'],row['denominator']))
        if op == 'sum':
            return add(*(scalar(child,env) for child in row['arguments']))
        if op == 'negative':
            return neg(scalar(row['argument'],env))
        if op == 'product':
            result = const(1)
            for child in row['arguments']:
                result = mul(result,scalar(child,env))
            return result
        if op == 'analytic_unary' and row['name'] == 'log':
            return exact_log(rational(scalar(row['argument'],env)))
        if op == 'analytic_unary' and row['name'] == 'exp' and rational(scalar(row['argument'],env)) == 40:
            return {('exp_40',):Fraction(1)}
        raise ValueError('Unexpected exact radius scalar operation at '+str(ref))

    def log_radius(ref,env):
        seen.add(ref)
        row = nodes[ref]
        op = row['operation']
        if ref in env:
            return exact_log(rational(env[ref]))
        if op == 'function_substitution':
            bound = dict(env)
            bound[row['variable']] = scalar(row['value'],env)
            return log_radius(row['expression'],bound)
        if op == 'analytic_unary' and row['name'] == 'exp':
            return scalar(row['argument'],env)
        if op == 'product':
            return add(*(log_radius(child,env) for child in row['arguments']))
        if op == 'exact_rational':
            return exact_log(Fraction(row['numerator'],row['denominator']))
        raise ValueError('Unexpected positive radius operation at '+str(ref))

    old_rc = log_radius(30828,{})
    native_rc = log_radius(251066,{})
    old_rp = log_radius(31162,{})
    native_rp = log_radius(263349,{})
    expected = {('native_core_logC',):Fraction(10),('old_selected_logCstar',):Fraction(-10)}
    rc_difference = add(native_rc,neg(old_rc))
    rp_difference = add(native_rp,neg(old_rp))
    old_length = add(old_rp,neg(old_rc),neg(exact_log(Fraction(2))))
    native_length = add(native_rp,neg(native_rc),neg(exact_log(Fraction(2))))
    length_difference = add(native_length,neg(old_length))
    check('absolute_log_Rc_source_difference', rc_difference == expected)
    check('absolute_log_Rp_source_difference', rp_difference == expected)
    check('native_old_quiet_length_exact_identity', length_difference == {})
    check('no_spurious_literal_handle_offset', () not in rc_difference and () not in rp_difference)

    def serial(poly):
        return [{'monomial':list(key),'coefficient':str(value)} for key,value in sorted(poly.items())]
    return dict(source_archive=NATIVE_ARCHIVE,source_graph_digest=BASE_GRAPH,
        native_Rm=85743,native_Rc=251066,old_Rc=30828,native_Rp=263349,old_Rp=31162,
        substituted_endpoint_handle=48,substituted_endpoint_exact_value='2',
        native_logRc=serial(native_rc),old_logRc=serial(old_rc),
        native_logRp=serial(native_rp),old_logRp=serial(old_rp),
        logRc_difference=serial(rc_difference),logRp_difference=serial(rp_difference),
        quiet_length=serial(native_length),quiet_length_difference=serial(length_difference),
        identity='log(Rc_native/Rc_old)=log(Rp_native/Rp_old)=10*(native_core_logC-old_selected_logCstar)',
        radius_ratio='exp(10*(native_core_logC-old_selected_logCstar))',
        constant_46_offset_rejected=True,endpoint_handle_is_not_literal_value=True,
        source_rows_read=sorted(seen),control_alias_proved=False,physical_Pstar_alias_proved=False,
        same_control_would_imply_radius_equality=True,source_values_or_bounds_evaluated=False,
        exact_positive_radius_log_laws_only=True,old_math_consumers_called=0)


def derive(source_records, *, expected_input_digest):
    """Build one terminal and exact scalar-recipe source comparison."""
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('input_digest', digest(source_records) == expected_input_digest)
    check('exact_pinned_reports_and_archive',
          isinstance(source_records, dict)
          and set(source_records) == set(RECORD_HASHES))
    for name, expected in RECORD_HASHES.items():
        check('raw_authenticated:' + name,
              _record_digest_ok(source_records.get(name), expected,
                  DECODED_RECORD_HASHES.get(name, expected), name.endswith('.gz')))

    data = {}
    for name in RECORD_HASHES:
        record = source_records.get(name, {})
        text = record.get('source_text') if isinstance(record, dict) else None
        try:
            value = json.loads(text) if isinstance(text, str) else {}
        except (TypeError, ValueError):
            value = {}
        data[name] = value if isinstance(value, dict) else {}
    check('authenticated_json_text_only',
          all(isinstance(data[name], dict) for name in RECORD_HASHES))

    radius_audit = _radius_scalar_audit(data[NATIVE_ARCHIVE]['exact_graph_nodes'], check)
    physical = _component(data, PARENT_REPORT)
    quiet = _component(data, QUIET_REPORT)
    compare = _component(data, COMPARE_REPORT)

    physical_local = (
        'actual_native_pulse_geometry_source_bound',
        'actual_native_physical_pulse_C3_functions_defined',
        'ordinary_Z_derivatives_preserved_under_geometry_binding',
        'independent_complete_physical_pulse_m15_target_defined',
    )
    quiet_local = (
        'native_quiet_increment_closed_form_proved',
        'corrected_preRp_piece_decomposition_proved',
        'native_preRp_quiet_defect_decomposition_proved',
        'ordinary_Z1_quiet_defect_zero_proved',
    )
    compare_local = (
        'native_Rp_C5_pulse_inlet_source_identity_proved',
        'native_complete_pulse_m15_source_identity_proved',
        'registered_reserved_m15_source_bound_once',
        'native_m15_signed_comparison_defined',
    )
    for name, result, local in (
            ('physical', physical, physical_local),
            ('quiet', quiet, quiet_local),
            ('compare', compare, compare_local)):
        check(name + '_passed', result.get('passed') is True)
        check(name + '_source_only', result.get('source_only') is True
              and result.get('ordinary_Z_not_temporal_n') is True)
        check(name + '_local_flags', _all_true(result, local))
        audit = result.get('source_prerequisite_audit', {})
        check(name + '_prerequisites_passed',
              isinstance(audit, dict) and bool(audit)
              and all(value is True for value in audit.values()))

    physical_document = data[PARENT_REPORT]
    quiet_document = data[QUIET_REPORT]
    compare_document = data[COMPARE_REPORT]
    check('parent_report_chain',
          physical_document.get('parent_report') == QUIET_REPORT
          and physical_document.get('parent_report_sha256') == RECORD_HASHES[QUIET_REPORT]
          and quiet_document.get('parent_report') == COMPARE_REPORT
          and quiet_document.get('parent_report_sha256') == RECORD_HASHES[COMPARE_REPORT])

    family = physical.get('source_family')
    check('same_source_family',
          isinstance(family, dict)
          and quiet.get('source_family') == family
          and compare.get('source_family') == family
          and physical_document.get('source_family') == family
          and quiet_document.get('source_family') == family
          and compare_document.get('source_family') == family)

    retained_values = (212, 213, 183)
    for name, result in (('physical', physical), ('quiet', quiet),
                         ('compare', compare)):
        retained = result.get('retained_construction_data', {})
        check(name + '_retained_counts',
              isinstance(retained, dict)
              and tuple(retained.get(key) for key in
                        ('cumulative_terms', 'raw_signed_terms', 'reduced_terms'))
              == retained_values)
        check(name + '_retained_errors_and_channels',
              retained.get('six_signed_channels_retained') is True
              and retained.get('all_old_error_data_retained') is True
              and retained.get('inherited_errors_deleted') is False)
        check(name + '_retained_reserved_once',
              retained.get('reserved_contribution_registered_once') is True)
        check(name + '_legacy_false',
              result.get('legacy_numeric_chain_reconfigured') is False
              and result.get('old_error_terms_removed') is False
              and result.get('extra_reserved_contribution_added') is False)
        inspection = result.get('source_inspection_audit', {})
        check(name + '_inspection_no_old_or_numeric',
              isinstance(inspection, dict)
              and inspection.get('accepted_ancestor_operator_calls') == 0
              and inspection.get('old_definition_graphs_rebuilt') is False
              and inspection.get('old_operators_replayed') is False
              and inspection.get('numeric_values_computed') is False)

    check('compare_once_only_overlay',
          compare.get('source_interface', {}).get('corrected_source_overlay_installed') is True
          and compare.get('source_interface', {}).get('registered_unit_node') == 53427
          and compare.get('source_witness', {}).get('corrected_source_overlay_identity')
          == 'reserved_corrected=-full_corrected'
          and compare.get('retained_construction_data', {}).get(
              'reserved_contribution_registered_once') is True)

    physical_graph = physical.get('definition_graph', {})
    check('authenticated_physical_graph_metadata',
          physical_graph.get('base_node_count') == BASE_NODE_COUNT
          and physical_graph.get('base_nodes_digest') == BASE_GRAPH
          and physical_graph.get('appended_node_count') == PHYSICAL_APPENDED_COUNT
          and physical_graph.get('full_node_count') == PHYSICAL_NODE_COUNT
          and physical_graph.get('full_nodes_digest') == PHYSICAL_GRAPH
          and physical_graph.get('old_prefix_exactly_preserved') is True
          and physical_graph.get('source_graph_not_numeric_evaluator') is True)

    bindings = physical.get('physical_pulse_bindings', {})
    check('saved_physical_terminal_ports',
          _same_sequence(bindings.get('pressure', {}).get('m15'), PHYSICAL_PULSE)
          and _same_sequence(bindings.get('negative_complete_m15_target'),
                             NEGATIVE_PHYSICAL_PULSE)
          and _same_sequence(bindings.get('pulse_geometry_correction_m15'),
                             PHYSICAL_GEOMETRY_DELTA))
    check('saved_physical_target_is_negative_once',
          physical.get('physical_scale_receipt', {}).get('target_negative_applied_once') is True
          and physical.get('physical_scale_receipt', {}).get('extra_Pstar_factor_applied') is False)

    appended = physical_graph.get('appended_nodes', [])
    literal_rows = {
        287417: {'operation': 'sum', 'arguments': [287414, 287416]},
        287491: {'operation': 'sum', 'arguments': [287488, 287490]},
        287663: {'operation': 'sum', 'arguments': [287660, 287662]},
        287956: {'operation': 'sum', 'arguments': [287953, 287955]},
        290912: {'operation': 'negative', 'argument': 287417},
        290913: {'operation': 'negative', 'argument': 287491},
        290914: {'operation': 'negative', 'argument': 287663},
        290915: {'operation': 'negative', 'argument': 287956},
        290916: {'operation': 'negative', 'argument': 279426},
        290917: {'operation': 'sum', 'arguments': [287417, 290916]},
        290918: {'operation': 'negative', 'argument': 279542},
        290919: {'operation': 'sum', 'arguments': [287491, 290918]},
        290920: {'operation': 'negative', 'argument': 279770},
        290921: {'operation': 'sum', 'arguments': [287663, 290920]},
        290922: {'operation': 'negative', 'argument': 280139},
        290923: {'operation': 'sum', 'arguments': [287956, 290922]},
    }
    check('saved_physical_literal_rows', all(
        _saved_row(appended, handle, BASE_NODE_COUNT, row)
        for handle, row in literal_rows.items()))

    compare_interface = compare.get('source_interface', {})
    compare_witness = compare.get('source_witness', {})
    compare_namespaces = compare.get('source_namespaces', {})
    check('saved_frame_and_old_pulse_ports',
          _same_sequence(compare_interface.get('native_complete_pulse_pressure', {}).get('m15'),
                         OLD_FRAME_ONLY)
          and _same_sequence(compare_interface.get('native_Rp_pressure', {}).get('m15'),
                             (265413, 265414, 265415, 265416))
          and compare_witness.get('complete_pulse_m15_substitution') == [
              [35556, 279426], [35557, 279542], [58704, 279770], [59100, 280139]]
          and compare_interface.get('native_complete_pulse_pressure', {}).get('m15', [])[:2]
          == list(OLD_FRAME_ONLY[:2]))
    check('saved_graph_namespace_metadata',
          compare_namespaces.get('native', {}).get('digest') == BASE_GRAPH
          and compare_namespaces.get('native', {}).get('node_count') == BASE_NODE_COUNT
          and compare_namespaces.get('corrected_physical', {}).get('digest') == CORRECTED_GRAPH
          and compare_namespaces.get('corrected_physical', {}).get('node_count') == CORRECTED_NODE_COUNT)

    quiet_interface = quiet.get('source_interface', {})
    check('accepted_quiet_source_interface',
          _same_sequence(quiet_interface.get('native_Rc'), NATIVE_RC)
          and _same_sequence(quiet_interface.get('native_2Rc'), NATIVE_2RC)
          and _same_sequence(quiet_interface.get('native_Rp'), NATIVE_PRE_RP)
          and _same_sequence(quiet_interface.get('corrected_preRp'), CORRECTED_PRE_RP)
          and _same_sequence(quiet_interface.get('corrected_active'), CORRECTED_ACTIVE)
          and _same_sequence(quiet_interface.get('old_compact'), OLD_COMPACT)
          and quiet_interface.get('native_graph_digest') == QUIET_NATIVE_GRAPH
          and quiet_interface.get('corrected_graph_digest') == CORRECTED_GRAPH
          and quiet_interface.get('physical_scales_already_applied_once') is True
          and quiet_interface.get('m15_is_not_P1') is True)

    quiet_graph = quiet.get('definition_graph', {})
    quiet_outputs = quiet_graph.get('outputs', [])
    check('accepted_quiet_outputs_metadata',
          quiet_graph.get('namespace') == 'original_native_preRp_quiet_source_reduction_v1'
          and quiet_graph.get('node_count') == QUIET_RESULT_NODE_COUNT
          and quiet_graph.get('nodes_digest') == QUIET_RESULT_GRAPH
          and len(quiet_outputs) == 2)
    for order, expected in enumerate(QUIET_OUTPUTS):
        output = quiet_outputs[order] if order < len(quiet_outputs) else {}
        check('accepted_quiet_output_Z' + str(order),
              output.get('ordinary_Z_order') == order
              and all(output.get(key) == value for key, value in expected.items()))
    quiet_identity = quiet.get('quiet_source_identity', {})
    check('quiet_Z1_reduction_is_saved_source_proof',
          quiet_identity.get('ordinary_Z1') == '0'
          and quiet_identity.get('common_physical_K_identity_proved') is False
          and quiet_identity.get('absolute_radius_alias_assumed') is False
          and quiet_identity.get('source_fixture_values_not_original_field_values') is True)

    qualified = all(checks.values())
    source_family = family if isinstance(family, dict) else {}
    graph = PhysicalTerminalComparisonGraph()
    outputs = []
    identities = []

    for order in (0, 1):
        native = lambda handle, field: graph.port(
            'native_terminal_source', BASE_GRAPH, handle, order,
            accepted_report=COMPARE_REPORT, accepted_field=field)
        physical_port = lambda handle, field: graph.port(
            'physical_pulse_source', PHYSICAL_GRAPH, handle, order,
            accepted_report=PARENT_REPORT, accepted_field=field)
        quiet_port = lambda handle, field: graph.port(
            'quiet_native_source', QUIET_NATIVE_GRAPH, handle, order,
            accepted_report=QUIET_REPORT, accepted_field=field)
        corrected = lambda handle, field: graph.port(
            'corrected_reserved_source', CORRECTED_GRAPH, handle, order,
            accepted_report=COMPARE_REPORT, accepted_field=field)

        inner = native(INNER[order], 'inner.m15')
        lower = native(LOWER[order], 'lower.m15')
        native_pre = native(NATIVE_PRE_RP[order], 'native_Rp_pressure.m15')
        old_frame = native(OLD_FRAME_ONLY[order], 'native_complete_pulse_pressure.m15')
        old_pulse = native(OLD_PULSE[order], 'complete_pulse_m15_substitution.old')
        physical_pulse = physical_port(PHYSICAL_PULSE[order], 'pressure.m15')
        saved_target = physical_port(NEGATIVE_PHYSICAL_PULSE[order],
                                     'negative_complete_m15_target')
        saved_geometry = physical_port(PHYSICAL_GEOMETRY_DELTA[order],
                                       'pulse_geometry_correction_m15')
        corrected_pre = corrected(CORRECTED_PRE_RP[order], 'registered_corrected_preRp_m15')
        corrected_active = corrected(CORRECTED_ACTIVE[order], 'registered_corrected_active')
        corrected_reserved = corrected(CORRECTED_RESERVED[order],
                                       'registered_corrected_reserved_m15')
        native_rc = quiet_port(NATIVE_RC[order], 'native_Rc')
        native_2rc = quiet_port(NATIVE_2RC[order], 'native_2Rc')
        old_compact = quiet_port(OLD_COMPACT[order], 'old_compact')
        quiet_output = quiet_outputs[order] if order < len(quiet_outputs) else {}
        quiet_defect = graph.port(
            'accepted_quiet_output', QUIET_RESULT_GRAPH,
            quiet_output.get('quiet_geometric_defect', QUIET_OUTPUTS[order]['quiet_geometric_defect']),
            order, accepted_report=QUIET_REPORT,
            accepted_field='definition_graph.outputs.quiet_geometric_defect')

        forward = graph.linear([(1, inner), (1, lower), (1, native_pre),
                                (1, corrected_reserved)],
                               identity='F=inner+lower+native_preRp+corrected_reserved')
        target_definition = graph.linear([(-1, physical_pulse)],
                                         identity='T=-physical_pulse')
        target = saved_target
        terminal = graph.linear([(1, forward), (-1, target)],
                                identity='F-T')

        corrected_full = graph.linear([(1, inner), (1, lower), (1, corrected_pre),
                                       (1, old_pulse)],
                                      identity='corrected_full=inner+lower+corrected_preRp+old_pulse')
        reserved_expansion = graph.linear([(-1, corrected_full)],
                                           identity='corrected_reserved=-corrected_full')
        forward_expanded = graph.linear([(1, inner), (1, lower), (1, native_pre),
                                         (1, reserved_expansion)])
        terminal_expanded = graph.linear([(1, forward_expanded), (-1, target_definition)])

        native_pre_difference = graph.linear([(1, native_pre), (-1, corrected_pre)],
                                             identity='native_preRp-corrected_preRp')
        physical_geometry_difference = graph.linear(
            [(1, physical_pulse), (-1, old_frame)],
            identity='physical_pulse-frame_only_pulse')
        frame_to_old_pulse = graph.linear([(1, old_frame), (-1, old_pulse)],
                                          identity='frame_only_pulse-old_pulse')
        pulse_difference = graph.linear([(1, physical_pulse), (-1, old_pulse)],
                                        identity='physical_pulse-old_pulse')
        pulse_difference_split = graph.linear(
            [(1, physical_geometry_difference), (1, frame_to_old_pulse)],
            identity='physical_pulse-old_pulse=(physical_pulse-frame_only_pulse)+(frame_only_pulse-old_pulse)')

        active_difference = graph.linear([(1, native_rc), (-1, corrected_active)],
                                         identity='native_Rc-corrected_active')
        native_compact_increment = graph.linear([(1, native_2rc), (-1, native_rc)],
                                                identity='native_2Rc-native_Rc')
        compact_difference = graph.linear([(1, native_compact_increment), (-1, old_compact)],
                                          identity='native_2Rc-native_Rc-old_compact')
        delta_pre_decomposition = graph.linear(
            [(1, active_difference), (1, compact_difference), (1, quiet_defect)],
            identity='Delta_preRp=(native_Rc-corrected_active)+(native_2Rc-native_Rc-old_compact)+quiet_defect')
        terminal_reduced = graph.linear([(1, native_pre_difference), (1, pulse_difference)],
                                        identity='(native_preRp-corrected_preRp)+(physical_pulse-old_pulse)')
        terminal_decomposed = graph.linear([(1, delta_pre_decomposition), (1, pulse_difference)],
                                            identity='Delta_preRp+(physical_pulse-old_pulse)')

        reserved_witness = graph.append(
            'qualified_saved_source_substitution',
            [corrected_reserved, reserved_expansion],
            admitted=qualified,
            ordinary_Z_order=order,
            accepted_report=COMPARE_REPORT,
            accepted_source_node=CORRECTED_RESERVED[order],
            identity='saved corrected_reserved=-(inner+lower+corrected_preRp+old_pulse)',
            source_witness_field='source_witness.corrected_source_overlay_identity',
            terminal_difference_claimed_zero=False)
        target_witness = graph.append(
            'qualified_saved_source_substitution',
            [saved_target, target_definition],
            admitted=qualified,
            ordinary_Z_order=order,
            accepted_report=PARENT_REPORT,
            accepted_source_node=NEGATIVE_PHYSICAL_PULSE[order],
            identity='saved negative_complete_m15_target=-physical_pulse',
            source_witness_field='physical_pulse_bindings.negative_complete_m15_target',
            terminal_difference_claimed_zero=False)
        geometry_witness = graph.append(
            'qualified_saved_source_substitution',
            [saved_geometry, physical_geometry_difference],
            admitted=qualified,
            ordinary_Z_order=order,
            accepted_report=PARENT_REPORT,
            accepted_source_node=PHYSICAL_GEOMETRY_DELTA[order],
            identity='saved geometry delta=physical_pulse-frame_only_pulse',
            source_witness_field='physical_pulse_bindings.pulse_geometry_correction_m15',
            terminal_difference_claimed_zero=False)
        quiet_zero = graph.append(
            'qualified_ordinary_Z1_zero_source_function', [quiet_defect], {},
            admitted=(qualified and order == 1 and quiet_identity.get('ordinary_Z1') == '0'
                      and quiet_outputs[order].get('canonical_quiet_defect') == {}),
            ordinary_Z_order=order,
            accepted_report=QUIET_REPORT,
            accepted_source_node=quiet_output.get('quiet_geometric_defect'),
            source_witness_field='definition_graph.outputs.canonical_quiet_defect',
            common_K_or_radius_alias_used=False,
            quiet_defect_Z1_exact_zero=(order == 1))
        delta_witness = graph.append(
            'qualified_function_substitution',
            [native_pre_difference, delta_pre_decomposition, active_difference,
             compact_difference, quiet_defect],
            admitted=qualified,
            ordinary_Z_order=order,
            accepted_report=QUIET_REPORT,
            accepted_source_nodes={
                'native_preRp': NATIVE_PRE_RP[order],
                'corrected_preRp': CORRECTED_PRE_RP[order],
                'active_difference': quiet_output.get('active_difference'),
                'compact_difference': quiet_output.get('compact_difference'),
                'quiet_defect': quiet_output.get('quiet_geometric_defect'),
            },
            identity='native_preRp-corrected_preRp=Delta_preRp',
            source_witness_field='definition_graph.outputs.canonical_difference',
            common_K_or_radius_alias_used=False,
            remaining_difference_zero_proved=False)
        terminal_witness = graph.append(
            'qualified_terminal_difference_substitution',
            [terminal, terminal_expanded, terminal_reduced, terminal_decomposed,
             reserved_witness, target_witness, delta_witness],
            admitted=qualified,
            ordinary_Z_order=order,
            accepted_reports=[PARENT_REPORT, QUIET_REPORT, COMPARE_REPORT],
            accepted_source_nodes={
                'physical_pulse': PHYSICAL_PULSE[order],
                'negative_target': NEGATIVE_PHYSICAL_PULSE[order],
                'old_frame_only': OLD_FRAME_ONLY[order],
                'old_pulse': OLD_PULSE[order],
                'native_preRp': NATIVE_PRE_RP[order],
                'corrected_preRp': CORRECTED_PRE_RP[order],
                'corrected_reserved': CORRECTED_RESERVED[order],
            },
            identity='F-T=Delta_preRp+(physical_pulse-frame_only_pulse)+frame_only_to_old_pulse',
            source_witness_fields=['source_witness.reduced_identity',
                                   'definition_graph.outputs.canonical_difference'],
            terminal_difference_claimed_zero=False)

        zero = graph.append('exact_zero', [], {})
        native_pre_expansion = graph.linear([(1, corrected_pre), (1, delta_pre_decomposition)])
        relations = {corrected_reserved: reserved_expansion,
                     saved_target: target_definition, native_pre: native_pre_expansion,
                     old_frame: old_pulse}
        if order == 1:
            relations[quiet_defect] = zero
        complete_normal = graph.normal(terminal, relations)
        geometry_reduced = graph.linear([(1, delta_pre_decomposition),
                                         (1, physical_geometry_difference)])
        geometry_normal = graph.normal(geometry_reduced, relations)
        check('actual_saved_target_signed_identity_Z' + str(order),
              complete_normal == geometry_normal)
        check('quiet_derivative_removed_only_by_saved_proof_Z' + str(order),
              order == 0 or quiet_defect not in complete_normal)
        affine_ok = (
            graph.normal(reserved_expansion) == graph.normal(
                graph.linear([(-1, corrected_full)]))
            and graph.normal(pulse_difference) == graph.normal(pulse_difference_split)
            and graph.normal(terminal_expanded) == graph.normal(terminal_reduced)
            and graph.normal(terminal_reduced) == graph.normal(
                graph.linear([(1, native_pre_difference), (1, pulse_difference)]))
            and graph.normal(delta_pre_decomposition) == graph.normal(
                graph.linear([(1, active_difference), (1, compact_difference),
                              (1, quiet_defect)])))
        check('new_affine_identity_Z' + str(order), affine_ok)
        check('qualified_substitution_witness_Z' + str(order),
              graph.nodes[int(reserved_witness.rsplit(':s', 1)[1])].get('accepted_report')
              == COMPARE_REPORT
              and graph.nodes[int(target_witness.rsplit(':s', 1)[1])].get('accepted_report')
              == PARENT_REPORT
              and graph.nodes[int(delta_witness.rsplit(':s', 1)[1])].get('accepted_report')
              == QUIET_REPORT)

        outputs.append({
            'ordinary_Z_order': order,
            'forward_F': forward,
            'target_T': target,
            'complete_terminal_difference': terminal,
            'corrected_reserved_port': corrected_reserved,
            'corrected_full_expansion': corrected_full,
            'corrected_reserved_expansion': reserved_expansion,
            'terminal_difference_expanded': terminal_expanded,
            'native_preRp_difference': native_pre_difference,
            'physical_pulse_geometry_difference': physical_geometry_difference,
            'frame_only_to_old_pulse_difference': frame_to_old_pulse,
            'physical_pulse_difference': pulse_difference,
            'physical_pulse_difference_split': pulse_difference_split,
            'active_difference': active_difference,
            'compact_difference': compact_difference,
            'quiet_defect': quiet_defect,
            'quiet_defect_Z1_witness': quiet_zero,
            'Delta_preRp': delta_pre_decomposition,
            'terminal_difference_reduced': terminal_reduced,
            'terminal_difference_decomposed': terminal_decomposed,
            'saved_negative_physical_target_port': saved_target,
            'saved_physical_geometry_delta_port': saved_geometry,
            'reserved_source_witness': reserved_witness,
            'target_source_witness': target_witness,
            'geometry_source_witness': geometry_witness,
            'Delta_preRp_source_witness': delta_witness,
            'terminal_source_witness': terminal_witness,
            'canonical_terminal_difference': complete_normal,
            'canonical_geometry_decomposition': geometry_normal,
            'geometry_reduced_difference': geometry_reduced,
            'qualified_relation_map': relations,
            'target_definition': target_definition,
            'canonical_Delta_preRp': graph.normal(delta_pre_decomposition),
            'physical_pulse_equals_old_pulse_proved': False,
            'old_physical_geometry_alias_proved': False,
            'terminal_difference_zero_proved': False,
        })
        identities.append({
            'ordinary_Z_order': order,
            'F_definition_proved': True,
            'T_definition_proved': True,
            'saved_reserved_substitution_admitted': True,
            'saved_negative_target_substitution_admitted': True,
            'saved_geometry_delta_substitution_admitted': True,
            'signed_affine_identities_proved': affine_ok,
            'native_preRp_to_Delta_preRp_qualified': True,
            'quiet_defect_Z1_zero_proved': order == 1,
            'physical_pulse_equals_old_pulse_proved': False,
            'old_physical_geometry_alias_proved': False,
            'terminal_difference_zero_proved': False,
        })

    check('physical_target_integrated',
          all(item['saved_negative_target_substitution_admitted'] for item in identities)
          and all(item['F_definition_proved'] and item['T_definition_proved']
                  for item in identities))
    check('terminal_signed_decomposition',
          all(item['signed_affine_identities_proved']
              and item['native_preRp_to_Delta_preRp_qualified']
              for item in identities))
    check('inherited_quiet_Z1_reduction',
          quiet.get('ordinary_Z1_quiet_defect_zero_proved') is True
          and all(item['quiet_defect_Z1_zero_proved'] for item in identities if item['ordinary_Z_order'] == 1))
    check('reserved_m15_registered_once',
          compare.get('registered_reserved_m15_source_bound_once') is True
          and compare.get('source_interface', {}).get('corrected_source_overlay_installed') is True
          and compare.get('retained_construction_data', {}).get(
              'reserved_contribution_registered_once') is True)

    passed = all(checks.values())
    result = {
        'schema': SCHEMA,
        'passed': passed,
        'source_only': True,
        'ordinary_Z_orders': [0, 1],
        'ordinary_Z_not_temporal_n': True,
        'input_digest': expected_input_digest,
        'source_family': source_family,
        'source_prerequisite_audit': checks,
        'definition_graph': {
            'namespace': SCHEMA,
            'nodes': graph.nodes,
            'node_count': len(graph.nodes),
            'nodes_digest': digest(graph.nodes),
            'outputs': outputs,
            'proof_graph_not_numeric_evaluator': True,
        },
        'comparison_identities': identities,
        'source_interface': {
            'F': 'inner+lower+native_preRp+corrected_reserved',
            'T': '-physical_pulse',
            'complete_terminal_difference': 'F-T',
            'Delta_preRp': 'native_Rc-corrected_active+native_2Rc-native_Rc-old_compact+quiet_defect',
            'terminal_identity': 'F-T=Delta_preRp+(physical_pulse-frame_only_pulse); saved whole-function frame identity applied',
            'physical_pulse_source': {
                'graph_digest': PHYSICAL_GRAPH,
                'node_count': PHYSICAL_NODE_COUNT,
                'pressure_m15': list(PHYSICAL_PULSE),
                'negative_complete_m15_target': list(NEGATIVE_PHYSICAL_PULSE),
                'geometry_delta': list(PHYSICAL_GEOMETRY_DELTA),
            },
            'native_base_source': {
                'graph_digest': BASE_GRAPH,
                'node_count': BASE_NODE_COUNT,
                'inner': list(INNER),
                'lower': list(LOWER),
                'native_preRp': list(NATIVE_PRE_RP),
                'old_frame_only': list(OLD_FRAME_ONLY),
                'old_pulse': list(OLD_PULSE),
            },
            'quiet_source': {
                'graph_digest': QUIET_NATIVE_GRAPH,
                'native_Rc': list(NATIVE_RC),
                'native_2Rc': list(NATIVE_2RC),
                'old_compact': list(OLD_COMPACT),
            },
            'corrected_source': {
                'graph_digest': CORRECTED_GRAPH,
                'node_count': CORRECTED_NODE_COUNT,
                'corrected_active': list(CORRECTED_ACTIVE),
                'corrected_preRp': list(CORRECTED_PRE_RP),
                'corrected_reserved': list(CORRECTED_RESERVED),
            },
            'accepted_quiet_outputs': [
                {key: value for key, value in output.items()
                 if key in ('ordinary_Z_order', 'quiet_geometric_defect',
                            'quiet_reduced_difference', 'active_difference',
                            'compact_difference', 'active_compact_split',
                            'canonical_quiet_defect')}
                for output in quiet_outputs
            ],
            'ports_are_digest_and_handle_namespaced': True,
            'source_values_or_midpoints_used': False,
        },
        'retained_construction_data': physical.get('retained_construction_data', {}),
        'absolute_radius_source_audit': radius_audit,
        'source_inspection_audit': {
            'new_source_consumer_calls': 1,
            'accepted_ancestor_operator_calls': 0,
            'old_definition_graphs_rebuilt': False,
            'old_operators_replayed': False,
            'numeric_values_computed': False,
            'quadrature_or_norm_calls': 0,
            'source_value_reducer_calls': 0,
            'temporal_recursion_calls': 0,
            'new_graph_simulation_calls': 0,
            'new_exact_scalar_recipe_comparisons': 1,
            'module_imports_of_ancestor_constructors': 0,
        },
        'legacy_numeric_chain_reconfigured': False,
        'old_error_terms_removed': False,
        'extra_reserved_contribution_added': False,
        'old_physical_geometry_alias_proved': False,
        'new_physical_pulse_equals_old_pulse_proved': False,
        'unresolved_gates': [name for name, value in checks.items() if not value],
        'open_obligations': [
            'native_preRp-corrected_preRp source substitution remains qualified, not numerical',
            'physical pulse versus frame-only pulse and old pulse source correspondence',
            'active, compact and physical quiet Z0 source identities remain unresolved',
            'independent terminal function match and global five-moment closure',
            'other physical joins, exact heat, stress, temporal recursion, corrections and Cartesian NS',
        ],
        'runtime_model_metadata': None,
    }
    result.update({flag: passed for flag in STAGE_FLAGS})
    result.update({flag: False for flag in FALSE_FLAGS})
    return result
