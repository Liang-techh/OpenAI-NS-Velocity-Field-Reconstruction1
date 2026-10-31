"""Instantiate archived Q4 endpoint source ports for a chosen positive sc.

The caller authenticates the receipts before passing them here.  Ports retain
the original physical U_y expression and carry explicit parameter and radial
coordinate substitutions; they are not numerical field values.  A resolver
must independently authenticate the named archive and honor both bindings.
"""
from collections.abc import Mapping

Q4 = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz'
LOOP = 'experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz'
SCHEMA = 'original_q4_sc_flat_instance_v1'


class Missing(ValueError):
    pass


class Mismatch(ValueError):
    pass


def field(obj, name):
    if not isinstance(obj, Mapping) or name not in obj:
        raise Missing(name)
    return obj[name]


def dyadic(value):
    if not isinstance(value, list) or len(value) != 4 or any(type(n) is not int for n in value):
        raise Mismatch('invalid dyadic tuple')
    sign, mantissa, exponent, bits = value
    if sign != 0 or mantissa <= 0 or mantissa % 2 != 1 or bits != mantissa.bit_length():
        raise Mismatch('noncanonical positive dyadic')
    if bits + exponent > -3:
        raise Mismatch('sc not certified less than one eighth')
    return list(value)


def port_valid(port):
    if not isinstance(port, Mapping) or port.get('schema') != 'q4_physical_Uy_endpoint_source_port_v1':
        raise Mismatch('invalid source port')
    archive = field(port, 'archive_sha256')
    if not isinstance(archive, str) or len(archive) != 64 or any(c not in '0123456789abcdef' for c in archive):
        raise Mismatch('invalid archive hash')
    for key in ('expression_node', 'parameter_node', 'coordinate_node', 'ordinary_Z_order'):
        if type(field(port, key)) is not int or port[key] < 0:
            raise Mismatch('invalid ' + key)
    if port['ordinary_Z_order'] not in (0, 1, 2):
        raise Mismatch('unsupported ordinary Z order')
    sc = dyadic(field(port, 'parameter_exact_mpf_tuple'))
    half = dyadic(field(port, 'coordinate_exact_mpf_tuple'))
    if half != [sc[0], sc[1], sc[2] - 1, sc[3]]:
        raise Mismatch('coordinate is not exact sc/2')
    return True


def evaluate_port(port, resolver):
    """Resolve an endpoint port once using an externally supplied source evaluator.

    The resolver receives a copy of the entire port, including both exact
    substitutions.  It owns archive authentication and field evaluation;
    this adapter does not manufacture values or collapse background history.
    """
    port_valid(port)
    if not callable(resolver):
        raise TypeError('source resolver must be callable')
    copied = dict(port)
    copied['parameter_exact_mpf_tuple'] = list(port['parameter_exact_mpf_tuple'])
    copied['coordinate_exact_mpf_tuple'] = list(port['coordinate_exact_mpf_tuple'])
    return resolver(copied)


def join(instance_bundle, flat_bundle):
    checks, errors = {}, []

    def check(label, condition):
        checks[label] = bool(condition)
        if not condition:
            raise Mismatch(label)

    try:
        ir = field(instance_bundle, 'result')
        fr = field(flat_bundle, 'result')
        check('both_executions_returned', field(ir, 'status') == 'returned' and field(fr, 'status') == 'returned')
        inst = field(ir, 'instance_result')
        flat = field(fr, 'archive_result')
        check('both_results_verified', field(inst, 'status') == 'verified' and field(inst, 'passed') is True
              and field(flat, 'status') == 'verified' and field(flat, 'passed') is True)
        iscope, fscope = field(inst, 'scope'), field(flat, 'scope')
        if not isinstance(iscope, Mapping) or not isinstance(fscope, Mapping):
            raise Missing('receipt scope')
        check('instance_scope_preserved', iscope.get('newly_chosen_binding') is True
              and iscope.get('binding_materialized') is True
              and iscope.get('existing_archive_assignment_proved') is False
              and iscope.get('corrected_owner_identity_proved') is False
              and iscope.get('old_errors_removed') is False
              and iscope.get('flat_conclusion_proved_here') is False)
        check('flat_scope_preserved', fscope.get('Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive') is True
              and fscope.get('sc_positive_independently_certified') is False
              and fscope.get('old_errors_removed') is False
              and fscope.get('corrected_source_identity_proved') is False
              and fscope.get('ordinary_Z_orders') == [0, 1, 2])
        binding = field(inst, 'binding')
        graph = field(binding, 'graph_id')
        hashes = field(fr, 'input_hashes')
        q4, loop = field(graph, 'q4'), field(graph, 'loop')
        check('archive_identity_matches', field(hashes, Q4) == q4 and field(hashes, LOOP) == loop)
        node = field(binding, 'sc_node_id')
        check('selected_parameter_node_is_valid', type(node) is int and node >= 0
              and field(binding, 'node_id') == node)
        sc = dyadic(field(binding, 'sc_exact_mpf_tuple'))
        half = [sc[0], sc[1], sc[2] - 1, sc[3]]
        derived = field(inst, 'derived_dyadics')
        check('derived_half_matches_binding', field(field(derived, 'sc_over_two'), 'exact_mpf_tuple') == half)
        comparison = field(derived, 'comparison')
        if not isinstance(comparison, Mapping):
            raise Missing('dyadic comparison')
        check('prior_exact_collar_comparison', comparison.get('sc_over_two_lt_three_sc_over_four') is True
              and comparison.get('three_sc_over_four_lt_sc') is True)
        evidence = field(flat, 'evidence')
        coordinate = field(evidence, 'coordinate')
        physical = field(evidence, 'physical_U_y_rows')
        background = field(evidence, 'background_V_y_rows')
        pstar = field(evidence, 'Pstar_node_from_prior_receipt')
        check('source_rows_are_three_exact_references',
              type(coordinate) is int and coordinate >= 0 and type(pstar) is int and pstar >= 0
              and isinstance(physical, list) and len(physical) == 3
              and isinstance(background, list) and len(background) == 3
              and all(type(x) is int and x >= 0 for x in physical + background))
        ports = []
        for k, ref in enumerate(physical):
            port = {'schema': 'q4_physical_Uy_endpoint_source_port_v1',
                    'archive_sha256': q4, 'expression_node': ref, 'ordinary_Z_order': k,
                    'parameter_node': node, 'parameter_exact_mpf_tuple': list(sc),
                    'coordinate_node': coordinate, 'coordinate_exact_mpf_tuple': list(half)}
            port_valid(port)
            ports.append(port)
        evidence_out = {'parameter_node': node, 'coordinate_node': coordinate,
                        'physical_U_y_rows': physical, 'background_V_y_rows': background,
                        'Pstar_node': pstar, 'sc_exact_mpf_tuple': sc,
                        'endpoint_exact_mpf_tuple': half}
        status = 'verified'
    except Missing as exc:
        status, errors, evidence_out, ports = 'unresolved', [str(exc)], {}, []
    except Mismatch as exc:
        status, errors, evidence_out, ports = 'mismatch', [str(exc)], {}, []
    return {
        'schema': SCHEMA, 'status': status, 'passed': status == 'verified',
        'checks': checks, 'errors': errors, 'evidence': evidence_out,
        'source_ports': {'U_y': ports},
        'scope': {
            'newly_chosen_binding': status == 'verified',
            'flat_Uy_correction_zero_for_new_instance': status == 'verified',
            'other_archived_source_assumptions_required': True,
            'background_Uy_retained': status == 'verified',
            'existing_archive_assignment_proved': False,
            'corrected_owner_identity_proved': False,
            'history_function_equality_claimed': False,
            'global_five_moment_closure_claimed': False,
            'numeric_field_values_computed': False,
            'old_errors_removed': False,
            'ordinary_Z_orders': [0, 1, 2],
            'source_evaluator_not_installed': True,
        },
    }
