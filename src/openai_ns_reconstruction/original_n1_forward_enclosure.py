"""Directed full raw n1 enclosures from the accepted coefficient archive.

No coefficient is replaced by its midpoint.  A polynomial interval contains
the exact source coefficient polynomial; a separate symmetric analytic tail
contains the omitted radial series.  The positive F0 scale can be too small
to materialize.  In that case its logarithm proves a sign-aware polynomial
bound, while the exact factored enclosure remains available to consumers.
These are coefficient-field enclosures, not a completed physical NS field.
"""
from __future__ import annotations

from fractions import Fraction
import ast
import hashlib
import math
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext

from .original_finite_N_signed_source import load, sha, digest, require
from .original_n1_complex_source_bounds import _interval, _record, _serialization_roundtrip
from .original_n1_forward_analytic import _digest as _analytic_digest

PARENT = 'lei_ren_part1_paper_compliant_current_original_n1_canonical_norm_tail.json'
PARENT_SHA256 = 'b71ee21cf0f1dd2f3fe1d4edf4f45037ced997f1d72e2b8e0c9f1c861aae6989'
FIELDS = ('F1', 'Uz1', 'K1', 'P1', 'F1x', 'Uz1x')
NORMALIZED = ('f', 'u', 'k', 'p', 'f', 'u')


def collar_dependency_scope(root, parent):
    """Inspect the actual consumer type before attempting a response bridge."""
    root=Path(root);base=root/'experiments/root_st073'
    prefix='lei_ren_part1_paper_compliant_current_original_n1_'
    names=dict(source=prefix+'lower_collar_functions.py',
        lower=prefix+'lower_collar_functions.json.gz',
        volterra=prefix+'inner_volterra_functions.json.gz',
        radial=prefix+'inner_radial_functions.json.gz')
    for name in names.values():
        require(sha(base/name)==parent['input_hashes'][name], 'Changed collar dependency source: '+name)
    tree=ast.parse((base/names['source']).read_text(encoding='utf-8'))
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='collar_fields')
    def keys(variable):
        return sorted({n.slice.value for n in ast.walk(function) if isinstance(n,ast.Subscript)
            and isinstance(n.value,ast.Name) and n.value.id==variable
            and isinstance(n.slice,ast.Constant) and isinstance(n.slice.value,str)})
    rows,functions=keys('core_rows'),keys('core_functions')
    require(rows==['Phi','Psi','mean_Psi'] and functions==['F0','P0','U0'],
        'The lower collar input type changed')
    lower=load(base/names['lower']);volterra=load(base/names['volterra']);radial=load(base/names['radial'])
    ln=lower['exact_graph_nodes'];vn=volterra['exact_graph_nodes'];rn=radial['exact_graph_nodes']
    raw_nodes=list(range(43300,43324));pressure_nodes=list(range(52128,52132))
    require(len(ln)==36594 and len(vn)==49024 and len(rn)==52487 and
        max(range(len(ln)))<min(raw_nodes+pressure_nodes), 'Saved graph stage ordering changed')
    leading=[]
    for quantity,node in (('Phi',35674),('Psi',35678),('mean_Psi',35682)):
        require(lower['lower_collar_bindings']['core_radial_rows'][quantity][0][0]==node and
            ln[node]['operation']=='current_original_unique_analytic_core_component' and
            ln[node]['quantity']==quantity and ln[node]['source_family']==parent['source_family'],
            'Actual leading collar profile source changed')
        leading.append(dict(quantity=quantity,node=node,operation=ln[node]['operation']))
    for node in raw_nodes:
        require(vn[node]['operation']=='current_original_unique_regular_n1_component',
            'Actual Volterra raw response nodes changed')
    require(radial['inner_radial_bindings']['pressure_rho']['core']==pressure_nodes,
        'Radial core pressure-rho binding changed')
    for offset,node in enumerate(pressure_nodes):
        require(rn[node]==dict(operation='product',arguments=[52124+offset,35565]) and
            rn[52124+offset]==dict(operation='sum',arguments=[52120+offset,(36809,37964,38138,38431)[offset]]),
            'Radial pressure-rho expression changed')
    source_hash=sha(base/names['source'])
    return dict(passed=True,source_pins={n:sha(base/n) for n in names.values()},
        source_AST_sha256=hashlib.sha256(ast.dump(function,include_attributes=False).encode()).hexdigest(),
        source_module_sha256=source_hash,actual_leading_profile_dependencies=rows,
        actual_core_function_dependencies=functions,leading_profile_nodes=leading,
        saved_graph_node_counts=dict(lower=len(ln),volterra=len(vn),radial=len(rn)),
        raw_response_nodes=raw_nodes,radial_pressure_rho_nodes=pressure_nodes,
        raw_response_nodes_absent_from_lower_graph=True,
        pressure_rho_nodes_absent_from_lower_graph=True,radial_pressure_rho_binding_and_shape_checked=True,
        lower_consumer_kind='leading profile continuation and its differential N1p forcing',
        raw_n1_response_is_not_a_leading_Phi_Psi_input=True,
        actual_raw_response_to_collar_inlet_identity=False,
        actual_q_dependent_leading_collar_cross_producer_identity=False,
        geometry_alone_proves_consumer_identity=False,
        six_signed_source_equality=False,source_or_owner_executed=False,
        next_action='Bind the actual forward leading collar to Phi/Psi continuation; separately use the infinite core n1 pressure equation for the core signed source identity. Do not feed raw n1 rows into leading profile inputs.')


def _rational(value):
    require(type(value) in (str, int, Fraction), 'Use exact rational data, not binary floats')
    return Fraction(value)


def _rho_box(ctx, value):
    values = value if isinstance(value, (list, tuple)) else (value, value)
    require(len(values) == 2, 'Two rational rho endpoints required')
    lo, hi = map(_rational, values)
    require(0 <= lo <= hi <= 4, 'Real rho box must lie in [0,4]')
    left = ctx.mpf(lo.numerator) / lo.denominator
    right = ctx.mpf(hi.numerator) / hi.denominator
    return ctx.make_mpf((left._mpi_[0], right._mpi_[1])), lo, hi


def _point_upper(ctx, value):
    return ctx.make_mpf((value._mpi_[1], value._mpi_[1]))


def _index(parent):
    prefix = parent['n1_actual_raw_prefix_binding']['actual_forward_prefix_bound']
    require(prefix['passed'] is True and
        prefix['ordinary_Z_domain']['max_rho_degree_with_orders_0_to_3'] == 14,
        'Accepted common ordinary-Z degree14 prefix required')
    result = {}
    for row in prefix['coefficient_enclosures']:
        key = (row['Z'], row['field'], row['rho_degree'], row['Z_degree'])
        require(key not in result, 'Duplicate source coefficient slot')
        result[key] = row['source_enclosure_exact_mpf_tuples']
    return result


def _coefficient(ctx, index, center, field, degree, order):
    key = (center, field, degree, order)
    require(key in index, 'Unavailable source coefficient: ' + str(key))
    endpoints = tuple(tuple(int(v) for v in endpoint) for endpoint in index[key])
    return ctx.make_mpf(endpoints) * math.factorial(order)


def _polynomial_jets(ctx, index, center, field, rho, odd):
    result = []
    for j in range(4):
        terms = []
        for m in range(1, 15):
            coefficient = _coefficient(ctx, index, center, field, m, j)
            terms.append(coefficient * (2*m if odd else 1) * rho**(m-1 if odd else m))
        result.append(sum(terms, ctx.mpf(0)))
    return result


def _scaled_polynomial(ctx, logscale, factor, tail_upper, materialize):
    """Enclose exp(logscale)*factor without expanding a huge MPF exponent."""
    if factor.a == 0 and factor.b == 0:
        return ctx.mpf(0), dict(method='exact_zero_factor', exponential_materialized=False)
    if materialize:
        return ctx.exp(logscale)*factor, dict(method='directed_exp_of_reduced_log',
            exponential_materialized=True)
    amplitude = _point_upper(ctx, abs(factor))
    require(amplitude.a > 0 and tail_upper.a > 0, 'Positive comparison scale required')
    log_amplitude = logscale + ctx.log(amplitude)
    log_tail = ctx.log(tail_upper)
    require(log_amplitude.b <= log_tail.a,
        'Unmaterialized polynomial is not dominated by the certified tail scale')
    if factor.b <= 0:
        box = ctx.make_mpf(((-tail_upper)._mpi_[0], ctx.mpf(0)._mpi_[1]))
    elif factor.a >= 0:
        box = ctx.make_mpf((ctx.mpf(0)._mpi_[0], tail_upper._mpi_[1]))
    else:
        box = ctx.make_mpf(((-tail_upper)._mpi_[0], tail_upper._mpi_[1]))
    return box, dict(method='logarithmic_dominance_with_original_factor_sign',
        exponential_materialized=False, log_polynomial_absolute_upper=_record(log_amplitude),
        log_comparison_radius=_record(log_tail), domination_proved=True,
        comparison_radius_is_bound_not_value=True)


def evaluate(parent, forward, scales, center, rho_box, *, precision=160):
    """Evaluate a full raw six-field ordinary-Z0..3 inclusion at an archive center.

    ``scales`` is the same-source scale-jet packet, never arbitrary input
    coefficients.  This evaluator intentionally has no arbitrary-Z oracle.
    """
    require(type(precision) is int and precision >= 100, 'At least100 decimal digits required')
    require(forward.get('passed') is True and scales.get('passed') is True,
        'Actual analytic forward and scale bindings required')
    require(forward.get('parent_report_sha256') == PARENT_SHA256 and
        scales.get('parent_report_sha256') == PARENT_SHA256, 'Different source parent')
    require(all(forward.get(key) is True for key in ('actual_infinite_forward_analytic_response',
        'new_infinite_forward_definition','actual_forward_analytic_definition_bound',
        'actual_forward_regular_uniqueness','actual_forward_prefix_identity','actual_source_scale_conjugacy')),
        'Actual infinite analytic forward response, rather than a finite prefix, required')
    accepted_digests = forward['accepted_packet_digests']
    require(accepted_digests['norm_tail'] == _analytic_digest(parent['n1_actual_canonical_norm_tail']) and
        accepted_digests['raw_prefix_binding'] == _analytic_digest(parent['n1_actual_raw_prefix_binding']) and
        scales['source_family'] == parent['source_family'],
        'Changed canonical raw tail, coefficient archive, or scale source family')
    tail = parent['n1_actual_canonical_norm_tail']
    require(tail['passed'] is True and tail['actual_canonical_raw_radial_Taylor_remainder_certified'] is True
        and tail['raw_base_prefix_degree'] == 14, 'Accepted raw analytic tail required')
    centers = parent['n1_actual_raw_prefix_binding']['actual_forward_prefix_bound'][
        'ordinary_Z_domain']['coefficient_centers']
    require(type(center) is str and center in centers, 'Use an exact archived real Z center')
    ctx = MPIntervalContext(); ctx.dps = precision
    rho, lo, hi = _rho_box(ctx, rho_box)
    rmax = ctx.mpf(hi.numerator)/hi.denominator
    q = rmax/(ctx.mpf(9)/2)
    norm = _interval(ctx, tail['weighted_solution_norm_upper_majorant'])
    gap = _interval(ctx, tail['axial_domain']['fresh_Cauchy_gap'])
    require(gap.a > 0 and q.a >= 0 and q.b < 1, 'Positive Cauchy gap and radial ratio below1 required')
    weighted_tail = norm*q**14/(1-q)
    index = _index(parent)
    center_scales = scales['centers'][center]['raw_scales']
    rows = []
    for i, (field, normalized) in enumerate(zip(FIELDS, NORMALIZED)):
        odd = i >= 4
        polynomials = _polynomial_jets(ctx, index, center, normalized, rho, odd)
        source_scale = center_scales[normalized]
        logscale = _interval(ctx, source_scale['log_positive_scale'])
        require(source_scale['relative_jet_convention'] == 'D_Z^j(scale)/scale; no factorial divisor',
            'Relative ordinary derivative jets, not relative Taylor coefficients, required')
        relative = [_interval(ctx, v) for v in source_scale['relative_derivative_Z_jets']]
        require(len(relative) == 4 and relative[0].a == 1 and relative[0].b == 1,
            'Actual ordinary scale Bell jets required')
        for j in range(4):
            factor = sum((math.comb(j,k)*relative[k]*polynomials[j-k]
                for k in range(j+1)), ctx.mpf(0))
            if odd:
                factor *= ctx.sqrt(rho)
            radius = weighted_tail * (ctx.sqrt(rmax) if odd else rmax) * math.factorial(j)/gap**j
            radius_upper = _point_upper(ctx, radius)
            polynomial_box, method = _scaled_polynomial(ctx, logscale, factor, radius_upper,
                source_scale['exponential_materialization_safe'])
            symmetric = ctx.make_mpf(((-radius_upper)._mpi_[0], radius_upper._mpi_[1]))
            full = polynomial_box + symmetric
            rows.append(dict(field=field, normalized_source_field=normalized, ordinary_Z_order=j,
                normalized_polynomial_ordinary_jet=_record(polynomials[j]),
                raw_polynomial_factored_enclosure=dict(log_positive_scale=_record(logscale),
                    signed_factor=_record(factor), expression='exp(log_positive_scale)*signed_factor'),
                raw_polynomial_enclosure=_record(polynomial_box),
                polynomial_enclosure_method=method,
                analytic_radial_tail_upper_majorant=_record(radius),
                full_raw_function_enclosure=_record(full),
                exact_source_coefficient_enclosure_errors_included=True,
                derivative_factor_2m_applied=odd))
    result = dict(passed=True, center=center, rho_box=[str(lo),str(hi)], common_raw_degree=14,
        ordinary_Z_orders=[0,1,2,3], rows=rows,
        coefficient_archive_center_only=True, coefficient_midpoint_selected=False,
        coefficient_decimal_approximation_selected=False, enclosure_is_function_value=False,
        interval_encloses_actual_infinite_raw_function=True,
        physical_R_derivatives_or_NS_velocity_evaluated=False,
        new_coefficient_recurrence_or_word_sum_evaluated=False,
        all_old_error_caps_and_six_signed_terms_retained=True)
    result['exact_interval_records_roundtripped'] = _serialization_roundtrip(ctx, result)
    return result


def build(root, parent, forward, scales):
    root = Path(root); path = root/'experiments/root_st073'/PARENT
    require(sha(path) == PARENT_SHA256 and digest(parent) == digest(load(path)),
        'Accepted canonical norm/tail parent changed')
    evaluations = [evaluate(parent, forward, scales, center, ['0','4'])
        for center in ('-0.3','0','0.3')]
    exits = [evaluate(parent, forward, scales, center, '4')
        for center in ('-0.3','0','0.3')]
    axes = [evaluate(parent, forward, scales, center, '0')
        for center in ('-0.3','0','0.3')]
    for axis in axes:
        for row in axis['rows']:
            box = row['full_raw_function_enclosure']
            require(all(int(v) == 0 for v in box['lower_exact_mpf_tuple']) and
                all(int(v) == 0 for v in box['upper_exact_mpf_tuple']), 'Axis must be exact zero')
    ctx = MPIntervalContext(); ctx.dps = 160
    for whole, endpoint in zip(evaluations, exits):
        for outer, inner in zip(whole['rows'], endpoint['rows']):
            a = _interval(ctx, outer['full_raw_function_enclosure'])
            b = _interval(ctx, inner['full_raw_function_enclosure'])
            require(a.a <= b.a and b.b <= a.b, 'Endpoint enclosure escaped whole radial core')
    return dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        forward_binding_sha256=digest(forward), scale_binding_sha256=digest(scales),
        accepted_coefficient_packet_sha256=digest(parent['n1_actual_raw_prefix_binding']),
        accepted_raw_tail_packet_sha256=digest(parent['n1_actual_canonical_norm_tail']),
        exact_coefficient_polynomial_and_analytic_tail_combined=True,
        whole_real_radial_core_at_archive_centers= evaluations,
        core_exit_at_archive_centers=exits,
        endpoint_interval_in_whole_radial_box=True,
        exact_axis_checks=axes,
        full_raw_function_enclosures=True, arbitrary_real_Z_query_or_global_matching=False,
        physical_R_derivative_enclosure=False, physical_NS_velocity_enclosure=False,
        collar_six_source_equality_or_global_moment_closure=False,
        old_error_caps_modified=False, temporal_n_ge_2=False,
        accepted_owner_builder_fixture_recurrence_or_quadrature_replayed=False,
        limitations=['The new infinite analytic forward definition is matched to the coefficient archive, not to the still-unbound collar producer',
            'These are interval enclosures at three exact archived centers, not all-Z moment identities',
            'Logarithmic polynomial dominance is an inclusion bound, never an actual signed function value',
            'Physical R derivatives and the final time-dependent Cartesian velocity remain outside this packet'])
