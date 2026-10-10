"""Exact-scale and ordinary numeric error budgets for original physical rows.

An unevaluated exact positive scale can multiply a certified signed interval.
Its factored relative diameter differs from that of an ordinary numeric box.
Both states remain explicit; a failed enclosure test is not actual field error.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rp_pressure_binding as pressure
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

signed, correlated = pressure.signed, pressure.correlated
HERE, PREFIX, sha = pressure.HERE, pressure.PREFIX, pressure.sha
NAME = PREFIX+'current_original_Rp_physical_accuracy.json.gz'
RECEIPT = PREFIX+'current_original_Rp_physical_accuracy_check.json'
GATES = ('current_original_Rp_exact_scale_factored_relative_error_installed',
         'current_original_Rp_ordinary_scale_log_width_error_budget_installed',
         'current_original_Rp_exact_affine_origin_and_source_width_ledger_installed')
OPEN = pressure.OPEN
ends = pressure.ends


def add_poly(left, right, multiplier=Fraction(1)):
    result = dict(left)
    for atoms, value in right.items():
        result[atoms] = result.get(atoms, Fraction(0))+multiplier*value
    return {atoms: value for atoms, value in result.items() if value}


def mul_poly(left, right):
    result = {}
    for a, av in left.items():
        for b, bv in right.items():
            atoms = tuple(sorted(a+b))
            result[atoms] = result.get(atoms, Fraction(0))+av*bv
    return {atoms: value for atoms, value in result.items() if value}


def affine_origin(polynomial, origin, pivot):
    """Exact change of basis around logRp, using its unique 10*logC term.

    This is algebra on source functions, never a regression on their bounds.
    The coefficient may depend on other exact source atoms, e.g. delta.
    """
    q = origin.get((pivot,))
    if q != Fraction(10) or any(pivot in atoms and atoms != (pivot,) for atoms in origin):
        raise ValueError('Actual source origin must have its unique exact 10*logC term')
    coefficient = {}
    for atoms, value in polynomial.items():
        count = atoms.count(pivot)
        if count > 1:
            raise ValueError('Original row must remain affine in logC')
        if count:
            reduced = list(atoms)
            reduced.remove(pivot)
            coefficient[tuple(reduced)] = value/q
    residual = add_poly(polynomial, mul_poly(coefficient, origin), Fraction(-1))
    if any(pivot in atoms for atoms in residual):
        raise ValueError('Exact logC dependency must leave the residual')
    if add_poly(add_poly(mul_poly(coefficient, origin), residual), polynomial, Fraction(-1)):
        raise ValueError('Exact affine source identity required')
    return coefficient, residual


@source_precision
def relative_error_budget(ctx, row, target):
    """Two sufficient interval-width tests; neither chooses a nominal value."""
    target = Fraction(target)
    if not 0 < target < 1:
        raise ValueError('Exact relative width target in (0,1) required')
    eps = ctx.mpf(target.numerator)/target.denominator
    out = dict(relative_width_target=signed.rational_record(target),
        diameter_convention='(maximum magnitude-minimum magnitude)/minimum magnitude',
        exact_reference_scale_retained=True, coefficient_or_scale_nominal_value_selected=False,
        factored_relative_width_satisfied=False, ordinary_numeric_relative_width_satisfied=False,
        ordinary_numeric_materialized=row['ordinary_numeric_materialized'],
        ordinary_numeric_delivery_target_satisfied=False,
        failure_does_not_prove_true_field_error=True)
    if row['exact_zero_enclosure']:
        return dict(out, exact_zero_enclosure=True, sign=0, factored_relative_width_satisfied=True,
            ordinary_numeric_relative_width_satisfied=True, coefficient_relative_width_upper=ctx.mpf(0),
            ordinary_numeric_delivery_target_satisfied=row['ordinary_numeric_materialized'],
            log_ordinary_magnitude_ratio_bound=ctx.mpf(0),
            absolute_error_representation=dict(kind='exact_zero', bound=ctx.mpf(0)))
    common = row.get('common_scale_coefficient_enclosure')
    if common is None:
        return dict(out, exact_zero_enclosure=False, sign=None,
            unresolved_reason='Distinct scale ratios are unresolved; keep original signed groups')
    lo, hi = ends(common)
    width = ctx.mpf(hi)-ctx.mpf(lo)
    scale = ctx.mpf(row['directed_reference_log_scale'])
    scale_width = ctx.mpf(ends(scale)[1])-ctx.mpf(ends(scale)[0])
    out.update(exact_zero_enclosure=False, common_coefficient_enclosure=common,
        coefficient_absolute_width_bound=width, reference_log_width_bound=scale_width,
        absolute_error_representation=dict(kind='exact_positive_scale_times_coefficient_width',
            exact_scale_log_function=row['exact_reference_log_scale_function'],
            coefficient_width_upper=ctx.mpf(ends(width)[1]),
            coefficient_width_is_not_absolute_physical_error=True),
        all_signed_ratio_and_positive_tail_errors_already_in_coefficient=True)
    if lo <= 0 <= hi:
        return dict(out, sign=None, unresolved_reason='Common signed coefficient includes zero; relative accuracy needs sign/nonzero recovery')
    minimum, maximum = min(abs(lo), abs(hi)), max(abs(lo), abs(hi))
    relative = width/ctx.mpf(minimum)
    coefficient_log_ratio = ctx.ln(ctx.mpf(maximum)/ctx.mpf(minimum))
    ordinary_log_ratio = scale_width+coefficient_log_ratio
    allowed_log_ratio = ctx.ln(1+eps)
    out.update(sign=1 if lo > 0 else -1, coefficient_relative_width_upper=relative,
        log_coefficient_magnitude_ratio_bound=coefficient_log_ratio,
        log_ordinary_magnitude_ratio_bound=ordinary_log_ratio,
        log_allowed_magnitude_ratio_bound=allowed_log_ratio,
        factored_relative_width_satisfied=ends(relative)[1] <= ends(eps)[0],
        ordinary_numeric_relative_width_satisfied=ends(ordinary_log_ratio)[1] <= ends(allowed_log_ratio)[0],
        ordinary_error_test_never_expands_common_scale=True,
        exact_factored_certificate_does_not_certify_numeric_scale_realization=True)
    if ends(ordinary_log_ratio)[1] <= signed.EXP_LOG_LIMIT:
        out['ordinary_relative_width_bound'] = ctx.exp(ordinary_log_ratio)-1
    else:
        out['ordinary_relative_width_upper_representation'] = dict(
            kind='expm1_of_directed_log_ratio_bound', argument=ordinary_log_ratio,
            physical_exponential_not_materialized=True)
    if not out['ordinary_numeric_relative_width_satisfied']:
        out['ordinary_unresolved_reason'] = 'Current scale/coefficient box does not certify requested ordinary relative width; no claim about true point error'
    if not out['factored_relative_width_satisfied']:
        out['factored_unresolved_reason'] = 'Current signed coefficient/ratio enclosure exceeds requested factored relative width'
    out['ordinary_numeric_delivery_target_satisfied'] = (out['ordinary_numeric_relative_width_satisfied']
        and row['ordinary_numeric_materialized'] and row.get('requested_relative_width_satisfied', False))
    return out


class CurrentOriginalRpPhysicalAccuracy:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not pressure.CurrentOriginalRpPressureBinding or not before.acceptance_loaded:
            raise ValueError('Accepted current physical/P0 source owner required')
        self.before, self.graph, self.family_record = before, before.graph, before.family_record
        self.ctx = MPIntervalContext()
        self.ctx.dps = 400
        self.source_ctx = before.ctx
        self.hashes = dict(before.hashes)
        for name in (pressure.NAME, pressure.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded = False
        self.origin = before.amplitude.radius.logRp
        self.pivot = before.amplitude.radius.parameters['logC']
        self._dependency_cache = {}
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family'] != self.family_record:
                raise ValueError('Current physical accuracy receipt or scope differs')
            correlated.box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True

    def assert_graph(self):
        self.before.assert_graph()
        checks = dict(accepted_current_P0_source=self.before.acceptance_loaded,
            same_source_graph_and_family=self.graph is self.before.graph and self.family_record == self.before.family_record,
            same_original_radius_origin=self.origin == self.before.amplitude.radius.logRp,
            same_exact_original_logC=self.pivot == self.before.amplitude.radius.parameters['logC'],
            unchanged_error_arithmetic_precision=self.ctx.dps == 400)
        if not all(checks.values()):
            raise ValueError('Current physical accuracy source differs')
        return checks

    def dependencies(self, reader, node):
        """Small bound-leaf ledger; never serialize the whole source recipe."""
        node = node.node if hasattr(node, 'node') else node
        key = (node, tuple(sorted(reader.bindings)), tuple(sorted(reader.aliases.items())))
        if key in self._dependency_cache:
            return self._dependency_cache[key]
        if node in reader.aliases:
            result = self.dependencies(reader, reader.aliases[node])
        elif node in reader.bindings:
            result = (node,)
        else:
            n = self.graph.nodes[node]
            op = n['operation']
            children = n['arguments'] if op in ('sum', 'product') else \
                (n['argument'],) if op in ('negative', 'analytic_unary') else \
                (n['numerator'], n['denominator']) if op == 'positive_quotient' else ()
            result = tuple(sorted(set(v for child in children for v in self.dependencies(reader, child))))
        self._dependency_cache[key] = result
        return result

    def ledger(self, row, reader):
        g, c = self.graph, self.ctx
        full = reader.polynomial(row['exact_reference_log_scale_function'])
        origin = reader.polynomial(self.origin)
        coefficient, residual = affine_origin(full, origin, self.pivot.node)
        af = signed.polynomial_function(g, coefficient)
        bf = signed.polynomial_function(g, residual)
        reconstructed = g.add(g.mul(af, self.origin), bf)
        original = correlated.box.pulse.radius.FunctionRef(g, row['exact_reference_log_scale_function'])
        if reader.polynomial(g.sub(reconstructed, original)):
            raise ValueError('Actual current full scale must equal the affine source form')
        proof = g.node('current_original_Rp_exact_affine_common_log_scale',
            original_full_scale_log=row['exact_reference_log_scale_function'],
            original_logRp=self.origin.node, exact_affine_coefficient=af.node,
            exact_residual_function=bf.node, exact_logC_pivot=self.pivot.node,
            source_family=self.family_record, identity='Lref=A*original_logRp+B',
            no_source_parameter_or_interval_endpoint_substituted=True)
        terms = []
        for atoms, value in sorted(full.items()):
            term = {(atoms): value}
            function = signed.polynomial_function(g, term)
            bound = c.mpf(reader.at(function))
            width = c.mpf(ends(bound)[1])-c.mpf(ends(bound)[0])
            leaves = sorted(set(v for atom in atoms for v in self.dependencies(reader, atom)))
            terms.append(dict(exact_term_function=function.node, exact_atoms=list(atoms),
                exact_rational_coefficient=signed.rational_record(value), directed_term=bound,
                directed_term_width=width, atoms=[dict(node=atom,
                    operation=g.nodes[atom]['operation'], name=g.nodes[atom].get('name')) for atom in atoms],
                original_bound_leaves=[dict(node=node, operation=g.nodes[node]['operation'],
                    name=g.nodes[node].get('name'), directed_source_bound=reader.bindings[node]) for node in leaves]))
        dominant = max(terms, key=lambda term: ends(term['directed_term_width'])[1]) if terms else None
        return dict(original_reference_log_scale=row['exact_reference_log_scale_function'],
            exact_affine_logRp_identity=proof.node, exact_affine_coefficient=af.node,
            exact_origin_function=self.origin.node, exact_residual_function=bf.node,
            directed_affine_coefficient=reader.at(af), directed_original_logRp=reader.at(self.origin),
            directed_exact_residual=reader.at(bf), exact_origin_logC_removed_from_residual=True,
            source_term_width_ledger=terms,
            largest_standalone_term_width_function=dominant['exact_term_function'] if dominant else None,
            widths_are_conservative_enclosure_contributions_not_measured_true_errors=True,
            shared_dependencies_and_reader_rounding_prevent_additive_error_attribution=True)

    def row_accuracy(self, row, reader, target):
        out = relative_error_budget(self.ctx, row, target)
        c = self.ctx
        contributions = []
        for term in row['ratio_terms']:
            ratio = term['directed_ratio_enclosure']
            coefficient = c.mpf(term['signed_coefficient'])
            entry = dict(source_scale_log_function=term['source_log_scale'],
                exact_ratio_function=term['exact_ratio_function'],
                exact_log_ratio_function=term['exact_log_ratio_function'],
                original_signed_source_coefficient=term['signed_coefficient'],
                ratio_enclosure=ratio, ratio_method=term['method'])
            if ratio is None:
                entry['variation_bound_unresolved'] = True
            else:
                ratio = c.mpf(ratio)
                a, b = ends(coefficient)
                u, v = ends(ratio)
                dc, dr = c.mpf(b)-c.mpf(a), c.mpf(v)-c.mpf(u)
                # For any two points in the boxes, |a*r-b*s| is at most
                # max|r|*width(a)+max|a|*width(r), without a nominal value.
                coefficient_error = c.mpf(max(abs(u), abs(v)))*dc
                ratio_error = c.mpf(max(abs(a), abs(b)))*dr
                contribution = coefficient*ratio
                entry.update(normalized_contribution_enclosure=contribution,
                    coefficient_variation_budget=coefficient_error,
                    ratio_variation_budget=ratio_error,
                    normalized_variation_budget=coefficient_error+ratio_error,
                    strictly_positive_tail_upper_retained=term['method']['kind'] == 'retained_positive_tail',
                    budget_is_enclosure_variation_not_measured_true_error=True)
            contributions.append(entry)
        out['signed_coefficient_and_ratio_variation_ledger'] = contributions
        out['coefficient_ledger_keeps_shared_source_dependencies'] = True
        out['normalized_error_does_not_replace_absolute_physical_error'] = True
        if not row['exact_zero_enclosure']:
            out['actual_common_scale_source_ledger'] = self.ledger(row, reader)
        return out

    @source_precision
    def evaluate(self, delivery, relative_width_target='1/100000000'):
        self.assert_graph()
        target = correlated.inverse.exact_scalar(relative_width_target)
        if not 0 < target < 1:
            raise ValueError('Exact relative width target in (0,1) required')
        binding = self.before.source_binding(delivery)
        result = self.before.before.evaluate(delivery, relative_width_target)
        reader = self.before.amplitude.reader(delivery)
        counts = dict(rows=0, exact_zero_rows=0, nonzero_factored_target_rows=0,
            ordinary_width_target_rows=0, ordinary_numeric_target_rows=0,
            coefficient_sign_unresolved_rows=0, unresolved_ratio_rows=0)
        for components in result['physical_value_rows'].values():
            for component, rows in components.items():
                for row in rows.values():
                    accuracy = self.row_accuracy(row, reader, target)
                    row['physical_accuracy'] = accuracy
                    if component == 'p':
                        row['current_P0_numeric_inclusion_pending'] = False
                        row['ordinary_pressure_materialization_not_a_current_P0_binding'] = False
                    counts['rows'] += 1
                    counts['exact_zero_rows'] += int(row['exact_zero_enclosure'])
                    counts['nonzero_factored_target_rows'] += int(not row['exact_zero_enclosure'] and accuracy['factored_relative_width_satisfied'])
                    counts['ordinary_width_target_rows'] += int(accuracy['ordinary_numeric_relative_width_satisfied'])
                    counts['ordinary_numeric_target_rows'] += int(accuracy['ordinary_numeric_delivery_target_satisfied'])
                    counts['coefficient_sign_unresolved_rows'] += int(accuracy.get('sign') is None and row.get('common_scale_coefficient_enclosure') is not None)
                    counts['unresolved_ratio_rows'] += int(row.get('common_scale_coefficient_enclosure') is None and not row['exact_zero_enclosure'])
        result.pop('full_physical_accuracy_and_P0_source_binding_still_open', None)
        result.update(current_P0_numeric_source_binding=binding, physical_accuracy_counts=counts,
            error_budget_arithmetic_precision=400, ordinary_full_physical_accuracy_still_open=True,
            exact_factored_relative_accuracy_is_separate_from_ordinary_numeric_accuracy=True,
            **dict.fromkeys(GATES, self.acceptance_loaded))
        return result

    @source_precision
    def velocity_pressure(self, delivery, relative_width_target='1/100000000'):
        self.assert_graph()
        target = correlated.inverse.exact_scalar(relative_width_target)
        packet = self.before.velocity_pressure(delivery, relative_width_target)
        reader = self.before.amplitude.reader(delivery)
        for row in packet['values'].values():
            row['physical_accuracy'] = self.row_accuracy(row, reader, target)
        packet.update(exact_factored_relative_accuracy_is_separate_from_ordinary_numeric_accuracy=True,
            error_budget_arithmetic_precision=400, **dict.fromkeys(GATES, self.acceptance_loaded))
        return packet


@source_precision
def run(before, deliveries, relative_width_target='1/100'):
    began = time.monotonic()
    owner = CurrentOriginalRpPhysicalAccuracy(before, require_checked=False)
    if not deliveries:
        raise ValueError('Actual live original physical source deliveries required')
    views = {name: owner.evaluate(delivery, relative_width_target) for name, delivery in deliveries.items()}
    result = dict(candidate_current_physical_accuracy_constructed=True, source_family=owner.family_record,
        actual_physical_accuracy_views=views, relative_width_target=relative_width_target,
        source_graph_assertions=owner.assert_graph(), **dict.fromkeys(GATES+OPEN, False),
        input_hashes=owner.hashes, execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result), separators=(',', ':'))+'\n').encode(), mtime=0))
    print('Current physical accuracy: exact-scale relative bounds and actual source-scale width ledger', flush=True)
    return owner, views
