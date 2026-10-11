"""Nonzero source-scaled physical values with controlled ratio operations.

An inverse-mu source expression is retained, not declared zero-width.
Ratios cancel identical source polynomials before directed arithmetic.
Absolute physical numeric materialization remains an independent open gate.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_source_product_arithmetic as product
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

correlated, signed, box = product.correlated, product.signed, product.box
HERE, PREFIX, sha, ends = product.HERE, product.PREFIX, product.sha, product.ends
NAME = PREFIX + 'current_original_Rp_scale_value_arithmetic.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_scale_value_arithmetic_check.json'
GATES = ('current_original_Rp_nonzero_source_scaled_value_operations_installed',
         'current_original_Rp_exact_inverse_mu_scale_partition_installed',
         'current_original_Rp_shared_scale_ratio_cancellation_installed')
OPEN = product.OPEN


@dataclass(frozen=True, eq=False)
class SourceScaledPhysicalValue:
    """Opaque handle: only its issuing owner supplies arithmetic and export."""
    component: str


class CurrentOriginalRpScaleValueArithmetic:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not product.CurrentOriginalRpSourceProductArithmetic or not before.acceptance_loaded:
            raise ValueError('Accepted original source-product arithmetic required')
        self.before, self.graph, self.ctx = before, before.graph, before.ctx
        self.family_record = before.family_record
        self.hashes = dict(before.hashes)
        for name in (product.NAME, product.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.acceptance_loaded = False
        self._values = {}
        self.inverse_mu = self.graph.quotient(self.graph.one, before.amplitude.mu,
            'current original positive denominator; exact reciprocal normalization')
        self._inverse_snapshot = copy.deepcopy(self.graph.nodes[self.inverse_mu.node])
        self.proof = self.graph.node('current_original_Rp_source_scaled_physical_value_arithmetic',
            exact_inverse_mu_source=self.inverse_mu.node,
            original_mu_source=before.amplitude.mu.node,
            absolute_inverse_mu_uncertainty_not_removed=True,
            scale_differences_cancel_by_exact_source_polynomials=True,
            coefficient_bounds_and_positive_tails_retained=True,
            issued_aliases_and_exponential_whitelist_sealed=True,
            existing_absolute_physical_width_flags_not_promoted=True,
            numeric_materialization_is_separate_from_source_scaled_delivery=True,
            original_physical_exponential_guard=signed.EXP_LOG_LIMIT)
        self._proof_snapshot = copy.deepcopy(self.graph.nodes[self.proof.node])
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or \
                    any(receipt[key] for key in OPEN) or receipt['source_family'] != self.family_record:
                raise ValueError('Source-scaled arithmetic receipt/scope differs')
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks = dict(accepted_original_source_product_owner=self.before.acceptance_loaded,
            same_source_context=self.graph is self.before.graph and self.ctx is self.before.ctx,
            same_source_family=self.family_record == self.before.family_record,
            exact_original_inverse_mu_unchanged=self.graph.nodes[self.inverse_mu.node] == self._inverse_snapshot,
            arithmetic_proof_unchanged=self.graph.nodes[self.proof.node] == self._proof_snapshot,
            physical_exponential_guard_unchanged=signed.EXP_LOG_LIMIT == 1000)
        if not all(checks.values()):
            raise ValueError('Source-scaled arithmetic source differs: ' + str(checks))
        return checks

    def _fresh_reader(self, record):
        # No caller memo/polynomial cache can act as an arithmetic oracle.
        reader = correlated.CorrelatedGraphBounds(correlated.locator.CancelledGraphBounds(
            self.graph, self.ctx, dict(record['bindings']), record['allowed']), dict(record['aliases']))
        if reader.polynomial(record['scale']) != record['scale_polynomial']:
            raise ValueError('Original issued physical log scale changed')
        return reader

    def _scale_dag(self, scale):
        pending, result = [scale.node], {}
        while pending:
            node = pending.pop()
            if node in result:
                continue
            data = self.graph.nodes[node]
            result[node] = copy.deepcopy(data)
            op = data['operation']
            if op in ('sum', 'product'):
                pending.extend(data['arguments'])
            elif op in ('negative', 'analytic_unary'):
                pending.append(data['argument'])
            elif op == 'positive_quotient':
                pending.extend((data['numerator'], data['denominator']))
        return result

    def _require(self, value):
        self.assert_graph()
        item = self._values.get(id(value))
        if type(value) is not SourceScaledPhysicalValue or item is None or item[0] is not value or \
                value.component != item[1]['component']:
            raise ValueError('Unchanged source-scaled value issued by this owner required')
        record = item[1]
        self.before.amplitude.before._validate_delivery(record['delivery'])
        if {node: bound._mpi_ for node, bound in record['bindings'].items()} != record['binding_snapshot']:
            raise ValueError('Issued source-scaled bindings changed')
        if record['aliases'] != record['alias_snapshot'] or frozenset(record['allowed']) != record['allowed_snapshot']:
            raise ValueError('Issued source aliases or exponential whitelist changed')
        if any(self.graph.nodes[node] != data for node, data in record['scale_dag'].items()):
            raise ValueError('Issued physical scale defining DAG changed')
        self._fresh_reader(record)
        return record

    def _issue(self, record):
        value = SourceScaledPhysicalValue(record['component'])
        self._values[id(value)] = (value, record)
        return value

    @source_precision
    def velocity_pressure(self, delivery, relative_width_target='1/1000'):
        self.assert_graph()
        target = Fraction(relative_width_target)
        if not 0 < target < 1:
            raise ValueError('Exact relative-width target in (0,1) required')
        packet = self.before.velocity_pressure(delivery, relative_width_target)
        reader = self.before.reader(delivery)
        values = {}
        for component, row in packet['values'].items():
            coefficient = row['common_scale_coefficient_enclosure']
            sign = row['physical_accuracy']['sign']
            if coefficient is None or sign not in (-1, 1):
                raise ValueError('A proved nonzero original component is required')
            scale = box.pulse.radius.FunctionRef(self.graph, row['exact_reference_log_scale_function'])
            poly = reader.polynomial(scale)
            record = dict(component=component, unit='pressure' if component == 'p' else 'velocity',
                scale=scale, scale_polynomial=poly, coefficient_tuple=coefficient._mpi_, sign=sign,
                scale_dag=self._scale_dag(scale),
                root=object(), multiplier=Fraction(1), delivery=delivery,
                bindings=dict(reader.bindings), aliases=dict(reader.aliases), allowed=reader.allowed_exponentials,
                alias_snapshot=dict(reader.aliases), allowed_snapshot=frozenset(reader.allowed_exponentials),
                binding_snapshot={node: bound._mpi_ for node, bound in reader.bindings.items()},
                absolute_width_satisfied=row['physical_accuracy']['ordinary_numeric_relative_width_satisfied'],
                original_materialized=row['ordinary_numeric_materialized'],
                original_row=row)
            values[component] = self._issue(record)
        return dict(values=values, source_family=self.family_record,
            source_scaled_nonzero_delivery=True, ordinary_numeric_materialized=False,
            full_certified_physical_accuracy=False, unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES, self.acceptance_loaded), **dict.fromkeys(OPEN, False))

    def _coefficient(self, record):
        # Restore exact endpoints, with no scalar-precision midpoint conversion.
        value = self.ctx.mpf(0)
        value._mpi_ = record['coefficient_tuple']
        return value

    @source_precision
    def export(self, value):
        record = self._require(value)
        reader = self._fresh_reader(record)
        center = self.before.center(self.before.reader(record['delivery']), record['scale'])
        poly = dict(record['scale_polynomial'])
        fixed_coefficient = poly.pop((self.before.pivot.node,), Fraction(0))
        inverse_coefficient = poly.pop((self.inverse_mu.node,), Fraction(0))
        fixed = self.graph.mul(self.graph.constant(fixed_coefficient), self.before.pivot)
        inverse = self.graph.mul(self.graph.constant(inverse_coefficient), self.inverse_mu)
        residual = signed.polynomial_function(self.graph, poly)
        if reader.polynomial(self.graph.sub(record['scale'], self.graph.add(fixed, inverse, residual))):
            raise ValueError('Exact original fixed/inverse/residual partition required')
        proof = self.graph.node('exact_original_Rp_source_scaled_inverse_mu_partition',
            original_scale_function=record['scale'].node, fixed_origin_function=fixed.node,
            exact_inverse_mu_term=inverse.node, exact_residual_function=residual.node,
            proof_method='exact rational polynomial under original source aliases',
            inverse_mu_defining_source_retained=True, inverse_mu_uncertainty_not_discarded=True)
        return dict(component=record['component'], physical_unit=record['unit'],
            exact_original_log_scale_function=record['scale'].node,
            exact_positive_scale_function=self.graph.unary('exp', record['scale']).node,
            signed_common_scale_coefficient_enclosure=self._coefficient(record), sign=record['sign'],
            proved_nonzero=True, exact_fixed_origin=center['exact_fixed_origin'],
            exact_source_inverse_mu=dict(source_function=self.inverse_mu.node,
                original_mu_function=self.before.amplitude.mu.node,
                coefficient=signed.rational_record(inverse_coefficient), term_function=inverse.node,
                directed_term_enclosure=reader.at(inverse),
                defining_identity='inverse_mu=1/mu, with mu=exp(-4*(exp(40)+11))/1000',
                zero_error_variation=False, directed_uncertainty_retained=True),
            exact_remaining_log_scale_function=residual.node,
            directed_remaining_log_scale=reader.at(residual),
            exact_partition_identity=proof.node,
            original_total_scale_width=center['residual_log_width_bound'],
            original_absolute_width_target_satisfied=record['absolute_width_satisfied'],
            ordinary_numeric_materialized=False, absolute_scale_accuracy_not_promoted=True,
            operation_precision=self.ctx.dps, source_family=self.family_record,
            arithmetic_proof=self.proof.node,
            root_original_component=record['original_row']['component'],
            exact_rational_multiplier=signed.rational_record(record['multiplier']))

    @source_precision
    def multiply_rational(self, value, factor):
        record = self._require(value)
        q = Fraction(factor)
        if not q:
            raise ValueError('This nonzero value API requires a nonzero rational factor')
        fresh = dict(record)
        fresh.update(coefficient_tuple=(self._coefficient(record) * (self.ctx.mpf(q.numerator) / q.denominator))._mpi_,
            sign=record['sign'] * (1 if q > 0 else -1), multiplier=record['multiplier'] * q)
        return self._issue(fresh)

    @source_precision
    def ratio(self, numerator, denominator):
        a, b = self._require(numerator), self._require(denominator)
        if a['delivery'] is not b['delivery'] or a['binding_snapshot'] != b['binding_snapshot'] or \
                a['aliases'] != b['aliases'] or a['allowed'] != b['allowed']:
            raise ValueError('Ratio requires the same source-issued physical delivery and bindings')
        reader = self._fresh_reader(a)
        difference = self.graph.sub(a['scale'], b['scale'])
        poly = reader.polynomial(difference)
        if poly != {key: value for key, value in self._subtract(a['scale_polynomial'], b['scale_polynomial']).items() if value}:
            raise ValueError('Exact issued source-scale subtraction required')
        bound = reader.at(difference)
        if not poly and a['root'] is b['root']:
            q = a['multiplier'] / b['multiplier']
            coefficient = self.ctx.mpf(q.numerator) / q.denominator
            same_value = True
        else:
            coefficient = self._coefficient(a) / self._coefficient(b)
            same_value = False
        ratio, method = (self.ctx.mpf(1), dict(kind='exact_common_scale', exact_ratio_is_one=True)) \
            if not poly else signed.ratio_enclosure(self.ctx, bound)
        numeric = None if ratio is None else coefficient * ratio
        sign = a['sign'] * b['sign']
        proof = self.graph.node('exact_original_Rp_source_scaled_ratio_cancellation',
            numerator_scale=a['scale'].node, denominator_scale=b['scale'].node,
            exact_log_scale_difference=difference.node,
            exact_polynomial_zero=not poly, same_original_value_rational_multiple=same_value,
            identical_source_delivery_required=True,
            defining_identity='(exp(La)*ca)/(exp(Lb)*cb)=exp(La-Lb)*(ca/cb)',
            coefficient_enclosure_is_not_an_exact_coefficient_substitution=True)
        return dict(sign=sign, exact_ratio_is_nonzero=True,
            numerator_unit=a['unit'], denominator_unit=b['unit'],
            exact_scale_log_difference=difference.node, directed_scale_log_difference=bound,
            exact_positive_scale_ratio_function=self.graph.unary('exp', difference).node,
            scale_polynomial_terms=len(poly), common_scale_cancelled_exactly=not poly,
            signed_coefficient_ratio_enclosure=coefficient, directed_scale_ratio_enclosure=ratio,
            ordinary_numeric_ratio_enclosure=numeric, ratio_method=method,
            same_original_value_rational_multiple=same_value,
            numerator_denominator_absolute_materialization_required=False,
            absolute_physical_accuracy_not_promoted=True, exact_ratio_identity=proof.node)

    @staticmethod
    def _subtract(a, b):
        result = dict(a)
        for key, value in b.items():
            result[key] = result.get(key, Fraction(0)) - value
        return result

    @source_precision
    def compare(self, left, right):
        a, b = self._require(left), self._require(right)
        if a['unit'] != b['unit']:
            raise ValueError('Comparison requires matching physical units')
        ratio = self.ratio(left, right)
        if a['sign'] != b['sign']:
            order = 1 if a['sign'] > b['sign'] else -1
            log_ratio = None
        else:
            magnitude = ratio['signed_coefficient_ratio_enclosure']
            log_ratio = ratio['directed_scale_log_difference'] + self.ctx.ln(magnitude)
            lo, hi = ends(log_ratio)
            order = (a['sign'] if lo > 0 else -a['sign'] if hi < 0 else
                     0 if lo == hi == 0 and ratio['same_original_value_rational_multiple'] else None)
        return dict(order=order, directed_log_magnitude_ratio=log_ratio,
            meaning='-1: left<right; 0: exact equality; 1: left>right; null: unresolved',
            exact_ratio_identity=ratio['exact_ratio_identity'],
            numeric_absolute_values_not_materialized=True)


@source_precision
def run(before, deliveries, relative_width_target='1/1000'):
    began = time.monotonic()
    owner = CurrentOriginalRpScaleValueArithmetic(before, require_checked=False)
    views = {}
    for name, delivery in deliveries.items():
        packet = owner.velocity_pressure(delivery, relative_width_target)
        u, v, w, p = (packet['values'][key] for key in ('u', 'v', 'w', 'p'))
        views[name] = dict(values={key: owner.export(value) for key, value in packet['values'].items()},
            u_over_v=owner.ratio(u, v), u_over_u=owner.ratio(u, u),
            w_over_u=owner.ratio(w, u), u_over_p=owner.ratio(u, p),
            u_compared_with_w=owner.compare(u, w),
            twice_u_over_u=owner.ratio(owner.multiply_rational(u, 2), u))
    result = dict(source_family=owner.family_record, actual_scale_value_operations=views,
        input_hashes=owner.hashes, source_assertions=owner.assert_graph(),
        execution_seconds=time.monotonic() - began, **dict.fromkeys(GATES + OPEN, False))
    (HERE / NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result), separators=(',', ':')) + '\n').encode(), mtime=0))
    return owner, views
