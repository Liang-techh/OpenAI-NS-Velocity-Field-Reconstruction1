"""Original Cartesian derivative values and source-correlated linear fields.

Cancel exact operator coefficients on the same source primitive before
enclosing curl/strain/Laplacian sums. Interval overlap is not a zero proof.
Nonlinear momentum, stress, flat remainder and recursion remain open.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_scale_value_arithmetic as scales
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

product, correlated, signed, box = scales.product, scales.correlated, scales.signed, scales.box
HERE, PREFIX, sha, ends = scales.HERE, scales.PREFIX, scales.sha, scales.ends
NAME = PREFIX + 'current_original_Rp_differential_field.json.gz'
RECEIPT = PREFIX + 'current_original_Rp_differential_field_check.json'
GATES = ('current_original_Rp_all_Cartesian4_time1_source_value_rows_installed',
         'current_original_Rp_source_primitive_linear_operator_cancellation_installed',
         'current_original_Rp_point_vorticity_strain_Laplacian_interface_installed')
OPEN = scales.OPEN
COMPONENTS = ('ux', 'uy', 'uz', 'p')
AXES = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def derivative_grid():
    indices = {(i, j, k) for i in range(5) for j in range(5-i) for k in range(5-i-j)}
    return {(component, index) for component in COMPONENTS for index in indices | {('t',)}}


@dataclass(frozen=True, eq=False)
class OriginalDifferentialField:
    """Live owner-bound physical derivative field; reports cannot hydrate it."""
    chart: str


def unit(component, derivative):
    base = 'pressure' if component == 'p' else 'velocity'
    if derivative == ('t',):
        return base + '/time'
    degree = sum(derivative)
    return base if not degree else base + '/length**' + str(degree)


def standard_forms():
    def term(component, derivative, weight=1):
        return (component, derivative, Fraction(weight))
    forms = {'divergence': (term('ux', AXES[0]), term('uy', AXES[1]), term('uz', AXES[2])),
        'vorticity_x': (term('uz', AXES[1]), term('uy', AXES[2], -1)),
        'vorticity_y': (term('ux', AXES[2]), term('uz', AXES[0], -1)),
        'vorticity_z': (term('uy', AXES[0]), term('ux', AXES[1], -1))}
    for i, component in enumerate(COMPONENTS[:3]):
        forms['Laplacian_' + component] = tuple(term(component, tuple(2 * v for v in axis)) for axis in AXES)
        for j in range(i, 3):
            forms['strain_' + str(i) + str(j)] = (term(component, AXES[j], Fraction(1, 2)),
                term(COMPONENTS[j], AXES[i], Fraction(1, 2)))
    return forms


def forms_definition():
    return {name: [dict(component=c, derivative=list(d), weight=signed.rational_record(q))
        for c, d, q in terms] for name, terms in standard_forms().items()}


class CurrentOriginalRpDifferentialField:
    @source_precision
    def __init__(self, before, require_checked=True):
        if type(before) is not scales.CurrentOriginalRpScaleValueArithmetic or not before.acceptance_loaded:
            raise ValueError('Accepted original source-scale arithmetic required')
        self.before, self.product = before, before.before
        self.graph, self.ctx, self.family_record = before.graph, before.ctx, before.family_record
        self.physical = self.product.amplitude.before.physical
        self.hashes = dict(before.hashes)
        for name in (scales.NAME, scales.RECEIPT, Path(__file__).name):
            self.hashes[name] = sha(name)
        self._fields = {}
        self.acceptance_loaded = False
        self._source_identities = copy.deepcopy(self.physical.canonical_source_identities)
        self._forms_definition = forms_definition()
        self._map_definition = correlated.physical.operator_definitions()
        self.proof = self.graph.node('current_original_Rp_differential_field_source_operations',
            actual_original_physical_map=self.physical.NAME if hasattr(self.physical, 'NAME') else self.physical.__class__.__name__,
            source_value_arithmetic_proof=before.proof.node,
            all_spatial_total_orders_through_four=True, fixed_physical_x_time_derivative=True,
            no_second_factorial_Jacobian_or_amplitude_derivative=True,
            exact_same_primitive_operator_cancellation_before_enclosure=True,
            independent_row_overlap_not_used_as_zero_identity=True,
            no_background_stress_flat_remainder_or_full_NS_claim=True)
        self._proof_snapshot = copy.deepcopy(self.graph.nodes[self.proof.node])
        self.assert_graph()
        if require_checked:
            receipt = json.loads((HERE / RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or \
                    any(receipt[k] for k in OPEN) or receipt['source_family'] != self.family_record or \
                    receipt['linear_operator_definitions'] != self._forms_definition or \
                    receipt['original_physical_map_definitions'] != self._map_definition:
                raise ValueError('Original differential field receipt/scope differs')
            for name in (Path(__file__).name, Path(__file__).stem + '_check.py', NAME):
                if receipt['input_hashes'].get(name) != sha(name):
                    raise ValueError('Unbound original differential source: ' + name)
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes, receipt['input_hashes'])
            self.acceptance_loaded = True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        checks = dict(accepted_original_value_arithmetic=self.before.acceptance_loaded,
            original_graph_context=self.graph is self.product.graph and self.ctx is self.product.ctx,
            original_source_family=self.family_record == self.product.family_record,
            original_operator_proof=self.graph.nodes[self.proof.node] == self._proof_snapshot,
            original_source_divergence_identity=self.physical.canonical_source_identities == self._source_identities,
            original_linear_operator_definitions=forms_definition() == self._forms_definition,
            original_physical_exponential_guard=signed.EXP_LOG_LIMIT == 1000)
        if not all(checks.values()):
            raise ValueError('Original differential field source differs: ' + str(checks))
        return checks

    @staticmethod
    def _signature(record):
        return (record['component'], record['unit'], record['scale'].node,
            tuple(sorted(record['scale_polynomial'].items())), record['coefficient_tuple'], record['sign'],
            id(record['root']), record['multiplier'], record['absolute_width_satisfied'],
            record['original_row']['component'], record['original_row']['physical_derivative'])

    def _entry_signature(self, entry):
        raw, row = entry['raw'], entry['row']
        if type(raw) is not correlated.physical.PhysicalSourceRow or (raw.component, raw.derivative) not in derivative_grid():
            raise ValueError('Typed original physical source row required')
        terms = []
        for term in raw.terms:
            if type(term) is not correlated.physical.PhysicalSourceTerm or \
                    type(term.source_row) is not correlated.physical.mixed.FactorizedMixedSourceRow or \
                    term.source_row.powers[-1] != 0:
                raise ValueError('Typed original physical source term required')
            polynomials = (term.source_row.coefficients, term.operator_coefficient, term.signed_coefficient)
            if any(type(value) is not correlated.physical.mixed.IntervalTaylor or value.order != 0 or
                    value.ctx is not self.physical.ctx for value in polynomials):
                raise ValueError('Ordinary original derivatives in their source context required')
            terms.append((term.source_label, term.source_row.name, term.source_row.derivative,
                term.source_row.source_units, term.source_row.powers,
                tuple((value.order, id(value.ctx), tuple(v._mpi_ for v in value.coefficients)) for value in polynomials),
                term.operator_function.node, term.radial_power, term.lambda_exponent.node,
                tuple((name, ref.node) for name, ref in term.log_scale_parts)))
        coefficient = row.get('common_scale_coefficient_enclosure')
        return (id(entry['handle']) if entry['handle'] is not None else None, entry['status'],
            raw.component, raw.derivative, tuple(terms), row['component'], row['physical_derivative'],
            row.get('exact_reference_log_scale_function'), None if coefficient is None else coefficient._mpi_,
            row['physical_accuracy'].get('sign'), row['physical_accuracy']['ordinary_numeric_relative_width_satisfied'],
            row['ordinary_numeric_materialized'])

    def _reader(self, record):
        return correlated.CorrelatedGraphBounds(correlated.locator.CancelledGraphBounds(self.graph, self.ctx,
            dict(record['bindings']), record['allowed']), dict(record['aliases']))

    def _require(self, field):
        self.assert_graph()
        item = self._fields.get(id(field))
        if type(field) is not OriginalDifferentialField or item is None or item[0] is not field or field.chart != item[1]['chart']:
            raise ValueError('Unchanged differential field issued by this owner required')
        record = item[1]
        self.product.amplitude.before._validate_delivery(record['delivery'])
        if {node: value._mpi_ for node, value in record['bindings'].items()} != record['binding_snapshot'] or \
                record['aliases'] != record['alias_snapshot'] or record['allowed'] != record['allowed_snapshot']:
            raise ValueError('Original differential field bindings/aliases/whitelist changed')
        if any(self.graph.nodes[node] != value for node, value in record['dag'].items()):
            raise ValueError('Original differential field defining DAG changed')
        if set(record['entries']) != derivative_grid() or set(record['entries']) != set(record['entry_snapshots']):
            raise ValueError('Original derivative row set changed')
        for key, entry in record['entries'].items():
            if self._entry_signature(entry) != record['entry_snapshots'][key][0]:
                raise ValueError('Original derivative entry or primitive changed')
            value = entry['handle']
            if value is None:
                continue
            issued = self.before._values.get(id(value))
            if issued is None or issued[0] is not value or self._signature(issued[1]) != record['entry_snapshots'][key][1]:
                raise ValueError('Original derivative value changed')
            if issued[1]['binding_snapshot'] != record['binding_snapshot'] or \
                    issued[1]['aliases'] != record['alias_snapshot'] or issued[1]['allowed'] != record['allowed_snapshot']:
                raise ValueError('Derivative source reader changed')
        return record

    @source_precision
    def evaluate(self, delivery, relative_width_target='1/1000'):
        self.assert_graph()
        target = Fraction(relative_width_target)
        if not 0 < target < 1:
            raise ValueError('Exact relative-width target in (0,1) required')
        original = self.product.evaluate(delivery, relative_width_target)
        mapped = self.product.before.physical_rows(delivery)
        reader = self.product.reader(delivery)
        entries, signatures, dag = {}, {}, {}
        bindings = dict(reader.bindings)
        snapshot = {node: value._mpi_ for node, value in bindings.items()}
        for section, components in original['physical_value_rows'].items():
            if section not in ('Cartesian_spatial_rows', 'fixed_x_time_rows'):
                raise ValueError('Original spatial/time section required')
            for component, rows in components.items():
                for label, row in rows.items():
                    derivative = tuple(row['physical_derivative'])
                    key = (component, derivative)
                    if key not in derivative_grid() or key in entries:
                        raise ValueError('Unique original component and ordinary derivative required')
                    raw = mapped[section][component][label] if section == 'Cartesian_spatial_rows' else mapped[section][component]
                    if type(raw) is not correlated.physical.PhysicalSourceRow or raw.component != component or raw.derivative != derivative:
                        raise ValueError('Same typed original physical derivative required')
                    coefficient = row.get('common_scale_coefficient_enclosure')
                    sign = row['physical_accuracy'].get('sign')
                    status = 'exact_zero' if row['exact_zero_enclosure'] else 'nonzero' if sign in (-1, 1) else 'sign_unresolved'
                    handle = None
                    if status == 'nonzero':
                        scale = box.pulse.radius.FunctionRef(self.graph, row['exact_reference_log_scale_function'])
                        record = dict(component=component + ':' + label, unit=unit(component, derivative),
                            scale=scale, scale_polynomial=reader.polynomial(scale), scale_dag=self.before._scale_dag(scale),
                            coefficient_tuple=coefficient._mpi_, sign=sign, root=object(), multiplier=Fraction(1),
                            delivery=delivery, bindings=bindings, binding_snapshot=snapshot,
                            aliases=dict(reader.aliases), alias_snapshot=dict(reader.aliases),
                            allowed=reader.allowed_exponentials, allowed_snapshot=reader.allowed_exponentials,
                            absolute_width_satisfied=row['physical_accuracy']['ordinary_numeric_relative_width_satisfied'],
                            original_materialized=row['ordinary_numeric_materialized'], original_row=row)
                        handle = self.before._issue(record)
                        issued_signature = self._signature(record)
                        dag.update(record['scale_dag'])
                    else:
                        issued_signature = None
                    entries[key] = dict(handle=handle, status=status, row=row, raw=raw)
                    signatures[key] = (self._entry_signature(entries[key]), issued_signature)
                    for term in raw.terms:
                        dag.update(self.before._scale_dag(term.operator_function))
        if set(entries) != derivative_grid():
            raise ValueError('Complete original spatial4/time1 grid required')
        field = OriginalDifferentialField(original['chart'])
        record = dict(chart=original['chart'], delivery=delivery, target=target, entries=entries,
            entry_snapshots=signatures, bindings=bindings, binding_snapshot=snapshot,
            aliases=dict(reader.aliases), alias_snapshot=dict(reader.aliases),
            allowed=reader.allowed_exponentials, allowed_snapshot=reader.allowed_exponentials,
            dag=dag, original=original, mapped=mapped, counts=tuple(sorted(original['delivery_counts'].items())))
        self._fields[id(field)] = (field, record)
        return field

    @source_precision
    def derivative(self, field, component, derivative=(0, 0, 0)):
        record = self._require(field)
        key = (component, tuple(derivative))
        if key not in record['entries']:
            raise ValueError('Original component and available derivative order required')
        entry = record['entries'][key]
        if entry['handle'] is None:
            raise ValueError('Derivative has ' + entry['status'] + '; inspect its reported enclosure')
        return entry['handle']

    def _linear(self, record, terms):
        reader = self._reader(record)
        basis = {}
        for component, derivative, weight in terms:
            key = (component, tuple(derivative))
            if key not in record['entries']:
                raise ValueError('Available original derivative required')
            basis[key] = basis.get(key, Fraction(0)) + Fraction(weight)
        basis = {key: value for key, value in basis.items() if value}
        units = {unit(*key) for key in basis}
        if len(units) > 1:
            raise ValueError('Linear combination requires matching physical units')
        primitives = {}
        for key, weight in basis.items():
            for term in record['entries'][key]['raw'].terms:
                logref = self.graph.add(*(ref for _, ref in term.log_scale_parts))
                poly = reader.polynomial(logref)
                primitive = (term.source_label, term.source_row.derivative, term.source_row.source_units,
                    term.source_row.powers, tuple(sorted(poly.items())))
                value = self.ctx.mpf(term.source_row.coefficients[0])
                if primitive not in primitives:
                    primitives[primitive] = dict(source=value, source_tuple=value._mpi_, log=logref,
                        poly=poly, exact=[], bounds=[], source_label=term.source_label,
                        source_derivative=term.source_row.derivative,
                        source_units=term.source_row.source_units, powers=term.source_row.powers)
                group = primitives[primitive]
                if value._mpi_ != group['source_tuple']:
                    raise ValueError('Identical primitive must keep its original coefficient enclosure')
                group['exact'].append(self.graph.mul(self.graph.constant(weight), term.operator_function))
                group['bounds'].append((weight, term.operator_coefficient[0]))
        groups, cancelled, ledger = {}, [], []
        for primitive, part in primitives.items():
            operator = self.graph.add(*part['exact'])
            exact_zero = not reader.polynomial(operator)
            item = dict(source_label=part['source_label'], ordinary_source_derivative=part['source_derivative'],
                exact_combined_operator=operator.node, source_log_scale=part['log'].node,
                exact_operator_zero=exact_zero, original_source_coefficient=part['source'])
            if exact_zero:
                cancelled.append(item)
                continue
            bound = sum((self.ctx.mpf(q.numerator) / q.denominator * self.ctx.mpf(value)
                for q, value in part['bounds']), self.ctx.mpf(0))
            coefficient = part['source'] * bound
            groupkey = (tuple(sorted(part['poly'].items())), part['source_units'], part['powers'])
            if groupkey not in groups:
                canonical = signed.polynomial_function(self.graph, part['poly'])
                groups[groupkey] = dict(log_function=canonical, log_bound=reader.at(canonical),
                    coefficient=self.ctx.mpf(0), proofs=[])
            groups[groupkey]['coefficient'] += coefficient
            item.update(directed_combined_operator=bound, signed_primitive_coefficient=coefficient)
            ledger.append(item)
        result = signed.enclose_factored_sum(self.graph, self.ctx, reader, list(groups.values()), record['target'])
        result.update(original_derivative_basis=[dict(component=key[0], physical_derivative=key[1],
                weight=signed.rational_record(weight)) for key, weight in sorted(basis.items())],
            physical_unit=next(iter(units), None), original_primitive_ledger=ledger,
            exact_cancelled_primitive_operators=cancelled,
            exact_primitive_cancellation_count=len(cancelled),
            interval_zero_overlap_not_claimed_as_exact_cancellation=True,
            all_absolute_and_NS_accuracy_gates_unpromoted=True)
        return result

    @source_precision
    def linear_combination(self, field, terms):
        return self._linear(self._require(field), terms)

    @source_precision
    def report(self, field):
        record = self._require(field)
        sections = {'Cartesian_spatial_rows': {}, 'fixed_x_time_rows': {}}
        for (component, derivative), entry in record['entries'].items():
            row = entry['row']
            output = dict(component=component, physical_derivative=derivative,
                physical_unit=unit(component, derivative), status=entry['status'],
                exact_reference_log_scale_function=row.get('exact_reference_log_scale_function'),
                signed_coefficient_enclosure=row.get('common_scale_coefficient_enclosure'),
                sign=row['physical_accuracy'].get('sign'),
                original_absolute_width_target_satisfied=row['physical_accuracy']['ordinary_numeric_relative_width_satisfied'],
                ordinary_numeric_materialized=row['ordinary_numeric_materialized'])
            if derivative == ('t',):
                sections['fixed_x_time_rows'][component] = output
            else:
                label = 'x%d_y%d_z%d' % derivative
                sections['Cartesian_spatial_rows'].setdefault(component, {})[label] = output
        linear = {name: self._linear(record, terms) for name, terms in standard_forms().items()}
        return dict(chart=record['chart'], source_family=self.family_record, derivative_rows=sections,
            derivative_counts=dict(record['counts']), actual_linear_differential_fields=linear,
            original_analytic_divergence_identity=copy.deepcopy(self._source_identities),
            numeric_divergence_sum_is_separate_from_analytic_identity=True,
            background_momentum_products_stress_and_flat_remainder_not_installed=True,
            unrestricted_physical_point_API=False, full_certified_physical_accuracy=False,
            differential_field_proof=self.proof.node,
            **dict.fromkeys(GATES, self.acceptance_loaded), **dict.fromkeys(OPEN, False))


@source_precision
def run(before, deliveries, relative_width_target='1/1000'):
    began = time.monotonic()
    owner = CurrentOriginalRpDifferentialField(before, require_checked=False)
    fields = {name: owner.evaluate(delivery, relative_width_target) for name, delivery in deliveries.items()}
    views = {name: owner.report(field) for name, field in fields.items()}
    result = dict(source_family=owner.family_record, actual_differential_fields=views,
        source_assertions=owner.assert_graph(), input_hashes=owner.hashes,
        execution_seconds=time.monotonic() - began, **dict.fromkeys(GATES + OPEN, False))
    (HERE / NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result), separators=(',', ':')) + '\n').encode(), mtime=0))
    return owner, fields, views
