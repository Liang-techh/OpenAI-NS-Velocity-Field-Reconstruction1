"""All-original-chart physical fields used by the n1 pressure continuation.

This is a small adapter over the accepted all-chart partial-mean owner.  It
binds the original leading V and own-m histories, applies the one common
physical unit Pstar, and retains signed value/Z1 pressure integrands.  Source
cell unions are exact native unions, including the symbolic sc_half and Rh
endpoints.  Magnitude rows are directed whole-cell enclosures; they are never
used as signed field values.
"""
import copy
import json
import math
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as all_chart
import lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields as physical


partial = all_chart.partial
c3, current, source = all_chart.c3, all_chart.current, all_chart.source
rg = all_chart.rg
HERE, PREFIX, sha, ep = all_chart.HERE, all_chart.PREFIX, all_chart.sha, all_chart.ep
source_precision = all_chart.source_precision
fingerprint = all_chart.fingerprint
OPEN = tuple(dict.fromkeys(all_chart.OPEN))

NAME = PREFIX + 'current_original_n1_all_chart_physical_pressure_fields.json'
RECEIPT = PREFIX + 'current_original_n1_all_chart_physical_pressure_fields_check.json'
GATE = 'current_original_n1_all_chart_physical_pressure_fields_installed'
I0_NAMES = ('q_Z', 'U_q_Z', 'U_squared', 'U_U_Z')
I1_NAMES = ('q', 'q_Z', 'q_squared', 'U_q_Z', 'U_squared', 'U_U_Z')


@dataclass(frozen=True, eq=False)
class OriginalAllChartPhysicalPressurePacket:
    leading_hierarchy_order: int = 0


def _as_ref(g, value):
    return value if isinstance(value, source.FunctionRef) else source.FunctionRef(g, value)


def _canonical_endpoint(value):
    if isinstance(value, list):
        return tuple(_canonical_endpoint(item) for item in value)
    return value


def _endpoint_key(value):
    return json.dumps(c3.encoded(_canonical_endpoint(value)),
                      sort_keys=True, separators=(',', ':'))


def _interval_key(left, right):
    return json.dumps([c3.encoded(_canonical_endpoint(left)),
                       c3.encoded(_canonical_endpoint(right))],
                      sort_keys=True, separators=(',', ':'))


def _same_endpoint(left, right):
    return _canonical_endpoint(left) == _canonical_endpoint(right)


def _same_source(value, c, flow, label):
    if value.ctx is not c or value.scale.bases is not flow.logs or value.ledger is not flow.ledger:
        raise ValueError('Same live source context/basis/ledger required for ' + label)


def _equivalent_rows(actual, expected):
    return partial.current.current.previous.equivalent_rows(actual, expected)


def source_rows(owner, ends, chart, left, right=None):
    """Return the actual original V/m/R/Pstar source objects for one union."""
    if type(owner) is not CurrentOriginalN1AllChartPhysicalPressureFields:
        raise ValueError('All-chart physical source provider requires retained owner')
    ends = tuple(ends)
    if chart not in owner.charts:
        raise ValueError('Accepted original chart required')
    packet = owner.modified.target.phase.leading_source(ends, chart, left, right)
    op = owner.modified.ranges.phase.outer.owner.owner(ends)
    if (packet.get('source_identity') != owner.family or
            packet.get('exact_common_P0_axial5') is not op.P0 or
            packet.get('actual_phase_Z_exact_zero') is not True):
        raise ValueError('Same physical source family/P0 and Z-independent phase required')
    generic = packet.get('original_generic_source')
    if not isinstance(generic, dict):
        raise ValueError('Original generic source packet missing')
    histories = generic.get('common_own_five_histories_axial5')
    required = ('common_velocity_V_axial5', 'original_physical_radius')
    if any(name not in generic for name in required) or not isinstance(histories, dict):
        raise ValueError('Required original V/m/radius source quantities absent')
    V, m, R = generic['common_velocity_V_axial5'], histories.get('m'), generic['original_physical_radius']
    if any(len(rows) < 4 for rows in (V, m)):
        raise ValueError('Original V/m Taylor rows through order 3 required')
    flow = op.flow
    if (len(op.P0) < 4 or
            not _equivalent_rows(generic.get('common_original_P0_axial5', ()), op.P0)):
        raise ValueError('Equivalent original independent P0 Taylor rows required')
    for k, value in enumerate(op.P0):
        _same_source(value, owner.c, flow, 'P0 Taylor row %d' % k)
    for label, values in (('V', V), ('m', m)):
        for k, value in enumerate(values):
            _same_source(value, owner.c, flow, '%s Taylor row %d' % (label, k))
    _same_source(R, owner.c, flow, 'original physical radius')
    _same_source(op.Pstar, owner.c, flow, 'Pstar')
    return dict(V=V, m=m, R=R, Pstar=op.Pstar, source_packet=packet,
                exact_common_P0=op.P0, source_family=owner.family)


def physical_source_scale(owner, ends):
    if type(owner) is not CurrentOriginalN1AllChartPhysicalPressureFields:
        raise ValueError('All-chart physical scale provider requires retained owner')
    op = owner.modified.ranges.phase.outer.owner.owner(tuple(ends))
    _same_source(op.Pstar, owner.c, op.flow, 'Pstar')
    return op.Pstar


def _query_substitution(g, body, chart, query):
    return g.node('substitute_free_query_at_native_coordinate', body=body.node,
                   free_query=query.node,
                   native_variable='pressure_moment_native_' + chart,
                   bound_partial_mean_native_variable='native_' + chart,
                   free_query_is_outer_pressure_endpoint=True,
                   defining_source=current.ast_binding(_query_substitution),
                   source_query_substitution='free physical query -> pressure moment native coordinate',
                   inner_Duhamel_native_variable_preserved=True,
                   same_original_source_and_global_phase=True)


def _physical_radius(g, log_radius, chart, query):
    return g.node('current_original_physical_radius', log_radius=log_radius.node,
                  free_query=query.node, native_chart=chart,
                  raw_quantity_path='original_generic_source.original_physical_radius',
                  radius_expression='exp(original physical log-radius)',
                  original_log_radius_is_absolute=True,
                  radius_Z_independent=True,
                  defining_source=current.ast_binding(_physical_radius),
                  range_enclosure_is_not_selected_function_value=True)


def _at_query(g, value, chart, query):
    return c3.C3Function(*(partial.at_native(g, _as_ref(g, row), chart, query)
                           for row in physical.rows(value)))


def _m_source_function(g, owner, chart):
    template = owner.modified.target.primitives[chart]['source_roots']['V']
    attrs = dict(g.nodes[template.value.node])
    coordinate, Z_variable = attrs.get('coordinate'), attrs.get('Z_variable')
    if Z_variable is None:
        raise ValueError('Original leading source Z variable metadata missing')
    out = []
    for k in range(4):
        out.append(g.node(
            'current_original_physical_source_quantity',
            source_family=owner.family, native_chart=chart, quantity='m', Z_order=k,
            coordinate=coordinate, Z_variable=Z_variable,
            phase_independent_leading_source=True,
            defining_quantity_not_a_range_value=True,
            ordinary_slow_Z_derivative_order=k,
            Taylor_coefficient_factorial=math.factorial(k),
            source_projection=('factorial(%d)*original_generic_source.'
                               'common_own_five_histories_axial5["m"][%d]') % (k, k),
            raw_quantity_path=('original_generic_source.'
                               'common_own_five_histories_axial5["m"][%d]') % k,
            source_rows_provider=current.ast_binding(source_rows),
            source_projection_binding=current.ast_binding(source_rows),
            source_row_is_Taylor_coefficient=True,
            factorial_conversion_applied_once=True,
            original_leading_source_recipe=attrs.get('recipe'),
            same_live_basis_context_and_ledger=True,
            range_enclosure_is_not_selected_function_value=True))
    return c3.C3Function(*out)


def _physical_scale_function(g, owner):
    node = g.node('current_original_physical_source_scale',
                  source_family=owner.family, name='Pstar',
                  raw_quantity_path='outer.owner.owner(ends).Pstar',
                  source_rows_provider=current.ast_binding(physical_source_scale),
                  source_projection_binding=current.ast_binding(physical_source_scale),
                  Z_independent=True, physical_common_unit=True,
                  defining_quantity_not_a_range_value=True,
                  same_live_basis_context_and_ledger=True)
    return c3.C3Function(node, g.zero, g.zero, g.zero)


def _original_log_r_minus(g, owner):
    params = owner.modified.target.phase.built['parameters']
    return g.add(g.unary('log', g.constant(4)),
                 g.neg(g.mul(g.constant(4), params['logP'])),
                 g.constant(-1000),
                 g.mul(params['hbB'], params['sc'], g.constant('1/2')))


def _ordered_cells(owner, chart):
    base = owner.modified.target.phase.windows[chart]
    cells = list(base['actual_current_source_cells'])
    partition = tuple(owner.partitions[chart])
    if len(cells) != len(partition) - 1:
        raise ValueError('Original source partition cell count changed')
    for cell, left, right in zip(cells, partition, partition[1:]):
        if (not _same_endpoint(cell['exact_left'], left) or
                not _same_endpoint(cell['exact_right'], right)):
            raise ValueError('Original source cells no longer match native partition')
    return cells


def _defining_integrals(g, owner, chart, integrands, radius, jacobian, query):
    # Geometry Jacobians are source functions in the original native
    # coordinate for charts such as actual_patch (J=1/native_actual_patch).
    # Bind that coordinate to the free query before the outer pressure
    # integral substitutes free query -> pressure_moment_native_<chart>.
    jacobian = partial.at_native(g, _as_ref(g, jacobian), chart, query)
    cells = _ordered_cells(owner, chart)
    by_cell = {name: {} for name in integrands}
    for name, value in integrands.items():
        power = 1 if name.startswith('I0_') else 2
        radius_power = radius if power == 1 else g.mul(radius, radius)
        for cell in cells:
            lower, upper = _as_ref(g, cell['lower']), _as_ref(g, cell['upper'])
            rows_integral = []
            for order, row in enumerate(physical.z1_rows(value)):
                body = g.mul(radius_power, g.mul(_as_ref(g, row), jacobian))
                body = _query_substitution(g, body, chart, query)
                rows_integral.append(current.controls.integral(
                    g, body, 'pressure_moment_native_' + chart, lower, upper,
                    measure=('physical dR=R*Jacobian*dnative exactly once'
                             if power == 1 else
                             'physical R*dR=R^2*Jacobian*dnative exactly once'),
                    ordinary_slow_Z_derivative_order=order,
                    free_query_endpoint=query.node,
                    free_query_substituted_at_native=True,
                    original_physical_radius_source='original_generic_source.original_physical_radius',
                    original_physical_Jacobian_applied_once=True,
                    endpoint_independent_of_Z=True))
            by_cell[name][_interval_key(cell['exact_left'], cell['exact_right'])] = physical.Z1Function(*rows_integral)
    unions = {name: {} for name in integrands}
    for first in range(len(cells)):
        for last in range(first, len(cells)):
            key = _interval_key(cells[first]['exact_left'], cells[last]['exact_right'])
            for name in integrands:
                value = by_cell[name][_interval_key(cells[first]['exact_left'], cells[first]['exact_right'])]
                for index in range(first + 1, last + 1):
                    value = physical._z1_add(g, value,
                        by_cell[name][_interval_key(cells[index]['exact_left'], cells[index]['exact_right'])])
                unions[name][key] = value
    return unions


def _field_functions(owner, chart):
    g = owner.g
    alg = c3.C3Algebra(g)
    query = owner.parent.functions[chart]['query']
    base = owner.modified.target.phase.windows[chart]
    leading = owner.modified.target.primitives[chart]['source_roots']
    scale = _physical_scale_function(g, owner)
    baseline_V = _at_query(g, leading['V'], chart, query)
    baseline_m = _at_query(g, _m_source_function(g, owner, chart), chart, query)
    delta_V = owner.parent.functions[chart]['delta_V']
    delta_m = owner.parent.functions[chart]['mean']
    U = alg.mul(scale, alg.add(baseline_V, delta_V))
    m = alg.mul(scale, alg.add(baseline_m, delta_m))
    q = partial.recover_Q_increment(g, U, m, owner.delta_node)
    integrands = physical.nonradial_integrands(g, U, q)
    offset = _as_ref(g, base['offset'])
    relative = partial.at_native(g, offset, chart, query)
    radius = _physical_radius(g, g.add(_original_log_r_minus(g, owner), relative), chart, query)
    unions = _defining_integrals(g, owner, chart, integrands, radius,
                                 _as_ref(g, base['Jacobian']), query)
    return dict(query=query, scale=scale, baseline_V=baseline_V,
                baseline_m=baseline_m, delta_V=delta_V, delta_m=delta_m,
                U=U, m=m, m_y=alg.add(U, alg.neg(m)), q=q,
                q_Z=physical._z1_shift(q), U_Z=physical._z1_shift(U),
                radius=radius, integrands=integrands,
                defining_integrals_by_cell=unions,
                source_geometry=dict(
                    raw_quantity_path='original_generic_source.original_physical_radius',
                    original_log_r_minus='log(4)-4*logP-1000+hbB*sc/2',
                    relative_log_radius_offset=base['offset'],
                    Jacobian=base['Jacobian'],
                    pressure_moment_native_variable='pressure_moment_native_' + chart,
                    bound_partial_mean_native_variable='native_' + chart,
                    physical_measure='R^power*Jacobian*dnative with power 1 or 2'))


def _raw_ordinary_jet(c, bd, values, label):
    if len(values) < 4:
        raise ValueError('Four source Taylor rows required for ' + label)
    return rg.JetBound(*(bd.row(values[k] * math.factorial(k)) for k in range(4)))


def _bound_record(value):
    if isinstance(value, physical.Z1Bound):
        return [value.value.record(), value.Z.record()]
    if isinstance(value, rg.JetBound):
        return rg.record(value)
    if isinstance(value, dict):
        return {key: _bound_record(item) for key, item in value.items()
                if key not in ('source_packet', 'source_family', 'source_view')}
    if isinstance(value, (list, tuple)):
        return [_bound_record(item) for item in value]
    if hasattr(value, 'record'):
        return value.record()
    return value


def _find_parent_packet(owner, ends, chart, left, right):
    for entry in owner.parent._packets.values():
        state = entry[1]
        view = state.get('view') if isinstance(state, dict) else None
        if (isinstance(view, dict) and tuple(view.get('exact_Z_cell', ())) == tuple(ends)
                and view.get('chart') == chart
                and len(view.get('exact_native_cell', ())) == 2
                and _same_endpoint(view['exact_native_cell'][0], left)
                and _same_endpoint(view['exact_native_cell'][1], right)):
            return entry[0]
    return None


def _source_and_physical_bounds(owner, ends, chart, left, right, partial_packet):
    c, bd = owner.c, rg.Bounds(owner.c)
    source_view = source_rows(owner, ends, chart, left, right)
    op = owner.modified.ranges.phase.outer.owner.owner(tuple(ends))
    if source_view['Pstar'] is not op.Pstar:
        raise ValueError('Pstar source object changed')
    entry = owner.parent._packets.get(id(partial_packet))
    if entry is None or entry[0] is not partial_packet:
        raise ValueError('Retained all-chart partial mean/Q packet required')
    state = entry[1]
    if not isinstance(state, dict):
        raise ValueError('All-chart partial packet state must be a dict')
    mean, delta_V = state.get('mean'), state.get('delta_V')
    if not isinstance(mean, rg.JetBound) or not isinstance(delta_V, rg.JetBound):
        raise ValueError('Parent partial state must retain JetBound rows')
    baseline_V = _raw_ordinary_jet(c, bd, source_view['V'], 'baseline V')
    baseline_m = _raw_ordinary_jet(c, bd, source_view['m'], 'baseline m')
    S = bd.row(source_view['Pstar'])
    U = bd.scaled(bd.add(baseline_V, delta_V), S)
    m = bd.scaled(bd.add(baseline_m, mean), S)
    q = physical._physical_q_bounds(c, ends, owner.delta, U, m)
    qjet, qz = physical.Z1Bound(q[0], q[1]), physical.Z1Bound(q[1], q[2])
    u1, uz = physical.Z1Bound(U.value, U.Z), physical._bound_z1_shift(U)
    products = dict(
        I0_q_Z=qz, I0_U_q_Z=physical._bound_z1_product(bd, u1, qz),
        I0_U_squared=physical._bound_z1_product(bd, u1, u1),
        I0_U_U_Z=physical._bound_z1_product(bd, u1, uz),
        I1_q=qjet, I1_q_Z=qz,
        I1_q_squared=physical._bound_z1_product(bd, qjet, qjet),
        I1_U_q_Z=physical._bound_z1_product(bd, u1, qz),
        I1_U_squared=physical._bound_z1_product(bd, u1, u1),
        I1_U_U_Z=physical._bound_z1_product(bd, u1, uz))
    weight = owner.modified.ranges.weights(tuple(ends), chart, left, right, Fraction(0))
    width, R = bd.row(weight['original_mass']), bd.row(source_view['R'])
    measure0, measure1 = R * width, R * R * width
    integrated = {name: physical._bound_z1_scale(
        bd, value, measure0 if name.startswith('I0_') else measure1)
        for name, value in products.items()}
    return dict(source_view=source_view,
                baseline_V=baseline_V, baseline_m=baseline_m, S=S, U=U, m=m,
                q=q, q_Z=qz, U_Z=uz, original_mass=width, original_R=R,
                physical_dR_measure=measure0, physical_R_dR_measure=measure1,
                products=products, integrated=integrated,
                no_extra_Jacobian=True, ordinary_baseline_factorials_once=True,
                physical_common_scale_applied_once=True)


class CurrentOriginalN1AllChartPhysicalPressureFields:
    @source_precision
    def __init__(self, checked_all_chart_owner, require_checked=True):
        if (type(checked_all_chart_owner) is not
                all_chart.CurrentOriginalN1AllChartPartialMeanQ or
                not checked_all_chart_owner.acceptance_loaded):
            raise ValueError('Accepted checked all-chart partial mean/Q owner required')
        checked_all_chart_owner.assert_graph()
        self.parent = checked_all_chart_owner
        self.modified = checked_all_chart_owner.modified
        self.g = checked_all_chart_owner.g
        self.c = checked_all_chart_owner.c
        self.family = copy.deepcopy(checked_all_chart_owner.family)
        self.charts = tuple(checked_all_chart_owner.charts)
        self.partitions = copy.deepcopy(checked_all_chart_owner.partitions)
        self.N, self.delta, self.delta_node = (checked_all_chart_owner.N,
                                               checked_all_chart_owner.delta,
                                               checked_all_chart_owner.delta_node)
        self.hashes = dict(checked_all_chart_owner.hashes)
        for name in (all_chart.NAME, all_chart.RECEIPT,
                     Path(all_chart.__file__).name, Path(__file__).name):
            self.hashes[name] = sha(name)
        self._delta = self.delta._mpi_
        self.functions = {chart: _field_functions(self, chart) for chart in self.charts}
        self._identities = (self.parent, self.modified, self.g, self.c,
                            self.N, self.delta_node)
        self._graph_length = len(self.g.nodes)
        self._graph = fingerprint(self.g.nodes)
        self._functions = fingerprint(physical._z1_encode(self.functions))
        self._partitions = fingerprint(self.partitions)
        self._family = copy.deepcopy(self.family)
        self._hashes = copy.deepcopy(self.hashes)
        self._packets = {}
        self.acceptance_loaded = False
        if require_checked:
            checked = json.loads((HERE / RECEIPT).read_bytes())
            if (not checked.get('all_passed') or not checked.get(GATE) or
                    checked.get('source_family') != self.family or
                    any(checked.get(key, False) for key in OPEN)):
                raise ValueError('All-chart physical pressure receipt/scope differs')
            for name, digest in checked['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('All-chart physical pressure source changed ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if (self._identities != (self.parent, self.modified, self.g, self.c,
                                 self.N, self.delta_node) or
                not self.parent.acceptance_loaded or
                self.family != self._family or self.hashes != self._hashes or
                self.delta._mpi_ != self._delta or
                fingerprint(self.g.nodes[:self._graph_length]) != self._graph or
                fingerprint(physical._z1_encode(self.functions)) != self._functions or
                fingerprint(self.partitions) != self._partitions):
            raise ValueError('All-chart physical pressure graph changed')
        self.parent.assert_graph()

    @source_precision
    def range_cell(self, ends, chart, left, right):
        self.assert_graph()
        began = time.monotonic()
        ends, left, right = tuple(ends), _canonical_endpoint(left), _canonical_endpoint(right)
        partition = tuple(self.partitions[chart])
        try:
            li, ri = next(i for i, value in enumerate(partition) if _same_endpoint(value, left)), \
                     next(i for i, value in enumerate(partition) if _same_endpoint(value, right))
        except StopIteration as exc:
            raise ValueError('Exact original native partition endpoints required') from exc
        if li >= ri:
            raise ValueError('Nonempty ordered native union required')
        partial_packet = _find_parent_packet(self, ends, chart, left, right)
        if partial_packet is None:
            partial_packet = self.parent.range_cell(ends, chart, left, right)
        bounds = _source_and_physical_bounds(self, ends, chart, left, right, partial_packet)
        functions = self.functions[chart]
        key = _interval_key(left, right)
        try:
            defining = {name: functions['defining_integrals_by_cell'][name][key]
                        for name in functions['integrands']}
        except KeyError as exc:
            raise ValueError('Requested union is not an exact source-cell union') from exc
        view = dict(
            source_family=self.family, exact_Z_cell=ends, chart=chart,
            exact_native_cell=[left, right],
            exact_native_query_union=[list(pair) for pair in zip(partition[li:ri], partition[li + 1:ri + 1])],
            exact_native_partition=list(partition),
            signed_physical_field_functions=physical._z1_encode(functions),
            signed_nonradial_integrands=physical._z1_encode(functions['integrands']),
            signed_defining_integrals=physical._z1_encode(defining),
            original_source_quantity_paths=dict(
                V='original_generic_source.common_velocity_V_axial5',
                m='original_generic_source.common_own_five_histories_axial5["m"]',
                R='original_generic_source.original_physical_radius',
                S='outer.owner.owner(ends).Pstar'),
            whole_cell_physical_U_Z3_magnitude=[v.record() for v in rg.rows(bounds['U'])],
            whole_cell_physical_m_Z3_magnitude=[v.record() for v in rg.rows(bounds['m'])],
            whole_cell_physical_q_Z2_magnitude=[v.record() for v in bounds['q']],
            whole_cell_value_Z1_integrand_magnitudes={
                name: [value.value.record(), value.Z.record()]
                for name, value in bounds['products'].items()},
            whole_cell_integrated_value_Z1_magnitudes={
                name: [value.value.record(), value.Z.record()]
                for name, value in bounds['integrated'].items()},
            original_positive_dlogR_width=bounds['original_mass'].record(),
            original_R_cap=bounds['original_R'].record(),
            physical_dR_measure=bounds['physical_dR_measure'].record(),
            physical_R_dR_measure=bounds['physical_R_dR_measure'].record(),
            retained_parent_partial_packet=True,
            baseline_rows_are_source_rows_not_range_values=True,
            ordinary_baseline_factorial_conversion_once=True,
            physical_common_scale_Pstar_applied_once=True,
            physical_radius_original_geometry=True,
            whole_cell_range_not_endpoint_value=True,
            exact_original_source_cell_union=True,
            numerical_evidence_absolute_whole_cell_magnitude_only=True,
            no_signed_numerical_centers_or_relative_width=True,
            no_finite_endpoint_q_R_or_support_claim=True,
            complete_pressure_tail_support_repair_n2_and_global_claims_not_installed=True,
            execution_seconds=time.monotonic() - began,
            **{GATE: self.acceptance_loaded}, **dict.fromkeys(OPEN, False))
        packet = OriginalAllChartPhysicalPressurePacket()
        state = dict(view=view, bounds=bounds, functions=functions,
                     partial_packet=partial_packet, defining_integrals=defining)
        self._packets[id(packet)] = (packet, state,
            fingerprint((view, _bound_record(bounds), physical._z1_encode(functions))))
        return packet

    @source_precision
    def range_chart(self, ends, chart):
        partition = self.partitions[chart]
        return self.range_cell(ends, chart, partition[0], partition[-1])

    @source_precision
    def report(self, packet):
        self.assert_graph()
        entry = self._packets.get(id(packet))
        if (type(packet) is not OriginalAllChartPhysicalPressurePacket or
                entry is None or entry[0] is not packet):
            raise ValueError('Live all-chart physical pressure packet required')
        state = entry[1]
        self.parent.report(state['partial_packet'])
        if fingerprint((state['view'], _bound_record(state['bounds']),
                        physical._z1_encode(state['functions']))) != entry[2]:
            raise ValueError('All-chart physical pressure packet/source ranges changed')
        return {**copy.deepcopy(state['view']), GATE: self.acceptance_loaded}


@source_precision
def run(checked_all_chart_owner=None):
    began = time.monotonic()
    with mp.workdps(540):
        if checked_all_chart_owner is None:
            modified = partial.parent.CurrentOriginalN1ModifiedLeadingC3(
                rg.CurrentC3TargetRanges(), require_checked=True)
            checked_partial = partial.CurrentOriginalN1PartialMeanQ(
                modified, require_checked=True)
            checked_all_chart_owner = all_chart.CurrentOriginalN1AllChartPartialMeanQ(
                checked_partial, require_checked=True)
        owner = CurrentOriginalN1AllChartPhysicalPressureFields(
            checked_all_chart_owner, require_checked=False)
        packets = {chart: owner.range_chart(('0', '.5'), chart)
                   for chart in owner.charts}
        samples = {chart: {key: value for key, value in owner.report(packet).items()
                           if key not in ('signed_physical_field_functions',
                                          'signed_nonradial_integrands',
                                          'signed_defining_integrals')}
                   for chart, packet in packets.items()}
        raw = dict(candidate_all_chart_physical_pressure_fields_constructed=True,
                   source_family=owner.family, exact_Z_cell=('0', '.5'),
                   charts=list(owner.charts), samples=samples,
                   signed_function_charts=list(owner.functions),
                   original_graph_prefix_length=owner._graph_length,
                   input_hashes=owner.hashes, **{GATE: False},
                   **dict.fromkeys(OPEN, False),
                   execution_seconds=time.monotonic() - began)
        (HERE / NAME).write_text(
            json.dumps(c3.encoded(raw), separators=(',', ':')) + '\n',
            encoding='utf-8', newline='\n')
        print('ALL_CHART_PHYSICAL_PRESSURE_FIELDS_COMPUTED', len(packets), flush=True)
        return owner, packets


if __name__ == '__main__':
    run()
