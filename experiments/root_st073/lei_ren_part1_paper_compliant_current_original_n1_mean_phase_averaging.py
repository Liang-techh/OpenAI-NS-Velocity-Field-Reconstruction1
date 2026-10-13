"""Signed phase averaging for the accepted finite-N own-mean source.

The accepted all-chart owner already contains the exact moving-query Duhamel
mean and directed fallback caps.  This adapter adds only signed graph
functions: the original B density is given a periodic phase primitive and the
own-mean integral is rewritten by endpoint-retaining integration by parts.
No slow-y cap is fabricated.  Since the current original source packets do
not expose a genuine mixed slow-y source bound, numerical averaging remains
disabled and the parent whole-cell caps remain the fallback.
"""
import copy
import json
import math
from dataclasses import dataclass
from pathlib import Path

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as atlas
import lei_ren_part1_paper_compliant_current_native_signed_averaging as averaging


partial = atlas.partial
c3, current, source = atlas.c3, atlas.current, atlas.source
rg = atlas.rg
HERE, PREFIX, sha = atlas.HERE, atlas.PREFIX, atlas.sha
source_precision = atlas.source_precision
fingerprint = atlas.fingerprint
OPEN = tuple(dict.fromkeys(atlas.OPEN))

NAME = PREFIX + 'current_original_n1_mean_phase_averaging.json'
RECEIPT = PREFIX + 'current_original_n1_mean_phase_averaging_check.json'
GATE = 'current_original_n1_mean_phase_averaging_signed_IBP_installed'
NUMERIC_GATE = 'current_original_n1_mean_phase_averaging_numerical_cap_installed'


@dataclass(frozen=True, eq=False)
class OriginalMeanPhaseAveragingPacket:
    leading_hierarchy_order: int = 0


def _as_ref(g, value):
    return value if isinstance(value, source.FunctionRef) else source.FunctionRef(g, value)


def _as_endpoint(value):
    """Normalize JSON list syntax without interpreting symbolic endpoints."""
    if isinstance(value, list):
        return tuple(_as_endpoint(item) for item in value)
    return value


def _endpoint_index(partition, endpoint):
    endpoint = _as_endpoint(endpoint)
    for index, candidate in enumerate(partition):
        if endpoint == _as_endpoint(candidate):
            return index
    raise ValueError('Exact original native partition endpoint required')


def _target_window(modified, chart):
    return next(value for value in modified.target.functions['windows']
                if value['chart'] == chart)


def _at_c3(g, value, chart, query):
    return c3.C3Function(*(partial.at_native(g, _as_ref(g, row), chart, query)
                           for row in c3.rows(value)))


def _periodic_primitive_rows(g, owner, chart, B):
    base = owner.modified.target.phase.windows[chart]
    phase_node = base.get('phase')
    rows = []
    for order, row in enumerate(c3.rows(B)):
        rows.append(g.node(
            'current_original_n1_periodic_phase_primitive',
            body=row.node, native_variable='native_' + chart,
            phase_variable='original_global_phase', source_phase_node=phase_node,
            period='1', phase_endpoint_zero_mean=True,
            primitive_endpoint_identity='G(0)=G(1)=0 from accepted original B parity',
            integrand='original B density row %d' % order,
            source_family=owner.family, chart=chart, Z_order=order,
            selected_N=owner.N.node, same_original_P0_and_global_phase=True,
            source_zero_mean_receipt=averaging.RECEIPT,
            defining_source=current.ast_binding(_periodic_primitive_rows),
            signed_function_not_a_range_endpoint=True))
    return c3.C3Function(*rows)


def _slow_y_rows(g, owner, chart, G):
    base = owner.modified.target.phase.windows[chart]
    phase_node = base.get('phase')
    rows = []
    for order, row in enumerate(c3.rows(G)):
        rows.append(g.node(
            'current_original_n1_periodic_phase_primitive_slow_y',
            body=row.node, native_variable='native_' + chart,
            phase_variable='original_global_phase', source_phase_node=phase_node,
            fixed_phase=True, ordinary_slow_y_derivative_order=1,
            ordinary_Z_order=order, total_y_chain_excluded=True,
            source_family=owner.family, chart=chart, selected_N=owner.N.node,
            genuine_source_cap_required=True,
            source_rows_provider=current.ast_binding(_slow_y_rows),
            signed_function_not_a_range_endpoint=True))
    return c3.C3Function(*rows)


def _ibp_row_from_handles(g, chart, query, offset, left_offset, G, Gy,
                          jacobian, lower, inverse_N, order, rate=1):
    """Production moving-query IBP row, with a distinct outer native variable."""
    outer_name = 'mean_phase_native_' + chart
    outer = g.symbol(outer_name)
    offset_query = partial.at_native(g, _as_ref(g, offset), chart, query)
    offset_outer = partial.at_native(g, _as_ref(g, offset), chart, outer)
    offset_left = partial.at_native(g, _as_ref(g, offset), chart, _as_ref(g, lower))
    G_query = partial.at_native(g, _as_ref(g, G), chart, query)
    G_outer = partial.at_native(g, _as_ref(g, G), chart, outer)
    G_left = partial.at_native(g, _as_ref(g, G), chart, _as_ref(g, lower))
    Gy_outer = partial.at_native(g, _as_ref(g, Gy), chart, outer)
    J_outer = partial.at_native(g, _as_ref(g, jacobian), chart, outer)
    kernel_left = g.unary('exp', g.neg(g.sub(offset_query, offset_left)))
    kernel_outer = g.unary('exp', g.neg(g.sub(offset_query, offset_outer)))
    rate_ref = g.constant(rate)
    endpoint = g.sub(G_query, g.mul(kernel_left, G_left))
    body = g.mul(kernel_outer, g.add(Gy_outer, g.mul(rate_ref, G_outer)), J_outer)
    integral = current.controls.integral(
        g, body, outer_name, _as_ref(g, lower), query,
        measure='original dlogR: Jacobian exactly once',
        outer_integration_variable=outer_name,
        ordinary_slow_Z_derivative_order=order,
        free_query_endpoint=query.node,
        free_query_substituted_at_outer_native=True,
        bound_inner_mean_native_variable='native_' + chart,
        bound_phase_variable='original_global_phase',
        original_phase_primitive_slow_y_not_total_y=True,
        endpoint_terms_retained=True,
        same_original_source_and_global_phase=True)
    # ∫ K B/N = ([K G] - ∫ K(G_y+lambda G) dy)/N^2.
    return g.mul(g.mul(inverse_N, inverse_N), g.sub(endpoint, integral))


def _phase_certificate(owner):
    checked = json.loads((HERE / averaging.RECEIPT).read_bytes())
    if (not checked.get('all_passed') or not checked.get(averaging.GATE) or
            checked.get('source_family') != owner.family):
        raise ValueError('Accepted same-family original zero-phase theorem required')
    theorem = averaging.exact_theorem()
    if not theorem.get('uniform_fractional_phi_mean_of_leading_and_fixed_phase_slow_jets_exact_zero'):
        raise ValueError('Original B zero-phase theorem is not available')
    return theorem


def _source_packet(owner, ends, chart, left, right=None):
    packet = owner.modified.target.phase.leading_source(tuple(ends), chart, left, right)
    op = owner.modified.ranges.phase.outer.owner.owner(tuple(ends))
    if (packet.get('source_identity') != owner.family or
            packet.get('exact_common_P0_axial5') is not op.P0 or
            packet.get('actual_phase_Z_exact_zero') is not True):
        raise ValueError('Same accepted source/P0 and global phase required')
    generic = packet.get('original_generic_source')
    if (not isinstance(generic, dict) or
            'common_velocity_V_axial5' not in generic or
            'common_own_five_histories_axial5' not in generic):
        raise ValueError('Actual original source packet missing')
    return packet, op


def _field_functions(owner, chart):
    g = owner.g
    base = owner.modified.target.phase.windows[chart]
    window = _target_window(owner.modified, chart)
    query = owner.parent.functions[chart]['query']
    B = owner.modified.target.primitives[chart]['B']
    G = _periodic_primitive_rows(g, owner, chart, B)
    Gy = _slow_y_rows(g, owner, chart, G)
    offset = _as_ref(g, base['offset'])
    left_offset = _as_ref(g, base['left_offset'])
    query_offset = partial.at_native(g, offset, chart, query)
    memory = g.unary('exp', g.neg(g.sub(query_offset, left_offset)))
    inverse_N = g.quotient(g.one, owner.N, 'same selected exact positive N')
    incoming = window['incoming']['m']
    lower = _as_ref(g, base['domain'][0])
    rows_mean = []
    rows_ibp = []
    for order, (inlet, primitive, slow) in enumerate(zip(
            c3.rows(incoming), c3.rows(G), c3.rows(Gy))):
        ibp = _ibp_row_from_handles(
            g, chart, query, offset, left_offset, primitive, slow,
            base['Jacobian'], lower, inverse_N, order, rate=1)
        rows_ibp.append(ibp)
        rows_mean.append(g.add(g.mul(inverse_N, memory, inlet), ibp))
    averaged = c3.C3Function(*rows_mean)
    ibp = c3.C3Function(*rows_ibp)
    direct = owner.parent.functions[chart]['mean']
    return dict(query=query, B=B, G=G, G_slow_y=Gy, mean_ibp=averaged,
                ibp_contribution=ibp, direct_mean=direct, incoming=incoming,
                memory=memory, selected_N=owner.N,
                source_geometry=dict(
                    original_offset=base['offset'],
                    original_left_offset=base['left_offset'],
                    original_Jacobian=base['Jacobian'],
                    Jacobian_bound_at_outer_native=True,
                    inner_mean_native_variable='native_' + chart,
                    outer_IBP_native_variable='mean_phase_native_' + chart,
                    phase_variable='original_global_phase',
                    physical_measure='Jacobian exactly once'))


def _bound_record(value):
    if isinstance(value, dict):
        return {key: _bound_record(item) for key, item in value.items()
                if key not in ('source_packet', 'source_family')}
    if isinstance(value, (list, tuple)):
        return [_bound_record(item) for item in value]
    if isinstance(value, rg.JetBound):
        return rg.record(value)
    if hasattr(value, 'record'):
        return value.record()
    return value


def _packet_bound_record(state):
    """Fingerprint only retained quantitative state, never live source handles."""
    return dict(
        mean=_bound_record(state['state']['mean']),
        mean_y=_bound_record(state['state']['mean_y']),
        delta_V=_bound_record(state['state']['delta_V']),
        Q=_bound_record(state['state']['Q']),
    )


class CurrentOriginalN1MeanPhaseAveraging:
    @source_precision
    def __init__(self, checked_all_chart_owner, require_checked=True):
        if (type(checked_all_chart_owner) is not
                atlas.CurrentOriginalN1AllChartPartialMeanQ or
                not checked_all_chart_owner.acceptance_loaded):
            raise ValueError('Accepted all-chart partial mean/Q owner required')
        checked_all_chart_owner.assert_graph()
        self.parent = checked_all_chart_owner
        self.modified = checked_all_chart_owner.modified
        self.g, self.c = checked_all_chart_owner.g, checked_all_chart_owner.c
        self.family = copy.deepcopy(checked_all_chart_owner.family)
        self.charts = tuple(checked_all_chart_owner.charts)
        self.partitions = copy.deepcopy(checked_all_chart_owner.partitions)
        self.N = checked_all_chart_owner.N
        self.delta = checked_all_chart_owner.delta
        self.theorem = _phase_certificate(self)
        self.hashes = dict(checked_all_chart_owner.hashes)
        for name in (atlas.NAME, atlas.RECEIPT, Path(atlas.__file__).name,
                     averaging.NAME, averaging.RECEIPT,
                     Path(averaging.__file__).name, Path(__file__).name):
            self.hashes[name] = sha(name)
        self.functions = {chart: _field_functions(self, chart)
                          for chart in self.charts}
        self._identities = (self.parent, self.modified, self.g, self.c, self.N)
        self._graph_length = len(self.g.nodes)
        self._graph = fingerprint(self.g.nodes)
        self._functions = fingerprint(c3.encoded(self.functions))
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
                raise ValueError('Mean phase averaging receipt/scope differs')
            for name, digest in checked['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('Mean phase averaging source changed ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if (self._identities != (self.parent, self.modified, self.g, self.c, self.N) or
                not self.parent.acceptance_loaded or self.family != self._family or
                self.hashes != self._hashes or
                fingerprint(self.g.nodes[:self._graph_length]) != self._graph or
                fingerprint(c3.encoded(self.functions)) != self._functions or
                fingerprint(self.partitions) != self._partitions):
            raise ValueError('Mean phase averaging graph changed')
        self.parent.assert_graph()

    @source_precision
    def range_cell(self, ends, chart, left, right):
        self.assert_graph()
        if chart not in self.partitions:
            raise ValueError('Accepted original chart required')
        partition = tuple(_as_endpoint(item) for item in self.partitions[chart])
        left = _as_endpoint(left)
        right = _as_endpoint(right)
        left_index = _endpoint_index(partition, left)
        right_index = _endpoint_index(partition, right)
        if left_index >= right_index:
            raise ValueError('Exact contiguous original source-cell union required')
        source_packet, op = _source_packet(self, ends, chart, left, right)
        parent_packet = self.parent.range_cell(ends, chart, left, right)
        entry = self.parent._packets.get(id(parent_packet))
        if entry is None or entry[0] is not parent_packet:
            raise ValueError('Retained parent all-chart packet required')
        state = entry[1]
        key = (left, right)
        view = dict(
            source_family=self.family, exact_Z_cell=tuple(ends), chart=chart,
            exact_native_cell=[left, right], exact_native_partition=list(partition),
            signed_phase_primitive_Z0_to_Z3=c3.encoded(self.functions[chart]['G']),
            signed_slow_y_phase_primitive_Z0_to_Z3=c3.encoded(self.functions[chart]['G_slow_y']),
            signed_moving_query_IBP_mean_Z0_to_Z3=c3.encoded(self.functions[chart]['mean_ibp']),
            signed_direct_parent_mean_Z0_to_Z3=c3.encoded(self.functions[chart]['direct_mean']),
            source_quantity_paths=dict(
                B='target.primitives[chart].B / original finite-N m density',
                incoming='target.functions.windows[chart].incoming.m',
                P0='outer.owner.owner(ends).P0'),
            selected_N_node=self.N.node,
            source_phase_zero_mean_theorem=averaging.RECEIPT,
            original_inlet_memory_retained=True,
            endpoint_terms_retained=True,
            same_N_squared_IBP_factor=True,
            original_Jacobian_applied_once=True,
            inner_mean_native_variable='native_' + chart,
            outer_IBP_native_variable='mean_phase_native_' + chart,
            original_global_phase_retained=True,
            genuine_slow_y_source_cap_available=False,
            numerical_phase_averaging_gate=False,
            numerical_phase_averaging_gate_reason=(
                'Accepted source packets expose signed B and ordinary Z3 rows, '
                'but no genuine fixed-phase mixed slow-y source cap; direct '
                'whole-cell mean/mean_y/Q bounds remain the only numerical fallback.'),
            numerical_fallback_direct_old_whole_mean_bounds=True,
            fallback_parent_mean_Z3_magnitude=rg.record(state['mean']),
            fallback_parent_mean_y_Z3_magnitude=rg.record(state['mean_y']),
            fallback_parent_Q_Z2_magnitude=[value.record() for value in state['Q']],
            source_packet_same_family_P0=True,
            complete_signed_defects_pressure_repair_and_global_claims_not_installed=True,
            **{GATE: self.acceptance_loaded},
            **{NUMERIC_GATE: False}, **dict.fromkeys(OPEN, False))
        packet = OriginalMeanPhaseAveragingPacket()
        packet_state = dict(view=view, parent_packet=parent_packet,
                            source_packet=source_packet, state=state)
        self._packets[id(packet)] = (packet, packet_state,
                                      fingerprint((view, c3.encoded(self.functions[chart]),
                                                   _packet_bound_record(packet_state))))
        return packet

    @source_precision
    def range_chart(self, ends, chart):
        partition = self.partitions[chart]
        return self.range_cell(ends, chart, partition[0], partition[-1])

    @source_precision
    def report(self, packet):
        self.assert_graph()
        entry = self._packets.get(id(packet))
        if (type(packet) is not OriginalMeanPhaseAveragingPacket or
                entry is None or entry[0] is not packet):
            raise ValueError('Live mean phase averaging packet required')
        state = entry[1]
        self.parent.report(state['parent_packet'])
        if fingerprint((state['view'], c3.encoded(self.functions[state['view']['chart']]),
                        _packet_bound_record(state))) != entry[2]:
            raise ValueError('Mean phase averaging packet changed')
        return {**copy.deepcopy(state['view']), GATE: self.acceptance_loaded,
                NUMERIC_GATE: False}


@source_precision
def run(checked_all_chart_owner=None):
    if checked_all_chart_owner is None:
        modified = partial.parent.CurrentOriginalN1ModifiedLeadingC3(
            rg.CurrentC3TargetRanges(), require_checked=True)
        checked_partial = partial.CurrentOriginalN1PartialMeanQ(
            modified, require_checked=True)
        checked_all_chart_owner = atlas.CurrentOriginalN1AllChartPartialMeanQ(
            checked_partial, require_checked=True)
    owner = CurrentOriginalN1MeanPhaseAveraging(
        checked_all_chart_owner, require_checked=False)
    ends = ('0', '.5')
    chart = 'first_micro'
    partition = owner.partitions[chart]
    packet = owner.range_cell(ends, chart, partition[0], partition[1])
    raw = dict(candidate_mean_phase_averaging_constructed=True,
               source_family=owner.family, exact_Z_cell=ends,
               actual_chart=chart, sample=owner.report(packet),
               signed_function_charts=list(owner.functions),
               original_graph_prefix_length=owner._graph_length,
               input_hashes=owner.hashes, **{GATE: False},
               **{NUMERIC_GATE: False}, **dict.fromkeys(OPEN, False))
    (HERE / NAME).write_text(json.dumps(c3.encoded(raw), separators=(',', ':')) + '\n',
                             encoding='utf-8', newline='\n')
    print('N1_MEAN_PHASE_AVERAGING_COMPUTED', chart, flush=True)
    return owner, {chart: packet}


if __name__ == '__main__':
    run()
