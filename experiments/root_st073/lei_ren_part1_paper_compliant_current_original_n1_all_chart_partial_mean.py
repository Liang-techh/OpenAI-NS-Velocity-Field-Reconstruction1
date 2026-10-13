"""All-chart finite-N own-mean and meridional-Q increment ranges.

This is a small adapter around the accepted finite-N partial-mean owner.  It
does not rebuild the source graph or manufacture fields for the outer charts.
It only projects the accepted signed functions onto exact unions of the
original native partition cells and applies the saved C3 whole-Z transport
caps to those cells.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q as partial


c3, current, source = partial.c3, partial.current, partial.source
rg = partial.rg
HERE, PREFIX, sha, ep = partial.HERE, partial.PREFIX, partial.sha, partial.ep
source_precision = partial.source_precision
fingerprint = partial.fingerprint
OPEN = partial.OPEN

NAME = PREFIX + 'current_original_n1_all_chart_partial_mean.json.gz'
RECEIPT = PREFIX + 'current_original_n1_all_chart_partial_mean_check.json'
GATE = 'current_original_finite_N_all_17_chart_partial_own_mean_Z3_mean_y_Z3_and_Q_increment_Z2_installed'


def _as_endpoint(value):
    """Normalize only JSON list syntax; never interpret an endpoint numerically."""
    return tuple(value) if isinstance(value, list) else value


def _endpoint_index(partition, endpoint):
    endpoint = _as_endpoint(endpoint)
    for index, candidate in enumerate(partition):
        if endpoint == candidate:
            return index
    raise ValueError('Exact original native partition endpoint required')


def _partition_manifest(partitions, charts):
    return [dict(chart=chart, exact_native_partition=list(partitions[chart]),
                 source_cell_count=len(partitions[chart]) - 1,
                 symbolic_endpoint_tags=[q for q in partitions[chart]
                                         if q == ('sc_half',) or q == 'Rh'])
            for chart in charts]


def _row_is_nonzero(row):
    return not bool(row.get('exact_zero', False))


@dataclass(frozen=True, eq=False)
class OriginalAllChartPartialMeanQPacket:
    leading_hierarchy_order: int = 0


class CurrentOriginalN1AllChartPartialMeanQ:
    """Checked all-chart adapter for the accepted partial own mean.

    ``range_cell`` accepts only exact adjacent partition endpoints or their
    exact contiguous unions.  The returned magnitude caps are for the query
    union; the own-mean prefix always starts at the original chart inlet.
    """

    @source_precision
    def __init__(self, checked_partial_owner, require_checked=True):
        if type(checked_partial_owner) is not partial.CurrentOriginalN1PartialMeanQ:
            raise ValueError('Checked CurrentOriginalN1PartialMeanQ owner required')
        if not checked_partial_owner.acceptance_loaded:
            raise ValueError('Accepted finite-N partial mean/Q owner required')
        checked_partial_owner.assert_graph()

        self.partial = checked_partial_owner
        self.parent = checked_partial_owner
        self.modified = checked_partial_owner.parent
        self.ranges = self.modified.ranges
        self.saved = self.modified.saved
        self.g = checked_partial_owner.g
        self.graph = self.g
        self.c = checked_partial_owner.c
        self.N = checked_partial_owner.N
        self.delta = checked_partial_owner.delta
        self.delta_node = checked_partial_owner.delta_node
        self.logN = checked_partial_owner.logN
        self.epsilon_N = checked_partial_owner.epsilon_N
        self.family = copy.deepcopy(checked_partial_owner.family)
        if self.modified.family != self.family or self.ranges.identity != self.family:
            raise ValueError('All-chart adapter source family changed')
        if self.saved.get('source_family') != self.family:
            raise ValueError('Saved C3 transport source family changed')

        # These are the accepted signed source functions.  Keep the same live
        # graph and function handles; all new work below is range arithmetic.
        self.functions = checked_partial_owner.functions
        self.charts = tuple(self.functions)
        self.partitions = current.native_partitions()
        if set(self.charts) != set(self.partitions) or len(self.charts) != 17:
            raise ValueError('All 17 accepted original charts are required')
        if tuple(self.partitions) != tuple(current.native_partitions()):
            raise ValueError('Original native partition ordering changed')

        self.hashes = dict(checked_partial_owner.hashes)
        for name in (
            partial.parent.NAME, partial.parent.RECEIPT,
            Path(partial.parent.__file__).name,
            partial.NAME, partial.RECEIPT, Path(partial.__file__).name,
            Path(__file__).name,
        ):
            self.hashes[name] = sha(name)

        self._identities = (self.partial, self.modified, self.ranges, self.g,
                            self.c, self.N, self.delta_node)
        self._family = copy.deepcopy(self.family)
        self._hashes = copy.deepcopy(self.hashes)
        self._graph_length = len(self.g.nodes)
        self._graph = fingerprint(self.g.nodes)
        self._functions = fingerprint(self.functions)
        self._saved = fingerprint(self.saved)
        self._partitions = fingerprint(self.partitions)
        self._delta = self.delta._mpi_
        self._logN = self.logN._mpi_
        self._epsilon = fingerprint(self.epsilon_N.record())
        self._packets = {}
        self.acceptance_loaded = False
        if require_checked:
            checked = json.loads((HERE / RECEIPT).read_bytes())
            if (not checked.get('all_passed') or not checked.get(GATE)
                    or checked.get('source_family') != self.family
                    or any(checked.get(key, False) for key in OPEN)):
                raise ValueError('All-chart partial mean/Q receipt or scope differs')
            for name, digest in checked['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('All-chart partial mean/Q source changed ' + name)
                self.hashes[name] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)

    @source_precision
    def assert_graph(self):
        if (
            (self.partial, self.modified, self.ranges, self.g, self.c,
             self.N, self.delta_node) != self._identities
            or not self.partial.acceptance_loaded
            or not self.modified.acceptance_loaded
            or not self.ranges.acceptance_loaded
            or self.family != self._family
            or self.family != self.partial.family
            or self.family != self.modified.family
            or self.ranges.identity != self.family
            or self.hashes != self._hashes
            or self.delta._mpi_ != self._delta
            or self.logN._mpi_ != self._logN
            or fingerprint(self.epsilon_N.record()) != self._epsilon
            or len(self.g.nodes) < self._graph_length
            or fingerprint(self.g.nodes[:self._graph_length]) != self._graph
            or fingerprint(self.functions) != self._functions
            or fingerprint(self.saved) != self._saved
            or fingerprint(self.partitions) != self._partitions
            or tuple(self.functions) != self.charts
            or set(self.charts) != set(current.native_partitions())
        ):
            raise ValueError('All-chart partial mean/Q source graph changed')
        self.partial.assert_graph()

    def _transport(self, ends, chart):
        ends = tuple(ends)
        transport = next((value for value in
                          self.saved['actual_four_Z_C3_target_transports']
                          if tuple(value['exact_Z_cell']) == ends), None)
        if transport is None:
            raise ValueError('Original saved whole-Z transport cell required')
        if transport.get('source_identity') != self.family:
            raise ValueError('Original saved transport family changed')
        stored = next((value for value in
                       transport['actual_17_chart_C3_transports']
                       if value['chart'] == chart), None)
        if stored is None:
            raise ValueError('Original all-chart transport window required')
        if stored.get('global_phase_cover') != [0, 1]:
            raise ValueError('Whole phase cover required')
        if not stored.get('third_row_predecessor_never_reset', False):
            raise ValueError('Incoming C3 history must not reset')
        incoming = stored.get('actual_N_scaled_incoming_C3', {}).get('m')
        if not isinstance(incoming, list) or len(incoming) != 4:
            raise ValueError('Incoming N-scaled C3 own-mean history required')
        return transport, stored

    def _live_cell(self, ends, chart, left, right):
        cell = self.ranges.cell(tuple(ends), chart, left, right)
        weight = self.ranges.weights(tuple(ends), chart, left, right, Fraction(0))
        return cell, weight

    def _incoming_bound(self, stored):
        bd = rg.Bounds(self.c)

        def read(row):
            if row.get('exact_zero', False):
                return bd.zero
            return rg.LogUpper(self.c,
                               current.packets.interval(self.c,
                                                        row['log_absolute_upper']))

        return rg.JetBound(*(read(row) for row in
                             stored['actual_N_scaled_incoming_C3']['m']))

    @staticmethod
    def _source_payload(source_cells):
        return [dict(chart=item['chart'], exact_left=item['exact_left'],
                     exact_right=item['exact_right'],
                     original_B_Z3_magnitude=rg.record(item['cell']['primitive']['B']),
                     original_positive_logR_width=rg.record(item['width']))
                for item in source_cells]

    def _packet_payload(self, state):
        return (
            state['view'], self._source_payload(state['source_cells']),
            rg.record(state['local_B']), rg.record(state['mean']),
            rg.record(state['mean_y']), [rg.record(value) for value in state['Q']]
        )

    @source_precision
    def range_cell(self, ends, chart, left, right):
        self.assert_graph()
        began = time.monotonic()
        ends = tuple(ends)
        if chart not in self.partitions:
            raise ValueError('Accepted original chart required')
        partition = self.partitions[chart]
        left = _as_endpoint(left)
        right = _as_endpoint(right)
        left_index = _endpoint_index(partition, left)
        right_index = _endpoint_index(partition, right)
        if left_index >= right_index:
            raise ValueError('Nonempty ordered native partition union required')
        _transport, stored = self._transport(ends, chart)

        bd = rg.Bounds(self.c)
        incoming = self._incoming_bound(stored)
        scaled_mean = incoming
        source_cells = []
        # The mean at the right endpoint sees every preceding original cell,
        # with the rate-one kernel and incoming memory bounded by one.
        for index, (cell_left, cell_right) in enumerate(
                zip(partition, partition[1:])):
            if index >= right_index:
                break
            cell, weight = self._live_cell(ends, chart, cell_left, cell_right)
            width = bd.row(weight['original_mass'])
            scaled_mean = bd.add(scaled_mean,
                                 bd.scaled(cell['primitive']['B'], width))
            source_cells.append(dict(chart=chart, exact_left=cell_left,
                                     exact_right=cell_right, cell=cell,
                                     width=width))

        mean = bd.scaled(scaled_mean, self.epsilon_N)
        query_cells = source_cells[left_index:right_index]
        if not query_cells:
            raise ValueError('Query union did not contain an original cell')
        local_B = bd.add(*(item['cell']['primitive']['B'] for item in query_cells))
        delta_V = bd.scaled(local_B, self.epsilon_N)
        # |delta_V - mean| is enclosed by the sum of the two directed caps.
        mean_y = bd.add(delta_V, mean)
        Q = partial.coefficient_Q_bounds(self.c, ends, self.delta,
                                         delta_V, mean)

        history_cells = [dict(exact_left=item['exact_left'],
                              exact_right=item['exact_right'],
                              original_positive_logR_width=rg.record(item['width']),
                              original_B_Z3_magnitude=rg.record(item['cell']['primitive']['B']),
                              rate_one_kernel_upper='<=1 for source before query',
                              incoming_decay_upper='<=1; predecessor history retained')
                         for item in source_cells]
        query_manifest = [dict(exact_left=item['exact_left'],
                               exact_right=item['exact_right'])
                          for item in query_cells]
        view = dict(
            source_family=self.family,
            exact_Z_cell=ends,
            chart=chart,
            exact_native_cell=[left, right],
            exact_native_query_union=query_manifest,
            exact_native_partition=list(partition),
            whole_chart_scope=dict(
                complete_original_partition=True,
                exact_native_partition=list(partition),
                source_cell_count=len(partition) - 1,
                query_is_contiguous_union_of_original_cells=True,
                own_mean_prefix_starts_at_original_inlet=True,
                incoming_history_is_from_saved_whole_Z_transport=True,
                original_phase_and_Jacobian_Z_independent=True,
            ),
            signed_partial_mean_query_functions=dict(
                mean=c3.encoded(self.functions[chart]['mean']),
                mean_y=c3.encoded(self.functions[chart]['mean_y']),
                delta_V=c3.encoded(self.functions[chart]['delta_V']),
                original_N_scaled_mean=c3.encoded(
                    self.functions[chart]['original_N_scaled_mean']),
            ),
            signed_Q_increment_Z2_function=c3.encoded(
                self.functions[chart]['delta_Q_Z2']),
            whole_cell_partial_own_mean_increment_Z3_magnitude=rg.record(mean),
            whole_cell_mean_y_increment_Z3_magnitude=rg.record(mean_y),
            whole_cell_delta_V_increment_Z3_magnitude=rg.record(delta_V),
            whole_cell_Q_increment_Z2_magnitude=[value.record() for value in Q],
            selected_N_log=self.logN,
            epsilon_N_log_bound=self.epsilon_N.record(),
            original_delta_recipe='exp(-4*(exp(40)+11)-30)',
            original_incoming_N_scaled_mean_Z3=copy.deepcopy(
                stored['actual_N_scaled_incoming_C3']['m']),
            partial_history_cells=history_cells,
            common_unit='physical increments are S*delta_m and S*delta_Q; original common S retained',
            same_selected_N_division_once_no_integer_materialization=True,
            partial_mean_free_query_is_distinct_from_bound_integration_coordinate=True,
            mean_y_is_delta_V_minus_mean_enclosed_by_sum_of_magnitudes=True,
            whole_cell_range_not_outgoing_endpoint_value=True,
            modified_Q_increment_only_not_complete_Q=True,
            Q_only_value_Z_ZZ_rows_no_fabricated_Q_ZZZ=True,
            internal_radial_phase_derivative_of_velocity_not_exported=True,
            complete_pressure_source_tail_support_and_repair_not_installed=True,
            magnitude_bounds_not_point_values_and_no_cap_differentiation=True,
            execution_seconds=time.monotonic() - began,
            **{GATE: True}, **dict.fromkeys(OPEN, False),
        )
        packet = OriginalAllChartPartialMeanQPacket()
        state = dict(view=view, source_cells=source_cells, local_B=local_B,
                     mean=mean, mean_y=mean_y, delta_V=delta_V, Q=Q)
        self._packets[id(packet)] = (packet, state,
                                     fingerprint(self._packet_payload(state)))
        return packet

    @source_precision
    def range_chart(self, ends, chart):
        partition = self.partitions[chart]
        return self.range_cell(ends, chart, partition[0], partition[-1])

    @source_precision
    def all_cells(self, ends, chart):
        partition = self.partitions[chart]
        return [self.range_cell(ends, chart, left, right)
                for left, right in zip(partition, partition[1:])]

    @source_precision
    def report(self, packet):
        self.assert_graph()
        entry = self._packets.get(id(packet))
        if (type(packet) is not OriginalAllChartPartialMeanQPacket
                or entry is None or entry[0] is not packet):
            raise ValueError('Live all-chart partial mean/Q packet required')
        if fingerprint(self._packet_payload(entry[1])) != entry[2]:
            raise ValueError('All-chart partial mean/Q packet or source range changed')
        return {**copy.deepcopy(entry[1]['view']), GATE: self.acceptance_loaded}


# Short alias retained for downstream adapters that name the object without Q.
CurrentOriginalN1AllChartPartialMean = CurrentOriginalN1AllChartPartialMeanQ


@source_precision
def run(checked_partial_owner=None):
    """Build one exact native cell per original chart at Z cell [0, 1/2]."""
    if checked_partial_owner is None:
        modified = partial.parent.CurrentOriginalN1ModifiedLeadingC3(
            rg.CurrentC3TargetRanges())
        checked_partial_owner = partial.CurrentOriginalN1PartialMeanQ(
            modified, require_checked=True)
    owner = CurrentOriginalN1AllChartPartialMeanQ(checked_partial_owner,
                                                  require_checked=False)
    packets = {}
    ends = ('0', '.5')
    for chart in owner.charts:
        partition = owner.partitions[chart]
        packets[chart] = owner.range_cell(ends, chart, partition[0], partition[1])
    samples = {}
    for chart, packet in packets.items():
        view = owner.report(packet)
        samples[chart] = {key: value for key, value in view.items()
                          if key not in ('signed_partial_mean_query_functions',
                                         'signed_Q_increment_Z2_function')}
    raw = dict(
        source_family=owner.family,
        exact_Z_cell=ends,
        charts=list(owner.charts),
        whole_chart_partition_manifest=_partition_manifest(owner.partitions,
                                                           owner.charts),
        samples=samples,
        signed_function_charts=list(owner.functions),
        input_hashes=owner.hashes,
        **{GATE: False}, **dict.fromkeys(OPEN, False),
    )
    (HERE / NAME).write_bytes(gzip.compress(
        (json.dumps(c3.encoded(raw), separators=(',', ':')) + '\n').encode(),
        mtime=0))
    print('ALL_CHART_PARTIAL_MEAN_COMPUTED', len(packets), flush=True)
    return owner, packets


if __name__ == '__main__':
    run()
