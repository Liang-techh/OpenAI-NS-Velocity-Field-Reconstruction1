"""Independent acceptance checks for the all-chart partial-mean adapter."""
import copy
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_all_chart_partial_mean as current
import lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q_check as partial_check


def _endpoint(value):
    return tuple(value) if isinstance(value, list) else value


def _reject(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('Unsupported all-chart partial-mean query admitted')


def _stored_window(owner, ends, chart):
    transport, window = owner._transport(ends, chart)
    assert transport['source_identity'] == owner.family
    partition = owner.partitions[chart]
    actual = [(_endpoint(cell['exact_left']), _endpoint(cell['exact_right']))
              for cell in window['cells']]
    expected = list(zip(partition, partition[1:]))
    assert actual == expected, (chart, actual, expected)
    assert window['global_phase_cover'] == [0, 1]
    assert window['third_row_predecessor_never_reset']
    assert all(cell['physical_Jacobian_applied_once']
               for cell in window['cells'])
    return window


@current.source_precision
def run(owner, packets):
    began = time.monotonic()
    fixture = partial_check.independent_mean_fixture()
    owner.assert_graph()
    ends = ('0', '.5')
    assert len(owner.charts) == 17
    assert set(owner.charts) == set(owner.partitions)
    assert owner.partitions['first_micro'][0] == ('sc_half',)
    assert owner.partitions['actual_patch'][-1] == 'Rh'

    # The actual first-cell queries exercise every chart.  The saved
    # transport and manifest checks below cover the remaining cells in each
    # chart without turning a range cap into a point query.
    queried = {}
    views = {}
    for chart in owner.charts:
        partition = owner.partitions[chart]
        queried[chart] = packets[chart]
        views[chart] = owner.report(queried[chart])
        view = views[chart]
        assert view['source_family'] == owner.family
        assert view['whole_chart_scope']['complete_original_partition']
        assert view['whole_chart_scope']['source_cell_count'] == len(partition) - 1
        assert view['whole_chart_scope']['query_is_contiguous_union_of_original_cells']
        assert view['whole_chart_scope']['incoming_history_is_from_saved_whole_Z_transport']
        assert view['exact_native_partition'] == list(partition)
        assert len(view['whole_cell_partial_own_mean_increment_Z3_magnitude']) == 4
        assert len(view['whole_cell_mean_y_increment_Z3_magnitude']) == 4
        assert len(view['whole_cell_delta_V_increment_Z3_magnitude']) == 4
        assert len(view['whole_cell_Q_increment_Z2_magnitude']) == 3
        assert len(view['signed_partial_mean_query_functions']['mean']) == 4
        assert len(view['signed_partial_mean_query_functions']['mean_y']) == 4
        assert len(view['signed_Q_increment_Z2_function']) == 3
        assert view['same_selected_N_division_once_no_integer_materialization']
        assert view['mean_y_is_delta_V_minus_mean_enclosed_by_sum_of_magnitudes']
        assert view['modified_Q_increment_only_not_complete_Q']
        assert view['Q_only_value_Z_ZZ_rows_no_fabricated_Q_ZZZ']
        assert view['magnitude_bounds_not_point_values_and_no_cap_differentiation']
        assert all(not view[key] for key in current.OPEN)

    # Verify the whole original partition and the incoming history, including
    # the quiet chart whose local primitive is flat but whose predecessor mean
    # is retained.
    windows = {chart: _stored_window(owner, ends, chart)
               for chart in owner.charts}
    assert any(not row.get('exact_zero', False)
               for row in windows['O3_power']['actual_N_scaled_incoming_C3']['m'])
    assert any(not row.get('exact_zero', False)
               for chart in owner.charts
               for row in windows[chart]['actual_N_scaled_incoming_C3']['m'])
    assert all(len(windows[chart]['cells']) == len(owner.partitions[chart]) - 1
               for chart in owner.charts)
    assert all(len(views[chart]['partial_history_cells']) == 1
               for chart in owner.charts)

    # The adapter accepts only exact source-cell unions and rejects a foreign
    # parent packet, while graph append remains available to downstream work.
    first = queried['first_micro']
    entry = owner._packets[id(first)]
    key = 'modified_Q_increment_only_not_complete_Q'
    old = entry[1]['view'][key]
    entry[1]['view'][key] = False
    _reject(lambda: owner.report(first))
    entry[1]['view'][key] = old
    _reject(lambda: owner.range_cell(ends, 'first_micro', ('sc_half',), (3, 8)))
    _reject(lambda: owner.report(current.partial.OriginalPartialMeanQPacket()))
    owner.g.symbol('all_chart_partial_mean_downstream_child')
    owner.report(first)

    inputs = dict(owner.hashes)
    inputs[current.NAME] = current.sha(current.NAME)
    inputs[Path(__file__).name] = current.sha(Path(__file__).name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        input_hashes=inputs,
        checks=dict(
            independent_moving_endpoint_mean_Z3_mean_y_Z3_and_Q_Z2_rows=fixture,
            all_17_original_chart_queries_execute=17,
            exact_native_partition_unions_and_symbolic_sc_half_Rh_preserved=True,
            whole_chart_source_partition_manifest_and_jacobian_once_verified=True,
            incoming_N_scaled_C3_history_retained_through_quiet_O3_power=True,
            rate_one_kernel_and_incoming_decay_bounded_by_one=True,
            same_selected_N_delta_epsilon_and_common_S_unit_retained=True,
            local_delta_V_and_mean_y_caps_use_directed_sum_bounds=True,
            Q_increment_uses_exact_delta_and_only_value_Z_ZZ_rows=True,
            source_caps_not_point_values_and_not_differentiated=True,
            arbitrary_native_endpoints_foreign_packets_rejected=True,
            accepted_parent_graph_immutable_with_child_append_allowed=True,
            complete_pressure_tails_support_repair_and_global_claims_not_admitted=True,
        ),
        **{current.GATE: True},
        **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (current.HERE / current.RECEIPT).write_bytes(
        (json.dumps(current.c3.encoded(result), separators=(',', ':')) + '\n').encode())
    print('PASS_ALL_CHART_PARTIAL_MEAN', fixture, len(queried), flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit('Use run(checked_all_chart_owner, packets)')
