"""Check current predecessor wiring and all original signed integral cells."""
import contextlib
import gzip
import io
import json
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_full_signed_transport as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions_check as backend_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.encode(current.serialized(value))


def independent_current_wiring_and_measure():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    op = SimpleNamespace(flow=f, c=c, P0=[], Rm_factor=f.scalar(5))
    first = current.previous.source_module.first.FirstSwitchFunctions.__new__(current.previous.source_module.first.FirstSwitchFunctions)
    first.flow, first.c = f, c
    stale = {name: [f.scalar(-91), f.scalar(37)] for name in current.RATES}
    seed = {name: [f.scalar('.73'), f.scalar('-.04')] for name in current.RATES}
    owner = current.WholeZFullSignedTransport.__new__(current.WholeZFullSignedTransport)
    owner.c, owner.N, owner.identity, owner.transports = c, 160, {'manufactured': True}, {}
    owner.owner = lambda ends: op
    owner.bridge_transport = lambda ends: dict(actual_R100_signed_correction_C0_Z=seed)
    owner.switch_owner = lambda ends: SimpleNamespace(first=first, incoming=stale)
    owner.long_source = SimpleNamespace(T=c.mpf('1.25'))
    owner.long_owner = lambda ends: SimpleNamespace(reference=SimpleNamespace(gap=c.mpf('2.25')))
    owner.patch_source = SimpleNamespace(inlet=forbidden)
    charts = (*current.switch_backend.CHARTS, *current.previous.previous.CHARTS,
        'actual_patch', *current.OUTER_CHARTS)
    drivers = {chart: p.mpf(i+1)/1000*(-1 if i % 2 else 1) for i, chart in enumerate(charts)}
    derivatives = {chart: p.mpf(i+2)/2000*(-1 if i % 3 else 1) for i, chart in enumerate(charts)}
    def query(ends, chart, *args):
        return dict(actual_signed_nonlinear_density_C0_Z=dict(
            kernels={name: f.scalar(str(drivers[chart])) for name in current.RATES},
            Z_derivatives={name: f.scalar(str(derivatives[chart])) for name in current.RATES}),
            actual_signed_primitive_C0_Z_phi={name: f.scalar(0) for name in current.backend.OUTPUTS},
            actual_source_geometry=dict(full_period=False), actual_signed_inverse_function_alternatives=[],
            actual_signed_E_V_C0_Z=dict(E=[f.scalar(2), f.scalar('.15')]),
            actual_original_source_packet=dict(original_roots=dict(E={current.C0:f.scalar(2), current.Z:f.scalar('.15')}),
                original_generic_source=dict(common_own_five_histories_axial5={name:
                    [f.scalar('-.2'), f.scalar('.03')] for name in current.RATES})))
    owner.query = query
    h = p.exp(-30)
    lengths = dict(zip(charts, (h,h,p.log(p.mpf(11)/10)-2*h,
        p.mpf('1.25'),p.mpf('2.25'),p.mpf(1),p.mpf(1),p.mpf(1),
        p.mpf(5),p.mpf(1),p.exp(40)-1,p.mpf(11),p.mpf(1),p.mpf(2))))
    comparisons = 0
    for chart in current.OUTER_CHARTS:
        partition = current.OUTER_PARTITIONS[chart]
        for left, right in zip(partition, partition[1:]):
            got = owner.local_integral(('0','0'), chart, left, right)
            l, r = p.mpf(left[0])/left[1], p.mpf(right[0])/right[1]
            width = p.exp(40*l)*p.expm1(40*(r-l)) if chart == 'O2_axial' else r-l
            for name, rate in current.RATES.items():
                lam = p.mpf(rate.numerator)/rate.denominator
                mass = -p.expm1(-lam*width)/lam if lam else width
                _, _, bound, _, _ = owner.outer_weights(('0','0'), chart, left, right, rate)
                lo, hi = current.ep(bound.finite_interval(max_log=2000))
                assert lo <= mass <= hi
                for row, want in zip(got['actual_signed_local_integral_C0_Z'][name],
                        (drivers[chart]*mass, derivatives[chart]*mass)):
                    lo, hi = current.ep(row.finite_interval(max_log=2000))
                    assert lo <= want <= hi
                    comparisons += 1
    with contextlib.redirect_stdout(io.StringIO()):
        live = owner.full_transport(('0','0'))
    stages = live['actual_source_owned_zero_inlet_to_Rc_signed_stages']
    assert encoded(stages['switch']['actual_R100_correction_incoming_C0_Z']) == encoded(seed)
    assert encoded(stages['switch']['actual_R100_correction_incoming_C0_Z']) != encoded(stale)
    assert encoded(stages['patch']['actual_Rm_correction_incoming_C0_Z']) == encoded(stages['long']['actual_Rm_signed_correction_C0_Z'])
    assert stages['patch']['actual_Rm_incoming_binding']['correction_only_not_complete_background_history']
    total = sum(lengths.values(), p.mpf(0))
    for name, rate in current.RATES.items():
        lam = p.mpf(rate.numerator)/rate.denominator
        wants = [p.mpf('.73')*p.exp(-lam*total), -p.mpf('.04')*p.exp(-lam*total)]
        end = p.mpf(0)
        for chart in charts:
            length = lengths[chart]
            end += length
            weight = p.exp(-lam*(total-end))*(-p.expm1(-lam*length)/lam if lam else length)
            wants[0] += drivers[chart]*weight
            wants[1] += derivatives[chart]*weight
        for row, want in zip(live['actual_Rc_signed_correction_C0_Z'][name], wants):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
        for row, want in zip(live['actual_Rc_signed_complete_own_history_C0_Z'][name],
                (wants[0]-p.mpf('.2'), wants[1]+p.mpf('.03'))):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
    return dict(passed=True, independent_original_outer_cell_C0_Z_mass_comparisons=220,
        independent_fourteen_piece_current_R100_to_Rc_signed_integral_and_background_comparisons=comparisons-220,
        stale_switch_seed_sentinel_and_forbidden_legacy_patch_inlet_checked=True,
        distinct_signed_drivers_nonzero_incoming_actual_axial_Jacobian_and_pressure_memory_checked=True)


def integral_key(ends, chart, left, right):
    coordinate = lambda value: value if isinstance(value, str) else tuple(value)
    return tuple(ends), chart, coordinate(left), coordinate(right)


def check_windows(owner, ends, windows, initial, counts, integral_snapshots, patch_window=False):
    prior = initial
    for window in windows:
        chart = 'actual_patch' if patch_window else window['actual_chart']
        assert encoded(window['actual_incoming_correction_C0_Z'] if not patch_window else
            window['actual_Rm_correction_incoming_C0_Z']) == encoded(prior)
        cells = window['actual_signed_patch_cells'] if patch_window else window['actual_signed_cells']
        for cell in cells:
            left = cell['exact_left'] if cell['exact_left'] == 'Rh' else tuple(cell['exact_left'])
            right = cell['exact_right'] if cell['exact_right'] == 'Rh' else tuple(cell['exact_right'])
            # This snapshot was produced by the fresh full-transport replay,
            # never hydrated from the saved report. Verify it exactly once.
            integral = integral_snapshots.pop(integral_key(ends, chart, left, right))
            if patch_window:
                digest_function = current.patch_backend.function_digest
            elif chart in current.switch_backend.CHARTS:
                digest_function = current.switch_backend.function_digest
            elif chart in current.previous.previous.CHARTS:
                digest_function = current.previous.previous.function_digest
            else:
                digest_function = current.function_digest
            digest = digest_function(integral)
            assert digest == cell['current_signed_function_and_integral_sha256']
            backend_check.verify(owner, integral['actual_source_function'], counts)
            assert integral['actual_source_function']['actual_source_bounded_exponent_callback_binding'] == current.previous.DENSITY_BINDING
            assert integral['physical_Jacobian_applied_once']
            counts['actual_signed_integral_cells'] += 1
            counts['actual_local_integral_rows'] += 10
        memory = window['original_incoming_own_rate_memory']
        inherited = window['actual_retained_incoming_C0_Z']
        local = window['actual_signed_local_driver_at_Rh_C0_Z'] if patch_window else window['actual_signed_local_driver_C0_Z']
        outgoing = window['actual_Rh_signed_correction_C0_Z'] if patch_window else window['actual_exit_signed_correction_C0_Z']
        for name in current.RATES:
            for i in range(2):
                retained = prior[name][i]*memory[name]
                assert encoded(retained) == encoded(inherited[name][i])
                assert encoded(retained+local[name][i]) == encoded(outgoing[name][i])
                counts['actual_incoming_memory_rows'] += 1
        assert current.ep(memory['p'].finite_interval(max_log=2000)) == (1,1)
        prior = outgoing
    return prior


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Source changed: '+name
    independent = independent_current_wiring_and_measure()
    counts = dict(actual_source_function_queries=0, nonflat_directed_inverse_alternatives=0,
        flat_inverse_alternatives=0, direct_signed_primitive_rows=0, nonlinear_signed_density_rows=0,
        exact_directed_interval_rows=0, actual_signed_integral_cells=0, actual_local_integral_rows=0,
        actual_incoming_memory_rows=0, actual_current_predecessor_joins=0,
        actual_Rc_correction_background_rows=0, original_quiet_power_inherited_correction_rows=0)
    with mp.workdps(540), patch.object(current.previous.r100.WholeZR100FiniteN, 'contribution', forbidden), \
            patch.object(current.switch_backend.switch.WholeZSwitchFiniteN, 'contribution', forbidden), \
            patch.object(current.previous.previous.source_module.WholeZLongRmFiniteN, 'contribution', forbidden), \
            patch.object(current.patch_backend.patch.WholeZRmPatchFiniteN, 'contribution', forbidden), \
            patch.object(current.rh.WholeZRhReferenceFiniteN, 'contribution', forbidden), \
            patch.object(current.slope.WholeZO2SlopeFiniteN, 'contribution', forbidden), \
            patch.object(current.axial.WholeZO2AxialBufferFiniteN, 'contribution', forbidden), \
            patch.object(current.o3.WholeZO3RcFiniteN, 'contribution', forbidden):
        owner = current.WholeZFullSignedTransport()
        fresh_local_integral, integral_snapshots = owner.local_integral, {}
        def captured_local_integral(ends, chart, left, right):
            result = fresh_local_integral(ends, chart, left, right)
            key = integral_key(result['exact_Z_cell'], result['actual_chart'], result['exact_left'], result['exact_right'])
            assert key not in integral_snapshots, 'Fresh transport evaluated a cell more than once'
            integral_snapshots[key] = result
            return result
        owner.local_integral = captured_local_integral
        assert saved['source_family'] == owner.identity and saved['candidate_N'] == owner.N == 2**3981
        assert saved['original_single_incoming_callback_bindings'] == dict(switch=current.SWITCH_BINDING, patch=current.PATCH_BINDING)
        for old in saved['actual_signed_full_transports']:
            ends = tuple(old['exact_Z_cell'])
            live = owner.full_transport(ends)
            same_source(encoded(live), old, counts)
            stages = live['actual_source_owned_zero_inlet_to_Rc_signed_stages']
            bridge, switch, long, patch_stage, outer = (stages[name] for name in ('bridge','switch','long','patch','outer'))
            assert switch['actual_current_bridge_inlet_transport_sha256'] == current.function_digest(bridge)
            assert long['actual_signed_R100_R110_inlet_transport_sha256'] == current.previous.previous.function_digest(switch)
            assert patch_stage['actual_Rm_incoming_binding']['actual_current_long_inlet_transport_sha256'] == current.function_digest(long)
            assert outer['actual_current_patch_inlet_transport_sha256'] == current.function_digest(patch_stage)
            joins = ((bridge['actual_R100_signed_correction_C0_Z'], switch['actual_R100_correction_incoming_C0_Z']),
                (switch['actual_R110_signed_correction_C0_Z'], long['actual_R110_current_signed_correction_incoming_C0_Z']),
                (long['actual_Rm_signed_correction_C0_Z'], patch_stage['actual_Rm_correction_incoming_C0_Z']),
                (patch_stage['actual_Rh_signed_correction_C0_Z'], outer['actual_Rh_current_signed_correction_incoming_C0_Z']))
            for first, second in joins:
                assert encoded(first) == encoded(second)
                counts['actual_current_predecessor_joins'] += 1
            initial = bridge['actual_current_inlet_correction_C0_Z']
            assert all(row.zero for rows in initial.values() for row in rows)
            check_windows(owner, ends, bridge['actual_signed_bridge_windows'], initial, counts, integral_snapshots)
            check_windows(owner, ends, switch['actual_signed_switch_windows'], switch['actual_R100_correction_incoming_C0_Z'], counts, integral_snapshots)
            check_windows(owner, ends, long['actual_signed_long_windows'], long['actual_R110_current_signed_correction_incoming_C0_Z'], counts, integral_snapshots)
            check_windows(owner, ends, [patch_stage], patch_stage['actual_Rm_correction_incoming_C0_Z'], counts, integral_snapshots, patch_window=True)
            check_windows(owner, ends, outer['actual_signed_outer_windows'], outer['actual_Rh_current_signed_correction_incoming_C0_Z'], counts, integral_snapshots)
            assert not integral_snapshots, 'Every cell from the fresh transport must be verified exactly once'
            quiet = outer['actual_signed_outer_windows'][-1]
            assert quiet['actual_chart'] == 'O3_power'
            assert all(row.zero for rows in quiet['actual_signed_local_driver_C0_Z'].values() for row in rows)
            for name in current.RATES:
                for i in range(2):
                    assert encoded(quiet['actual_exit_signed_correction_C0_Z'][name][i]) == encoded(quiet['actual_retained_incoming_C0_Z'][name][i])
                    counts['original_quiet_power_inherited_correction_rows'] += 1
                    complete = live['actual_original_Rc_background_C0_Z'][name][i]+live['actual_Rc_signed_correction_C0_Z'][name][i]
                    assert encoded(complete) == encoded(live['actual_Rc_signed_complete_own_history_C0_Z'][name][i])
                    counts['actual_Rc_correction_background_rows'] += 1
            source = outer['actual_original_Rc_signed_source_function']
            roots = source['actual_original_source_packet']['original_roots']['E']
            assert all(encoded(row) == encoded(roots[part]) for row, part in
                zip(outer['actual_positive_Rc_normalized_amplitude_C0_Z'], (current.C0,current.Z)))
            assert current.ep(roots[current.C0].coefficient)[0] > 0
            assert all(live[key] is False for key in current.OPEN)
            assert not live['actual_full_prefix_C1_defect_functions_installed'] and not live['actual_terminal_controls_installed']
        expected = len(current.CELLS)*(10+6+8+len(current.patch_backend.PARTITION)-1+22)
        assert counts['actual_signed_integral_cells'] == expected
    receipt = dict(all_passed=True, **{current.GATE: True}, source_family=owner.identity, candidate_N=owner.N,
        current_signed_bridge_switch_long_patch_outer_predecessor_wiring_checked=True,
        all_original_signed_integral_cells_and_true_own_memories_checked=True,
        every_cell_verified_from_fresh_replay_snapshot_without_second_inverse_evaluation=True,
        independent_current_predecessor_and_original_measure=independent, replay_counts=counts,
        original_single_incoming_callback_bindings=dict(switch=current.SWITCH_BINDING, patch=current.PATCH_BINDING),
        actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
        actual_global_frequency_admitted=False, **dict.fromkeys(current.OPEN, False),
        input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
            Path(__file__).name: current.sha(Path(__file__).name),
            Path(backend_check.__file__).name: current.sha(Path(backend_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('Current whole-Z true signed zero-inlet-to-Rc integral graph checks passed', flush=True)
    return receipt


if __name__ == '__main__':
    run()
