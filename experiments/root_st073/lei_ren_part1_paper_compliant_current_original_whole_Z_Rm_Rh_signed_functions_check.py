"""Check actual patch inverse functions, dx/x integrals and incoming memory."""
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
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_Rh_signed_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions_check as previous_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.encode(current.serialized(value))


def coordinate(value):
    return 'Rh' if value == 'Rh' else tuple(value)


def independent_measure_and_memory():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    owner = current.WholeZRmRhSignedFunctions.__new__(current.WholeZRmRhSignedFunctions)
    owner.c, owner.N, owner.identity = c, 160, {'manufactured': True}
    op = SimpleNamespace(flow=f, c=c, P0=[], Rm_factor=f.scalar(5))
    owner.owner = lambda ends: op
    constant = False

    def log_value(value):
        return c.mpf(1) if value == 'Rh' else c.ln(c.mpf(value[0])/value[1])

    def query(ends, chart, left, right=None):
        right = left if right is None else right
        yl, yr = log_value(left), log_value(right)
        y = c.mpf((current.ep(yl)[0], current.ep(yr)[1]))
        density = dict(kernels={name: f.scalar(c.mpf('-.37')+(0 if constant else c.mpf('.08')*y))
            for name in current.RATES},
            Z_derivatives={name: f.scalar(c.mpf('.021')+(0 if constant else c.mpf('.006')*y))
            for name in current.RATES})
        return dict(actual_signed_nonlinear_density_C0_Z=density,
            actual_signed_primitive_C0_Z_phi={name: f.scalar(0) for name in current.previous.OUTPUTS},
            actual_signed_inverse_function_alternatives=[{'manufactured': True}],
            actual_source_geometry=dict(full_period=False),
            actual_original_source_packet=dict(original_generic_source=dict(
                common_own_five_histories_axial5={name: [f.scalar('-.2'), f.scalar('.03')]
                    for name in current.RATES})))

    owner.query = query
    comparisons = 0
    for left, right in (((1, 1), (2, 1)), ((3, 2), (5, 2)), ((1, 1), 'Rh')):
        got = owner.local_integral(('0', '0'), 'actual_patch', left, right)
        yl = p.log(p.mpf(left[0])/left[1])
        yr = p.mpf(1) if right == 'Rh' else p.log(p.mpf(right[0])/right[1])
        W = yr-yl
        lo, hi = current.ep(got['actual_physical_log_width'])
        assert lo <= W <= hi
        for name, rate in current.RATES.items():
            lam = p.mpf(rate.numerator)/rate.denominator
            mass = -p.expm1(-lam*W)/lam if lam else W
            first = (W-mass)/lam if lam else W*W/2
            targets = (-p.mpf('.37')*mass+p.mpf('.08')*(yl*mass+first),
                p.mpf('.021')*mass+p.mpf('.006')*(yl*mass+first))
            for row, target in zip(got['actual_signed_local_integral_C0_Z'][name], targets):
                lo, hi = current.ep(row.finite_interval(max_log=2000))
                assert lo <= target <= hi
                comparisons += 1
    constant = True
    incoming = {name: [f.scalar('.73'), f.scalar('-.04')] for name in current.RATES}
    owner.patch_source = SimpleNamespace(inlet=lambda ends: dict(actual_Rm_correction_C0_Z=incoming))
    with contextlib.redirect_stdout(io.StringIO()):
        result = owner.transport(('0', '0'))
    for name, rate in current.RATES.items():
        lam = p.mpf(rate.numerator)/rate.denominator
        memory = p.exp(-lam)
        mass = -p.expm1(-lam)/lam if lam else p.mpf(1)
        for row, target in zip(result['actual_Rh_signed_correction_C0_Z'][name],
                (p.mpf('.73')*memory-p.mpf('.37')*mass,
                    -p.mpf('.04')*memory+p.mpf('.021')*mass)):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= target <= hi
            comparisons += 1
        for row, target in zip(result['actual_Rh_signed_complete_own_history_C0_Z'][name],
                (-p.mpf('.2')+p.mpf('.73')*memory-p.mpf('.37')*mass,
                    p.mpf('.03')-p.mpf('.04')*memory+p.mpf('.021')*mass)):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= target <= hi
            comparisons += 1
    assert current.ep(result['original_incoming_own_rate_memory']['p'].finite_interval(max_log=2000)) == (1, 1)
    return dict(passed=True, independent_affine_log_radius_signed_C0_Z_integral_comparisons=30,
        independent_nonzero_incoming_suffix_transport_and_background_separation_comparisons=comparisons-30,
        original_log_Jacobian_and_true_length_one_checked=True, rate_zero_pressure_memory_exactly_one=True)


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Source changed: '+name
    independent = independent_measure_and_memory()
    counts = dict(actual_source_function_queries=0, nonflat_directed_inverse_alternatives=0,
        flat_inverse_alternatives=0, direct_signed_primitive_rows=0, nonlinear_signed_density_rows=0,
        exact_directed_interval_rows=0, actual_patch_integral_cells=0,
        actual_signed_local_integral_rows=0, actual_retained_incoming_rows=0,
        actual_Rh_correction_and_background_rows=0, typed_rejections=0)
    with mp.workdps(540), patch.object(current.patch.WholeZRmPatchFiniteN, 'contribution', forbidden), \
            patch.object(current.current.WholeZO3RcFiniteN, 'contribution', forbidden), \
            patch.object(current.previous.phase.first, 'NativePhaseFirstJets', forbidden), \
            patch.object(current.previous.signed, 'NativeSignedUDensityCover', forbidden):
        owner = current.WholeZRmRhSignedFunctions()
        assert saved['source_family'] == owner.identity and saved['candidate_N'] == owner.N == 2**3981
        assert encoded(current.PARTITION) == saved['exact_patch_partition']
        for old in saved['actual_signed_patch_point_functions']:
            live = owner.query(tuple(old['exact_Z_cell']), 'actual_patch', coordinate(old['exact_left']))
            same_source(encoded(live), old, counts)
            previous_check.verify(owner, live, counts)
            assert not live['actual_source_geometry']['full_period']
            for box in live['actual_source_geometry']['actual_phase_boxes']:
                lo, hi = current.ep(box)
                assert 0 <= lo <= hi <= 1 and hi-lo < mp.mpf('1e-70')
        for old in saved['actual_signed_patch_transports']:
            ends = tuple(old['exact_Z_cell'])
            live = owner.transport(ends)
            same_source(encoded(live), old, counts)
            op = owner.owner(ends)
            assert live['exact_common_P0_axial5'] is op.P0
            assert live['exact_total_log_width'] == 1
            assert live['actual_nonzero_incoming_memory_retained']
            assert live['correction_only_incoming_distinct_from_background']
            assert live['current_source_signed_inverse_used_in_every_cell']
            assert live['interval_caps_not_substituted_for_defining_functions']
            assert all(live[key] is False for key in current.OPEN)
            for cell in live['actual_signed_patch_cells']:
                integral = owner.local_integral(ends, 'actual_patch', coordinate(cell['exact_left']), coordinate(cell['exact_right']))
                assert current.function_digest(integral) == cell['current_signed_function_and_integral_sha256']
                previous_check.verify(owner, integral['actual_source_function'], counts)
                assert integral['physical_Jacobian_applied_once']
                assert integral['incoming_correction_not_supplied_or_reset']
                assert integral['live_inverse_function_extension_under_Z_independent_integral']
                assert integral['directed_rectangle_function_enclosure_not_selected_integral_value']
                counts['actual_patch_integral_cells'] += 1
                counts['actual_signed_local_integral_rows'] += 10
            for name in current.RATES:
                for i in range(2):
                    incoming = live['actual_Rm_correction_incoming_C0_Z'][name][i]
                    memory = live['original_incoming_own_rate_memory'][name]
                    retained = live['actual_retained_incoming_C0_Z'][name][i]
                    assert encoded(incoming*memory) == encoded(retained)
                    outgoing = retained + live['actual_signed_local_driver_at_Rh_C0_Z'][name][i]
                    assert encoded(outgoing) == encoded(live['actual_Rh_signed_correction_C0_Z'][name][i])
                    complete = live['actual_original_Rh_background_C0_Z'][name][i] + outgoing
                    assert encoded(complete) == encoded(live['actual_Rh_signed_complete_own_history_C0_Z'][name][i])
                    counts['actual_retained_incoming_rows'] += 1
                    counts['actual_Rh_correction_and_background_rows'] += 2
            assert live['original_incoming_own_rate_memory']['p'].scale.powers == (0, 0, 0, 0, 0)
            assert current.ep(live['original_incoming_own_rate_memory']['p'].coefficient) == (1, 1)
        assert counts['actual_patch_integral_cells'] == 4*(len(current.PARTITION)-1)
        assert counts['actual_source_function_queries'] == 12+counts['actual_patch_integral_cells']
        ends = current.CELLS[0]
        for call in (lambda: owner.query(ends, 'actual_patch', (0, 1)),
                lambda: owner.query(ends, 'actual_patch', (3, 1)),
                lambda: owner.radius_query(ends, 'actual_patch', (2, 1), (1, 1)),
                lambda: owner.local_integral(ends, 'actual_patch', (1, 1), (1, 1)),
                lambda: owner.radius_query(ends, 'actual_patch', 'Rh', (2, 1)),
                lambda: owner.radius_query(ends, 'actual_patch', (1, 0))):
            try:
                call()
            except (ValueError, TypeError):
                counts['typed_rejections'] += 1
            else:
                raise AssertionError('Invalid actual patch function request admitted')
    receipt = dict(all_passed=True, **{current.GATE: True}, source_family=owner.identity, candidate_N=owner.N,
        actual_current_patch_inverse_C0_Z_phase_functions_checked=True,
        all_original_patch_cells_signed_C1_integrals_and_true_suffix_transport_checked=True,
        nonzero_same_N_incoming_and_independent_P0_background_separation_checked=True,
        independent_actual_log_measure_and_memory=independent, replay_counts=counts,
        actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
        actual_global_frequency_admitted=False, **dict.fromkeys(current.OPEN, False),
        input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
            Path(__file__).name: current.sha(Path(__file__).name),
            Path(previous_check.__file__).name: current.sha(Path(previous_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('Current whole-Z signed Rm/Rh functions, C1 integrals and incoming-memory checks passed', flush=True)
    return receipt


if __name__ == '__main__':
    run()
