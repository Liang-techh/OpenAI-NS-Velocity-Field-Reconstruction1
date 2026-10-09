"""Check signed bridge functions, source collar and true zero-inlet transport."""
import contextlib
from fractions import Fraction
import gzip
import io
import json
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_signed_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions_check as backend_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.encode(current.serialized(value))


def independent_bounded_signed_exponent():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    prior = current.backend.signed.density.local.prior
    H = c.mpf(2)**10000
    scale = prior.FormalScale(f.logs, offset=c.mpf((-current.ep(H)[1], 0)))
    comparisons = 0
    for coefficient in (c.mpf((-.4,-.2)), c.mpf((.2,.4)), c.mpf((-.4,.3))):
        value = prior.ScaledEnclosure(scale, coefficient, f.ledger)
        before = encoded(value)
        got = current.bounded_signed_exponent(value)
        assert encoded(value) == before
        for offset in (p.mpf(0), p.mpf(-5), p.mpf(-1000)):
            for coefficient_endpoint in current.ep(coefficient):
                want = p.mpf(coefficient_endpoint)*p.exp(offset)
                assert current.ep(got)[0] <= want <= current.ep(got)[1]
                comparisons += 1
        if current.ep(coefficient)[1] < 0:
            assert current.ep(got)[1] <= 0
        if current.ep(coefficient)[0] > 0:
            assert current.ep(got)[0] >= 0
    A = prior.ScaledEnclosure(scale, c.mpf(('-.4','-.2')), f.ledger)
    primitives = {name: f.scalar(0) for name in current.backend.OUTPUTS}
    primitives.update(A=A, A_Z=f.scalar('.1'), B_over_Pstar=f.scalar('-.3'), B_Z_over_Pstar=f.scalar('.02'))
    result = current.SIGNED_DENSITY(f.scalar(2), f.scalar('.3'), f.scalar('-.2'), f.scalar('.04'), primitives, 160)
    for logscale in (p.mpf(0), p.mpf(-5), p.mpf(-1000)):
        a = -p.mpf('.3')*p.exp(logscale)/160
        de, dez = 2*p.expm1(a), p.mpf('.3')*p.expm1(a)+2*p.exp(a)*p.mpf('.1')/160
        dv, dvz = -p.mpf('.3')/160, p.mpf('.02')/160
        crossz = p.mpf('.3')*de+2*dez
        squarez = de*dez
        kernels = dict(m=dv,h=de,k=-p.mpf('.2')*de+2*dv+de*dv,
            e=-p.mpf('.4')*dv+dv*dv-2*de-de*de/2,p=2*de+de*de/2)
        derivatives = dict(m=dvz,h=dez,
            k=p.mpf('.04')*de-p.mpf('.2')*dez+p.mpf('.3')*dv+2*dvz+dez*dv+de*dvz,
            e=p.mpf('.08')*dv-p.mpf('.4')*dvz+2*dv*dvz-crossz-squarez,p=crossz+squarez)
        for part, wants in (('kernels', kernels), ('Z_derivatives', derivatives)):
            for name, want in wants.items():
                got = current.bounded_signed_exponent(result[part][name])
                assert current.ep(got)[0] <= want <= current.ep(got)[1]
                comparisons += 1
    module = current.backend.current.reference.phase.densities
    assert current.SIGNED_DENSITY.__code__ is module.density_Z_kernels.__code__
    assert current.FACTORED_INCREMENT.__code__ is module.density.factored_expm1.__code__
    try:
        current.bounded_signed_exponent(f.scalar(2))
    except ArithmeticError:
        pass
    else:
        raise AssertionError('Unbounded exponent admitted')
    return dict(passed=True, independent_signed_tail_and_nonlinear_C0_Z_comparisons=comparisons,
        actual_formal_row_unchanged_and_signed_tail_not_zeroed=True,
        original_density_and_expm1_code_bodies_identical=True, unbounded_input_rejected=True)


def independent_phase():
    c = MPIntervalContext()
    c.dps = 160
    p = mp.mp.clone()
    p.dps = 1600
    N = 2**3981
    sc = c.mpf(1)/8
    P = c.exp(40)+11
    flow = SimpleNamespace(logs=(c.mpf(-10000), 2*P, c.mpf(0), c.ln(4)-1000-4*P, c.mpf(0)), h=c.mpf(1))
    flow.Y = c.ln(100)-flow.logs[3]
    op = SimpleNamespace(flow=flow, c=c, Rm_factor=object())
    stub = SimpleNamespace(sc=sc, parameter_binding={'manufactured_source_binding': True})
    with patch.object(current.backend.phase, 'RmRadiusPhase', lambda *args: stub):
        phase = current.CurrentBridgeRadiusPhase(op, {'actual_five_defect_family_sha256': 'fixture'}, N, sc)
    Y, h = p.log(25)+1000+4*(p.exp(40)+11), p.exp(-10000)
    cases = (('first_micro', current.INLET), ('first_micro', (1,2)), ('first_micro', (1,1)),
        ('second_micro', (1,1)), ('second_micro', (3,2)), ('second_micro', (2,1)),
        ('frozen_macro', (0,1)), ('frozen_macro', (1,3)), ('frozen_macro', (1,1)))
    comparisons = 0
    for chart, coordinate in cases:
        q = p.mpf(1)/16 if coordinate == current.INLET else p.mpf(coordinate[0])/coordinate[1]
        value = N*(Y*q+h*(2*(1-q)-p.mpf(1)/16)) if chart == 'frozen_macro' else N*h*(q-p.mpf(1)/16)
        want = value-p.floor(value)
        got = phase.query(chart, coordinate)
        assert not got['full_period']
        assert any(current.ep(box)[0] <= want <= current.ep(box)[1] for box in got['actual_phase_boxes'])
        assert got['actual_positive_width_retained_not_zeroed']
        comparisons += 1
    for a, qa, b, qb in (('first_micro', (1,1), 'second_micro', (1,1)),
            ('second_micro', (2,1), 'frozen_macro', (0,1))):
        assert encoded(phase.query(a, qa)['actual_phase_boxes']) == encoded(phase.query(b, qb)['actual_phase_boxes'])
    inlet = phase.query('first_micro', current.INLET)
    assert inlet['exact_inlet_phase_zero'] and current.ep(inlet['actual_phase_boxes'][0]) == (0,0)
    assert current.ep(inlet['signed_actual_positive_bridge_width_error']) == (0,0)
    for left, right in (((0,1), (1,2)), ((1,2), (1,1))):
        got = phase.query('frozen_macro', left, right)
        assert got['full_period'] and got['full_period_proved_from_actual_N_physical_width']
    q = Fraction(1,3)
    r = q+Fraction(1, 2**4050)
    narrow = phase.query('frozen_macro', current.rational(q), current.rational(r))
    assert not narrow['full_period']
    for coordinate in (q, r):
        t = p.mpf(coordinate.numerator)/coordinate.denominator
        value = N*(Y*t+h*(2*(1-t)-p.mpf(1)/16))
        want = value-p.floor(value)
        assert any(current.ep(box)[0] <= want <= current.ep(box)[1] for box in narrow['actual_phase_boxes'])
        comparisons += 1
    # The selected production origin has an astronomical exponent: retain
    # its MPF tuple and symbolic coordinate rather than allocating 2**(-e).
    tiny = c.make_mpf(((0,1,-10**12,1), (0,1,-10**12,1)))
    tiny_stub = SimpleNamespace(sc=tiny, parameter_binding={})
    with patch.object(current.backend.phase, 'RmRadiusPhase', lambda *args: tiny_stub):
        tiny_phase = current.CurrentBridgeRadiusPhase(op, {'actual_five_defect_family_sha256': 'fixture'}, N, tiny)
    point = tiny_phase.query('first_micro', current.INLET)
    assert point['selected_source_sc']._mpi_ == tiny._mpi_
    assert point['exact_inlet_phase_zero'] and current.ep(point['actual_phase_boxes'][0]) == (0,0)
    return dict(passed=True, independent_actual_radius_phase_comparisons=comparisons,
        exact_micro_micro_and_micro_macro_joins=2, actual_macro_full_period_images=2,
        exact_zero_inlet_and_astronomical_sc_symbolic_coordinate_checked=True,
        positive_micro_width_phase_not_zeroed=True, narrow_macro_cell_retained=True)


def independent_measure_and_transport():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    series = current.source_module.first.FirstSwitchFunctions.__new__(current.source_module.first.FirstSwitchFunctions)
    series.flow, series.c = f, c
    op = SimpleNamespace(flow=f, c=c, P0=[], Rm_factor=f.scalar(5))
    owner = current.WholeZBridgeSignedFunctions.__new__(current.WholeZBridgeSignedFunctions)
    owner.c, owner.N, owner.identity, owner.sc = c, 160, {'manufactured': True}, c.mpf(1)/8
    owner.start, owner.inlet_proof = owner.sc/2, {'source_owned': True}
    owner.owner = lambda ends: op
    owner.bridge_owner = lambda ends: SimpleNamespace(series=series)
    owner.bridge_background = SimpleNamespace(seam=lambda ends: {'same_source': True})
    drivers = dict(zip(current.CHARTS, ('-.37', '.13', '-.21')))
    derivatives = dict(zip(current.CHARTS, ('.021', '-.015', '.04')))
    def query(ends, chart, left, right=None):
        flat = chart == 'first_micro' and left == 'inlet'
        return dict(actual_signed_nonlinear_density_C0_Z=dict(
            kernels={name: f.scalar(0 if flat else drivers[chart]) for name in current.RATES},
            Z_derivatives={name: f.scalar(0 if flat else derivatives[chart]) for name in current.RATES}),
            actual_signed_primitive_C0_Z_phi={name: f.scalar(0) for name in current.backend.OUTPUTS},
            actual_source_geometry=dict(full_period=False), actual_original_source_packet=dict(
                source_owned_collar_q_C0_Z_exact_zero=flat, original_generic_source=dict(
                    common_own_five_histories_axial5={name: [f.scalar('-.2'), f.scalar('.03')] for name in current.RATES})))
    owner.query = query
    comparisons = 0
    lengths = dict(zip(current.CHARTS, (p.exp(-30)*p.mpf(15)/16, p.exp(-30), p.log(20)-2*p.exp(-30))))
    for chart, left, right in (('first_micro', current.INLET, (3,4)),
            ('second_micro', (5,4), (7,4)), ('frozen_macro', (1,4), (3,4))):
        got = owner.local_integral(('0','0'), chart, left, right)
        W = p.exp(-30)*p.mpf(11)/16 if chart == 'first_micro' else p.exp(-30)/2 \
            if chart == 'second_micro' else lengths[chart]/2
        lo, hi = current.ep(got['actual_physical_log_width'].finite_interval(max_log=2000))
        assert 0 < lo <= W <= hi
        for name, rate in current.RATES.items():
            lam = p.mpf(rate.numerator)/rate.denominator
            mass = -p.expm1(-lam*W)/lam if lam else W
            for row, want in zip(got['actual_signed_local_integral_C0_Z'][name],
                    (p.mpf(drivers[chart])*mass, p.mpf(derivatives[chart])*mass)):
                lo, hi = current.ep(row.finite_interval(max_log=2000))
                assert lo <= want <= hi
                comparisons += 1
    with contextlib.redirect_stdout(io.StringIO()):
        got = owner.bridge_transport(('0','0'))
    total = sum(lengths.values(), p.mpf(0))
    assert abs(total-(p.log(20)-p.exp(-30)/16)) < p.mpf('1e-230')
    for name, rate in current.RATES.items():
        lam = p.mpf(rate.numerator)/rate.denominator
        wants = [p.mpf(0), p.mpf(0)]
        end = p.mpf(0)
        for chart in current.CHARTS:
            length = lengths[chart]
            end += length
            weight = p.exp(-lam*(total-end))*(-p.expm1(-lam*length)/lam if lam else length)
            wants[0] += p.mpf(drivers[chart])*weight
            wants[1] += p.mpf(derivatives[chart])*weight
        for row, want in zip(got['actual_R100_signed_correction_C0_Z'][name], wants):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
        for row, want in zip(got['actual_R100_signed_complete_own_history_C0_Z'][name],
                (wants[0]-p.mpf('.2'), wants[1]+p.mpf('.03'))):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
    return dict(passed=True, independent_original_micro_macro_C0_Z_kernel_mass_comparisons=30,
        independent_piecewise_signed_global_integral_and_background_comparisons=comparisons-30,
        true_total_length_Y_minus_hb_sc_half_and_zero_inlet_checked=True,
        physical_measures_suffixes_and_memories_applied_once=True)


def verify(owner, live, counts):
    backend_check.verify(owner, live, counts)
    assert live['actual_source_bounded_exponent_callback_binding'] == current.DENSITY_BINDING
    packet, geometry = live['actual_original_source_packet'], live['actual_source_geometry']
    assert packet['identical_source_P0_tuples_rebound_to_live_descendant_object']
    assert packet['original_generic_roots_preserved_without_all_u_support_evaluation']
    assert geometry['actual_positive_width_retained_not_zeroed']
    upper = geometry['actual_N_bridge_width_phase_log_upper']
    assert upper is None or current.ep(upper)[1] < -100*mp.log(10)
    for before, after in zip(packet['original_bridge_independent_P0_axial5'], live['exact_common_P0_axial5']):
        assert encoded(before) == encoded(after)
        assert before.ctx is after.ctx and before.scale.bases is after.scale.bases and before.ledger is after.ledger
    if packet['source_owned_collar_q_C0_Z_exact_zero']:
        assert all(row.zero for row in packet['original_q_C0_Z'].values())
        assert all(row.zero for row in live['actual_signed_primitive_C0_Z_phi'].values())
        assert all(row['original_inverse_and_Z_function_proof']['geometry'] == 'flat'
            for row in live['actual_signed_inverse_function_alternatives'])
        counts['source_defined_exact_flat_collar_queries'] += 1


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Source changed: '+name
    exponent, phase, measure = independent_bounded_signed_exponent(), independent_phase(), independent_measure_and_transport()
    counts = dict(actual_source_function_queries=0, nonflat_directed_inverse_alternatives=0,
        flat_inverse_alternatives=0, direct_signed_primitive_rows=0, nonlinear_signed_density_rows=0,
        exact_directed_interval_rows=0, actual_signed_bridge_integral_cells=0,
        actual_micro_narrow_phase_cells=0, actual_macro_full_period_cells=0, actual_local_integral_rows=0,
        actual_incoming_memory_rows=0, actual_R100_correction_background_rows=0,
        source_defined_exact_flat_collar_queries=0, exact_phase_joins=0, typed_rejections=0)
    with mp.workdps(540), patch.object(current.r100.WholeZR100FiniteN, 'prepare', forbidden), \
            patch.object(current.r100.WholeZR100FiniteN, 'select_N', forbidden), \
            patch.object(current.r100.WholeZR100FiniteN, 'contribution', forbidden), \
            patch.object(current.original, 'OriginalMicroFiniteN', forbidden), \
            patch.object(current.macro, 'OriginalMacroFiniteN', forbidden), \
            patch.object(current.backend.phase.first, 'NativePhaseFirstJets', forbidden), \
            patch.object(current.backend.signed, 'NativeSignedUDensityCover', forbidden):
        owner = current.WholeZBridgeSignedFunctions()
        assert saved['source_family'] == owner.identity and saved['candidate_N'] == owner.N == 2**3981
        for old in saved['actual_signed_bridge_point_functions']:
            live = owner.query(tuple(old['exact_Z_cell']), old['actual_chart'], tuple(old['exact_left']))
            same_source(encoded(live), old, counts)
            verify(owner, live, counts)
            assert not live['actual_source_geometry']['full_period']
        for old in saved['actual_signed_bridge_transports']:
            ends = tuple(old['exact_Z_cell'])
            live = owner.bridge_transport(ends)
            same_source(encoded(live), old, counts)
            assert live['zero_inlet_follows_from_source_definition_not_arbitrary_incoming_reset']
            assert live['original_background_nonzero_histories_remain_separate']
            assert all(live[key] is False for key in current.OPEN)
            inlet = live['actual_source_owned_zero_inlet_signed_function']
            verify(owner, inlet, counts)
            assert any(not row.zero for rows in inlet['actual_original_source_packet']
                ['original_generic_source']['common_own_five_histories_axial5'].values() for row in rows[:2])
            assert any(not row.zero for rows in live['actual_original_R100_background_C0_Z'].values() for row in rows)
            assert inlet['actual_source_geometry']['exact_inlet_phase_zero']
            assert current.ep(inlet['actual_source_geometry']['actual_phase_boxes'][0]) == (0,0)
            prior = live['actual_current_inlet_correction_C0_Z']
            assert all(row.zero for pair in prior.values() for row in pair)
            for window in live['actual_signed_bridge_windows']:
                chart = window['actual_chart']
                assert encoded(window['actual_incoming_correction_C0_Z']) == encoded(prior)
                for cell in window['actual_signed_cells']:
                    integral = owner.local_integral(ends, chart, tuple(cell['exact_left']), tuple(cell['exact_right']))
                    assert current.function_digest(integral) == cell['current_signed_function_and_integral_sha256']
                    verify(owner, integral['actual_source_function'], counts)
                    assert integral['physical_Jacobian_applied_once'] and integral['original_hb_hb_Y_minus_2hb_measure_and_cutoff_mass_preserved']
                    geometry = integral['actual_source_function']['actual_source_geometry']
                    if chart == 'frozen_macro':
                        assert geometry['full_period'] and geometry['full_period_proved_from_actual_N_physical_width']
                        counts['actual_macro_full_period_cells'] += 1
                    else:
                        assert not geometry['full_period']
                        counts['actual_micro_narrow_phase_cells'] += 1
                    counts['actual_signed_bridge_integral_cells'] += 1
                    counts['actual_local_integral_rows'] += 10
                for name in current.RATES:
                    for i in range(2):
                        retained = prior[name][i]*window['original_incoming_own_rate_memory'][name]
                        assert encoded(retained) == encoded(window['actual_retained_incoming_C0_Z'][name][i])
                        outgoing = retained+window['actual_signed_local_driver_C0_Z'][name][i]
                        assert encoded(outgoing) == encoded(window['actual_exit_signed_correction_C0_Z'][name][i])
                        counts['actual_incoming_memory_rows'] += 1
                assert current.ep(window['original_incoming_own_rate_memory']['p'].finite_interval(max_log=2000)) == (1,1)
                prior = window['actual_exit_signed_correction_C0_Z']
            assert encoded(prior) == encoded(live['actual_R100_signed_correction_C0_Z'])
            for name in current.RATES:
                for i in range(2):
                    complete = live['actual_original_R100_background_C0_Z'][name][i]+prior[name][i]
                    assert encoded(complete) == encoded(live['actual_R100_signed_complete_own_history_C0_Z'][name][i])
                    counts['actual_R100_correction_background_rows'] += 1
            assert live['exact_micro_macro_source_join']['actual_micro_histories_are_the_live_macro_incoming']
            for a, qa, b, qb in (('first_micro', (1,1), 'second_micro', (1,1)),
                    ('second_micro', (2,1), 'frozen_macro', (0,1)),
                    ('frozen_macro', (1,1), 'first_switch', (0,1))):
                first, second = owner.radius_query(ends, a, qa), owner.radius_query(ends, b, qb)
                assert encoded(first['actual_phase_boxes']) == encoded(second['actual_phase_boxes'])
                counts['exact_phase_joins'] += 1
        assert counts['actual_source_function_queries'] == 56 and counts['actual_signed_bridge_integral_cells'] == 40
        ends = current.CELLS[0]
        for call in (lambda: owner.radius_query(ends, 'first_micro', (0,1)),
                lambda: owner.radius_query(ends, 'second_micro', (0,1)),
                lambda: owner.radius_query(ends, 'frozen_macro', (2,1)),
                lambda: owner.radius_query(ends, 'first_micro', (1,1), (1,2)),
                lambda: owner.local_integral(ends, 'first_micro', 'inlet', 'inlet'),
                lambda: owner.local_integral(ends, 'frozen_macro', (1,1), (1,1)),
                lambda: owner.source_query(ends, 'first_micro', (1,3))):
            try:
                call()
            except ValueError:
                counts['typed_rejections'] += 1
            else:
                raise AssertionError('Invalid original bridge source/measure request admitted')
    receipt = dict(all_passed=True, **{current.GATE: True}, source_family=owner.identity, candidate_N=owner.N,
        all_three_original_bridge_signed_inverse_C1_source_and_integral_functions_checked=True,
        exact_source_defined_zero_correction_inlet_and_nonzero_background_separation_checked=True,
        original_positive_hb_hb_Y_minus2hb_widths_cutoff_masses_and_memories_checked=True,
        exact_source_sc_symbolic_inlet_P0_tuple_context_basis_ledger_binding_checked=True,
        independent_actual_bridge_radius_phase=phase, independent_actual_bridge_measure_and_transport=measure,
        independent_original_bounded_signed_exponent=exponent,
        replay_counts=counts, actual_full_prefix_C1_defect_functions_installed=False,
        actual_terminal_controls_installed=False, actual_global_frequency_admitted=False,
        **dict.fromkeys(current.OPEN, False),
        input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
            Path(__file__).name: current.sha(Path(__file__).name),
            Path(backend_check.__file__).name: current.sha(Path(backend_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('Current whole-Z signed bridge functions, source collar and true R100 C1 transport checks passed', flush=True)
    return receipt


if __name__ == '__main__':
    run()
