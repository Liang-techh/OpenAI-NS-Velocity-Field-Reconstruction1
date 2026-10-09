"""Check real switch phase/widths, inverse functions and C1 incoming transport."""
import contextlib
import gzip
import io
import json
from pathlib import Path
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_switch_signed_functions as current
import lei_ren_part1_paper_compliant_current_original_whole_Z_signed_loop_functions_check as backend_check
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow


def encoded(value):
    return current.encode(current.serialized(value))


def independent_phase():
    c = MPIntervalContext()
    c.dps = 220
    p = mp.mp.clone()
    p.dps = 360
    flow = MacroFlow(c, c.mpf(-600), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    op = SimpleNamespace(flow=flow, c=c, Rm_factor=flow.scalar(5))
    mapper = current.CurrentSwitchRadiusPhase.__new__(current.CurrentSwitchRadiusPhase)
    mapper.op, mapper.identity, mapper.N = op, {'manufactured': True}, 257
    mapper.base = SimpleNamespace(sc=c.mpf('.25'), parameter_binding={'manufactured': True})
    mapper.cache = {}
    comparisons = 0
    for chart in current.CHARTS:
        for coordinate in ((0,1), (1,3), (1,1)):
            got = mapper.query(chart, coordinate)
            t = p.mpf(coordinate[0])/coordinate[1]
            regular = p.log(100)-p.log(4)+1000+4*(p.exp(40)+11)
            micro = t-p.mpf('.125') if chart == 'first_switch' else 1+t-p.mpf('.125')
            if chart == 'post_power':
                regular += p.log(p.mpf(11)/10)*t
                micro = 2*(1-t)-p.mpf('.125')
            wanted = mapper.N*(regular+p.exp(-600)*micro)
            wanted -= p.floor(wanted)
            assert any(current.ep(box)[0] <= wanted <= current.ep(box)[1] for box in got['actual_phase_boxes'])
            assert got['actual_positive_width_retained_not_zeroed'] and got['actual_phase_Z_exact_zero']
            assert not got['full_period']
            comparisons += 1
    for first, second in (('first_switch','second_switch'), ('second_switch','post_power')):
        a, b = mapper.query(first, (1,1)), mapper.query(second, (0,1))
        assert encoded(a['actual_phase_boxes']) == encoded(b['actual_phase_boxes'])
        assert encoded(a['original_microscopic_width_coefficient']) == encoded(b['original_microscopic_width_coefficient'])
    return dict(passed=True, independent_actual_radius_phase_scalar_comparisons=comparisons,
        exact_original_first_second_post_radius_join_identities=2,
        second_switch_uses_one_plus_local_coordinate=True,
        positive_microscopic_phase_term_is_bounded_not_set_zero=True)


def independent_measure_and_transport():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    first = current.original.first.FirstSwitchFunctions.__new__(current.original.first.FirstSwitchFunctions)
    first.flow, first.c = f, c
    incoming = {name: [f.scalar('.73'), f.scalar('-.04')] for name in current.RATES}
    source = SimpleNamespace(first=first, incoming=incoming)
    op = SimpleNamespace(flow=f, c=c, P0=[], Rm_factor=f.scalar(5))
    owner = current.WholeZSwitchSignedFunctions.__new__(current.WholeZSwitchSignedFunctions)
    owner.c, owner.N, owner.identity = c, 160, {'manufactured': True}
    owner.owner = lambda ends: op
    owner.switch_owner = lambda ends: source
    density = dict(kernels={name: f.scalar('-.37') for name in current.RATES},
        Z_derivatives={name: f.scalar('.021') for name in current.RATES})
    owner.query = lambda *args: dict(actual_signed_nonlinear_density_C0_Z=density,
        actual_signed_primitive_C0_Z_phi={name: f.scalar(0) for name in current.backend.OUTPUTS},
        actual_source_geometry=dict(full_period=False),
        actual_original_source_packet=dict(original_generic_source=dict(
            common_own_five_histories_axial5={name: [f.scalar('-.2'), f.scalar('.03')]
                for name in current.RATES})))
    comparisons = 0
    Y, h = p.log(p.mpf(11)/10), p.exp(-30)
    for chart in current.CHARTS:
        length = h if chart != 'post_power' else Y-2*h
        got = owner.local_integral(('0','0'), chart, (1,4), (3,4))
        W = length/2
        lo, hi = current.ep(got['actual_physical_log_width'].finite_interval(max_log=2000))
        assert 0 < lo <= W <= hi
        for name, rate in current.RATES.items():
            lam = p.mpf(rate.numerator)/rate.denominator
            mass = -p.expm1(-lam*W)/lam if lam else W
            for row, target in zip(got['actual_signed_local_integral_C0_Z'][name],
                    (-p.mpf('.37')*mass, p.mpf('.021')*mass)):
                lo, hi = current.ep(row.finite_interval(max_log=2000))
                assert lo <= target <= hi
                comparisons += 1
    with contextlib.redirect_stdout(io.StringIO()):
        got = owner.switch_transport(('0','0'))
    for name, rate in current.RATES.items():
        lam = p.mpf(rate.numerator)/rate.denominator
        memory = p.exp(-lam*Y)
        mass = -p.expm1(-lam*Y)/lam if lam else Y
        want = (p.mpf('.73')*memory-p.mpf('.37')*mass,
            -p.mpf('.04')*memory+p.mpf('.021')*mass)
        for row, target in zip(got['actual_R110_signed_correction_C0_Z'][name], want):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= target <= hi
            comparisons += 1
        for row, target in zip(got['actual_R110_signed_complete_own_history_C0_Z'][name],
                (-p.mpf('.2')+want[0], p.mpf('.03')+want[1])):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= target <= hi
            comparisons += 1
    return dict(passed=True, independent_positive_micro_post_C0_Z_mass_comparisons=30,
        independent_three_window_nonzero_incoming_and_background_comparisons=comparisons-30,
        exact_total_length_h_plus_h_plus_Y_minus2h_is_Y=True,
        original_Jacobians_applied_once_and_positive_micro_width_retained=True)


def independent_inverse_identity():
    a, t, q, h, psi, angle, r, s, sine = sy.symbols('a t q h psi angle r s sine')
    nu = 1+t*t+2*q*q
    F = ((1+t*t)*psi+2*t*q*h*(angle-psi)/r
        +q*q*((2-3*s)*angle+s*psi+r*sine/sy.pi)/(r*r))/nu
    A = a*(2*t*q*h*(angle-psi)/r
        +q*q*((2-3*s)*(angle-psi)+r*sine/sy.pi)/(r*r))/(2*nu)
    assert sy.simplify((A-a*(F-psi)/2).subs(s, 1-r*r)) == 0
    assert current.switch.primitives.THEOREM['exact_inverse_identity'] == 'A=a/2*(phi-psi_fraction)'
    def tested(source, qr, dstar, phi):
        _, branches, _ = current.backend.signed.signed_u_branches(source, dstar)
        rows = [current.inverse_identity_branch(source, qr, dstar, phi, branch) for branch in branches]
        assert all(row['values'] is not None for row in rows)
        return dict(values={name: current.directed_dominating_union([row['values'][name] for row in rows])
            for name in current.backend.OUTPUTS}, record=dict(status='enclosed', geometry=rows[0]['record']['geometry']))
    scalar = backend_check.scalar_check
    proxy = SimpleNamespace(**{**vars(scalar.current), 'Z_FIRST': tested})
    with patch.object(scalar, 'current', proxy), mp.workdps(170):
        cases = scalar.scalar_fixtures()
    return dict(passed=True, independent_original_inverse_Z_phi_and_nonlinear_density_cases=cases,
        exact_original_Mobius_and_bound_small_r_A_inverse_identities_checked=True,
        actual_selected_inverse_identity_backend_used=True,
        original_Z_phi_derivatives_retained_not_zeroed=True)


def independent_formal_union():
    c = MPIntervalContext()
    c.dps = 180
    prior = current.backend.signed.density.local.prior
    H = c.mpf(2)**10000
    bases = (c.mpf((-current.ep(H)[1], current.ep(H)[1])),) + tuple(c.mpf(0) for _ in range(4))
    ledger = {'directed_independent_log_rescalings': 0, 'directed_small_exponential_tails': 0}
    rows = [prior.ScaledEnclosure(prior.FormalScale(bases, (1,0,0,0,0), H), c.mpf((-3,-1)), ledger),
        prior.ScaledEnclosure(prior.FormalScale(bases, (2,0,0,0,0), c.mpf((-current.ep(H)[1], current.ep(H)[1]))),
            c.mpf((1,4)), ledger)]
    got = current.directed_dominating_union(rows)
    assert all(current.ep((row.scale-got.scale).evaluate())[1] <= 0 for row in rows)
    assert current.ep(got.coefficient)[0] < 0 < current.ep(got.coefficient)[1]
    assert ledger['directed_small_exponential_tails'] > 0
    assert got.scale.powers == (0,0,0,0,0)
    finite = [prior.ScaledEnclosure(prior.FormalScale(bases), c.mpf((-3,-2)), ledger),
        prior.ScaledEnclosure(prior.FormalScale(bases, offset=c.mpf('.7')), c.mpf((1,4)), ledger)]
    hull = current.directed_dominating_union(finite).finite_interval(max_log=2000)
    p = mp.mp.clone()
    p.dps = 240
    lo, hi = current.ep(hull)
    for value in (-p.mpf(3), -p.mpf(2), p.exp(p.mpf('.7')), 4*p.exp(p.mpf('.7'))):
        assert lo <= value <= hi
    return dict(passed=True, signed_finite_union_scalar_comparisons=4,
        independent_wide_formal_logs_have_checked_nonpositive_relative_exponents=True,
        enormous_physical_exponentials_not_materialized=True,
        nonzero_directed_tail_and_signed_union_retained=True,
        dominating_numeric_log_is_arithmetic_coordinate_not_source_function_value=True)


def verify(owner, live, counts):
    backend_check.verify(owner, live, counts)
    geometry, packet = live['actual_source_geometry'], live['actual_original_source_packet']
    assert geometry['actual_positive_width_retained_not_zeroed']
    assert packet['identical_source_P0_tuples_rebound_to_live_descendant_object']
    assert packet['empty_partial_cutoff_integral_exact_zero_nonempty_recipe_unchanged']
    assert geometry['actual_source_logRa_and_analytic_logP_collection_bound_by_original_Rm_binder']
    for before, after in zip(packet['original_switch_independent_P0_axial5'], live['exact_common_P0_axial5']):
        assert encoded(before) == encoded(after)
        assert before.ctx is after.ctx and before.scale.bases is after.scale.bases and before.ledger is after.ledger
    chart = live['actual_chart']
    if chart != 'post_power':
        assert not geometry['full_period']
    if geometry['actual_N_width_coefficient_log_upper'] is not None:
        assert current.ep(geometry['actual_N_width_coefficient_log_upper'])[1] < -100*mp.log(10)
    for alternative in live['actual_signed_inverse_function_alternatives']:
        proof = alternative['original_inverse_and_Z_function_proof']
        if proof['geometry'] == 'flat':
            continue
        assert proof['exact_original_inverse_identity'] == 'A=a/2*(phi-psi_fraction)'
        assert proof['actual_selected_inverse_psi_used_not_support_cap']
        assert proof['original_Z_phi_derivative_enclosures_retained']
        assert proof['original_defining_inverse_primitive_function_unchanged']
        assert encoded(alternative['actual_signed_primitive_C0_Z_phi']['A'].record()) == encoded(proof['actual_inverse_identity_A_C0'])


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Source changed: '+name
    radius, measure, inverse = independent_phase(), independent_measure_and_transport(), independent_inverse_identity()
    formal_union = independent_formal_union()
    counts = dict(actual_source_function_queries=0, nonflat_directed_inverse_alternatives=0,
        flat_inverse_alternatives=0, direct_signed_primitive_rows=0, nonlinear_signed_density_rows=0,
        exact_directed_interval_rows=0, actual_signed_switch_integral_cells=0,
        positive_micro_phase_cells=0, actual_full_period_post_cells=0,
        actual_local_integral_rows=0, actual_incoming_memory_rows=0,
        actual_R110_correction_background_rows=0, exact_phase_joins=0, typed_rejections=0)
    with mp.workdps(540), patch.object(current.switch.WholeZSwitchFiniteN, 'contribution', forbidden), \
            patch.object(current.original, 'OriginalSwitchFiniteN', forbidden), \
            patch.object(current.backend.phase.first, 'NativePhaseFirstJets', forbidden), \
            patch.object(current.backend.signed, 'NativeSignedUDensityCover', forbidden):
        owner = current.WholeZSwitchSignedFunctions()
        assert saved['source_family'] == owner.identity and saved['candidate_N'] == owner.N == 2**3981
        for old in saved['actual_signed_switch_point_functions']:
            live = owner.query(tuple(old['exact_Z_cell']), old['actual_chart'], tuple(old['exact_left']))
            same_source(encoded(live), old, counts)
            verify(owner, live, counts)
            assert not live['actual_source_geometry']['full_period']
            for box in live['actual_source_geometry']['actual_phase_boxes']:
                lo, hi = current.ep(box)
                assert 0 <= lo <= hi <= 1 and hi-lo < mp.mpf('1e-70')
        for old in saved['actual_signed_switch_transports']:
            ends = tuple(old['exact_Z_cell'])
            live = owner.switch_transport(ends)
            same_source(encoded(live), old, counts)
            assert live['actual_R100_incoming_not_zeroed_or_replaced_by_background']
            assert live['interval_caps_not_substituted_for_defining_functions']
            assert live['actual_all_three_local_coordinates_in_zero_one']
            assert all(live[key] is False for key in current.OPEN)
            prior = live['actual_R100_correction_incoming_C0_Z']
            for window in live['actual_signed_switch_windows']:
                chart = window['actual_chart']
                assert encoded(window['actual_incoming_correction_C0_Z']) == encoded(prior)
                for cell in window['actual_signed_cells']:
                    integral = owner.local_integral(ends, chart, tuple(cell['exact_left']), tuple(cell['exact_right']))
                    assert current.function_digest(integral) == cell['current_signed_function_and_integral_sha256']
                    verify(owner, integral['actual_source_function'], counts)
                    assert integral['physical_Jacobian_applied_once']
                    assert integral['actual_positive_micro_width_not_zeroed']
                    assert integral['live_inverse_function_extension_under_Z_independent_integral']
                    if chart == 'post_power':
                        assert integral['actual_source_function']['actual_source_geometry']['full_period']
                        counts['actual_full_period_post_cells'] += 1
                    else:
                        assert not integral['actual_source_function']['actual_source_geometry']['full_period']
                        counts['positive_micro_phase_cells'] += 1
                    counts['actual_signed_switch_integral_cells'] += 1
                    counts['actual_local_integral_rows'] += 10
                for name in current.RATES:
                    for i in range(2):
                        retained = prior[name][i]*window['original_incoming_own_rate_memory'][name]
                        assert encoded(retained) == encoded(window['actual_retained_incoming_C0_Z'][name][i])
                        outgoing = retained+window['actual_signed_local_driver_C0_Z'][name][i]
                        assert encoded(outgoing) == encoded(window['actual_exit_signed_correction_C0_Z'][name][i])
                        counts['actual_incoming_memory_rows'] += 1
                pressure_memory = window['original_incoming_own_rate_memory']['p']
                assert current.ep(pressure_memory.finite_interval(max_log=2000)) == (1,1)
                prior = window['actual_exit_signed_correction_C0_Z']
            assert encoded(prior) == encoded(live['actual_R110_signed_correction_C0_Z'])
            for name in current.RATES:
                for i in range(2):
                    complete = live['actual_original_R110_background_C0_Z'][name][i]+prior[name][i]
                    assert encoded(complete) == encoded(live['actual_R110_signed_complete_own_history_C0_Z'][name][i])
                    counts['actual_R110_correction_background_rows'] += 1
            for a, b in (('first_switch','second_switch'), ('second_switch','post_power')):
                first, second = owner.radius_query(ends, a, (1,1)), owner.radius_query(ends, b, (0,1))
                assert encoded(first['actual_phase_boxes']) == encoded(second['actual_phase_boxes'])
                assert encoded(first['original_microscopic_width_coefficient']) == encoded(second['original_microscopic_width_coefficient'])
                counts['exact_phase_joins'] += 1
        assert counts['actual_signed_switch_integral_cells'] == 24
        assert counts['actual_source_function_queries'] == 36
        assert counts['positive_micro_phase_cells'] == 16 and counts['actual_full_period_post_cells'] == 8
        ends = current.CELLS[0]
        for call in (lambda: owner.radius_query(ends, 'first_switch', (-1,1)),
                lambda: owner.radius_query(ends, 'second_switch', (2,1)),
                lambda: owner.radius_query(ends, 'post_power', (1,1), (0,1)),
                lambda: owner.local_integral(ends, 'first_switch', (1,1), (1,1)),
                lambda: owner.radius_query(ends, 'first_switch', (1,0))):
            try:
                call()
            except ValueError:
                counts['typed_rejections'] += 1
            else:
                raise AssertionError('Invalid original switch request admitted')
    receipt = dict(all_passed=True, **{current.GATE: True}, source_family=owner.identity, candidate_N=owner.N,
        all_three_actual_switch_inverse_Z_phase_and_C1_integral_functions_checked=True,
        positive_micro_width_global_phase_and_original_coordinate_maps_checked=True,
        genuine_R100_incoming_three_window_suffix_transport_and_R110_background_separation_checked=True,
        exact_source_P0_rebinding_not_selected_or_recomputed=True,
        exact_conditioned_inverse_A_identity_retains_actual_psi_and_genuine_derivatives=True,
        independent_actual_radius_phase=radius, independent_actual_measure_and_transport=measure,
        independent_actual_inverse_identity=inverse,
        independent_signed_formal_branch_union=formal_union,
        replay_counts=counts, actual_full_prefix_C1_defect_functions_installed=False,
        actual_terminal_controls_installed=False, actual_global_frequency_admitted=False,
        **dict.fromkeys(current.OPEN, False), input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
            Path(__file__).name: current.sha(Path(__file__).name),
            Path(backend_check.__file__).name: current.sha(Path(backend_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('Current whole-Z switch signed functions, positive widths and C1 transport checks passed', flush=True)
    return receipt


if __name__ == '__main__':
    run()
