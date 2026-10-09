"""Check actual long source frontends, phase, signed C1 integrals and memory."""
import ast
import contextlib
from fractions import Fraction
import gzip
import hashlib
import inspect
import io
import json
from pathlib import Path
import textwrap
from types import SimpleNamespace
import time
from unittest.mock import patch
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_whole_Z_long_Rm_signed_functions as current
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
    P = c.exp(40)+11
    logRa = c.ln(4)-4*P-1000
    f = MacroFlow(c, c.mpf(-600), P, c.ln(4), logRa, c.ln(100)-logRa, c.mpf('.01'))
    op = SimpleNamespace(flow=f, c=c, Rm_factor=f.scalar(5))
    mapper = current.CurrentLongRadiusPhase.__new__(current.CurrentLongRadiusPhase)
    mapper.op, mapper.identity, mapper.N = op, {'manufactured': True}, 257
    mapper.T, mapper.logC = c.mpf(12), c.mpf(5)
    mapper.gap = 10*(mapper.logC+P)-mapper.T-8
    mapper.base = SimpleNamespace(sc=c.mpf('.25'), parameter_binding={'manufactured': True})
    mapper.cache = {}
    comparisons = 0
    def scalar(chart, q):
        t = p.mpf(q.numerator)/q.denominator
        Ps, Ts, Cs = p.exp(40)+11, p.mpf(12), p.mpf(5)
        radius = p.log(110)+Ts*t if chart == 'long_reshape' else \
            p.log(110)+Ts+(10*(Cs+Ps)-Ts-8)*t if chart == 'reference' else \
            p.log(110)+10*(Cs+Ps)+(-8 if chart == 'restoration' else -7)+t
        wanted = mapper.N*(radius-(p.log(4)-4*Ps-1000)-p.exp(-600)/8)
        return wanted-p.floor(wanted)
    for chart in current.CHARTS:
        for point in ((0,1), (1,3), (1,1)):
            got = mapper.query(chart, point)
            wanted = scalar(chart, Fraction(*point))
            assert any(current.ep(box)[0] <= wanted <= current.ep(box)[1] for box in got['actual_phase_boxes'])
            assert not got['full_period'] and got['actual_positive_width_retained_not_zeroed']
            assert got['actual_phase_Z_exact_zero'] and got['global_phase_not_restarted']
            for proof in got['exact_selected_dyadic_modulus_proofs']:
                assert not proof['huge_integer_materialized']
            comparisons += 1
        cell = mapper.query(chart, (0,1), (1,2))
        assert cell['full_period'] and cell['full_period_proved_from_actual_N_physical_width']
    for first, second in zip(current.CHARTS, current.CHARTS[1:]):
        a, b = mapper.query(first, (1,1)), mapper.query(second, (0,1))
        assert encoded(a['actual_phase_boxes']) == encoded(b['actual_phase_boxes'])
    narrow = mapper.query('long_reshape', (0,1), (1,10000))
    assert not narrow['full_period'] and not narrow['full_period_proved_from_actual_N_physical_width']
    for q in (Fraction(0), Fraction(1,20000), Fraction(1,10000)):
        wanted = scalar('long_reshape', q)
        assert any(current.ep(box)[0] <= wanted <= current.ep(box)[1] for box in narrow['actual_phase_boxes'])
        comparisons += 1
    giant = c.mpf(2)**10000
    box, proof = current.backend.signed.first.spatial.binary_mod_one(c, giant, Fraction(mapper.N, 3))
    exact = Fraction((pow(2,10000,3)*mapper.N) % 3, 3)
    wanted = p.mpf(exact.numerator)/exact.denominator
    assert current.ep(box)[0] <= wanted <= current.ep(box)[1]
    assert proof['exact_fraction'] == dict(numerator=exact.numerator, denominator=exact.denominator)
    assert not proof['huge_integer_materialized']
    return dict(passed=True, independent_scalar_original_radius_phase_comparisons=comparisons,
        exact_three_long_reference_restoration_post_joins=3,
        full_period_images_proved_from_actual_positive_lengths=4,
        narrow_monotone_cell_samples=3, exact_enormous_selected_T_modulus_checked=True,
        original_signed_positive_inlet_width_retained=True)


def independent_measure_and_transport():
    c = MPIntervalContext()
    c.dps = 180
    p = mp.mp.clone()
    p.dps = 240
    f = MacroFlow(c, c.mpf(-30), c.ln(9), c.ln(4), c.ln(5), c.ln(100)-c.ln(5), c.mpf('.01'))
    lengths = dict(zip(current.CHARTS, ('1.25', '2.25', '1', '1')))
    drivers = dict(zip(current.CHARTS, ('-.37', '.13', '-.21', '.09')))
    derivatives = dict(zip(current.CHARTS, ('.021', '-.015', '.04', '-.006')))
    op = SimpleNamespace(flow=f, c=c, P0=[], Rm_factor=f.scalar(5))
    source = SimpleNamespace(reference=SimpleNamespace(gap=c.mpf('2.25')))
    owner = current.WholeZLongRmSignedFunctions.__new__(current.WholeZLongRmSignedFunctions)
    owner.c, owner.N, owner.identity = c, 160, {'manufactured': True}
    owner.owner = lambda ends: op
    owner.long_owner = lambda ends: source
    owner.long_source = SimpleNamespace(T=c.mpf('1.25'))
    incoming = {name: [f.scalar('.73'), f.scalar('-.04')] for name in current.RATES}
    owner.switch_transport = lambda ends: dict(actual_R110_signed_correction_C0_Z=incoming)
    def query(ends, chart, *args):
        return dict(actual_signed_nonlinear_density_C0_Z=dict(
            kernels={name: f.scalar(drivers[chart]) for name in current.RATES},
            Z_derivatives={name: f.scalar(derivatives[chart]) for name in current.RATES}),
            actual_signed_primitive_C0_Z_phi={name: f.scalar(0) for name in current.backend.OUTPUTS},
            actual_source_geometry=dict(full_period=False), actual_original_source_packet=dict(original_generic_source=dict(
                common_own_five_histories_axial5={name: [f.scalar('-.2'), f.scalar('.03')] for name in current.RATES})))
    owner.query = query
    comparisons = 0
    for chart in current.CHARTS:
        got = owner.local_integral(('0','0'), chart, (1,4), (3,4))
        W = p.mpf(lengths[chart])/2
        assert current.ep(got['actual_physical_log_width'])[0] <= W <= current.ep(got['actual_physical_log_width'])[1]
        for name, rate in current.RATES.items():
            lam = p.mpf(rate.numerator)/rate.denominator
            mass = -p.expm1(-lam*W)/lam if lam else W
            for row, want in zip(got['actual_signed_local_integral_C0_Z'][name],
                    (p.mpf(drivers[chart])*mass, p.mpf(derivatives[chart])*mass)):
                lo, hi = current.ep(row.finite_interval(max_log=2000))
                assert lo <= want <= hi
                comparisons += 1
    with contextlib.redirect_stdout(io.StringIO()):
        got = owner.long_transport(('0','0'))
    total = sum((p.mpf(lengths[chart]) for chart in current.CHARTS), p.mpf(0))
    for name, rate in current.RATES.items():
        lam = p.mpf(rate.numerator)/rate.denominator
        wants = [p.mpf('.73')*p.exp(-lam*total), -p.mpf('.04')*p.exp(-lam*total)]
        end = p.mpf(0)
        for chart in current.CHARTS:
            length = p.mpf(lengths[chart])
            end += length
            weight = p.exp(-lam*(total-end))*(-p.expm1(-lam*length)/lam if lam else length)
            wants[0] += p.mpf(drivers[chart])*weight
            wants[1] += p.mpf(derivatives[chart])*weight
        for row, want in zip(got['actual_Rm_signed_correction_C0_Z'][name], wants):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
        for row, want in zip(got['actual_Rm_signed_complete_own_history_C0_Z'][name],
                (wants[0]-p.mpf('.2'), wants[1]+p.mpf('.03'))):
            lo, hi = current.ep(row.finite_interval(max_log=2000))
            assert lo <= want <= hi
            comparisons += 1
    return dict(passed=True, independent_T_gap_one_one_C0_Z_kernel_mass_comparisons=40,
        independent_four_piece_signed_global_integrals_incoming_and_background_comparisons=comparisons-40,
        genuine_nonzero_incoming_memory_and_distinct_signed_drivers_checked=True,
        true_positive_Jacobians_and_suffixes_applied_once=True)


def independent_frontend_binding():
    node = ast.parse(textwrap.dedent(inspect.getsource(current.source_module.WholeZLongRmFiniteN.query))).body[0]
    stop = next(i for i, statement in enumerate(node.body) if isinstance(statement, ast.Assign)
        and isinstance(statement.value, ast.Call) and ast.unparse(statement.value.func) == 'primitives.all_u_primitive_bounds')
    prefix = ast.Module(body=node.body[:stop], type_ignores=[])
    digest = hashlib.sha256(ast.dump(prefix).encode()).hexdigest()
    assert digest == current.FRONTEND_BINDING['exact_original_frontend_AST_sha256']
    calls = [ast.unparse(n.func) for n in ast.walk(prefix) if isinstance(n, ast.Call)]
    for name in ('long.source_packet', 'long.source_quotients', 'reference.background_cell',
            'reference.recover_cell', 'reference.quotient_cell'):
        assert name in calls
    assert not any('primitive_bounds' in name or 'contribution' in name for name in calls)
    return dict(passed=True, defining_long_and_reference_source_frontend_AST_unchanged=True,
        same_original_pure_source_callbacks_and_roots_retained=True,
        no_support_value_or_ancestor_contribution_in_new_frontend=True)


def verify(owner, live, counts):
    backend_check.verify(owner, live, counts)
    packet, geometry = live['actual_original_source_packet'], live['actual_source_geometry']
    assert packet['identical_source_P0_tuples_rebound_to_live_descendant_object']
    assert packet['actual_original_source_frontend_binding'] == current.FRONTEND_BINDING
    assert geometry['exact_defining_T_400_A_and_live_long_reference_tuples_bound']
    assert geometry['actual_positive_width_retained_not_zeroed']
    assert current.ep(geometry['actual_N_inlet_width_phase_log_upper'])[1] < -100*mp.log(10)
    for before, after in zip(packet['original_long_independent_P0_axial5'], live['exact_common_P0_axial5']):
        assert encoded(before) == encoded(after)
        assert before.ctx is after.ctx and before.scale.bases is after.scale.bases and before.ledger is after.ledger
    for proof in geometry['exact_selected_dyadic_modulus_proofs']:
        assert not proof['huge_integer_materialized']
        selected = geometry['exact_live_source_T'] if proof['component'] == 'T' else geometry['exact_source_logCstar']
        assert proof['source_exact_mpf_tuple'] == list(current.ep(selected)[0]._mpf_)
    for alternative in live['actual_signed_inverse_function_alternatives']:
        proof = alternative['original_inverse_and_Z_function_proof']
        if proof['geometry'] != 'flat':
            assert proof['exact_original_inverse_identity'] == 'A=a/2*(phi-psi_fraction)'
            assert proof['actual_selected_inverse_psi_used_not_support_cap']
            assert proof['original_Z_phi_derivative_enclosures_retained']
            assert encoded(alternative['actual_signed_primitive_C0_Z_phi']['A'].record()) == encoded(proof['actual_inverse_identity_A_C0'])


def run():
    began = time.monotonic()
    saved = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE]
    checked_switch = json.loads(gzip.decompress((current.HERE/current.previous.NAME).read_bytes()))
    checked_inlets = {tuple(row['exact_Z_cell']): row for row in checked_switch['actual_signed_switch_transports']}
    for name, digest in saved['input_hashes'].items():
        assert current.sha(name) == digest, 'Source changed: '+name
    phase, measure, frontend = independent_phase(), independent_measure_and_transport(), independent_frontend_binding()
    counts = dict(actual_source_function_queries=0, nonflat_directed_inverse_alternatives=0, flat_inverse_alternatives=0,
        direct_signed_primitive_rows=0, nonlinear_signed_density_rows=0, exact_directed_interval_rows=0,
        actual_signed_long_integral_cells=0, actual_full_period_cells=0, actual_local_integral_rows=0,
        actual_incoming_memory_rows=0, actual_Rm_correction_background_rows=0, exact_phase_joins=0,
        exact_accepted_Rm_phase_joins=0, typed_rejections=0)
    with mp.workdps(540), patch.object(current.source_module.WholeZLongRmFiniteN, 'query', forbidden), \
            patch.object(current.source_module.WholeZLongRmFiniteN, 'contribution', forbidden), \
            patch.object(current.source_module.long, 'OriginalLongReshapeFiniteN', forbidden), \
            patch.object(current.source_module.reference, 'OriginalReferenceRestoreFiniteN', forbidden), \
            patch.object(current.previous.switch.WholeZSwitchFiniteN, 'contribution', forbidden), \
            patch.object(current.backend.phase.first, 'NativePhaseFirstJets', forbidden), \
            patch.object(current.backend.signed, 'NativeSignedUDensityCover', forbidden):
        owner = current.WholeZLongRmSignedFunctions()
        assert saved['source_family'] == owner.identity and saved['candidate_N'] == owner.N == 2**3981
        for old in saved['actual_signed_long_point_functions']:
            live = owner.query(tuple(old['exact_Z_cell']), old['actual_chart'], tuple(old['exact_left']))
            same_source(encoded(live), old, counts)
            verify(owner, live, counts)
            assert not live['actual_source_geometry']['full_period']
            assert all(current.ep(box)[1]-current.ep(box)[0] < mp.mpf('1e-70')
                for box in live['actual_source_geometry']['actual_phase_boxes'])
        for old in saved['actual_signed_long_transports']:
            ends = tuple(old['exact_Z_cell'])
            live = owner.long_transport(ends)
            same_source(encoded(live), old, counts)
            inlet = checked_inlets[ends]
            assert current.function_digest(inlet) == live['actual_signed_R100_R110_inlet_transport_sha256']
            assert encoded(inlet['actual_R110_signed_correction_C0_Z']) == encoded(live['actual_R110_current_signed_correction_incoming_C0_Z'])
            assert live['current_signed_switch_output_not_ancestor_support_incoming']
            assert live['original_R110_background_remains_separate_from_signed_correction']
            assert all(live[key] is False for key in current.OPEN)
            prior = live['actual_R110_current_signed_correction_incoming_C0_Z']
            for window in live['actual_signed_long_windows']:
                chart = window['actual_chart']
                assert encoded(window['actual_incoming_correction_C0_Z']) == encoded(prior)
                for cell in window['actual_signed_cells']:
                    integral = owner.local_integral(ends, chart, tuple(cell['exact_left']), tuple(cell['exact_right']))
                    assert current.function_digest(integral) == cell['current_signed_function_and_integral_sha256']
                    verify(owner, integral['actual_source_function'], counts)
                    assert integral['physical_Jacobian_applied_once'] and integral['original_T_gap_one_one_windows_not_shortened']
                    geometry = integral['actual_source_function']['actual_source_geometry']
                    assert geometry['full_period'] and geometry['full_period_proved_from_actual_N_physical_width']
                    counts['actual_signed_long_integral_cells'] += 1
                    counts['actual_full_period_cells'] += 1
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
            assert encoded(prior) == encoded(live['actual_Rm_signed_correction_C0_Z'])
            for name in current.RATES:
                for i in range(2):
                    complete = live['actual_original_Rm_background_C0_Z'][name][i]+prior[name][i]
                    assert encoded(complete) == encoded(live['actual_Rm_signed_complete_own_history_C0_Z'][name][i])
                    counts['actual_Rm_correction_background_rows'] += 1
            for a, b in zip(current.CHARTS, current.CHARTS[1:]):
                first, second = owner.radius_query(ends, a, (1,1)), owner.radius_query(ends, b, (0,1))
                assert encoded(first['actual_phase_boxes']) == encoded(second['actual_phase_boxes'])
                counts['exact_phase_joins'] += 1
            point = owner.radius_query(ends, 'postrestore', (1,1))
            old_rm = owner.long_radius[ends].base.point((1,1), owner.N, decimal_digits=80)
            assert encoded(point['actual_phase_boxes']) == encoded(old_rm['phase_boxes'])
            counts['exact_accepted_Rm_phase_joins'] += 1
        assert counts['actual_source_function_queries'] == 48 and counts['actual_signed_long_integral_cells'] == 32
        ends = current.CELLS[0]
        for call in (lambda: owner.radius_query(ends, 'long_reshape', (-1,1)),
                lambda: owner.radius_query(ends, 'reference', (2,1)),
                lambda: owner.radius_query(ends, 'restoration', (1,1), (0,1)),
                lambda: owner.local_integral(ends, 'postrestore', (1,1), (1,1)),
                lambda: owner.source_query(ends, 'postrestore', (1,0))):
            try:
                call()
            except ValueError:
                counts['typed_rejections'] += 1
            else:
                raise AssertionError('Invalid current long/reference request admitted')
    receipt = dict(all_passed=True, **{current.GATE: True}, source_family=owner.identity, candidate_N=owner.N,
        all_four_original_long_reference_signed_inverse_C1_source_and_integral_functions_checked=True,
        source_selected_singleton_T_logC_and_adaptive_global_phase_checked=True,
        original_T_gap_one_one_lengths_incoming_memories_and_Rm_background_separation_checked=True,
        actual_signed_switch_inlet_used_without_ancestor_support_substitution=True,
        complete_source_P0_tuple_context_basis_ledger_binding_checked=True,
        independent_actual_long_radius_phase=phase, independent_actual_long_measure_and_transport=measure,
        independent_original_source_frontend_binding=frontend, replay_counts=counts,
        actual_full_prefix_C1_defect_functions_installed=False, actual_terminal_controls_installed=False,
        actual_global_frequency_admitted=False, **dict.fromkeys(current.OPEN, False),
        input_hashes={**owner.hashes, current.NAME: current.sha(current.NAME),
            Path(__file__).name: current.sha(Path(__file__).name),
            Path(backend_check.__file__).name: current.sha(Path(backend_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(receipt, indent=2)+'\n').encode())
    print('Current whole-Z long/reference signed source, phase and C1 transport checks passed', flush=True)
    return receipt


if __name__ == '__main__':
    run()
