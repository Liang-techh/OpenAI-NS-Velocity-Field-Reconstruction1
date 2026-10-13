"""Independent calculus and source checks for the C3 repair-band pressure adapter."""
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source as current
import lei_ren_part1_paper_compliant_current_original_n1_radial_free_pressure_integrals as radial_free


def _B(value, R, left, right):
    return s.cancel(value.subs(R, right) - value.subs(R, left))


def _direct_fields(m, R, Z, delta):
    U = s.cancel(m + R * s.diff(m, R))
    d = 1 - Z**2
    L = 1 - delta * Z**2
    q = s.cancel((2 * Z * U - (1 - delta) * Z * m - d * s.diff(m, Z)) / L)
    q_Z = s.diff(q, Z)
    q_R = s.diff(q, R)
    q_y = R * q_R
    S = s.cancel((q + (1 - delta) * Z * q_Z / 2 + q_y) / L
                 + q * (q / 2 + q_y)
                 + U * (-2 * Z * q + d * q_Z - 2 * Z * q_y) / L)
    N1p = s.cancel(-S / 2 + 2 * q_R + R * s.diff(q, R, 2))
    return dict(U=U, q=q, q_Z=q_Z, q_R=q_R, q_y=q_y, S=S, N1p=N1p)


def independent_radial_free_fixture():
    """Use a nonlinear own mean and compare production identities to direct R integrals."""
    R, Z = s.symbols('R Z', positive=True, real=True)
    delta = s.Rational(2, 11)
    k = s.Rational(3, 7)
    eps = s.Rational(1, 13)
    # The perturbation is smooth on [1,2] and gives nonzero q_R, q_Z and U.
    m = k * Z / R + eps * Z**2 * (R - 1)**2 / R
    fields = _direct_fields(m, R, Z, delta)
    assert s.simplify(fields['U'] - (m + R * s.diff(m, R))) == 0
    i0 = {}
    i1 = {}
    for name, value in dict(q_Z=fields['q_Z'],
                            U_q_Z=fields['U'] * fields['q_Z'],
                            U_squared=fields['U']**2,
                            U_U_Z=fields['U'] * s.diff(fields['U'], Z),
                            q=fields['q'], q_squared=fields['q']**2).items():
        i0[name] = s.integrate(value, (R, 1, 2))
        i1[name] = s.integrate(R * value, (R, 1, 2))
    moments = dict(unweighted={key: i0[key] for key in ('q_Z', 'U_q_Z', 'U_squared', 'U_U_Z')},
                   weighted={key: i1[key] for key in ('q', 'q_Z', 'q_squared',
                                                        'U_q_Z', 'U_squared', 'U_U_Z')})
    boundaries = dict(
        R_q=_B(R * fields['q'], R, 1, 2),
        R_q_squared=_B(R * fields['q']**2, R, 1, 2),
        R_squared_q=_B(R**2 * fields['q'], R, 1, 2),
        R_squared_q_squared=_B(R**2 * fields['q']**2, R, 1, 2),
        R_U_squared=_B(R * fields['U']**2, R, 1, 2),
        R_squared_U_squared=_B(R**2 * fields['U']**2, R, 1, 2),
        q_plus_R_q_R=_B(fields['q'] + R * fields['q_R'], R, 1, 2),
        R_squared_q_R=_B(R**2 * fields['q_R'], R, 1, 2))
    reduced = radial_free.fixed_pressure_integrals(Z, delta, moments, boundaries)
    direct0 = s.integrate(fields['N1p'], (R, 1, 2))
    direct1 = s.integrate(R * fields['N1p'], (R, 1, 2))
    assert s.simplify(reduced['N1p'] - direct0) == 0
    assert s.simplify(reduced['R_N1p'] - direct1) == 0

    # Independently recover the compact-band baseline and compare the
    # production symbolic constants without using archived value records.
    m0 = k * Z / R
    base = _direct_fields(m0, R, Z, delta)
    base0 = s.integrate(base['N1p'], (R, 1, 2))
    base1 = s.integrate(R * base['N1p'], (R, 1, 2))
    expected0 = k**2 / 8
    expected1 = -k**2 * s.log(2) / 4
    assert s.simplify(base0 - expected0) == 0
    assert s.simplify(-base1 - expected1) == 0
    generated = current.baseline_pressure_coefficients(k)
    assert s.simplify(generated['m13'] - expected0) == 0
    assert s.simplify(generated['m15'] - expected1) == 0
    return dict(nonlin_own_mean_identity=True,
                direct_pressure_identity=True,
                finite_endpoint_q_R_terms_checked=True,
                compact_baseline_Rc_coefficients_checked=True,
                fixture_delta=str(delta), fixture_k=str(k), fixture_epsilon=str(eps))


def _upper(value):
    if isinstance(value, dict):
        if 'upper_exact_mpf_tuple' in value:
            return mp.make_mpf(tuple(value['upper_exact_mpf_tuple']))
        if 'upper' in value:
            return mp.mpf(value['upper'])
    if hasattr(value, '_mpi_'):
        return current.ep(value)[1]
    return mp.mpf(str(value))


def _expansion_bounds(c, baseline_lower, baseline_upper, error_upper):
    """Recompute a safe log expansion around the independent baseline.

    The production source uses a logarithmic padding when the error is many
    orders below the baseline.  This independently derives the corresponding
    enclosure inequality, avoiding exponentiation of the actual astronomical
    source scales.
    """
    b0 = c.mpf(baseline_lower)
    b1 = c.mpf(baseline_upper)
    err = c.mpf(error_upper)
    gap = current.ep(b0 - err)[0]
    if gap > mp.mpf(100):
        # exp(err-b) <= exp(-100), log(1-r) >= -2r and log(1+r) <= r.
        tiny = c.exp(c.mpf(-100))
        lower = current.ep(b0 - 2 * tiny)[0]
        upper = current.ep(b1 + tiny)[1]
        mode = 'relative_gap_gt_100'
    else:
        ratio0 = c.exp(err - b0)
        ratio1 = c.exp(err - b1)
        lower = current.ep(b0 + c.log1p(-ratio0))[0]
        upper = current.ep(b1 + c.log1p(ratio1))[1]
        mode = 'direct_log1p_expansion'
    return lower, upper, mode


def _recorded_budget_fixture(owner, views):
    """Check recorded logarithmic intervals independently of production flags."""
    c = owner.c
    modes = {}
    for key, view in sorted(views.items()):
        modes['|'.join(key)] = {}
        for name in ('signed_pressure_m13', 'signed_pressure_m15'):
            row = view[name]
            coefficient = current.quiet.read_interval(c, row['log_absolute_coefficient'])
            baseline = current.quiet.read_interval(c, row['baseline_log'])
            error_record = row['perturbation_log_upper']
            assert not error_record['exact_zero']
            error_upper = _upper(error_record['log_absolute_upper'])
            coefficient_lower, coefficient_upper = current.ep(coefficient)
            baseline_lower, baseline_upper = current.ep(baseline)
            assert error_upper < baseline_lower
            expansion_lower, expansion_upper, mode = _expansion_bounds(
                c, baseline_lower, baseline_upper, error_upper)
            # The independently recomputed signed expansion must lie inside
            # the archived coefficient interval, including every cell.
            assert coefficient_lower < baseline_lower <= baseline_upper < coefficient_upper
            assert coefficient_lower <= expansion_lower <= baseline_lower
            assert baseline_upper <= expansion_upper <= coefficient_upper
            assert all(mp.isfinite(q) for q in
                       (coefficient_lower, coefficient_upper, baseline_lower,
                        baseline_upper, error_upper, expansion_lower,
                        expansion_upper))
            modes['|'.join(key)][name] = mode
    return dict(recorded_log_intervals_finite=True,
                recorded_error_strictly_below_baseline=True,
                recorded_signed_enclosure_contains_baseline=True,
                recorded_error_expansion_recomputed=True,
                recorded_error_expansion_modes=modes)


def focused_production_fixtures(owner, view, views=None):
    """Exercise source geometry, endpoint scaling, derivative selection, and delta."""
    c = owner.c
    bd = current.rg.Bounds(c)

    # The accepted native frame is the absolute radius source, not the phase
    # offset alone.  Its singleton node identity is checked before using the
    # logarithmic Rc interval for endpoint normalization.
    nodes = owner.native_radius_binding['expected_nodes']
    assert nodes == dict(original_Rc_phase_offset=4893,
                         original_log_r_minus=26684, original_logRc=26685)
    assert current.ep(owner.log_Rc)[0] > 0

    # Nonunit Rc catches an extra radius in the defining I0/I1 measures and
    # independently checks the exact dR and R dR powers in x=R/Rc.
    x, Rc, f = s.symbols('x Rc f', positive=True)
    assert s.simplify(s.integrate(f * Rc, (x, 1, 2)) - f * Rc) == 0
    assert s.simplify(s.integrate((Rc * x) * f * Rc, (x, 1, 2)) -
                      s.Rational(3, 2) * f * Rc**2) == 0
    for name, rows in owner.functions['defining_integrals'].items():
        expected = ('physical dR=Rc*dx exactly once' if name.startswith('I0_')
                    else 'physical R*dR=R*Rc*dx exactly once')
        for row in physical_rows(rows):
            assert owner.g.nodes[row.node]['measure'] == expected

    dq = current.rg.LogUpper.constant(c, 3)
    dqy = current.rg.LogUpper.constant(c, 5)
    scaled0, scaled1 = current._normalized_endpoint_errors(bd, dq, dqy, c.ln(3))
    unscaled0, unscaled1 = current._normalized_endpoint_errors(bd, dq, dqy, c.mpf(0))
    assert current.ep(unscaled0.log)[0] - current.ep(scaled0.log)[1] > c.mpf('.9')
    assert current.ep(unscaled1.log)[0] - current.ep(scaled1.log)[1] > c.mpf('.9')

    # Give U_Z and U_{yZ} deliberately different caps.  Only U.Z is valid.
    U = current.rg.JetBound(current.rg.LogUpper.constant(c, 2),
                            current.rg.LogUpper.constant(c, 3),
                            bd.zero, bd.zero)
    Uy = current.rg.JetBound(current.rg.LogUpper.constant(c, 2),
                             current.rg.LogUpper.constant(c, 7),
                             bd.zero, bd.zero)
    got_uuz = current._uuz_bound(U)
    wrong_uuz = U.value * Uy.Z
    assert current.ep(got_uuz.log)[1] < current.ep(wrong_uuz.log)[0]

    # A tiny but nonzero delta must survive the arithmetic path.
    tiny_delta = c.mpf('1e-1000')
    tiny_term = current._delta_u2_bound(c, bd, bd.one, tiny_delta)
    assert tiny_term.log is not None
    assert float('1e-1000') == 0.0
    assert current.ep(owner.delta)[1] > 0
    assert current.ep(owner.delta)[1] < c.mpf('.5')
    if views is None:
        views = {tuple(view['exact_Z_cell']): view}
    recorded = _recorded_budget_fixture(owner, views)
    return dict(absolute_native_Rc_binding=True,
                nonunit_Rc_measure_and_endpoint_scaling=True,
                derivative_mismatch_rejected=True,
                nonzero_tiny_delta_retained=True, **recorded)


def physical_rows(value):
    return (value.value, value.Z)


def _check_sample(owner, view):
    assert tuple(view['exact_Z_cell']) in owner.rows
    for name, sign in (('signed_pressure_m13', 'positive'),
                       ('signed_pressure_m15', 'negative')):
        value = view[name]
        assert value['sign'] == sign
        assert value['perturbation_smaller_than_baseline']
        assert value['signed_value_enclosure']
        assert value['log_absolute_coefficient']['lower'] != '0'
    assert view['finite_N_velocity_and_mean_retained']
    assert view['original_delta_upper_lt_half']
    assert view['inverse_L_majorant'] == '1/L <= 2 from 0 <= delta < 1/2'
    assert view['physical_common_scale_Pstar_applied_once']
    assert view['original_Rc_scale_retained']
    assert view['absolute_Rc_source']['native_frame_nodes'] == {
        'original_Rc_phase_offset': 4893,
        'original_log_r_minus': 26684,
        'original_logRc': 26685}
    assert view['absolute_Rc_source']['radius_not_materialized']
    assert view['finite_endpoint_q_R_retained']
    assert view['no_radius_or_N_materialization']
    assert view['no_outer_support_or_full_tail_claim']


def run(owner, packets):
    began = time.monotonic()
    raw = json.loads(gzip.decompress((current.HERE / current.NAME).read_bytes()))
    assert raw['candidate_actual_repair_band_pressure_source_constructed']
    assert not raw[current.GATE]
    for name, digest in raw['input_hashes'].items():
        assert current.sha(name) == digest, name
    owner.assert_graph()
    fixture = independent_radial_free_fixture()
    assert set(packets) == set(owner.rows)
    checked = {}
    for key, packet in packets.items():
        checked[key] = owner.report(packet)
        _check_sample(owner, checked[key])
    assert len(checked) == 4
    assert {tuple(view['exact_Z_cell']) for view in checked.values()} == set(owner.rows)
    default_dps = mp.mp.dps
    public = current.CurrentOriginalN1RepairBandPressureSource(
        owner.repair, owner.power, owner.quiet, require_checked=False)
    public_packet = public.evaluate(next(iter(owner.rows)))
    public_view = public.report(public_packet)
    assert mp.mp.dps == default_dps
    _check_sample(public, public_view)
    focused = focused_production_fixtures(owner, next(iter(checked.values())), checked)
    inputs = dict(owner.hashes)
    inputs[current.NAME] = current.sha(current.NAME)
    inputs[Path(__file__).name] = current.sha(Path(__file__).name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            accepted_repair_band_power_and_quiet_sources_bound=True,
            exact_four_Z_cell_cover=True,
            signed_m13_positive_and_m15_negative=True,
            finite_N_velocity_and_own_mean_corrections_retained=True,
            physical_Pstar_applied_once_and_Rc_scale_formal=True,
            radial_free_identity_and_finite_endpoint_q_R_terms=True,
            independent_nonlinear_physical_mean_fixture=fixture,
            focused_production_geometry_and_error_fixtures=focused,
            public_default_precision_evaluated_query=True,
            public_default_precision_dps=default_dps,
            no_radius_or_N_materialization=True,
            compact_repair_band_only=True,
            outer_support_full_tails_n1_higher_repair_and_global_false=True),
        exact_Z_cells=[list(key) for key in sorted(checked)],
        input_hashes=inputs,
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began)
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(result, separators=(',', ':')) + '\n', encoding='utf8', newline='\n')
    print('PASS_REPAIR_BAND_SIGNED_PRESSURE_SOURCE', len(checked), flush=True)
    return result


if __name__ == '__main__':
    raise SystemExit('run(owner, packets) requires the retained checked owner')
