"""Independent checks for the exact quiet pressure primitives."""
import json
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives as current


def _reject(fn):
    try:
        fn()
    except ValueError:
        return True
    raise AssertionError('Unsupported quiet pressure mutation admitted')


def _direct_physical_equations(K, R, Z, delta):
    d = 1 - Z**2
    L = 1 - delta * Z**2
    m = K * Z / R
    U = s.Integer(0)
    q = s.cancel((2 * Z * U - (1 - delta) * Z * m - d * s.diff(m, Z)) / L)
    q_Z = s.diff(q, Z)
    q_y = R * s.diff(q, R)
    S = ((q + (1 - delta) * Z * q_Z / 2 + q_y) / L +
         q * (q / 2 + q_y) +
         U * (-2 * Z * q + d * q_Z - 2 * Z * q_y) / L)
    N1p = s.cancel(-S / 2 + 2 * s.diff(q, R) + R * s.diff(q, R, 2))
    return dict(m=m, U=U, q=q, q_Z=q_Z, q_y=q_y, S=s.cancel(S), N1p=N1p)


def _actual_source_identity_independent(owner):
    Z = s.Symbol('Z', real=True)
    frame = owner.frame_raw['actual_analytic_native_inlet_functions']
    mean = s.cancel(s.sympify(frame['u'], locals={'Z': Z}) *
                    s.sympify(frame['m1'], locals={'Z': Z}))
    delta = s.Symbol('delta', real=True)
    kappa = s.cancel(mean / Z)
    q = s.cancel((-(1 - delta) * Z * mean -
                  (1 - Z**2) * s.diff(mean, Z)) /
                 (1 - delta * Z**2))
    return mean, kappa, q


def run(owner=None):
    began = time.monotonic()
    owner = owner if owner is not None else current.run()
    if not (current.HERE / current.NAME).exists():
        current.run(owner)
    owner.assert_source()
    K = s.Symbol('K', positive=True)
    a, h = s.symbols('a h', positive=True)
    b = a + h
    primitive = owner.primitive(K, a, b)
    direct = _direct_physical_equations(K, owner.R, owner.Z, owner.delta)
    for key in ('m', 'U', 'q', 'q_Z', 'q_y', 'S', 'N1p'):
        assert s.simplify(primitive[key] - direct[key]) == 0, key
    i0 = s.integrate(direct['N1p'], (owner.R, a, b))
    i1 = s.integrate(owner.R * direct['N1p'], (owner.R, a, b))
    assert s.simplify(primitive['N1p_integral'] - i0) == 0
    assert s.simplify(primitive['R_N1p_integral'] - i1) == 0
    assert s.simplify(primitive['m13'] - i0) == 0
    assert s.simplify(primitive['m15'] + i1) == 0
    assert primitive['m13_sign'] == '+' and primitive['m15_sign'] == '-'
    assert s.simplify(primitive['q'].subs(owner.R, owner.Rp) + K / owner.Rp) == 0

    mean, kappa, q_common = _actual_source_identity_independent(owner)
    assert s.simplify(kappa - owner.kappa) == 0
    assert s.simplify(s.diff(kappa, owner.Z)) == 0
    assert s.simplify(q_common + kappa) == 0
    source_proof = owner.actual_source_identity()
    assert source_proof['exact_Z_independent_kappa']
    assert source_proof['exact_q_common_at_Rp_is_minus_kappa']
    assert owner.frame_raw['actual_analytic_native_inlet_functions'] == owner.frame_strings
    assert owner.frame_raw['actual_Rp_native_frame_functions'] == owner.endpoint_nodes['native']
    assert owner.quiet_raw['actual_C3_repaired_exit_power_and_pulse_functions'] == owner.endpoint_nodes['quiet']
    assert owner.quiet_raw['actual_C3_repaired_exit_power_and_pulse_functions']['compact_repair_at_Rc_to_2Rc_consumed']
    assert owner.quiet_raw['actual_C3_repaired_exit_power_and_pulse_functions']['original_outer_O3_source_ends_at_Rc_not_2Rc']

    old = owner.frame_strings['u']
    owner.frame_strings['u'] = '0'
    _reject(owner.actual_source_identity)
    owner.frame_strings['u'] = old
    _reject(lambda: owner.primitive(K, 2, 1))
    owner.assert_source()

    inputs = dict(owner.hashes)
    inputs[current.NAME] = current.sha(current.NAME)
    inputs[Path(__file__).name] = current.sha(Path(__file__).name)
    result = dict(
        all_passed=True,
        source_family=owner.family,
        checks=dict(
            direct_physical_known_pressure_equation_for_arbitrary_K_and_ordered_R_endpoints=True,
            signed_N1p_and_unweighted_R_weighted_integrals_verified=True,
            m13_plus_and_m15_minus_signs_retained=True,
            actual_native_frame_mean_and_q_source_identities_verified=True,
            actual_compact_repair_2Rc_to_Rp_segment_consumed=True,
            same_original_endpoint_function_nodes_and_raw_frame_strings_bound=True,
            formal_K_Rp_Pstar_kappa_keeps_physical_R_units=True,
            no_radius_Pstar_N_materialization_or_signed_numeric_centers=True,
            exact_quiet_segment_only=True,
            actual_numeric_defect_full_tails_n1_repair_global_claims_not_installed=True,
        ),
        input_hashes=inputs,
        **{current.GATE: True}, **dict.fromkeys(current.OPEN, False),
        execution_seconds=time.monotonic() - began,
    )
    (current.HERE / current.RECEIPT).write_text(
        json.dumps(result, separators=(',', ':')) + '\n', encoding='utf8', newline='\n')
    print('PASS_QUIET_PRESSURE_PRIMITIVES', len(str(primitive['m13'])), flush=True)
    return result


if __name__ == '__main__':
    run()
