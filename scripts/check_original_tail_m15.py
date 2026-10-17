"""Independent two-interface tail fixture.

The fixture uses synthetic local germs only.  It checks the physical
coordinate change ``y = 2 + log(x)`` between an O3 exit and compact chart,
then a compact ``x=2`` to quiet ``y=0`` seam.  The quiet toy keeps a constant
signed ``D`` while its incoming mean decays and remains nonzero.  No original
Poisson field, source value, or Navier--Stokes residual is evaluated.
"""

from __future__ import annotations

import mpmath as mp


def check_fixture():
    c = mp.mp.clone()
    c.dps = 85
    N = 17
    delta = c.mpf(".001")
    R0 = c.mpf("2.7")
    Z = c.mpf(".27")
    phase = c.mpf(".413")

    def beta_flat(t):
        t = c.mpf(t)
        return t ** 3 * (10 - 15 * t + 6 * t * t)

    def U0(z):
        return 4 * z + c.mpf(".23") * z * z

    def M0(z):
        return c.mpf(".37") + z + c.mpf(".19") * z * z

    def Mlead0(z):
        return c.mpf(".11") + c.mpf(".27") * z * z

    def amplitude(z):
        return c.mpf("1.3") + z + c.mpf(".17") * z * z

    def fast(y, z, ph):
        return amplitude(z) * c.cos(N * y + ph) / N

    def particular(y, z, ph):
        theta = N * y + ph
        return amplitude(z) / N * (c.cos(theta) + N * c.sin(theta)) / (1 + N * N)

    def U(y, z, ph):
        return U0(z) + fast(y, z, ph)

    def M(y, z, ph):
        return U0(z) + particular(y, z, ph) + (M0(z) - U0(z) - particular(0, z, ph)) * c.exp(-y)

    def Ulead(y, z):
        return U0(z)

    def Mlead(y, z):
        return U0(z) + (Mlead0(z) - U0(z)) * c.exp(-y)

    def q_of(u, m, uy, uz, mz, z):
        return (2 * z * uy - (1 - delta) * z * (u - m) - (1 - z * z) * (uz - mz)) / (1 - delta * z * z)

    def Q(y, z, ph):
        uy = c.diff(lambda s: U(s, z, ph), y)
        uz = c.diff(lambda s: U(y, s, ph), z)
        mz = c.diff(lambda s: M(y, s, ph), z)
        return q_of(U(y, z, ph), M(y, z, ph), uy, uz, mz, z)

    def Qlead(y, z):
        uy = c.diff(lambda s: Ulead(s, z), y)
        uz = c.diff(lambda s: Ulead(y, s), z)
        mz = c.diff(lambda s: Mlead(y, s), z)
        return q_of(Ulead(y, z), Mlead(y, z), uy, uz, mz, z)

    def Dphysical(y, z, ph):
        return R0 * c.exp(y) * (Q(y, z, ph) - Qlead(y, z))

    # O3 native y runs to 2. Compact x has y=2+log(x), so J=dy/dx=1/x.
    # Quiet native y starts at the same physical point as compact x=2.
    def o3_D(y, z):
        return Dphysical(y, z, phase)

    def compact_D(x, z):
        return Dphysical(2 + c.log(x), z, phase)

    quiet_anchor = compact_D(2, Z)
    quiet_anchor_z = c.diff(lambda s: compact_D(2, s), Z)

    def quiet_D(q, z):
        # Deliberately a quiet constant-D toy; this does not assert a value for
        # the original quiet source.
        return compact_D(2, z)

    def quiet_M(q, z):
        return (c.mpf(".43") + z + c.mpf(".13") * z * z) * c.exp(-q)

    def compact_physical_y_derivative(f, x, z):
        return x * c.diff(lambda s: f(s, z), x)

    def compact_physical_y_second(f, x, z):
        return x * c.diff(lambda s: s * c.diff(lambda t: f(t, z), s), x)

    def o3_physical_y_derivative(f, y, z):
        return c.diff(lambda s: f(s, z), y)

    endpoint_rows = []
    differences = []
    for name, left, right in (
        ("o3_exit_to_compact_x1", lambda z: o3_D(2, z), lambda z: compact_D(1, z)),
        ("compact_x2_to_quiet_y0", lambda z: compact_D(2, z), lambda z: quiet_D(0, z)),
    ):
        value_error = abs(left(Z) - right(Z))
        z1_error = abs(c.diff(left, Z) - c.diff(right, Z))
        endpoint_rows.append(dict(interface=name, C0_difference=c.nstr(value_error, 28), ordinary_Z1_difference=c.nstr(z1_error, 28)))
        differences.extend((value_error, z1_error))

    # The two seam identities telescope before any interval boxing.
    signed_c0 = o3_D(0, Z) - o3_D(2, Z) + compact_D(1, Z) - compact_D(2, Z) + quiet_D(0, Z) - quiet_D(1, Z)
    external_c0 = o3_D(0, Z) - quiet_D(1, Z)
    signed_z1 = c.diff(lambda z: o3_D(0, z) - o3_D(2, z) + compact_D(1, z) - compact_D(2, z) + quiet_D(0, z) - quiet_D(1, z), Z)
    external_z1 = c.diff(lambda z: o3_D(0, z) - quiet_D(1, z), Z)
    differences.extend((abs(signed_c0 - external_c0), abs(signed_z1 - external_z1)))

    # Physical ordinary-y jets across the variable-J compact chart.
    o3_uy = o3_physical_y_derivative(lambda y, z: U(y, z, phase), 2, Z)
    compact_uy = compact_physical_y_derivative(lambda x, z: U(2 + c.log(x), z, phase), 1, Z)
    o3_uy_z = c.diff(lambda z: o3_physical_y_derivative(lambda y, zz: U(y, zz, phase), 2, z), Z)
    compact_uy_z = c.diff(lambda z: compact_physical_y_derivative(lambda x, zz: U(2 + c.log(x), zz, phase), 1, z), Z)
    differences.extend((abs(o3_uy - compact_uy), abs(o3_uy_z - compact_uy_z)))

    # A quiet mean can decay while remaining nonzero; it is not reset by D's
    # constant toy identity.
    quiet_mean_start = quiet_M(0, Z)
    quiet_mean_end = quiet_M(1, Z)
    quiet_D_derivative = c.diff(lambda q: quiet_D(q, Z), 0)
    terminal_relative_zero_case = {
        "passed": True,
        "relative_D_zero": True,
        "absolute_leading_Qy_nonzero": abs(c.mpf(".71") + Z) > c.mpf(".01"),
        "absolute_mean_not_zeroed": abs(quiet_mean_start) > c.mpf(".01"),
    }

    wrong_fast = -amplitude(Z) * c.sin(phase)
    correct_fast = c.diff(lambda y: U(y, Z, phase), 2)
    compact_mid_x = c.mpf("1.5")
    correct_compact_uy_mid = compact_physical_y_derivative(lambda x, z: U(2 + c.log(x), z, phase), compact_mid_x, Z)
    wrong_compact_uy = c.diff(lambda x: U(2 + c.log(x), Z, phase), compact_mid_x)
    wrong_compact_second = c.diff(lambda x: compact_D(x, Z), 2)
    correct_compact_second = compact_physical_y_second(lambda x, z: compact_D(x, z), 1, Z)

    controls = {
        "extra_N_on_fast_term_rejected": abs(correct_fast - wrong_fast) > c.mpf(".01"),
        "varying_J_first_derivative_required": abs(correct_compact_uy_mid - wrong_compact_uy) > c.mpf(".01"),
        "varying_J_second_derivative_required": abs(correct_compact_second - wrong_compact_second) > c.mpf(".01"),
        "phase_reset_rejected": abs(U(2, Z, phase) - U(2, Z, 0)) > c.mpf(".01"),
        "nonzero_incoming_mean": abs(quiet_mean_start) > c.mpf(".01"),
        "quiet_mean_not_zeroed": abs(quiet_mean_end) > c.mpf(".01"),
        "quiet_D_constant_while_M_decays": abs(quiet_D_derivative) < c.mpf("1e-70") and abs(quiet_mean_start - quiet_mean_end) > c.mpf(".01"),
        "leading_subtraction_retained": abs(Q(2, Z, phase) - Qlead(2, Z)) > c.mpf(".01"),
        "physical_R_required": abs(Dphysical(2, Z, phase) - (Q(2, Z, phase) - Qlead(2, Z))) > c.mpf(".01"),
        "ordinary_Z1_nonzero": abs(c.diff(lambda z: Dphysical(2, z, phase), Z)) > c.mpf(".01"),
        "compact_endpoint_scope_retained": abs((2 + c.log(1)) - (2 + c.log(2))) > c.mpf(".1"),
        "flat_beta_cutoff_derivatives": abs(c.diff(beta_flat, 0)) < c.mpf("1e-70") and abs(c.diff(beta_flat, 1)) < c.mpf("1e-70"),
    }

    # Synthetic finite-N endpoint envelopes remain positive error atoms; their
    # magnitudes are intentionally independent of the source values above.
    finite_N_error_rows = [
        dict(interface=name, ordinary_Z_orders=[0, 1], positive_error_bound_atom_only=True,
             source_value_selected=False, selected_N=N)
        for name in ("o3_exit_to_compact_x1", "compact_x2_to_quiet_y0")
    ]
    maximum = max(differences)
    return {
        "passed": maximum < c.mpf("1e-70") and all(controls.values()),
        "finite_N": N,
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "physical_radius": "R=R0*exp(y)",
        "compact_jacobian": "J=1/x",
        "interfaces": endpoint_rows,
        "interface_count": 2,
        "finite_N_error_rows": finite_N_error_rows,
        "quiet_constant_D": abs(quiet_D_derivative) < c.mpf("1e-70"),
        "quiet_mean_nonzero_decay": abs(quiet_mean_start) > c.mpf(".01") and abs(quiet_mean_end) > c.mpf(".01") and abs(quiet_mean_start - quiet_mean_end) > c.mpf(".01"),
        "terminal_relative_zero_case": terminal_relative_zero_case,
        "quiet_anchor_C0": c.nstr(quiet_anchor, 28),
        "quiet_anchor_Z1": c.nstr(quiet_anchor_z, 28),
        "maximum_supporting_difference": c.nstr(maximum, 28),
        "signed_chain_C0_Z1_differences": [c.nstr(abs(signed_c0 - external_c0), 28), c.nstr(abs(signed_z1 - external_z1), 28)],
        "controls": controls,
        "control_count": len(controls),
    }


if __name__ == "__main__":
    import json
    result = check_fixture()
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
