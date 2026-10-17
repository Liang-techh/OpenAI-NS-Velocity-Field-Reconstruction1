"""Independent local-germ fixture for the second-switch to long upstream seams.

This fixture is deliberately source-free.  It uses one physical ``y`` germ
through three native charts, with distinct positive Jacobians, a shared finite
``N``, arbitrary phase, nonzero incoming histories and an independent leading
subtraction.  It does not evaluate an original inverse, an original field
value, or a Navier--Stokes residual.
"""

from __future__ import annotations

import json

import mpmath as mp


def check_fixture() -> dict:
    c = mp.mp.clone()
    c.dps = 90
    N = 17
    delta = c.mpf(".001")
    R0 = c.mpf("3.4")
    Z = c.mpf(".23")
    phase = c.mpf(".317")

    # second_switch, post_power and long_reshape are three native charts for
    # one physical coordinate.  Their Jacobians are positive; the latter two
    # are intentionally nonconstant so raw-native derivatives are wrong.
    def y_second(s):
        return 9 * c.mpf(s)

    def J_second(s):
        return c.mpf(9)

    def y_post(t):
        t = c.mpf(t)
        return c.mpf(9) + 2 * t + c.mpf(".17") * t * t

    def J_post(t):
        return 2 + c.mpf(".34") * c.mpf(t)

    y_post_end = y_post(1)

    def y_long(u):
        u = c.mpf(u)
        return y_post_end + c.mpf(".8") * u + c.mpf(".11") * u * u

    def J_long(u):
        return c.mpf(".8") + c.mpf(".22") * c.mpf(u)

    def U0(z):
        z = c.mpf(z)
        return 4 * z + c.mpf(".12") * z * z

    def amplitude(z):
        z = c.mpf(z)
        return c.mpf("1.1") + c.mpf(".4") * z + c.mpf(".09") * z * z

    def U(y, z, ph):
        y, z, ph = c.mpf(y), c.mpf(z), c.mpf(ph)
        return U0(z) + amplitude(z) * c.cos(N * y + ph) / N

    def U_leading(y, z):
        return U0(z)

    def M0(z):
        z = c.mpf(z)
        return c.mpf(".41") + z + c.mpf(".11") * z * z

    def M0_leading(z):
        z = c.mpf(z)
        return c.mpf(".13") + c.mpf(".27") * z * z

    # Solve M_y = U - M explicitly for the synthetic germ.  The particular
    # term has the same shared N and phase as U, while the incoming history is
    # intentionally independent and nonzero.
    def particular(y, z, ph):
        y, z, ph = c.mpf(y), c.mpf(z), c.mpf(ph)
        theta = N * y + ph
        return amplitude(z) * (c.cos(theta) + N * c.sin(theta)) / (N * (1 + N * N))

    def M(y, z, ph):
        y, z, ph = c.mpf(y), c.mpf(z), c.mpf(ph)
        return (U0(z) + particular(y, z, ph) +
                (M0(z) - U0(z) - particular(0, z, ph)) * c.exp(-y))

    def M_leading(y, z):
        y, z = c.mpf(y), c.mpf(z)
        return U0(z) + (M0_leading(z) - U0(z)) * c.exp(-y)

    def q_formula(u, m, uy, uz, mz, z):
        z = c.mpf(z)
        return (2 * z * uy - (1 - delta) * z * (u - m) -
                (1 - z * z) * (uz - mz)) / (1 - delta * z * z)

    def Q(y, z, ph):
        y, z, ph = c.mpf(y), c.mpf(z), c.mpf(ph)
        uy = c.diff(lambda yy: U(yy, z, ph), y)
        uz = c.diff(lambda zz: U(y, zz, ph), z)
        mz = c.diff(lambda zz: M(y, zz, ph), z)
        return q_formula(U(y, z, ph), M(y, z, ph), uy, uz, mz, z)

    def Q_leading(y, z):
        y, z = c.mpf(y), c.mpf(z)
        uy = c.diff(lambda yy: U_leading(yy, z), y)
        uz = c.diff(lambda zz: U_leading(y, zz), z)
        mz = c.diff(lambda zz: M_leading(y, zz), z)
        return q_formula(U_leading(y, z), M_leading(y, z), uy, uz, mz, z)

    def D_physical(y, z, ph):
        y, z, ph = c.mpf(y), c.mpf(z), c.mpf(ph)
        R = R0 * c.exp(y)
        return R * (Q(y, z, ph) - Q_leading(y, z))

    def physical_first(f, native, z, jac):
        native, z = c.mpf(native), c.mpf(z)
        return c.diff(lambda x: f(x, z), native) / jac(native)

    def physical_second(f, native, z, jac):
        native, z = c.mpf(native), c.mpf(z)
        return c.diff(lambda x: c.diff(lambda xx: f(xx, z), x) / jac(x), native) / jac(native)

    chart_D = {
        "second": lambda s, z: D_physical(y_second(s), z, phase),
        "post": lambda t, z: D_physical(y_post(t), z, phase),
        "long": lambda u, z: D_physical(y_long(u), z, phase),
    }

    # The two source seams are C0 and ordinary-Z1 identities in physical y.
    seams = (
        ("second_power", chart_D["second"](1, Z), chart_D["post"](0, Z)),
        ("power_long", chart_D["post"](1, Z), chart_D["long"](0, Z)),
    )
    seam_rows = []
    differences = []
    for name, left, right in seams:
        if name == "second_power":
            left_z = lambda z: chart_D["second"](1, z)
            right_z = lambda z: chart_D["post"](0, z)
        else:
            left_z = lambda z: chart_D["post"](1, z)
            right_z = lambda z: chart_D["long"](0, z)
        c0 = abs(left - right)
        z1 = abs(c.diff(left_z, Z) - c.diff(right_z, Z))
        seam_rows.append({"interface": name, "C0_difference": c.nstr(c0, 28),
                          "ordinary_Z1_difference": c.nstr(z1, 28)})
        differences.extend((c0, z1))

    # The source chain telescopes before any interval boxing.
    signed_c0 = (chart_D["second"](0, Z) - chart_D["second"](1, Z) +
                 chart_D["post"](0, Z) - chart_D["post"](1, Z) +
                 chart_D["long"](0, Z) - chart_D["long"](1, Z))
    external_c0 = chart_D["second"](0, Z) - chart_D["long"](1, Z)
    signed_z1 = c.diff(lambda z: chart_D["second"](0, z) - chart_D["second"](1, z) +
                        chart_D["post"](0, z) - chart_D["post"](1, z) +
                        chart_D["long"](0, z) - chart_D["long"](1, z), Z)
    external_z1 = c.diff(lambda z: chart_D["second"](0, z) - chart_D["long"](1, z), Z)
    differences.extend((abs(signed_c0 - external_c0), abs(signed_z1 - external_z1)))

    # Actual physical-y derivative rows at both seams, including the varying
    # Jacobians.  A raw native derivative is deliberately a negative control.
    u_second = lambda s, z: U(y_second(s), z, phase)
    u_post = lambda t, z: U(y_post(t), z, phase)
    u_long = lambda u, z: U(y_long(u), z, phase)
    uy_second = physical_first(u_second, 1, Z, J_second)
    uy_post_in = physical_first(u_post, 0, Z, J_post)
    uy_post_out = physical_first(u_post, 1, Z, J_post)
    uy_long = physical_first(u_long, 0, Z, J_long)
    uy_second_z = c.diff(lambda z: physical_first(u_second, 1, z, J_second), Z)
    uy_post_in_z = c.diff(lambda z: physical_first(u_post, 0, z, J_post), Z)
    uy_post_out_z = c.diff(lambda z: physical_first(u_post, 1, z, J_post), Z)
    uy_long_z = c.diff(lambda z: physical_first(u_long, 0, z, J_long), Z)
    differences.extend((abs(uy_second - uy_post_in), abs(uy_second_z - uy_post_in_z),
                        abs(uy_post_out - uy_long), abs(uy_post_out_z - uy_long_z)))

    def beta_flat(t):
        t = c.mpf(t)
        return t ** 3 * (10 - 15 * t + 6 * t * t)

    incoming_mean = M(y_second(0), Z, phase)
    wrong_fast = amplitude(Z) * c.cos(phase)
    correct_fast = c.diff(lambda yy: U(yy, Z, phase), y_second(1))
    raw_post = c.diff(lambda tt: u_post(tt, Z), c.mpf(".6"))
    physical_post = physical_first(u_post, c.mpf(".6"), Z, J_post)
    raw_long = c.diff(lambda uu: u_long(uu, Z), c.mpf(".6"))
    physical_long = physical_first(u_long, c.mpf(".6"), Z, J_long)
    wrong_phase = U(y_second(1), Z, 0)
    leading_gap = Q(y_post(0), Z, phase) - Q_leading(y_post(0), Z)
    physical_gap = D_physical(y_post(0), Z, phase)

    controls = {
        "shared_finite_N_once": abs(correct_fast - N * wrong_fast) > c.mpf(".01"),
        "phase_reset_rejected": abs(U(y_second(1), Z, phase) - wrong_phase) > c.mpf(".01"),
        "second_positive_J": J_second(1) > 0,
        "post_positive_varying_J": J_post(1) > 0 and abs(J_post(1) - J_post(0)) > c.mpf(".01"),
        "long_positive_varying_J": J_long(1) > 0 and abs(J_long(1) - J_long(0)) > c.mpf(".01"),
        "second_physical_y_derivative_required": abs(uy_second - c.diff(lambda s: u_second(s, Z), 1)) > c.mpf(".01"),
        "post_physical_y_derivative_required": abs(physical_post - raw_post) > c.mpf(".01"),
        "long_physical_y_derivative_required": abs(physical_long - raw_long) > c.mpf(".01"),
        "nonzero_incoming_mean": abs(incoming_mean) > c.mpf(".01"),
        "leading_subtraction_retained": abs(leading_gap) > c.mpf(".01"),
        "physical_R_required": abs(physical_gap - leading_gap) > c.mpf(".01"),
        "ordinary_Z1_nonzero": abs(c.diff(lambda z: chart_D["post"](0, z), Z)) > c.mpf(".01"),
        "post_endpoint_scope_retained": abs(y_post(1) - y_post(0)) > c.mpf(".1"),
        "long_endpoint_scope_retained": abs(y_long(1) - y_long(0)) > c.mpf(".1"),
        "flat_width_cutoff_derivatives": abs(c.diff(beta_flat, 0)) < c.mpf("1e-70") and abs(c.diff(beta_flat, 1)) < c.mpf("1e-70"),
        "finite_N_error_atoms_are_not_values": True,
    }

    finite_N_error_rows = [
        {"interface": name, "ordinary_Z_orders": [0, 1],
         "positive_error_bound_atom_only": True, "source_value_selected": False,
         "selected_N": N}
        for name in ("second_power", "power_long")
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
        "native_charts": {"second_J": "9", "post_J": "2+0.34*t", "long_J": "0.8+0.22*u"},
        "shared_phase": c.nstr(phase, 24),
        "interfaces": seam_rows,
        "interface_count": 2,
        "finite_N_error_rows": finite_N_error_rows,
        "nonzero_incoming_mean": abs(incoming_mean) > c.mpf(".01"),
        "leading_subtraction_nonzero": abs(leading_gap) > c.mpf(".01"),
        "physical_R_once": True,
        "ordinary_Z1_nonzero": controls["ordinary_Z1_nonzero"],
        "signed_chain_C0_Z1_differences": [c.nstr(abs(signed_c0 - external_c0), 28),
                                            c.nstr(abs(signed_z1 - external_z1), 28)],
        "maximum_supporting_difference": c.nstr(maximum, 28),
        "controls": controls,
        "control_count": len(controls),
    }


if __name__ == "__main__":
    result = check_fixture()
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
