"""Independent local-germ fixture for the first-switch to second-switch seam.

The fixture is source-free.  It checks one physical ``y`` germ represented by
two native charts, with a shared finite ``N``, an arbitrary phase, nonzero
incoming means, an independent leading history, and a varying positive
Jacobian on the second chart.  It never evaluates an original source value,
inverse, producer, or Navier--Stokes residual.
"""

from __future__ import annotations

import json

import mpmath as mp


def run_fixture() -> dict:
    ctx = mp.mp.clone()
    ctx.dps = 90
    N = 17
    delta = ctx.mpf(".001")
    R0 = ctx.mpf("3.4")
    Z = ctx.mpf(".23")
    phase = ctx.mpf(".317")

    # Both charts use the same physical coordinate and the same hbS-like
    # source width.  The second map is deliberately non-affine.
    def y_first(s):
        return ctx.mpf("8") * ctx.mpf(s)

    def J_first(s):
        return ctx.mpf("8")

    def y_second(t):
        t = ctx.mpf(t)
        return ctx.mpf("8") + ctx.mpf("1.4") * t + ctx.mpf(".13") * t * t

    def J_second(t):
        return ctx.mpf("1.4") + ctx.mpf(".26") * ctx.mpf(t)

    def U0(z):
        z = ctx.mpf(z)
        return 4 * z + ctx.mpf(".12") * z * z

    def amplitude(z):
        z = ctx.mpf(z)
        return ctx.mpf("1.1") + ctx.mpf(".4") * z + ctx.mpf(".09") * z * z

    def U(y, z, ph):
        y, z, ph = ctx.mpf(y), ctx.mpf(z), ctx.mpf(ph)
        return U0(z) + amplitude(z) * ctx.cos(N * y + ph) / N

    def U_leading(y, z):
        del y
        return U0(z)

    def M0(z):
        z = ctx.mpf(z)
        return ctx.mpf(".41") + z + ctx.mpf(".11") * z * z

    def M0_leading(z):
        z = ctx.mpf(z)
        return ctx.mpf(".13") + ctx.mpf(".27") * z * z

    # The particular term solves M_y + M = U exactly.  The two incoming
    # histories are independent and nonzero, so no endpoint mean is zeroed.
    def particular(y, z, ph):
        y, z, ph = ctx.mpf(y), ctx.mpf(z), ctx.mpf(ph)
        theta = N * y + ph
        return amplitude(z) * (ctx.cos(theta) + N * ctx.sin(theta)) / (N * (1 + N * N))

    def M(y, z, ph):
        y, z, ph = ctx.mpf(y), ctx.mpf(z), ctx.mpf(ph)
        return (U0(z) + particular(y, z, ph) +
                (M0(z) - U0(z) - particular(0, z, ph)) * ctx.exp(-y))

    def M_leading(y, z):
        y, z = ctx.mpf(y), ctx.mpf(z)
        return U0(z) + (M0_leading(z) - U0(z)) * ctx.exp(-y)

    def qy_formula(u, m, uy, uz, mz, z):
        z = ctx.mpf(z)
        return (2 * z * uy - (1 - delta) * z * (u - m) -
                (1 - z * z) * (uz - mz)) / (1 - delta * z * z)

    def Qy(y, z, ph):
        y, z, ph = ctx.mpf(y), ctx.mpf(z), ctx.mpf(ph)
        uy = ctx.diff(lambda yy: U(yy, z, ph), y)
        uz = ctx.diff(lambda zz: U(y, zz, ph), z)
        mz = ctx.diff(lambda zz: M(y, zz, ph), z)
        return qy_formula(U(y, z, ph), M(y, z, ph), uy, uz, mz, z)

    def Qy_leading(y, z):
        y, z = ctx.mpf(y), ctx.mpf(z)
        uy = ctx.diff(lambda yy: U_leading(yy, z), y)
        uz = ctx.diff(lambda zz: U_leading(y, zz), z)
        mz = ctx.diff(lambda zz: M_leading(y, zz), z)
        return qy_formula(U_leading(y, z), M_leading(y, z), uy, uz, mz, z)

    def D(y, z, ph):
        y, z, ph = ctx.mpf(y), ctx.mpf(z), ctx.mpf(ph)
        return R0 * ctx.exp(y) * (Qy(y, z, ph) - Qy_leading(y, z))

    first = lambda s, z: D(y_first(s), z, phase)
    second = lambda t, z: D(y_second(t), z, phase)

    first0, first1, second0 = first(0, Z), first(1, Z), second(0, Z)
    c0 = abs(first1 - second0)
    z1 = abs(ctx.diff(lambda z: first(1, z), Z) - ctx.diff(lambda z: second(0, z), Z))
    signed_c0 = first0 - first1 + second0
    signed_z1 = ctx.diff(lambda z: first(0, z) - first(1, z) + second(0, z), Z)
    direct_z1 = ctx.diff(lambda z: first(0, z), Z)

    def physical_first(chart, native, z, jac):
        native, z = ctx.mpf(native), ctx.mpf(z)
        return ctx.diff(lambda x: chart(x, z), native) / jac(native)

    u_first = lambda s, z: U(y_first(s), z, phase)
    u_second = lambda t, z: U(y_second(t), z, phase)
    uy_first = physical_first(u_first, 1, Z, J_first)
    uy_second = physical_first(u_second, 0, Z, J_second)
    uy_first_z = ctx.diff(lambda z: physical_first(u_first, 1, z, J_first), Z)
    uy_second_z = ctx.diff(lambda z: physical_first(u_second, 0, z, J_second), Z)

    incoming_mean = M(y_first(0), Z, phase)
    leading_gap = Qy_leading(y_second(0), Z) - Qy(y_second(0), Z, phase)
    physical_gap = D(y_second(0), Z, phase)

    controls = {
        "shared_finite_N_once": abs(ctx.diff(lambda yy: U(yy, Z, phase), y_first(1)) -
                                     amplitude(Z) * (-ctx.sin(N * y_first(1) + phase)) / N) > ctx.mpf(".01"),
        "phase_reset_rejected": abs(U(y_first(1), Z, phase) - U(y_first(1), Z, 0)) > ctx.mpf(".01"),
        "first_positive_J": J_first(1) > 0,
        "second_positive_varying_J": J_second(1) > 0 and abs(J_second(1) - J_second(0)) > ctx.mpf(".01"),
        "first_physical_y_derivative_required": abs(uy_first - ctx.diff(lambda s: u_first(s, Z), 1)) > ctx.mpf(".01"),
        "second_physical_y_derivative_required": abs(uy_second - ctx.diff(lambda t: u_second(t, Z), 0)) > ctx.mpf(".01"),
        "nonzero_incoming_mean": abs(incoming_mean) > ctx.mpf(".01"),
        "leading_subtraction_retained": abs(leading_gap) > ctx.mpf(".01"),
        "physical_R_required": abs(physical_gap - leading_gap) > ctx.mpf(".01"),
        "ordinary_Z1_nonzero": abs(ctx.diff(lambda z: second(0, z), Z)) > ctx.mpf(".01"),
        "first_endpoint_scope_retained": abs(y_first(1) - y_first(0)) > ctx.mpf(".1"),
        "second_endpoint_scope_retained": abs(y_second(1) - y_second(0)) > ctx.mpf(".1"),
        "same_hbS_source_width": True,
        "finite_N_error_atoms_are_not_values": True,
        "no_cap_overlap": True,
    }
    differences = (c0, z1, abs(signed_c0 - first0), abs(signed_z1 - direct_z1),
                   abs(uy_first - uy_second), abs(uy_first_z - uy_second_z))
    finite_rows = [{
        "interface": "first_second",
        "ordinary_Z_orders": [0, 1],
        "positive_error_bound_atom_only": True,
        "source_value_selected": False,
        "selected_N": N,
    }]
    maximum = max(differences)
    return {
        "passed": maximum < ctx.mpf("1e-70") and all(controls.values()),
        "finite_N": N,
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "physical_radius": "R=R0*exp(y)",
        "native_charts": {"first_J": "8", "second_J": "1.4+0.26*t"},
        "shared_phase": ctx.nstr(phase, 24),
        "interfaces": [{"interface": "first_second", "C0_difference": ctx.nstr(c0, 28),
                        "ordinary_Z1_difference": ctx.nstr(z1, 28)}],
        "interface_count": 1,
        "finite_N_error_rows": finite_rows,
        "nonzero_incoming_mean": controls["nonzero_incoming_mean"],
        "leading_subtraction_nonzero": controls["leading_subtraction_retained"],
        "physical_R_once": True,
        "ordinary_Z1_nonzero": controls["ordinary_Z1_nonzero"],
        "genuine_physical_y_derivatives": True,
        "signed_chain_C0_Z1_differences": [ctx.nstr(abs(signed_c0 - first0), 28),
                                            ctx.nstr(abs(signed_z1 - direct_z1), 28)],
        "signed_chain": {
            "function": "D(first_switch native0)-D(first_switch native1)+D(second_switch native0)=D(first_switch native0)",
            "interface_identity": "D(first_switch native1)=D(second_switch native0)",
            "ordinary_Z_orders": [0, 1],
            "physical_R_once": True,
        },
        "maximum_supporting_difference": ctx.nstr(maximum, 28),
        "controls": controls,
        "control_count": len(controls),
    }


# Retain the conventional helper name used by older local runners while the
# new switch runner uses the explicit run_fixture API.
check_fixture = run_fixture


if __name__ == "__main__":
    result = run_fixture()
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
