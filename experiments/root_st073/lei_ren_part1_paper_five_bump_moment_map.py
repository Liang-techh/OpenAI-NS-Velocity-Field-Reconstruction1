"""Exact linear/quadratic change of the five moments under 3+2 bumps.

An integration backend supplies directed integrals from the actual base field
and the chosen compact bumps. This module assembles all cross terms; it does
not choose arbitrary replacement profiles, infer missing integrals or claim a
solution. Axial Taylor jets can be used for every integral and coefficient.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_candidate_shared_inlet import radial_product, integral, finite_moments
from lei_ren_part1_paper_interval_taylor import IntervalTaylor, constant
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAMES = ("theta", "z", "theta_z", "z_theta", "p")
PHYSICAL_CONTROL_ORDER = ("xi1", "xi2", "xi3", "c1", "c2")
PAPER_CONTROL_ORDER = ("c1", "c2", "xi1", "xi2", "xi3")


def physical_defects_to_paper_rows(ctx, defects, logRm, Am, Z):
    """Convert actual moment defects to the existing (10.8) row convention.

    Am and Z carry their axial jets. Am must not be frozen to a scalar when
    derivatives are requested. Physical radial factors are applied through
    logarithms, and all five defects are required from the common profile.
    Output row order: z, z_weighted, theta, z_theta, p.
    """
    if set(defects) != set(NAMES):
        raise ValueError("Require all five actual physical moment defects")
    invAm = Am.reciprocal()
    invRm = ctx.exp(-ctx.mpf(logRm))
    invAngular = invAm*ctx.exp(-ctx.mpf(logRm)*ctx.mpf(3)/2)/ctx.sqrt(ctx.mpf(2))
    z = defects["z"]*invRm
    theta = defects["theta"]*invAngular
    weighted = defects["theta_z"]*invAngular-theta*Z*4
    ztheta = (defects["z_theta"]-defects["z"]*Z*8)*invRm*invAm*invAm
    pressure = defects["p"]*invAm*invAm
    return (z, weighted, theta, ztheta, pressure)


def paper_rows_to_physical_changes(ctx, rows, logRm, Am, Z):
    if len(rows) != 5:
        raise ValueError("Require all five paper rows")
    a1, a2, a3, a4, a5 = rows
    Rm = ctx.exp(ctx.mpf(logRm))
    angular = Am*ctx.exp(ctx.mpf(logRm)*ctx.mpf(3)/2)*ctx.sqrt(ctx.mpf(2))
    return dict(z=a1*Rm, theta=a3*angular,
                theta_z=(a2+a3*Z*4)*angular,
                z_theta=a4*Am*Am*Rm+a1*Z*Rm*8,
                p=a5*Am*Am)


def assemble_map(F, Uz, angular_bumps, axial_bumps, integrate_products, zero):
    """Assemble A*h+Q(h,h) in physical R units.

    integrate_products(weight, factors) returns int R^weight prod(factors)dR
    over the common compact support band. The backend must retain the actual
    base F/Uz, bump supports, and verified error/jet data. Q is symmetric and
    evaluated as a full double sum, so each angular-axial entry in theta_z
    contains half of the cross coefficient.
    """
    if len(angular_bumps) != 3 or len(axial_bumps) != 2:
        raise ValueError("The paper repair requires three angular and two axial controls")
    A = {name: [zero for _ in range(5)] for name in NAMES}
    Q = {name: [[zero for _ in range(5)] for _ in range(5)] for name in NAMES}
    I = integrate_products
    for i, f in enumerate(angular_bumps):
        A["theta"][i] = I(1, (f,))*2
        A["theta_z"][i] = I(1, (f, Uz))*2
        A["z_theta"][i] = I(1, (F, f))*(-2)
        A["p"][i] = I(0, (F, f))*2
        for k, other in enumerate(angular_bumps):
            Q["z_theta"][i][k] = -I(1, (f, other))
            Q["p"][i][k] = I(0, (f, other))
        for j, g in enumerate(axial_bumps, start=3):
            Q["theta_z"][i][j] = I(1, (f, g))
            Q["theta_z"][j][i] = Q["theta_z"][i][j]
    for j, g in enumerate(axial_bumps, start=3):
        A["z"][j] = I(0, (g,))
        A["theta_z"][j] = I(1, (F, g))*2
        A["z_theta"][j] = I(0, (Uz, g))*2
        for k, other in enumerate(axial_bumps, start=3):
            Q["z_theta"][j][k] = I(0, (g, other))
    return dict(moment_names=NAMES, control_order=PHYSICAL_CONTROL_ORDER, A=A, Q=Q,
                radial_units="physical R", symmetric_double_sum=True,
                actual_support_and_integral_certificate_required=True)


def evaluate_map(mapping, controls):
    if len(controls) != 5:
        raise ValueError("Require five repair controls")
    A, Q = mapping["A"], mapping["Q"]
    return {name: sum((A[name][i]*controls[i] for i in range(5)), A[name][0]*0)
                  +sum((Q[name][i][j]*controls[i]*controls[j]
                        for i in range(5) for j in range(5)), A[name][0]*0)
            for name in NAMES}


def control_jacobian(mapping, controls):
    """Exact control derivative, preserving each axial coefficient/interval."""
    if len(controls) != 5 or not mapping["symmetric_double_sum"]:
        raise ValueError("Require five controls and the symmetric quadratic map")
    A, Q = mapping["A"], mapping["Q"]
    return {name: [A[name][i]+sum((Q[name][i][j]*controls[j]*2
                                  for j in range(5)), A[name][0]*0)
                   for i in range(5)] for name in NAMES}


def run():
    ctx = MPIntervalContext()
    ctx.dps = 100
    with mp.workdps(140):
        order = 3
        z = IntervalTaylor.variable(ctx, ctx.mpf(".3"), order)
        c = lambda x: constant(ctx, x, order)
        F = [c(1), z]
        U = [z, c(0), c(2)-z]
        angular = [[c(1)], [c(0), c(1)], [c(0), c(0), c(1)]]
        axial = [[c(1)], [c(0), c(1)]]
        def integrate_products(weight, factors):
            product = [c(1)]
            for factor in factors:
                product = radial_product(product, factor)
            return integral(product, ctx.mpf(2), weight)-integral(product, ctx.mpf(1), weight)
        mapping = assemble_map(F, U, angular, axial, integrate_products, c(0))
        controls = [c(".01")*(1+z), c("-.02"), c(".005")*z,
                    c(".03"), c("-.01")*(1-z)]
        changes = evaluate_map(mapping, controls)
        def changed(poly, bumps, h):
            result = list(poly)+[c(0)]*max(0, max(map(len, bumps))-len(poly))
            for bump, value in zip(bumps, h):
                for n, coefficient in enumerate(bump):
                    result[n] = result[n]+coefficient*value
            return result
        Fc = changed(F, angular, controls[:3])
        Uc = changed(U, axial, controls[3:])
        def moments(f, u):
            def at(r):
                m = finite_moments(f, u, ctx.mpf(r))
                m["z_theta"] = m["u_squared"]-m["weighted_phi_squared"]
                return m
            upper, lower = at(2), at(1)
            return {name: upper[name]-lower[name] for name in NAMES}
        before, after = moments(F, U), moments(Fc, Uc)
        checks = 0
        differences = {}
        for name in NAMES:
            difference = after[name]-before[name]-changes[name]
            differences[name] = list(difference.coefficients)
            for value in difference.coefficients:
                lo, hi = endpoints(value)
                if not lo <= 0 <= hi:
                    raise AssertionError("Five-moment quadratic identity failed: " + name)
                checks += 1
        jacobian = control_jacobian(mapping, controls)
        direct_jacobian = assemble_map(Fc, Uc, angular, axial, integrate_products, c(0))["A"]
        jacobian_checks = 0
        for name in NAMES:
            for i in range(5):
                difference = jacobian[name][i]-direct_jacobian[name][i]
                for value in difference.coefficients:
                    lo, hi = endpoints(value)
                    if not lo <= 0 <= hi:
                        raise AssertionError("Control Jacobian disagrees with the changed-field variation")
                    jacobian_checks += 1
        Am = (c(3)+z).reciprocal()
        logRm = ctx.log(ctx.mpf(2))
        normalized = physical_defects_to_paper_rows(ctx, changes, logRm, Am, z)
        restored = paper_rows_to_physical_changes(ctx, normalized, logRm, Am, z)
        conversion_checks = 0
        for name in NAMES:
            for value in (restored[name]-changes[name]).coefficients:
                lo, hi = endpoints(value)
                if not lo <= 0 <= hi:
                    raise AssertionError("Physical/paper row transformation failed")
                conversion_checks += 1
        # The independent direct integral contains every cross/quadratic term;
        # angular-only or axial-only fits cannot satisfy this comparison.
        here = Path(__file__).parent
        report = dict(input_hashes={name: hashlib.sha256((here/name).read_bytes()).hexdigest()
                                    for name in ("lei_ren_part1_paper_candidate_shared_inlet.py",
                                                 "lei_ren_part1_paper_interval_taylor.py")},
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      paper_repair_structure="(10.8): A*h+Q_Z(h,h)=-d(Z)",
                      moment_definitions="(2.22)",
                      controlled_component_definitions=dict(F="F+sum_0^2 h_i*f_i",
                                                            Uz="Uz+sum_3^4 h_i*g_i"),
                      nonlinear_map_assembled=True,
                      direct_changed_field_integral_comparison_count=checks,
                      changed_field_control_jacobian_comparison_count=jacobian_checks,
                      physical_paper_row_conversion_comparison_count=conversion_checks,
                      existing_paper_map_module="lei_ren_part1_paper_five_bump_map.py",
                      existing_paper_control_order=PAPER_CONTROL_ORDER,
                      physical_control_order=PHYSICAL_CONTROL_ORDER,
                      independent_physical_identity_backend=True,
                      comparison_ordinary_axial_Taylor_order=order,
                      fixture_interval_identity_differences=differences,
                      angular_axial_cross_terms_retained=True,
                      pressure_F_squared_terms_retained=True,
                      actual_paper_bump_supports_generated=False,
                      actual_transition_defect_data_supplied=False,
                      nonlinear_repair_solved=False,
                      functional_five_moment_closure=False)
        Path(__file__).with_suffix(".json").write_text(json.dumps(encode(report), indent=2)+"\n", encoding="utf-8")
        print("Five-bump exact quadratic map: direct physical moment comparisons", checks,
              "repair solution and actual supports remain open", flush=True)
        return report


if __name__ == "__main__":
    run()
