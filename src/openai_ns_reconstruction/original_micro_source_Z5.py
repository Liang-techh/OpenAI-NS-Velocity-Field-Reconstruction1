"""Extend the accepted ordinary micro source jets by one axial order.

This module consumes the signed E/a packet as an interval receipt.  It reads
the frozen analytic archives only to form the two new all-coefficient Xh
rows, then appends the fifth ordinary Z row to the already accepted jets.
No accepted producer, writer, recurrence, norm owner, or quadrature is run.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from .original_core_signed_direction import read, record, require, sha


BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_"
PARENT = PREFIX + "compliant_current_original_micro_Ea_mixed_bounds.json"
PARENT_SHA256 = "a48e8ba9de2f2865598d78b6ba7cfcc94f95ec4b94cdda45b7acb2c168413576"
PARENT_MODULE = "src/openai_ns_reconstruction/original_micro_Ea_mixed_bounds.py"
HELPER = "src/openai_ns_reconstruction/original_ordinary_jet_bounds.py"

TRANSFER = PREFIX + "compliant_core_transfer.json"
LINEAR = PREFIX + "shared_linear_resolvent.json"
TUBE = PREFIX + "shared_analytic_tube.json"
UNIFORM = PREFIX + "compliant_core_uniform_bounds.json"
C0_PARENT = PREFIX + "compliant_current_original_micro_direction_C0.json"


def _load(path):
    return json.loads(Path(path).read_bytes())


def _pin(root, path, expected):
    require(isinstance(expected, str) and sha(root / path) == expected,
            "Changed pinned source: " + path)
    return dict(path=path, expected_sha256=expected)


def _normalize_parent(root, parent):
    """Accept either the full signed report or its inner E/a packet."""
    disk = _load(root / PARENT)
    require(sha(root / PARENT) == PARENT_SHA256 and disk.get("checks", {}).get("passed") is True,
            "Accepted micro E/a parent report required")
    if isinstance(parent, dict) and "actual_micro_Ea_mixed_bounds" in parent:
        require(parent == disk, "Parent report differs from the accepted E/a report")
        return parent, parent["actual_micro_Ea_mixed_bounds"]
    require(isinstance(parent, dict) and parent == disk["actual_micro_Ea_mixed_bounds"],
            "Inner E/a packet differs from the accepted E/a report")
    return disk, parent


def _consume(c, packet, name):
    rows = packet["comparison_source_Z0_Z4_enclosures"][name]
    require(isinstance(rows, list) and len(rows) == 5,
            "Accepted source jet must contain exactly Z0..Z4: " + name)
    return [read(c, row) for row in rows]


def _sym(c, value):
    return c.mpf([-value.b, value.b])


def _append(old, value):
    require(len(old) == 5, "Only an accepted Z0..Z4 jet may be extended")
    return list(old) + [value]


def _validate_cached_rows(packet, uniform):
    """Check the accepted lower rows without recomputing any of them."""
    fresh = packet["analytic_derivative_provider"]["new_rows"]
    require(
        {(row["scaled_radial_order"], row["ordinary_Z_order"]) for row in fresh}
        == {(0, 4), (1, 3), (1, 4)},
        "Accepted analytic rows do not have the required fresh prefix",
    )
    cached = {(row["scaled_radial_order"], row["axial_order"]): row
              for row in uniform["derivative_bounds"]}
    required = {(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 2)}
    require(required <= set(cached), "Uniform lower derivative rows are incomplete")
    return dict(
        accepted_analytic_rows_consumed=True,
        accepted_fresh_prefix=sorted((row["scaled_radial_order"], row["ordinary_Z_order"])
                                     for row in fresh),
        uniform_cached_lower_rows=sorted(required),
    )


def _next_log(f, g):
    order = 5
    rest = sum((math.comb(order - 1, k) * f[k] * g[order - k]
                for k in range(1, order)), f[0] * 0)
    return (f[order] - rest) / f[0]


def _quotient_jet(c, numerator, denominator):
    """Ordinary quotient recurrence through the appended order."""
    q = [numerator[0] / denominator[0]]
    for order in range(1, 6):
        rest = sum((math.comb(order, k) * denominator[k] * q[order - k]
                    for k in range(1, order + 1)), denominator[0] * 0)
        q.append((numerator[order] - rest) / denominator[0])
    return q


def _next_exp(logjet, valuejet):
    order = 5
    return sum((math.comb(order - 1, k - 1) * logjet[k] * valuejet[order - k]
                for k in range(1, order + 1)), valuejet[0] * 0)


def _new_rows(c, transfer, linear, tube, uniform):
    """Form only (rho_i,Z_k)=(0,5),(1,5) from the frozen embedding."""
    h = read(c, tube["Xh_parameter"])
    correction = read(c, transfer["scaled_map_size_upper"])
    phi_norm = read(c, linear["Phi_model_Xh_norm_upper"]) + correction
    psi_norm = read(c, transfer["fresh_Psi_model_Xh_norm_upper"]) + correction
    epsilon = read(c, transfer["epsilon"])
    radius = c.mpf(41) / 10
    require((radius / 20).b < 1, "Xh embedding disk must remain inside its pole radius")
    rows = {}
    proof_rows = []
    for radial in (0, 1):
        order = 5
        factor = (c.mpf(math.factorial(radial + order)) /
                  ((order + 1) ** 2 * h ** order * 20 ** radial) /
                  (1 - radius / 20) ** (radial + order + 1))
        phi = phi_norm * factor
        uz = epsilon * psi_norm * factor
        rows[radial] = (phi.b, uz.b)
        proof_rows.append(dict(
            scaled_radial_order=radial,
            ordinary_Z_order=order,
            embedding=record(factor),
            Phi_absolute_bound=record(phi),
            Uz_minus_U0_absolute_bound=record(uz),
            factorial_normalization=False,
        ))
    return rows, dict(
        passed=True,
        new_rows=proof_rows,
        new_rows_only=[(0, 5), (1, 5)],
        old_uniform_and_analytic_rows_consumed=True,
        embedding="(i+k)!/[20^i*h^k*(k+1)^2*(1-rho_max/20)^(i+k+1)]",
        model_norm_plus_uniform_self_map_once=True,
        epsilon_applied_once_to_Psi_norm=True,
        rho_max="41/10",
        no_new_coefficient_or_temporal_recursion=True,
        ordinary_derivatives_no_factorial_discount=True,
    )


def _new_gaxis(c, z0, delta, j, old_g):
    from .original_ordinary_jet_bounds import jet_add, jet_mul

    zero = c.mpf(0)
    const = lambda value: [c.mpf(value)] + [zero] * 5
    z = [z0, c.mpf(1)] + [zero] * 4
    one = const(1)
    zz = jet_mul(c, z, z)
    d = jet_add(c, one, [-value for value in zz])
    L = jet_add(c, one, [-delta * value for value in zz])
    U0 = jet_add(c, [4 * value for value in z], const(j))
    H0 = jet_add(c, [(1 - delta) * value / 2 for value in z],
                 jet_mul(c, d, U0))
    denominator = jet_add(c, jet_mul(c, H0, H0), const((j / 500) ** 2))
    numerator = jet_mul(c, L, H0)
    # The first five gaxis rows are consumed from the accepted packet.  Only
    # the order-five quotient recurrence is newly formed.
    order = 5
    rest = sum((math.comb(order, k) * denominator[k] * old_g[order - k]
                for k in range(1, order + 1)), denominator[0] * 0)
    g5 = (numerator[order] - rest) / denominator[0]
    return d[5], L[5], g5, z[5]


def _pair_add(c, a, b):
    from .original_ordinary_jet_bounds import jet_add
    return jet_add(c, a[0], b[0]), jet_add(c, a[1], b[1])


def _pair_neg(a):
    return [-x for x in a[0]], [-x for x in a[1]]


def _pair_mul(c, a, b):
    from .original_ordinary_jet_bounds import jet_add, jet_mul
    return (jet_mul(c, a[0], b[0]),
            jet_add(c, jet_mul(c, a[1], b[0]), jet_mul(c, a[0], b[1])))


def _pair_inv(c, a):
    from .original_ordinary_jet_bounds import jet_mul, jet_reciprocal
    reciprocal = jet_reciprocal(c, a[0])
    return reciprocal, [-x for x in jet_mul(c, jet_mul(c, reciprocal, reciprocal), a[1])]


def _direction(c, values, scalar, packet):
    """Build order-four directions from order-five source jets."""
    from .original_ordinary_jet_bounds import jet_add, jet_derivative, jet_mul

    zero = [c.mpf(0)] * 5
    const = lambda value: ([c.mpf(value)] + zero[1:], zero)

    def val(name):
        if name in ("Z", "d", "L", "gaxis"):
            return values[name][:5], zero
        return values[name][:5], values[name + "_y"][:5]

    add = lambda a, b: _pair_add(c, a, b)
    scale = lambda a, value: _pair_mul(c, a, const(value))

    def mul(a, b):
        return _pair_mul(c, a, b)

    def inv(a):
        return _pair_inv(c, a)

    def derivative(name):
        return (jet_derivative(c, values[name]),
                jet_derivative(c, values[name + "_y"]))

    z, d, L, g = [val(name) for name in ("Z", "d", "L", "gaxis")]
    H, K, M = [val(name) for name in ("JH", "JK", "M")]
    HZ, KZ, MZ = [derivative(name) for name in ("JH", "JK", "M")]
    phi = val("comparison_Phi")
    delta, epsilon = scalar["delta"], scalar["epsilon"]
    A = add(add(scale(H, 1 - delta / 2),
                _pair_neg(scale(mul(z, HZ), (1 - delta) / 2))),
            add(_pair_neg(mul(d, KZ)), scale(mul(z, K), 2 * delta - 1)))
    B = add(scale(mul(z, H), (1 - delta) / 2), mul(d, K))
    W = add(const(1), _pair_neg(add(scale(mul(z, M), 1 - delta), mul(d, MZ))))
    T = add(scale(mul(mul(g, B), inv(phi)), c.mpf(1) / 2),
            scale(add(_pair_neg(W), scale(mul(A, inv(phi)), c.mpf(1) / 2)), epsilon))
    rho = ([scalar["rho"]] + zero[1:], [scalar["rho"]] + zero[1:])
    result = _pair_mul(c, _pair_mul(c, rho, inv(L)), T)
    old_D = [read(c, row) for row in packet["comparison_Dbar_Z0_Z3"]]
    old_Dy = [read(c, row) for row in packet["comparison_Dbar_yZ0_Z3"]]
    require(len(old_D) == 4 and len(old_Dy) == 4, "Accepted Dbar rows must be Z0..Z3")
    result[0][:4] = old_D
    result[1][:4] = old_Dy
    return result, B


def _source_pins(root, full):
    hashes = full["input_hashes"]
    implementation = full["implementation_hashes"]
    paths = [
        C0_PARENT,
        TRANSFER,
        TRANSFER[:-5] + ".py",
        TRANSFER[:-5] + "_check.json",
        LINEAR,
        TUBE,
        TUBE[:-5] + ".py",
        UNIFORM,
        UNIFORM[:-5] + ".py",
    ]
    pins = {PARENT: _pin(root, PARENT, PARENT_SHA256)}
    for path in paths:
        name = Path(path).name
        require(name in hashes, "Accepted parent input hash missing: " + name)
        pins[path] = _pin(root, path, hashes[name])
    require(PARENT_MODULE in implementation, "Accepted E/a module implementation pin missing")
    pins[PARENT_MODULE] = _pin(root, PARENT_MODULE, implementation[PARENT_MODULE])
    helper_expected = full["actual_micro_Ea_mixed_bounds"]["source_pins"][HELPER]["expected_sha256"]
    pins[HELPER] = _pin(root, HELPER, helper_expected)
    return pins


def _identity(full, packet, transfer, linear, tube, uniform):
    require(packet["passed"] is True and packet["source_family"] == full["source_family"],
            "Accepted E/a source family required")
    provenance = packet["source_provenance"]
    require(provenance["passed"] is True and provenance["source_producers_replayed"] is False,
            "Accepted direct source provenance required")
    require(full["actual_core_signed_direction"]["passed"] is True and
            full["actual_micro_direction_C0"]["passed"] is True,
            "Accepted signed core and C0 identities required")
    require(transfer["contraction_proved"] and transfer["analytic_core_exists_for_compliant_source"],
            "Frozen analytic transfer required")
    require(linear["analytic_linear_inverse_bound_certified"],
            "Frozen analytic linear inverse required")
    require(transfer["analytic_core_family_sha256"] == tube["analytic_core_family_sha256"] ==
            uniform["analytic_core_family_sha256"],
            "Analytic core family identity differs")
    require(uniform["implicit_source_sha256"] == full["source_family"]["implicit_source_sha256"],
            "Uniform implicit source identity differs")
    require(uniform["unique_analytic_field_bounded_not_sampled"] is True and
            uniform["leading_model"] == "sum_n (-chi*r/2)^n/[n!(n+1)!], 0<=chi<=1",
            "Uniform analytic model identity missing")
    return dict(
        passed=True,
        parent_report_sha256=PARENT_SHA256,
        source_family=full["source_family"],
        accepted_source_producers_replayed=False,
        analytic_core_family_sha256=uniform["analytic_core_family_sha256"],
        direct_parent_Ea_packet=True,
        signed_core_and_C0_identities=True,
    )


def extend(c, root, parent):
    """Return interval values, (Dbar,Dbar_y), proof, and source pins."""
    from mpmath.ctx_iv import MPIntervalContext
    from .original_ordinary_jet_bounds import absolute_majorant, bound_product, jet_add, jet_mul

    require(isinstance(c, MPIntervalContext), "An MPIntervalContext is required")
    root = Path(root)
    full, packet = _normalize_parent(root, parent)
    pins = _source_pins(root, full)
    transfer, linear, tube, uniform = [_load(root / path) for path in (TRANSFER, LINEAR, TUBE, UNIFORM)]
    identity = _identity(full, packet, transfer, linear, tube, uniform)

    names = ("core_Phi", "core_Phi_rho", "core_U", "core_U_rho", "core_logPhi",
             "comparison_logPhi", "comparison_Phi", "comparison_U", "JH", "JK", "M",
             "JH_y", "JK_y", "M_y", "comparison_Phi_y", "comparison_U_y", "gaxis", "L", "d", "Z")
    source = {name: _consume(c, packet, name) for name in names}
    scalar_parent = full["actual_core_signed_direction"]["parameters"]
    scalar = dict(
        epsilon=read(c, full["actual_micro_direction_C0"]["comparison_direction"]["bounds"]["epsilon"]),
        rho=c.mpf([4, c.mpf(41) / 10]),
        delta=read(c, scalar_parent["delta"]),
    )
    require(scalar["epsilon"].a > 0, "Accepted epsilon must be positive")
    cached_proof = _validate_cached_rows(packet, uniform)
    new_rows, row_proof = _new_rows(c, transfer, linear, tube, uniform)
    phi5 = _sym(c, c.mpf(new_rows[0][0]))
    pr5 = _sym(c, c.mpf(new_rows[1][0]))
    u5 = _sym(c, c.mpf(new_rows[0][1]))
    ur5 = _sym(c, c.mpf(new_rows[1][1]))

    Phi = _append(source["core_Phi"], phi5)
    Pr = _append(source["core_Phi_rho"], pr5)
    U = _append(source["core_U"], u5)
    Ur = _append(source["core_U_rho"], ur5)
    coreLog = _append(source["core_logPhi"], _next_log(Phi, source["core_logPhi"] + [c.mpf(0)]))
    ratio = _quotient_jet(c, Pr, Phi)
    length = 2 * c.mpf(1) / 1000
    rho_max = scalar["rho"].b
    barLog5 = coreLog[5] + _sym(c, length * rho_max * c.fabs(ratio[5]).b)
    barLog = _append(source["comparison_logPhi"], barLog5)
    barPhi5 = _next_exp(barLog, source["comparison_Phi"] + [c.mpf(0)])
    barPhi = _append(source["comparison_Phi"], barPhi5)
    barU5 = U[5] + _sym(c, length * rho_max * c.fabs(Ur[5]).b)
    barU = _append(source["comparison_U"], barU5)

    global_phi = [read(c, uniform["actual_Phi_upper"])] + source["core_Phi"][1:] + [c.mpf([0, phi5.b])]
    product_in = bound_product(c, global_phi, absolute_majorant(c, U))
    product_bar = jet_mul(c, barPhi, barU)
    H5 = _sym(c, max(global_phi[5].b, c.fabs(barPhi[5]).b))
    K5 = _sym(c, max(product_in[5].b, c.fabs(product_bar[5]).b))
    M5 = _sym(c, max(c.fabs(U[5]).b, c.fabs(barU[5]).b))
    H = _append(source["JH"], H5)
    K = _append(source["JK"], K5)
    M = _append(source["M"], M5)

    alpha = c.mpf([0, 1])
    phi_y6 = [alpha * scalar["rho"] * value for value in jet_mul(c, barPhi, ratio)]
    u_y6 = [alpha * scalar["rho"] * value for value in Ur]
    HY = _append(source["JH_y"], 2 * (barPhi[5] - H5))
    KY = _append(source["JK_y"], 2 * (product_bar[5] - K5))
    MY = _append(source["M_y"], barU5 - M5)
    comparison_phi_y = _append(source["comparison_Phi_y"], phi_y6[5])
    comparison_u_y = _append(source["comparison_U_y"], u_y6[5])

    z0, delta, j = source["Z"][0], scalar["delta"], read(c, scalar_parent["j"])
    d5, L5, g5, z5 = _new_gaxis(c, z0, delta, j, source["gaxis"])
    d = _append(source["d"], d5)
    L = _append(source["L"], L5)
    z = _append(source["Z"], z5)
    gaxis = _append(source["gaxis"], g5)

    values = dict(core_Phi=Phi, core_Phi_rho=Pr, core_U=U, core_U_rho=Ur,
                  core_logPhi=coreLog, comparison_logPhi=barLog,
                  comparison_Phi=barPhi, comparison_U=barU,
                  JH=H, JK=K, M=M, JH_y=HY, JK_y=KY, M_y=MY,
                  comparison_Phi_y=comparison_phi_y, comparison_U_y=comparison_u_y,
                  gaxis=gaxis, L=L, d=d, Z=z)
    directions, _ = _direction(c, values, scalar, packet)
    B0 = jet_add(c, [((1 - scalar["delta"]) / 2) * value
                     for value in jet_mul(c, z, H)], jet_mul(c, d, K))
    values["B0"] = B0
    proof = dict(
        passed=True,
        source_identity=identity,
        ordinary_Z_orders=list(range(6)),
        direction_ordinary_Z_orders=list(range(5)),
        old_source_Z0_Z4_consumed_without_recomputation=True,
        old_Dbar_Z0_Z3_consumed_byte_exact=True,
        new_Dbar_Z4_and_yZ4_only=True,
        new_core_rows=row_proof,
        cached_prefix=cached_proof,
        core_row_keys=["core_Phi", "core_Phi_rho", "core_U", "core_U_rho"],
        comparison_row5_formula="coreLog5, ratio5, barLog5, exp recurrence, and barU5",
        normalized_history_max_bounds=True,
        history_bounds_formula="H=max(globalPhi,abs(barPhi)); K=max(abs(Phi*U),abs(barPhi*barU)); M=max(abs(U),abs(barU))",
        cached_inlet_rows_0_to_4=True,
        cached_inlet_U0_formula="4*Z+j+accepted_Uz_minus_U0_error",
        physical_y_ODE_rows=True,
        geometry_rho_y_equals_rho=True,
        B0_formula="(1-delta)*Z*JH/2+(1-Z^2)*JK",
        gaxis_exact_rational_recurrence=True,
        ordinary_derivatives_no_factorial_discount=True,
        no_zero_padding=True,
        no_accepted_writer_or_producer_replayed=True,
    )
    return values, directions, proof, pins


__all__ = ["extend"]
