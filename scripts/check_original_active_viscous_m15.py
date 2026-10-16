"""Independent active-seam finite-N viscous trace checks.

This fixture uses an affine native-coordinate seam and direct quadrature of
``-d_y(R Q_y)``.  It does not construct a production owner, load an archive,
or invoke an accepted ancestor producer.  The source-ID and endpoint
fingerprint checks exercise only the small exact APIs used by the owner.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_micro_viscous_m15 import (  # noqa: E402
    EndpointFingerprint,
    contract_trace_chain,
)
from openai_ns_reconstruction.original_active_viscous_m15 import (  # noqa: E402
    ActiveEndpointFingerprint,
)
from check_original_micro_viscous_m15 import integrate, q_bundle  # noqa: E402


Z_TEXT = "0.37"
EPS_TEXT = "0.17"
R0_TEXT = "2.3"
TOL_TEXT = "2e-48"


def amplitude_perturbation(ctx, z):
    """Nonconstant slow amplitude used only for the modified-minus-leading trace."""
    return ctx.mpf("0.15") + ctx.mpf("0.11") * z + ctx.mpf("0.07") * z * z


def amplitude_perturbation_z(ctx, z):
    return ctx.mpf("0.11") + ctx.mpf("0.14") * z


def radius_affine(ctx, y):
    """Physical radius with a convenient nonzero scale near y=-8."""
    return ctx.mpf(R0_TEXT) * ctx.exp(y + 8)


def delta_q_bundle(ctx, y, z):
    """Return the modified-minus-leading Qy jet, including slow-Z products."""
    base = q_bundle(ctx, y, z)
    eps = ctx.mpf(EPS_TEXT)
    b = amplitude_perturbation(ctx, z)
    bz = amplitude_perturbation_z(ctx, z)
    scale = eps * b
    scale_z = eps * bz
    return {
        "Qy": scale * base["Qy"],
        "QyZ": scale * base["QyZ"] + scale_z * base["Qy"],
        "Qy_y": scale * base["Qy_y"],
        "Qy_yZ": scale * base["Qy_yZ"] + scale_z * base["Qy_y"],
    }


def delta_rq(ctx, y, z, *, ordinary_z=False, include_radius=True):
    q = delta_q_bundle(ctx, y, z)
    if ordinary_z:
        value = q["QyZ"]
        derivative = q["Qy_yZ"]
    else:
        value = q["Qy"]
        derivative = q["Qy_y"]
    radius = radius_affine(ctx, y) if include_radius else ctx.mpf(1)
    return radius * value, radius * (value + derivative)


def native_integral(ctx, native_to_y, z, *, ordinary_z=False, jacobian=1, include_radius=True):
    """Integrate the physical signed derivative over a native interval [0,1]."""
    return integrate(
        ctx,
        lambda native: -jacobian
        * delta_rq(
            ctx,
            native_to_y(native),
            z,
            ordinary_z=ordinary_z,
            include_radius=include_radius,
        )[1],
        0,
        1,
        segments=12,
    )


def endpoint_orientation(ctx, left_y, right_y, z, *, ordinary_z=False):
    left = delta_rq(ctx, left_y, z, ordinary_z=ordinary_z)[0]
    right = delta_rq(ctx, right_y, z, ordinary_z=ordinary_z)[0]
    return left - right


def close(ctx, left, right, tolerance):
    return abs(left - right) <= tolerance * max(ctx.mpf(1), abs(left), abs(right))


def fingerprint_adversarial_fixture():
    """Ensure endpoint scope distinguishes radius/phase and rejects 0/0."""
    nodes = [
        {"operation": "bound_variable", "name": "coordinate"},
        {"operation": "bound_variable", "name": "phase"},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "analytic_unary", "name": "exp", "argument": 0},
        {
            "operation": "current_original_unique_analytic_core_component",
            "rho": 1,
            "Z_variable": 2,
            "source_family": "fixture-family",
            "unique_analytic_Banach_solution": True,
            "actual_twenty_term_operator": True,
            "admitted_solution_not_a_model_or_enclosure_center": True,
        },
        {
            "operation": "positive_quotient",
            "numerator": 3,
            "denominator": 0,
            "source_positive_certificate": "fixture positive denominator",
        },
        {
            "operation": "definite_integral",
            "variable": "coordinate",
            "lower": 2,
            "upper": 2,
            "integrand": 5,
        },
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
    ]
    f = EndpointFingerprint(nodes, {"coordinate": 0})
    radius_zero = f.at(4, {0: f.zero})
    radius_one = f.at(4, {0: f.one})
    if radius_zero == radius_one:
        raise AssertionError(("radius endpoint alias", radius_zero, radius_one))
    phase_zero = f.at(5, {0: f.zero, 1: f.zero})
    phase_one = f.at(5, {0: f.zero, 1: f.one})
    if phase_zero == phase_one:
        raise AssertionError(("phase endpoint alias", phase_zero, phase_one))
    empty = f.at(7, {0: f.one})
    if empty[1] != 0 or 7 not in f.empty:
        raise AssertionError(("empty endpoint integral", empty, f.empty))
    quotient_rejected = False
    try:
        f.at(6, {0: f.zero})
    except ValueError:
        quotient_rejected = True
    if not quotient_rejected:
        raise AssertionError("undefined endpoint quotient was accepted")
    return dict(
        passed=True,
        radius_endpoint_mismatch=True,
        phase_endpoint_mismatch=True,
        exact_empty_integral=True,
        undefined_quotient_rejected=True,
    )


def active_fingerprint_adversarial_fixture():
    """Exercise the active endpoint algebra without constructing production state."""
    # The two opaque leaves have different coordinate arguments.  Their
    # distinct identities must survive endpoint substitution, while the exact
    # additive identity a + ((b - a) - 8) = b - 8 must normalize.
    nodes = [
        {"operation": "bound_variable", "name": "coordinate"},
        {"operation": "bound_variable", "name": "Z"},
        {"operation": "bound_variable", "name": "phase"},
        {"operation": "exact_rational", "numerator": -8, "denominator": 1},
        {
            "operation": "current_original_leading_function_recipe",
            "coordinate": 0,
            "Z_variable": 1,
            "recipe": "opaque-a",
            "source_family": "active-fixture",
        },
        {
            "operation": "current_original_leading_function_recipe",
            "coordinate": 2,
            "Z_variable": 1,
            "recipe": "opaque-b",
            "source_family": "active-fixture",
        },
        {"operation": "negative", "argument": 4},
        {"operation": "sum", "arguments": [5, 6]},
        {"operation": "sum", "arguments": [7, 3]},
        {"operation": "sum", "arguments": [4, 8]},
        {"operation": "sum", "arguments": [5, 3]},
        {"operation": "sum", "arguments": [4, 6]},
    ]
    f = ActiveEndpointFingerprint(nodes, {"coordinate": 0})
    scope = {
        0: f.rational(0),
        1: f.rational(1),
        2: f.rational(1),
    }
    normalized_left = f.at(9, scope)
    normalized_right = f.at(10, scope)
    if normalized_left != normalized_right:
        raise AssertionError(("active additive normalization", normalized_left, normalized_right))
    cancelled = f.at(11, scope)
    if cancelled != f.zero:
        raise AssertionError(("nonrational plus negative did not cancel", cancelled, f.zero))
    opaque_a = f.at(4, scope)
    opaque_b = f.at(5, scope)
    if opaque_a == opaque_b:
        raise AssertionError(("opaque source arguments aliased", opaque_a, opaque_b))
    radius_phase_nodes = [
        {"operation": "analytic_unary", "name": "exp", "argument": 0},
        {"operation": "analytic_unary", "name": "exp", "argument": 2},
    ]
    # These nodes are independent of the algebra graph above but use the same
    # endpoint scope, so a coordinate/phase substitution mismatch remains
    # visible to the active fingerprint.
    offset = len(nodes)
    f.nodes.extend(radius_phase_nodes)
    radius_at_zero = f.at(offset, scope)
    radius_at_one = f.at(offset, {**scope, 0: f.rational(1)})
    phase_at_zero = f.at(offset + 1, {**scope, 2: f.rational(0)})
    phase_at_one = f.at(offset + 1, {**scope, 2: f.rational(1)})
    if radius_at_zero == radius_at_one or phase_at_zero == phase_at_one:
        raise AssertionError("active radius/phase mismatch was erased")
    return dict(
        passed=True,
        exact_additive_normalization=True,
        nonrational_negative_cancellation=True,
        opaque_source_argument_distinction=True,
        radius_endpoint_mismatch=True,
        phase_endpoint_mismatch=True,
    )


def source_id_fixture():
    """Cancel only the proved source ID and retain unknown IDs."""
    proved = contract_trace_chain(
        [
            {"name": "restoration", "lower": "reference_end", "upper": "restoration_end"},
            {"name": "postrestore", "lower": "postrestore_start", "upper": "Rm"},
        ],
        {"postrestore_start": "restoration_end"},
    )
    expected = {"reference_end": 1, "Rm": -1}
    if proved.get("signed_coefficients") != expected:
        raise AssertionError(("signed source coefficients", proved))
    if proved.get("cancelled_interfaces") != ["restoration_end"]:
        raise AssertionError(("signed source cancellation", proved))
    unknown = contract_trace_chain(
        [{"name": "unknown", "lower": "unknown_left", "upper": "unknown_right"}],
        {},
    )
    if unknown.get("signed_coefficients") != {"unknown_left": 1, "unknown_right": -1}:
        raise AssertionError(("unknown source IDs cancelled", unknown))
    return dict(
        passed=True,
        signed_coefficients=proved["signed_coefficients"],
        cancelled_interfaces=proved["cancelled_interfaces"],
        unknown_ids_retained=unknown["signed_coefficients"],
        bounds_never_cancelled_as_functions=True,
    )


def physical_fixture():
    started = time.perf_counter()
    ctx = mp.mp.clone()
    ctx.dps = 70
    z = ctx.mpf(Z_TEXT)
    tolerance = ctx.mpf(TOL_TEXT)

    # The two native coordinates meet at the same physical point:
    # restoration y=-8+1=-7 and postrestore y=-7+0=-7.
    restoration_y = lambda native: -ctx.mpf(8) + native
    postrestore_y = lambda native: -ctx.mpf(7) + native
    outer_left, seam, outer_right = -ctx.mpf(8), -ctx.mpf(7), -ctx.mpf(6)
    pieces = []
    integral_values = []
    for name, native_to_y, left_y, right_y in (
        ("restoration", restoration_y, outer_left, seam),
        ("postrestore", postrestore_y, seam, outer_right),
    ):
        values = {}
        for ordinary_z in (False, True):
            integral = native_integral(ctx, native_to_y, z, ordinary_z=ordinary_z)
            endpoint = endpoint_orientation(ctx, left_y, right_y, z, ordinary_z=ordinary_z)
            if not close(ctx, integral, endpoint, tolerance):
                raise AssertionError(("physical endpoint orientation", name, ordinary_z, integral, endpoint))
            values["C1" if ordinary_z else "C0"] = dict(
                integral=ctx.nstr(integral, 20),
                endpoint=ctx.nstr(endpoint, 20),
            )
            integral_values.append(integral)
        pieces.append(dict(name=name, native_interval=["0", "1"], physical_interval=[str(left_y), str(right_y)], **values))

    seam_c0_left = delta_rq(ctx, seam, z)[0]
    seam_c0_right = delta_rq(ctx, seam, z)[0]
    seam_c1_left = delta_rq(ctx, seam, z, ordinary_z=True)[0]
    seam_c1_right = delta_rq(ctx, seam, z, ordinary_z=True)[0]
    if seam_c0_left != seam_c0_right or seam_c1_left != seam_c1_right:
        raise AssertionError("affine native seam did not use the same physical y")

    # Both external traces are nonzero; only the internal seam cancels.
    outer_values = {
        "left_C0": delta_rq(ctx, outer_left, z)[0],
        "right_C0": delta_rq(ctx, outer_right, z)[0],
        "left_C1": delta_rq(ctx, outer_left, z, ordinary_z=True)[0],
        "right_C1": delta_rq(ctx, outer_right, z, ordinary_z=True)[0],
    }
    if any(abs(value) <= ctx.mpf("1e-12") for value in outer_values.values()):
        raise AssertionError(("external endpoint unexpectedly zero", outer_values))

    first_c0, first_c1, second_c0, second_c1 = integral_values
    total_c0 = first_c0 + second_c0
    total_c1 = first_c1 + second_c1
    if not close(ctx, total_c0, endpoint_orientation(ctx, outer_left, outer_right, z), tolerance):
        raise AssertionError("internal C0 seam did not cancel")
    if not close(ctx, total_c1, endpoint_orientation(ctx, outer_left, outer_right, z, ordinary_z=True), tolerance):
        raise AssertionError("internal C1 seam did not cancel")

    # Negative controls: each uses the same source but the wrong physical rule.
    correct_first_c0 = first_c0
    wrong_radius = native_integral(ctx, restoration_y, z, include_radius=False)
    wrong_phase = integrate(
        ctx,
        lambda native: -2 * delta_rq(ctx, -7 + 2 * native, z)[1],
        0,
        1,
        segments=12,
    )
    correct_stretched = integrate(
        ctx,
        lambda native: -2 * delta_rq(ctx, -7 + 2 * native, z)[1],
        0,
        1,
        segments=12,
    )
    wrong_jacobian = integrate(
        ctx,
        lambda native: -delta_rq(ctx, -7 + 2 * native, z)[1],
        0,
        1,
        segments=12,
    )
    wrong_sign = -correct_first_c0
    endpoint_omission = delta_rq(ctx, outer_left, z)[0]
    controls = {
        "missing_physical_radius_gap": abs(correct_first_c0 - wrong_radius),
        "wrong_phase_gap": abs(second_c0 - wrong_phase),
        "missing_physical_jacobian_gap": abs(correct_stretched - wrong_jacobian),
        "incorrect_sign_gap": abs(correct_first_c0 - wrong_sign),
        "omitted_outer_endpoint_gap": abs(second_c0 - endpoint_omission),
    }
    for name, gap in controls.items():
        if gap <= ctx.mpf("1e-8"):
            raise AssertionError((name, gap))

    return dict(
        passed=True,
        execution_seconds=round(time.perf_counter() - started, 3),
        fixture=dict(
            finite_N=17,
            Z=Z_TEXT,
            epsilon=EPS_TEXT,
            physical_radius="2.3*exp(y+8)",
            affine_native_maps=["restoration y=-8+s", "postrestore y=-7+s"],
            seam_identity="-8+1=-7+0",
            quadrature="mpmath tanh-sinh over 12 fixed subintervals per native interval",
            source_midpoints_used=False,
            accepted_ancestor_constructors_loaded=False,
        ),
        checks=dict(
            C0_endpoint_orientation=True,
            ordinary_Z1_endpoint_orientation=True,
            ordinary_Z1_nonconstant_slow_amplitude=True,
            internal_seam_cancellation=True,
            both_external_endpoints_nonzero=True,
            physical_RQy_integrated=True,
            physical_R_included=True,
            physical_jacobian_included=True,
            total_C0=ctx.nstr(total_c0, 20),
            total_C1=ctx.nstr(total_c1, 20),
            seam_C0=ctx.nstr(seam_c0_left, 20),
            seam_C1=ctx.nstr(seam_c1_left, 20),
            amplitude_Z_derivative=ctx.nstr(amplitude_perturbation_z(ctx, z), 12),
        ),
        pieces=pieces,
        outer_endpoints={key: ctx.nstr(value, 20) for key, value in outer_values.items()},
        negative_controls={key: ctx.nstr(value, 12) for key, value in controls.items()},
        limits=dict(
            independent_physics_only=True,
            production_archives_not_loaded=True,
            full_NS_validation=False,
        ),
    )


def run():
    physical = physical_fixture()
    fingerprint = fingerprint_adversarial_fixture()
    active_fingerprint = active_fingerprint_adversarial_fixture()
    source_ids = source_id_fixture()
    return dict(
        passed=physical["passed"] and fingerprint["passed"] and active_fingerprint["passed"] and source_ids["passed"],
        physical=physical,
        fingerprint=fingerprint,
        active_fingerprint=active_fingerprint,
        source_ids=source_ids,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2), flush=True)
