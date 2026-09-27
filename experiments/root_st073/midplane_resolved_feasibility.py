"""Bounded moment-only feasibility diagnostic for the resolved midplane repair.

The earlier resolved joint SLSQP run imposed split-quadrature moment equations
and sampled cone inequalities simultaneously.  This experiment keeps the
same resolved split quadrature and direct finite-difference delta jets, but
fits the moment equations alone.  Since affine jets make momentum quadratic in
the mode amplitudes, the moment map is cached as a quadratic polynomial and
the least-squares Jacobian is evaluated analytically.

This is a numerical feasibility diagnostic.  A stopped optimizer, or a
moment-only candidate, is not an infeasibility proof and does not establish
the scale recursion.
"""

import argparse
import json
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares

from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from affine_momentum import jets, momentum
from joined_field import coordinates
from midplane_physical_covariance_pairs import support as physical_support
from radial_continuation import ROOT
from separated_moment_modes import RADIAL_WINDOWS_THREE, SeparatedMomentModes


KNOTS = (11.0, 15.0, 19.0)
RADIAL_BREAKS = tuple(sorted({v for lo, hi in RADIAL_WINDOWS_THREE
                              for v in (lo, (lo + hi) / 2.0, hi)}))


class ZeroBackground:
    """Metadata-compatible zero field for exact affine basis jets."""

    def __init__(self, base):
        for name in ("inner", "nu", "join_X", "ratio"):
            setattr(self, name, getattr(base, name))

    def fields(self, points, tau):
        points = np.asarray(points, float)
        return np.zeros_like(points), np.zeros(len(points))


def _integrated_moment_coefficients(slice_data):
    """Return offset, linear, quadratic coefficients for one scale.

    If ``a`` is the delta amplitude vector, the returned arrays satisfy
    ``outer_moments(slice_data, a) == offset + linear @ a +
    einsum('mij,i,j->m', quadratic, a, a)`` up to roundoff.
    """
    u0, g0, l0 = slice_data["baseline"]
    mu, mg, ml = slice_data["modes"]
    mode_count = mu.shape[0]

    # Momentum is l + grad @ u.  Symmetrising the quadratic coefficient keeps
    # the polynomial exact while making the analytic derivative transparent.
    r0 = momentum((u0, g0, l0))
    linear_r = ml + np.einsum("pnij,nj->pni", mg, u0)
    linear_r += np.einsum("nij,pnj->pni", g0, mu)
    quadratic_r = np.einsum("pnij,qnj->pqni", mg, mu)
    quadratic_r = 0.5 * (quadratic_r + quadratic_r.swapaxes(0, 1))

    offsets = []
    linears = []
    quadratics = []
    for panel in slice_data["panels"]:
        sl = panel["slice"]
        rr = panel["radii"]
        ww = panel["weights"]
        ro = panel["outer_radius"]
        theta_weights = ww * rr**2 / ro**2
        axial_weights = ww * rr / ro
        offsets.append(np.array([
            -np.dot(theta_weights, r0[sl, 1]),
            -np.dot(axial_weights, r0[sl, 2]),
        ]))
        linears.append(np.stack([
            -np.einsum("n,pn->p", theta_weights, linear_r[:, sl, 1]),
            -np.einsum("n,pn->p", axial_weights, linear_r[:, sl, 2]),
        ], axis=1))
        q = np.empty((mode_count, mode_count, 2), float)
        q[:, :, 0] = -np.einsum(
            "n,pqn->pq", theta_weights, quadratic_r[:, :, sl, 1])
        q[:, :, 1] = -np.einsum(
            "n,pqn->pq", axial_weights, quadratic_r[:, :, sl, 2])
        quadratics.append(q)
    return (np.concatenate(offsets), np.concatenate(linears, axis=1).T,
            np.moveaxis(np.concatenate(quadratics, axis=2), 2, 0))


def _assemble_coefficients(moment_data):
    pieces = [_integrated_moment_coefficients(item) for item in moment_data]
    return (np.concatenate([piece[0] for piece in pieces]),
            np.vstack([piece[1] for piece in pieces]),
            np.concatenate([piece[2] for piece in pieces], axis=0))


def _evaluate(coefficients, amplitudes):
    offset, linear, quadratic = coefficients
    amplitudes = np.asarray(amplitudes, float)
    return (offset + linear @ amplitudes
            + np.einsum("mij,i,j->m", quadratic, amplitudes, amplitudes))


def _jacobian(coefficients, amplitudes):
    _, linear, quadratic = coefficients
    amplitudes = np.asarray(amplitudes, float)
    return (linear
            + np.einsum("mij,j->mi", quadratic, amplitudes)
            + np.einsum("mji,j->mi", quadratic, amplitudes))


def _make_sources(inner, supports, k):
    if k in (11.0, 19.0):
        source = json.loads((ROOT / "compact_potential" /
                             f"midplane_wave_source_k{int(k)}.json").read_text())
        support = next(row["support"] for row in supports["scales"]
                       if int(row["k"]) == int(k))
    else:
        tau = 0.5 * 2.0**(-k)
        X = inner.p.X_max * (1.0 + 15.0 * 0.325)**2
        source = dict(tau=tau,
                      point=inner.from_similarity([X], [-0.0125], tau)[0])
        support = physical_support(inner, source, inner.nu)
    return source, support


def _cone_replay(inner, field, supports, k, offsets=(-0.99, 0.0, 0.99),
                 order=48):
    """Evaluate the sampled physical cone directly from candidate jets.

    The radial primitive uses a fresh split Gauss rule and one field jet per
    quadrature point.  This evaluates the complete candidate field independently of the
    earlier support-local four-mode cone cache.
    """
    source, support = _make_sources(inner, supports, k)
    tau = float(source["tau"])
    radius, _, height = map(float, source["point"])
    centers = [(radius + x * support["radial_halfwidth"], 0.0,
                height + z * 1.4 * support["axial_halfwidth"])
               for x in offsets for z in offsets]
    points = list(centers)
    panels = []
    nodes, weights = leggauss(order)
    for r, _, z in centers:
        q = float(coordinates(0.0, z / np.sqrt(inner.nu), tau,
                              inner.h)["q"])
        ri = np.sqrt(2.0 * inner.nu * q * inner.p.X_max)
        cuts = sorted(set(np.clip([0.0, ri, r], 0.0, r)))
        start = len(points)
        radii = []
        radial_weights = []
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            rr = (hi + lo) / 2.0 + (hi - lo) / 2.0 * nodes
            ww = (hi - lo) / 2.0 * weights
            points.extend((float(x), 0.0, float(z)) for x in rr)
            radii.extend(rr)
            radial_weights.extend(ww)
        panels.append((slice(start, len(points)), np.asarray(radii),
                       np.asarray(radial_weights), float(r)))
    points = np.asarray(points)
    args = (points, tau, 0.0005 * np.sqrt(inner.nu * tau), 0.0001 * tau)
    jet = jets(field, *args)
    residual = momentum(jet)
    rows = []
    for index, (panel, rr, ww, r) in enumerate(panels):
        u = jet[0][index]
        gradient = jet[1][index]
        F = float(u[1] / r)
        shear = np.array([gradient[1, 0] - F, gradient[2, 0]])
        shear_norm = float(np.linalg.norm(shear))
        row = dict(index=index, radial_offset=float(offsets[index // 3]),
                   axial_offset=float(offsets[index % 3]),
                   momentum_norm=float(np.linalg.norm(residual[index])),
                   lambda_squared=None, target_dot_N=None,
                   target_dot_K=None, cone_ratio=None, cone_pass=False)
        if shear_norm == 0.0:
            rows.append(row)
            continue
        normal = shear / shear_norm
        tangent = np.array([-normal[1], normal[0]])
        lambda_squared = float(-2.0 * F * normal[0]
                               * (2.0 * F * normal[0] + shear_norm))
        row["lambda_squared"] = lambda_squared
        local = residual[panel]
        target = np.array([
            -np.dot(ww * rr**2, local[:, 1]) / r**2,
            -np.dot(ww * rr, local[:, 2]) / r,
        ])
        dot_n = float(target @ normal)
        dot_k = float(target @ tangent)
        row["target_dot_N"] = dot_n
        row["target_dot_K"] = dot_k
        denominator = abs(2.0 * F * normal[0] * dot_n)
        if lambda_squared > 0.0 and denominator > 1e-30:
            ratio = float(np.sqrt(lambda_squared) * abs(dot_k) / denominator)
            row["cone_ratio"] = ratio
            row["cone_pass"] = bool(dot_n < 0.0 and ratio < 1.0)
        rows.append(row)
    return dict(k=int(k), tau=tau, order=int(order), rows=rows,
                pass_count=int(sum(row["cone_pass"] for row in rows)),
                max_ratio=float(max((row["cone_ratio"] for row in rows
                                     if row["cone_ratio"] is not None),
                                    default=float("inf"))))


def _direct_replay(inner, base, field, n):
    """Recompute moments at a higher Gauss order using the field directly."""
    # Passing the candidate as both baseline and a zero-difference unit field
    # reuses the repository's split-panel construction without the fitted
    # quadratic cache.  The resulting values are direct field jets.
    slices = moment_slices(
        inner, base, field, orders=KNOTS, n=int(n),
        unit_fields=[field], unit_fields_are_deltas=False,
        radial_breaks=RADIAL_BREAKS)
    return np.concatenate([outer_moments(item, np.zeros(1))
                           for item in slices])


def run(fit_max_nfev=180, cone_order=48,
        replay_orders=(48, 96)):
    started = time.perf_counter()
    inner, fields = build_fields()
    start_amplitudes = np.asarray(json.loads(
        (ROOT / "midplane_axial_cone_all_knots_repair.json").read_text()
    )["amplitudes"], float)
    base_field = fields["two_sided_cone"]

    def make_field(amplitudes):
        return SeparatedMomentModes(
            base_field, amplitudes, windows=RADIAL_WINDOWS_THREE,
            knots=KNOTS)

    zero = ZeroBackground(base_field)
    eye = np.eye(len(start_amplitudes))
    unit_fields = [SeparatedMomentModes(
        zero, eye[index], windows=RADIAL_WINDOWS_THREE, knots=KNOTS)
                    for index in range(len(start_amplitudes))]
    current = make_field(start_amplitudes)

    # This is the resolved split fit grid used by the failed joint repair.
    moment_data = moment_slices(
        inner, base_field, current, orders=KNOTS, unit_fields=unit_fields,
        unit_fields_are_deltas=True, radial_breaks=RADIAL_BREAKS, n=24)
    coefficients = _assemble_coefficients(moment_data)
    probe = np.random.default_rng(73).normal(size=len(eye)) * 0.01
    reference = np.concatenate([outer_moments(item, probe) for item in moment_data])
    np.testing.assert_allclose(_evaluate(coefficients, probe), reference, rtol=1e-10, atol=1e-7)
    print("Resolved quadratic moment cache verified", flush=True)
    initial_moments = _evaluate(coefficients, np.zeros(len(eye)))
    mscale = np.maximum(np.max(np.abs(np.asarray([
        _evaluate(coefficients, row) - initial_moments for row in eye
    ])), axis=0), 1e-8)

    def normalized(delta):
        return _evaluate(coefficients, delta) / mscale

    def normalized_jacobian(delta):
        return _jacobian(coefficients, delta) / mscale[:, None]

    initial_jacobian = normalized_jacobian(np.zeros(len(eye)))
    singular_values = np.linalg.svd(initial_jacobian, compute_uv=False)
    rank_cutoff = float(singular_values[0] * 1e-10) if len(singular_values) else 0.0
    numerical_rank = int(np.sum(singular_values > rank_cutoff))
    rng = np.random.default_rng(73073)
    direction = rng.normal(size=len(eye))
    direction /= np.linalg.norm(direction)
    jac_step = 1e-6
    finite_jacobian_error = float(np.max(np.abs(
        (normalized(jac_step * direction)
         - normalized(-jac_step * direction)) / (2.0 * jac_step)
        - initial_jacobian @ direction)))

    fit = least_squares(
        normalized, np.zeros(len(eye)), jac=normalized_jacobian,
        bounds=(-40.0, 40.0), max_nfev=int(fit_max_nfev),
        ftol=1e-11, xtol=1e-11, gtol=1e-11, x_scale="jac")
    correction = np.asarray(fit.x)
    fitted_moments = _evaluate(coefficients, correction)
    fitted_amplitudes = start_amplitudes + correction
    fitted_field = make_field(fitted_amplitudes)

    # Independent higher-order replay is deliberately outside the quadratic
    # cache and uses fresh field evaluations at each requested order.
    replay = {}
    for order in replay_orders:
        values = _direct_replay(inner, base_field, fitted_field, int(order))
        replay[str(int(order))] = dict(
            moments=values.tolist(),
            normalized_max=float(np.max(np.abs(values / mscale))),
            absolute_max=float(np.max(np.abs(values))),
        )
    if len(replay_orders) >= 2:
        first = replay[str(int(replay_orders[0]))]["moments"]
        second = replay[str(int(replay_orders[1]))]["moments"]
        replay["cross_order_max_difference"] = float(
            np.max(np.abs(np.asarray(first) - np.asarray(second))))

    supports = json.loads((ROOT / "midplane_physical_covariance_pairs.json").read_text())
    cone = {}
    for label, field in (("initial", current), ("moment_only", fitted_field)):
        rows = [_cone_replay(inner, field, supports, k, order=cone_order)
                for k in KNOTS]
        cone[label] = dict(
            scales=rows,
            pass_counts={str(row["k"]): row["pass_count"] for row in rows},
            max_ratios={str(row["k"]): row["max_ratio"] for row in rows},
        )

    prior = json.loads((ROOT / "midplane_connected_cone_resolved_repair.json").read_text())
    fit_normalized_max = float(np.max(np.abs(fitted_moments / mscale)))
    next_action = (
        "Use the moment-only candidate as the initializer for a separate cone "
        "feasibility solve, first identifying the sampled rows that still fail "
        "after the independent replay; joint SLSQP failure is not evidence of "
        "moment infeasibility."
        if fit_normalized_max <= 1e-6 else
        "Inspect conditioning, scaling and bounds before deciding whether to redesign "
        "the resolved moment basis; this moment-only run did not reach numerical feasibility "
        "residual, and neither optimizer outcome proves infeasibility."
    )
    report = dict(
        source=("Resolved split-quadrature direct-delta-jet diagnostic with "
                "quadratic moment cache and moment-only least squares"),
        fit_quadrature_order=24,
        replay_quadrature_orders=[int(order) for order in replay_orders],
        radial_breaks=list(RADIAL_BREAKS),
        mode_count=int(len(start_amplitudes)),
        moment_count=int(len(initial_moments)),
        quadratic_cache=dict(
            polynomial_identity="moment(delta)=offset+linear@delta+delta^T quadratic delta",
            coefficient_shapes=[list(np.shape(item)) for item in coefficients],
            finite_difference_directional_error=finite_jacobian_error,
        ),
        initial=dict(
            amplitudes=start_amplitudes.tolist(),
            moment_values=initial_moments.tolist(),
            normalized_moment_max=float(np.max(np.abs(initial_moments / mscale))),
            absolute_moment_max=float(np.max(np.abs(initial_moments))),
            moment_scales=mscale.tolist(),
            jacobian_singular_values=singular_values.tolist(),
            jacobian_rank_cutoff=rank_cutoff,
            jacobian_numerical_rank=numerical_rank,
            jacobian_condition_number=(
                float(singular_values[0] / singular_values[-1])
                if singular_values[-1] > 0 else float("inf")),
            jacobian_row_norms=np.linalg.norm(initial_jacobian, axis=1).tolist(),
            jacobian_column_norms=np.linalg.norm(initial_jacobian, axis=0).tolist(),
        ),
        moment_only=dict(
            optimizer_success=bool(fit.success), message=str(fit.message),
            status=int(fit.status), nfev=int(fit.nfev), njev=int(fit.njev),
            cost=float(fit.cost), optimality=float(fit.optimality),
            correction=correction.tolist(), amplitudes=fitted_amplitudes.tolist(),
            correction_l2=float(np.linalg.norm(correction)),
            max_abs_correction=float(np.max(np.abs(correction))),
            normalized_moment_max=fit_normalized_max,
            absolute_moment_max=float(np.max(np.abs(fitted_moments))),
            moments=fitted_moments.tolist(),
            jacobian_singular_values=_jacobian_singular_values(
                coefficients, correction, mscale),
        ),
        independent_higher_order_replay=replay,
        sampled_cone_replay=dict(order=int(cone_order), **cone),
        prior_resolved_joint=dict(
            optimizer_success=bool(prior["optimizer_success"]),
            message=prior["message"],
            iterations=int(prior["iterations"]),
            normalized_moment_max=float(prior["normalized_moment_max"]),
            minimum_constraint=float(prior["minimum_constraint"]),
        ),
        interpretation=(
            "The previous resolved joint SLSQP iteration-limit stop is an "
            "algorithmic stopping event. It does not prove that the joint "
            "moment/cone constraints are infeasible. The moment-only result "
            "tests only the resolved sampled moment equations; cone rows remain "
            "a finite physical analogue and do not establish continuous support, "
            "a PDE solution, or scale recursion."),
        next_action=next_action,
        elapsed_seconds=float(time.perf_counter() - started),
        accepted=False,
        scale_recursion_established=False,
    )
    output = ROOT / "midplane_resolved_feasibility.json"
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode())
    print(json.dumps({
        "output": str(output),
        "optimizer_success": report["moment_only"]["optimizer_success"],
        "normalized_moment_max": report["moment_only"]["normalized_moment_max"],
        "jacobian_rank": report["initial"]["jacobian_numerical_rank"],
        "jacobian_singular_values": report["initial"]["jacobian_singular_values"],
        "replay": {key: value.get("normalized_max")
                    for key, value in replay.items() if isinstance(value, dict)},
        "cone_pass_counts": report["sampled_cone_replay"]["moment_only"]["pass_counts"],
        "elapsed_seconds": report["elapsed_seconds"],
    }, indent=2), flush=True)
    return report


def _jacobian_singular_values(coefficients, amplitudes, scales):
    values = _jacobian(coefficients, amplitudes) / np.asarray(scales)[:, None]
    return np.linalg.svd(values, compute_uv=False).tolist()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fit-max-nfev", type=int, default=180)
    parser.add_argument("--cone-order", type=int, default=48)
    parser.add_argument("--replay-orders", type=int, nargs="+", default=[48, 96])
    args = parser.parse_args()
    run(fit_max_nfev=args.fit_max_nfev, cone_order=args.cone_order,
        replay_orders=tuple(args.replay_orders))
