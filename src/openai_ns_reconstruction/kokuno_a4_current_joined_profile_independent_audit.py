"""Implementation-distinct audit of the current Kokuno joined profile through X_h.

This Agent-4 validator is intentionally downstream of Agent-1 PR #959.  The
scientific reference path first save/reloads ``KokunoPA16CurrentJoinedProfile``
and then reads only its public ``profile_values(X, eta)`` scalar surface.
Production ``radial_derivatives`` values are treated only as quantities under
audit.

Derivatives are reconstructed in ``s = log X`` with a centered fourth-order
finite difference and converted by ``d/dX = (1/X) d/ds``.  This is independent
of Agent-1's source-form scaled-x derivative implementation.  One-sided fourth
order log-X stencils independently probe the Xi and X_h endpoints.

This is a source-coordinate prerequisite audit.  It does not construct a
Cartesian velocity, pressure, forcing or Navier--Stokes residual, and it cannot
promote ``pde_validated``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import tempfile
from typing import Any

import numpy as np

from .kokuno_pa16_current_joined_profile import KokunoPA16CurrentJoinedProfile


SCHEMA = "kokuno-a4-current-joined-profile-independent-audit-v1"
SEED = 9173651
LOG_STEPS = (4.0e-3, 2.0e-3, 1.0e-3)
OFFGRID_PER_REGION = 6
ETA_CENTER_PROBES = (0.0, -1.0e-8, 1.0e-8, -1.0e-6, 1.0e-6)
DERIVATIVE_RMS_GATE = 5.0e-3
DERIVATIVE_MAX_GATE = 2.0e-2
REFINEMENT_GATE = 8.0
REFINEMENT_FLOOR = 2.0e-9
H_CLOSURE_GATE = 5.0e-12
DERIVATIVE_CHAIN_GATE = 5.0e-5
XI_JET_GATE = 2.0e-2
XH_E_SLOPE_GATE = 2.0e-4
XH_F_SLOPE_GATE = 2.0e-4
XH_U_SLOPE_GATE = 5.0e-3
NONTRIVIALITY_FLOOR = 1.0e-14
FINAL_MOMENTUM_GATE = 1.0e-3
FINAL_DIVERGENCE_GATE = 1.0e-5

_TRUTH_BOUNDARY = {
    "current_joined_profile_independently_audited": True,
    "public_values_only_numerical_reference": True,
    "implementation_distinct_logX_operator": True,
    "physical_cartesian_axis_near_assessed": False,
    "global_cartesian_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "leading_only_ns_residual_assessed": False,
    "leading_plus_oscillatory_ns_residual_assessed": False,
    "after_correction_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
}


def _rms(values: Any) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _public_matrix(field: Any, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
    values = field.profile_values(X, eta)
    return np.stack(
        [
            np.asarray(values["F_current_joined"], dtype=float),
            np.asarray(values["U_current_joined"], dtype=float),
            np.asarray(values["E_current_joined"], dtype=float),
            np.asarray(values["H_current_joined"], dtype=float),
        ],
        axis=-1,
    )


def _production_derivatives(field: Any, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
    values = field.radial_derivatives(X, eta)
    return np.stack(
        [
            np.asarray(values["F_current_joined_X"], dtype=float),
            np.asarray(values["U_current_joined_X"], dtype=float),
            np.asarray(values["E_current_joined_X"], dtype=float),
        ],
        axis=-1,
    )


def _centered_log_fd4(field: Any, X: np.ndarray, eta: np.ndarray, h: float) -> np.ndarray:
    """Reconstruct d(F,U,E)/dX from public values only."""
    X = np.asarray(X, dtype=float).reshape(-1)
    eta = np.asarray(eta, dtype=float).reshape(-1)
    offsets = np.asarray([-2.0, -1.0, 1.0, 2.0])
    coeff = np.asarray([1.0, -8.0, 8.0, -1.0]) / 12.0
    Xs = X[:, None] * np.exp(h * offsets[None, :])
    etas = np.broadcast_to(eta[:, None], Xs.shape)
    public = _public_matrix(field, Xs, etas)[..., :3]
    d_ds = np.einsum("nkj,k->nj", public, coeff) / float(h)
    return d_ds / X[:, None]


def _one_sided_log_fd4(
    field: Any,
    X0: float,
    eta: np.ndarray,
    h: float,
    *,
    direction: str,
) -> np.ndarray:
    """Fourth-order one-sided log-X derivative from public values only."""
    eta = np.asarray(eta, dtype=float).reshape(-1)
    if direction == "forward":
        nodes = np.arange(5.0)
        coeff = np.asarray([-25.0, 48.0, -36.0, 16.0, -3.0]) / 12.0
    elif direction == "backward":
        nodes = -np.arange(5.0)
        coeff = np.asarray([25.0, -48.0, 36.0, -16.0, 3.0]) / 12.0
    else:
        raise ValueError("direction must be 'forward' or 'backward'")
    Xs = float(X0) * np.exp(float(h) * nodes[None, :])
    Xs = np.broadcast_to(Xs, (eta.size, nodes.size))
    etas = np.broadcast_to(eta[:, None], Xs.shape)
    public = _public_matrix(field, Xs, etas)[..., :3]
    d_ds = np.einsum("nkj,k->nj", public, coeff) / float(h)
    return d_ds / float(X0)


def _component_metrics(
    production: np.ndarray,
    independent: np.ndarray,
    X: np.ndarray,
    eta: np.ndarray,
    regions: np.ndarray,
) -> dict[str, Any]:
    names = ("F_X", "U_X", "E_X")
    out: dict[str, Any] = {}
    for j, name in enumerate(names):
        p = np.asarray(production[:, j], dtype=float)
        r = np.asarray(independent[:, j], dtype=float)
        error = p - r
        scale_rms = max(_rms(p), _rms(r), NONTRIVIALITY_FLOOR)
        scale_max = max(
            float(np.max(np.abs(p))),
            float(np.max(np.abs(r))),
            scale_rms,
            NONTRIVIALITY_FLOOR,
        )
        idx = int(np.argmax(np.abs(error)))
        out[name] = {
            "absolute_rms": _rms(error),
            "absolute_max": float(np.max(np.abs(error))),
            "scale_rms": scale_rms,
            "relative_rms": _rms(error) / scale_rms,
            "relative_max": float(np.max(np.abs(error))) / scale_max,
            "worst": {
                "X": float(X[idx]),
                "eta": float(eta[idx]),
                "region": str(regions[idx]),
                "production": float(p[idx]),
                "independent": float(r[idx]),
                "absolute_error": float(abs(error[idx])),
            },
        }
    return out


def _sample_outer_points(field: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(SEED)
    log_x_i = math.log(float(field.X_i)) - float(field.log_X_R)
    T_sh = float(field.certificate.selected_T_sh)
    log_x_sep = log_x_i + T_sh
    max_h = max(LOG_STEPS)
    shaping_lo = log_x_i + max(0.20 * T_sh, 4.0 * max_h)
    shaping_hi = log_x_sep - max(0.20 * T_sh, 4.0 * max_h)
    if not shaping_lo < shaping_hi:
        raise RuntimeError("insufficient angular-shaping collar for frozen log-X stencil")

    windows = {
        "angular_shaping": (shaping_lo, shaping_hi),
        "axial_restoration": (-7.85, -7.15),
        "pa16_repair": (-5.85, -5.15),
    }
    log_x_parts: list[np.ndarray] = []
    eta_parts: list[np.ndarray] = []
    region_parts: list[np.ndarray] = []
    for name, (lo, hi) in windows.items():
        # Stratification prevents a random draw from clustering near one edge.
        edges = np.linspace(lo, hi, OFFGRID_PER_REGION + 1)
        u = rng.random(OFFGRID_PER_REGION)
        log_x = edges[:-1] + u * (edges[1:] - edges[:-1])
        eta = rng.uniform(-0.70, 0.70, OFFGRID_PER_REGION)
        log_x_parts.append(log_x)
        eta_parts.append(eta)
        region_parts.append(np.asarray([name] * OFFGRID_PER_REGION, dtype=object))

    # eta-center/near-center probes are source-coordinate symmetry probes; they
    # are not relabelled as physical Cartesian axis-near validation.
    log_x_parts.append(np.full(len(ETA_CENTER_PROBES), -7.50))
    eta_parts.append(np.asarray(ETA_CENTER_PROBES, dtype=float))
    region_parts.append(np.asarray(["eta_center_restore"] * len(ETA_CENTER_PROBES), dtype=object))

    log_x = np.concatenate(log_x_parts)
    eta = np.concatenate(eta_parts)
    regions = np.concatenate(region_parts)
    X = np.exp(float(field.log_X_R) + log_x)
    if np.any(~np.isfinite(X)) or np.any(X <= float(field.X_i)) or np.any(X >= float(field.X_h)):
        raise RuntimeError("frozen outer audit points are outside (X_i,X_h)")
    return X, eta, regions


def _refinement_assessment(level_metrics: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for name in ("F_X", "U_X", "E_X"):
        e0 = float(level_metrics[0][name]["absolute_rms"])
        e1 = float(level_metrics[1][name]["absolute_rms"])
        e2 = float(level_metrics[2][name]["absolute_rms"])
        scale = max(float(level_metrics[2][name]["scale_rms"]), NONTRIVIALITY_FLOOR)
        r01 = math.inf if e1 == 0.0 else e0 / e1
        r12 = math.inf if e2 == 0.0 else e1 / e2
        normalized_change = abs(e1 - e2) / scale
        floor = normalized_change <= REFINEMENT_FLOOR
        passed = floor or min(r01, r12) >= REFINEMENT_GATE
        out[name] = {
            "coarse_to_medium_ratio": r01,
            "medium_to_fine_ratio": r12,
            "medium_to_fine_normalized_change": normalized_change,
            "numerical_floor_exception": floor,
            "pass": passed,
        }
    return out


def _endpoint_refinement(sequence: list[np.ndarray], value_scale: np.ndarray) -> dict[str, Any]:
    """Assess convergence of an endpoint derivative sequence componentwise."""
    out: dict[str, Any] = {}
    names = ("F_X", "U_X", "E_X")
    for j, name in enumerate(names):
        d01 = _rms(sequence[0][:, j] - sequence[1][:, j])
        d12 = _rms(sequence[1][:, j] - sequence[2][:, j])
        scale = max(_rms(sequence[2][:, j]), _rms(value_scale[:, j]), NONTRIVIALITY_FLOOR)
        ratio = math.inf if d12 == 0.0 else d01 / d12
        normalized_change = d12 / scale
        floor = normalized_change <= REFINEMENT_FLOOR
        out[name] = {
            "coarse_medium_change": d01,
            "medium_fine_change": d12,
            "change_ratio": ratio,
            "medium_fine_normalized_change": normalized_change,
            "numerical_floor_exception": floor,
            "pass": bool(floor or ratio >= REFINEMENT_GATE),
        }
    return out


def audit_rebound_public_profile(field: Any) -> dict[str, Any]:
    """Run the frozen scientific audit on an already rebound public field."""
    if not bool(field.route_ready):
        return {
            "schema": SCHEMA,
            "status": "blocked_by_current_tsh_geometry",
            "audit_pass": False,
            "route_ready": False,
            "protocol": protocol(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    X, eta, regions = _sample_outer_points(field)
    production = _production_derivatives(field, X, eta)
    level_refs: list[np.ndarray] = []
    level_metrics: list[dict[str, Any]] = []
    for h in LOG_STEPS:
        independent = _centered_log_fd4(field, X, eta, h)
        level_refs.append(independent)
        level_metrics.append(_component_metrics(production, independent, X, eta, regions))

    fine = level_refs[-1]
    public = _public_matrix(field, X, eta)
    F, U, E, H = (public[:, k] for k in range(4))
    root = np.sqrt(2.0 * X)
    h_scale = max(_rms(H), _rms(root * E), NONTRIVIALITY_FLOOR)
    H_closure = _rms(H - root * E) / h_scale
    chain_target = fine[:, 2] / root - E / root**3
    chain_scale = max(_rms(fine[:, 0]), _rms(chain_target), NONTRIVIALITY_FLOOR)
    derivative_chain = _rms(fine[:, 0] - chain_target) / chain_scale

    refinement = _refinement_assessment(level_metrics)
    finest_pass = all(
        level_metrics[-1][name]["relative_rms"] <= DERIVATIVE_RMS_GATE
        and level_metrics[-1][name]["relative_max"] <= DERIVATIVE_MAX_GATE
        for name in ("F_X", "U_X", "E_X")
    )
    refinement_pass = all(refinement[name]["pass"] for name in refinement)

    endpoint_eta = np.asarray([-0.60, -0.25, 0.0, 0.25, 0.60], dtype=float)
    xi_left: list[np.ndarray] = []
    xi_right: list[np.ndarray] = []
    xh_back: list[np.ndarray] = []
    for h in LOG_STEPS:
        xi_left.append(_one_sided_log_fd4(field, float(field.X_i), endpoint_eta, h, direction="backward"))
        xi_right.append(_one_sided_log_fd4(field, float(field.X_i), endpoint_eta, h, direction="forward"))
        xh_back.append(_one_sided_log_fd4(field, float(field.X_h), endpoint_eta, h, direction="backward"))

    xi_values = _public_matrix(
        field,
        np.full(endpoint_eta.size, float(field.X_i)),
        endpoint_eta,
    )[:, :3]
    xi_value_derivative_scale = np.abs(xi_values) / float(field.X_i)
    xi_component_mismatch: dict[str, float] = {}
    for j, name in enumerate(("F_X", "U_X", "E_X")):
        scale = max(
            _rms(xi_left[-1][:, j]),
            _rms(xi_right[-1][:, j]),
            _rms(xi_value_derivative_scale[:, j]),
            NONTRIVIALITY_FLOOR,
        )
        xi_component_mismatch[name] = float(
            np.max(np.abs(xi_left[-1][:, j] - xi_right[-1][:, j])) / scale
        )
    xi_mismatch_max = max(xi_component_mismatch.values())
    xi_ref_left = _endpoint_refinement(xi_left, xi_value_derivative_scale)
    xi_ref_right = _endpoint_refinement(xi_right, xi_value_derivative_scale)
    xi_refinement_pass = all(v["pass"] for v in xi_ref_left.values()) and all(
        v["pass"] for v in xi_ref_right.values()
    )

    xh_values = _public_matrix(
        field,
        np.full(endpoint_eta.size, float(field.X_h)),
        endpoint_eta,
    )[:, :3]
    xh_scale = np.abs(xh_values) / float(field.X_h)
    xh_refinement = _endpoint_refinement(xh_back, xh_scale)
    xh_refinement_pass = all(v["pass"] for v in xh_refinement.values())
    xh_fine = xh_back[-1]
    Fh, Uh, Eh = (xh_values[:, k] for k in range(3))
    xh_E_slope_error = float(np.max(np.abs(float(field.X_h) * xh_fine[:, 2] / Eh - 0.1)))
    xh_F_slope_error = float(np.max(np.abs(float(field.X_h) * xh_fine[:, 0] / Fh + 0.4)))
    u_scale = max(_rms(Uh), 1.0)
    xh_U_slope_normalized = float(np.max(np.abs(float(field.X_h) * xh_fine[:, 1])) / u_scale)

    nontrivial = bool(
        max(_rms(F), _rms(U), _rms(E)) > NONTRIVIALITY_FLOOR
        and max(_rms(fine[:, 0]), _rms(fine[:, 1]), _rms(fine[:, 2])) > NONTRIVIALITY_FLOOR
    )
    endpoint_pass = bool(
        xi_mismatch_max <= XI_JET_GATE
        and xi_refinement_pass
        and xh_refinement_pass
        and xh_E_slope_error <= XH_E_SLOPE_GATE
        and xh_F_slope_error <= XH_F_SLOPE_GATE
        and xh_U_slope_normalized <= XH_U_SLOPE_GATE
    )
    audit_pass = bool(
        finest_pass
        and refinement_pass
        and H_closure <= H_CLOSURE_GATE
        and derivative_chain <= DERIVATIVE_CHAIN_GATE
        and endpoint_pass
        and nontrivial
    )

    return {
        "schema": SCHEMA,
        "status": "pass" if audit_pass else "fail",
        "audit_pass": audit_pass,
        "route_ready": True,
        "protocol": protocol(),
        "samples": {
            "count": int(X.size),
            "X_min": float(np.min(X)),
            "X_max": float(np.max(X)),
            "eta_min": float(np.min(eta)),
            "eta_max": float(np.max(eta)),
            "regions": {name: int(np.sum(regions == name)) for name in np.unique(regions)},
        },
        "derivative_levels": [
            {"log_step": h, "metrics": metrics}
            for h, metrics in zip(LOG_STEPS, level_metrics)
        ],
        "refinement": refinement,
        "algebraic": {
            "H_relative_rms_closure": H_closure,
            "independent_F_E_derivative_chain_relative_rms": derivative_chain,
        },
        "endpoints": {
            "eta": endpoint_eta.tolist(),
            "Xi": {
                "left_right_component_scale_normalized_max": xi_component_mismatch,
                "left_right_scale_normalized_max": xi_mismatch_max,
                "left_refinement": xi_ref_left,
                "right_refinement": xi_ref_right,
            },
            "X_h": {
                "E_dimensionless_slope_max_abs_error": xh_E_slope_error,
                "F_dimensionless_slope_max_abs_error": xh_F_slope_error,
                "U_dimensionless_slope_scale_normalized_max": xh_U_slope_normalized,
                "refinement": xh_refinement,
            },
        },
        "nontrivial": nontrivial,
        "gates": {
            "finest_derivative_pass": finest_pass,
            "refinement_pass": refinement_pass,
            "H_closure_pass": H_closure <= H_CLOSURE_GATE,
            "derivative_chain_pass": derivative_chain <= DERIVATIVE_CHAIN_GATE,
            "endpoint_pass": endpoint_pass,
            "nontriviality_pass": nontrivial,
        },
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }


def protocol() -> dict[str, Any]:
    return {
        "seed": SEED,
        "offgrid_per_region": OFFGRID_PER_REGION,
        "eta_center_probes": list(ETA_CENTER_PROBES),
        "independent_operator": "centered-FD4-in-log-X; one-sided-FD4-at-Xi/X_h",
        "log_step_ladder": list(LOG_STEPS),
        "derivative_relative_rms_gate": DERIVATIVE_RMS_GATE,
        "derivative_relative_max_gate": DERIVATIVE_MAX_GATE,
        "refinement_ratio_gate": REFINEMENT_GATE,
        "refinement_normalized_change_floor": REFINEMENT_FLOOR,
        "H_closure_gate": H_CLOSURE_GATE,
        "derivative_chain_gate": DERIVATIVE_CHAIN_GATE,
        "Xi_left_right_jet_gate": XI_JET_GATE,
        "Xh_E_slope_gate": XH_E_SLOPE_GATE,
        "Xh_F_slope_gate": XH_F_SLOPE_GATE,
        "Xh_U_slope_gate": XH_U_SLOPE_GATE,
        "final_momentum_gate_unassessed": FINAL_MOMENTUM_GATE,
        "final_divergence_gate_unassessed": FINAL_DIVERGENCE_GATE,
        "residual_defined_forcing_forbidden": True,
    }


def materialize_current_joined_profile_independent_audit(
    field: KokunoPA16CurrentJoinedProfile | None = None,
) -> dict[str, Any]:
    """Save/reload the current A1 artifact and run the frozen public-surface audit."""
    original = KokunoPA16CurrentJoinedProfile() if field is None else field
    if not isinstance(original, KokunoPA16CurrentJoinedProfile):
        raise TypeError("field must be KokunoPA16CurrentJoinedProfile")
    with tempfile.TemporaryDirectory(prefix="kokuno-a4-joined-") as tmp:
        path = Path(tmp) / "current_joined_profile.json"
        original.save_configuration(path)
        rebound = KokunoPA16CurrentJoinedProfile.load_configuration(path)
        if rebound.semantic_sha256 != original.semantic_sha256:
            raise RuntimeError("save/reload semantic identity changed")
        receipt = audit_rebound_public_profile(rebound)
        receipt["upstream_semantic_sha256"] = original.semantic_sha256
        receipt["rebound_semantic_sha256"] = rebound.semantic_sha256
        receipt["save_reload_identity_pass"] = True
        return receipt


def enforce_current_joined_profile_independent_audit(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA:
        raise RuntimeError("unexpected A4 joined-profile audit schema")
    if not receipt.get("audit_pass", False):
        raise RuntimeError(f"A4 joined-profile audit did not pass: {receipt.get('status')}")
    if receipt.get("truth_boundary", {}).get("pde_validated") is not False:
        raise RuntimeError("joined-profile audit cannot promote pde_validated")


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--enforce", action="store_true")
    args = parser.parse_args()
    receipt = materialize_current_joined_profile_independent_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False))
    if args.enforce:
        enforce_current_joined_profile_independent_audit(receipt)


if __name__ == "__main__":
    _main()
