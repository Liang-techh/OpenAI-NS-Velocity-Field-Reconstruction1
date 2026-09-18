"""Agent-7 render-space capacity transfer of one axial Piola degree onto ST051-B.

This module consumes the immutable 33^3 ST051-B `.025` redistribution grid used
by CR-A9-053/054. It transfers exactly the divergence-preserving axial
coordinate map previously screened on ST048-S in Agent-7 PR #407:

    h_beta(z) = z * (1 - beta * (1 - (z/2)^2)^4),  |z| < 2,

with the contravariant Piola pullback

    u_beta = (h'(z) u_x(x,y,h(z)), h'(z) u_y(x,y,h(z)), u_z(x,y,h(z))).

The only tested nonzero value is the inherited beta=.075. A positive common
scale restores the sampled t=.25 grid kinetic energy. The result is a
visualization/expression-capacity diagnostic only: pressure/forcing/PDE
receipts do not transfer through this velocity transform.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from openai_ns_reconstruction import constrained_st051b_streamline_vorticity_render as render

TASK_ID = "CR003-ST051B-PIOLA-TRANSFER-075"
CLAIM_ISSUE = 498
SOURCE_PIOLA_PR = 407
SOURCE_PIOLA_HEAD = "2d7cb2fe05696c6303e09ad20b91c61aad5f977e"
SOURCE_RENDER_PR = 490
SOURCE_RENDER_HEAD = "09682c9a4c70148113639220387344ea6c9b9538"
BETA = 0.075
PHYSICAL_Z_CUT = 2.0
COARSE_STRIDE = 2

CRITERIA = {
    "minimum_axial_span_gain": 0.005,
    "minimum_mean_turn_gain": 0.0,
    "maximum_radial_span_gain": 0.05,
    "minimum_response_rank": 2,
    "maximum_response_condition": 5.0,
    "maximum_abs_response_cosine": 0.9,
}

TRUTH = {
    "canonical_velocity_changed": False,
    "production_beta_selected": False,
    "production_candidate_selected": False,
    "pressure_or_force_rebuilt": False,
    "held_out_pde_residual_evaluated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def warp_z_and_jacobian(z: np.ndarray, beta: float = BETA) -> tuple[np.ndarray, np.ndarray]:
    beta = float(beta)
    if not np.isfinite(beta) or not 0.0 <= beta <= 0.20:
        raise ValueError("beta must lie in [0,0.20]")
    z = np.asarray(z, dtype=np.float64)
    if not np.isfinite(z).all():
        raise ValueError("z must be finite")
    s = z / PHYSICAL_Z_CUT
    inside = np.abs(s) < 1.0
    one_minus = np.maximum(1.0 - s * s, 0.0)
    bump = one_minus**4
    mapped = np.where(inside, z * (1.0 - beta * bump), z)
    jac_inside = 1.0 - beta * bump + 8.0 * beta * s * s * one_minus**3
    jac = np.where(inside, jac_inside, 1.0)
    if not np.isfinite(mapped).all() or not np.isfinite(jac).all() or np.any(jac <= 0.0):
        raise ValueError("Piola map lost orientation")
    return mapped, jac


def sampled_kinetic_energy(axis: np.ndarray, field: np.ndarray) -> float:
    axis = np.asarray(axis, dtype=np.float64)
    field = np.asarray(field, dtype=np.float64)
    if field.shape != (len(axis), len(axis), len(axis), 3):
        raise ValueError("field shape mismatch")
    if len(axis) < 3 or not np.allclose(np.diff(axis), axis[1] - axis[0]):
        raise ValueError("axis must be uniform")
    h = float(axis[1] - axis[0])
    weights = np.ones(len(axis), dtype=np.float64)
    weights[[0, -1]] = 0.5
    volume_weights = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]
    return float(0.5 * h**3 * np.sum(volume_weights * np.sum(field * field, axis=-1)))


def apply_piola_grid(axis: np.ndarray, fields: np.ndarray, *, beta: float = BETA, restore_reference_energy: bool = True) -> tuple[np.ndarray, float, dict]:
    axis = np.asarray(axis, dtype=np.float64)
    fields = np.asarray(fields, dtype=np.float64)
    expected = (len(render.REFERENCE_TIMES), len(axis), len(axis), len(axis), 3)
    if fields.shape != expected or not np.isfinite(fields).all():
        raise ValueError("unexpected source grid")
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    mapped_z, jac = warp_z_and_jacobian(points[:, 2], beta)
    mapped = points.copy()
    mapped[:, 2] = mapped_z
    out = []
    for field in fields:
        evaluate = render.grid_interpolator(axis, field)
        values = np.asarray(evaluate(mapped), dtype=np.float64)
        values[:, :2] *= jac[:, None]
        out.append(values.reshape(len(axis), len(axis), len(axis), 3))
    out = np.asarray(out, dtype=np.float64)
    parent_energy = sampled_kinetic_energy(axis, fields[0])
    raw_energy = sampled_kinetic_energy(axis, out[0])
    scale = 1.0
    if restore_reference_energy:
        if not (parent_energy > 0.0 and raw_energy > 0.0):
            raise ValueError("sampled energy must stay positive")
        scale = float(np.sqrt(parent_energy / raw_energy))
        out *= scale
    dense_z = np.linspace(-PHYSICAL_Z_CUT, PHYSICAL_Z_CUT, 4001)
    dense_mapped, dense_jac = warp_z_and_jacobian(dense_z, beta)
    structure = {
        "minimum_coordinate_jacobian": float(np.min(dense_jac)),
        "maximum_coordinate_jacobian": float(np.max(dense_jac)),
        "map_monotone": bool(np.all(np.diff(dense_mapped) > 0.0)),
        "endpoint_map_error": float(max(abs(dense_mapped[0] + 2.0), abs(dense_mapped[-1] - 2.0))),
        "reference_energy": parent_energy,
        "raw_warped_reference_energy": raw_energy,
        "common_reference_energy_scale": scale,
        "scaled_reference_energy": sampled_kinetic_energy(axis, out[0]),
    }
    return out, scale, structure


def _metrics(axis: np.ndarray, field: np.ndarray) -> dict:
    evaluate = render.grid_interpolator(axis, field)
    lines = render.integrate_streamlines(evaluate, render.seed_points())
    streamline = render.streamline_metrics(lines)
    spacing = float(axis[1] - axis[0])
    _omega, magnitude = render.vorticity(field, spacing)
    vort = render._vorticity_metrics(magnitude, axis)
    vort["full_grid_rms"] = float(np.sqrt(np.mean(magnitude * magnitude)))
    return {"streamlines": streamline, "vorticity": vort}


def _relative(child: float, parent: float) -> float:
    return float(child / max(abs(parent), 1.0e-300) - 1.0)


def morphology_table(axis: np.ndarray, baseline: np.ndarray, warped: np.ndarray) -> list[dict]:
    rows = []
    for index, time in enumerate(render.REFERENCE_TIMES):
        base = _metrics(axis, baseline[index])
        child = _metrics(axis, warped[index])
        rows.append({
            "time": float(time),
            "baseline": base,
            "warped": child,
            "relative": {
                "mean_turns": _relative(child["streamlines"]["mean_absolute_turns"], base["streamlines"]["mean_absolute_turns"]),
                "max_turns": _relative(child["streamlines"]["max_absolute_turns"], base["streamlines"]["max_absolute_turns"]),
                "mean_axial_span": _relative(child["streamlines"]["mean_axial_span"], base["streamlines"]["mean_axial_span"]),
                "mean_radial_span": _relative(child["streamlines"]["mean_radial_span"], base["streamlines"]["mean_radial_span"]),
                "selected_axial_rms": _relative(child["vorticity"]["selected_axial_rms"], base["vorticity"]["selected_axial_rms"]),
                "selected_radial_rms": _relative(child["vorticity"]["selected_radial_rms"], base["vorticity"]["selected_radial_rms"]),
                "full_grid_vorticity_rms": _relative(child["vorticity"]["full_grid_rms"], base["vorticity"]["full_grid_rms"]),
            },
        })
    return rows


def _region_masks(axis: np.ndarray) -> dict[str, np.ndarray]:
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    az = np.abs(zz)
    return {
        "core": (radius <= 0.75) & (az <= 0.75),
        "radial_collar": (radius >= 1.40) & (radius < 1.90) & (az < 1.40),
        "axial_collar": (radius < 1.40) & (az >= 1.20) & (az < 1.90),
        "corner": (radius >= 1.40) & (radius < 1.90) & (az >= 1.20) & (az < 1.90),
    }


def response_diagnostics(axis: np.ndarray, parent_fields: np.ndarray, radial_child_fields: np.ndarray, piola_fields: np.ndarray) -> dict:
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    interior = (np.hypot(xx, yy) < 1.90) & (np.abs(zz) < 1.90)
    radial_delta = radial_child_fields - parent_fields
    piola_delta = piola_fields - radial_child_fields
    columns = [radial_delta[:, interior, :].ravel(), piola_delta[:, interior, :].ravel()]
    raw = np.column_stack(columns)
    norms = np.linalg.norm(raw, axis=0)
    if np.any(norms <= np.finfo(float).tiny):
        raise ValueError("degenerate response column")
    normalized = raw / norms[None, :]
    singular = np.linalg.svd(normalized, compute_uv=False)
    rank = int(np.sum(singular > 1.0e-8))
    condition = float(singular[0] / singular[-1])
    cosine = float(np.dot(normalized[:, 0], normalized[:, 1]))
    projection = np.dot(columns[1], columns[0]) / np.dot(columns[0], columns[0]) * columns[0]
    novelty = float(np.linalg.norm(columns[1] - projection) / np.linalg.norm(columns[1]))
    regions = {}
    for name, mask in _region_masks(axis).items():
        base_rms = float(np.sqrt(np.mean(radial_child_fields[:, mask, :] ** 2)))
        if base_rms <= np.finfo(float).tiny:
            raise ValueError(f"zero baseline region {name}")
        regions[name] = {
            "radial_redistribution_response_over_baseline_rms": float(np.sqrt(np.mean(radial_delta[:, mask, :] ** 2)) / base_rms),
            "piola_response_over_baseline_rms": float(np.sqrt(np.mean(piola_delta[:, mask, :] ** 2)) / base_rms),
        }
    return {
        "column_order": ["frozen_radial_redistribution_parent_to_025", "piola_beta_0075_on_025"],
        "raw_column_norms": [float(x) for x in norms],
        "normalized_rank": rank,
        "normalized_singular_values": [float(x) for x in singular],
        "normalized_condition_number": condition,
        "normalized_column_cosine": cosine,
        "piola_novelty_outside_radial_span": novelty,
        "regions": regions,
    }


def _nested17(axis: np.ndarray, fields: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if len(axis) != 33:
        raise ValueError("nested refinement expects 33 source nodes")
    index = np.arange(0, 33, COARSE_STRIDE)
    return axis[index], fields[:, index][:, :, index][:, :, :, index]


def decision(fine_rows: list[dict], response: dict, structure: dict) -> dict:
    axial = [row["relative"]["mean_axial_span"] for row in fine_rows]
    turns = [row["relative"]["mean_turns"] for row in fine_rows]
    radial = [row["relative"]["mean_radial_span"] for row in fine_rows]
    passed = bool(
        structure["map_monotone"]
        and structure["minimum_coordinate_jacobian"] > 0.0
        and min(axial) >= CRITERIA["minimum_axial_span_gain"]
        and min(turns) >= CRITERIA["minimum_mean_turn_gain"]
        and max(radial) <= CRITERIA["maximum_radial_span_gain"]
        and response["normalized_rank"] >= CRITERIA["minimum_response_rank"]
        and response["normalized_condition_number"] <= CRITERIA["maximum_response_condition"]
        and abs(response["normalized_column_cosine"]) <= CRITERIA["maximum_abs_response_cosine"]
    )
    return {
        "useful_axial_capacity_transfer": passed,
        "criteria": CRITERIA,
        "fine_min_axial_span_gain": float(min(axial)),
        "fine_min_mean_turn_gain": float(min(turns)),
        "fine_max_radial_span_gain": float(max(radial)),
    }


def run(grid_path: str | Path, receipt_path: str | Path | None = None) -> dict:
    axis, child = render.load_child_grid_netcdf(grid_path)
    parent = render.recover_parent_grid(axis, child)
    roundtrip = render.apply_frozen_child_transform(axis, parent)
    roundtrip_error = float(np.max(np.abs(roundtrip - child)))
    if roundtrip_error > 5.0e-13:
        raise ValueError("frozen radial child round-trip drift")
    warped, scale, structure = apply_piola_grid(axis, child, beta=BETA)
    fine = morphology_table(axis, child, warped)
    coarse_axis, coarse_child = _nested17(axis, child)
    coarse_warped, coarse_scale, coarse_structure = apply_piola_grid(coarse_axis, coarse_child, beta=BETA)
    coarse = morphology_table(coarse_axis, coarse_child, coarse_warped)
    response = response_diagnostics(axis, parent, child, warped)
    verdict = decision(fine, response, structure)
    boundary = np.zeros(child.shape[1:4], dtype=bool)
    boundary[[0, -1], :, :] = True
    boundary[:, [0, -1], :] = True
    boundary[:, :, [0, -1]] = True
    boundary_max = float(np.max(np.abs(warped[:, boundary, :])))
    receipt = None
    if receipt_path is not None:
        receipt = json.loads(Path(receipt_path).read_text())
        if receipt.get("receipt_sha256") != render.GRID_SOURCE_RECEIPT_SHA256:
            raise ValueError("source receipt drift")
        if receipt.get("grid", {}).get("grid_sha256") != render.GRID_SOURCE_SHA256:
            raise ValueError("source grid drift")
    return {
        "task_id": TASK_ID,
        "claim_issue": CLAIM_ISSUE,
        "sources": {
            "integration_head": SOURCE_RENDER_HEAD,
            "piola_source_pr": SOURCE_PIOLA_PR,
            "piola_source_head": SOURCE_PIOLA_HEAD,
            "grid_source_pr": render.GRID_SOURCE_PR,
            "grid_source_head": render.GRID_SOURCE_HEAD,
            "grid_source_run": render.GRID_SOURCE_RUN,
            "grid_source_artifact_id": render.GRID_SOURCE_ARTIFACT_ID,
        },
        "contract": {
            "beta": BETA,
            "beta_scan_performed": False,
            "input_child_gain": render.GAIN,
            "times": list(render.REFERENCE_TIMES),
            "fine_grid_resolution": len(axis),
            "coarse_nested_resolution": len(coarse_axis),
            "image_derived_numeric_target_used": False,
            "pressure_or_force_fit_performed": False,
            "held_out_pde_residual_evaluated": False,
        },
        "source_roundtrip_max_abs": roundtrip_error,
        "fine_reference_energy_scale": scale,
        "coarse_reference_energy_scale": coarse_scale,
        "fine_structure": structure,
        "coarse_structure": coarse_structure,
        "fine_boundary_max_abs_velocity": boundary_max,
        "fine_morphology": fine,
        "coarse_morphology": coarse,
        "response_diagnostics": response,
        "decision": verdict,
        "truth": TRUTH,
        "source_receipt_present": receipt is not None,
    }


def save_report(report: dict, path: str | Path) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    path.write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.grid, args.receipt)
    digest = save_report(report, args.out)
    print(json.dumps({"report_sha256": digest, "decision": report["decision"], "response": report["response_diagnostics"]}, indent=2, sort_keys=True))
    return 0 if report["decision"]["useful_axial_capacity_transfer"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
