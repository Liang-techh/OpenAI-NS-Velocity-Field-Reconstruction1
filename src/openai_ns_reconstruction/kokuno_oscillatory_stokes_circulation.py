"""Independent Stokes-circulation audit for the frozen Kokuno A2 oscillation.

This module changes no velocity coefficient, support, phase, scale, or orientation.
It samples only the public ``velocity_osc(x,y,z,t)`` and compares

    integral_{boundary S} u_osc . dl

against

    integral_S curl(u_osc) . n dS

on fixed strictly interior off-axis rectangles.  The surface curl is reconstructed
from public velocity values with a centered FD4 stencil; no production Jacobian,
divergence, or vorticity API participates in the independent path.

The corrected 2026-09-09 Kokuno reconstruction is structural provenance for the
localized complete-curl organization.  The public-z pullback, concrete support,
autonomous mode/background data, quadrature, FD steps, and this audit are
repository realizations, not paper-exact hidden data.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-OSCILLATORY-STOKES-CIRCULATION-077"
SCHEMA = "kokuno-a2-oscillatory-stokes-circulation-v1"
PARENT_AGENT2_PR = 941
PARENT_AGENT2_HEAD = "34e0fd5efab6604f9719ed54a49dd58741ed2308"
PUBLIC_Z_VELOCITY_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"

# Each plane uses oriented local axes (a,b), with normal e_a x e_b.
# The rectangles are strictly interior to the registered radial/axial support and
# intentionally off-axis/non-symmetric so the circulation signal is nontrivial.
SURFACES = (
    {
        "label": "xy_offaxis",
        "center": (0.62, 0.26, -0.55),
        "axis_a": 0,
        "axis_b": 1,
        "normal_axis": 2,
        "halfwidth_a": 0.16,
        "halfwidth_b": 0.13,
    },
    {
        "label": "yz_offaxis",
        "center": (0.58, -0.28, 0.42),
        "axis_a": 1,
        "axis_b": 2,
        "normal_axis": 0,
        "halfwidth_a": 0.14,
        "halfwidth_b": 0.18,
    },
    {
        "label": "zx_offaxis",
        "center": (-0.52, 0.47, -0.36),
        "axis_a": 2,
        "axis_b": 0,
        "normal_axis": 1,
        "halfwidth_a": 0.18,
        "halfwidth_b": 0.15,
    },
)
EVALUATION_TIMES = (0.43, 0.61)
QUADRATURE_ORDER = 24
FD_STEPS = (0.004, 0.002, 0.001)
SHIFT_VECTOR = (0.011, -0.013, 0.017)

# Frozen before exact-head Actions output.  These are scoped integral-consistency
# guards, not CR001 PDE acceptance thresholds.
COARSE_RELATIVE_MISMATCH_MAX = 5.0e-2
MEDIUM_RELATIVE_MISMATCH_MAX = 2.0e-2
FINE_RELATIVE_MISMATCH_MAX = 1.0e-2
SHIFTED_FINE_RELATIVE_MISMATCH_MAX = 1.5e-2
STOKES_SIGNAL_FLOOR = 1.0e-10


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(payload: Mapping[str, Any]) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _as_velocity(value: Any) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.ndim < 1 or out.shape[-1] != 3 or not np.all(np.isfinite(out)):
        raise RuntimeError("public oscillatory velocity returned an invalid Cartesian array")
    return out


def _surface_spec(raw: Mapping[str, Any], shift: np.ndarray | None = None) -> dict[str, Any]:
    center = np.asarray(raw["center"], dtype=float)
    if shift is not None:
        center = center + np.asarray(shift, dtype=float)
    spec = {
        "label": str(raw["label"]),
        "center": center,
        "axis_a": int(raw["axis_a"]),
        "axis_b": int(raw["axis_b"]),
        "normal_axis": int(raw["normal_axis"]),
        "halfwidth_a": float(raw["halfwidth_a"]),
        "halfwidth_b": float(raw["halfwidth_b"]),
    }
    axes = (spec["axis_a"], spec["axis_b"], spec["normal_axis"])
    if sorted(axes) != [0, 1, 2]:
        raise RuntimeError("Stokes surface axes must be one oriented Cartesian permutation")
    ea = np.eye(3)[spec["axis_a"]]
    eb = np.eye(3)[spec["axis_b"]]
    en = np.eye(3)[spec["normal_axis"]]
    if not np.allclose(np.cross(ea, eb), en, rtol=0.0, atol=0.0):
        raise RuntimeError("Stokes surface orientation is not right-handed")
    return spec


def _surface_corners(spec: Mapping[str, Any]) -> np.ndarray:
    center = np.asarray(spec["center"], dtype=float)
    ea = np.eye(3)[int(spec["axis_a"])]
    eb = np.eye(3)[int(spec["axis_b"])]
    ha = float(spec["halfwidth_a"])
    hb = float(spec["halfwidth_b"])
    return np.asarray(
        [center + sa * ha * ea + sb * hb * eb for sa in (-1.0, 1.0) for sb in (-1.0, 1.0)],
        dtype=float,
    )


def _assert_surface_inside_support(spec: Mapping[str, Any], margin: float) -> None:
    field = default_field()
    r0 = float(field.radial_inner)
    r1 = float(field.radial_outer)
    z0 = float(field.axial_lower)
    z1 = float(field.axial_upper)
    corners = _surface_corners(spec)
    # FD points can extend by two steps along either tangent direction.
    pad = 2.0 * float(margin)
    samples = [corners]
    for axis in (int(spec["axis_a"]), int(spec["axis_b"])):
        direction = np.eye(3)[axis]
        samples.append(corners + pad * direction)
        samples.append(corners - pad * direction)
    cloud = np.concatenate(samples, axis=0)
    radii = np.hypot(cloud[:, 0], cloud[:, 1])
    if not (np.min(radii) > r0 and np.max(radii) < r1):
        raise RuntimeError("registered Stokes surface or FD collar left radial support")
    if not (np.min(cloud[:, 2]) > z0 and np.max(cloud[:, 2]) < z1):
        raise RuntimeError("registered Stokes surface or FD collar left axial support")


def _gauss_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if int(order) < 4:
        raise ValueError("Gauss-Legendre order must be at least four")
    nodes, weights = np.polynomial.legendre.leggauss(int(order))
    return np.asarray(nodes, dtype=float), np.asarray(weights, dtype=float)


def _evaluate_points(points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.shape[-1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be one finite Cartesian [...,3] array")
    tt = np.full(points.shape[:-1], float(time), dtype=float)
    return _as_velocity(velocity_osc(points[..., 0], points[..., 1], points[..., 2], tt))


def _boundary_circulation(spec: Mapping[str, Any], time: float, order: int) -> float:
    nodes, weights = _gauss_rule(order)
    center = np.asarray(spec["center"], dtype=float)
    ia = int(spec["axis_a"])
    ib = int(spec["axis_b"])
    ea = np.eye(3)[ia]
    eb = np.eye(3)[ib]
    ha = float(spec["halfwidth_a"])
    hb = float(spec["halfwidth_b"])

    # Positive boundary orientation induced by n=e_a x e_b:
    # bottom (+a), right (+b), top (-a), left (-b).
    edges = (
        (center - hb * eb, ha * nodes, ea, ia, +ha),
        (center + ha * ea, hb * nodes, eb, ib, +hb),
        (center + hb * eb, -ha * nodes, ea, ia, -ha),
        (center - ha * ea, -hb * nodes, eb, ib, -hb),
    )
    circulation = 0.0
    for base, coordinate, direction, component, jacobian in edges:
        points = base[None, :] + coordinate[:, None] * direction[None, :]
        velocity = _evaluate_points(points, time)
        circulation += float(np.sum(weights * velocity[:, component] * jacobian))
    return circulation


def _fd4_component_derivative(
    points: np.ndarray,
    time: float,
    *,
    component: int,
    axis: int,
    step: float,
) -> np.ndarray:
    h = float(step)
    if not np.isfinite(h) or h <= 0.0:
        raise ValueError("FD step must be positive and finite")
    direction = np.eye(3)[int(axis)]
    fm2 = _evaluate_points(points - 2.0 * h * direction, time)[..., int(component)]
    fm1 = _evaluate_points(points - h * direction, time)[..., int(component)]
    fp1 = _evaluate_points(points + h * direction, time)[..., int(component)]
    fp2 = _evaluate_points(points + 2.0 * h * direction, time)[..., int(component)]
    return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)


def _surface_curl_flux(
    spec: Mapping[str, Any],
    time: float,
    order: int,
    step: float,
) -> float:
    nodes, weights = _gauss_rule(order)
    center = np.asarray(spec["center"], dtype=float)
    ia = int(spec["axis_a"])
    ib = int(spec["axis_b"])
    ea = np.eye(3)[ia]
    eb = np.eye(3)[ib]
    ha = float(spec["halfwidth_a"])
    hb = float(spec["halfwidth_b"])

    aa, bb = np.meshgrid(ha * nodes, hb * nodes, indexing="ij")
    points = center + aa[..., None] * ea + bb[..., None] * eb
    # For n=e_a x e_b, curl(u).n = partial_a u_b - partial_b u_a.
    d_a_u_b = _fd4_component_derivative(points, time, component=ib, axis=ia, step=step)
    d_b_u_a = _fd4_component_derivative(points, time, component=ia, axis=ib, step=step)
    curl_normal = d_a_u_b - d_b_u_a
    w2 = weights[:, None] * weights[None, :]
    return float(ha * hb * np.sum(w2 * curl_normal))


def _stokes_case(
    spec: Mapping[str, Any],
    time: float,
    *,
    order: int,
    step: float,
) -> dict[str, float]:
    circulation = _boundary_circulation(spec, time, order)
    curl_flux = _surface_curl_flux(spec, time, order, step)
    signal = max(abs(circulation), abs(curl_flux))
    denominator = max(abs(circulation) + abs(curl_flux), STOKES_SIGNAL_FLOOR)
    relative = abs(circulation - curl_flux) / denominator
    return {
        "boundary_circulation": float(circulation),
        "surface_curl_flux": float(curl_flux),
        "absolute_mismatch": float(abs(circulation - curl_flux)),
        "relative_mismatch": float(relative),
        "signal_scale": float(signal),
    }


@dataclass(frozen=True)
class OscillatoryStokesCirculationResult:
    surfaces: tuple[dict[str, Any], ...]
    times: tuple[float, ...]
    levels: tuple[dict[str, Any], ...]
    shifted_fine: dict[str, Any]

    def summary_payload(self) -> dict[str, Any]:
        def serial_surface(spec: Mapping[str, Any]) -> dict[str, Any]:
            return {
                "label": str(spec["label"]),
                "center": np.asarray(spec["center"], dtype=float).tolist(),
                "axis_a": int(spec["axis_a"]),
                "axis_b": int(spec["axis_b"]),
                "normal_axis": int(spec["normal_axis"]),
                "halfwidth_a": float(spec["halfwidth_a"]),
                "halfwidth_b": float(spec["halfwidth_b"]),
            }
        return {
            "surfaces": [serial_surface(spec) for spec in self.surfaces],
            "times": list(self.times),
            "levels": [dict(level) for level in self.levels],
            "shifted_fine": dict(self.shifted_fine),
        }


class KokunoOscillatoryStokesCirculation:
    """Materialize the frozen implementation-distinct Stokes receipt."""

    def semantic_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "task": TASK,
            "parent_agent2_pr": PARENT_AGENT2_PR,
            "parent_agent2_head": PARENT_AGENT2_HEAD,
            "public_z_velocity_blob": PUBLIC_Z_VELOCITY_BLOB,
            "source_provenance": {
                "corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
                "corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
                "scope": "localized oscillatory fields and complete-curl organization",
            },
            "surfaces": [
                {
                    **{k: v for k, v in raw.items() if k != "center"},
                    "center": list(raw["center"]),
                }
                for raw in SURFACES
            ],
            "times": list(EVALUATION_TIMES),
            "quadrature_order": QUADRATURE_ORDER,
            "fd_steps": list(FD_STEPS),
            "shift_vector": list(SHIFT_VECTOR),
            "numerical_contract": {
                "boundary_path": "public velocity_osc values + oriented Gauss-Legendre line integral",
                "surface_path": "public velocity_osc values + independent centered-FD4 curl + Gauss-Legendre surface integral",
                "stokes_identity": "integral_boundary u.dl = integral_surface curl(u).n dS",
                "uses_production_jacobian": False,
                "uses_production_divergence": False,
                "uses_production_vorticity": False,
                "strictly_interior_off_axis_surfaces": True,
            },
            "frozen_engineering_gates": {
                "coarse_relative_mismatch_max": COARSE_RELATIVE_MISMATCH_MAX,
                "medium_relative_mismatch_max": MEDIUM_RELATIVE_MISMATCH_MAX,
                "fine_relative_mismatch_max": FINE_RELATIVE_MISMATCH_MAX,
                "shifted_fine_relative_mismatch_max": SHIFTED_FINE_RELATIVE_MISMATCH_MAX,
                "stokes_signal_floor": STOKES_SIGNAL_FLOOR,
            },
            "scientific_boundary": {
                "oscillatory_velocity_changed": False,
                "oscillatory_support_changed": False,
                "new_oscillatory_parameters_added": False,
                "stokes_circulation_consistency_assessed": True,
                "whole_domain_vorticity_assessed": False,
                "whole_domain_divergence_assessed": False,
                "mean_projection_performed": False,
                "radial_inverse_performed": False,
                "correction_velocity_constructed": False,
                "pressure_included": False,
                "forcing_included": False,
                "complete_ns_residual": False,
                "same_protocol_st006_comparison_valid": False,
                "residual_reduction_claimed": False,
                "paper_exact": False,
                "pde_validated": False,
            },
        }

    @property
    def stokes_circulation_sha256(self) -> str:
        return _sha256(self.semantic_payload())

    def materialize(self) -> OscillatoryStokesCirculationResult:
        surfaces = tuple(_surface_spec(raw) for raw in SURFACES)
        max_step = max(FD_STEPS)
        for spec in surfaces:
            _assert_surface_inside_support(spec, max_step)

        levels: list[dict[str, Any]] = []
        for step in FD_STEPS:
            cases: list[dict[str, Any]] = []
            for surface_index, spec in enumerate(surfaces):
                for time in EVALUATION_TIMES:
                    case = _stokes_case(
                        spec,
                        time,
                        order=QUADRATURE_ORDER,
                        step=step,
                    )
                    cases.append({
                        "surface_index": surface_index,
                        "surface_label": spec["label"],
                        "time": float(time),
                        **case,
                    })
            levels.append({
                "fd_step": float(step),
                "cases": cases,
                "max_relative_mismatch": float(max(case["relative_mismatch"] for case in cases)),
                "min_signal_scale": float(min(case["signal_scale"] for case in cases)),
            })

        shift = np.asarray(SHIFT_VECTOR, dtype=float)
        shifted_surfaces = tuple(_surface_spec(raw, shift=shift) for raw in SURFACES)
        for spec in shifted_surfaces:
            _assert_surface_inside_support(spec, max_step)
        shifted_cases: list[dict[str, Any]] = []
        for surface_index, spec in enumerate(shifted_surfaces):
            for time in EVALUATION_TIMES:
                case = _stokes_case(
                    spec,
                    time,
                    order=QUADRATURE_ORDER,
                    step=FD_STEPS[-1],
                )
                shifted_cases.append({
                    "surface_index": surface_index,
                    "surface_label": spec["label"],
                    "time": float(time),
                    **case,
                })
        shifted_fine = {
            "fd_step": float(FD_STEPS[-1]),
            "shift_vector": list(SHIFT_VECTOR),
            "cases": shifted_cases,
            "max_relative_mismatch": float(max(case["relative_mismatch"] for case in shifted_cases)),
            "min_signal_scale": float(min(case["signal_scale"] for case in shifted_cases)),
        }
        return OscillatoryStokesCirculationResult(
            surfaces=surfaces,
            times=tuple(float(t) for t in EVALUATION_TIMES),
            levels=tuple(levels),
            shifted_fine=shifted_fine,
        )

    def manifest(self) -> dict[str, Any]:
        payload = self.semantic_payload()
        return {"payload": payload, "stokes_circulation_sha256": _sha256(payload)}

    def save_manifest(self, path: str | Path) -> dict[str, Any]:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        manifest = self.manifest()
        target.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    @classmethod
    def load_manifest(cls, path: str | Path) -> "KokunoOscillatoryStokesCirculation":
        raw = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or set(raw) != {"payload", "stokes_circulation_sha256"}:
            raise ValueError("invalid oscillatory Stokes manifest envelope")
        payload = raw["payload"]
        if not isinstance(payload, dict) or _sha256(payload) != raw["stokes_circulation_sha256"]:
            raise ValueError("oscillatory Stokes manifest checksum mismatch")
        canonical = cls()
        if payload != canonical.semantic_payload():
            raise ValueError("oscillatory Stokes manifest semantics mismatch")
        return canonical


def public_contract() -> dict[str, Any]:
    forbidden = {
        "amplitude", "phase", "scale", "orientation", "support",
        "pressure", "forcing", "viscosity", "nu", "residual", "target",
        "gain", "damping", "threshold", "mean", "stress", "inverse",
    }
    params = set(inspect.signature(KokunoOscillatoryStokesCirculation.materialize).parameters)
    return {
        "task": TASK,
        "stokes_circulation_sha256": KokunoOscillatoryStokesCirculation().stokes_circulation_sha256,
        "forbidden_materialize_inputs_present": sorted(forbidden.intersection(params)),
        "uses_only_public_velocity_values": True,
        "production_jacobian_used": False,
        "production_divergence_used": False,
        "production_vorticity_used": False,
        **KokunoOscillatoryStokesCirculation().semantic_payload()["scientific_boundary"],
    }
