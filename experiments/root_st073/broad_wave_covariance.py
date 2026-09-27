"""Exact-curl mode-1 wave covariance versus the current mean stress target.

The wave is reconstructed directly from the selected degree-2 potential
coefficients in ``broad_shear_growth.json`` and ``supported_fourier_basis``.
That basis evaluates all compact radial/axial cutoff derivatives and the
cylindrical curl connection terms.  The physical realization used here is
``Re(V c)`` with an explicit factor of one; no eigenvector renormalization is
hidden in the covariance.

The saved broad dynamic-control cone rows provide the raw mean target.  The
constrained initial-mean target is reconstructed from its saved ``H``, ``A``,
``b`` and selected coefficients, which is exact for the saved linear
zero-instantaneous-velocity directions.  This remains an instantaneous
finite-node covariance diagnostic and makes no PDE, finite-time, or recursion
acceptance claim.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from broad_shear_dynamic_control import (  # noqa: E402
    load_saved_field as load_dynamic_field,
    wave_locations,
)
from affine_momentum import momentum  # noqa: E402
from grouped_outer_cache import outer_cache  # noqa: E402
from supported_fourier_basis import basis_data  # noqa: E402


GROWTH_PATH = ROOT / "broad_shear_growth.json"
DYNAMIC_PATH = ROOT / "broad_shear_dynamic_control.json"
CONSTRAINED_PATH = ROOT / "broad_meridional_constrained.json"
OUTPUT_PATH = ROOT / "broad_wave_covariance.json"
ANGLES = (8, 16)
MARGIN = np.array([0.02, 0.0, 0.0], dtype=float)


def _save(report):
    OUTPUT_PATH.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _decode_complex(packed):
    values = np.asarray(packed, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("packed complex coefficients must have shape (N, 2)")
    return values[:, 0] + 1j * values[:, 1]


class ExactCurlWave:
    """Callable real realization of the saved mode-1 exact-curl field."""

    def __init__(self, center, widths, degree, carrier, coefficients):
        self.center = tuple(float(value) for value in center)
        self.widths = tuple(float(value) for value in widths)
        self.degree = int(degree)
        self.carrier = tuple(float(value) for value in carrier)
        self.coefficients = np.asarray(coefficients, dtype=complex)
        expected = 3 * (self.degree + 1) ** 2
        if self.coefficients.shape != (expected,):
            raise ValueError(f"expected {expected} potential coefficients")

    def complex_velocity(self, points):
        velocity, _, _ = basis_data(
            points,
            center=self.center,
            widths=self.widths,
            mode=1,
            degree=self.degree,
            carrier=self.carrier,
        )
        return np.einsum("niq,q->ni", velocity, self.coefficients)

    def fields(self, points, tau):
        points = np.asarray(points, dtype=float)
        return np.real(self.complex_velocity(points)), np.zeros(len(points))


def load_saved_wave(path=OUTPUT_PATH):
    """Load the serialized exact-curl wave and its covariance report."""

    report = json.loads(Path(path).read_text(encoding="utf-8"))
    metadata = report["wave_field"]
    coefficients = _decode_complex(metadata["potential_coefficients"])
    return ExactCurlWave(
        metadata["center"],
        metadata["widths"],
        metadata["degree"],
        metadata["carrier"],
        coefficients,
    ), report


def _cylindrical_velocity(complex_velocity, theta, realization_factor=1.0):
    """Convert the explicit real wave realization to cylindrical components."""

    theta = np.asarray(theta, dtype=float)
    velocity = realization_factor * np.real(complex_velocity)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    return np.column_stack((
        cos_theta * velocity[:, 0] + sin_theta * velocity[:, 1],
        -sin_theta * velocity[:, 0] + cos_theta * velocity[:, 1],
        velocity[:, 2],
    ))


def _covariance_at(wave, radius, z, angle_count):
    angles = 2.0 * np.pi * np.arange(int(angle_count)) / float(angle_count)
    points = np.column_stack((
        float(radius) * np.cos(angles),
        float(radius) * np.sin(angles),
        np.full(len(angles), float(z)),
    ))
    complex_velocity = wave.complex_velocity(points)
    cylindrical = _cylindrical_velocity(complex_velocity, angles, 1.0)
    covariance = cylindrical.T @ cylindrical / float(len(angles))
    return covariance, cylindrical


def _physical_wave_nodes(dynamic, wave_report, k):
    tau = 0.5 * 2.0 ** (-float(k))
    locations = wave_locations(wave_report)
    nodes = []
    for eta, y in locations:
        point = dynamic.inner.from_similarity(
            [dynamic.inner.p.X_max * (1.0 + 15.0 * y) ** 2],
            [eta],
            tau,
        )[0]
        nodes.append({
            "eta": float(eta),
            "radial_fraction": float(y),
            "point": np.asarray(point, dtype=float).tolist(),
        })
    return nodes


def _saved_dynamic_targets(dynamic_report, nodes):
    """Read the actual grouped-cache targets saved by dynamic-control."""

    diagnostics = [
        row for row in dynamic_report["cone_problem"]["diagnostics"]
        if row.get("kind") == "wave_center"
    ]
    result = []
    for node in nodes:
        matches = [
            row for row in diagnostics
            if abs(float(row["label"]["eta"]) - node["eta"]) < 1.0e-12
            and abs(float(row["label"]["y"]) - node["radial_fraction"]) < 1.0e-12
        ]
        if len(matches) != 1:
            raise ValueError("dynamic-control wave target node is missing or duplicated")
        row = matches[0]
        source_target = np.asarray(row["source_T"], dtype=float)
        control_contribution = (
            np.asarray(row["unit_target"], dtype=float)
            @ np.asarray(dynamic_report["control"], dtype=float)
        )
        target = source_target + control_contribution
        normal = np.asarray(row["normal"], dtype=float)
        tangent = np.array([-normal[1], normal[0]])
        result.append({
            "T": target.tolist(),
            "baseline_source_T": source_target.tolist(),
            "control_contribution": control_contribution.tolist(),
            "target_dot_N": float(target @ normal),
            "target_dot_K": float(target @ tangent),
            "source": "broad_shear_dynamic_control.json cone_problem diagnostics; grouped outer cache/full jets",
        })
    return result


def _saved_constrained_targets(constrained_report, nodes):
    """Recover constrained mean targets from saved H/A/b rows."""

    diagnostics = constrained_report["cone_problem"]["diagnostics"]
    A = np.asarray(constrained_report["cone_problem"]["A"], dtype=float)
    b = np.asarray(constrained_report["cone_problem"]["b"], dtype=float)
    coefficients = np.asarray(constrained_report["coefficients"], dtype=float)
    result = []
    for node in nodes:
        matches = [
            (index, row) for index, row in enumerate(diagnostics)
            if row.get("kind") == "wave_center"
            and abs(float(row["label"]["eta"]) - node["eta"]) < 1.0e-12
            and abs(float(row["label"]["y"]) - node["radial_fraction"]) < 1.0e-12
        ]
        if len(matches) != 1:
            raise ValueError("constrained wave target node is missing or duplicated")
        index, row = matches[0]
        H = np.asarray(row["H"], dtype=float)
        transformed = (
            A[3 * index:3 * index + 3] @ coefficients
            + b[3 * index:3 * index + 3]
            + MARGIN
        )
        target, _, _, _ = np.linalg.lstsq(H, transformed, rcond=None)
        reconstruction_error = float(np.max(np.abs(H @ target - transformed)))
        result.append({
            "T": target,
            "H": H,
            "transformed_target": transformed,
            "reconstruction_error_max": reconstruction_error,
            "source": "broad_meridional_constrained.json saved H/A/b and selected coefficients",
        })
    return result


def _fit_nonnegative_scalar(W_flux, targets):
    W_flux = np.asarray(W_flux, dtype=float)
    targets = np.asarray(targets, dtype=float)
    numerator = float(np.sum(W_flux * targets))
    denominator = float(np.sum(W_flux * W_flux))
    unconstrained = numerator / denominator if denominator > 0.0 else 0.0
    alpha = max(0.0, unconstrained)
    predicted = alpha * W_flux
    target_norm = float(np.linalg.norm(targets))
    predicted_norm = float(np.linalg.norm(predicted))
    residual = predicted - targets
    cosine = float(numerator / max(
        np.linalg.norm(W_flux) * target_norm, np.finfo(float).tiny
    ))
    node_cosines = np.sum(W_flux * targets, axis=1) / np.maximum(
        np.linalg.norm(W_flux, axis=1) * np.linalg.norm(targets, axis=1),
        np.finfo(float).tiny,
    )
    return {
        "unconstrained_alpha": unconstrained,
        "alpha_nonnegative": alpha,
        "amplitude_nonnegative": float(np.sqrt(alpha)),
        "predicted_flux": predicted.tolist(),
        "residual_flux": residual.tolist(),
        "relative_l2_mismatch": float(np.linalg.norm(residual)
                                       / max(target_norm, np.finfo(float).tiny)),
        "global_direction_cosine": cosine,
        "node_direction_cosines": node_cosines.tolist(),
        "positive_alpha_directional_match": bool(alpha > 0.0 and cosine > 0.0),
        "target_norm": target_norm,
        "predicted_norm": predicted_norm,
        "fit_definition": "minimize ||alpha W_flux - T||_2 over alpha >= 0; alpha is amplitude squared relative to saved coefficient vector",
    }


def _direct_dynamic_target_check(dynamic, dynamic_report, node, saved_target):
    """Verify one saved target with a fresh grouped-cache full-jet evaluation."""

    cache = outer_cache(
        dynamic,
        [],
        k=float(dynamic_report["k"]),
        order=8,
        locations=[(node["eta"], node["radial_fraction"])],
        radial_breaks=dynamic_report["radial_breaks"],
    )
    residual = momentum(cache["baseline"])
    sl, radii, weights, radius = cache["panels"][0]
    direct_target = np.array([
        -np.dot(weights * radii**2, residual[sl, 1]) / radius**2,
        -np.dot(weights * radii, residual[sl, 2]) / radius,
    ])
    saved_target = np.asarray(saved_target, dtype=float)
    return {
        "order": 8,
        "point_count": int(len(cache["points"])),
        "direct_T": direct_target.tolist(),
        "saved_T": saved_target.tolist(),
        "difference": (direct_target - saved_target).tolist(),
        "maximum_absolute_difference": float(np.max(np.abs(direct_target - saved_target))),
        "relative_l2_difference": float(
            np.linalg.norm(direct_target - saved_target)
            / max(np.linalg.norm(saved_target), np.finfo(float).tiny)
        ),
        "source": "fresh grouped_outer_cache full finite-difference jets at first wave node",
    }


def run():
    growth = json.loads(GROWTH_PATH.read_text(encoding="utf-8"))
    dynamic_report = json.loads(DYNAMIC_PATH.read_text(encoding="utf-8"))
    constrained_report = json.loads(CONSTRAINED_PATH.read_text(encoding="utf-8"))
    selected = growth["selected_candidate"]
    if int(selected["mode"]) != 1:
        raise ValueError("broad_shear_growth selected candidate is not mode 1")
    dynamic, loaded_dynamic_report = load_dynamic_field(DYNAMIC_PATH)
    k = float(loaded_dynamic_report["k"])
    wave_report = json.loads((ROOT / "broad_shear_wave_cone.json").read_text(
        encoding="utf-8"
    ))
    nodes = _physical_wave_nodes(dynamic, wave_report, k)
    center = tuple(float(value) for value in growth["center"])
    widths = tuple(float(value) for value in growth["widths"])
    degree = int(growth["degree"])
    carrier = tuple(float(value) for value in growth["carriers"]["1"])
    coefficients = _decode_complex(selected["potential_coefficients"])
    wave = ExactCurlWave(center, widths, degree, carrier, coefficients)
    report = {
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "status": "initialized",
        "k": k,
        "tau": float(0.5 * 2.0 ** (-k)),
        "source_growth": GROWTH_PATH.name,
        "source_dynamic_mean": DYNAMIC_PATH.name,
        "source_constrained_mean": CONSTRAINED_PATH.name,
        "wave_field": {
            "mode": 1,
            "center": list(center),
            "widths": list(widths),
            "degree": degree,
            "carrier": list(carrier),
            "potential_coefficients": selected["potential_coefficients"],
            "potential_coefficient_dimension": int(len(coefficients)),
            "coefficient_l2_norm": float(np.linalg.norm(coefficients)),
            "coefficient_packing": "[real, imaginary] per coefficient; radial/theta/axial blocks and a-major,b-minor polynomial order from supported_fourier_basis",
            "curl_definition": "supported_fourier_basis.basis_data velocity columns = -curl(A), including cylindrical 1/r terms and compact bump cutoff derivatives",
            "realization": "w = 1.0 * Re(sum_q basis_data.velocity_q * potential_coefficient_q); no conjugate duplicate and no hidden eigenvector normalization",
            "growth_amplitude_metadata": float(selected["amplitude"]),
            "loader": "load_saved_wave",
        },
        "wave_nodes": nodes,
        "angular_quadrature": {
            "counts": list(ANGLES),
            "weights": "uniform 2*pi/N angular average",
            "cylindrical_components": ["r", "theta", "z"],
            "covariance_definition": "W_ab = mean_theta w_a(theta) w_b(theta), with w = Re(V c)",
        },
        "target_sign_convention": {
            "definition": "T = [-integral(r^2 R_theta) / R^2, -integral(r R_z) / R]",
            "paper_relation": "E contains div(W) - div(Sigma), so the required wave radial flux pair is compared directly to T",
            "pair_order": ["W_rtheta", "W_rz"],
        },
        "covariance_nodes": [],
        "target_comparison": {},
        "scope": "Instantaneous exact-curl mode-1 covariance and finite wave-node radial-flux direction diagnostic. No PDE, finite-time, continuum, total-energy, or scale-recursion acceptance claim.",
    }
    _save(report)

    raw_targets = _saved_dynamic_targets(dynamic_report, nodes)
    constrained_targets = _saved_constrained_targets(constrained_report, nodes)
    report["raw_target_direct_check"] = _direct_dynamic_target_check(
        dynamic,
        dynamic_report,
        nodes[0],
        raw_targets[0]["T"],
    )
    covariance_rows = []
    for node in nodes:
        point = np.asarray(node["point"], dtype=float)
        covariance = {}
        cylindrical_peak = {}
        for count in ANGLES:
            W, velocities = _covariance_at(wave, point[0], point[2], count)
            covariance[str(count)] = W.tolist()
            cylindrical_peak[str(count)] = float(np.max(
                np.linalg.norm(velocities, axis=1)
            ))
        W8 = np.asarray(covariance[str(ANGLES[0])], dtype=float)
        W16 = np.asarray(covariance[str(ANGLES[1])], dtype=float)
        covariance_rows.append({
            **node,
            "covariance": covariance,
            "flux_pair_8": [float(W8[0, 1]), float(W8[0, 2])],
            "flux_pair_16": [float(W16[0, 1]), float(W16[0, 2])],
            "angular_max_abs_difference_8_vs_16": float(np.max(np.abs(W8 - W16))),
            "cylindrical_velocity_peak": cylindrical_peak,
        })
    report["covariance_nodes"] = covariance_rows
    W_flux = np.asarray([row["flux_pair_16"] for row in covariance_rows])
    raw_T = np.asarray([row["T"] for row in raw_targets])
    constrained_T = np.asarray([row["T"] for row in constrained_targets])
    report["target_comparison"] = {
        "raw_dynamic_control": {
            "targets": raw_T.tolist(),
            "metadata": raw_targets,
            "fit": _fit_nonnegative_scalar(W_flux, raw_T),
        },
        "constrained_initial_mean": {
            "targets": constrained_T.tolist(),
            "metadata": [
                {
                    "source": item["source"],
                    "transformed_target": item["transformed_target"].tolist(),
                    "reconstruction_error_max": item["reconstruction_error_max"],
                }
                for item in constrained_targets
            ],
            "fit": _fit_nonnegative_scalar(W_flux, constrained_T),
        },
    }
    report["angular_convergence_max_abs"] = float(max(
        row["angular_max_abs_difference_8_vs_16"] for row in covariance_rows
    ))
    report["status"] = "completed"
    _save(report)
    print(json.dumps({
        "stage": "completed",
        "status": report["status"],
        "angular_convergence_max_abs": report["angular_convergence_max_abs"],
        "raw_alpha": report["target_comparison"]["raw_dynamic_control"]["fit"]["alpha_nonnegative"],
        "raw_cosine": report["target_comparison"]["raw_dynamic_control"]["fit"]["global_direction_cosine"],
        "constrained_alpha": report["target_comparison"]["constrained_initial_mean"]["fit"]["alpha_nonnegative"],
        "constrained_cosine": report["target_comparison"]["constrained_initial_mean"]["fit"]["global_direction_cosine"],
    }), flush=True)
    return report


if __name__ == "__main__":
    run()
