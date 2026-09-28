"""Rescale-only analysis of the frozen multiscale momentum report.

This module does no field evaluation.  It pulls the completed physical
residual metrics back to the reference coordinates using

    Rhat(y) = s**(3/2) R_s(S y),

and reports the corresponding L2 and maximum norms.  Because the source
report stores aggregate vector metrics rather than every residual component,
the anisotropic acceleration normalization is reported with rigorous scalar
lower and upper bounds instead of an invented exact component norm.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE_PATH = ROOT / "scale_reference_multiscale.json"
OUTPUT_PATH = ROOT / "scale_reference_rescaled_defect.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _save(path: Path, report: dict) -> None:
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def _ratio(value: float, reference: float) -> float:
    return float(value / max(abs(float(reference)), 1.0e-300))


def run(source_path: Path | str = SOURCE_PATH,
        output_path: Path | str = OUTPUT_PATH) -> dict:
    source_path = Path(source_path)
    output_path = Path(output_path)
    source = json.loads(source_path.read_text(encoding="utf-8"))
    if source.get("status") != "completed":
        raise ValueError(f"Multiscale source is not completed: {source.get('status')}")
    inputs = source["inputs"]
    h = float(inputs["similarity_exponent_h"])
    scale_items = []
    for key, entry in source["scales"].items():
        scale = float(key)
        residual = entry["residual"]
        physical_l2 = float(residual["volume_L2"])
        physical_max = float(residual["max_norm"])
        determinant = float(entry["determinant"])
        isotropic_factor = scale ** 1.5
        isotropic_l2_factor = isotropic_factor / math.sqrt(determinant)
        # D_s = diag(s^(3/2), s^(3/2), s^(3/2+h)).  For 0 < s <= 1,
        # the axial factor is the lower scalar singular value.
        horizontal_factor = isotropic_factor
        axial_factor = scale ** (1.5 + h)
        physical_component_l2_lower = axial_factor * physical_l2
        physical_component_l2_upper = horizontal_factor * physical_l2
        physical_component_max_lower = axial_factor * physical_max
        physical_component_max_upper = horizontal_factor * physical_max
        pulled_component_l2_lower = physical_component_l2_lower / math.sqrt(determinant)
        pulled_component_l2_upper = physical_component_l2_upper / math.sqrt(determinant)
        scale_items.append({
            "scale": scale,
            "determinant": determinant,
            "physical": {
                "volume_L2": physical_l2,
                "max_norm": physical_max,
            },
            "isotropic_pulled_back": {
                "definition": "Rhat(y)=s^(3/2) R_s(S y)",
                "L2": isotropic_l2_factor * physical_l2,
                "max_norm": isotropic_factor * physical_max,
                "L2_factor_from_physical": isotropic_l2_factor,
                "max_factor_from_physical": isotropic_factor,
            },
            "component_anisotropic_bounds": {
                "definition": "D_s=diag(s^(3/2),s^(3/2),s^(3/2+h)) applied to R_s",
                "horizontal_factor": horizontal_factor,
                "axial_factor": axial_factor,
                "physical_L2_lower_axial_factor": physical_component_l2_lower,
                "physical_L2_upper_horizontal_factor": physical_component_l2_upper,
                "physical_max_lower_axial_factor": physical_component_max_lower,
                "physical_max_upper_horizontal_factor": physical_component_max_upper,
                "pulled_back_L2_lower_axial_factor": pulled_component_l2_lower,
                "pulled_back_L2_upper_horizontal_factor": pulled_component_l2_upper,
                "bound_reason": "For every vector v, axial_factor*||v|| <= ||D_s v|| <= horizontal_factor*||v||; exact component aggregation is unavailable in the scalar source report.",
            },
        })
    scale_items.sort(key=lambda item: item["scale"], reverse=True)
    reference = scale_items[0]
    reference_iso_l2 = reference["isotropic_pulled_back"]["L2"]
    reference_iso_max = reference["isotropic_pulled_back"]["max_norm"]
    reference_comp_l2 = reference["component_anisotropic_bounds"]["physical_L2_lower_axial_factor"]
    reference_comp_max = reference["component_anisotropic_bounds"]["physical_max_lower_axial_factor"]
    for index, item in enumerate(scale_items):
        iso = item["isotropic_pulled_back"]
        bounds = item["component_anisotropic_bounds"]
        iso["L2_ratio_to_s1"] = _ratio(iso["L2"], reference_iso_l2)
        iso["max_ratio_to_s1"] = _ratio(iso["max_norm"], reference_iso_max)
        bounds["physical_L2_lower_ratio_to_s1"] = _ratio(bounds["physical_L2_lower_axial_factor"], reference_comp_l2)
        bounds["physical_L2_upper_ratio_to_s1"] = _ratio(bounds["physical_L2_upper_horizontal_factor"], reference_comp_l2)
        bounds["physical_max_lower_ratio_to_s1"] = _ratio(bounds["physical_max_lower_axial_factor"], reference_comp_max)
        bounds["physical_max_upper_ratio_to_s1"] = _ratio(bounds["physical_max_upper_horizontal_factor"], reference_comp_max)
        bounds["pulled_L2_lower_ratio_to_s1"] = _ratio(bounds["pulled_back_L2_lower_axial_factor"], reference_iso_l2)
        bounds["pulled_L2_upper_ratio_to_s1"] = _ratio(bounds["pulled_back_L2_upper_horizontal_factor"], reference_iso_l2)
        if index > 0:
            previous = scale_items[index - 1]
            iso_prev = previous["isotropic_pulled_back"]
            iso["L2_ratio_to_previous_larger_scale"] = _ratio(iso["L2"], iso_prev["L2"])
            iso["max_ratio_to_previous_larger_scale"] = _ratio(iso["max_norm"], iso_prev["max_norm"])

    iso_l2_values = [item["isotropic_pulled_back"]["L2"] for item in scale_items]
    iso_max_values = [item["isotropic_pulled_back"]["max_norm"] for item in scale_items]
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "scope": (
            "Metrics-only pullback of the completed multiscale residual report. "
            "The isotropic s^(3/2) normalization is computed exactly from stored "
            "physical L2/max metrics. The anisotropic component normalization is "
            "bounded using its two scalar factors because component-resolved residual "
            "arrays were not saved. No PDE or scale-recursion conclusion is made."
        ),
        "sources": {
            "multiscale_report": {
                "path": source_path.name,
                "sha256": _sha256(source_path),
            },
        },
        "inputs": {
            "similarity_exponent_h": h,
            "scales": [item["scale"] for item in scale_items],
            "source_point_count": int(inputs["point_count"]),
            "source_pressure_rule": inputs.get("pressure_rule"),
            "source_quadrature_rule": inputs.get("quadrature_rule"),
        },
        "normalizations": {
            "isotropic": {
                "residual_definition": "Rhat(y)=s^(3/2) R_s(S y)",
                "L2_definition": "||Rhat||_L2(dy)=s^(3/2)/sqrt(det(S))*||R_s||_L2(dx)",
                "max_definition": "max_y ||Rhat(y)||=s^(3/2)*max_x ||R_s(x)||",
            },
            "component_anisotropic": {
                "matrix": "diag(s^(3/2), s^(3/2), s^(3/2+h))",
                "available_evidence": "scalar physical L2/max metrics only",
                "reported_evidence": "rigorous lower/upper bounds from minimum axial and maximum horizontal factors",
                "exact_component_norm_claim": False,
            },
        },
        "scales": scale_items,
        "summary": {
            "isotropic_pulled_L2_values_descending_scale": iso_l2_values,
            "isotropic_pulled_max_values_descending_scale": iso_max_values,
            "isotropic_L2_contracts_with_halving": bool(all(item["isotropic_pulled_back"]["L2_ratio_to_previous_larger_scale"] <= 1.0 + 1.0e-12 for item in scale_items[1:])),
            "isotropic_max_contracts_with_halving": bool(all(item["isotropic_pulled_back"]["max_ratio_to_previous_larger_scale"] <= 1.0 + 1.0e-12 for item in scale_items[1:])),
            "isotropic_L2_ratio_s025_to_s1": _ratio(iso_l2_values[-1], iso_l2_values[0]),
            "isotropic_max_ratio_s025_to_s1": _ratio(iso_max_values[-1], iso_max_values[0]),
            "interpretation": (
                "The isotropic pulled-back norms remain approximately constant as the "
                "scale halves (small sampled decreases). These decreases are not evidence "
                "of a uniform contraction sufficient for dynamical closure, and no "
                "discretization uncertainty estimate is supplied. Absolute physical "
                "residual L2 and maximum increase toward smaller scales, while that "
                "absolute increase alone does not establish worsening relative dynamical "
                "balance or refute the kinematic recursion mechanism. Absolute residuals "
                "remain far above the 1e-3 target."
            ),
        },
    }
    _save(output_path, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args()
    run(args.source, args.output)
