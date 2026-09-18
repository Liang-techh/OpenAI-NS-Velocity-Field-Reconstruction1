"""Independent black-box audit of Agent 2's concrete slow partition.

Consumes only the public ``KokunoSourceCompatiblePartitionRealization.evaluate``
value map.  It does not call the construction's raw bump, normalization, or
analytic derivative helpers.  Derivatives are rechecked with a centered
seven-point sixth-order finite difference on fresh held-out inputs, aligning
finite active beta sets by public label identity.

This is a local structural preflight, not a Navier--Stokes PDE validation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_source_compatible_partition import KokunoSourceCompatiblePartitionRealization

SCHEMA = "kokuno-agent4-source-compatible-partition-independent-audit-v1"
TASK_ID = "KOKUNO-A4-SOURCE-COMPATIBLE-PARTITION-INDEPENDENT-AUDIT-021"
BASE_PR = 392
BASE_HEAD = "e2f73ccdc4eccaebf03fd20e3b3d9c3d9df53977"
SEED = 9173121
N_HELD_OUT = 24
FD6_STEPS = (0.04, 0.02, 0.01)
ORIGIN_CASES = ((0.0, 0.0, 0.0), (1.3e-5, -0.7e-5, 2.1e-5))

PARTITION_MAX_ABS_GUARD = 5.0e-12
DERIVATIVE_CLOSURE_MAX_ABS_GUARD = 5.0e-10
FD6_FINE_RELATIVE_RMS_GUARD = 1.0e-6
FD6_MIN_REFINEMENT_RATIO_GUARD = 8.0
DROPPED_LABEL_PARTITION_ERROR_FLOOR = 1.0e-4
DROPPED_DERIVATIVE_RELATIVE_ERROR_FLOOR = 1.0e-3
FORMAL_MOMENTUM_GATE = 1.0e-3
FORMAL_DIVERGENCE_GATE = 1.0e-5


def _fresh_inputs() -> dict[str, np.ndarray]:
    rng = np.random.default_rng(SEED)
    y = rng.uniform(5.15, 7.85, N_HELD_OUT)
    q = np.power(2.0, -y)
    D_r_y = rng.uniform(-0.12, 0.12, N_HELD_OUT)
    D_z_y = rng.uniform(-0.12, 0.12, N_HELD_OUT)
    return {
        "q": q,
        "D_r_q": -q * math.log(2.0) * D_r_y,
        "D_z_q": -q * math.log(2.0) * D_z_y,
        "slow": rng.uniform(-1.5e-4, 1.5e-4, (N_HELD_OUT, 3)),
        "D_r_slow": rng.uniform(-2.0e-6, 2.0e-6, (N_HELD_OUT, 3)),
        "D_z_slow": rng.uniform(-2.0e-6, 2.0e-6, (N_HELD_OUT, 3)),
    }


def _evaluate(realization, inputs, shift=0.0, direction="r"):
    if direction == "r":
        q = inputs["q"] + shift * inputs["D_r_q"]
        slow = inputs["slow"] + shift * inputs["D_r_slow"]
    elif direction == "z":
        q = inputs["q"] + shift * inputs["D_z_q"]
        slow = inputs["slow"] + shift * inputs["D_z_slow"]
    else:
        raise ValueError(direction)
    return realization.evaluate(
        q, inputs["D_r_q"], inputs["D_z_q"], slow,
        inputs["D_r_slow"], inputs["D_z_slow"]
    )


def _union_labels(records):
    labels = set()
    for record in records:
        labels.update(record["beta_labels"])
    return tuple(sorted(labels))


def _align(values, labels, union):
    values = np.asarray(values, dtype=float)
    if values.shape[-1] != len(labels):
        raise ValueError("beta-array width does not match labels")
    index = {label: i for i, label in enumerate(labels)}
    if len(index) != len(labels):
        raise ValueError("beta labels must be unique")
    out = np.zeros(values.shape[:-1] + (len(union),), dtype=float)
    for j, label in enumerate(union):
        i = index.get(label)
        if i is not None:
            out[..., j] = values[..., i]
    return out


def _rms(value):
    value = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(value * value)))


def _fd6(realization, inputs, base, direction, step):
    ks = (-3, -2, -1, 1, 2, 3)
    shifted = {k: _evaluate(realization, inputs, k * step, direction) for k in ks}
    union = _union_labels([base, *shifted.values()])
    f = {k: _align(v["eta"], v["beta_labels"], union) for k, v in shifted.items()}
    derivative = (-f[-3] + 9*f[-2] - 45*f[-1] + 45*f[1] - 9*f[2] + f[3]) / (60*step)
    key = "D_r_eta" if direction == "r" else "D_z_eta"
    analytic = _align(base[key], base["beta_labels"], union)
    error = derivative - analytic
    return {
        "step": float(step),
        "aligned_beta_count": len(union),
        "analytic_rms": _rms(analytic),
        "fd6_rms": _rms(derivative),
        "error_rms": _rms(error),
        "error_max_abs": float(np.max(np.abs(error))),
        "relative_rms": _rms(error) / max(_rms(analytic), 1.0e-300),
    }


def _active_ell_span(record):
    eta = np.asarray(record["eta"], dtype=float)
    labels = record["beta_labels"]
    worst = 0
    for row in eta:
        active = np.flatnonzero(np.abs(row) > 1.0e-14)
        if active.size:
            ells = [int(labels[i][0]) for i in active]
            worst = max(worst, max(ells) - min(ells))
    return worst


def _origin_report(origin):
    inputs = _fresh_inputs()
    realization = KokunoSourceCompatiblePartitionRealization(grid_origin=origin)
    base = _evaluate(realization, inputs)
    ladders = {
        d: [_fd6(realization, inputs, base, d, h) for h in FD6_STEPS]
        for d in ("r", "z")
    }
    refinement = {
        d: [ladders[d][i]["relative_rms"] / max(ladders[d][i+1]["relative_rms"], 1e-300)
            for i in range(2)]
        for d in ("r", "z")
    }

    eta = np.asarray(base["eta"], dtype=float)
    j_weight = int(np.argmax(np.mean(eta * eta, axis=0)))
    eta_mut = eta.copy(); eta_mut[..., j_weight] = 0.0
    dropped_partition = float(np.max(np.abs(np.sum(eta_mut * eta_mut, axis=-1) - 1.0)))

    # Independent finest radial FD6 reused only for mutation calibration.
    h = FD6_STEPS[-1]
    shifted = {k: _evaluate(realization, inputs, k*h, "r") for k in (-3,-2,-1,1,2,3)}
    union = _union_labels([base, *shifted.values()])
    f = {k: _align(v["eta"], v["beta_labels"], union) for k, v in shifted.items()}
    fd6 = (-f[-3] + 9*f[-2] - 45*f[-1] + 45*f[1] - 9*f[2] + f[3]) / (60*h)
    Dr = np.asarray(base["D_r_eta"], dtype=float)
    j_derivative = int(np.argmax(np.mean(Dr * Dr, axis=0)))
    Dr_mut = Dr.copy(); Dr_mut[..., j_derivative] = 0.0
    Dr_ref = _align(Dr, base["beta_labels"], union)
    dropped_derivative = _rms(fd6 - _align(Dr_mut, base["beta_labels"], union)) / max(_rms(Dr_ref), 1e-300)

    partition_max = float(np.max(np.abs(base["partition_error"])))
    radial_closure = float(np.max(np.abs(base["radial_closure"])))
    axial_closure = float(np.max(np.abs(base["axial_closure"])))
    mesh_error = max(abs(float(v) - float(k)**-6) for k, v in base["mesh_by_ell"].items())
    checks = {
        "squared_partition_closure": partition_max <= PARTITION_MAX_ABS_GUARD,
        "radial_differentiated_closure": radial_closure <= DERIVATIVE_CLOSURE_MAX_ABS_GUARD,
        "axial_differentiated_closure": axial_closure <= DERIVATIVE_CLOSURE_MAX_ABS_GUARD,
        "mesh_formula_exact_to_float": mesh_error == 0.0,
        "public_beta_labels_unique": len(set(base["beta_labels"])) == len(base["beta_labels"]),
        "source_overlap_bound": _active_ell_span(base) <= 2,
        "radial_fd6_fine_relative_rms": ladders["r"][-1]["relative_rms"] <= FD6_FINE_RELATIVE_RMS_GUARD,
        "axial_fd6_fine_relative_rms": ladders["z"][-1]["relative_rms"] <= FD6_FINE_RELATIVE_RMS_GUARD,
        "radial_fd6_refinement": min(refinement["r"]) >= FD6_MIN_REFINEMENT_RATIO_GUARD,
        "axial_fd6_refinement": min(refinement["z"]) >= FD6_MIN_REFINEMENT_RATIO_GUARD,
        "dropped_label_detected": dropped_partition >= DROPPED_LABEL_PARTITION_ERROR_FLOOR,
        "dropped_derivative_detected": dropped_derivative >= DROPPED_DERIVATIVE_RELATIVE_ERROR_FLOOR,
    }
    return {
        "grid_origin": list(origin),
        "held_out_points": N_HELD_OUT,
        "base_beta_count": len(base["beta_labels"]),
        "ell_labels": [int(x) for x in base["ell_labels"]],
        "partition_max_abs": partition_max,
        "radial_closure_max_abs": radial_closure,
        "axial_closure_max_abs": axial_closure,
        "mesh_formula_max_abs_error": float(mesh_error),
        "maximum_active_ell_span": _active_ell_span(base),
        "fd6_ladders": ladders,
        "fd6_relative_rms_refinement_ratios": refinement,
        "mutation": {
            "dropped_weight_label": repr(base["beta_labels"][j_weight]),
            "dropped_label_partition_error_max_abs": dropped_partition,
            "dropped_derivative_label": repr(base["beta_labels"][j_derivative]),
            "dropped_derivative_relative_error_vs_finest_fd6": float(dropped_derivative),
        },
        "checks": checks,
        "local_preflight_passed": bool(all(checks.values())),
    }


def build_report() -> dict[str, Any]:
    reports = [_origin_report(origin) for origin in ORIGIN_CASES]
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "dependency": {"agent2_parent_pr": BASE_PR, "agent2_exact_head": BASE_HEAD,
                       "public_interface": "KokunoSourceCompatiblePartitionRealization.evaluate"},
        "preregistration": {
            "seed": SEED, "held_out_points": N_HELD_OUT, "fd6_steps": list(FD6_STEPS),
            "origin_cases": [list(x) for x in ORIGIN_CASES],
            "thresholds": {
                "partition_max_abs": PARTITION_MAX_ABS_GUARD,
                "differentiated_closure_max_abs": DERIVATIVE_CLOSURE_MAX_ABS_GUARD,
                "fd6_finest_relative_rms": FD6_FINE_RELATIVE_RMS_GUARD,
                "minimum_fd6_refinement_ratio": FD6_MIN_REFINEMENT_RATIO_GUARD,
                "dropped_label_partition_error_min": DROPPED_LABEL_PARTITION_ERROR_FLOOR,
                "dropped_derivative_relative_error_min": DROPPED_DERIVATIVE_RELATIVE_ERROR_FLOOR,
            },
            "thresholds_changed_after_results": False,
        },
        "independent_operator": {
            "value_path": "public evaluate(...)[eta] only at shifted caller coordinates",
            "derivative_oracle": "centered seven-point sixth-order finite difference",
            "beta_alignment": "public label identity; absent shifted supports are zero",
            "construction_raw_bump_helper_used": False,
            "construction_normalization_helper_used": False,
            "construction_analytic_derivative_helper_used": False,
            "agent2_fd4_test_operator_reused": False,
            "training_tensor_or_loss_read": False,
            "free_forcing_used": False,
        },
        "origin_reports": reports,
        "local_guards": {
            "source_compatible_partition_independent_preflight_passed": bool(
                all(v for report in reports for v in report["checks"].values())
            ),
            "scope": "autonomous slow-partition implementation seam only",
        },
        "formal_project_gates": {
            "normalized_momentum_max_l2": FORMAL_MOMENTUM_GATE,
            "divergence_max_l2": FORMAL_DIVERGENCE_GATE,
            "formal_full_domain_pde_gate_assessed": False,
        },
        "truth_boundary": {
            "source_actual_partition_recovered": False,
            "source_hidden_grid_origin_recovered": False,
            "actual_positive_order_background_instantiated": False,
            "actual_source_public_oscillatory_xyz_t_velocity_available": False,
            "genuinely_independent_second_covariance_column_available": False,
            "global_leading_velocity_pressure_available": False,
            "after_correction_global_velocity_available": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def write_report(path: str | Path) -> Path:
    target = Path(path)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing audit report: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(build_report(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/kokuno_agent4/source_compatible_partition_independent_audit.json")
    args = parser.parse_args(argv)
    path = write_report(args.output)
    print(path.read_text(encoding="utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
