"""Endpoint replay for the frozen disjoint global collar correction.

The collar coefficients are frozen from ``global_collar_tangent.json``.  This
script evaluates the same order-13 disjoint-box grid at ``k0 + 1e-6`` with
the frozen Cartesian FD steps, compares the parent global field with the
collar-corrected field, and checks that the reference-time compact wave is
unchanged.  It is a bounded exterior replay, not a whole-space or recursion
test.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

for _name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from affine_momentum import jets, momentum  # noqa: E402
from global_axial_extension import build_candidate  # noqa: E402
from global_collar_tangent import BOXES, CollarCorrection, grid, metrics  # noqa: E402


COLLAR_PATH = ROOT / "global_collar_tangent.json"
OUTPUT_PATH = ROOT / "global_collar_time_replay.json"
SNAPSHOT_PATH = ROOT / "full_wave_frozen_cache.json"
GLOBAL_SOURCE_PATH = ROOT / "global_axial_extension.py"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _wave_support(snapshot):
    wave = snapshot["inputs"]["wave"]
    center = np.asarray(wave["center"], dtype=float)
    widths = np.asarray(wave["widths"], dtype=float)
    radial = np.asarray([center[0] - widths[0], center[0] + widths[0]])
    axial = np.asarray([center[1] - widths[1], center[1] + widths[1]])
    points = np.asarray(
        [[r, 0.0, z] for r in radial for z in axial]
        + [[center[0], 0.0, center[1]]],
        dtype=float,
    )
    return wave, center, widths, radial, axial, points


def _support_disjointness(boxes, radial, axial):
    rows = []
    wave_r0, wave_r1 = map(float, radial)
    wave_z0, wave_z1 = map(float, axial)
    for index, box in enumerate(boxes):
        r0, r1, z0, z1 = map(float, box)
        radial_gap = max(wave_r0 - r1, r0 - wave_r1, 0.0)
        axial_gap = max(wave_z0 - z1, z0 - wave_z1, 0.0)
        separated_radially = bool(r1 <= wave_r0 or r0 >= wave_r1)
        separated_axially = bool(z1 <= wave_z0 or z0 >= wave_z1)
        disjoint = bool(separated_radially or separated_axially)
        rows.append(
            {
                "box_index": int(index),
                "box": [r0, r1, z0, z1],
                "separated_radially": separated_radially,
                "separated_axially": separated_axially,
                "radial_gap": float(radial_gap),
                "axial_gap": float(axial_gap),
                "disjoint": disjoint,
            }
        )
    return {
        "wave_radial_interval": [wave_r0, wave_r1],
        "wave_axial_interval": [wave_z0, wave_z1],
        "boxes": rows,
        "all_boxes_disjoint": bool(all(row["disjoint"] for row in rows)),
    }


def run(output_path=OUTPUT_PATH):
    started = time.perf_counter()
    output_path = Path(output_path)
    frozen = json.loads(COLLAR_PATH.read_text(encoding="utf-8"))
    if frozen.get("status") != "completed":
        raise ValueError("global_collar_tangent.json is not completed")
    control = np.asarray(frozen["control"], dtype=float)
    if control.shape != (108,):
        raise ValueError(f"Expected 108 collar controls, got {control.shape}")
    full, localized, mean, candidate, snapshot, mean_report = build_candidate()
    tau0 = float(snapshot["inputs"]["mean"]["tau"])
    k0 = float(snapshot["inputs"]["mean"]["k"])
    delta_k = 1.0e-6
    k_endpoint = k0 + delta_k
    tau_endpoint = 0.5 * 2.0 ** (-k_endpoint)
    hspace = float(snapshot["timesteps"]["hspace"])
    htime = float(snapshot["timesteps"]["htime"])
    points, weights = grid(13)
    corrected = CollarCorrection(full, control, tau0)

    parent_endpoint = momentum(jets(full, points, tau_endpoint, hspace, htime))
    corrected_endpoint = momentum(jets(corrected, points, tau_endpoint, hspace, htime))
    wave, center, widths, radial, axial, wave_points = _wave_support(snapshot)
    disjointness = _support_disjointness(BOXES, radial, axial)

    # At tau0 the collar velocity term is zero.  Because each box is disjoint
    # from the compact wave rectangle, both velocity and pressure remain the
    # original frozen wave values at these support probes.
    reference_parent, reference_pressure = full.fields(wave_points, tau0)
    reference_corrected, reference_corrected_pressure = corrected.fields(wave_points, tau0)
    reference_wave_check = {
        "tau": tau0,
        "point_count": int(len(wave_points)),
        "points": wave_points.tolist(),
        "velocity_max_abs_difference": float(
            np.max(np.abs(reference_corrected - reference_parent))
        ),
        "pressure_max_abs_difference": float(
            np.max(np.abs(reference_corrected_pressure - reference_pressure))
        ),
        "unchanged": bool(
            np.max(np.abs(reference_corrected - reference_parent)) <= 1.0e-12
            and np.max(np.abs(reference_corrected_pressure - reference_pressure)) <= 1.0e-12
        ),
    }
    report = {
        "status": "completed",
        "accepted": False,
        "pde_validated": False,
        "scale_recursion_established": False,
        "finite_time_evolution_validated": False,
        "scope": (
            "Endpoint replay of the frozen three-box global collar correction on "
            "the same order-13 collider grid and finite-difference steps used by "
            "global_collar_tangent. This is an exterior sampled diagnostic only; "
            "no whole-space, PDE, or scale-recursion acceptance is claimed."
        ),
        "sources": {
            "collar_candidate": {
                "path": COLLAR_PATH.name,
                "sha256": _sha256(COLLAR_PATH),
            },
            "global_source": {
                "path": GLOBAL_SOURCE_PATH.name,
                "sha256": _sha256(GLOBAL_SOURCE_PATH),
            },
            "frozen_cache": {
                "path": SNAPSHOT_PATH.name,
                "sha256": _sha256(SNAPSHOT_PATH),
            },
            "mean_source": {
                "path": mean_report.get("source_dynamic_report", "broad_shear_dynamic_control.json"),
                "coefficients": mean_report.get("coefficients"),
            },
        },
        "reference": {
            "k": k0,
            "tau": tau0,
            "frozen_report_metrics": frozen.get("independent_replay"),
            "training_metrics": frozen.get("training"),
        },
        "endpoint": {
            "k": k_endpoint,
            "tau": tau_endpoint,
            "delta_k": delta_k,
            "physical_time": "t=-tau",
            "grid_order": 13,
            "grid_point_count": int(len(points)),
            "finite_difference_steps": {"hspace": hspace, "htime": htime},
            "parent_global": metrics(parent_endpoint, weights),
            "collar_corrected": metrics(corrected_endpoint, weights),
            "difference_metrics": metrics(corrected_endpoint - parent_endpoint, weights),
        },
        "support_disjointness": disjointness,
        "reference_wave_unchanged": reference_wave_check,
        "limitations": {
            "domain": "Only the three registered disjoint collar boxes are sampled; no whole-R3 estimate is made.",
            "time": "The endpoint is a single instantaneous replay at k0 + 1e-6, not a time-integrated trajectory.",
            "acceptance": "No Navier--Stokes closure, finite-energy theorem, PDE, or scale-recursion acceptance is asserted.",
        },
        "elapsed_seconds": float(time.perf_counter() - started),
    }
    output_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "output": str(output_path),
                "endpoint_parent_l2": report["endpoint"]["parent_global"]["momentum_volume_L2"],
                "endpoint_corrected_l2": report["endpoint"]["collar_corrected"]["momentum_volume_L2"],
                "support_disjoint": disjointness["all_boxes_disjoint"],
                "reference_wave_unchanged": reference_wave_check["unchanged"],
                "elapsed_seconds": report["elapsed_seconds"],
            }
        ),
        flush=True,
    )
    return report


if __name__ == "__main__":
    run()
