"""Independent local parameter-sensitivity audit for frozen velocity artifacts.

This module intentionally evaluates only the public ``VelocityField.at_points``
interface after loading a candidate artifact from disk.  It does not import any
optimizer/training loss and it does not classify Navier--Stokes validity.
"""
from __future__ import annotations

import copy
import json
import tempfile
from pathlib import Path
from typing import Iterable

import numpy as np

from .velocity_components import DEFAULT_CANDIDATE, VelocityField


def _validate_inputs(delta, indices, seed, held_out_points, times, box):
    if isinstance(delta, (bool, np.bool_)):
        raise TypeError("delta must be a positive finite real number")
    delta = float(delta)
    if not np.isfinite(delta) or delta <= 0.0 or delta > 0.25:
        raise ValueError("delta must satisfy 0 < delta <= 0.25")

    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise TypeError("seed must be an integer")
    seed = int(seed)

    if isinstance(held_out_points, (bool, np.bool_)) or int(held_out_points) != held_out_points:
        raise TypeError("held_out_points must be an integer")
    held_out_points = int(held_out_points)
    if held_out_points < 8:
        raise ValueError("held_out_points must be >= 8")

    indices = tuple(indices)
    if not indices:
        raise ValueError("at least one coefficient index is required")
    if any(isinstance(i, (bool, np.bool_)) or int(i) != i for i in indices):
        raise TypeError("coefficient indices must be integers")
    indices = tuple(int(i) for i in indices)
    if len(set(indices)) != len(indices):
        raise ValueError("coefficient indices must be unique")

    times = np.asarray(tuple(times), dtype=float)
    if times.ndim != 1 or not times.size or not np.all(np.isfinite(times)):
        raise ValueError("times must be a nonempty finite 1D sequence")
    if np.any((times < 0.25) | (times > 0.75)):
        raise ValueError("times must stay inside the packaged velocity interval [.25,.75]")

    box = np.asarray(box, dtype=float)
    if box.shape != (3, 2) or not np.all(np.isfinite(box)):
        raise ValueError("box must have finite shape (3,2)")
    if np.any(box[:, 1] <= box[:, 0]):
        raise ValueError("box upper bounds must exceed lower bounds")
    return delta, indices, seed, held_out_points, times, box


def _field_values(field, points, times):
    rows = []
    for time in times:
        value = np.asarray(field.at_points(points, float(time)), dtype=float)
        if value.shape != points.shape or not np.all(np.isfinite(value)):
            raise ValueError("public velocity interface returned malformed/nonfinite data")
        rows.append(value)
    return np.stack(rows, axis=0)


def _relative_change(reference, perturbed):
    diff = perturbed - reference
    vector_error = np.linalg.norm(diff, axis=-1)
    reference_speed = np.linalg.norm(reference, axis=-1)
    base_rms = float(np.sqrt(np.mean(reference_speed * reference_speed)))
    if not np.isfinite(base_rms) or base_rms <= np.finfo(float).eps:
        raise ValueError("reference velocity is numerically trivial on the held-out sample")
    rms = float(np.sqrt(np.mean(vector_error * vector_error)))
    return {
        "vector_error_rms": rms,
        "vector_error_max": float(np.max(vector_error)),
        "relative_rms_to_base": float(rms / base_rms),
        "base_velocity_rms": base_rms,
    }


def audit_velocity_parameter_perturbations(
    candidate_path=DEFAULT_CANDIDATE,
    *,
    delta: float = 0.01,
    coefficient_indices: Iterable[int] = (0, 58, 116),
    seed: int = 781223,
    held_out_points: int = 256,
    times: Iterable[float] = (0.30, 0.50, 0.70),
    box=((-1.75, 1.75), (-1.75, 1.75), (-1.75, 1.75)),
):
    """Measure local velocity sensitivity to bounded coefficient perturbations.

    Each variant is serialized to a temporary candidate artifact and then loaded
    through ``VelocityField``.  The same fixed held-out points are reused for the
    reference and every perturbation.  No optimization/training sample is used.

    The report deliberately contains no pass/fail threshold: local parameter
    sensitivity is evidence about robustness of the visualization artifact, not
    evidence that the PDE is satisfied.
    """
    delta, indices, seed, held_out_points, times, box = _validate_inputs(
        delta, coefficient_indices, seed, held_out_points, times, box
    )
    path = Path(candidate_path)
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("family") != "coupled_velocity_v1":
        raise ValueError("parameter perturbation audit currently requires coupled_velocity_v1")
    coefficients = np.asarray(payload.get("coefficients"), dtype=float)
    if coefficients.shape != (117,) or not np.all(np.isfinite(coefficients)):
        raise ValueError("candidate must contain 117 finite coupled coefficients")
    if np.any(np.abs(coefficients) > 1.0):
        raise ValueError("candidate coefficient bounds are violated")
    if any(i < 0 or i >= coefficients.size for i in indices):
        raise ValueError("coefficient index out of range")

    rng = np.random.default_rng(seed)
    points = rng.uniform(box[:, 0], box[:, 1], size=(held_out_points, 3))
    reference_field = VelocityField(path)
    reference = _field_values(reference_field, points, times)

    variants = []
    with tempfile.TemporaryDirectory(prefix="ns-velocity-perturb-") as tmp:
        root = Path(tmp)
        for index in indices:
            before = float(coefficients[index])
            direction = 1.0 if before + delta <= 1.0 else -1.0
            after = before + direction * delta
            if abs(after) > 1.0:
                raise ValueError("requested perturbation cannot remain inside coefficient bounds")

            variant_payload = copy.deepcopy(payload)
            variant_coefficients = list(variant_payload["coefficients"])
            variant_coefficients[index] = after
            variant_payload["coefficients"] = variant_coefficients
            variant_path = root / f"candidate_{index}.json"
            variant_path.write_text(json.dumps(variant_payload, indent=2) + "\n", encoding="utf-8")
            variant_field = VelocityField(variant_path)
            perturbed = _field_values(variant_field, points, times)
            variants.append(
                {
                    "coefficient_index": index,
                    "coefficient_before": before,
                    "coefficient_after": after,
                    "signed_delta": float(after - before),
                    **_relative_change(reference, perturbed),
                }
            )

    return {
        "candidate_path": str(path),
        "candidate_sha256": reference_field.sha256,
        "seed": seed,
        "held_out_points": held_out_points,
        "times": [float(t) for t in times],
        "box": box.tolist(),
        "coefficient_delta": delta,
        "variants": variants,
        "truth_boundary": (
            "held-out public-velocity parameter sensitivity only; not training loss, "
            "PDE validation, visual correspondence proof, paper-exact recovery, or blow-up evidence"
        ),
    }


def main():
    print(json.dumps(audit_velocity_parameter_perturbations(), indent=2))


if __name__ == "__main__":
    main()
