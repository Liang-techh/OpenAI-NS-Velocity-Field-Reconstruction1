"""Frozen behavioral binding for the exact ST052-M parent evaluator.

The package-level ST052 morphology transforms accept an externally supplied
``base_velocity(points, time)`` callable.  Caller-provided candidate/head strings
are not enough to authenticate that callable.  This module migrates the frozen
reference probes already embedded in the exact ST052-M reconstruction recipe and
uses them as a fail-closed runtime behavioral binding.

This is intentionally narrower than full parent materialization: passing ten
source-native probes does not prove equality of two functions on the whole
space-time domain and is not a parent artifact checksum.  It is a concrete
runtime guard that rejects mismatched injected callables and can be composed
with the package transform identity while full parent save/load remains pending.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np

Array = np.ndarray
BaseVelocity = Callable[[Array, float], Array]


class St052ParentBindingError(ValueError):
    """Raised when an injected parent callable fails the frozen source probes."""


@dataclass(frozen=True)
class St052ParentReferenceSpec:
    """Source-native ST052-M velocity reference migrated from PR #508 recipe."""

    schema_version: int = 1
    binding_kind: str = "source_native_behavioral_reference_v1"
    parent_candidate_id: str = "ST052-M"
    parent_source_pr: int = 508
    parent_source_head: str = "b3b8bfdbe1077f9ec967d158602951997d81e17d"
    source_recipe_path: str = "experiments/root_st052/recipe.json"
    source_recipe_git_blob_sha1: str = "e30c769052379f72afeee46ca264482884cc5ac7"
    atol: float = 1.0e-8
    rtol: float = 1.0e-8
    points: tuple[tuple[float, float, float], ...] = (
        (0.1, 0.0, 0.1),
        (0.0, 0.0, 0.2),
        (0.0, 0.0, 0.0),
        (0.7, 0.0, 1.55),
        (1.4, 0.0, -0.1),
        (0.02, 0.0, 1.8),
        (-0.3, 0.5, 0.9),
        (0.0, 0.3, 0.0),
        (0.57, 0.0, 1.92),
        (0.3, 0.0, -1.91),
    )
    times: tuple[float, ...] = (
        0.25,
        0.371,
        0.5,
        0.633,
        0.75,
        0.25,
        0.687,
        0.52,
        0.75,
        0.25,
    )
    velocity: tuple[tuple[float, float, float], ...] = (
        (-0.0032184221683559306, 0.006423453107365135, 0.008068000014928186),
        (-0.0, 0.0, 0.01690318903930133),
        (-0.0, 0.0, 0.001458334118599456),
        (0.24170754974178593, 0.11057825847128837, 0.060790620984724056),
        (-0.03897196815729013, 0.0644574398842995, 0.02742232167082357),
        (-0.00013463236030421128, -6.478183237001746e-05, -0.0010794508897124186),
        (-0.22716143336195221, -0.2207036714880389, 0.3745305256308158),
        (-0.05644376828773303, -0.026354784624148798, 0.0012078791358357771),
        (0.001138406029274419, 5.086012311515379e-05, 1.4364739520658454e-05),
        (0.001151910654119985, -0.00033354306377833624, -2.4675171172199502e-05),
    )

    def payload(self) -> dict:
        return asdict(self)

    def sha256(self) -> str:
        raw = json.dumps(
            self.payload(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        return hashlib.sha256(raw).hexdigest()

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(
                {"spec": self.payload(), "spec_sha256": self.sha256()},
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load(cls, path: str | Path) -> "St052ParentReferenceSpec":
        obj = json.loads(Path(path).read_text(encoding="utf-8"))
        payload = dict(obj["spec"])
        payload["points"] = tuple(tuple(row) for row in payload["points"])
        payload["times"] = tuple(payload["times"])
        payload["velocity"] = tuple(tuple(row) for row in payload["velocity"])
        spec = cls(**payload)
        if obj.get("spec_sha256") != spec.sha256():
            raise St052ParentBindingError("ST052 parent reference checksum mismatch")
        return spec


DEFAULT_PARENT_REFERENCE = St052ParentReferenceSpec()


def verify_parent_callable(
    base_velocity: BaseVelocity,
    spec: St052ParentReferenceSpec = DEFAULT_PARENT_REFERENCE,
) -> dict:
    """Verify an injected parent callable on the frozen source-native probes.

    The tolerance exactly matches the ``np.allclose(..., atol=1e-8, rtol=1e-8)``
    velocity check used by the exact ST052-M recipe replay.  Failure raises
    :class:`St052ParentBindingError`; success returns a serializable receipt.
    """
    points = np.asarray(spec.points, dtype=float)
    times = np.asarray(spec.times, dtype=float)
    expected = np.asarray(spec.velocity, dtype=float)
    actual = np.empty_like(expected)

    for time in np.unique(times):
        mask = times == time
        value = np.asarray(base_velocity(points[mask], float(time)), dtype=float)
        if value.shape != expected[mask].shape:
            raise St052ParentBindingError(
                f"parent callable returned shape {value.shape}, expected {expected[mask].shape}"
            )
        if not np.isfinite(value).all():
            raise St052ParentBindingError("parent callable returned non-finite values")
        actual[mask] = value

    diff = np.abs(actual - expected)
    tolerance = spec.atol + spec.rtol * np.abs(expected)
    passed = bool(np.all(diff <= tolerance))
    max_abs = float(np.max(diff))
    max_scaled = float(np.max(diff / tolerance))
    if not passed:
        index = np.unravel_index(int(np.argmax(diff / tolerance)), diff.shape)
        raise St052ParentBindingError(
            "ST052 parent behavioral reference mismatch: "
            f"index={index}, max_abs={max_abs:.17g}, max_scaled={max_scaled:.17g}"
        )

    return {
        "binding_kind": spec.binding_kind,
        "parent_candidate_id": spec.parent_candidate_id,
        "parent_source_head": spec.parent_source_head,
        "source_recipe_git_blob_sha1": spec.source_recipe_git_blob_sha1,
        "reference_spec_sha256": spec.sha256(),
        "probe_count": int(len(points)),
        "atol": spec.atol,
        "rtol": spec.rtol,
        "max_abs": max_abs,
        "max_scaled_error": max_scaled,
        "passed": True,
        "full_parent_function_equality_proved": False,
        "full_parent_artifact_materialized": False,
    }


def behavioral_binding_id(
    spec: St052ParentReferenceSpec = DEFAULT_PARENT_REFERENCE,
) -> str:
    """Return an immutable ID for this reference contract, not for the full field."""
    payload = {
        "kind": spec.binding_kind,
        "parent_candidate_id": spec.parent_candidate_id,
        "parent_source_head": spec.parent_source_head,
        "reference_spec_sha256": spec.sha256(),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


class St052ReferenceBoundParent:
    """Callable wrapper that verifies parent behavior before allowing evaluation."""

    def __init__(
        self,
        base_velocity: BaseVelocity,
        spec: St052ParentReferenceSpec = DEFAULT_PARENT_REFERENCE,
    ) -> None:
        self.spec = spec
        self.receipt = verify_parent_callable(base_velocity, spec)
        self.binding_id = behavioral_binding_id(spec)
        self._base_velocity = base_velocity

    def __call__(self, points: Array, time: float) -> Array:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be finite with shape (n,3)")
        time = float(time)
        if not np.isfinite(time) or not (0.25 - 1.0e-12 <= time <= 0.75 + 1.0e-12):
            raise ValueError("time outside frozen ST052 interval")
        value = np.asarray(self._base_velocity(points, time), dtype=float)
        if value.shape != points.shape or not np.isfinite(value).all():
            raise ValueError("parent callable must return finite shape-(n,3) values")
        return value

    def save_receipt(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "behavioral_binding_id": self.binding_id,
            "verification": self.receipt,
            "truth_boundary": TRUTH_BOUNDARY,
        }
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )


TRUTH_BOUNDARY = {
    "runtime_parent_behavior_bound_on_frozen_reference": True,
    "full_st052_parent_materialized_here": False,
    "parent_artifact_sha256_assigned": False,
    "whole_domain_parent_equivalence_proved": False,
    "complete_child_candidate_materialized": False,
    "complete_candidate_save_load_ready": False,
    "production_candidate_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}
