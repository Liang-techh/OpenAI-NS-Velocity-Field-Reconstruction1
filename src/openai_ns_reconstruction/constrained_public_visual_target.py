"""Truth-bounded public visual target for the OpenAI Navier--Stokes figure.

This module records only qualitative facts that are observable in OpenAI's
public 2026-09-08 Navier--Stokes article/caption.  It deliberately contains no
hidden frame time, camera calibration, streamline seeds, numerical field data,
or fitted visual threshold.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from types import MappingProxyType
from typing import Any, Iterable, Mapping


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class PublicObservable:
    observable_id: str
    source_scope: str
    kind: str
    statement: str


@dataclass(frozen=True)
class PublicVisualTarget:
    target_id: str
    publisher: str
    published_date: str
    source_url: str
    retrieved_date: str
    source_type: str
    source_license_status: str
    observables: tuple[PublicObservable, ...]
    explicitly_unobserved: tuple[str, ...]

    @property
    def observable_ids(self) -> tuple[str, ...]:
        return tuple(item.observable_id for item in self.observables)


OPENAI_NS_PUBLIC_VISUAL_TARGET = PublicVisualTarget(
    target_id="openai_ns_public_visual_2026-09-08_v1",
    publisher="OpenAI",
    published_date="2026-09-08",
    source_url="https://openai.com/index/navier-stokes-solution/",
    retrieved_date="2026-09-18",
    source_type="public_webpage_text_figure_and_caption",
    source_license_status=(
        "No redistribution license is asserted here; only short factual/qualitative "
        "observations are recorded, with no copied image or source implementation."
    ),
    observables=(
        PublicObservable(
            "inward_spiraling",
            "figure_caption_and_article_text",
            "qualitative_geometry",
            "The displayed trajectories spiral inward.",
        ),
        PublicObservable(
            "axial_stretching",
            "figure_caption",
            "qualitative_geometry",
            "The displayed trajectories exhibit axial stretching.",
        ),
        PublicObservable(
            "vortex_swirl",
            "article_text",
            "qualitative_geometry",
            "The public description identifies the flow as a spinning vortex/swirl.",
        ),
        PublicObservable(
            "increasing_elongation",
            "article_text",
            "qualitative_time_trend",
            "The vortex is described as becoming increasingly elongated over the dynamics.",
        ),
        PublicObservable(
            "central_region_shrinks",
            "article_text",
            "qualitative_time_trend",
            "The central region is described as shrinking over the dynamics.",
        ),
        PublicObservable(
            "speed_increases",
            "article_text",
            "qualitative_time_trend",
            "The flow is described as speeding up while the central region shrinks.",
        ),
        PublicObservable(
            "relative_angular_rotation_color_encoding",
            "figure_caption",
            "visual_encoding",
            "Orange marks faster angular rotation and teal marks slower rotation.",
        ),
        PublicObservable(
            "circulating_speed_depends_on_radius",
            "figure_caption",
            "qualitative_kinematics",
            "Circulating speed is stated to depend on radius.",
        ),
    ),
    explicitly_unobserved=(
        "numerical_velocity_samples",
        "physical_coordinate_scale",
        "numerical_frame_time",
        "camera_pose_or_projection",
        "streamline_seed_positions",
        "hidden_model_parameters",
        "quantitative_visual_acceptance_thresholds",
        "paper_exact_coefficients",
    ),
)


_FORBIDDEN_EVIDENCE_KEYS = frozenset(
    {
        "target_value",
        "target_threshold",
        "acceptance_threshold",
        "hidden_frame_time",
        "camera_pose",
        "camera_registration",
        "streamline_seed_target",
        "recovered_velocity",
        "recovered_parameters",
    }
)
_TRUTH_KEYS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def public_visual_target() -> PublicVisualTarget:
    """Return the frozen public-observable target contract."""

    return OPENAI_NS_PUBLIC_VISUAL_TARGET


def _nonempty_text(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def audit_visual_evidence_receipt(
    *,
    candidate_sha256: str,
    evidence: Iterable[Mapping[str, Any]],
    target: PublicVisualTarget = OPENAI_NS_PUBLIC_VISUAL_TARGET,
) -> Mapping[str, Any]:
    """Audit coverage of public observables without creating a visual pass/fail gate.

    Evidence may report candidate measurements, but it cannot introduce numerical
    target values, hidden-time/camera recovery, or truth-state promotion.  The
    returned mapping is immutable at the top level and is explicitly diagnostic.
    """

    if not isinstance(candidate_sha256, str) or not _SHA256_RE.fullmatch(candidate_sha256):
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")

    allowed_ids = set(target.observable_ids)
    seen: set[str] = set()
    measured: list[str] = []
    inconclusive: list[str] = []
    not_measured: list[str] = []

    for index, raw in enumerate(evidence):
        if not isinstance(raw, Mapping):
            raise ValueError(f"evidence[{index}] must be a mapping")
        forbidden = _FORBIDDEN_EVIDENCE_KEYS.intersection(raw)
        if forbidden:
            raise ValueError(
                "public-observable evidence may not contain hidden/numerical target keys: "
                + ", ".join(sorted(forbidden))
            )
        for truth_key in _TRUTH_KEYS:
            if raw.get(truth_key) is True:
                raise ValueError(f"evidence may not promote {truth_key}")

        observable_id = _nonempty_text(raw.get("observable_id"), name="observable_id")
        if observable_id not in allowed_ids:
            raise ValueError(f"unknown public observable: {observable_id}")
        if observable_id in seen:
            raise ValueError(f"duplicate public observable evidence: {observable_id}")
        seen.add(observable_id)

        status = _nonempty_text(raw.get("status"), name="status")
        if status not in {"measured", "inconclusive", "not_measured"}:
            raise ValueError("status must be measured, inconclusive, or not_measured")
        if status != "not_measured":
            _nonempty_text(raw.get("method"), name="method")
            _nonempty_text(raw.get("provenance"), name="provenance")

        if status == "measured":
            measured.append(observable_id)
        elif status == "inconclusive":
            inconclusive.append(observable_id)
        else:
            not_measured.append(observable_id)

    missing = [observable_id for observable_id in target.observable_ids if observable_id not in seen]
    coverage_fraction = len(measured) / len(target.observable_ids)
    result = {
        "schema_version": 1,
        "task_id": "CR-A9-037",
        "claim_scope": "public_observable_coverage_only",
        "candidate_sha256": candidate_sha256,
        "target_id": target.target_id,
        "source_url": target.source_url,
        "measured_observables": tuple(measured),
        "inconclusive_observables": tuple(inconclusive),
        "declared_not_measured": tuple(not_measured),
        "missing_observables": tuple(missing),
        "measured_coverage_fraction": coverage_fraction,
        "posthoc_visual_threshold_defined": False,
        "hidden_data_recovery_attempted": False,
        **{key: False for key in _TRUTH_KEYS},
    }
    return MappingProxyType(result)
