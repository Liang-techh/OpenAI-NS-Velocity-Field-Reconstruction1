"""Precision-qualified total Cartesian composition of the late relative-swirl edit.

This module stacks exactly on A1 #1171, which preserves Kokuno's public
post-flattening relative-swirl correction as representable binary64 base and
delta channels because the current O(lambda^28) relative edit is below ordinary
binary64 relative epsilon.

The public reconstruction structure remains the same as #1171:

    E = E_unedited * (1 + c1 beta1 + c2 beta2),

with the public angular endpoint condition and exact pressure-neutrality row.
No source coefficient, pointwise bump, pressure, forcing, residual or threshold
is changed here.  The only new operation is numerical representation: each
binary64 base/delta channel is embedded into a fixed, non-user-tunable Decimal
context and added there.  This yields one public

    velocity(x, y, z, t) -> [..., 3]

surface that retains the mathematically nonzero #1171 correction instead of
silently rounding it away.  The returned ndarray has dtype=object and Decimal
entries.  It is intentionally a precision-qualified stage representation, not
a claim that ordinary binary64 delivery, terminal/global leading, or PDE
validation is complete.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_relative_swirl_split import (
    KokunoPA16CurrentCartesianRelativeSwirlSplit,
)

SCHEMA = "kokuno-pa16-current-cartesian-relative-swirl-decimal-composed-v1"
PARENT_EXACT_HEAD = "657818dd83e119e6091bd2490d3804c04c3ef723"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
DECIMAL_DIGITS = 96
PRECISION_MARGIN_DIGITS = 24

_SOURCE_FORMULAS = {
    "relative_edit": "E=E_unedited*(1+c1*beta1+c2*beta2)",
    "postpulse_U": "U=0 permanently after the pulse endpoint",
    "angular_target": "I=XH/(1-lambda) after the two relative-swirl bumps",
    "pressure_neutrality": "integral(E^2-E_unedited^2)dy=0",
    "cartesian_velocity": (
        "u_r=(v0/(2q))r; u_theta=q^(-A-1/2) r F; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1171 binary64 base+delta split channels",
    "total_representation": (
        "fixed 96-significant-decimal-digit object-array composition; "
        "precision is not a public/caller tuning knob"
    ),
    "binary64_embedding": (
        "Decimal.from_float followed by the frozen 96-digit context before addition"
    ),
    "rounding": "ROUND_HALF_EVEN",
    "pointwise_bump_shape": (
        "repository-autonomous standard C-infinity bump inherited unchanged from #1154/#1171"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "public_relative_swirl_two_bump_algebra_materialized": True,
    "absolute_current_r_I_at_flattening_endpoint_materialized": True,
    "current_relative_swirl_total_target_materialized": True,
    "relative_swirl_profile_delta_channel_materialized": True,
    "relative_swirl_cartesian_delta_channel_materialized": True,
    "precision_qualified_relative_swirl_total_materialized": True,
    "precision_qualified_public_velocity_materialized": True,
    "precision_qualified_decimal_json_export_materialized": True,
    "current_cartesian_relative_swirl_composed": True,
    "binary64_total_relative_swirl_sum_is_resolved": False,
    "binary64_total_velocity_export_ready": False,
    "source_exact_pointwise_bump_shape_recovered": False,
    "source_hidden_parameters_recovered": False,
    "source_exterior_heat_replacement_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _configure_context(ctx) -> None:
    ctx.prec = DECIMAL_DIGITS
    ctx.rounding = ROUND_HALF_EVEN
    ctx.Emax = 999999999
    ctx.Emin = -999999999


def _embed_binary64(value: float) -> Decimal:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("only finite binary64 values may enter Decimal composition")
    with localcontext() as ctx:
        _configure_context(ctx)
        return +Decimal.from_float(value)


def _compose_decimal(base: float, delta: float) -> tuple[Decimal, Decimal, Decimal]:
    """Return context-embedded base, delta, and their precision-qualified total."""
    with localcontext() as ctx:
        _configure_context(ctx)
        base_d = +Decimal.from_float(float(base))
        delta_d = +Decimal.from_float(float(delta))
        total_d = +(base_d + delta_d)
    return base_d, delta_d, total_d


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianRelativeSwirlDecimalComposed:
    """One public Decimal-valued velocity surface retaining exact #1171 delta."""

    split: KokunoPA16CurrentCartesianRelativeSwirlSplit = field(
        default_factory=KokunoPA16CurrentCartesianRelativeSwirlSplit,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.split, KokunoPA16CurrentCartesianRelativeSwirlSplit):
            raise TypeError("split must be KokunoPA16CurrentCartesianRelativeSwirlSplit")
        truth = self.split.truth_boundary
        if not truth["relative_swirl_cartesian_delta_channel_materialized"]:
            raise ValueError("exact #1171 Cartesian delta channel is required")
        if truth["current_cartesian_relative_swirl_composed"]:
            raise ValueError("parent unexpectedly already composes the total field")
        report = self.split.representation_report()
        max_edit = float(report["max_abs_relative_edit"])
        if not (0.0 < max_edit < np.finfo(float).eps):
            raise ValueError("#1171 sub-epsilon representation blocker is not present")
        required = int(math.ceil(-math.log10(max_edit))) + PRECISION_MARGIN_DIGITS
        if DECIMAL_DIGITS < required:
            raise ValueError(
                f"frozen Decimal precision {DECIMAL_DIGITS} is insufficient; need >= {required}"
            )

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.split.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    @property
    def geometry(self):
        return self.split.geometry

    @property
    def compensator(self):
        return self.split.compensator

    @property
    def log_X_flatten_end(self) -> float:
        return float(self.split.log_X_flatten_end)

    def velocity_decimal_split(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Embed exact #1171 public base/delta channels in the frozen Decimal context."""
        base, delta = self.split.velocity_split(x, y, z, t)
        base = np.asarray(base, dtype=float)
        delta = np.asarray(delta, dtype=float)
        if base.shape != delta.shape or base.shape[-1] != 3:
            raise RuntimeError("parent split velocity shape contract drifted")
        base_d = np.empty(base.shape, dtype=object)
        delta_d = np.empty(delta.shape, dtype=object)
        for idx in np.ndindex(base.shape):
            b, d, _ = _compose_decimal(base[idx], delta[idx])
            base_d[idx] = b
            delta_d[idx] = d
        return base_d, delta_d

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return one precision-qualified Decimal-valued Cartesian velocity array."""
        base, delta = self.split.velocity_split(x, y, z, t)
        base = np.asarray(base, dtype=float)
        delta = np.asarray(delta, dtype=float)
        total = np.empty(base.shape, dtype=object)
        for idx in np.ndindex(base.shape):
            _, _, total[idx] = _compose_decimal(base[idx], delta[idx])
        return total

    def profile_total_decimal_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Compose low-dimensional E/F and analytic D_logX derivatives in Decimal."""
        p = self.split.split_profile_logX(log_X, eta)
        shape = np.asarray(p["E_base"]).shape
        keys = {
            "E": (p["E_base"], p["delta_E"]),
            "F": (p["F_base"], p["delta_F"]),
            "U": (p["U_base"], p["delta_U"]),
            "DlogX_E": (
                -(0.5 + self.split.lambda_value) * np.asarray(p["E_base"], dtype=float),
                p["delta_E_DlogX"],
            ),
            "DlogX_F": (
                -(1.0 + self.split.lambda_value) * np.asarray(p["F_base"], dtype=float),
                p["delta_F_DlogX"],
            ),
            "DlogX_U": (np.zeros(shape, dtype=float), p["delta_U_DlogX"]),
        }
        out: dict[str, np.ndarray] = {}
        for name, (base, delta) in keys.items():
            ba = np.broadcast_to(np.asarray(base, dtype=float), shape)
            da = np.broadcast_to(np.asarray(delta, dtype=float), shape)
            arr = np.empty(shape, dtype=object)
            for idx in np.ndindex(shape):
                _, _, arr[idx] = _compose_decimal(ba[idx], da[idx])
            out[name] = arr
        return out

    def export_velocity_decimal_json(
        self, path: str | Path, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, Any]:
        """Deterministic stage export using exact decimal strings; not a global export claim."""
        velocity = self.velocity(x, y, z, t)
        payload = {
            "schema": "kokuno-relative-swirl-decimal-velocity-export-v1",
            "candidate_semantic_sha256": self.semantic_sha256,
            "decimal_digits": DECIMAL_DIGITS,
            "shape": list(velocity.shape),
            "values_row_major": [str(v) for v in velocity.reshape(-1)],
            "truth_boundary": self.truth_boundary,
        }
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    def representation_report(self) -> dict[str, Any]:
        parent_report = self.split.representation_report()
        log_X = self.log_X_flatten_end + self.compensator.y1
        radius = math.exp(0.5 * (math.log(2.0) + log_X))
        base, delta = self.split.velocity_split(radius, 0.0, 0.0, 0.0)
        total = self.velocity(radius, 0.0, 0.0, 0.0)
        base_d, delta_d = self.velocity_decimal_split(radius, 0.0, 0.0, 0.0)
        active = [idx for idx in np.ndindex(delta.shape) if float(delta[idx]) != 0.0]
        if not active:
            raise RuntimeError("representative relative-swirl Cartesian delta is unexpectedly zero")
        decimal_survival = []
        binary64_erasure = []
        exact_delta_replay = []
        with localcontext() as ctx:
            _configure_context(ctx)
            for idx in active:
                binary64_erasure.append(bool((base + delta)[idx] == base[idx]))
                decimal_survival.append(bool(total[idx] != base_d[idx]))
                exact_delta_replay.append(bool(+(total[idx] - base_d[idx]) == delta_d[idx]))
        max_edit = float(parent_report["max_abs_relative_edit"])
        required = int(math.ceil(-math.log10(max_edit))) + PRECISION_MARGIN_DIGITS
        return {
            "schema": "kokuno-agent1-relative-swirl-decimal-composition-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.split),
            "decimal_digits": DECIMAL_DIGITS,
            "required_decimal_digits_by_frozen_guard": required,
            "precision_margin_available": DECIMAL_DIGITS - required,
            "active_cartesian_components": len(active),
            "binary64_erases_all_active_components": bool(all(binary64_erasure)),
            "decimal_total_retains_all_active_components": bool(all(decimal_survival)),
            "decimal_total_minus_base_replays_delta": bool(all(exact_delta_replay)),
            "parent_max_abs_relative_edit": max_edit,
            "truth_boundary": self.truth_boundary,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "source_path": SOURCE_PATH,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "split_semantic_sha256": _semantic_value(self.split),
            "decimal_digits": DECIMAL_DIGITS,
            "precision_margin_digits": PRECISION_MARGIN_DIGITS,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]):
        obj = cls()
        if dict(payload) != obj.configuration():
            raise ValueError("configuration/provenance drift")
        return obj

    @classmethod
    def load_configuration(cls, path: str | Path):
        return cls.from_configuration(json.loads(Path(path).read_text()))
