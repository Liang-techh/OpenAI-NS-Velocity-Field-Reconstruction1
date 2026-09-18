"""Log-space bridge from the selected Appendix-B normalization to outer scales.

KokunoYumeto's corrected 2026-09-09 reconstruction uses one fixed constant
``C>1`` in both

    F = phi / C

and the PA.11/PA.15 outer scale

    X_R = 110 (C P_*)^10.

The source chooses ``C`` only after several finite lower bounds are available;
in particular it requires ``C >= sup |phi_*|``, the PA.10 separation bound
``(C P_*)^10 > exp(T_sh+8)``, and additional source bounds.  The current
source-rescaled core/reference path has an executable *real-axis* normalization
``log_C_real_axis`` while the older finite outer schedule still stores the
unrelated diagnostic value ``C=2``.  Reusing that old ``log_X_R`` with the new
rescaled Appendix-B path would therefore mix two different C normalizations.

This module makes the shared-C dependency executable without materializing the
enormous selected ``C``.  It computes the conditional source-form outer scales
entirely in log space and exposes a hard compatibility guard against the old
finite ``C``.  The selected real-axis normalization is still an autonomous
numerical realization: the source complex-neighborhood C certificate and the
remaining PA.11 lower bounds are *not* proved here.  Consequently this bridge
is a prerequisite for selected-pressure PA.15 moments/PA.16 matching, not a
completed global leading profile or PDE validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_rescaled_appendix_b_boundary import KokunoSourceRescaledAppendixBBoundary
from .kokuno_rescaled_appendix_b_selected import (
    make_source_scale_aware_appendix_b_boundary,
)


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-rescaled-outer-scale-binding-v1"
LOG_X_RESTORE_START = -8.0

_SOURCE_FORMULAS = {
    "shared_C": "one fixed C>1; F=phi/C",
    "PA11_C_choice": (
        "choose C once with C>=sup|phi_*|, C>=C_0, "
        "C>=(R_*/110)^(1/10)/P_*, (C P_*)^10>exp(T_sh+8), "
        "and X_R=110(C P_*)^10>14 exp(8)/(3 Q_min)"
    ),
    "PA15_outer_scale": "X_R=110(C P_*)^10; x=X/X_R",
    "outer_power_start": "X_w=X_R*exp(T_d+2)",
    "reserved_intervals": (
        "I1=X_w*exp((T_w-25,T_w-20)); "
        "I2=X_w*exp((T_w-20,T_w-15)); "
        "I3=X_w*exp((T_w-14,T_w-9)); "
        "I4=X_w*exp((T_w-8,T_w-3))"
    ),
    "PA15_moment_scaling": (
        "(M,I,J,S,C_p)=(X_R*Mhat,X_R^(3/2)*Ihat,"
        "X_R^(3/2)*Jhat,X_R*Shat,Cphat)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "shared_C_dependency_recorded": True,
    "source_outer_X_R_formula_executable_in_log_space": True,
    "selected_real_axis_log_C_executable": True,
    "old_finite_outer_C_compatibility_guard_executable": True,
    "PA15_log_scaling_factors_executable": True,
    "selected_real_axis_C_is_autonomous_normalization": True,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
    "source_hidden_C_recovered": False,
    "actual_source_incoming_five_moment_discrepancy_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoRescaledOuterScaleBinding:
    """Bind the selected rescaled ``log C`` to the public outer scale formulas.

    The nested :class:`KokunoOuterReservedPatchSchedule` supplies only the
    already-explicit autonomous choices ``M_d``, ``P_*`` margin, ``lambda``,
    ``h`` and reserved-patch offset.  Its finite ``C`` is *not* reused for the
    rescaled path unless it exactly matches ``boundary.log_C``.
    """

    boundary: KokunoSourceRescaledAppendixBBoundary = field(
        default_factory=make_source_scale_aware_appendix_b_boundary
    )

    def __post_init__(self) -> None:
        if not isinstance(self.boundary, KokunoSourceRescaledAppendixBBoundary):
            raise TypeError(
                "boundary must be a KokunoSourceRescaledAppendixBBoundary"
            )
        schedule = self.template_outer_schedule
        if not math.isclose(
            float(schedule.h), float(self.boundary.h), rel_tol=0.0, abs_tol=1.0e-14
        ):
            raise ValueError("boundary and outer schedule must use the same h")
        if not math.isfinite(self.selected_log_C) or self.selected_log_C <= 0.0:
            raise ValueError("selected log C must be positive and finite")

    @property
    def template_outer_schedule(self) -> KokunoOuterReservedPatchSchedule:
        return self.boundary.reference.seed.binding.outer_schedule

    @property
    def selected_log_C(self) -> float:
        return float(self.boundary.log_C)

    @property
    def template_log_C(self) -> float:
        return math.log(float(self.template_outer_schedule.C))

    @property
    def finite_template_C_matches_selected_C(self) -> bool:
        return math.isclose(
            self.selected_log_C,
            self.template_log_C,
            rel_tol=0.0,
            abs_tol=8.0 * np.finfo(float).eps * max(1.0, abs(self.selected_log_C)),
        )

    @property
    def selected_C_materializable_binary64(self) -> bool:
        return self.selected_log_C <= math.log(np.finfo(float).max)

    @property
    def log_X_R(self) -> float:
        schedule = self.template_outer_schedule
        return math.log(110.0) + 10.0 * (
            self.selected_log_C + float(schedule.log_P_star)
        )

    @property
    def template_log_X_R(self) -> float:
        return float(self.template_outer_schedule.log_X_R)

    @property
    def delta_log_X_R_from_finite_template(self) -> float:
        # Written in this form to expose the exact source of the scale shift.
        return 10.0 * (self.selected_log_C - self.template_log_C)

    @property
    def log_x_i(self) -> float:
        return math.log(110.0) - self.log_X_R

    @property
    def log_X_w(self) -> float:
        schedule = self.template_outer_schedule
        return self.log_X_R + float(schedule.T_d) + 2.0

    @property
    def log_X_star(self) -> float:
        schedule = self.template_outer_schedule
        return self.log_X_w + float(schedule.T_w) + float(schedule.repair_log_offset)

    @property
    def log_e_w(self) -> float:
        # This amplitude depends on P_*, T_d and lambda, but not on C.
        return float(self.template_outer_schedule.log_e_w)

    @property
    def log_c_patch(self) -> float:
        schedule = self.template_outer_schedule
        return self.log_e_w + (0.5 + float(schedule.lambda_outer)) * self.log_X_w

    @property
    def log_e_star(self) -> float:
        # Evaluate after cancellation instead of subtracting two O(log X_R)
        # values.  This stays accurate even when selected log C is enormous.
        schedule = self.template_outer_schedule
        return self.log_e_w - (0.5 + float(schedule.lambda_outer)) * (
            float(schedule.T_w) + float(schedule.repair_log_offset)
        )

    def reserved_log_intervals(self) -> dict[str, tuple[float, float]]:
        base = self.log_X_w + float(self.template_outer_schedule.T_w)
        return {
            "I1": (base - 25.0, base - 20.0),
            "I2": (base - 20.0, base - 15.0),
            "I3": (base - 14.0, base - 9.0),
            "I4": (base - 8.0, base - 3.0),
        }

    def pa15_log_inverse_scales(self) -> dict[str, float]:
        """Return logs of factors mapping physical moments to normalized rows."""

        return {
            "M": -self.log_X_R,
            "I": -1.5 * self.log_X_R,
            "J": -1.5 * self.log_X_R,
            "S": -self.log_X_R,
            "C_p": 0.0,
        }

    def log_x_sep(self, T_sh: float) -> float:
        value = float(T_sh)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("T_sh must be finite and nonnegative")
        return math.log(110.0) + value - self.log_X_R

    def separation_geometry_feasible(self, T_sh: float) -> bool:
        return self.log_x_sep(T_sh) < LOG_X_RESTORE_START

    def require_finite_template_C_match(self) -> None:
        """Fail closed before old finite-C PA.15 scaling is reused by mistake."""

        if not self.finite_template_C_matches_selected_C:
            raise ValueError(
                "the selected rescaled C does not match the finite outer schedule C; "
                "use this log-space binding for X_R/PA.15 scales"
            )

    def scale_report(self) -> dict[str, Any]:
        template_intervals = self.template_outer_schedule.reserved_log_intervals()
        selected_intervals = self.reserved_log_intervals()
        return {
            "selected_log_C": self.selected_log_C,
            "template_log_C": self.template_log_C,
            "finite_template_C_matches_selected_C": self.finite_template_C_matches_selected_C,
            "selected_C_materializable_binary64": self.selected_C_materializable_binary64,
            "log_P_star": float(self.template_outer_schedule.log_P_star),
            "selected_log_X_R": self.log_X_R,
            "finite_template_log_X_R": self.template_log_X_R,
            "delta_log_X_R": self.delta_log_X_R_from_finite_template,
            "selected_log_x_i": self.log_x_i,
            "selected_log_X_w": self.log_X_w,
            "selected_reserved_log_intervals": {
                key: list(value) for key, value in selected_intervals.items()
            },
            "finite_template_reserved_log_intervals": {
                key: list(value) for key, value in template_intervals.items()
            },
            "selected_log_X_star": self.log_X_star,
            "selected_log_e_w": self.log_e_w,
            "selected_log_c_patch": self.log_c_patch,
            "selected_log_e_star_stable": self.log_e_star,
            "pa15_log_inverse_scales": self.pa15_log_inverse_scales(),
            "source_complex_C_bound_verified": False,
            "source_all_PA11_C_bounds_verified": False,
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "boundary": self.boundary.to_payload(),
            "selected_scale_report": self.scale_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoRescaledOuterScaleBinding":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected rescaled outer-scale binding schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("rescaled outer-scale source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("rescaled outer-scale truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("rescaled outer-scale payload SHA-256 mismatch")
        obj = cls(
            boundary=KokunoSourceRescaledAppendixBBoundary.from_payload(
                payload.get("boundary")
            )
        )
        if obj.to_payload() != payload:
            raise ValueError("rescaled outer-scale replay changed payload")
        return obj

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoRescaledOuterScaleBinding":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit the selected rescaled-C outer log-scale binding"
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    binding = KokunoRescaledOuterScaleBinding()
    report = {
        "schema": SCHEMA,
        "binding_sha256": binding.sha256,
        "scale_report": binding.scale_report(),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
