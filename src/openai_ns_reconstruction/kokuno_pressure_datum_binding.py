"""Bind the public Appendix-B axis-pressure datum to the selected outer scale.

This module is a narrow Kokuno Agent-1 bridge.  The corrected 2026-09-09
reconstruction requires, on the axis,

    Pi_0(eta) <= -(5/2) P_*^2 f(eta)^2,
    f(eta)=(1+eta^2)^(-1).

The repository's earlier finite B.13/core reference instead used an autonomous
``pressure_scale=1`` in

    Pi_0(eta)=-pressure_scale^2 f(eta)^2.

For that restricted reference ansatz, satisfying the public inequality is
equivalent to

    pressure_scale >= sqrt(5/2) P_*.

The source does not publish a unique hidden value.  This file therefore makes
one explicit autonomous source-compatible choice: the lower bound (optionally
multiplied by a user-visible margin), rounded one floating-point step upward.
It exposes a callable B.13 reference velocity and a fail-closed probe of whether
the existing finite nonlinear core/reference continuation can carry that scale.

A successful datum check is not global pressure matching.  A failed finite-core
probe is a numerical representation barrier, not evidence that the source
construction is inconsistent.  No forcing, PDE threshold, ST006 evidence, or
paper-exact state is changed here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import warnings

import numpy as np

from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate
from .kokuno_appendix_b_boundary import KokunoAppendixBBoundary
from .paper_core_reference import PaperCoreReference


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pressure-datum-binding-v1"

_SOURCE_FORMULAS = {
    "outer_f": "f(eta)=(1+eta^2)^(-1)",
    "source_axis_pressure_bound": "Pi_0 <= -(5/2) P_*^2 f(eta)^2",
    "restricted_reference_datum": "Pi_0=-pressure_scale^2 f(eta)^2",
    "restricted_reference_lower_bound": "pressure_scale >= sqrt(5/2) P_*",
    "outer_P_star": "log(P_*)=T_d+log_p_star_margin in the selected outer schedule",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_axis_pressure_bound_formula_executable": True,
    "selected_source_compatible_pressure_datum_executable": True,
    "selected_B13_reference_velocity_executable": True,
    "selected_pressure_scale_is_autonomous_existence_choice": True,
    "source_hidden_pressure_scale_recovered": False,
    "source_outer_pressure_datum_bound_recovered": False,
    "actual_source_appendix_B_trajectory_reconstructed": False,
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


def _finite_array(value: Any, name: str) -> np.ndarray:
    arr = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def _log10_abs_max(values: Any) -> float | None:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or np.any(~np.isfinite(arr)):
        return None
    maximum = float(np.max(np.abs(arr)))
    if maximum == 0.0:
        return None
    return math.log10(maximum)


@dataclass(frozen=True)
class KokunoPressureDatumBinding:
    """Selected source-compatible realization of the public axis-pressure bound."""

    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    pressure_margin: float = 1.0

    # Keep the existing executable B.13/core-reference choices explicit.  Only
    # pressure_scale is replaced by the source-compatible selected value.
    j0: float = 0.02
    sigma: float = 0.5
    Lambda: float = 10.0
    log_transition_width: float = 0.005
    quadrature_points: int = 24

    def __post_init__(self) -> None:
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        margin = float(self.pressure_margin)
        if not math.isfinite(margin) or not 1.0 <= margin <= 10.0:
            raise ValueError("pressure_margin must lie in [1,10]")
        j0 = float(self.j0)
        sigma = float(self.sigma)
        Lambda = float(self.Lambda)
        width = float(self.log_transition_width)
        if not math.isfinite(j0) or not 0.0 < j0 <= 0.05:
            raise ValueError("j0 must lie in (0,0.05]")
        if not math.isfinite(sigma) or sigma <= 0.0:
            raise ValueError("sigma must be positive")
        if not math.isfinite(Lambda) or Lambda < 1.0:
            raise ValueError("Lambda must be at least 1")
        if not math.isfinite(width) or not 0.0 < width < 0.02:
            raise ValueError("log_transition_width must lie in (0,0.02)")
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        q = int(self.quadrature_points)
        if not 8 <= q <= 64:
            raise ValueError("quadrature_points must lie in [8,64]")
        object.__setattr__(self, "pressure_margin", margin)
        object.__setattr__(self, "j0", j0)
        object.__setattr__(self, "sigma", sigma)
        object.__setattr__(self, "Lambda", Lambda)
        object.__setattr__(self, "log_transition_width", width)
        object.__setattr__(self, "quadrature_points", q)

    @property
    def P_star(self) -> float:
        value = math.exp(float(self.outer_schedule.log_P_star))
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("selected outer P_* is outside positive float range")
        return value

    @property
    def required_pressure_scale(self) -> float:
        return math.sqrt(2.5) * self.P_star

    @property
    def selected_pressure_scale(self) -> float:
        # The source inequality is non-strict.  One upward float step avoids a
        # false negative from squaring/rounding at equality without inventing a
        # scientifically meaningful margin.
        selected = self.pressure_margin * self.required_pressure_scale
        return math.nextafter(selected, math.inf)

    @property
    def pressure_ratio_to_lower_bound(self) -> float:
        required = self.required_pressure_scale
        return self.selected_pressure_scale / required

    @property
    def pressure_square_ratio_to_source_bound(self) -> float:
        return (self.selected_pressure_scale / self.required_pressure_scale) ** 2

    @staticmethod
    def f(eta: Any) -> np.ndarray:
        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        return 1.0 / (1.0 + values * values)

    def selected_axis_pressure(self, eta: Any) -> np.ndarray:
        f = self.f(eta)
        return -(self.selected_pressure_scale**2) * f * f

    def source_axis_pressure_upper_bound(self, eta: Any) -> np.ndarray:
        f = self.f(eta)
        return -2.5 * (self.P_star**2) * f * f

    def pressure_bound_satisfied(self, eta: Any) -> np.ndarray:
        return self.selected_axis_pressure(eta) <= self.source_axis_pressure_upper_bound(eta)

    def make_reference(self) -> PaperCoreReference:
        """Return the executable B.13 reference with the selected pressure datum."""
        return PaperCoreReference(
            h=float(self.outer_schedule.h),
            j=self.j0,
            sigma=self.sigma,
            Lambda=self.Lambda,
            C=float(self.outer_schedule.C),
            pressure_scale=self.selected_pressure_scale,
        )

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Vectorized B.13 ``F,U,Pi`` for the selected pressure datum."""
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        ref = self.make_reference()
        F = np.empty_like(X_array)
        U = np.empty_like(X_array)
        Pi = np.empty_like(X_array)
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            F.reshape(-1)[index] = ref.F(float(xx), float(ee))
            U.reshape(-1)[index] = ref.U(float(xx), float(ee))
            Pi.reshape(-1)[index] = ref.Pi(float(xx), float(ee))
        return {"F": F, "U": U, "Pi": Pi}

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Vectorized Cartesian B.13 reference velocity ``[...,3]``.

        This is only the near-axis finite reference.  It is not the nonlinear
        source core, Appendix-B join, or global Kokuno leading field.
        """
        x_arr, y_arr, z_arr, t_arr = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            _finite_array(t, "t"),
        )
        if np.any((t_arr < 0.0) | (t_arr >= 1.0)):
            raise ValueError("time must satisfy 0<=t<1")
        ref = self.make_reference()
        out = np.empty(x_arr.shape + (3,), dtype=float)
        flat = out.reshape((-1, 3))
        for index, (xx, yy, zz, tt) in enumerate(
            zip(
                x_arr.reshape(-1),
                y_arr.reshape(-1),
                z_arr.reshape(-1),
                t_arr.reshape(-1),
                strict=True,
            )
        ):
            value = np.asarray(
                ref.velocity(float(xx), float(yy), float(zz), float(tt)),
                dtype=float,
            )
            if value.shape != (3,) or np.any(~np.isfinite(value)):
                raise RuntimeError("selected B.13 reference velocity became non-finite")
            flat[index] = value
        return out

    __call__ = velocity

    def make_core_candidate(
        self, *, maxdegree: int = 14, eta_nodes: int = 257
    ) -> KokunoLeadingCoreSeriesCandidate:
        """Build the existing finite nonlinear core at the selected pressure scale.

        This method intentionally propagates the datum through the real existing
        recurrence.  Any dynamic-range failure is surfaced to the caller.
        """
        return KokunoLeadingCoreSeriesCandidate(
            h=float(self.outer_schedule.h),
            j0=self.j0,
            sigma=self.sigma,
            Lambda=self.Lambda,
            C=float(self.outer_schedule.C),
            pressure_scale=self.selected_pressure_scale,
            maxdegree=maxdegree,
            eta_nodes=eta_nodes,
            quadrature_points=self.quadrature_points,
        )

    def series_probe(self, *, maxdegree: int = 14, eta_nodes: int = 65) -> dict[str, Any]:
        """Fail-closed binary64 feasibility probe for the existing nonlinear core."""
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                candidate = self.make_core_candidate(
                    maxdegree=int(maxdegree), eta_nodes=int(eta_nodes)
                )
                series = candidate.series
                report = {
                    "status": "finite",
                    "maxdegree": int(maxdegree),
                    "eta_nodes": int(eta_nodes),
                    "log10_max_abs_u_coefficient": _log10_abs_max(series.u_coefficients),
                    "log10_max_abs_pi_coefficient": _log10_abs_max(series.pi_coefficients),
                    "candidate_sha256": candidate.sha256,
                }
        except Exception as exc:  # numerical barrier is part of the receipt
            report = {
                "status": "nonfinite_or_failed",
                "maxdegree": int(maxdegree),
                "eta_nodes": int(eta_nodes),
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            }
        report["scientific_interpretation"] = (
            "finite binary64 replay of the existing truncated core"
            if report["status"] == "finite"
            else "binary64/current-recurrence barrier only; not source infeasibility"
        )
        return report

    def make_reference_continuation(
        self, *, maxdegree: int = 14, eta_nodes: int = 257
    ) -> KokunoReferenceContinuationCandidate:
        return KokunoReferenceContinuationCandidate(
            h=float(self.outer_schedule.h),
            j0=self.j0,
            sigma=self.sigma,
            Lambda=self.Lambda,
            C=float(self.outer_schedule.C),
            pressure_scale=self.selected_pressure_scale,
            maxdegree=maxdegree,
            eta_nodes=eta_nodes,
            quadrature_points=self.quadrature_points,
            log_transition_width=self.log_transition_width,
        )

    def appendix_b_probe(
        self, *, eta: float = 0.0, maxdegree: int = 14, eta_nodes: int = 65
    ) -> dict[str, Any]:
        """Attempt one actual selected Appendix-B boundary solve and report honestly."""
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                continuation = self.make_reference_continuation(
                    maxdegree=int(maxdegree), eta_nodes=int(eta_nodes)
                )
                boundary = KokunoAppendixBBoundary(reference=continuation)
                values = boundary.boundary_values(np.asarray([eta], dtype=float))
                result = {
                    key: float(np.asarray(value).reshape(-1)[0])
                    for key, value in values.items()
                }
            status = {
                "status": "finite",
                "eta": eta,
                "maxdegree": int(maxdegree),
                "eta_nodes": int(eta_nodes),
                "boundary_values": result,
                "reference_sha256": continuation.sha256,
            }
        except Exception as exc:
            status = {
                "status": "nonfinite_or_failed",
                "eta": eta,
                "maxdegree": int(maxdegree),
                "eta_nodes": int(eta_nodes),
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            }
        status["scientific_interpretation"] = (
            "selected source-compatible datum propagated through the existing finite Appendix-B path"
            if status["status"] == "finite"
            else "numerical/current-finite-model barrier; global source join remains unclaimed"
        )
        return status

    def datum_report(self) -> dict[str, Any]:
        sample_eta = np.asarray([-1.0, -0.5, 0.0, 0.5, 1.0])
        checks = self.pressure_bound_satisfied(sample_eta)
        return {
            "log_P_star": float(self.outer_schedule.log_P_star),
            "P_star": self.P_star,
            "required_pressure_scale": self.required_pressure_scale,
            "selected_pressure_scale": self.selected_pressure_scale,
            "pressure_margin": self.pressure_margin,
            "pressure_ratio_to_lower_bound": self.pressure_ratio_to_lower_bound,
            "pressure_square_ratio_to_source_bound": self.pressure_square_ratio_to_source_bound,
            "sample_eta": sample_eta.tolist(),
            "sample_bound_satisfied": np.asarray(checks, dtype=bool).tolist(),
            "all_sample_bounds_satisfied": bool(np.all(checks)),
            "algebraic_bound_satisfied_by_construction": bool(
                self.selected_pressure_scale >= self.required_pressure_scale
            ),
            "baseline_pressure_scale_before_binding": 1.0,
            "selected_to_previous_baseline_ratio": self.selected_pressure_scale,
            "hidden_source_value_recovered": False,
        }

    def report(
        self,
        *,
        probe_maxdegree: int | None = None,
        probe_eta_nodes: int = 65,
        appendix_b_eta: float | None = None,
    ) -> dict[str, Any]:
        result: dict[str, Any] = {
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "datum": self.datum_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
        if probe_maxdegree is not None:
            result["finite_core_probe"] = self.series_probe(
                maxdegree=int(probe_maxdegree), eta_nodes=int(probe_eta_nodes)
            )
        if appendix_b_eta is not None:
            result["appendix_B_probe"] = self.appendix_b_probe(
                eta=float(appendix_b_eta),
                maxdegree=int(probe_maxdegree if probe_maxdegree is not None else 14),
                eta_nodes=int(probe_eta_nodes),
            )
        return result

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
            "parameters": {
                "pressure_margin": self.pressure_margin,
                "j0": self.j0,
                "sigma": self.sigma,
                "Lambda": self.Lambda,
                "log_transition_width": self.log_transition_width,
                "quadrature_points": self.quadrature_points,
            },
            "dependencies": {"outer_schedule": self.outer_schedule.to_payload()},
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode()).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPressureDatumBinding":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno pressure-datum binding schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("pressure-datum source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("pressure-datum truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode()).hexdigest()
        if claimed != expected:
            raise ValueError("pressure-datum payload SHA mismatch")
        params = dict(payload.get("parameters", {}))
        outer = KokunoOuterReservedPatchSchedule.from_payload(
            payload.get("dependencies", {}).get("outer_schedule")
        )
        obj = cls(outer_schedule=outer, **params)
        if obj.to_payload() != payload:
            raise ValueError("pressure-datum payload does not replay exactly")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoPressureDatumBinding":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pressure-margin", type=float, default=1.0)
    parser.add_argument("--probe-maxdegree", type=int, default=14)
    parser.add_argument("--probe-eta-nodes", type=int, default=65)
    parser.add_argument("--appendix-b-eta", type=float, default=None)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    obj = KokunoPressureDatumBinding(pressure_margin=args.pressure_margin)
    report = obj.report(
        probe_maxdegree=args.probe_maxdegree,
        probe_eta_nodes=args.probe_eta_nodes,
        appendix_b_eta=args.appendix_b_eta,
    )
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
