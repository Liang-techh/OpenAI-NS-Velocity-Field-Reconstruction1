"""Selected Appendix-B prefix moments at ``X_i=110`` in PA.15 coordinates.

This module is a narrow bridge between two already executable reconstruction
stages: ``KokunoAppendixBBoundary`` supplies one explicitly selected
Appendix-B continuation through ``X_i=110`` and ``KokunoInnerJoinExit``
consumes five incoming moment discrepancies for the later PA.16 repair.

The corrected KokunoYumeto 2026-09-09 reconstruction retains exactly

    M   = int U dX,
    I   = int H dX,             H = sqrt(2X) E,
    J   = int U H dX,
    S   = int (U^2-E^2/2) dX,
    C_p = int E^2/(2X) dX,

and, for ``x=X/X_R``, the PA.15 scaling

    (M,I,J,S,C_p)
      = (X_R Mhat, X_R^(3/2) Ihat, X_R^(3/2) Jhat,
         X_R Shat, Cphat).

The selected discrepancy below is the executable #429 Appendix-B prefix at
``X_i`` minus the source ideal pair ``U_0=4 eta`` and
``E_0=P_* f(eta) x^(1/10)`` at the same endpoint.  The ideal prefix is
integrated analytically; the selected actual prefix is propagated with the
same Appendix-B ODE used by the boundary object.  This avoids materializing
the enormous ``X_R`` and avoids subtracting the singular-looking ideal
``C_p`` density pointwise near ``x=0``.

A critical truth boundary is explicit here.  The public Appendix-B axis datum
requires ``Pi_0 <= -(5/2) P_*^2 f^2``.  The current executable #429 boundary
inherits an autonomous finite ``pressure_scale`` from the earlier core
reference and is *not* yet bound to the outer ``P_*`` schedule.  Therefore the
quantity computed here is a selected autonomous PA.15 discrepancy, useful for
wiring/testing the source matching map, but it is not claimed to be the actual
source incoming discrepancy.  ``T_sh`` certification is also still separate.
No hidden Kokuno/OpenAI parameter is inferred and no completed global join or
PDE validation is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .kokuno_appendix_b_boundary import KokunoAppendixBBoundary, X_I
from .kokuno_inner_join_exit import KokunoInnerJoinExit, InnerJoinSolveResult
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-appendix-b-incoming-moments-v1"

_SOURCE_FORMULAS = {
    "five_prefix_moments": (
        "M=int U dX; I=int H dX; J=int U H dX; "
        "S=int(U^2-E^2/2)dX; C_p=int E^2/(2X)dX; H=sqrt(2X)E"
    ),
    "PA15_scaling": (
        "(M,I,J,S,C_p)=(X_R*Mhat,X_R^(3/2)*Ihat,"
        "X_R^(3/2)*Jhat,X_R*Shat,Cphat)"
    ),
    "ideal_pair": "U_0=4 eta; E_0=P_* f(eta) x^(1/10); f=(1+eta^2)^(-1)",
    "axis_pressure_datum_bound": "Pi_0 <= -(5/2) P_*^2 f(eta)^2",
    "incoming_discrepancy": (
        "selected Appendix-B five-prefix moments at X_i minus the ideal-pair "
        "five-prefix moments at x_i=X_i/X_R"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_appendix_B_prefix_moments_executable": True,
    "selected_PA15_scaled_incoming_discrepancy_executable": True,
    "selected_upstream_boundary_and_moment_data_bindable_to_PA16": True,
    "ideal_prefix_integrated_analytically": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_outer_pressure_datum_bound": False,
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
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


@dataclass(frozen=True)
class AppendixBIncomingMomentResult:
    """One scalar-eta receipt in PA.15 row order ``(M,I,J,S,C_p)``."""

    eta: float
    ell_i: float
    G_i: float
    actual_scaled_moments: tuple[float, float, float, float, float]
    ideal_scaled_moments: tuple[float, float, float, float, float]
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    max_abs_incoming_discrepancy: float


@dataclass(frozen=True)
class KokunoAppendixBIncomingMoments:
    """Propagate the selected Appendix-B prefix moments to ``X_i=110``."""

    boundary: KokunoAppendixBBoundary = field(default_factory=KokunoAppendixBBoundary)
    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    axis_quadrature_points: int = 96
    _axis_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _axis_weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.boundary, KokunoAppendixBBoundary):
            raise TypeError("boundary must be a KokunoAppendixBBoundary")
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if isinstance(self.axis_quadrature_points, bool) or not isinstance(
            self.axis_quadrature_points, (int, np.integer)
        ):
            raise TypeError("axis_quadrature_points must be an integer")
        order = int(self.axis_quadrature_points)
        if not 32 <= order <= 256:
            raise ValueError("axis_quadrature_points must lie in [32,256]")
        if not math.isclose(
            float(self.boundary.reference.C),
            float(self.outer_schedule.C),
            rel_tol=0.0,
            abs_tol=1e-14,
        ):
            raise ValueError("boundary C and outer-schedule C must agree")
        if not math.isclose(
            float(self.boundary.reference.h),
            float(self.outer_schedule.h),
            rel_tol=0.0,
            abs_tol=1e-14,
        ):
            raise ValueError("boundary h and outer-schedule h must agree")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "axis_quadrature_points", order)
        object.__setattr__(self, "_axis_nodes", nodes)
        object.__setattr__(self, "_axis_weights", weights)

    @property
    def log_x_i(self) -> float:
        return math.log(X_I) - float(self.outer_schedule.log_X_R)

    @property
    def x_i(self) -> float:
        value = math.exp(self.log_x_i)
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("x_i is outside positive float range; use log_x_i")
        return value

    @property
    def P_star(self) -> float:
        value = math.exp(float(self.outer_schedule.log_P_star))
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("P_* is outside positive float range")
        return value

    @property
    def required_source_pressure_scale_lower_bound(self) -> float:
        """Minimum scale implied by ``Pi_0<=-(5/2)P_*^2 f^2``.

        The earlier finite reference uses ``Pi_0=-pressure_scale^2 f^2`` at
        the axis, so the public source inequality would require
        ``pressure_scale >= sqrt(5/2) P_*`` for that restricted ansatz.
        """

        return math.sqrt(2.5) * self.P_star

    @property
    def selected_pressure_scale(self) -> float:
        return float(self.boundary.reference.pressure_scale)

    @property
    def source_outer_pressure_datum_bound(self) -> bool:
        return self.selected_pressure_scale >= self.required_source_pressure_scale_lower_bound

    def pressure_datum_report(self) -> dict[str, Any]:
        required = self.required_source_pressure_scale_lower_bound
        selected = self.selected_pressure_scale
        return {
            "selected_reference_pressure_scale": selected,
            "required_source_lower_bound_for_reference_ansatz": required,
            "selected_to_required_ratio": selected / required,
            "source_outer_pressure_datum_bound": self.source_outer_pressure_datum_bound,
            "actual_source_incoming_five_moment_discrepancy_bound": False,
        }

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X_i": X_I,
            "log_X_R": float(self.outer_schedule.log_X_R),
            "log_x_i": self.log_x_i,
            "x_i": self.x_i,
            "axis_quadrature_points": self.axis_quadrature_points,
            "pressure_datum": self.pressure_datum_report(),
            "source_T_sh_lower_bound_verified": False,
        }

    def _axis_actual_physical_moments(self, eta: float) -> np.ndarray:
        """Integrate the smooth selected/reference prefix on ``0<=X<=X0``."""

        X0 = float(self.boundary.X0)
        X = 0.5 * X0 * (self._axis_nodes + 1.0)
        weights = 0.5 * X0 * self._axis_weights
        F = np.asarray(self.boundary.reference.F(X, eta), dtype=float)
        U = np.asarray(self.boundary.reference.U(X, eta), dtype=float)
        if np.any(~np.isfinite(F)) or np.any(~np.isfinite(U)) or np.any(F <= 0.0):
            raise RuntimeError("reference axis profile produced invalid values")
        H = 2.0 * X * F
        densities = np.stack(
            [
                U,
                H,
                U * H,
                U * U - X * F * F,
                F * F,
            ],
            axis=0,
        )
        return np.asarray(densities @ weights, dtype=float)

    def _augmented_rhs(self, y: float, state: np.ndarray, eta: float, sl: Any) -> np.ndarray:
        base = self.boundary._rhs(float(y), np.asarray(state[:4], dtype=float), eta, sl)
        X = float(self.boundary.X0) * math.exp(float(y))
        F = math.exp(float(state[2]))
        U = float(state[3])
        H = 2.0 * X * F
        moment_rhs = np.asarray(
            [
                X * U,
                X * H,
                X * U * H,
                X * (U * U - X * F * F),
                X * F * F,
            ],
            dtype=float,
        )
        return np.concatenate((np.asarray(base, dtype=float), moment_rhs))

    def _selected_actual_physical_moments(self, eta: float) -> tuple[np.ndarray, float, float]:
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        sl = self.boundary._make_slice(eta)
        X0 = float(self.boundary.X0)
        F0 = float(self.boundary.reference.F(X0, eta))
        U0 = float(self.boundary.reference.U(X0, eta))
        p10 = -2.0 * X0 * float(self.boundary.reference.F_X(X0, eta)) / F0
        n0 = -2.0 * float(self.boundary.reference.U_X(X0, eta))
        initial_moments = self._axis_actual_physical_moments(eta)
        initial = np.concatenate(
            (
                np.asarray([p10, n0, math.log(F0), U0], dtype=float),
                initial_moments,
            )
        )
        solved = solve_ivp(
            lambda yy, zz: self._augmented_rhs(yy, zz, eta, sl),
            (0.0, float(self.boundary.y_i)),
            initial,
            method="DOP853",
            rtol=float(self.boundary.rtol),
            atol=float(self.boundary.atol),
            max_step=float(self.boundary.max_step),
        )
        if not solved.success or solved.y.size == 0:
            raise RuntimeError(f"Appendix-B moment propagation failed: {solved.message}")
        final = np.asarray(solved.y[:, -1], dtype=float)
        if final.shape != (9,) or np.any(~np.isfinite(final)):
            raise RuntimeError("Appendix-B moment propagation produced invalid state")
        F_i = math.exp(float(final[2]))
        ell_i = math.log(
            float(self.boundary.reference.C) * math.sqrt(2.0 * X_I) * F_i
        )
        G_i = float(final[3])
        return np.asarray(final[4:], dtype=float), ell_i, G_i

    def _actual_scaled_moments(self, physical: np.ndarray) -> np.ndarray:
        log_X_R = float(self.outer_schedule.log_X_R)
        scale_1 = math.exp(-log_X_R)
        scale_3_2 = math.exp(-1.5 * log_X_R)
        scaled = np.asarray(
            [
                physical[0] * scale_1,
                physical[1] * scale_3_2,
                physical[2] * scale_3_2,
                physical[3] * scale_1,
                physical[4],
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(scaled)):
            raise RuntimeError("PA.15 scaling produced nonfinite moments")
        return scaled

    def ideal_scaled_moments(self, eta: float) -> np.ndarray:
        """Exact ideal-pair prefix moments on ``0<=x<=x_i``."""

        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        x = self.x_i
        P = self.P_star
        f = 1.0 / (1.0 + eta * eta)
        I = math.sqrt(2.0) * P * f * x**1.6 / 1.6
        values = np.asarray(
            [
                4.0 * eta * x,
                I,
                4.0 * eta * I,
                16.0 * eta * eta * x
                - 0.5 * P * P * f * f * x**1.2 / 1.2,
                0.5 * P * P * f * f * x**0.2 / 0.2,
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(values)):
            raise RuntimeError("analytic ideal prefix produced nonfinite moments")
        return values

    def result_at_eta(self, eta: float) -> AppendixBIncomingMomentResult:
        physical, ell_i, G_i = self._selected_actual_physical_moments(float(eta))
        actual = self._actual_scaled_moments(physical)
        ideal = self.ideal_scaled_moments(float(eta))
        incoming = actual - ideal
        if np.any(~np.isfinite(incoming)):
            raise RuntimeError("incoming PA.15 discrepancy is nonfinite")
        return AppendixBIncomingMomentResult(
            eta=float(eta),
            ell_i=float(ell_i),
            G_i=float(G_i),
            actual_scaled_moments=tuple(map(float, actual)),
            ideal_scaled_moments=tuple(map(float, ideal)),
            incoming_scaled_discrepancy=tuple(map(float, incoming)),
            max_abs_incoming_discrepancy=float(np.max(np.abs(incoming))),
        )

    def incoming_scaled_discrepancy(self, eta: Any) -> np.ndarray:
        """Vectorized selected PA.15 discrepancy in ``(M,I,J,S,C_p)`` order."""

        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        out = np.empty(eta_array.shape + (5,), dtype=float)
        for index, value in enumerate(eta_array.reshape(-1)):
            out.reshape(-1, 5)[index, :] = np.asarray(
                self.result_at_eta(float(value)).incoming_scaled_discrepancy,
                dtype=float,
            )
        return out

    def solve_inner_join_at_eta(
        self,
        inner_join: KokunoInnerJoinExit,
        *,
        eta: float,
        absolute_tolerance: float = 2.0e-11,
    ) -> InnerJoinSolveResult:
        """Bind selected executable upstream data into the PA.16 inverse.

        This is a selected-realization software path only. ``inner_join.T_sh``
        remains caller supplied, the source pressure datum is not yet bound to
        outer ``P_*``, and neither condition is promoted by this method.
        """

        if not isinstance(inner_join, KokunoInnerJoinExit):
            raise TypeError("inner_join must be a KokunoInnerJoinExit")
        if inner_join.outer_schedule.to_payload() != self.outer_schedule.to_payload():
            raise ValueError("inner_join and incoming moments must use the same outer schedule")
        receipt = self.result_at_eta(float(eta))
        return inner_join.solve_at_eta(
            eta=float(eta),
            ell_i=receipt.ell_i,
            G_i=receipt.G_i,
            incoming_scaled_discrepancy=np.asarray(
                receipt.incoming_scaled_discrepancy, dtype=float
            ),
            absolute_tolerance=float(absolute_tolerance),
        )

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
            "parameters": {"axis_quadrature_points": self.axis_quadrature_points},
            "dependencies": {
                "boundary": self.boundary.to_payload(),
                "outer_schedule": self.outer_schedule.to_payload(),
            },
            "geometry": self.geometry_report(),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoAppendixBIncomingMoments":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Appendix-B incoming-moments schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("Appendix-B incoming-moment source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("Appendix-B incoming-moment truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
        if claimed != expected:
            raise ValueError("Appendix-B incoming-moment payload SHA mismatch")
        dependencies = payload.get("dependencies", {})
        parameters = payload.get("parameters", {})
        obj = cls(
            boundary=KokunoAppendixBBoundary.from_payload(dependencies.get("boundary")),
            outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(
                dependencies.get("outer_schedule")
            ),
            axis_quadrature_points=parameters.get("axis_quadrature_points"),
        )
        if obj.to_payload() != payload:
            raise ValueError("Appendix-B incoming-moment payload does not replay exactly")
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
    def load_json(cls, path: str | Path) -> "KokunoAppendixBIncomingMoments":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
