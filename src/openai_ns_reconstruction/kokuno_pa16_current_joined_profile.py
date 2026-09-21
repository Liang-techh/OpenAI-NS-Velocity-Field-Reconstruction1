"""Current-lineage executable Kokuno source-coordinate join through ``X_h``.

Pinned public provenance is the corrected 2026-09-09 KokunoYumeto
reconstruction at commit ``143f6773feb424ad9ed3a8d116653200f20346b7``.
The existing Agent-1 lineage already materializes the actual candidate-side
PA.10 profile through ``X_i=110`` and the current numerical PA.10 ``T_sh``
certificate.  ``KokunoInnerJoinExit`` already implements the public later join
in scaled radius ``x=X/X_R``:

* angular shaping from the Xi data to the ideal angular factor;
* axial restoration on ``-8 < log x < -7``;
* the two-U / three-E PA.16 five-moment repair on ``-6 < log x < -5``;
* the ideal-pair handoff at ``X_h=X_R exp(-5)``.

This module performs one missing binding only: it exposes that existing join as
one callable physical-X ``F(X,eta), U(X,eta), E(X,eta)`` profile on the current
Agent-1 identity, including the already-executable prefix ``0<=X<=X_i`` and
source-form first radial derivatives.  Inside the PA.16 repair interval it uses
the coefficients returned by the existing ``solve_selected_at_eta`` route; it
does not implement a second repair solver.

The selected ``T_sh`` is still based on finite eta sampling rather than the
source analytic C0 theorem, and the upstream pressure/stress/kappa choices are
still explicitly repository-autonomous.  Therefore this is a candidate-side
source-coordinate profile, not a source-admitted/global Cartesian velocity,
not matched pressure/forcing, and not Navier--Stokes validation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_inner_join_exit import LOG_JOIN_EXIT, LOG_REPAIR_START
from .kokuno_pa16_current_tsh_certificate import KokunoPA16CurrentTshCertificate
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)


SCHEMA = "kokuno-pa16-current-joined-profile-v1"

_SOURCE_FORMULAS = {
    "scaled_radius": "x=X/X_R; X_i=110; X_h=X_R*exp(-5)",
    "current_Xi_handoff": "G_i=U(X_i,eta); ell_i=log(C E(X_i,eta))",
    "angular_shaping": (
        "log E=-log C+y_i/10+(1-sigma(y_i/T_sh))*ell_i+"
        "sigma(y_i/T_sh)*log f"
    ),
    "axial_restoration": "restore G_i to 4 eta on -8<log x<-7",
    "PA16_repair": "two U bumps + three E bumps on -6<log x<-5",
    "ideal_exit": "U=4 eta; E=P_* f x^(1/10); f=(1+eta^2)^(-1)",
    "physical_F": "F(X,eta)=E(X,eta)/sqrt(2X)",
    "physical_radial_derivative": (
        "d/dX=(1/X_R)d/dx; F_X=E_X/sqrt(2X)-E/(2X)^(3/2)"
    ),
}

_NUMERICAL_REALIZATION = {
    "upstream": "exact current Agent-1 #953 ancestry",
    "T_sh": (
        "selected numerical lower bound from current finite-sample B0 envelope; "
        "not source analytic admission"
    ),
    "repair_coefficients": (
        "existing KokunoInnerJoinExit.solve_at_eta coefficients; no duplicate solver"
    ),
    "vectorization": (
        "broadcast physical X/eta; one existing PA.16 solve per distinct eta only "
        "when a requested point lies strictly inside the repair annulus"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_lineage_prefix_0_to_Xi_bound": True,
    "current_candidate_side_join_Xi_to_Xh_materialized": True,
    "current_candidate_side_PA16_coefficients_consumed": True,
    "physical_X_F_U_E_profile_executable_through_Xh": True,
    "physical_X_radial_derivatives_executable_through_Xh": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_prepared_upstream_five_moment_discrepancy_materialized": False,
    "source_global_inner_to_outer_join_admitted": False,
    "outer_global_leading_velocity_materialized": False,
    "cartesian_velocity_materialized": False,
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


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


@dataclass(frozen=True)
class KokunoPA16CurrentJoinedProfile:
    """Bind the current PA.10 prefix to the existing public PA.16 join."""

    certificate: KokunoPA16CurrentTshCertificate = field(
        default_factory=KokunoPA16CurrentTshCertificate,
        repr=False,
        compare=False,
    )
    quadrature_points: int = 96
    absolute_tolerance: float = 2.0e-11

    def __post_init__(self) -> None:
        if not isinstance(self.certificate, KokunoPA16CurrentTshCertificate):
            raise TypeError("certificate must be KokunoPA16CurrentTshCertificate")
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        q = int(self.quadrature_points)
        if not 32 <= q <= 256:
            raise ValueError("quadrature_points must lie in [32,256]")
        tol = float(self.absolute_tolerance)
        if not math.isfinite(tol) or not 0.0 < tol <= 1.0e-6:
            raise ValueError("absolute_tolerance must lie in (0,1e-6]")
        object.__setattr__(self, "quadrature_points", q)
        object.__setattr__(self, "absolute_tolerance", tol)

    @property
    def moments(self):
        return self.certificate.moments

    @property
    def X_i(self) -> float:
        return float(self.moments.X_i)

    @property
    def log_X_R(self) -> float:
        return float(self.moments.outer_schedule.log_X_R)

    @property
    def X_R(self) -> float:
        value = math.exp(self.log_X_R)
        if not math.isfinite(value):
            raise RuntimeError("current X_R is outside binary64 physical-X range")
        return value

    @property
    def log_X_h(self) -> float:
        return self.log_X_R + LOG_JOIN_EXIT

    @property
    def X_h(self) -> float:
        value = math.exp(self.log_X_h)
        if not math.isfinite(value):
            raise RuntimeError("current X_h is outside binary64 physical-X range")
        return value

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.moments.eta_interval)

    @property
    def route_ready(self) -> bool:
        return bool(self.certificate.separation_geometry_feasible)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_h)):
            raise ValueError(f"X must lie in [0,X_h], X_h={self.X_h:.17e}")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _selected_join(self):
        if not self.route_ready:
            raise ValueError(
                "current numerical T_sh route does not fit the source separation geometry"
            )
        return self.certificate.build_selected_inner_join(
            quadrature_points=self.quadrature_points
        )

    def _join_inputs(self, eta: float) -> dict[str, Any]:
        eta_f = float(eta)
        lo, hi = self.eta_interval
        if not math.isfinite(eta_f) or not lo <= eta_f <= hi:
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return self.moments.pa16_input_at_eta(eta_f)

    def _join_values_for_group(
        self, X: np.ndarray, *, eta: float
    ) -> dict[str, np.ndarray]:
        """Evaluate the existing later join for one eta, solving PA.16 only if needed."""
        join = self._selected_join()
        data = self._join_inputs(float(eta))
        x = np.asarray(X, dtype=float) / self.X_R
        log_x = np.log(x)
        needs_repair = bool(
            np.any((log_x > LOG_REPAIR_START) & (log_x < LOG_JOIN_EXIT))
        )
        coefficients = None
        if needs_repair:
            solved = self.certificate.solve_selected_at_eta(
                float(eta),
                quadrature_points=self.quadrature_points,
                absolute_tolerance=self.absolute_tolerance,
            )
            if not solved.repair.success:
                raise RuntimeError("existing PA.16 solve did not return success")
            coefficients = np.asarray(solved.repair.coefficients, dtype=float)
            if coefficients.shape != (5,) or np.any(~np.isfinite(coefficients)):
                raise RuntimeError("existing PA.16 solve returned invalid coefficients")

        scaled = join.profile_values_scaled(
            x,
            eta=float(eta),
            ell_i=float(data["ell_i"]),
            G_i=float(data["G_i"]),
            coefficients=coefficients,
        )
        U = np.asarray(scaled["U"], dtype=float)
        E = np.asarray(scaled["E"], dtype=float)
        U_X = np.asarray(scaled["U_x"], dtype=float) / self.X_R
        E_X = np.asarray(scaled["E_x"], dtype=float) / self.X_R
        root = np.sqrt(2.0 * np.asarray(X, dtype=float))
        F = E / root
        F_X = E_X / root - E / root**3
        H = root * E
        return {
            "F": F,
            "U": U,
            "E": E,
            "H": H,
            "F_X": F_X,
            "U_X": U_X,
            "E_X": E_X,
        }

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return broadcast current-lineage ``F,U,E,H`` on ``0<=X<=X_h``.

        The existing prefix evaluator is used on ``X<=X_i``.  For ``X>X_i``
        the existing later-join formulas are evaluated.  Distinct eta values
        are grouped so a PA.16 nonlinear solve, when required, is performed at
        most once per eta in this call.
        """
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F = np.empty_like(xf)
        U = np.empty_like(xf)
        E = np.empty_like(xf)
        H = np.empty_like(xf)

        inner = xf <= self.X_i
        if np.any(inner):
            vals = self.moments.profile_values(xf[inner], ef[inner])
            F[inner] = np.asarray(vals["F_actual_prefix"])
            U[inner] = np.asarray(vals["U_actual_prefix"])
            E[inner] = np.asarray(vals["E_actual_prefix"])
            H[inner] = np.sqrt(2.0 * xf[inner]) * E[inner]

        outer = ~inner
        if np.any(outer):
            if not self.route_ready:
                raise ValueError(
                    "X>X_i requested but current numerical T_sh geometry is infeasible"
                )
            outer_indices = np.flatnonzero(outer)
            for eta_value in np.unique(ef[outer]):
                local = outer_indices[ef[outer] == eta_value]
                vals = self._join_values_for_group(xf[local], eta=float(eta_value))
                F[local] = vals["F"]
                U[local] = vals["U"]
                E[local] = vals["E"]
                H[local] = vals["H"]

        if (
            np.any(~np.isfinite(F))
            or np.any(~np.isfinite(U))
            or np.any(~np.isfinite(E))
            or np.any(~np.isfinite(H))
        ):
            raise RuntimeError("current joined profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current joined profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_joined": F.reshape(shape),
            "U_current_joined": U.reshape(shape),
            "E_current_joined": E.reshape(shape),
            "H_current_joined": H.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return source-form first physical-X derivatives for ``X>0``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)

        inner = xf <= self.X_i
        if np.any(inner):
            vals = self.moments.radial_derivatives(xf[inner], ef[inner])
            F_X[inner] = np.asarray(vals["F_actual_prefix_X"])
            U_X[inner] = np.asarray(vals["U_actual_prefix_X"])
            E_X[inner] = np.asarray(vals["E_actual_prefix_X"])

        outer = ~inner
        if np.any(outer):
            if not self.route_ready:
                raise ValueError(
                    "X>X_i requested but current numerical T_sh geometry is infeasible"
                )
            outer_indices = np.flatnonzero(outer)
            for eta_value in np.unique(ef[outer]):
                local = outer_indices[ef[outer] == eta_value]
                vals = self._join_values_for_group(xf[local], eta=float(eta_value))
                F_X[local] = vals["F_X"]
                U_X[local] = vals["U_X"]
                E_X[local] = vals["E_X"]

        if np.any(~np.isfinite(F_X)) or np.any(~np.isfinite(U_X)) or np.any(~np.isfinite(E_X)):
            raise RuntimeError("current joined radial derivatives are non-finite")
        return {
            "F_current_joined_X": F_X.reshape(shape),
            "U_current_joined_X": U_X.reshape(shape),
            "E_current_joined_X": E_X.reshape(shape),
        }

    def Xi_handoff_report(self, eta: float) -> dict[str, Any]:
        """Compare the current prefix value/jet with the right-hand public join at Xi."""
        eta_f = float(eta)
        left_v = self.moments.profile_values(self.X_i, eta_f)
        left_d = self.moments.radial_derivatives(self.X_i, eta_f)
        if not self.route_ready:
            return {
                "eta": eta_f,
                "route_ready": False,
                "reason": "selected numerical T_sh does not fit source separation geometry",
            }
        right = self._join_values_for_group(np.asarray([self.X_i]), eta=eta_f)
        left = np.asarray(
            [
                float(left_v["F_actual_prefix"]),
                float(left_v["U_actual_prefix"]),
                float(left_v["E_actual_prefix"]),
            ]
        )
        right_values = np.asarray([right["F"][0], right["U"][0], right["E"][0]])
        left_jet = np.asarray(
            [
                float(left_d["F_actual_prefix_X"]),
                float(left_d["U_actual_prefix_X"]),
                float(left_d["E_actual_prefix_X"]),
            ]
        )
        right_jet = np.asarray([right["F_X"][0], right["U_X"][0], right["E_X"][0]])
        return {
            "eta": eta_f,
            "route_ready": True,
            "left_F_U_E": left.tolist(),
            "right_F_U_E": right_values.tolist(),
            "max_abs_value_jump": float(np.max(np.abs(left - right_values))),
            "left_FX_UX_EX": left_jet.tolist(),
            "right_FX_UX_EX": right_jet.tolist(),
            "max_abs_jet_jump": float(np.max(np.abs(left_jet - right_jet))),
        }

    def ideal_exit_report(self, eta: float) -> dict[str, Any]:
        eta_f = float(eta)
        if not self.route_ready:
            return {
                "eta": eta_f,
                "route_ready": False,
                "reason": "selected numerical T_sh does not fit source separation geometry",
            }
        vals = self._join_values_for_group(np.asarray([self.X_h]), eta=eta_f)
        f = 1.0 / (1.0 + eta_f * eta_f)
        E_ideal = math.exp(float(self.moments.outer_schedule.log_P_star)) * f * math.exp(-0.5)
        F_ideal = E_ideal / math.sqrt(2.0 * self.X_h)
        return {
            "eta": eta_f,
            "route_ready": True,
            "X_h": self.X_h,
            "U_exit": float(vals["U"][0]),
            "E_exit": float(vals["E"][0]),
            "F_exit": float(vals["F"][0]),
            "U_ideal": 4.0 * eta_f,
            "E_ideal": E_ideal,
            "F_ideal": F_ideal,
            "Xh_U_X": float(self.X_h * vals["U_X"][0]),
            "Xh_E_X_over_E": float(self.X_h * vals["E_X"][0] / vals["E"][0]),
            "Xh_F_X_over_F": float(self.X_h * vals["F_X"][0] / vals["F"][0]),
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "certificate": self.certificate.configuration(),
            "quadrature_points": self.quadrature_points,
            "absolute_tolerance": self.absolute_tolerance,
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA16CurrentJoinedProfile":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current joined-profile schema")
        return cls(
            certificate=KokunoPA16CurrentTshCertificate.from_configuration(
                payload.get("certificate")
            ),
            quadrature_points=int(payload.get("quadrature_points")),
            absolute_tolerance=float(payload.get("absolute_tolerance")),
        )

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(_canonical_json(self.configuration()) + "\n", encoding="utf-8")
        return target

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA16CurrentJoinedProfile":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    def report(self, *, eta_probe: float = 0.25) -> dict[str, Any]:
        eta_probe = float(eta_probe)
        return {
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "numerical_realization": dict(_NUMERICAL_REALIZATION),
            "parent_certificate_semantic_sha256": self.certificate.semantic_sha256,
            "geometry": {
                "X_i": self.X_i,
                "log_X_R": self.log_X_R,
                "log_X_h": self.log_X_h,
                "X_h": self.X_h,
                "selected_T_sh": self.certificate.selected_T_sh,
                "route_ready": self.route_ready,
            },
            "Xi_handoff_probe": self.Xi_handoff_report(eta_probe),
            "ideal_exit_probe": self.ideal_exit_report(eta_probe),
            "truth_boundary": self.truth_boundary,
            "semantic_sha256": self.semantic_sha256,
        }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta", type=float, default=0.25)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    profile = KokunoPA16CurrentJoinedProfile()
    text = json.dumps(profile.report(eta_probe=args.eta), indent=2, sort_keys=True, allow_nan=False)
    if args.output is None:
        print(text)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    _main()
