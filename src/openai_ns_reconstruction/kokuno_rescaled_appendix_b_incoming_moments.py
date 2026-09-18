"""PA.15 incoming five-moment target on the selected source-rescaled path.

This module is the scale-correct successor of the earlier finite-``C``
``KokunoAppendixBIncomingMoments`` bridge.  KokunoYumeto's corrected
2026-09-09 reconstruction retains the five prefix moments

    M   = int U dX,
    I   = int H dX,             H = sqrt(2X) E,
    J   = int U H dX,
    S   = int (U^2-E^2/2) dX,
    C_p = int E^2/(2X) dX,

and, for ``x=X/X_R``, PA.15 rescales them as

    (M,I,J,S,C_p)
      = (X_R Mhat, X_R^(3/2) Ihat, X_R^(3/2) Jhat,
         X_R Shat, Cphat).

The current selected-pressure Appendix-B trajectory is executable through
``X_i=110`` in :mod:`kokuno_rescaled_appendix_b_boundary`.  PR #481 also made
an important source dependency explicit: the same fixed ``C`` used in
``F=phi/C`` must be used in ``X_R=110(C P_*)^10``.  On the selected real-axis
normalization this makes ``log X_R`` far outside a range where ``X_R`` or
``x_i=X_i/X_R`` can be materialized in binary64.

The implementation below therefore:

* propagates the selected Appendix-B five physical prefix moments using the
  same source ODE as the executable boundary;
* applies PA.15 row scaling with signed log arithmetic, without materializing
  ``X_R``;
* evaluates the ideal-pair prefix entirely from ``log x_i``;
* exposes endpoint ``ell_i,G_i`` and a vectorized incoming discrepancy in the
  PA.16 row order ``(M,I,J,S,C_p)``;
* records when a mathematically nonzero normalized row is below binary64
  materialization range rather than silently treating scale underflow as a
  source identity.

This remains one explicitly selected/autonomous realization.  The upstream
source fixed point and complex-domain ``C`` certificate are not reconstructed,
``T_sh`` is not certified here, and the result is not a completed PA.16 join,
PDE validation, paper-exact field, or OpenAI-field identification.
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
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .kokuno_rescaled_outer_scale_binding import KokunoRescaledOuterScaleBinding
from .kokuno_rescaled_appendix_b_boundary import X_I


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-rescaled-appendix-b-incoming-moments-v1"
ROW_NAMES = ("M", "I", "J", "S", "C_p")

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
    "incoming_discrepancy": (
        "selected source-rescaled Appendix-B five-prefix moments at X_i=110 "
        "minus ideal-pair moments at x_i=X_i/X_R"
    ),
    "shared_C": "the same fixed C in F=phi/C is used in X_R=110(C P_*)^10",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_source_rescaled_appendix_B_prefix_moments_executable": True,
    "shared_C_PA15_scaling_applied": True,
    "rescaled_PA15_incoming_discrepancy_executable": True,
    "PA16_input_tuple_executable": True,
    "ideal_prefix_evaluated_in_log_space": True,
    "binary64_scale_underflow_reported_explicitly": True,
    "selected_real_axis_C_is_autonomous_normalization": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "source_complex_C_bound_verified": False,
    "source_all_PA11_C_bounds_verified": False,
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

_LOG_MAX = math.log(np.finfo(float).max)
_LOG_TINY = math.log(np.finfo(float).tiny)


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite_array(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _exp_array_from_log(log_values: Any, name: str) -> np.ndarray:
    values = np.asarray(log_values, dtype=float)
    if np.any(~np.isfinite(values)):
        raise OverflowError(f"{name} logarithm became non-finite")
    if np.any(values > _LOG_MAX):
        raise OverflowError(f"{name} exceeds binary64 range")
    with np.errstate(under="ignore"):
        return np.exp(values)


def _positive_from_log(log_value: float, name: str) -> tuple[float, bool]:
    value = float(log_value)
    if not math.isfinite(value):
        raise OverflowError(f"{name} logarithm became non-finite")
    if value > _LOG_MAX:
        raise OverflowError(f"{name} exceeds binary64 range")
    if value < _LOG_TINY:
        return 0.0, True
    return math.exp(value), False


def _signed_scale(value: float, log_factor: float, name: str) -> tuple[float, bool]:
    raw = float(value)
    factor = float(log_factor)
    if not math.isfinite(raw) or not math.isfinite(factor):
        raise OverflowError(f"{name} scale input became non-finite")
    if raw == 0.0:
        return 0.0, False
    log_abs = math.log(abs(raw)) + factor
    scaled_abs, underflow = _positive_from_log(log_abs, name)
    return math.copysign(scaled_abs, raw), underflow


@dataclass(frozen=True)
class RescaledAppendixBIncomingMomentResult:
    """One scalar-eta selected receipt in PA.15 row order."""

    eta: float
    ell_i: float
    G_i: float
    actual_physical_moments: tuple[float, float, float, float, float]
    actual_scaled_moments: tuple[float, float, float, float, float]
    ideal_scaled_moments: tuple[float, float, float, float, float]
    incoming_scaled_discrepancy: tuple[float, float, float, float, float]
    actual_scaled_underflow_rows: tuple[bool, bool, bool, bool, bool]
    ideal_scaled_underflow_rows: tuple[bool, bool, bool, bool, bool]
    max_abs_incoming_discrepancy: float


@dataclass(frozen=True)
class KokunoRescaledAppendixBIncomingMoments:
    """Compute the selected source-rescaled PA.15 incoming moment target."""

    binding: KokunoRescaledOuterScaleBinding = field(
        default_factory=KokunoRescaledOuterScaleBinding
    )
    axis_quadrature_points: int = 64
    _axis_nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _axis_weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.binding, KokunoRescaledOuterScaleBinding):
            raise TypeError("binding must be a KokunoRescaledOuterScaleBinding")
        if isinstance(self.axis_quadrature_points, bool) or not isinstance(
            self.axis_quadrature_points, (int, np.integer)
        ):
            raise TypeError("axis_quadrature_points must be an integer")
        order = int(self.axis_quadrature_points)
        if not 32 <= order <= 192:
            raise ValueError("axis_quadrature_points must lie in [32,192]")
        if self.binding.finite_template_C_matches_selected_C:
            raise ValueError(
                "rescaled incoming moments require the selected shared-C scale, "
                "not the old finite-C template"
            )
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "axis_quadrature_points", order)
        object.__setattr__(self, "_axis_nodes", nodes)
        object.__setattr__(self, "_axis_weights", weights)

    @property
    def boundary(self):
        return self.binding.boundary

    @property
    def outer_schedule(self):
        return self.binding.template_outer_schedule

    @property
    def log_x_i(self) -> float:
        return float(self.binding.log_x_i)

    @property
    def log_P_star(self) -> float:
        return float(self.outer_schedule.log_P_star)

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X_i": X_I,
            "selected_log_C": self.binding.selected_log_C,
            "selected_log_X_R": self.binding.log_X_R,
            "selected_log_x_i": self.log_x_i,
            "log_P_star": self.log_P_star,
            "axis_quadrature_points": self.axis_quadrature_points,
            "finite_template_C_matches_selected_C": self.binding.finite_template_C_matches_selected_C,
            "source_T_sh_lower_bound_verified": False,
        }

    def _axis_extra_physical_moments(self, eta: float) -> np.ndarray:
        """Integrate ``(I,J,S,C_p)`` over the smooth ``0<=X<=X0`` prefix."""

        X0 = float(self.boundary.X0)
        X = 0.5 * X0 * (self._axis_nodes + 1.0)
        weights = 0.5 * X0 * self._axis_weights
        values = self.boundary.reference.profile_values(X, eta)
        log_F = np.asarray(values["log_F"], dtype=float)
        U = np.asarray(values["U"], dtype=float)
        if np.any(~np.isfinite(log_F)) or np.any(~np.isfinite(U)):
            raise RuntimeError("rescaled reference axis profile produced invalid values")
        F = _exp_array_from_log(log_F, "axis F")
        F_sq = _exp_array_from_log(2.0 * log_F, "axis F^2")
        H = 2.0 * X * F
        densities = np.stack(
            (
                H,
                U * H,
                U * U - X * F_sq,
                F_sq,
            ),
            axis=0,
        )
        out = np.asarray(densities @ weights, dtype=float)
        if np.any(~np.isfinite(out)):
            raise RuntimeError("rescaled axis moment quadrature produced invalid values")
        return out

    def _augmented_rhs(
        self, y: float, state: np.ndarray, eta: float, reference_slice: Any
    ) -> np.ndarray:
        base = np.asarray(
            self.boundary._rhs(float(y), np.asarray(state[:6], dtype=float), eta, reference_slice),
            dtype=float,
        )
        X = float(self.boundary.X0) * math.exp(float(y))
        log_F = float(state[2])
        F, _ = _positive_from_log(log_F, "Appendix-B F")
        F_sq, _ = _positive_from_log(2.0 * log_F, "Appendix-B F^2")
        U = float(state[3])
        H = 2.0 * X * F
        # d/dy = X d/dX.
        extra = np.asarray(
            (
                X * H,
                X * U * H,
                X * (U * U - X * F_sq),
                X * F_sq,
            ),
            dtype=float,
        )
        out = np.concatenate((base, extra))
        if np.any(~np.isfinite(out)):
            raise OverflowError("rescaled Appendix-B moment derivative became non-finite")
        return out

    def _selected_actual_physical_moments(
        self, eta: float
    ) -> tuple[np.ndarray, float, float]:
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        initial_base = np.asarray(self.boundary._initial_state(eta), dtype=float)
        initial_extra = self._axis_extra_physical_moments(eta)
        initial = np.concatenate((initial_base, initial_extra))
        reference_slice = self.boundary._make_slice(eta)
        solved = solve_ivp(
            lambda yy, zz: self._augmented_rhs(yy, zz, eta, reference_slice),
            (0.0, float(self.boundary.y_i)),
            initial,
            method="DOP853",
            rtol=float(self.boundary.rtol),
            atol=float(self.boundary.atol),
            max_step=float(self.boundary.max_step),
        )
        if not solved.success or solved.y.size == 0:
            raise RuntimeError(
                f"source-rescaled Appendix-B moment propagation failed: {solved.message}"
            )
        final = np.asarray(solved.y[:, -1], dtype=float)
        if final.shape != (10,) or np.any(~np.isfinite(final)):
            raise RuntimeError("source-rescaled moment propagation produced invalid state")
        log_E_i = 0.5 * math.log(2.0 * X_I) + float(final[2])
        ell_i = self.binding.selected_log_C + log_E_i
        G_i = float(final[3])
        physical = np.asarray(
            (final[4], final[6], final[7], final[8], final[9]), dtype=float
        )
        return physical, float(ell_i), G_i

    def _actual_scaled_moments(
        self, physical: np.ndarray
    ) -> tuple[np.ndarray, tuple[bool, bool, bool, bool, bool]]:
        logs = self.binding.pa15_log_inverse_scales()
        values: list[float] = []
        underflow: list[bool] = []
        for row, raw in zip(ROW_NAMES, np.asarray(physical, dtype=float), strict=True):
            scaled, hit = _signed_scale(float(raw), float(logs[row]), f"scaled {row}")
            values.append(scaled)
            underflow.append(bool(hit))
        return np.asarray(values, dtype=float), tuple(underflow)  # type: ignore[return-value]

    def ideal_scaled_moments(
        self, eta: float
    ) -> tuple[np.ndarray, tuple[bool, bool, bool, bool, bool]]:
        """Evaluate the ideal-pair PA.15 prefix on ``0<=x<=x_i`` in log space."""

        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        log_f = -math.log1p(eta * eta)
        lx = self.log_x_i
        lp = self.log_P_star

        x, ux = _positive_from_log(lx, "ideal x_i")
        i_value, ui = _positive_from_log(
            math.log(math.sqrt(2.0) / 1.6) + lp + log_f + 1.6 * lx,
            "ideal I",
        )
        s_energy, us_energy = _positive_from_log(
            math.log(0.5 / 1.2) + 2.0 * lp + 2.0 * log_f + 1.2 * lx,
            "ideal S energy term",
        )
        cp, ucp = _positive_from_log(
            math.log(2.5) + 2.0 * lp + 2.0 * log_f + 0.2 * lx,
            "ideal C_p",
        )
        M = 4.0 * eta * x
        J = 4.0 * eta * i_value
        S = 16.0 * eta * eta * x - s_energy
        values = np.asarray((M, i_value, J, S, cp), dtype=float)
        if np.any(~np.isfinite(values)):
            raise RuntimeError("log-space ideal PA.15 prefix produced non-finite values")
        underflow = (
            bool(ux and eta != 0.0),
            bool(ui),
            bool(ui and eta != 0.0),
            bool((ux and eta != 0.0) or us_energy),
            bool(ucp),
        )
        return values, underflow

    def result_at_eta(self, eta: float) -> RescaledAppendixBIncomingMomentResult:
        physical, ell_i, G_i = self._selected_actual_physical_moments(float(eta))
        actual, actual_underflow = self._actual_scaled_moments(physical)
        ideal, ideal_underflow = self.ideal_scaled_moments(float(eta))
        incoming = actual - ideal
        if np.any(~np.isfinite(incoming)):
            raise RuntimeError("rescaled incoming PA.15 discrepancy is non-finite")
        return RescaledAppendixBIncomingMomentResult(
            eta=float(eta),
            ell_i=ell_i,
            G_i=G_i,
            actual_physical_moments=tuple(map(float, physical)),
            actual_scaled_moments=tuple(map(float, actual)),
            ideal_scaled_moments=tuple(map(float, ideal)),
            incoming_scaled_discrepancy=tuple(map(float, incoming)),
            actual_scaled_underflow_rows=actual_underflow,
            ideal_scaled_underflow_rows=ideal_underflow,
            max_abs_incoming_discrepancy=float(np.max(np.abs(incoming))),
        )

    def incoming_scaled_discrepancy(self, eta: Any) -> np.ndarray:
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        out = np.empty(eta_array.shape + (5,), dtype=float)
        flat = out.reshape(-1, 5)
        for index, value in enumerate(eta_array.reshape(-1)):
            flat[index, :] = np.asarray(
                self.result_at_eta(float(value)).incoming_scaled_discrepancy,
                dtype=float,
            )
        return out

    def pa16_inputs_at_eta(self, eta: float) -> dict[str, Any]:
        receipt = self.result_at_eta(float(eta))
        return {
            "eta": receipt.eta,
            "ell_i": receipt.ell_i,
            "G_i": receipt.G_i,
            "incoming_scaled_discrepancy": list(receipt.incoming_scaled_discrepancy),
            "row_order": list(ROW_NAMES),
            "source_T_sh_lower_bound_verified": False,
        }

    def report(self) -> dict[str, Any]:
        stationary = float(self.boundary.reference.seed.phase_stationary_eta)
        sample_etas = (stationary, 0.0)
        receipts = [self.result_at_eta(value) for value in sample_etas]
        return {
            "schema": SCHEMA,
            "binding_sha256": self.binding.sha256,
            "geometry": self.geometry_report(),
            "sample_receipts": [
                {
                    "eta": item.eta,
                    "ell_i": item.ell_i,
                    "G_i": item.G_i,
                    "actual_physical_moments": list(item.actual_physical_moments),
                    "actual_scaled_moments": list(item.actual_scaled_moments),
                    "ideal_scaled_moments": list(item.ideal_scaled_moments),
                    "incoming_scaled_discrepancy": list(item.incoming_scaled_discrepancy),
                    "actual_scaled_underflow_rows": list(item.actual_scaled_underflow_rows),
                    "ideal_scaled_underflow_rows": list(item.ideal_scaled_underflow_rows),
                    "max_abs_incoming_discrepancy": item.max_abs_incoming_discrepancy,
                }
                for item in receipts
            ],
            "truth_boundary": dict(_TRUTH_BOUNDARY),
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
            "binding": self.binding.to_payload(),
            "parameters": {"axis_quadrature_points": self.axis_quadrature_points},
            "autonomous_numerics": {
                "axis_quadrature": "Gauss-Legendre on 0<=X<=X0",
                "trajectory": "same DOP853 source-y ODE/tolerances as selected Appendix-B boundary",
                "PA15_scaling": "signed log-space scaling; materialization underflow reported",
                "ideal_prefix": "analytic source formula evaluated from log(x_i)",
            },
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoRescaledAppendixBIncomingMoments":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected rescaled Appendix-B incoming-moment schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("rescaled Appendix-B incoming source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("rescaled Appendix-B incoming truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("rescaled Appendix-B incoming payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing rescaled Appendix-B incoming parameters")
        obj = cls(
            binding=KokunoRescaledOuterScaleBinding.from_payload(payload.get("binding")),
            axis_quadrature_points=int(params["axis_quadrature_points"]),
        )
        if obj.to_payload() != payload:
            raise ValueError("rescaled Appendix-B incoming replay changed payload")
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
    def load_json(cls, path: str | Path) -> "KokunoRescaledAppendixBIncomingMoments":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Emit selected source-rescaled PA.15 incoming five-moment receipt"
    )
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--axis-quadrature-points", type=int, default=64)
    args = parser.parse_args(argv)
    moments = KokunoRescaledAppendixBIncomingMoments(
        axis_quadrature_points=args.axis_quadrature_points
    )
    text = json.dumps(moments.report(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
