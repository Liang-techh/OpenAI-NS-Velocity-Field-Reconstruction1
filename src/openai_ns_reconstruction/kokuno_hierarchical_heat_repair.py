"""Hierarchical high-precision realization of the Kokuno three-bump heat repair.

This module closes one *numerical representation* gap between the signed-log
heat-replacement target and the executable I2 reserved-patch correction.

The corrected public reconstruction supplies the three normalized heat-moment
rows and the source requires three compact bumps on an earlier reserved patch.
The repository already implements those rows and one autonomous three-bump
basis in :mod:`kokuno_heat_discrepancy_repair`.  For the default existence-style
outer schedule, however, the corrected target spans more than 400 decades, so a
single float64 coefficient vector erases two nonzero constraints.

Here the existing repository discrete moment map (the same fixed Gauss-Legendre
quadrature and autonomous bump basis) is solved with Python ``decimal`` at
adaptive precision.  The coefficient is kept as a hierarchy of source-channel
linear pieces plus nonlinear Newton counterterms.  This preserves all three
signed-log target channels and also permits a high-precision local I2 velocity
evaluation.

Important scope boundary: this is a high-precision solve of the repository's
*frozen discrete realization* of the source moment map.  It is not a proof that
the continuous source integrals are closed to hundreds of digits, and it does
not create a global Kokuno leading field or pass the Navier--Stokes PDE gate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_heat_discrepancy_repair import KokunoHeatDiscrepancyRepair
from .kokuno_log_rescaled_heat_repair_target import KokunoLogRescaledHeatRepairTarget


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-hierarchical-heat-repair-v1"
CHANNELS = ("C_p", "S", "I_sub")

_SOURCE_FORMULAS = {
    "target": (
        "target=-(Delta C_p/e_*^2, Delta S/(X_*e_*^2), "
        "Delta I_sub/(X_*^(3/2)e_*))"
    ),
    "repair": "delta E=e_* sum_i c_i(eta) beta_i(X/X_*)",
    "rows": (
        "delta C_p/e_*^2=int delta(E^2/e_*^2)/(2x_*) dx_*; "
        "delta S/(X_*e_*^2)=-int delta(E^2/e_*^2)/2 dx_*; "
        "delta I/(X_*^(3/2)e_*)=int sqrt(2x_*) delta E/e_* dx_*"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "corrected_v2_signed_log_target_consumed": True,
    "repository_discrete_three_bump_map_solved_hierarchically": True,
    "actual_heat_discrepancy_applied_to_I2_hierarchical_representation": True,
    "float64_preserves_all_repair_channels": False,
    "continuous_source_moment_compensation_certified": False,
    "heat_compensation_completed": False,
    "core_to_heat_matching_completed": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _D(value: Any) -> Decimal:
    if isinstance(value, Decimal):
        return value
    if isinstance(value, (np.integer, int)):
        return Decimal(int(value))
    value_float = float(value)
    if not math.isfinite(value_float):
        raise ValueError("non-finite value cannot be converted to Decimal")
    return Decimal(repr(value_float))


def _decimal_linear_solve(
    matrix: list[list[Decimal]], rhs: list[Decimal]
) -> list[Decimal]:
    n = len(rhs)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for pivot in range(n):
        row = max(range(pivot, n), key=lambda r: abs(augmented[r][pivot]))
        if augmented[row][pivot] == 0:
            raise ValueError("hierarchical repair matrix is singular")
        if row != pivot:
            augmented[pivot], augmented[row] = augmented[row], augmented[pivot]
        scale = augmented[pivot][pivot]
        augmented[pivot] = [value / scale for value in augmented[pivot]]
        for row_index in range(n):
            if row_index == pivot:
                continue
            factor = augmented[row_index][pivot]
            if factor == 0:
                continue
            augmented[row_index] = [
                augmented[row_index][column]
                - factor * augmented[pivot][column]
                for column in range(n + 1)
            ]
    return [augmented[index][-1] for index in range(n)]


def _decimal_abs_log(value: Decimal) -> float:
    if value == 0:
        return -math.inf
    return float(abs(value).ln())


@dataclass(frozen=True)
class HierarchicalHeatRepairSolution:
    eta: float
    precision_digits: int
    term_labels: tuple[str, ...]
    coefficient_terms: tuple[tuple[str, str, str], ...]
    coefficients: tuple[str, str, str]
    target: tuple[str, str, str]
    achieved: tuple[str, str, str]
    residual: tuple[str, str, str]
    target_sign: tuple[int, int, int]
    target_log_abs: tuple[float, float, float]
    max_relative_residual: float
    newton_steps: int
    sha256: str

    def coefficient_decimals(self) -> tuple[Decimal, Decimal, Decimal]:
        return tuple(Decimal(value) for value in self.coefficients)  # type: ignore[return-value]

    def to_payload(self) -> dict[str, Any]:
        return {
            "eta": self.eta,
            "precision_digits": self.precision_digits,
            "term_labels": list(self.term_labels),
            "coefficient_terms": [list(row) for row in self.coefficient_terms],
            "coefficients": list(self.coefficients),
            "target": list(self.target),
            "achieved": list(self.achieved),
            "residual": list(self.residual),
            "target_sign": list(self.target_sign),
            "target_log_abs": list(self.target_log_abs),
            "max_relative_residual": self.max_relative_residual,
            "newton_steps": self.newton_steps,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class KokunoHierarchicalHeatRepair:
    """Precision-preserving local I2 realization of the source heat repair."""

    target: KokunoLogRescaledHeatRepairTarget = field(
        default_factory=KokunoLogRescaledHeatRepairTarget
    )
    quadrature_points: int = 192
    coefficient_limit: float = 0.05
    guard_digits: int = 64
    relative_residual_digits: int = 40
    max_newton_steps: int = 8

    _repair: KokunoHeatDiscrepancyRepair = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.target, KokunoLogRescaledHeatRepairTarget):
            raise TypeError("target must be a KokunoLogRescaledHeatRepairTarget")
        if self.target.outer_schedule.lambda_outer > 0.5:
            raise ValueError(
                "hierarchical repair currently requires lambda_outer<=0.5, "
                "matching the existing three-bump repair contract"
            )
        for name, value, lo, hi in (
            ("quadrature_points", self.quadrature_points, 64, 256),
            ("guard_digits", self.guard_digits, 32, 160),
            ("relative_residual_digits", self.relative_residual_digits, 20, 80),
            ("max_newton_steps", self.max_newton_steps, 2, 16),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer")
            if not lo <= int(value) <= hi:
                raise ValueError(f"{name} must lie in [{lo},{hi}]")
            object.__setattr__(self, name, int(value))
        limit = float(self.coefficient_limit)
        if not math.isfinite(limit) or not (0.0 < limit <= 0.2):
            raise ValueError("coefficient_limit must lie in (0,0.2]")
        object.__setattr__(self, "coefficient_limit", limit)
        object.__setattr__(
            self,
            "_repair",
            KokunoHeatDiscrepancyRepair(
                lambda_outer=self.target.outer_schedule.lambda_outer,
                quadrature_points=self.quadrature_points,
                coefficient_limit=limit,
            ),
        )

    @property
    def outer_schedule(self):
        return self.target.outer_schedule

    def _discrete_tensors(
        self,
        f_eta: float,
    ) -> tuple[list[list[Decimal]], list[list[list[Decimal]]]]:
        zero = np.zeros(3, dtype=float)
        linear_np = self._repair.coefficient_jacobian(zero, f_eta=f_eta)
        x, weights = self._repair._quadrature()
        basis = self._repair.basis_values(x)
        quadratic_np = np.zeros((3, 3, 3), dtype=float)
        for i in range(3):
            for j in range(3):
                product = basis[i] * basis[j]
                quadratic_np[0, i, j] = np.sum(
                    weights * 0.5 * product / x
                )
                quadratic_np[1, i, j] = np.sum(
                    weights * -0.5 * product
                )
        linear = [[_D(linear_np[r, c]) for c in range(3)] for r in range(3)]
        quadratic = [
            [
                [_D(quadratic_np[r, i, j]) for j in range(3)]
                for i in range(3)
            ]
            for r in range(3)
        ]
        return linear, quadratic

    @staticmethod
    def _map(
        coefficients: list[Decimal],
        linear: list[list[Decimal]],
        quadratic: list[list[list[Decimal]]],
    ) -> list[Decimal]:
        out: list[Decimal] = []
        for row in range(3):
            value = sum(
                (linear[row][i] * coefficients[i] for i in range(3)),
                Decimal(0),
            )
            for i in range(3):
                for j in range(3):
                    value += (
                        quadratic[row][i][j]
                        * coefficients[i]
                        * coefficients[j]
                    )
            out.append(value)
        return out

    @staticmethod
    def _jacobian(
        coefficients: list[Decimal],
        linear: list[list[Decimal]],
        quadratic: list[list[list[Decimal]]],
    ) -> list[list[Decimal]]:
        jacobian = [[linear[r][c] for c in range(3)] for r in range(3)]
        for row in range(3):
            for column in range(3):
                correction = Decimal(0)
                for j in range(3):
                    correction += (
                        quadratic[row][column][j]
                        + quadratic[row][j][column]
                    ) * coefficients[j]
                jacobian[row][column] += correction
        return jacobian

    def _target_decimals(
        self, eta: float, precision_digits: int
    ) -> tuple[list[Decimal], tuple[int, int, int], tuple[float, float, float]]:
        encoded = self.target.signed_log_target(float(eta))
        signs = tuple(int(value) for value in np.asarray(encoded["sign"]).tolist())
        logs = tuple(float(value) for value in np.asarray(encoded["log_abs"]).tolist())
        values: list[Decimal] = []
        with localcontext() as ctx:
            ctx.prec = precision_digits
            for sign, log_abs in zip(signs, logs):
                if sign == 0:
                    values.append(Decimal(0))
                else:
                    values.append(Decimal(sign) * _D(log_abs).exp())
        return values, signs, logs

    def solve(self, eta: float) -> HierarchicalHeatRepairSolution:
        eta_value = float(eta)
        if not math.isfinite(eta_value) or abs(eta_value) > 1.0:
            raise ValueError("eta must lie in [-1,1]")
        report = self.target.precision_report(eta_value)
        precision_digits = max(
            128,
            int(report["estimated_required_decimal_digits"]) + self.guard_digits,
        )
        if precision_digits > 900:
            raise ValueError("required hierarchical precision exceeds the declared bound")

        with localcontext() as ctx:
            ctx.prec = precision_digits
            target, signs, logs = self._target_decimals(
                eta_value, precision_digits
            )
            if not any(signs):
                coefficients = [Decimal(0), Decimal(0), Decimal(0)]
                achieved = [Decimal(0), Decimal(0), Decimal(0)]
                residual = [Decimal(0), Decimal(0), Decimal(0)]
                term_labels: list[str] = []
                terms: list[list[Decimal]] = []
                newton_steps = 0
                max_relative = 0.0
            else:
                f_eta = float(self.outer_schedule.source_f(eta_value))
                linear, quadratic = self._discrete_tensors(f_eta)
                coefficients = [Decimal(0), Decimal(0), Decimal(0)]
                term_labels = []
                terms = []

                for channel_index, channel_name in enumerate(CHANNELS):
                    if target[channel_index] == 0:
                        continue
                    rhs = [Decimal(0), Decimal(0), Decimal(0)]
                    rhs[channel_index] = target[channel_index]
                    term = _decimal_linear_solve(linear, rhs)
                    terms.append(term)
                    term_labels.append(f"linear_{channel_name}")
                    coefficients = [
                        coefficients[i] + term[i] for i in range(3)
                    ]

                tolerance = Decimal(10) ** (-self.relative_residual_digits)
                newton_steps = 0
                for step in range(1, self.max_newton_steps + 1):
                    achieved = self._map(coefficients, linear, quadratic)
                    residual = [
                        achieved[row] - target[row] for row in range(3)
                    ]
                    relative = [
                        abs(residual[row]) / abs(target[row])
                        for row in range(3)
                        if target[row] != 0
                    ]
                    if relative and max(relative) <= tolerance:
                        break
                    jacobian = self._jacobian(
                        coefficients, linear, quadratic
                    )
                    delta = _decimal_linear_solve(
                        jacobian, [-value for value in residual]
                    )
                    coefficients = [
                        coefficients[i] + delta[i] for i in range(3)
                    ]
                    terms.append(delta)
                    term_labels.append(f"newton_{step}")
                    newton_steps = step
                else:
                    raise ValueError(
                        "hierarchical three-bump solve did not reach the declared "
                        "relative residual"
                    )

                achieved = self._map(coefficients, linear, quadratic)
                residual = [
                    achieved[row] - target[row] for row in range(3)
                ]
                relative = [
                    abs(residual[row]) / abs(target[row])
                    for row in range(3)
                    if target[row] != 0
                ]
                max_relative = float(max(relative)) if relative else 0.0

            limit = _D(self.coefficient_limit)
            if any(abs(value) > limit for value in coefficients):
                raise ValueError(
                    "hierarchical coefficient exceeds the declared local repair bound"
                )

            coefficient_strings = tuple(str(value) for value in coefficients)
            term_strings = tuple(
                tuple(str(value) for value in term) for term in terms
            )
            target_strings = tuple(str(value) for value in target)
            achieved_strings = tuple(str(value) for value in achieved)
            residual_strings = tuple(str(value) for value in residual)
            unsigned = {
                "schema": "kokuno-hierarchical-heat-repair-solution-v1",
                "eta": eta_value,
                "precision_digits": precision_digits,
                "term_labels": term_labels,
                "coefficient_terms": [list(row) for row in term_strings],
                "coefficients": list(coefficient_strings),
                "target": list(target_strings),
                "achieved": list(achieved_strings),
                "residual": list(residual_strings),
                "target_sign": list(signs),
                "target_log_abs": list(logs),
                "max_relative_residual": max_relative,
                "newton_steps": newton_steps,
            }
            digest = hashlib.sha256(
                _canonical_json(unsigned).encode("utf-8")
            ).hexdigest()
            return HierarchicalHeatRepairSolution(
                eta=eta_value,
                precision_digits=precision_digits,
                term_labels=tuple(term_labels),
                coefficient_terms=term_strings,
                coefficients=coefficient_strings,  # type: ignore[arg-type]
                target=target_strings,  # type: ignore[arg-type]
                achieved=achieved_strings,  # type: ignore[arg-type]
                residual=residual_strings,  # type: ignore[arg-type]
                target_sign=signs,
                target_log_abs=logs,
                max_relative_residual=max_relative,
                newton_steps=newton_steps,
                sha256=digest,
            )

    @staticmethod
    def _decimal_bump(x_star: Decimal, center: float) -> Decimal:
        half_width = Decimal("0.105")
        s = (x_star - _D(center)) / half_width
        if abs(s) >= 1:
            return Decimal(0)
        return (Decimal(1) - Decimal(1) / (Decimal(1) - s * s)).exp()

    def correction_over_e_star_decimal(
        self, x_star: float, eta: float
    ) -> Decimal:
        x_value = float(x_star)
        if not math.isfinite(x_value) or x_value <= 0.0:
            raise ValueError("x_star must be positive and finite")
        solution = self.solve(float(eta))
        with localcontext() as ctx:
            ctx.prec = solution.precision_digits
            x_decimal = _D(x_value)
            coefficients = list(solution.coefficient_decimals())
            centers = (0.96, 1.24, 1.53)
            return sum(
                (
                    coefficients[index]
                    * self._decimal_bump(x_decimal, centers[index])
                    for index in range(3)
                ),
                Decimal(0),
            )

    def velocity_decimal(
        self, x: float, y: float, z: float, t: float
    ) -> tuple[Decimal, Decimal, Decimal]:
        """Return the local I2 corrected velocity without collapsing its hierarchy.

        Similarity coordinates are taken from the already tested repository
        float evaluator.  The profile correction and the final base+correction
        addition are then carried in Decimal at the precision required by the
        signed-log target.
        """

        coordinates = self.outer_schedule._coordinates.evaluate(x, y, z, t)
        for name in ("X", "eta", "q"):
            if np.asarray(coordinates[name]).shape != ():
                raise ValueError("velocity_decimal currently requires scalar inputs")
        X = float(coordinates["X"])
        eta = float(coordinates["eta"])
        q = float(coordinates["q"])
        profile = self.outer_schedule.patch_profile_values(X, eta)
        x_star = float(profile["x_star"])
        base_velocity = np.asarray(
            self.outer_schedule.patch_velocity(x, y, z, t), dtype=float
        )
        solution = self.solve(eta)
        with localcontext() as ctx:
            ctx.prec = solution.precision_digits
            correction_ratio = self.correction_over_e_star_decimal(
                x_star, eta
            )
            e_star = _D(self.outer_schedule.log_e_star).exp()
            delta_E = e_star * correction_ratio
            exponent = -Decimal(1) - _D(self.outer_schedule.h)
            q_factor = (exponent * _D(q).ln()).exp()
            geom = q_factor / (Decimal(2) * _D(X)).sqrt()
            delta = (
                -_D(y) * geom * delta_E,
                _D(x) * geom * delta_E,
                Decimal(0),
            )
            return tuple(
                _D(base_velocity[index]) + delta[index]
                for index in range(3)
            )  # type: ignore[return-value]

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Float64 compatibility view of the hierarchical local I2 velocity.

        The Decimal path is authoritative for the tiny repair.  This method is
        supplied for candidate-evaluator compatibility and may round the repair
        away when it is below one ulp of the base velocity.
        """

        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        output = np.empty(x_array.shape + (3,), dtype=float)
        for index in np.ndindex(x_array.shape):
            values = self.velocity_decimal(
                float(x_array[index]),
                float(y_array[index]),
                float(z_array[index]),
                float(t_array[index]),
            )
            output[index] = [float(value) for value in values]
        return output

    def precision_report(self, eta: float) -> dict[str, Any]:
        target_report = self.target.precision_report(float(eta))
        solution = self.solve(float(eta))
        coefficient_logs = [
            _decimal_abs_log(value)
            for value in solution.coefficient_decimals()
        ]
        return {
            "eta": float(eta),
            "target": target_report,
            "solve_precision_digits": solution.precision_digits,
            "solution_sha256": solution.sha256,
            "newton_steps": solution.newton_steps,
            "max_relative_residual": solution.max_relative_residual,
            "coefficient_log10_abs": [
                value / math.log(10.0) if math.isfinite(value) else -math.inf
                for value in coefficient_logs
            ],
            "float64_velocity_preserves_all_repair_channels": False,
            "continuous_moment_compensation_certified": False,
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "provenance": {
                "source_repository": SOURCE_REPOSITORY,
                "source_commit": SOURCE_COMMIT,
                "source_path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "target": self.target.to_payload(),
            "parameters": {
                "quadrature_points": self.quadrature_points,
                "coefficient_limit": self.coefficient_limit,
                "guard_digits": self.guard_digits,
                "relative_residual_digits": self.relative_residual_digits,
                "max_newton_steps": self.max_newton_steps,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "autonomous_numerics": {
                "high_precision_engine": "python decimal",
                "discrete_map": (
                    "existing autonomous three-bump basis and fixed "
                    "Gauss-Legendre quadrature"
                ),
                "hierarchy": (
                    "per-target-channel linear pieces followed by exact "
                    "discrete-polynomial Newton counterterms"
                ),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoHierarchicalHeatRepair":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected hierarchical heat-repair schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source formula contract mismatch")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("truth boundary mismatch")
        target_payload = payload.get("target")
        params = payload.get("parameters")
        if not isinstance(target_payload, dict) or not isinstance(params, dict):
            raise ValueError("serialized target/parameters are missing")
        obj = cls(
            target=KokunoLogRescaledHeatRepairTarget.from_payload(target_payload),
            **params,
        )
        if obj.to_payload() != payload:
            raise ValueError("hierarchical heat-repair payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoHierarchicalHeatRepair":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
