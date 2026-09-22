"""Materialize the current A1 pulse-entry M/J data required by Kokuno's end solve.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

The public reconstruction defines

    H = sqrt(2X) E,
    M(X,eta) = int_0^X U dx,
    J(X,eta) = int_0^X U H dx,

and the pulse-end 2x2 solve consumes

    m_p = M(X_p)/(X_p E_p),
    j_p = J(X_p)/(X_p H_p E_p).

The current repository lineage is not the hidden source construction.  It uses
an explicitly autonomous finite-N modulation followed by the frozen PA.17 I1
five-moment repair.  That repair is admitted with a nonzero held-out closure
tolerance, so this adapter does *not* set J(X_p)=0 by fiat.  Instead it measures
the actual post-I1 normalized residual in row order (M,J,I,S,C_p) and transports
that residual to X_p.  Since the current leading U is exactly zero from I1 exit
to pulse entry, physical M and J remain constant on that interval.

With PA.17 scales X_1,e_1, E_p=e_b f and f=(1+eta^2)^(-1), the exact conversion
from the post-I1 normalized residuals r_M,r_J is

    m_p = r_M exp[s1(log X_1-log X_p)] / f,
    j_p = r_J exp[s2(log X_1-log X_p)] / (sqrt(2) f^2),

where s1=1/2-lambda and s2=1/2-2 lambda.  Analytic eta derivatives are obtained
from the frozen Chebyshev PA.17 coefficient jet plus the exact f-dependence of
the autonomous discrepancy.  The M result is independently cross-checked
against the already executable #1107 current M/X primitive at X_p.

This is a current-candidate bookkeeping seam.  It does not recover the source
Amp(eta), source-hidden bump shapes, or a paper/OpenAI-exact field; it does not
compose c1/c2 into the Cartesian velocity and is not PDE validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_i1_frozen_gate import (
    I1_CLOSURE_TOLERANCE,
    KokunoPA16CurrentCartesianI1FrozenGate,
)
from .kokuno_pa16_current_cartesian_main_pulse_logx import (
    KokunoPA16CurrentCartesianMainPulseLogX,
)
from .kokuno_public_pulse_end_compensator import (
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_RELEASE,
    SOURCE_RELEASE_DATE,
    SOURCE_REPOSITORY,
    KokunoPublicPulseEndCompensator,
)

SCHEMA = "kokuno-current-pulse-entry-moments-v1"
PARENT_EXACT_HEAD = "c144d00fd81a938e497ca5d841fed6d1fe448a8b"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
MOMENT_ORDER = ("M", "J", "I", "S", "C_p")

_SOURCE_FORMULAS = {
    "moments": (
        "M=int_0^X U dx; H=sqrt(2X)E; J=int_0^X U H dx; "
        "I=int_0^X H dx; S=int_0^X(U^2-E^2/2)dx; C_p=int_0^X E^2/(2x)dx"
    ),
    "pulse_entry_normalization": (
        "m_p=M(X_p)/(X_p E_p); j_p=J(X_p)/(X_p H_p E_p), H_p=sqrt(2X_p)E_p"
    ),
    "row_slopes": "s1=1/2-lambda; s2=1/2-2lambda",
    "PA17_J_scale": "delta J/(X_1^(3/2)e_1^2)",
    "zero_U_transport": (
        "after the completed I1 repair the current leading U=0 through I2/I3/I4 "
        "and the pre-pulse bridge, so cumulative physical M and J are constant"
    ),
    "M_conversion": "m_p=r_M exp[s1(log X_1-log X_p)]/f",
    "J_conversion": "j_p=r_J exp[s2(log X_1-log X_p)]/(sqrt(2)f^2)",
}

_NUMERICAL_REALIZATION = {
    "parent": "stack exactly on A1 #1116 public end-compensator algebra",
    "current_leading": "consume exact #1107 overflow-safe current principal main-pulse identity",
    "I1_history": (
        "read the exact frozen 5e-7 current PA.17 family already embedded in the #1107 ancestry"
    ),
    "J_policy": (
        "use the actual frozen-family post-I1 J closure residual; never replace it by an assumed zero"
    ),
    "autonomous_eta_scaling": (
        "for the selected autonomous modulation A and B/E are eta-independent; its normalized "
        "M/I rows scale as f and J/S/C_p rows as f^2"
    ),
    "M_crosscheck": "independent current primitive path: (M/X)_p/E_p from #1107",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_lineage_M_entry_materialized": True,
    "current_lineage_J_entry_materialized": True,
    "current_J_from_frozen_I1_closure_residual": True,
    "current_J_assumed_zero": False,
    "analytic_eta_jet_materialized": True,
    "current_M_independent_primitive_crosscheck_materialized": True,
    "public_end_compensator_input_contract_ready": True,
    "source_exact_amplitude_root_materialized": False,
    "source_exact_bump_shape_recovered": False,
    "current_cartesian_end_compensation_composed": False,
    "source_hidden_parameters_recovered": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
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


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _eta_array(value: Any) -> np.ndarray:
    eta = _finite(value, "eta")
    if np.any(np.abs(eta) > 1.0):
        raise ValueError("eta must lie in [-1,1]")
    return eta


@dataclass(frozen=True)
class KokunoCurrentPulseEntryMoments:
    """Current-lineage pulse-entry M/J adapter for the public c1/c2 solver."""

    leading: KokunoPA16CurrentCartesianMainPulseLogX = field(
        default_factory=KokunoPA16CurrentCartesianMainPulseLogX,
        repr=False,
        compare=False,
    )
    compensator: KokunoPublicPulseEndCompensator = field(
        default_factory=KokunoPublicPulseEndCompensator,
        repr=False,
        compare=False,
    )
    _i1_gate: KokunoPA16CurrentCartesianI1FrozenGate = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.leading, KokunoPA16CurrentCartesianMainPulseLogX):
            raise TypeError("leading must be KokunoPA16CurrentCartesianMainPulseLogX")
        if not isinstance(self.compensator, KokunoPublicPulseEndCompensator):
            raise TypeError("compensator must be KokunoPublicPulseEndCompensator")
        if self.leading.lambda_value != self.compensator.lambda_value:
            raise ValueError("current leading and public end compensator use different lambda")

        node: Any = self.leading
        gate = None
        for _ in range(10):
            if isinstance(node, KokunoPA16CurrentCartesianI1FrozenGate):
                gate = node
                break
            node = getattr(node, "parent", None)
            if node is None:
                break
        if gate is None:
            raise ValueError("#1107 ancestry no longer contains the frozen current I1 gate")
        if gate.closure_tolerance != I1_CLOSURE_TOLERANCE:
            raise ValueError("current pulse-entry moments require the exact frozen 5e-7 I1 gate")
        if (
            gate.i1_family.modulation.outer_schedule.to_payload()
            != self.leading.outer_schedule.to_payload()
        ):
            raise ValueError("current I1 family and pulse-entry leading use different outer schedules")
        if not gate.i1_family.repair.log_X_1 < self.leading.log_X_p:
            raise ValueError("PA.17 repair scale must precede the pulse entry")
        object.__setattr__(self, "_i1_gate", gate)

    @property
    def i1_gate(self) -> KokunoPA16CurrentCartesianI1FrozenGate:
        return self._i1_gate

    @property
    def i1_family(self):
        return self.i1_gate.i1_family

    @property
    def repair(self):
        return self.i1_family.repair

    @property
    def modulation(self):
        return self.i1_family.modulation

    @property
    def lambda_value(self) -> float:
        return float(self.leading.lambda_value)

    @staticmethod
    def source_f(eta: Any) -> np.ndarray:
        values = _eta_array(eta)
        return 1.0 / (1.0 + values * values)

    @staticmethod
    def source_f_eta(eta: Any) -> np.ndarray:
        values = _eta_array(eta)
        denominator = 1.0 + values * values
        return -2.0 * values / (denominator * denominator)

    def _increment_partial_f(self, coefficients: np.ndarray) -> np.ndarray:
        """Exact partial derivative of PA.17 normalized rows with c held fixed."""

        coeff = np.asarray(coefficients, dtype=float)
        if coeff.shape != (5,):
            raise ValueError("coefficients must have shape (5,)")
        xi, weights = self.repair._quadrature()
        basis = self.repair.basis_values(xi)
        u = np.tensordot(coeff[:2], basis[:2], axes=(0, 0))
        e = np.tensordot(coeff[2:], basis[2:], axes=(0, 0))
        power = np.power(xi, -0.5 - self.lambda_value)
        root = np.sqrt(2.0 * xi)
        return np.asarray(
            [
                0.0,
                float(np.dot(root * u * power, weights)),
                0.0,
                float(np.dot(-power * e, weights)),
                float(np.dot(power * e / xi, weights)),
            ],
            dtype=float,
        )

    def post_i1_residual_normalized_with_eta(
        self, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return actual frozen-family post-I1 residual and analytic eta jet.

        The residual is ``autonomous discrepancy + PA.17 achieved increment``.
        It is *not* forced to zero; the frozen Chebyshev family is admitted with
        the existing 5e-7 held-out construction gate.
        """

        eta_array = _eta_array(eta)
        shape = eta_array.shape
        flat = eta_array.reshape(-1)
        residual = np.empty((flat.size, 5), dtype=float)
        residual_eta = np.empty_like(residual)
        powers = np.asarray([1.0, 2.0, 1.0, 2.0, 2.0], dtype=float)

        for index, eta_value in enumerate(flat):
            f = float(self.source_f(float(eta_value)))
            f_eta = float(self.source_f_eta(float(eta_value)))
            coefficients = np.asarray(
                self.i1_family.coefficients(float(eta_value)), dtype=float
            ).reshape(5)
            coefficients_eta = np.asarray(
                self.i1_family.coefficient_eta(float(eta_value)), dtype=float
            ).reshape(5)

            achieved = np.asarray(
                self.repair.normalized_increment(coefficients, f_eta=f), dtype=float
            )
            discrepancy = np.asarray(
                self.modulation.normalized_discrepancy(float(eta_value)), dtype=float
            )
            jacobian = np.asarray(
                self.repair.coefficient_jacobian(coefficients, f_eta=f), dtype=float
            )
            partial_f = self._increment_partial_f(coefficients)

            achieved_eta = jacobian @ coefficients_eta + partial_f * f_eta
            discrepancy_eta = discrepancy * powers * (f_eta / f)
            residual[index] = achieved + discrepancy
            residual_eta[index] = achieved_eta + discrepancy_eta

        return (
            residual.reshape(shape + (5,)),
            residual_eta.reshape(shape + (5,)),
        )

    def entry_ratios(self, eta: Any) -> dict[str, np.ndarray]:
        """Return current M/J pulse-entry data in the public end-solve normalization."""

        eta_array = _eta_array(eta)
        residual, residual_eta = self.post_i1_residual_normalized_with_eta(eta_array)
        f = self.source_f(eta_array)
        f_eta = self.source_f_eta(eta_array)

        delta_log = float(self.repair.log_X_1 - self.leading.log_X_p)
        scale_m = math.exp(self.compensator.s1 * delta_log)
        scale_j = math.exp(self.compensator.s2 * delta_log) / math.sqrt(2.0)

        r_m = np.asarray(residual[..., 0], dtype=float)
        r_j = np.asarray(residual[..., 1], dtype=float)
        r_m_eta = np.asarray(residual_eta[..., 0], dtype=float)
        r_j_eta = np.asarray(residual_eta[..., 1], dtype=float)

        m_entry = scale_m * r_m / f
        j_entry = scale_j * r_j / (f * f)
        m_entry_eta = scale_m * (r_m_eta / f - r_m * f_eta / (f * f))
        j_entry_eta = scale_j * (
            r_j_eta / (f * f) - 2.0 * r_j * f_eta / (f * f * f)
        )

        E_entry, E_entry_eta, m_over_X, m_eta_over_X = self.leading._pulse_entry_state(
            eta_array
        )
        E_entry = np.asarray(E_entry, dtype=float)
        E_entry_eta = np.asarray(E_entry_eta, dtype=float)
        m_over_X = np.asarray(m_over_X, dtype=float)
        m_eta_over_X = np.asarray(m_eta_over_X, dtype=float)
        m_direct = m_over_X / E_entry
        m_direct_eta = (
            m_eta_over_X * E_entry - m_over_X * E_entry_eta
        ) / (E_entry * E_entry)
        m_difference = m_entry - m_direct
        denominator = np.maximum.reduce(
            [np.abs(m_entry), np.abs(m_direct), np.full_like(m_entry, 1.0e-300)]
        )
        m_relative_difference = np.abs(m_difference) / denominator

        scale_logs = self.repair.log_scale_report()["physical_moment_scale_logs"]
        M_scale = math.exp(float(scale_logs["M"]))
        J_scale = math.exp(float(scale_logs["J"]))

        arrays = (
            m_entry,
            j_entry,
            m_entry_eta,
            j_entry_eta,
            m_direct,
            m_direct_eta,
            m_difference,
            m_relative_difference,
        )
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("current pulse-entry moment conversion produced non-finite values")

        return {
            "eta": np.asarray(eta_array, dtype=float),
            "moment_order": np.asarray(MOMENT_ORDER, dtype=object),
            "post_I1_residual_normalized": np.asarray(residual, dtype=float),
            "post_I1_residual_normalized_eta": np.asarray(residual_eta, dtype=float),
            "M_physical_at_pulse_entry": M_scale * r_m,
            "J_physical_at_pulse_entry": J_scale * r_j,
            "M_physical_at_pulse_entry_eta": M_scale * r_m_eta,
            "J_physical_at_pulse_entry_eta": J_scale * r_j_eta,
            "m_entry_ratio": np.asarray(m_entry, dtype=float),
            "j_entry_ratio": np.asarray(j_entry, dtype=float),
            "m_entry_ratio_eta": np.asarray(m_entry_eta, dtype=float),
            "j_entry_ratio_eta": np.asarray(j_entry_eta, dtype=float),
            "m_entry_ratio_direct_current_primitive": np.asarray(m_direct, dtype=float),
            "m_entry_ratio_direct_current_primitive_eta": np.asarray(
                m_direct_eta, dtype=float
            ),
            "m_entry_crosscheck_difference": np.asarray(m_difference, dtype=float),
            "m_entry_crosscheck_relative_difference": np.asarray(
                m_relative_difference, dtype=float
            ),
        }

    def compensator_inputs(self, eta: Any) -> dict[str, np.ndarray]:
        """Return #1116-ready inputs using the existing autonomous principal Amp proxy.

        This helper prepares low-dimensional inputs only.  The amplitude remains
        repository-autonomous and no c1/c2 correction is composed into velocity.
        """

        entry = self.entry_ratios(eta)
        m_entry = np.asarray(entry["m_entry_ratio"], dtype=float)
        shape = m_entry.shape
        amplitude = np.full(shape, float(self.leading.kernel.A_principal), dtype=float)
        return {
            "amplitude": amplitude,
            "m_entry_ratio": m_entry,
            "j_entry_ratio": np.asarray(entry["j_entry_ratio"], dtype=float),
            "amplitude_eta": np.zeros(shape, dtype=float),
            "m_entry_ratio_eta": np.asarray(entry["m_entry_ratio_eta"], dtype=float),
            "j_entry_ratio_eta": np.asarray(entry["j_entry_ratio_eta"], dtype=float),
            "amplitude_role": np.asarray(
                "repository_autonomous_principal_not_source_exact", dtype=object
            ),
        }

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "blob": SOURCE_BLOB,
                "path": SOURCE_PATH,
                "release": SOURCE_RELEASE,
                "release_date": SOURCE_RELEASE_DATE,
            },
            "current_leading": self.leading.configuration(),
            "public_end_compensator": self.compensator.configuration(),
            "binding": {
                "i1_frozen_semantic_sha256": self.i1_gate.semantic_sha256,
                "i1_closure_tolerance": I1_CLOSURE_TOLERANCE,
                "moment_order": list(MOMENT_ORDER),
                "s1": self.compensator.s1,
                "s2": self.compensator.s2,
                "J_policy": "propagate_actual_frozen_I1_closure_residual_not_assumed_zero",
                "mutable": False,
            },
            "truth_boundary": self.truth_boundary,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoCurrentPulseEntryMoments":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current pulse-entry moments schema")
        leading_payload = payload.get("current_leading")
        compensator_payload = payload.get("public_end_compensator")
        if not isinstance(leading_payload, Mapping) or not isinstance(
            compensator_payload, Mapping
        ):
            raise ValueError("missing current leading or end-compensator configuration")
        candidate = cls(
            leading=KokunoPA16CurrentCartesianMainPulseLogX.from_configuration(
                leading_payload
            ),
            compensator=KokunoPublicPulseEndCompensator.from_configuration(
                compensator_payload
            ),
        )
        if _canonical_json(dict(payload)) != _canonical_json(candidate.configuration()):
            raise ValueError("serialized current pulse-entry source/truth binding changed")
        return candidate

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoCurrentPulseEntryMoments":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def report(self, eta: Any = (-0.83, -0.61, -0.37, 0.19, 0.73)) -> dict[str, Any]:
        eta_array = _eta_array(eta)
        entry = self.entry_ratios(eta_array)
        inputs = self.compensator_inputs(eta_array)
        solution = self.compensator.solve(
            inputs["amplitude"],
            inputs["m_entry_ratio"],
            inputs["j_entry_ratio"],
        )
        closure = np.stack(
            (
                np.asarray(solution.residual_row1, dtype=float),
                np.asarray(solution.residual_row2, dtype=float),
            ),
            axis=-1,
        )
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "eta": eta_array.tolist(),
            "lambda_value": self.lambda_value,
            "log_X_1": float(self.repair.log_X_1),
            "log_X_p": float(self.leading.log_X_p),
            "max_abs_post_I1_normalized_residual": float(
                np.max(np.abs(entry["post_I1_residual_normalized"]))
            ),
            "max_abs_m_entry_ratio": float(np.max(np.abs(entry["m_entry_ratio"]))),
            "max_abs_j_entry_ratio": float(np.max(np.abs(entry["j_entry_ratio"]))),
            "max_m_entry_crosscheck_relative_difference": float(
                np.max(entry["m_entry_crosscheck_relative_difference"])
            ),
            "autonomous_principal_probe_max_abs_c": float(
                np.max(np.abs(np.stack((solution.c1, solution.c2), axis=-1)))
            ),
            "autonomous_principal_probe_max_abs_algebraic_closure_residual": float(
                np.max(np.abs(closure))
            ),
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "J is the current autonomous-candidate value implied by the frozen I1 closure residual, not a recovered source-hidden J",
                "the principal amplitude probe remains repository-autonomous and is not source-exact Amp(eta)",
                "c1/c2 are not composed into the Cartesian velocity in this increment",
                "no pressure, forcing, held-out complete NS residual, or PDE validation is supplied",
            ],
        }
