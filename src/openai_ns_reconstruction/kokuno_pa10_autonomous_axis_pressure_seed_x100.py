"""Source-admissible autonomous axis-pressure seed composed with PA.10 reference pressure.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

The corrected reader states that the *source-prepared* Appendix-A axis datum is
holomorphic near ``[-1,1]``, real/even there, and obeys

    Pi_0(eta) <= -(5/2) P_*^2 (1+eta^2)^-2,
    eta Pi_0'(eta) > 0  for eta != 0.

The numerical Appendix-A datum itself is not materialized on this Agent-1
branch.  This module therefore does **not** pretend to recover it.  Instead it
freezes the explicit repository-autonomous seed

    Pi0_seed(eta) = -(5/2) P_*^2 (1+eta^2)^-2,

which saturates the displayed envelope and satisfies the displayed derivative
sign condition.  ``P_*=1`` is the same already-registered autonomous positive
realization used by the Agent-1 ideal-join contract; it is not hidden source
data and is not residual-fitted.  The source's later large-``P_*`` regime is
not claimed.

Composing this autonomous axis seed with the source-determined Agent-1 pressure
increment ``C_p=Pi_r-Pi_0`` gives an executable *autonomous reference pressure*
through ``X=100``.  This is useful as an explicit pressure dependency for the
next ``S_n -> n_{s,r}`` candidate step, but it is not the source-prepared
Appendix-A pressure, not globally matched pressure, and not PDE validation.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa10_ideal_join_profile_contract import DEFAULT_P_STAR
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_pressure_increment_x100 import (
    KokunoPA10ReferencePressureIncrementToX100,
)


SCHEMA = "kokuno-pa10-autonomous-axis-pressure-seed-x100-v1"
FROZEN_P_STAR = float(DEFAULT_P_STAR)

_SOURCE_FORMULAS = {
    "axis_pressure_dependency": (
        "the exact Pi_0(eta) used by the analytic-axis construction is prepared "
        "by the public Appendix-A outer construction"
    ),
    "axis_pressure_envelope": "Pi_0 <= -(5/2) P_*^2 (1+eta^2)^-2",
    "axis_pressure_derivative_sign": "eta Pi_0'(eta)>0 for eta!=0",
    "reference_pressure_composition": "Pi_r=Pi_0+C_p",
    "reference_pressure_increment": "C_p(X,eta)=integral_0^X F_r(s,eta)^2 ds",
    "reference_pressure_radial_derivative": "partial_X Pi_r=partial_X C_p=F_r^2",
}

_AUTONOMOUS_REALIZATION = {
    "axis_pressure_seed": "Pi0_seed=-(5/2) P_*^2 (1+eta^2)^-2",
    "axis_pressure_seed_derivative": "Pi0_seed'=10 P_*^2 eta (1+eta^2)^-3",
    "P_star": FROZEN_P_STAR,
    "selection": "fixed before residual work; not inferred from figures or hidden data",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "displayed_axis_pressure_admissibility_conditions_encoded": True,
    "repository_autonomous_axis_pressure_seed_materialized": True,
    "autonomous_reference_pressure_to_X100_executable": True,
    "source_prepared_appendixA_Pi0_materialized": False,
    "absolute_axis_pressure_Pi0_materialized": False,
    "absolute_reference_pressure_materialized": False,
    "source_large_P_star_regime_admitted": False,
    "outer_tail_pressure_compatibility_verified": False,
    "matched_global_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
    "reference_nsr_materialized": False,
    "complete_reference_stress_pair_materialized": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "kappa0_continuation_to_X100_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
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
class KokunoPA10AutonomousAxisPressureSeedToX100:
    """Autonomous admissible ``Pi0`` seed plus source-determined ``C_p`` to X=100."""

    P_star: float = FROZEN_P_STAR
    pressure_increment: KokunoPA10ReferencePressureIncrementToX100 = field(
        default_factory=KokunoPA10ReferencePressureIncrementToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        p = float(self.P_star)
        if not math.isfinite(p) or p <= 0.0:
            raise ValueError("P_star must be finite and strictly positive")
        if p != FROZEN_P_STAR:
            raise ValueError(
                f"P_star is frozen at the prior Agent-1 autonomous realization {FROZEN_P_STAR:g}; "
                "this pressure seam is not a residual-tuning surface"
            )
        if not isinstance(
            self.pressure_increment, KokunoPA10ReferencePressureIncrementToX100
        ):
            raise TypeError(
                "pressure_increment must be KokunoPA10ReferencePressureIncrementToX100"
            )
        object.__setattr__(self, "P_star", p)

    @property
    def X_b_ref(self) -> float:
        return float(self.pressure_increment.X_b_ref)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.pressure_increment.eta_interval

    def _axis_eta(self, eta: Any) -> np.ndarray:
        eta_arr = _finite(eta, "eta")
        if np.any(np.abs(eta_arr) > math.nextafter(1.0, math.inf)):
            raise ValueError("eta must lie in the displayed analytic-axis interval [-1,1]")
        return np.clip(eta_arr, -1.0, 1.0)

    @staticmethod
    def f(eta: Any) -> np.ndarray:
        eta_arr = _finite(eta, "eta")
        return 1.0 / (1.0 + eta_arr * eta_arr)

    def axis_pressure(self, eta: Any) -> np.ndarray:
        """Return the frozen autonomous axis seed ``Pi0_seed(eta)``."""
        eta_arr = self._axis_eta(eta)
        f = self.f(eta_arr)
        return -2.5 * self.P_star * self.P_star * f * f

    def axis_pressure_eta(self, eta: Any) -> np.ndarray:
        """Return the analytic derivative of the autonomous axis seed."""
        eta_arr = self._axis_eta(eta)
        denom = 1.0 + eta_arr * eta_arr
        return 10.0 * self.P_star * self.P_star * eta_arr / (denom * denom * denom)

    def pressure(self, X: Any, eta: Any) -> np.ndarray:
        """Return autonomous ``Pi_ref_seed=Pi0_seed+C_p_reference`` vectorized."""
        cp = self.pressure_increment.pressure_increment(X, eta)
        return self.axis_pressure(eta) + cp

    def radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return source-exact ``partial_X Pi_ref_seed=F_r^2``."""
        return self.pressure_increment.radial_derivative(X, eta)

    def eta_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return ``Pi0_seed'(eta) + partial_eta C_p``.

        The first term is analytic.  The second inherits the separately-labelled
        fixed FD4 numerical realization from the parent pressure-increment API.
        """
        return self.axis_pressure_eta(eta) + self.pressure_increment.eta_derivative(X, eta)

    def admissibility(self, eta: Any) -> dict[str, np.ndarray]:
        """Expose the two displayed Appendix-B axis admissibility checks."""
        eta_arr = self._axis_eta(eta)
        envelope = -2.5 * self.P_star * self.P_star * self.f(eta_arr) ** 2
        pi0 = self.axis_pressure(eta_arr)
        derivative_product = eta_arr * self.axis_pressure_eta(eta_arr)
        return {
            "eta": eta_arr,
            "Pi0_seed": pi0,
            "source_envelope": envelope,
            "envelope_margin": envelope - pi0,
            "eta_times_Pi0_seed_eta": derivative_product,
        }

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        parent = self.pressure_increment.values(X, eta)
        eta_arr = parent["eta"]
        pi0 = self.axis_pressure(eta_arr)
        return {
            **parent,
            "Pi0_autonomous_seed": pi0,
            "Pi0_autonomous_seed_eta": self.axis_pressure_eta(eta_arr),
            "Pi_reference_autonomous": pi0 + parent["C_p_reference"],
            "Pi_reference_autonomous_X": parent["C_p_reference_X"],
        }

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        parent = self.pressure_increment.handoff_at_X100(eta)
        eta_arr = parent["eta"]
        pi0 = self.axis_pressure(eta_arr)
        return {
            **parent,
            "Pi0_autonomous_seed": pi0,
            "Pi0_autonomous_seed_eta": self.axis_pressure_eta(eta_arr),
            "Pi_reference_autonomous": pi0 + parent["C_p_reference"],
            "Pi_reference_autonomous_X": parent["C_p_reference_X"],
            "Pi_reference_autonomous_eta": (
                self.axis_pressure_eta(eta_arr) + parent["C_p_reference_eta"]
            ),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "P_star": self.P_star,
            "parent_pressure_increment_semantic_sha256": self.pressure_increment.semantic_sha256,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10AutonomousAxisPressureSeedToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "autonomous-axis-pressure configuration is frozen; serialized data differ "
                "from the registered realization"
            )
        return cls()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10AutonomousAxisPressureSeedToX100":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "autonomous_realization": _AUTONOMOUS_REALIZATION,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
        handoff = self.handoff_at_X100(eta_probe)
        admissibility = self.admissibility(eta_probe)
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "autonomous_realization": copy.deepcopy(_AUTONOMOUS_REALIZATION),
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "axis_admissibility_probe": {
                key: value.tolist() for key, value in admissibility.items()
            },
            "X100_handoff_probe": {
                "eta": eta_probe.tolist(),
                "Pi0_autonomous_seed": handoff["Pi0_autonomous_seed"].tolist(),
                "Pi0_autonomous_seed_eta": handoff["Pi0_autonomous_seed_eta"].tolist(),
                "Pi_reference_autonomous": handoff["Pi_reference_autonomous"].tolist(),
                "Pi_reference_autonomous_X": handoff["Pi_reference_autonomous_X"].tolist(),
                "Pi_reference_autonomous_eta": handoff["Pi_reference_autonomous_eta"].tolist(),
            },
            "truth_boundary": self.truth_boundary,
        }

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    pressure = KokunoPA10AutonomousAxisPressureSeedToX100()
    payload = pressure.save_report(args.output)
    if args.config_output is not None:
        pressure.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("X100_handoff_probe=", payload["X100_handoff_probe"])


if __name__ == "__main__":
    _main()
