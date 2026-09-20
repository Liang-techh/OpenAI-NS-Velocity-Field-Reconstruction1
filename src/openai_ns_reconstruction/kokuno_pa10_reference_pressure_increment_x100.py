"""Executable public PA.10 reference pressure *increment* through ``X=100``.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

The public reconstruction gives, with ``F=phi/C`` and ``E=sqrt(2X)F``,

    Pi_X = E^2/(2X) = F^2,
    Pi_r(X,eta) - Pi_0(eta) = integral_0^X F_r(x,eta)^2 dx.

The outer construction prepares the axis-pressure datum ``Pi_0(eta)``.  That
absolute datum is not yet materialized on the Agent-1 continuation branch, so
this module deliberately exposes only the source-determined pressure increment
``C_p = Pi_r-Pi_0``.  It must not be relabelled as matched/global pressure.

Production evaluation uses a fixed Gauss--Legendre quadrature for the radial
integral.  The radial derivative is the exact source identity ``C_{p,X}=F^2``.
A fixed fourth-order eta finite difference is exposed separately for the next
``S_n``/``n_s`` stage and is explicitly a repository numerical realization,
not a hidden source formula or fitted parameter.
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

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_plateau_x100 import KokunoPA10ReferencePlateauToX100


SCHEMA = "kokuno-pa10-reference-pressure-increment-x100-v1"
RADIAL_QUADRATURE_ORDER = 96
ETA_FD_STEP = 2.0e-5

_SOURCE_FORMULAS = {
    "pressure_radial_identity": "Pi_X=E^2/(2X)=F^2",
    "reference_pressure_increment": (
        "Pi_r(X,eta)-Pi_0(eta)=C^{-2} integral_0^X phi_r(x,eta)^2 dx "
        "= integral_0^X F_r(x,eta)^2 dx"
    ),
    "pressure_datum_dependency": (
        "Pi_0(eta) is prepared by the public outer construction; it is not "
        "identified by the local radial identity alone"
    ),
}

_NUMERICAL_REALIZATION = {
    "radial_integral": f"fixed Gauss-Legendre order {RADIAL_QUADRATURE_ORDER}",
    "eta_derivative": f"fixed centered FD4 with h={ETA_FD_STEP:.17g}",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "reference_pressure_increment_to_X100_executable": True,
    "reference_pressure_increment_radial_derivative_source_exact": True,
    "reference_pressure_increment_eta_derivative_numerical_realization": True,
    "absolute_axis_pressure_Pi0_materialized": False,
    "absolute_reference_pressure_materialized": False,
    "reference_nsr_materialized": False,
    "complete_reference_stress_pair_materialized": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "kappa0_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_inner_to_outer_bridge_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
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


def _gl_nodes_weights(order: int = RADIAL_QUADRATURE_ORDER) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return np.asarray(nodes, dtype=float), np.asarray(weights, dtype=float)


_GL_NODES, _GL_WEIGHTS = _gl_nodes_weights()


@dataclass(frozen=True)
class KokunoPA10ReferencePressureIncrementToX100:
    """Source-determined ``C_p=Pi_r-Pi_0`` on the reference profile to X=100."""

    reference: KokunoPA10ReferencePlateauToX100 = field(
        default_factory=KokunoPA10ReferencePlateauToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoPA10ReferencePlateauToX100):
            raise TypeError("reference must be KokunoPA10ReferencePlateauToX100")

    @property
    def X_b_ref(self) -> float:
        return float(self.reference.X_b_ref)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.reference.eta_interval

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_hi = math.nextafter(self.X_b_ref, math.inf)
        if np.any(X_arr < 0.0) or np.any(X_arr > x_hi):
            raise ValueError(f"X must lie in the public reference domain [0,{self.X_b_ref}]")
        eta_lo, eta_hi = self.eta_interval
        if np.any((eta_arr < eta_lo) | (eta_arr > eta_hi)):
            raise ValueError(
                f"eta must lie in the executable reference interval [{eta_lo},{eta_hi}]"
            )
        return np.minimum(X_arr, self.X_b_ref), eta_arr

    def pressure_increment(self, X: Any, eta: Any) -> np.ndarray:
        """Return ``C_p(X,eta)=integral_0^X F_r(s,eta)^2 ds`` vectorized."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        xs = X_arr.reshape(-1)
        etas = eta_arr.reshape(-1)
        out = np.zeros_like(xs)

        positive = xs > 0.0
        if np.any(positive):
            xp = xs[positive]
            ep = etas[positive]
            # One deterministic tensor batch: node axis first, sample axis second.
            nodes = 0.5 * xp[None, :] * (_GL_NODES[:, None] + 1.0)
            weights = 0.5 * xp[None, :] * _GL_WEIGHTS[:, None]
            eta_grid = np.broadcast_to(ep[None, :], nodes.shape)
            F = self.reference.values(nodes, eta_grid)["F_reference"]
            out[positive] = np.sum(weights * F * F, axis=0)
        return out.reshape(shape)

    def radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return the exact source derivative ``partial_X C_p = F_r^2``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        F = self.reference.values(X_arr, eta_arr)["F_reference"]
        return F * F

    def eta_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return a frozen FD4 realization of ``partial_eta C_p``.

        This is intentionally separate from the source-exact radial derivative.
        It fails closed when the centered stencil would leave the executable
        eta interval.
        """
        X_arr, eta_arr = self._broadcast(X, eta)
        h = ETA_FD_STEP
        eta_lo, eta_hi = self.eta_interval
        if np.any(eta_arr - 2.0 * h < eta_lo) or np.any(eta_arr + 2.0 * h > eta_hi):
            raise ValueError("eta FD4 stencil leaves the executable reference interval")
        fm2 = self.pressure_increment(X_arr, eta_arr - 2.0 * h)
        fm1 = self.pressure_increment(X_arr, eta_arr - h)
        fp1 = self.pressure_increment(X_arr, eta_arr + h)
        fp2 = self.pressure_increment(X_arr, eta_arr + 2.0 * h)
        return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) / (12.0 * h)

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast(X, eta)
        return {
            "X": X_arr,
            "eta": eta_arr,
            "C_p_reference": self.pressure_increment(X_arr, eta_arr),
            "C_p_reference_X": self.radial_derivative(X_arr, eta_arr),
        }

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_b_ref, dtype=float)
        return {
            **self.values(X, eta_arr),
            "C_p_reference_eta": self.eta_derivative(X, eta_arr),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "radial_quadrature_order": RADIAL_QUADRATURE_ORDER,
            "eta_fd_step": ETA_FD_STEP,
            "parent_reference_semantic_sha256": self.reference.semantic_sha256,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10ReferencePressureIncrementToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "reference-pressure-increment configuration is frozen; serialized data "
                "differ from the registered source/numerical realization"
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
    ) -> "KokunoPA10ReferencePressureIncrementToX100":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
        handoff = self.handoff_at_X100(eta_probe)
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
            "numerical_realization": copy.deepcopy(_NUMERICAL_REALIZATION),
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "X100_handoff_probe": {
                "eta": eta_probe.tolist(),
                "C_p_reference": handoff["C_p_reference"].tolist(),
                "C_p_reference_X": handoff["C_p_reference_X"].tolist(),
                "C_p_reference_eta": handoff["C_p_reference_eta"].tolist(),
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
    pressure = KokunoPA10ReferencePressureIncrementToX100()
    payload = pressure.save_report(args.output)
    if args.config_output is not None:
        pressure.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("X100_handoff_probe=", payload["X100_handoff_probe"])


if __name__ == "__main__":
    _main()
