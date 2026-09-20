"""Executable public PA.10 reference plateau through ``X_b^ref=100``.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

The corrected continuation derivation fixes

    X_0 = 4/Lambda,   X_b^ref = 100,   X_i = 110,
    y = log(X/X_0),

uses the first reference collar already materialized by Agent 1, and then states
that after ``y=2 t_1`` the reference pair ``phi_r,U_r`` is held constant in
``y`` until later source stages modify the *actual* profile.  Since the source
normalization ``C`` is constant, ``F_r=phi_r/C`` is constant on the same
plateau, while ``E_r=sqrt(2X)F_r`` keeps its geometric ``sqrt(X)`` factor.

This module closes only that reference-profile backbone.  It deliberately does
*not* invent the reference stress primitives ``p_{1,r},n_{s,r}``; those are
forward-integrated source quantities and are needed before the subsequent
fixed-``kappa_0`` actual continuation can be materialized truthfully.  No
pressure fit, forcing fit, residual optimization, arbitrary taper, or hidden
OpenAI parameter is introduced here.
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
from .kokuno_pa10_reference_continuation_collar import (
    SOURCE_X_B_REF,
    SOURCE_X_I,
    KokunoPA10ReferenceContinuationCollar,
)


SCHEMA = "kokuno-pa10-reference-plateau-x100-v1"

_SOURCE_FORMULAS = {
    "reference_geometry": "X_0=4/Lambda; X_b^ref=100; X_i=110; y=log(X/X_0)",
    "reference_first_collar": (
        "retain natural profile through y=t_1; flatten its radial slope on "
        "t_1<y<2t_1 with the public flat step"
    ),
    "reference_plateau": (
        "after 2t_1 keep phi_r,U_r constant in y until later source stages"
    ),
    "normalized_plateau": (
        "C is constant, hence F_r=phi_r/C is constant in y after 2t_1"
    ),
    "azimuthal_profile": "E_r=sqrt(2X) F_r",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_reference_plateau_to_X100_executable": True,
    "reference_F_U_vectorized_to_X100": True,
    "reference_radial_jets_executable_to_X100": True,
    "reference_X100_handoff_executable": True,
    "reference_plateau_is_actual_stressed_profile": False,
    "reference_stress_primitives_p1r_nsr_materialized": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "kappa0_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
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


@dataclass(frozen=True)
class KokunoPA10ReferencePlateauToX100:
    """Source reference pair continued exactly flat in ``y`` through ``X=100``."""

    collar: KokunoPA10ReferenceContinuationCollar = field(
        default_factory=KokunoPA10ReferenceContinuationCollar,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.collar, KokunoPA10ReferenceContinuationCollar):
            raise TypeError("collar must be KokunoPA10ReferenceContinuationCollar")
        if not self.collar.X_2 < SOURCE_X_B_REF < SOURCE_X_I:
            raise ValueError("source reference scales are not ordered X_2 < 100 < 110")

    @property
    def X_0(self) -> float:
        return float(self.collar.X_0)

    @property
    def X_2(self) -> float:
        return float(self.collar.X_2)

    @property
    def X_b_ref(self) -> float:
        return float(SOURCE_X_B_REF)

    @property
    def X_i(self) -> float:
        return float(SOURCE_X_I)

    @property
    def Lambda(self) -> float:
        return float(self.collar.Lambda)

    @property
    def C(self) -> float:
        return float(self.collar.C)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.collar.eta_interval

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

    def _plateau_endpoint(self, eta: np.ndarray) -> dict[str, np.ndarray]:
        X2 = np.full_like(eta, self.X_2, dtype=float)
        return self.collar.values(X2, eta)

    def values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return vectorized ``F_r,U_r,E_r`` on ``0<=X<=100``.

        The parent collar is used byte-for-byte on ``X<=X_2``.  For larger X,
        ``F_r`` and ``U_r`` are the exact endpoint values at ``X_2`` and only
        the public geometric factor in ``E_r=sqrt(2X)F_r`` changes.
        """
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F = np.empty_like(x)
        U = np.empty_like(x)

        collar_mask = x <= self.X_2
        if np.any(collar_mask):
            inherited = self.collar.values(x[collar_mask], e[collar_mask])
            F[collar_mask] = inherited["F_reference"]
            U[collar_mask] = inherited["U_reference"]

        plateau_mask = ~collar_mask
        if np.any(plateau_mask):
            endpoint = self._plateau_endpoint(e[plateau_mask])
            F[plateau_mask] = endpoint["F_reference"]
            U[plateau_mask] = endpoint["U_reference"]

        if np.any(F <= 0.0):
            raise RuntimeError("reference F must remain strictly positive through X=100")
        E = np.sqrt(2.0 * x) * F
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_reference": F.reshape(shape),
            "U_reference": U.reshape(shape),
            "E_reference": E.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return exact first X-derivatives of the materialized reference pair."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F_X = np.zeros_like(x)
        U_X = np.zeros_like(x)

        collar_mask = x <= self.X_2
        if np.any(collar_mask):
            inherited = self.collar.radial_derivatives(x[collar_mask], e[collar_mask])
            F_X[collar_mask] = inherited["F_reference_X"]
            U_X[collar_mask] = inherited["U_reference_X"]

        values = self.values(x, e)
        E_X = np.empty_like(x)
        positive = x > 0.0
        E_X[positive] = (
            values["F_reference"][positive] / np.sqrt(2.0 * x[positive])
            + np.sqrt(2.0 * x[positive]) * F_X[positive]
        )
        # E=sqrt(2X)F has the standard one-sided sqrt singular derivative at
        # X=0 as a scalar E-profile.  Cartesian regularity is carried by F,
        # exactly as in the public source.  Do not manufacture a finite E_X(0).
        E_X[~positive] = np.inf
        return {
            "F_reference_X": F_X.reshape(shape),
            "U_reference_X": U_X.reshape(shape),
            "E_reference_X": E_X.reshape(shape),
        }

    def radial_log_slopes(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Expose source-useful radial slope quantities without new parameters."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.values(X_arr, eta_arr)
        derivatives = self.radial_derivatives(X_arr, eta_arr)
        F = values["F_reference"]
        D_X_log_F = X_arr * derivatives["F_reference_X"] / F
        D_X_U = X_arr * derivatives["U_reference_X"]
        # H=sqrt(2X)E=2 X F, hence l=X d_X log H=1+X F_X/F.
        l_reference = 1.0 + D_X_log_F
        return {
            "D_X_log_F_reference": D_X_log_F,
            "D_X_U_reference": D_X_U,
            "l_reference": l_reference,
        }

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        """Return the exact reference value/jet handoff at ``X_b^ref=100``."""
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_b_ref, dtype=float)
        return {
            **self.values(X, eta_arr),
            **self.radial_derivatives(X, eta_arr),
            **self.radial_log_slopes(X, eta_arr),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source_X_b_ref": self.X_b_ref,
            "source_X_i": self.X_i,
            "parent_reference_collar_semantic_sha256": self.collar.semantic_sha256,
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA10ReferencePlateauToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "reference-plateau configuration is frozen; serialized data differ "
                "from the registered source/autonomous realization"
            )
        return cls()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10ReferencePlateauToX100":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
        handoff = self.handoff_at_X100(eta_probe)
        endpoint = self.collar.values(
            np.full_like(eta_probe, self.X_2), eta_probe
        )
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
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "geometry": {
                "Lambda": self.Lambda,
                "X_0": self.X_0,
                "X_2": self.X_2,
                "X_b_ref": self.X_b_ref,
                "X_i": self.X_i,
            },
            "X100_handoff_probe": {
                "eta": eta_probe.tolist(),
                "F_reference": handoff["F_reference"].tolist(),
                "U_reference": handoff["U_reference"].tolist(),
                "E_reference": handoff["E_reference"].tolist(),
                "F_X": handoff["F_reference_X"].tolist(),
                "U_X": handoff["U_reference_X"].tolist(),
                "l_reference": handoff["l_reference"].tolist(),
                "max_abs_F_plateau_value_drift": float(
                    np.max(np.abs(handoff["F_reference"] - endpoint["F_reference"]))
                ),
                "max_abs_U_plateau_value_drift": float(
                    np.max(np.abs(handoff["U_reference"] - endpoint["U_reference"]))
                ),
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
    profile = KokunoPA10ReferencePlateauToX100()
    payload = profile.save_report(args.output)
    if args.config_output is not None:
        profile.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("geometry=", payload["geometry"])
    print("X100_handoff_probe=", payload["X100_handoff_probe"])


if __name__ == "__main__":
    _main()
