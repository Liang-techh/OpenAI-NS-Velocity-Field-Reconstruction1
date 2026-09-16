"""Bounded compact 2D profile basis for the Eq. (4.5) candidate lane.

This module intentionally stops at ``Phi/F`` profile jets. The open
streamfunction-structure lane owns the algebraic map from ``Phi`` to ``v0/U``,
and the Eq. (4.5) velocity-backbone lane owns the Cartesian ``[u,v,w]`` mixing.
Keeping those pieces separate avoids duplicating still-open work while providing
the missing bounded, serializable profile parameterization.

All basis-shape choices in this module are autonomous project design choices,
not paper-sourced hidden coefficients.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np


SCHEMA = "eq45_compact_profile_basis_v1"
DEFAULT_COEFFICIENT_LIMIT = 4.0
DEFAULT_X_CUT = 4.0
DEFAULT_ETA_CUT = 1.0
DEFAULT_CUTOFF_POWER = 5


@dataclass(frozen=True)
class Eq45ProfileJets:
    """Vectorized profile values required by the downstream Phi/F adapter."""

    phi: np.ndarray
    phi_x: np.ndarray
    phi_eta: np.ndarray
    swirl: np.ndarray

    def stacked(self) -> np.ndarray:
        return np.stack((self.phi, self.phi_x, self.phi_eta, self.swirl), axis=-1)


@dataclass(frozen=True)
class Eq45CompactProfileBasis:
    """Tensor-product compact basis in similarity coordinates ``(X, eta)``.

    The same deterministic mode ordering is used for ``phi_coefficients`` and
    ``swirl_coefficients``::

        (i, j) for i=0..radial_degree, j=0..eta_degree.

    Each raw monomial is multiplied by autonomous C^4 compact envelopes

        (1 - X/x_cut)^p_+  and  (1 - (eta/eta_cut)^2)^p_+,

    with ``p >= 5`` by default. There is no division by ``X``; therefore finite
    coefficients give finite profile values and first derivatives on the axis.
    """

    radial_degree: int
    eta_degree: int
    phi_coefficients: tuple[float, ...]
    swirl_coefficients: tuple[float, ...]
    x_cut: float = DEFAULT_X_CUT
    eta_cut: float = DEFAULT_ETA_CUT
    cutoff_power: int = DEFAULT_CUTOFF_POWER
    coefficient_limit: float = DEFAULT_COEFFICIENT_LIMIT

    def __post_init__(self) -> None:
        if not isinstance(self.radial_degree, int) or self.radial_degree < 0:
            raise ValueError("radial_degree must be a nonnegative integer")
        if not isinstance(self.eta_degree, int) or self.eta_degree < 0:
            raise ValueError("eta_degree must be a nonnegative integer")
        if not np.isfinite(self.x_cut) or self.x_cut <= 0.0:
            raise ValueError("x_cut must be positive and finite")
        if not np.isfinite(self.eta_cut) or not 0.0 < self.eta_cut <= 1.0:
            raise ValueError("eta_cut must lie in (0, 1]")
        if not isinstance(self.cutoff_power, int) or self.cutoff_power < 5:
            raise ValueError("cutoff_power must be an integer >= 5")
        if not np.isfinite(self.coefficient_limit) or self.coefficient_limit <= 0.0:
            raise ValueError("coefficient_limit must be positive and finite")

        expected = self.mode_count
        for name, coefficients in (
            ("phi_coefficients", self.phi_coefficients),
            ("swirl_coefficients", self.swirl_coefficients),
        ):
            if len(coefficients) != expected:
                raise ValueError(f"{name} must have length {expected}")
            values = np.asarray(coefficients, dtype=float)
            if not np.all(np.isfinite(values)):
                raise ValueError(f"{name} must be finite")
            if np.any(np.abs(values) > self.coefficient_limit):
                raise ValueError(
                    f"{name} exceed declared coefficient_limit={self.coefficient_limit}"
                )

    @property
    def mode_count(self) -> int:
        return (self.radial_degree + 1) * (self.eta_degree + 1)

    @property
    def parameter_count(self) -> int:
        return 2 * self.mode_count

    @property
    def mode_indices(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (i, j)
            for i in range(self.radial_degree + 1)
            for j in range(self.eta_degree + 1)
        )

    @classmethod
    def seed(cls) -> "Eq45CompactProfileBasis":
        """Return a deterministic nonzero visualization-starting profile block.

        This is an autonomous initialization only. It is not a recovered OpenAI
        profile and is not claimed to satisfy the Navier--Stokes equations.
        """
        radial_degree = 1
        eta_degree = 2
        # Ordering: (0,0),(0,1),(0,2),(1,0),(1,1),(1,2).
        phi = (1.0, 0.0, -0.35, -0.30, 0.0, 0.10)
        swirl = (0.80, 0.0, -0.25, -0.20, 0.0, 0.05)
        return cls(
            radial_degree=radial_degree,
            eta_degree=eta_degree,
            phi_coefficients=phi,
            swirl_coefficients=swirl,
        )

    def _envelopes(
        self, x: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        sx = x / self.x_cut
        se = eta / self.eta_cut

        inside_x = (x >= 0.0) & (sx < 1.0)
        inside_eta = np.abs(se) < 1.0

        one_minus_x = np.where(inside_x, 1.0 - sx, 0.0)
        one_minus_eta2 = np.where(inside_eta, 1.0 - se * se, 0.0)

        bx = one_minus_x**self.cutoff_power
        be = one_minus_eta2**self.cutoff_power

        dbx = np.where(
            inside_x,
            -(self.cutoff_power / self.x_cut)
            * one_minus_x ** (self.cutoff_power - 1),
            0.0,
        )
        dbe = np.where(
            inside_eta,
            -(2.0 * self.cutoff_power * eta / (self.eta_cut**2))
            * one_minus_eta2 ** (self.cutoff_power - 1),
            0.0,
        )
        return bx, dbx, be, dbe

    def evaluate(self, x, eta) -> Eq45ProfileJets:
        """Evaluate ``Phi, Phi_X, Phi_eta, F`` on broadcast-compatible arrays."""
        x_arr = np.asarray(x, dtype=float)
        eta_arr = np.asarray(eta, dtype=float)
        try:
            x_arr, eta_arr = np.broadcast_arrays(x_arr, eta_arr)
        except ValueError as exc:
            raise ValueError("X and eta must be broadcast-compatible") from exc

        if not np.all(np.isfinite(x_arr)) or not np.all(np.isfinite(eta_arr)):
            raise ValueError("X and eta must be finite")
        if np.any(x_arr < 0.0):
            raise ValueError("X must be nonnegative")

        bx, dbx, be, dbe = self._envelopes(x_arr, eta_arr)
        sx = x_arr / self.x_cut
        se = eta_arr / self.eta_cut

        phi = np.zeros_like(x_arr, dtype=float)
        phi_x = np.zeros_like(x_arr, dtype=float)
        phi_eta = np.zeros_like(x_arr, dtype=float)
        swirl = np.zeros_like(x_arr, dtype=float)

        for index, (i, j) in enumerate(self.mode_indices):
            sx_i = sx**i
            se_j = se**j
            core = sx_i * se_j

            if i == 0:
                dcore_x = np.zeros_like(x_arr, dtype=float)
            else:
                dcore_x = (i / self.x_cut) * sx ** (i - 1) * se_j

            if j == 0:
                dcore_eta = np.zeros_like(x_arr, dtype=float)
            else:
                dcore_eta = (j / self.eta_cut) * sx_i * se ** (j - 1)

            basis = bx * be * core
            basis_x = be * (dbx * core + bx * dcore_x)
            basis_eta = bx * (dbe * core + be * dcore_eta)

            phi += self.phi_coefficients[index] * basis
            phi_x += self.phi_coefficients[index] * basis_x
            phi_eta += self.phi_coefficients[index] * basis_eta
            swirl += self.swirl_coefficients[index] * basis

        if not all(
            np.all(np.isfinite(value))
            for value in (phi, phi_x, phi_eta, swirl)
        ):
            raise ValueError("profile evaluation produced non-finite values")

        return Eq45ProfileJets(
            phi=phi,
            phi_x=phi_x,
            phi_eta=phi_eta,
            swirl=swirl,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "classification": {
                "basis_shape": "autonomous_design",
                "coefficient_values": "autonomous_design",
                "paper_exact": False,
            },
            "radial_degree": self.radial_degree,
            "eta_degree": self.eta_degree,
            "x_cut": self.x_cut,
            "eta_cut": self.eta_cut,
            "cutoff_power": self.cutoff_power,
            "coefficient_limit": self.coefficient_limit,
            "mode_indices": [list(item) for item in self.mode_indices],
            "phi_coefficients": list(self.phi_coefficients),
            "swirl_coefficients": list(self.swirl_coefficients),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Eq45CompactProfileBasis":
        if not isinstance(data, Mapping):
            raise ValueError("profile basis payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")

        required = (
            "radial_degree",
            "eta_degree",
            "x_cut",
            "eta_cut",
            "cutoff_power",
            "coefficient_limit",
            "phi_coefficients",
            "swirl_coefficients",
        )
        missing = [name for name in required if name not in data]
        if missing:
            raise ValueError(f"missing profile basis field(s): {missing}")

        return cls(
            radial_degree=int(data["radial_degree"]),
            eta_degree=int(data["eta_degree"]),
            x_cut=float(data["x_cut"]),
            eta_cut=float(data["eta_cut"]),
            cutoff_power=int(data["cutoff_power"]),
            coefficient_limit=float(data["coefficient_limit"]),
            phi_coefficients=tuple(float(v) for v in data["phi_coefficients"]),
            swirl_coefficients=tuple(float(v) for v in data["swirl_coefficients"]),
        )

    def save_json(self, path: str | Path) -> None:
        target = Path(path)
        target.write_text(json.dumps(self.to_dict(), indent=2) + "\n", encoding="utf-8")

    @classmethod
    def load_json(cls, path: str | Path) -> "Eq45CompactProfileBasis":
        target = Path(path)
        return cls.from_dict(json.loads(target.read_text(encoding="utf-8")))
