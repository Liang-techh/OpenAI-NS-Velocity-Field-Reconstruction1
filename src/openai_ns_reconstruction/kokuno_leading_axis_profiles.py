"""Executable axis data from Kokuno's public Navier--Stokes reconstruction.

This module intentionally implements only the prescribed leading *axis* data that
are explicit in the pinned public workbench source.  It does not reconstruct the
full leading profile, does not identify OpenAI's hidden field, and does not imply
Navier--Stokes validation.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
CORRECTED_PDF = "released_ns_reader_corrected_20260909.pdf"
SCHEMA = "kokuno-leading-axis-profile-v1"

_SOURCE_FORMULAS = {
    "A": "1/2 + h",
    "D": "1/2 - h",
    "d": "1 - eta^2",
    "U_star": "4*eta + j0",
    "H_star": "D*eta + d*U_star",
    "W_star": "1 - d*dU_star_deta - 2*D*eta*U_star",
    "minus_W_star_simplified": "3 - 8*h*eta^2 + (1 - 2*h)*j0*eta",
    "dH_star_deta": "D - 2*eta*U_star + 4*d",
    "dW_star_deta": "16*h*eta - (1 - 2*h)*j0",
}

_COORDINATE_CONTRACT = {
    "native_similarity_coordinate": "eta = z / q^(1/2-h)",
    "native_time_relation": "tau=1-t=q(1-eta^2)=q-z^2*q^(2h)",
    "current_repo_eq45_user_relation": "q-z^2*q^(-2h)=1-t",
    "mapping_status": "pending_explicit_convention_reconciliation",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "full_leading_profile_reconstructed": False,
    "full_3d_velocity_candidate": False,
    "velocity_api_compatible": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoLeadingAxisProfile:
    """Low-dimensional executable seed for Kokuno's prescribed axis profiles.

    ``h`` and ``j0`` are bounded by the public source.  The defaults are autonomous
    demo choices inside those intervals; they are not recovered OpenAI parameters.
    There is deliberately no amplitude parameter, so this layer cannot collapse to
    the zero profile by scaling.
    """

    h: float = 0.005
    j0: float = 0.025

    def __post_init__(self) -> None:
        h = float(self.h)
        j0 = float(self.j0)
        if not np.isfinite(h) or not (0.0 < h <= 1.0e-2):
            raise ValueError("h must be finite and satisfy 0 < h <= 1e-2")
        if not np.isfinite(j0) or not (0.0 < j0 <= 5.0e-2):
            raise ValueError("j0 must be finite and satisfy 0 < j0 <= 0.05")
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "j0", j0)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @staticmethod
    def _eta_array(eta: Any) -> np.ndarray:
        values = np.asarray(eta, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("eta must contain only finite values")
        return values

    def d(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        return 1.0 - eta_array * eta_array

    def u_star(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        return 4.0 * eta_array + self.j0

    def du_star_deta(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        return np.full_like(eta_array, 4.0, dtype=float)

    def h_star(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        d = 1.0 - eta_array * eta_array
        u = 4.0 * eta_array + self.j0
        return self.D * eta_array + d * u

    def dh_star_deta(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        d = 1.0 - eta_array * eta_array
        u = 4.0 * eta_array + self.j0
        return self.D - 2.0 * eta_array * u + 4.0 * d

    def w_star(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        d = 1.0 - eta_array * eta_array
        u = 4.0 * eta_array + self.j0
        return 1.0 - 4.0 * d - 2.0 * self.D * eta_array * u

    def w_star_simplified(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        return -3.0 + 8.0 * self.h * eta_array**2 - (1.0 - 2.0 * self.h) * self.j0 * eta_array

    def dw_star_deta(self, eta: Any) -> np.ndarray:
        eta_array = self._eta_array(eta)
        return 16.0 * self.h * eta_array - (1.0 - 2.0 * self.h) * self.j0

    def evaluate(self, eta: Any) -> dict[str, np.ndarray]:
        eta_array = self._eta_array(eta)
        return {
            "eta": eta_array.copy(),
            "U_star": self.u_star(eta_array),
            "dU_star_deta": self.du_star_deta(eta_array),
            "H_star": self.h_star(eta_array),
            "dH_star_deta": self.dh_star_deta(eta_array),
            "W_star": self.w_star(eta_array),
            "dW_star_deta": self.dw_star_deta(eta_array),
        }

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
                "corrected_pdf": CORRECTED_PDF,
                "formula_evidence": "pinned_public_workbench_commit",
                "pdf_independently_parsed": False,
            },
            "parameters": {
                "h": self.h,
                "j0": self.j0,
                "origin": "autonomous_demo_within_public_source_bounds",
                "source_bounds": {"h": "0<h<=1e-2", "j0": "0<j0<=0.05"},
            },
            "formulas": dict(_SOURCE_FORMULAS),
            "coordinates": dict(_COORDINATE_CONTRACT),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoLeadingAxisProfile":
        if not isinstance(payload, dict):
            raise ValueError("profile payload must be a JSON object")
        expected_keys = {"schema", "source", "parameters", "formulas", "coordinates", "truth_boundary"}
        keys_without_sha = set(payload) - {"sha256"}
        if keys_without_sha != expected_keys:
            raise ValueError("profile payload schema keys do not match")
        if payload["schema"] != SCHEMA:
            raise ValueError("unsupported Kokuno leading-axis schema")

        expected_source = cls().to_payload()["source"]
        if payload["source"] != expected_source:
            raise ValueError("Kokuno source/provenance metadata changed")
        if payload["formulas"] != _SOURCE_FORMULAS:
            raise ValueError("Kokuno source formula metadata changed")
        if payload["coordinates"] != _COORDINATE_CONTRACT:
            raise ValueError("Kokuno coordinate-contract metadata changed")
        if payload["truth_boundary"] != _TRUTH_BOUNDARY:
            raise ValueError("Kokuno truth-boundary metadata changed")

        parameters = payload["parameters"]
        if set(parameters) != {"h", "j0", "origin", "source_bounds"}:
            raise ValueError("Kokuno parameter metadata changed")
        expected_origin = "autonomous_demo_within_public_source_bounds"
        expected_bounds = {"h": "0<h<=1e-2", "j0": "0<j0<=0.05"}
        if parameters["origin"] != expected_origin or parameters["source_bounds"] != expected_bounds:
            raise ValueError("Kokuno parameter provenance metadata changed")

        profile = cls(h=float(parameters["h"]), j0=float(parameters["j0"]))
        if "sha256" in payload and payload["sha256"] != profile.sha256:
            raise ValueError("Kokuno leading-axis payload SHA mismatch")
        return profile

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoLeadingAxisProfile":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_payload(payload)
