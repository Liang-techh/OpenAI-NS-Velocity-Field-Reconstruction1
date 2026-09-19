"""Bind the source-defined PA.10 logarithmic length L_0 to executable Lambda.

Pinned public provenance is KokunoYumeto/yang-mills-interacting-workbench at
``143f6773feb424ad9ed3a8d116653200f20346b7``,
``navier-stokes/navier_stokes_workbench.tex``, corrected reader dated
2026-09-09 (Zenodo 22678406).

The corrected reconstruction fixes the rescaling ``Lambda`` before the PA.10
bounds and defines

    X_0 = 4 / Lambda,
    X_i = 110,
    L_0 = log(X_i / X_0) = log(X_i Lambda / 4).

This module binds that exact dependency to the repository's selected executable
``KokunoSourceRescaledCoreSeed``.  The seed's default Lambda is an explicitly
autonomous conditioning choice, not a recovered source threshold, so binding
L_0 does *not* certify the full source B_0 or T_sh estimate.  The remaining
pre-C B_0 dependencies

    ||log phi_* + log Phi(4,.)||_{C^0_eta},  M_0

remain fail closed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from .kokuno_rescaled_core_seed import KokunoSourceRescaledCoreSeed


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-pa10-l0-binding-v1"
X_I = 110.0

_SOURCE_FORMULAS = {
    "inner_radius": "X_0=4/Lambda",
    "join_radius": "X_i=110",
    "log_length": "L_0=log(X_i/X_0)=log(X_i*Lambda/4)",
    "B0_context": (
        "||log phi(X_i)||_{C^0_eta} <= "
        "||log phi_*+log Phi(4,.)||_{C^0_eta} + 0.5(M_0+0.8)L_0"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_X0_equals_4_over_Lambda_executable": True,
    "source_L0_identity_executable": True,
    "selected_executable_Lambda_bound_into_L0": True,
    "selected_Lambda_is_autonomous_conditioning_choice": True,
    "source_hidden_Lambda_recovered": False,
    "combined_profile_C0_norm_machine_bound": False,
    "source_M0_machine_bound": False,
    "source_L0_machine_bound": True,
    "source_B0_dependencies_machine_bound": False,
    "source_B0_analytic_bound_proved": False,
    "source_T_sh_lower_bound_verified": False,
    "selected_pa16_handoff_allowed": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
    "formal_full_domain_pde_gate_assessed": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class KokunoPA10L0Binding:
    """Exact PA.10 L_0 identity evaluated on an executable selected Lambda."""

    seed: KokunoSourceRescaledCoreSeed = field(default_factory=KokunoSourceRescaledCoreSeed)
    X_i: float = X_I

    def __post_init__(self) -> None:
        if not isinstance(self.seed, KokunoSourceRescaledCoreSeed):
            raise TypeError("seed must be a KokunoSourceRescaledCoreSeed")
        X_i = float(self.X_i)
        if not math.isfinite(X_i) or X_i <= 0.0:
            raise ValueError("X_i must be finite and positive")
        object.__setattr__(self, "X_i", X_i)

    @property
    def selected_Lambda(self) -> float:
        value = float(self.seed.rescaling_lambda)
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("selected Lambda must be finite and positive")
        return value

    @property
    def log_Lambda(self) -> float:
        return math.log(self.selected_Lambda)

    @property
    def X_0(self) -> float:
        value = 4.0 / self.selected_Lambda
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("materialized X_0 is outside positive binary64 range")
        return value

    @property
    def log_X_0(self) -> float:
        # Evaluate from logs rather than materializing a ratio inside log().
        return math.log(4.0) - self.log_Lambda

    @property
    def L_0(self) -> float:
        # Source identity, evaluated in the stable form log(X_i/4)+log Lambda.
        value = math.log(self.X_i / 4.0) + self.log_Lambda
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("selected source L_0 must be finite and nonnegative")
        return value

    def verify_identity(self) -> dict[str, float | bool]:
        materialized = math.log(self.X_i / self.X_0)
        stable = self.L_0
        error = abs(materialized - stable)
        tolerance = 8.0 * max(math.ulp(materialized), math.ulp(stable))
        return {
            "selected_Lambda": self.selected_Lambda,
            "log_Lambda": self.log_Lambda,
            "X_0": self.X_0,
            "log_X_0": self.log_X_0,
            "X_i": self.X_i,
            "L_0_stable": stable,
            "L_0_materialized_replay": materialized,
            "identity_abs_error": error,
            "identity_within_8ulp": bool(error <= tolerance),
        }

    def report(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "selected_seed_sha256": self.seed.sha256,
            "binding": self.verify_identity(),
            "resolved_source_B0_dependencies": ["L_0"],
            "unresolved_source_B0_dependencies": [
                "||log phi_*+log Phi(4,.)||_{C^0_eta}",
                "M_0",
            ],
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
        }

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": self.report()["source"],
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "seed": self.seed.to_payload(),
            "parameters": {"X_i": self.X_i},
            "truth_boundary": copy.deepcopy(_TRUTH_BOUNDARY),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoPA10L0Binding":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected PA.10 L0 binding schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("PA.10 L0 source formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("PA.10 L0 truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("PA.10 L0 payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing PA.10 L0 parameters")
        seed_payload = payload.get("seed")
        seed = KokunoSourceRescaledCoreSeed.from_payload(seed_payload)
        obj = cls(seed=seed, **params)
        if obj.to_payload() != payload:
            raise ValueError("PA.10 L0 replay changed payload")
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
    def load_json(cls, path: str | Path) -> "KokunoPA10L0Binding":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Bind source PA.10 L0 to the selected executable Lambda"
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)

    binding = KokunoPA10L0Binding()
    report = binding.report()
    report["binding_sha256"] = binding.sha256
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
