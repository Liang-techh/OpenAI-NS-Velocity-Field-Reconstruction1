"""Freeze the admitted closure gate for the current A1 PA.17 I1 candidate.

This is a narrow governance/execution wrapper around
``KokunoPA16CurrentCartesianI1Repair`` from A1 PR #1043.  It does not change
any profile, repair coefficient, similarity-coordinate map, or Cartesian
velocity formula.  It closes one concrete post-selection freedom identified by
A9 CR-A9-101: the parent schema serialized ``closure_tolerance`` and the
underlying PA.17 family accepts values as loose as 1e-5 even though the current
lineage was constructed and declared with the held-out midpoint closure gate
5e-7.

The admitted child therefore has exactly one closure gate, 5e-7.  Both direct
construction and configuration load fail closed if the wrapped parent or the
serialized frozen-gate record carries any other value.  Deliberately looser or
tighter gates require a different experiment/schema identity rather than a
post-hoc edit of this current lineage.

No statement in this module upgrades the public Kokuno reconstruction to a
paper-exact/OpenAI field or turns replay/closure checks into PDE validation.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_i1_repair import (
    KokunoPA16CurrentCartesianI1Repair,
)

SCHEMA = "kokuno-pa16-current-cartesian-i1-frozen-gate-v1"
PARENT_EXACT_HEAD = "5dc51d8d7ca48fed2a68d597f6e1a80e23becac2"
I1_CLOSURE_TOLERANCE = 5.0e-7

_TRUTH_BOUNDARY = {
    "parent_current_I1_candidate_consumed": True,
    "current_I1_closure_gate_frozen": True,
    "configuration_mutation_of_I1_closure_gate_fails_closed": True,
    "velocity_formula_changed_by_this_increment": False,
    "PA17_coefficients_changed_by_this_increment": False,
    "source_schedule_changed_by_this_increment": False,
    "current_I2_overlay_applied": False,
    "current_I3_overlay_applied": False,
    "current_I4_overlay_applied": False,
    "outer_global_leading_velocity_materialized": False,
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


def _require_frozen_gate(value: Any, *, location: str) -> float:
    try:
        gate = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{location} must equal the frozen I1 closure gate") from exc
    if not np.isfinite(gate) or gate != I1_CLOSURE_TOLERANCE:
        raise ValueError(
            f"{location} must equal the frozen I1 closure gate "
            f"{I1_CLOSURE_TOLERANCE:.17g}; got {gate!r}"
        )
    return gate


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianI1FrozenGate:
    """Current I1 Cartesian candidate with an immutable 5e-7 closure gate."""

    candidate: KokunoPA16CurrentCartesianI1Repair = field(
        default_factory=KokunoPA16CurrentCartesianI1Repair,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.candidate, KokunoPA16CurrentCartesianI1Repair):
            raise TypeError("candidate must be KokunoPA16CurrentCartesianI1Repair")
        _require_frozen_gate(
            self.candidate.closure_tolerance,
            location="candidate.closure_tolerance",
        )
        config = self.candidate.configuration()
        repair = config.get("i1_repair")
        if not isinstance(repair, Mapping):
            raise ValueError("wrapped candidate is missing i1_repair configuration")
        _require_frozen_gate(
            repair.get("closure_tolerance"),
            location="candidate configuration closure_tolerance",
        )

    @property
    def closure_tolerance(self) -> float:
        return I1_CLOSURE_TOLERANCE

    @property
    def pre_i1(self):
        return self.candidate.pre_i1

    @property
    def current(self):
        return self.candidate.current

    @property
    def geometry(self):
        return self.candidate.geometry

    @property
    def i1_family(self):
        return self.candidate.i1_family

    @property
    def A(self) -> float:
        return self.candidate.A

    @property
    def D(self) -> float:
        return self.candidate.D

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.candidate.eta_interval

    @property
    def X_I1_start(self) -> float:
        return self.candidate.X_I1_start

    @property
    def X_I1_end(self) -> float:
        return self.candidate.X_I1_end

    @property
    def log_X_I1_start(self) -> float:
        return self.candidate.log_X_I1_start

    @property
    def log_X_I1_end(self) -> float:
        return self.candidate.log_X_I1_end

    def similarity_coordinates(self, x: Any, y: Any, z: Any, t: Any):
        return self.candidate.similarity_coordinates(x, y, z, t)

    def similarity_profile_values(self, X: Any, eta: Any):
        return self.candidate.similarity_profile_values(X, eta)

    def similarity_radial_derivatives(self, X: Any, eta: Any):
        return self.candidate.similarity_radial_derivatives(X, eta)

    def values(self, x: Any, y: Any, z: Any, t: Any):
        return self.candidate.values(x, y, z, t)

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        return self.candidate.velocity(x, y, z, t)

    def i1_exit_report(self, eta: Any):
        return self.candidate.i1_exit_report(eta)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_i1_candidate": self.candidate.configuration(),
            "frozen_i1_closure_gate": {
                "closure_tolerance": I1_CLOSURE_TOLERANCE,
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianI1FrozenGate":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current I1 frozen-gate schema")
        frozen = payload.get("frozen_i1_closure_gate")
        if not isinstance(frozen, Mapping):
            raise ValueError("missing frozen_i1_closure_gate configuration")
        if frozen.get("mutable") is not False:
            raise ValueError("current I1 closure gate must be serialized as immutable")
        _require_frozen_gate(
            frozen.get("closure_tolerance"),
            location="serialized frozen I1 closure_tolerance",
        )
        parent = payload.get("parent_i1_candidate")
        if not isinstance(parent, Mapping):
            raise ValueError("missing parent_i1_candidate configuration")
        repair = parent.get("i1_repair")
        if not isinstance(repair, Mapping):
            raise ValueError("parent candidate is missing i1_repair configuration")
        _require_frozen_gate(
            repair.get("closure_tolerance"),
            location="serialized parent I1 closure_tolerance",
        )
        candidate = KokunoPA16CurrentCartesianI1Repair.from_configuration(parent)
        return cls(candidate=candidate)

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianI1FrozenGate":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "configuration": self.configuration(),
            "dependencies": {
                "parent_exact_head": PARENT_EXACT_HEAD,
                "parent_semantic_sha256": self.candidate.semantic_sha256,
            },
            "frozen_construction_gate": {
                "I1_closure_tolerance": I1_CLOSURE_TOLERANCE,
                "post_selection_mutation_allowed": False,
            },
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        base_report = self.candidate.report()
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.candidate.semantic_sha256,
            "semantic_sha256": self.semantic_sha256,
            "frozen_i1_closure_gate": {
                "closure_tolerance": I1_CLOSURE_TOLERANCE,
                "mutable": False,
                "parent_runtime_value": self.candidate.closure_tolerance,
            },
            "velocity_equivalence_to_parent": True,
            "parent_i1_exit_report": base_report["exit_report"],
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "this hardens the current I1 construction gate but changes no velocity formula",
                "the parent #1035 finite-N loop and PA.17 interpolation remain repository-autonomous",
                "I2/I3/I4 and terminal/global outer completion are not applied",
                "no matched pressure, restricted forcing, held-out complete NS residual, or PDE validation is supplied",
            ],
        }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()
    candidate = KokunoPA16CurrentCartesianI1FrozenGate()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
