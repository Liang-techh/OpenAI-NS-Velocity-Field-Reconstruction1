"""Fail-closed CR002 audit for the #868 axis-completed first-cell quadrature.

This audit keeps two statements distinct:

1. PR #868 now includes a discrete first cell from the virtual axis node to the
   first positive radius; and
2. that trapezoidal first cell is not, by itself, a certificate of the formal
   integral from zero to the first positive radius.

The distinction matters most for the e=2 tangential channel.  For the mechanics
source P_e(s)=1, the formal first-cell primitive is r_min^(e+1)/(e+1), whereas
#868's virtual-endpoint trapezoid returns 0.5*r_min^(e+1).  These agree for e=1
but the e=2 trapezoid is 50 percent high.  This is a quadrature-scope witness,
not candidate evidence and not a claim that the real typed source is constant.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_strict_inner_transport_axis_completed_radial_stress import (
    _axis_completed_cumulative_trapezoid,
)

TASK = "CR002-KOKUNO-AXIS-COMPLETED-FIRST-CELL-QUADRATURE-SCOPE-118"
SCHEMA = "cr002-kokuno-axis-completed-first-cell-quadrature-scope-v1"
UPSTREAM_PR = 868
UPSTREAM_HEAD = "c8268100a36942de9fa966d76f22ad6148ceaa02"
UPSTREAM_BLOB = "041d4298b5350b25a5699492ddf39c35fae4aaab"
CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = REPO_ROOT / "configs" / "kokuno_axis_completed_first_cell_quadrature_scope.json"
UPSTREAM_PATH = (
    REPO_ROOT
    / "src"
    / "openai_ns_reconstruction"
    / "kokuno_strict_inner_transport_axis_completed_radial_stress.py"
)
CONSTRAINTS_PATH = REPO_ROOT / "configs" / "constraints.json"

_AUTONOMOUS_MARKER = (
    "PR #868 prepends a virtual axis node with weighted endpoint value zero and applies a trapezoidal first cell"
)
_PUBLIC_FORMAL_MARKER = (
    "the registered compact radial inverse has formal lower limit r=0: sigma_e(r)=-r^(-e) integral_0^r s^e P_e(F)(s) ds"
)


def _git_blob_sha1(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object at {path}")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _close(actual: float, expected: float, *, atol: float = 2.0e-15) -> bool:
    return math.isclose(float(actual), float(expected), rel_tol=2.0e-14, abs_tol=atol)


def _mechanics_witness(r_min: float, exponent: int) -> dict[str, float]:
    if exponent not in (1, 2):
        raise ValueError("scope witness is frozen to e=1,2")
    # The upstream helper requires at least three positive nodes.  Only its
    # first returned primitive is used here; later nodes cannot affect it.
    radii = np.asarray([r_min, 1.5 * r_min, 2.0 * r_min], dtype=float)
    weighted = radii**exponent
    discrete = float(_axis_completed_cumulative_trapezoid(weighted, radii)[0])
    formal = float(r_min ** (exponent + 1) / (exponent + 1))
    discrete_stress = float(-discrete / (r_min**exponent))
    formal_stress = float(-formal / (r_min**exponent))
    relative_overestimate = float((discrete - formal) / formal)
    return {
        "formal_primitive": formal,
        "trapezoid_primitive": discrete,
        "formal_stress_at_r_min": formal_stress,
        "trapezoid_stress_at_r_min": discrete_stress,
        "relative_primitive_overestimate": relative_overestimate,
    }


def _audit_provenance(contract: dict[str, Any]) -> None:
    provenance = contract.get("provenance")
    _require(isinstance(provenance, dict), "missing four-way provenance")
    required = {"user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"}
    _require(set(provenance) == required, "provenance buckets must remain exactly four-way")
    for key in required:
        _require(isinstance(provenance[key], list) and provenance[key], f"empty provenance bucket: {key}")
    all_entries: list[str] = []
    for key in required:
        for entry in provenance[key]:
            _require(isinstance(entry, str) and entry.strip() != "", f"invalid provenance entry in {key}")
            all_entries.append(entry)
    _require(len(all_entries) == len(set(all_entries)), "provenance entries must not be duplicated across buckets")
    _require(_AUTONOMOUS_MARKER in provenance["autonomous_design"], "#868 first-cell realization must remain autonomous design")
    _require(_AUTONOMOUS_MARKER not in provenance["public_source_fact"], "repository trapezoid realization cannot be laundered into public-source fact")
    _require(_PUBLIC_FORMAL_MARKER in provenance["public_source_fact"], "formal radial-inverse source statement drifted")


def _audit_cr001(contract: dict[str, Any]) -> None:
    snapshot = contract.get("cr001_snapshot")
    _require(isinstance(snapshot, dict), "missing CR001 snapshot")
    _require(snapshot.get("constraints_blob_sha1") == CONSTRAINTS_BLOB, "canonical constraints blob drift")
    _require(_git_blob_sha1(CONSTRAINTS_PATH) == CONSTRAINTS_BLOB, "working-tree canonical constraints are not the frozen CR001 blob")
    canonical = _load_json(CONSTRAINTS_PATH)
    _require(canonical["nu"] == 0.01 == snapshot.get("nu"), "nu drift")
    _require(canonical["domain"]["physical"] == snapshot.get("physical_domain") == "R^3", "physical-domain drift")
    _require(canonical["domain"]["evaluation_box"] == snapshot.get("evaluation_box"), "evaluation-box drift")
    _require(canonical["domain"]["support"] == snapshot.get("support"), "support drift")
    _require(canonical["domain"]["time_interval"] == snapshot.get("time_interval"), "time-interval drift")
    _require(canonical["forcing"]["mode"] == snapshot.get("forcing_mode") == "restricted_two_parameter_family", "forcing-family drift")
    _require(canonical["nontriviality"]["reference_energy"] == snapshot.get("reference_energy") == 1.0, "energy normalization drift")
    _require(canonical["nontriviality"]["reference_energy_abs_tolerance"] == snapshot.get("reference_energy_abs_tolerance") == 0.001, "energy tolerance drift")
    validation = canonical["validation"]
    _require(validation["seed"] == snapshot.get("validation_seed") == 914027, "validation seed drift")
    _require(validation["held_out_points"] == snapshot.get("held_out_points") == 4096, "held-out count drift")
    _require(validation["derivative_steps"] == snapshot.get("derivative_steps"), "derivative ladder drift")
    _require(validation["quadrature_orders_per_axis"] == snapshot.get("quadrature_orders_per_axis"), "quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key, expected in (
        ("divergence_max", 1.0e-5),
        ("divergence_L2", 1.0e-5),
        ("pde_residual_max", 1.0e-3),
        ("pde_residual_L2", 1.0e-3),
    ):
        _require(thresholds[key] == expected == snapshot.get(key), f"CR001 threshold drift: {key}")
    _require(snapshot.get("residual_defined_pointwise_free_force_allowed") is False, "free residual-defined forcing must remain forbidden")
    _require(snapshot.get("amplitude_collapse_allowed") is False, "amplitude collapse must remain forbidden")
    _require(snapshot.get("post_hoc_threshold_relaxation_allowed") is False, "post-hoc threshold relaxation must remain forbidden")


def audit_contract(contract: dict[str, Any] | None = None) -> dict[str, Any]:
    """Audit the CR002 scope contract and return a deterministic mechanics receipt."""
    data = copy.deepcopy(_load_json(CONTRACT_PATH) if contract is None else contract)
    _require(data.get("schema") == SCHEMA, "schema drift")
    _require(data.get("task") == TASK, "task drift")

    upstream = data.get("audited_upstream")
    _require(isinstance(upstream, dict), "missing upstream binding")
    _require(upstream.get("pr") == UPSTREAM_PR, "upstream PR drift")
    _require(upstream.get("head") == UPSTREAM_HEAD, "upstream exact-head drift")
    _require(upstream.get("git_blob_sha1") == UPSTREAM_BLOB, "upstream blob registration drift")
    _require(_git_blob_sha1(UPSTREAM_PATH) == UPSTREAM_BLOB, "audited upstream source blob drifted")

    source_text = UPSTREAM_PATH.read_text()
    for snippet in (
        "augmented_r = np.concatenate((np.array([0.0]), r))",
        "augmented_values = np.concatenate((np.array([0.0]), values))",
        "increments = 0.5 * (augmented_values[1:] + augmented_values[:-1])",
        '"formal_axis_based_inverse_verified": False',
        '"formal_full_domain_moment_verified": False',
        '"axis_regularity_verified": False',
    ):
        _require(snippet in source_text, f"upstream first-cell/truth-boundary source drift: {snippet}")

    _audit_provenance(data)
    scope = data.get("representation_scope")
    _require(isinstance(scope, dict), "missing representation scope")
    _require(scope.get("virtual_axis_node_is_quadrature_realization_not_formal_integral_certificate") is True, "virtual endpoint must remain a numerical realization")
    _require(scope.get("first_cell_convergence_required_before_formal_lower_limit_accuracy_claim") is True, "formal first-cell accuracy must require separate convergence/error evidence")
    _require(scope.get("axis_completed_discrete_lower_limit_may_be_claimed") is True, "#868 scoped discrete repair must not be denied")
    for false_key in (
        "formal_axis_based_inverse_verified",
        "formal_first_cell_accuracy_verified",
        "formal_full_domain_moment_verified",
        "axis_regularity_verified",
    ):
        _require(scope.get(false_key) is False, f"unsupported promotion: {false_key}")

    witness = data.get("mechanics_witness")
    _require(isinstance(witness, dict), "missing mechanics witness")
    r_min = float(witness.get("r_min"))
    _require(_close(r_min, 0.2), "mechanics witness r_min drift")
    _require(witness.get("candidate_evidence") is False, "mechanics witness cannot become candidate evidence")
    _require(witness.get("public_source_numeric_target") is False, "mechanics witness cannot become a public numerical target")

    computed: dict[str, Any] = {}
    for exponent, key in ((1, "e1"), (2, "e2")):
        actual = _mechanics_witness(r_min, exponent)
        registered = witness.get(key)
        _require(isinstance(registered, dict), f"missing {key} witness")
        for metric, value in actual.items():
            _require(_close(float(registered.get(metric)), value), f"{key} mechanics witness drift: {metric}")
        computed[key] = actual

    _require(_close(computed["e1"]["relative_primitive_overestimate"], 0.0), "e=1 witness should be trapezoid-exact for constant source")
    _require(_close(computed["e2"]["relative_primitive_overestimate"], 0.5), "e=2 witness must expose 50% first-cell overestimate")
    _require(not _close(computed["e2"]["formal_stress_at_r_min"], computed["e2"]["trapezoid_stress_at_r_min"]), "e=2 formal/discrete stresses must remain distinguishable")

    _audit_cr001(data)

    truth = data.get("truth_state")
    _require(isinstance(truth, dict), "missing truth state")
    for false_key in (
        "changes_velocity",
        "changes_pressure",
        "changes_forcing",
        "changes_candidate_bytes",
        "changes_scientific_thresholds",
        "velocity_export_ready_promoted",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
        "pde_pending_blocks_existing_callable_velocity_delivery",
    ):
        _require(truth.get(false_key) is False, f"unsupported truth-state promotion: {false_key}")

    return {
        "schema": SCHEMA,
        "task": TASK,
        "audited_upstream_pr": UPSTREAM_PR,
        "audited_upstream_head": UPSTREAM_HEAD,
        "audited_upstream_blob": UPSTREAM_BLOB,
        "e1": computed["e1"],
        "e2": computed["e2"],
        "conclusion": "axis_completed_trapezoid_is_scoped_discrete_realization_not_formal_first_cell_accuracy_certificate",
        "formal_axis_based_inverse_verified": False,
        "formal_first_cell_accuracy_verified": False,
        "pde_validated": False,
    }


def main() -> int:
    print(json.dumps(audit_contract(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
