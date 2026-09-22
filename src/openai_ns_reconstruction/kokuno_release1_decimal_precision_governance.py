"""Fail-closed CR002 audit for the Kokuno A1 #1188 precision boundary.

This module governs representation provenance only.  It does not alter the
candidate, source schedule, forcing, validation thresholds, or scientific
readiness states.
"""

from __future__ import annotations

import copy
import hashlib
import json
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = ROOT / "configs" / "kokuno_release1_decimal_precision_scope.json"
IMPLEMENTATION_PATH = (
    ROOT
    / "src"
    / "openai_ns_reconstruction"
    / "kokuno_pa16_current_cartesian_postswirl_release1.py"
)
CONSTRAINTS_PATH = ROOT / "configs" / "constraints.json"

CANONICAL_SOURCE_CLASSES = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}


class GovernanceError(RuntimeError):
    """Raised when the pinned precision/source boundary drifts."""


def _git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GovernanceError(message)


def load_contract(path: Path = CONTRACT_PATH) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def mechanics_only_precision_witness() -> dict[str, Any]:
    """Show that Decimal embedding cannot recover a sub-ULP distinction.

    The numbers are synthetic representation mechanics only.  They are not
    Kokuno/OpenAI candidate data and are not PDE evidence.
    """

    with localcontext() as ctx:
        ctx.prec = 80
        exact_a = Decimal("0.123456789012345678901234567890123456789")
        exact_b = exact_a + Decimal("1e-30")
        f_a = float(exact_a)
        f_b = float(exact_b)
        embedded_a = Decimal.from_float(f_a)
        embedded_b = Decimal.from_float(f_b)
        return {
            "scope": "autonomous_representation_mechanics_only",
            "exact_delta": str(exact_b - exact_a),
            "binary64_values_equal": f_a == f_b,
            "decimal_embeddings_equal": embedded_a == embedded_b,
            "embedded_a": str(embedded_a),
            "embedding_error_nonzero": embedded_a != exact_a,
        }


def audit_payload(
    contract: Mapping[str, Any],
    implementation_text: str,
    implementation_bytes: bytes,
    constraints: Mapping[str, Any],
    constraints_bytes: bytes,
) -> dict[str, Any]:
    """Audit the exact representation boundary and frozen CR001 contract."""

    pinned = contract["pinned_lineage"]
    _require(contract["schema"] == "cr002-kokuno-release1-decimal-precision-scope-v1", "schema drift")
    _require(contract["task_id"] == "CR002-KOKUNO-RELEASE1-PRECISION-SCOPE", "task id drift")
    _require(
        set(contract["canonical_source_classes"]) == CANONICAL_SOURCE_CLASSES,
        "canonical source vocabulary drift",
    )
    _require(
        set(contract["source_classification"]) == CANONICAL_SOURCE_CLASSES,
        "source classification keys must be exactly the four canonical classes",
    )
    _require(
        all(contract["source_classification"][key] for key in CANONICAL_SOURCE_CLASSES),
        "each source class must remain nonempty",
    )

    _require(
        _git_blob_sha1(implementation_bytes) == pinned["implementation_git_blob_sha1"],
        "A1 #1188 implementation blob drifted; re-audit before transferring this contract",
    )
    _require(
        _git_blob_sha1(constraints_bytes) == pinned["canonical_constraints_git_blob_sha1"],
        "canonical CR001 constraints blob drifted",
    )

    required_fragments = (
        "SIGMA_QUADRATURE_ORDER = 96",
        "np.polynomial.legendre.leggauss(SIGMA_QUADRATURE_ORDER)",
        "float(np.dot(_GL_WEIGHTS, _sigma(x)))",
        "E_rel = np.asarray(seam[\"E_base\"], dtype=float)",
        "out_float = np.zeros((r_f.size, 3), dtype=float)",
        "out[idx] = _embed_decimal(out_float[idx])",
        '"source_l_minus1_hold_materialized": False',
        '"source_l_minus1_to_minus_h_transition_materialized": False',
        '"outer_global_leading_velocity_materialized": False',
        '"unified_global_cartesian_velocity_export_ready": False',
        '"restricted_forcing_materialized": False',
        '"heldout_ns_residual_assessed": False',
        '"pde_validated": False',
        '"paper_exact": False',
        '"openai_field_identified": False',
    )
    missing = [fragment for fragment in required_fragments if fragment not in implementation_text]
    _require(not missing, f"pinned A1 precision/truth mechanics drifted: {missing}")

    representation = contract["representation_identity"]
    _require(representation["declared_decimal_significant_digits"] == 96, "Decimal digit scope drift")
    _require(representation["sigma_quadrature_order"] == 96, "quadrature order scope drift")
    _require(
        representation["new_release1_internal_profile_arithmetic"] == "binary64"
        and representation["new_release1_sigma_quadrature_arithmetic"] == "binary64"
        and representation["new_release1_cartesian_assembly_arithmetic"] == "binary64",
        "current release1 arithmetic must remain classified as binary64",
    )
    _require(
        representation["decimal_embedding_stage"] == "after binary64 Cartesian component assembly",
        "Decimal embedding stage drift",
    )

    state = contract["machine_state"]
    _require(state["release1_callable_velocity_materialized"] is True, "callable stage unexpectedly demoted")
    _require(state["release1_decimal_output_encoding_materialized"] is True, "Decimal output encoding unexpectedly demoted")
    _require(state["release1_internal_binary64_arithmetic_materialized"] is True, "binary64 arithmetic fact unexpectedly demoted")
    for key in (
        "release1_end_to_end_96digit_arithmetic_materialized",
        "parent_1179_sub_epsilon_survival_evidence_transfers_to_new_release1_arithmetic",
        "binary64_to_decimal_embedding_recovers_discarded_sub_ulp_information",
        "source_exact_numerical_precision_identified",
        "unified_global_cartesian_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require(state[key] is False, f"forbidden precision/scientific promotion: {key}")

    # Bind the unchanged CR001 contract without making PDE validation a delivery gate.
    c = contract["cr001_unchanged"]
    _require(constraints["nu"] == c["nu"] == 0.01, "nu drift")
    _require(constraints["domain"]["physical"] == c["physical_domain"] == "R^3", "physical domain drift")
    _require(constraints["domain"]["evaluation_box"] == c["evaluation_box"], "evaluation box drift")
    _require(constraints["domain"]["support"] == c["support"], "support drift")
    _require(constraints["domain"]["time_interval"] == c["time_interval"], "time interval drift")
    _require(constraints["forcing"]["mode"] == c["forcing_mode"], "forcing family drift")
    _require(
        "No residual-dependent basis or pointwise free force" in constraints["forcing"]["restriction"],
        "residual-defined/free forcing prohibition drift",
    )
    _require(constraints["nontriviality"]["reference_energy"] == c["reference_energy"], "reference energy drift")
    _require(
        constraints["nontriviality"]["reference_energy_abs_tolerance"] == c["reference_energy_abs_tolerance"],
        "energy tolerance drift",
    )
    validation = constraints["validation"]
    _require(validation["seed"] == c["validation_seed"] == 914027, "validation seed drift")
    _require(validation["held_out_points"] == c["held_out_points"] == 4096, "held-out count drift")
    _require(validation["derivative_steps"] == c["derivative_steps"], "FD ladder drift")
    _require(validation["quadrature_orders_per_axis"] == c["quadrature_orders_per_axis"], "quadrature ladder drift")
    thresholds = validation["thresholds"]
    for key in ("pde_residual_max", "pde_residual_L2", "divergence_max", "divergence_L2"):
        _require(thresholds[key] == c[key], f"CR001 threshold drift: {key}")

    witness = mechanics_only_precision_witness()
    _require(witness["binary64_values_equal"], "mechanics witness no longer demonstrates binary64 collision")
    _require(witness["decimal_embeddings_equal"], "Decimal embedding unexpectedly recovered discarded distinction")
    _require(witness["embedding_error_nonzero"], "mechanics witness lost its precision distinction")

    return {
        "schema": "cr002-kokuno-release1-decimal-precision-audit-v1",
        "status": "pass",
        "a1_exact_head": pinned["a1_exact_head"],
        "implementation_git_blob_sha1": pinned["implementation_git_blob_sha1"],
        "canonical_constraints_git_blob_sha1": pinned["canonical_constraints_git_blob_sha1"],
        "decimal_output_encoding_materialized": True,
        "end_to_end_96digit_arithmetic_materialized": False,
        "global_velocity_export_promoted": False,
        "pde_promoted": False,
        "mechanics_witness": witness,
    }


def audit() -> dict[str, Any]:
    contract = load_contract()
    implementation_bytes = IMPLEMENTATION_PATH.read_bytes()
    constraints_bytes = CONSTRAINTS_PATH.read_bytes()
    return audit_payload(
        contract,
        implementation_bytes.decode("utf-8"),
        implementation_bytes,
        json.loads(constraints_bytes.decode("utf-8")),
        constraints_bytes,
    )


def mutated_contract(**state_updates: bool) -> dict[str, Any]:
    """Test helper for fail-closed state-promotion regressions."""

    payload = copy.deepcopy(load_contract())
    payload["machine_state"].update(state_updates)
    return payload


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
