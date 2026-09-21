"""Fail-closed CR002 audit for partial composite oscillatory runtime identity.

The current A2 #970 semantic identity binds the A1 leading runtime identity plus
A2 adapter/differential identities.  The exact upstream sources inspected for
this scope do not bind a digest of the complete concrete oscillatory runtime
``to_payload()`` realization into that composite semantic identity.  This audit
prevents adapter/source-wrapper identity from being promoted to a complete
save/load velocity identity.

This is representation/provenance governance only.  It changes no velocity,
pressure, forcing, candidate coefficient, validation sample, or CR001 gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-partial-composite-osc-runtime-identity-scope-v1"
TASK = "CR002-KOKUNO-PARTIAL-COMPOSITE-OSC-RUNTIME-IDENTITY-098"
SCOPE_REL = Path("configs/kokuno_partial_composite_osc_runtime_identity_scope.json")
CONSTRAINTS_REL = Path("configs/constraints.json")
CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
BASE_HEAD = "f989bd99888e6263f9acb7008a932aa1db08d0aa"
A2_HEAD = "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
A2_SOURCE_BLOB = "532435705bd5db341e82371180da8c04c3dd6c14"
BATCH_ADAPTER_BLOB = "599baa190ec742d0e532e5517078534124e68d6b"
PUBLIC_Z_BLOB = "4a3abc1a11f7e054f651dd7bace8ada6d5a6b653"
A4_HEAD = "c45ee161014ec138deb14ed9085464855b2f88ee"


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def _git_blob_sha1(path: Path) -> str:
    raw = path.read_bytes()
    return hashlib.sha1(f"blob {len(raw)}\0".encode("ascii") + raw).hexdigest()


def _repo_root(start: Path | None = None) -> Path:
    here = (start or Path(__file__)).resolve()
    for parent in (here, *here.parents):
        if (parent / "configs" / "constraints.json").is_file():
            return parent
    raise RuntimeError("repository root with configs/constraints.json not found")


def load_scope(path: Path | None = None) -> dict[str, Any]:
    root = _repo_root() if path is None else None
    target = (root / SCOPE_REL) if path is None else Path(path)
    payload = json.loads(target.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("scope contract must be one JSON object")
    return payload


def _witness_errors(witness: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    if witness.get("classification") != "autonomous_design" or witness.get("not_source_evidence") is not True:
        errors.append("mechanics_witness_source_laundering")
        return errors
    adapter = witness.get("adapter_metadata")
    a = witness.get("runtime_a")
    b = witness.get("runtime_b")
    if not all(isinstance(v, Mapping) for v in (adapter, a, b)):
        return errors + ["mechanics_witness_missing_payload"]
    adapter_a = _sha256(adapter)
    adapter_b = _sha256(adapter)
    if adapter_a != adapter_b:
        errors.append("mechanics_witness_adapter_hash_not_equal")
    if a.get("velocity_at_probe") == b.get("velocity_at_probe"):
        errors.append("mechanics_witness_velocity_not_distinct")
    bound_a = _sha256({"adapter": adapter, "runtime": dict(a)})
    bound_b = _sha256({"adapter": adapter, "runtime": dict(b)})
    if bound_a == bound_b:
        errors.append("mechanics_witness_runtime_hash_not_distinct")
    return errors


def _constraints_errors(root: Path, binding: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    path = root / CONSTRAINTS_REL
    if _git_blob_sha1(path) != CONSTRAINTS_BLOB:
        errors.append("canonical_constraints_blob_drift")
        return errors
    c = json.loads(path.read_text(encoding="utf-8"))
    expected = {
        "constraints_blob_sha1": CONSTRAINTS_BLOB,
        "nu": c["nu"],
        "physical_domain": c["domain"]["physical"],
        "evaluation_box": c["domain"]["evaluation_box"],
        "support": c["domain"]["support"],
        "time_interval": c["domain"]["time_interval"],
        "forcing_mode": c["forcing"]["mode"],
        "reference_energy": c["nontriviality"]["reference_energy"],
        "reference_energy_abs_tolerance": c["nontriviality"]["reference_energy_abs_tolerance"],
        "validation_seed": c["validation"]["seed"],
        "held_out_points": c["validation"]["held_out_points"],
        "derivative_steps": c["validation"]["derivative_steps"],
        "quadrature_orders_per_axis": c["validation"]["quadrature_orders_per_axis"],
        "momentum_max": c["validation"]["thresholds"]["pde_residual_max"],
        "momentum_L2": c["validation"]["thresholds"]["pde_residual_L2"],
        "divergence_max": c["validation"]["thresholds"]["divergence_max"],
        "divergence_L2": c["validation"]["thresholds"]["divergence_L2"],
        "residual_defined_free_forcing_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    if dict(binding) != expected:
        errors.append("cr001_binding_drift")
    restriction = str(c["forcing"]["restriction"])
    if "No residual-dependent basis or pointwise free force" not in restriction:
        errors.append("canonical_free_force_firewall_missing")
    if "reject collapsed candidates" not in str(c["nontriviality"]["enforcement"]):
        errors.append("canonical_collapse_firewall_missing")
    if "changing thresholds requires a new experiment version" not in str(c["validation"]["failure_policy"]):
        errors.append("canonical_threshold_firewall_missing")
    return errors


def audit_scope(payload: Mapping[str, Any], *, root: Path | None = None) -> list[str]:
    errors: list[str] = []
    if payload.get("schema") != SCHEMA:
        errors.append("schema_drift")
    if payload.get("task_id") != TASK:
        errors.append("task_id_drift")

    base = payload.get("exact_base")
    if not isinstance(base, Mapping) or base.get("pr") != 973 or base.get("head") != BASE_HEAD:
        errors.append("exact_base_drift")

    upstream = payload.get("observed_upstream")
    if not isinstance(upstream, Mapping):
        return errors + ["observed_upstream_missing"]
    a2 = upstream.get("partial_composite")
    adapter = upstream.get("oscillatory_batch_adapter")
    public_z = upstream.get("public_z_field")
    a4 = upstream.get("a4_partial_composite_audit")
    if not isinstance(a2, Mapping) or a2.get("pr") != 970 or a2.get("head") != A2_HEAD or a2.get("source_blob_sha1") != A2_SOURCE_BLOB:
        errors.append("a2_identity_drift")
    if not isinstance(adapter, Mapping) or adapter.get("source_blob_sha1") != BATCH_ADAPTER_BLOB or adapter.get("public_z_velocity_source_blob_sha1") != PUBLIC_Z_BLOB:
        errors.append("osc_adapter_identity_drift")
    if not isinstance(public_z, Mapping) or public_z.get("source_blob_sha1") != PUBLIC_Z_BLOB:
        errors.append("public_z_identity_drift")
    if not isinstance(a4, Mapping) or a4.get("pr") != 972 or a4.get("head") != A4_HEAD:
        errors.append("a4_identity_drift")

    if isinstance(a2, Mapping):
        for key in (
            "semantic_payload_binds_agent1_runtime_semantic_sha256",
            "semantic_payload_binds_agent1_runtime_configuration_sha256",
            "semantic_payload_binds_agent2_batch_adapter_sha256",
            "semantic_payload_binds_agent2_differentials_sha256",
        ):
            if a2.get(key) is not True:
                errors.append(f"known_a2_binding_removed:{key}")
        if a2.get("semantic_payload_binds_full_oscillatory_runtime_payload_or_configuration_sha256") is not False:
            errors.append("full_osc_runtime_binding_laundered")
        if a2.get("composite_direct_save_load_available") is not False:
            errors.append("direct_composite_save_load_laundered")
    if isinstance(adapter, Mapping):
        if adapter.get("adapter_semantic_binds_support_and_time_interval") is not True:
            errors.append("adapter_support_time_binding_removed")
        if adapter.get("adapter_semantic_binds_public_z_wrapper_source_blob") is not True:
            errors.append("adapter_wrapper_blob_binding_removed")
        if adapter.get("adapter_semantic_binds_full_public_z_field_to_payload_digest") is not False:
            errors.append("adapter_full_runtime_binding_laundered")
    if isinstance(public_z, Mapping):
        if public_z.get("runtime_exposes_to_payload") is not True or public_z.get("to_payload_extends_inherited_runtime_realization") is not True:
            errors.append("public_z_runtime_payload_surface_removed")
        if public_z.get("runtime_payload_digest_is_bound_by_pr970_composite_semantic_identity") is not False:
            errors.append("public_z_runtime_digest_laundered")
    if isinstance(a4, Mapping):
        if a4.get("records_upstream_composite_direct_save_load_available") is not False:
            errors.append("a4_save_load_scope_laundered")
        if "not a serialization" not in str(a4.get("scientific_scope", "")):
            errors.append("a4_scope_laundered")

    state = payload.get("machine_locked_state")
    if not isinstance(state, Mapping):
        errors.append("machine_locked_state_missing")
    else:
        for key in (
            "partial_composite_callable_through_current_Xh",
            "partial_composite_semantic_sha_available",
            "oscillatory_adapter_semantic_sha_available",
        ):
            if state.get(key) is not True:
                errors.append(f"known_delivery_state_removed:{key}")
        for key in (
            "full_oscillatory_runtime_realization_identity_bound_into_composite_semantic_sha",
            "direct_composite_save_load_available",
            "composite_semantic_sha_sufficient_for_saved_velocity_identity",
            "velocity_export_ready_for_kokuno_route",
            "visual_correspondence_verified",
            "pde_validated",
            "paper_exact",
            "openai_field_identified",
        ):
            if state.get(key) is not False:
                errors.append(f"premature_promotion:{key}")

    provenance = payload.get("provenance_classes")
    if not isinstance(provenance, Mapping) or set(provenance) != {
        "user_requirement", "public_source_fact", "autonomous_design", "pending_unknown"
    }:
        errors.append("four_way_provenance_missing")
    else:
        public = " ".join(map(str, provenance["public_source_fact"]))
        if "hidden" in public.lower() or "paper-exact" in public.lower():
            errors.append("public_source_hidden_identity_laundering")
        pending = " ".join(map(str, provenance["pending_unknown"]))
        for token in ("runtime realization digest", "save/load", "PDE validation", "OpenAI-field identity"):
            if token.lower() not in pending.lower():
                errors.append(f"pending_boundary_missing:{token}")

    requirements = payload.get("future_promotion_requirements")
    text = " ".join(map(str, requirements)) if isinstance(requirements, list) else ""
    for token in ("complete concrete oscillatory runtime", "bind the digest", "save/load", "recompute"):
        if token.lower() not in text.lower():
            errors.append(f"future_identity_requirement_missing:{token}")

    witness = payload.get("mechanics_only_witness")
    if not isinstance(witness, Mapping):
        errors.append("mechanics_witness_missing")
    else:
        errors.extend(_witness_errors(witness))

    root = root or _repo_root()
    binding = payload.get("cr001_binding")
    if not isinstance(binding, Mapping):
        errors.append("cr001_binding_missing")
    else:
        errors.extend(_constraints_errors(root, binding))
    return errors


def assert_scope(payload: Mapping[str, Any], *, root: Path | None = None) -> None:
    errors = audit_scope(payload, root=root)
    if errors:
        raise RuntimeError("CR002 partial-composite runtime identity audit failed: " + ", ".join(errors))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", type=Path, default=None)
    args = parser.parse_args(argv)
    payload = load_scope(args.scope)
    assert_scope(payload)
    print("CR002 partial-composite oscillatory runtime identity scope: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
