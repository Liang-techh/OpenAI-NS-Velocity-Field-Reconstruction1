"""Fail-closed CR002 audit for Kokuno source complete-curl delivery input scope.

This module governs a representation boundary only.  A2 #1018 makes the
corrected-source normalized complete-curl algebra executable, but its public
reference still consumes source directional jets supplied by the caller and
fails closed at the axis.  Those facts must not be promoted into a self-
contained Cartesian ``velocity(x,y,z,t)`` delivery claim.

No scientific threshold, field coefficient, pressure, forcing or source
formula is changed here.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any, Mapping

from openai_ns_reconstruction import kokuno_source_normalized_complete_curl_reference as source_ref


SCHEMA = "cr002-kokuno-source-complete-curl-delivery-input-scope-v1"
TASK_ID = "CR002-KOKUNO-SOURCE-CURL-DELIVERY-INPUT-SCOPE-106"
EXPECTED_BASE_HEAD = "ecc2be96eee441308cc0216bf187b43bbd7fd3a0"
EXPECTED_SOURCE_MODULE_BLOB = "42809fc1b10937b4246c65121ac1dbfdc2986822"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
CONTRACT_RELATIVE_PATH = Path("configs/kokuno_source_complete_curl_delivery_input_scope.json")
SOURCE_MODULE_RELATIVE_PATH = Path(
    "src/openai_ns_reconstruction/kokuno_source_normalized_complete_curl_reference.py"
)


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def load_contract(path: str | Path | None = None) -> dict[str, Any]:
    target = Path(path) if path is not None else _repository_root() / CONTRACT_RELATIVE_PATH
    return json.loads(target.read_text(encoding="utf-8"))


def _require(errors: list[str], condition: bool, message: str) -> None:
    if not condition:
        errors.append(message)


def _audit_contract_shape(contract: Mapping[str, Any], errors: list[str]) -> None:
    _require(errors, contract.get("schema") == SCHEMA, "scope schema drift")
    _require(errors, contract.get("task_id") == TASK_ID, "task id drift")

    base = contract.get("exact_base", {})
    _require(errors, base.get("pr") == 1018, "exact base PR must remain #1018")
    _require(errors, base.get("head") == EXPECTED_BASE_HEAD, "exact base head drift")
    _require(
        errors,
        base.get("source_module_blob") == EXPECTED_SOURCE_MODULE_BLOB,
        "exact #1018 source-reference blob drift",
    )

    provenance = contract.get("four_way_provenance", {})
    _require(
        errors,
        set(provenance)
        == {
            "user_requirements",
            "public_source_facts",
            "autonomous_repository_facts",
            "pending_or_unknown",
        },
        "four-way provenance classes must remain separate and complete",
    )
    for key in (
        "user_requirements",
        "public_source_facts",
        "autonomous_repository_facts",
        "pending_or_unknown",
    ):
        value = provenance.get(key)
        _require(errors, isinstance(value, list) and bool(value), f"{key} must be nonempty")

    upstream = contract.get("upstream_source_reference", {})
    expected_upstream = {
        "source_repository": "KokunoYumeto/yang-mills-interacting-workbench",
        "source_commit": "143f6773feb424ad9ed3a8d116653200f20346b7",
        "source_workbench_blob": "205a99807302e21a51c5eaf223390c0dfc42bcd0",
        "source_normalized_curl_formula_executable": True,
        "source_complete_harmonic_amplitude_algebra_executable": True,
        "source_annular_domain_requires_positive_radius": True,
        "source_directional_derivatives_supplied_by_caller": True,
        "current_runtime_phase_source_exact": False,
        "current_runtime_support_source_exact": False,
        "current_runtime_complete_curl_source_equivalence_verified": False,
        "source_to_runtime_parameter_map_complete": False,
    }
    for key, expected in expected_upstream.items():
        _require(errors, upstream.get(key) == expected, f"upstream source scope drift: {key}")

    boundary = contract.get("delivery_input_boundary", {})
    _require(
        errors,
        boundary.get("project_delivery_signature") == "velocity(x,y,z,t)->[u,v,w]",
        "project delivery signature drift",
    )
    required_inputs = {
        "radius",
        "normalized potential directional jet",
        "n_phi",
        "t_m",
        "D_r C_m",
        "D_z C_m",
        "k",
        "m",
    }
    _require(
        errors,
        set(boundary.get("source_reference_primary_inputs", [])) == required_inputs,
        "source-reference primary input surface drift",
    )
    for key in (
        "source_reference_self_contained_cartesian_velocity_provider",
        "source_reference_accepts_x_y_z_t_only",
        "source_reference_defines_global_axis_extension",
        "source_reference_establishes_axis_safe_global_cartesian_runtime",
        "source_reference_establishes_project_domain_totality",
        "source_reference_establishes_identity_preserving_velocity_save_load",
        "source_reference_velocity_export_ready",
    ):
        _require(errors, boundary.get(key) is False, f"forbidden delivery promotion: {key}")
    _require(
        errors,
        boundary.get("source_formula_executable_can_be_true_while_velocity_delivery_remains_false")
        is True,
        "formula-execution versus delivery independence must remain explicit",
    )

    locked = contract.get("machine_locked_nonimplications", {})
    for key in (
        "executable_source_formula_implies_self_contained_velocity_provider",
        "caller_supplied_directional_jets_imply_source_background_realized",
        "manufactured_curl_or_divergence_identity_implies_runtime_source_equivalence",
        "annular_fail_closed_reference_implies_project_velocity_cannot_be_axis_safe",
        "source_formula_reference_implies_visual_correspondence",
        "source_formula_reference_implies_pde_validation",
        "source_formula_reference_implies_paper_exact",
        "source_formula_reference_implies_openai_field_identity",
        "green_or_queued_ci_implies_scientific_promotion",
    ):
        _require(errors, locked.get(key) is False, f"nonimplication promoted: {key}")

    witness = contract.get("scope_logic_witness", {})
    _require(
        errors,
        witness.get("classification") == "autonomous_mechanics_only",
        "mechanics witness must remain autonomous-only",
    )
    _require(
        errors,
        witness.get("not_a_public_source_parameter_claim") is True,
        "mechanics witness must not be public-source parameter evidence",
    )
    _require(
        errors,
        witness.get("not_candidate_numerical_evidence") is True,
        "mechanics witness must not be candidate numerical evidence",
    )

    cr001 = contract.get("cr001_snapshot", {})
    expected_cr001 = {
        "constraints_blob_sha1": EXPECTED_CONSTRAINTS_BLOB,
        "nu": 0.01,
        "physical_domain": "R^3",
        "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
        "support": "r < 2 and abs(z) < 2",
        "time_interval": [0.25, 0.75],
        "forcing_mode": "restricted_two_parameter_family",
        "reference_energy": 1.0,
        "reference_energy_abs_tolerance": 0.001,
        "validation_seed": 914027,
        "held_out_points": 4096,
        "derivative_steps": [0.02, 0.01, 0.005],
        "quadrature_orders_per_axis": [24, 48, 96],
        "momentum_max": 0.001,
        "momentum_l2": 0.001,
        "divergence_max": 1e-5,
        "divergence_l2": 1e-5,
        "residual_defined_free_forcing_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    for key, expected in expected_cr001.items():
        _require(errors, cr001.get(key) == expected, f"CR001 snapshot drift: {key}")

    canonical = contract.get("canonical_delivery_independence", {})
    expected_canonical = {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_api": "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, expected in expected_canonical.items():
        _require(errors, canonical.get(key) == expected, f"canonical Eq45 state drift: {key}")

    kokuno = contract.get("kokuno_source_reference_state", {})
    for key in (
        "self_contained_cartesian_velocity_provider",
        "velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(errors, kokuno.get(key) is False, f"Kokuno source-reference promotion: {key}")

    promotions = contract.get("forbidden_promotions", [])
    required_fragments = (
        "self-contained velocity(x,y,z,t)",
        "source background/phase realization complete",
        "global axis-safe Cartesian realization",
        "current autonomous runtime is source-equivalent",
        "visual correspondence",
        "PDE validation or exact OpenAI field",
        "downgrade canonical Eq45 velocity_export_ready",
        "unrestricted pointwise forcing",
        "post-hoc threshold relaxation",
    )
    joined = "\n".join(str(value) for value in promotions)
    for fragment in required_fragments:
        _require(errors, fragment in joined, f"missing forbidden-promotion firewall: {fragment}")


def _audit_live_source_surface(errors: list[str]) -> None:
    upstream = source_ref.source_reference_contract()
    expected = {
        "source_normalized_curl_formula_executable": True,
        "source_complete_harmonic_amplitude_algebra_executable": True,
        "source_annular_domain_requires_positive_radius": True,
        "source_directional_derivatives_supplied_by_caller": True,
        "current_runtime_phase_source_exact": False,
        "current_runtime_support_source_exact": False,
        "current_runtime_complete_curl_source_equivalence_verified": False,
        "source_to_runtime_parameter_map_complete": False,
        "paper_exact": False,
        "matched_pressure_materialized": False,
        "restricted_forcing_materialized": False,
        "complete_ns_residual_assessed": False,
        "pde_validated": False,
    }
    for key, value in expected.items():
        _require(errors, upstream.get(key) == value, f"live #1018 source contract drift: {key}")

    complete_signature = inspect.signature(
        source_ref.complete_harmonic_amplitude_from_coefficient_jet
    )
    complete_parameters = set(complete_signature.parameters)
    _require(
        errors,
        {"radius", "n_phi", "t_m", "coefficient_dr", "coefficient_dz", "k", "m"}
        <= complete_parameters,
        "live complete-harmonic reference no longer exposes the governed caller-supplied jets",
    )
    _require(
        errors,
        not hasattr(source_ref, "velocity"),
        "#1018 source-reference module unexpectedly exposes a velocity provider",
    )


def _audit_live_repository_files(root: Path, errors: list[str]) -> None:
    constraints_path = root / "configs/constraints.json"
    project_status_path = root / "project_status.json"
    source_module_path = root / SOURCE_MODULE_RELATIVE_PATH

    _require(errors, constraints_path.is_file(), "canonical configs/constraints.json missing")
    _require(errors, project_status_path.is_file(), "project_status.json missing")
    _require(errors, source_module_path.is_file(), "exact #1018 source-reference module missing")
    if errors:
        return

    _require(
        errors,
        _git_blob_sha1(constraints_path) == EXPECTED_CONSTRAINTS_BLOB,
        "canonical constraints blob changed",
    )
    _require(
        errors,
        _git_blob_sha1(source_module_path) == EXPECTED_SOURCE_MODULE_BLOB,
        "exact #1018 source-reference module bytes changed",
    )

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    _require(errors, constraints.get("nu") == 0.01, "live nu drift")
    domain = constraints.get("domain", {})
    _require(errors, domain.get("physical") == "R^3", "live physical domain drift")
    _require(
        errors,
        domain.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]],
        "live evaluation box drift",
    )
    _require(errors, domain.get("support") == "r < 2 and abs(z) < 2", "live support drift")
    _require(errors, domain.get("time_interval") == [0.25, 0.75], "live time interval drift")

    forcing = constraints.get("forcing", {})
    _require(
        errors,
        forcing.get("mode") == "restricted_two_parameter_family",
        "live forcing family drift",
    )
    restriction = str(forcing.get("restriction", ""))
    _require(
        errors,
        "No residual-dependent basis or pointwise free force" in restriction,
        "live free-force prohibition missing",
    )
    _require(
        errors,
        set(forcing.get("parameters", {})) == {"a", "c"},
        "live forcing parameter set must remain exactly a,c",
    )

    nontriviality = constraints.get("nontriviality", {})
    _require(errors, nontriviality.get("reference_time") == 0.25, "live energy time drift")
    _require(errors, nontriviality.get("reference_energy") == 1.0, "live energy target drift")
    _require(
        errors,
        nontriviality.get("reference_energy_abs_tolerance") == 0.001,
        "live energy tolerance drift",
    )

    validation = constraints.get("validation", {})
    _require(errors, validation.get("seed") == 914027, "live validation seed drift")
    _require(errors, validation.get("held_out_points") == 4096, "live held-out count drift")
    _require(
        errors,
        validation.get("derivative_steps") == [0.02, 0.01, 0.005],
        "live derivative ladder drift",
    )
    _require(
        errors,
        validation.get("quadrature_orders_per_axis") == [24, 48, 96],
        "live quadrature ladder drift",
    )
    thresholds = validation.get("thresholds", {})
    _require(errors, thresholds.get("pde_residual_max") == 0.001, "live momentum max drift")
    _require(errors, thresholds.get("pde_residual_L2") == 0.001, "live momentum L2 drift")
    _require(errors, thresholds.get("divergence_max") == 1e-5, "live divergence max drift")
    _require(errors, thresholds.get("divergence_L2") == 1e-5, "live divergence L2 drift")

    project_status = json.loads(project_status_path.read_text(encoding="utf-8"))
    _require(
        errors,
        project_status.get("candidate_family") == "eq45_supported_velocity_candidate_v1",
        "canonical candidate family drift",
    )
    _require(
        errors,
        project_status.get("velocity_api")
        == "openai_ns_reconstruction.eq45_supported_delivery:velocity",
        "canonical velocity API drift",
    )
    states = project_status.get("states", {})
    expected_states = {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, value in expected_states.items():
        _require(errors, states.get(key) == value, f"live canonical state drift: {key}")


def audit_scope(
    contract: Mapping[str, Any] | None = None,
    *,
    repository_root: str | Path | None = None,
    verify_live_files: bool = True,
) -> list[str]:
    """Return fail-closed scope violations; an empty list means this audit passes."""

    payload: Mapping[str, Any] = load_contract() if contract is None else deepcopy(dict(contract))
    errors: list[str] = []
    _audit_contract_shape(payload, errors)
    _audit_live_source_surface(errors)
    if verify_live_files:
        root = Path(repository_root) if repository_root is not None else _repository_root()
        _audit_live_repository_files(root, errors)
    return errors


def assert_scope(
    contract: Mapping[str, Any] | None = None,
    *,
    repository_root: str | Path | None = None,
    verify_live_files: bool = True,
) -> None:
    errors = audit_scope(
        contract,
        repository_root=repository_root,
        verify_live_files=verify_live_files,
    )
    if errors:
        raise RuntimeError("CR002 source complete-curl delivery-input scope failed:\n- " + "\n- ".join(errors))


def main() -> int:
    errors = audit_scope()
    if errors:
        print(json.dumps({"ok": False, "errors": errors}, indent=2))
        return 1
    print(
        json.dumps(
            {
                "ok": True,
                "task_id": TASK_ID,
                "exact_base_head": EXPECTED_BASE_HEAD,
                "source_reference_self_contained_cartesian_velocity_provider": False,
                "canonical_eq45_velocity_export_ready": True,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
