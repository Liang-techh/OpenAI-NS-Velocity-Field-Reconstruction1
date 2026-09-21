"""Fail-closed CR002 audit for A2 #1052 source-sensitivity upstream inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


CONTRACT_RELATIVE = Path(
    "configs/kokuno_source_sensitivity_upstream_input_scope.json"
)
PARENT_MODULE_RELATIVE = Path(
    "src/openai_ns_reconstruction/kokuno_source_zero_data_amplitude_sensitivity.py"
)
CONSTRAINTS_RELATIVE = Path("configs/constraints.json")
PROJECT_STATUS_RELATIVE = Path("project_status.json")

EXPECTED_PARENT_HEAD = "d58cee2bf38bdca8da13a796215a079a569ad7d7"
EXPECTED_PARENT_BLOB = "eba3c00703bd5763b12e6f48a8c34eb918dec82e"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def load_contract(root: Path | None = None) -> dict[str, Any]:
    root = _repo_root() if root is None else Path(root)
    return json.loads((root / CONTRACT_RELATIVE).read_text(encoding="utf-8"))


def _append_if(violations: list[str], condition: bool, message: str) -> None:
    if condition:
        violations.append(message)


def audit_contract(
    contract: dict[str, Any],
    *,
    root: Path | None = None,
    check_live_files: bool = True,
) -> list[str]:
    """Return all provenance/representation violations; empty means coherent."""

    violations: list[str] = []
    parent = contract.get("exact_parent", {})
    truth = contract.get("parent_interface_truth", {})
    retirement = contract.get("input_retirement", {})
    boundary = contract.get("representation_boundary", {})
    forbidden = contract.get("forbidden_promotions", {})
    cr001 = contract.get("canonical_cr001", {})
    eq45 = contract.get("canonical_eq45_state_must_remain_independent", {})
    provenance = contract.get("four_way_provenance", {})

    _append_if(
        violations,
        contract.get("contract_id")
        != "CR002-KOKUNO-SOURCE-SENSITIVITY-UPSTREAM-INPUT-SCOPE-112",
        "contract id drift",
    )
    _append_if(violations, parent.get("pr") != 1052, "parent PR drift")
    _append_if(
        violations, parent.get("head") != EXPECTED_PARENT_HEAD, "parent head drift"
    )
    _append_if(
        violations,
        parent.get("module_git_blob") != EXPECTED_PARENT_BLOB,
        "parent module blob declaration drift",
    )

    required_provenance = {
        "user_requirement",
        "public_source_fact",
        "autonomous_design",
        "pending_unknown",
    }
    _append_if(
        violations,
        set(provenance) != required_provenance,
        "four-way provenance classes must remain exact and separate",
    )
    for key in sorted(required_provenance):
        _append_if(
            violations,
            not isinstance(provenance.get(key), list) or not provenance.get(key),
            f"provenance class {key} must remain nonempty",
        )

    required_true = {
        "source_zero_data_ivp_numerical_solver_materialized",
        "source_amplitude_directional_sensitivity_solver_materialized",
        "source_amplitude_state_t_m_generated_internally",
        "source_amplitude_dr_t_m_generated_internally",
        "source_amplitude_dz_t_m_generated_internally",
        "caller_supplies_corrected_background",
        "caller_supplies_mode_forcing_value",
        "caller_supplies_full_normalized_dr_forcing",
        "caller_supplies_full_normalized_dz_forcing",
        "finest_directional_jet_feeds_complete_curl_path",
    }
    required_false = {
        "caller_supplies_mode_state_t_m",
        "caller_supplies_amplitude_dr_t_m",
        "caller_supplies_amplitude_dz_t_m",
        "caller_can_tune_integrator_resolution",
        "source_exact_duhamel_frame_B_materialized",
        "source_exact_duhamel_propagator_Vm_materialized",
        "deterministic_corrected_background_provider_materialized",
        "deterministic_source_forcing_provider_materialized",
        "source_amplitude_mode_provider_complete",
        "source_support_cutoff_runtime_materialized",
        "self_contained_velocity_xyzt_provider",
        "source_to_runtime_parameter_map_complete",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    }
    for key in sorted(required_true):
        _append_if(violations, truth.get(key) is not True, f"{key} must remain true")
    for key in sorted(required_false):
        _append_if(violations, truth.get(key) is not False, f"{key} must remain false")

    _append_if(
        violations,
        retirement.get("retired_caller_inputs") != ["t_m", "D_r t_m", "D_z t_m"],
        "#1052 caller-input retirement set drift",
    )
    expected_upstream = [
        "corrected_background_jet",
        "forcing_m(v).value",
        "forcing_m(v).D_r_value_full_normalized",
        "forcing_m(v).D_z_value_full_normalized",
    ]
    _append_if(
        violations,
        retirement.get("still_caller_supplied_output_determining_inputs")
        != expected_upstream,
        "remaining upstream output-determining inputs drift",
    )
    for key in (
        "generated_directional_sensitivities_are_conditional_on_upstream_forcing_jets",
        "generated_directional_sensitivities_are_conditional_on_corrected_background",
        "retiring_amplitude_directional_inputs_does_not_make_source_amplitude_provider_complete",
        "retiring_amplitude_directional_inputs_does_not_make_source_route_self_contained",
    ):
        _append_if(
            violations,
            retirement.get(key) is not True,
            f"input-retirement boundary {key} must remain true",
        )

    for key in (
        "generated_amplitude_directional_jet_is_not_upstream_input_free",
        "internal_sensitivity_integration_is_not_deterministic_source_forcing_materialization",
        "fixed_rk4_refinement_evidence_is_not_exact_source_duhamel_propagator_evidence",
        "complete_curl_consumption_is_not_unified_cartesian_velocity_delivery",
        "source_route_incompleteness_does_not_downgrade_canonical_eq45_delivery",
    ):
        _append_if(
            violations,
            boundary.get(key) is not True,
            f"representation boundary {key} must remain true",
        )

    for key, value in forbidden.items():
        _append_if(
            violations, value is not False, f"forbidden promotion {key} must remain false"
        )

    expected_cr001 = {
        "constraints_git_blob": EXPECTED_CONSTRAINTS_BLOB,
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
        "momentum_L2": 0.001,
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
        "residual_defined_free_force_forbidden": True,
        "candidate_collapse_forbidden": True,
        "post_hoc_threshold_relaxation_forbidden": True,
    }
    for key, expected in expected_cr001.items():
        _append_if(
            violations, cr001.get(key) != expected, f"canonical CR001 field {key} drift"
        )

    expected_eq45 = {
        "candidate_family": "eq45_supported_velocity_candidate_v1",
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
    }
    for key, expected in expected_eq45.items():
        _append_if(
            violations, eq45.get(key) != expected, f"canonical Eq45 state {key} drift"
        )

    if not check_live_files:
        return violations

    root = _repo_root() if root is None else Path(root)
    parent_path = root / PARENT_MODULE_RELATIVE
    constraints_path = root / CONSTRAINTS_RELATIVE
    status_path = root / PROJECT_STATUS_RELATIVE
    for path in (parent_path, constraints_path, status_path):
        if not path.exists():
            violations.append(f"required live file missing: {path.relative_to(root)}")
    if any(item.startswith("required live file missing") for item in violations):
        return violations

    _append_if(
        violations,
        _git_blob_sha1(parent_path.read_bytes()) != EXPECTED_PARENT_BLOB,
        "exact #1052 sensitivity module blob mismatch",
    )
    _append_if(
        violations,
        _git_blob_sha1(constraints_path.read_bytes()) != EXPECTED_CONSTRAINTS_BLOB,
        "canonical configs/constraints.json blob mismatch",
    )

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
    domain = constraints.get("domain", {})
    _append_if(violations, constraints.get("nu") != 0.01, "live nu drift")
    _append_if(violations, domain.get("physical") != "R^3", "live physical domain drift")
    _append_if(
        violations,
        domain.get("evaluation_box") != [[-2, 2], [-2, 2], [-2, 2]],
        "live evaluation box drift",
    )
    _append_if(
        violations,
        domain.get("support") != "r < 2 and abs(z) < 2",
        "live support drift",
    )
    _append_if(
        violations,
        domain.get("time_interval") != [0.25, 0.75],
        "live time interval drift",
    )
    _append_if(
        violations,
        constraints.get("forcing", {}).get("mode") != "restricted_two_parameter_family",
        "live forcing mode drift",
    )
    _append_if(
        violations,
        validation.get("seed") != 914027 or validation.get("held_out_points") != 4096,
        "live held-out validation contract drift",
    )
    _append_if(
        violations,
        validation.get("derivative_steps") != [0.02, 0.01, 0.005],
        "live derivative ladder drift",
    )
    _append_if(
        violations,
        validation.get("quadrature_orders_per_axis") != [24, 48, 96],
        "live quadrature ladder drift",
    )
    _append_if(
        violations,
        thresholds.get("pde_residual_max") != 0.001
        or thresholds.get("pde_residual_L2") != 0.001,
        "live momentum threshold drift",
    )
    _append_if(
        violations,
        thresholds.get("divergence_max") != 1e-5
        or thresholds.get("divergence_L2") != 1e-5,
        "live divergence threshold drift",
    )

    status = json.loads(status_path.read_text(encoding="utf-8"))
    states = status.get("states", {})
    _append_if(
        violations,
        status.get("candidate_family") != "eq45_supported_velocity_candidate_v1",
        "live canonical candidate family drift",
    )
    for key, expected in expected_eq45.items():
        if key == "candidate_family":
            continue
        _append_if(
            violations,
            states.get(key) != expected,
            f"live canonical Eq45 state {key} drift",
        )

    return violations


def mechanics_only_upstream_forcing_jet_witness() -> dict[str, Any]:
    """Toy zero-data sensitivity: generated D_r t still depends on upstream D_r f."""

    pulse_v = 0.4
    forcing_value = 1.25
    dr_forcing_a = 0.25
    dr_forcing_b = 2.0

    state_a = pulse_v * forcing_value
    state_b = pulse_v * forcing_value
    dr_state_a = pulse_v * dr_forcing_a
    dr_state_b = pulse_v * dr_forcing_b

    if state_a != state_b:
        raise RuntimeError("base-state witness must share the same forcing value")
    if dr_state_a == dr_state_b:
        raise RuntimeError("directional sensitivity witness unexpectedly collapsed")

    return {
        "classification": "autonomous mechanics only; not Kokuno/OpenAI candidate data",
        "toy_system": "dt/dv=f, D_r t derivative=(D_r f), common zero data",
        "pulse_v": pulse_v,
        "shared_generated_state": state_a,
        "caller_dr_forcing_a": dr_forcing_a,
        "caller_dr_forcing_b": dr_forcing_b,
        "generated_dr_state_a": dr_state_a,
        "generated_dr_state_b": dr_state_b,
        "generated_sensitivity_remains_conditional_on_upstream_forcing_jet": True,
    }


def assert_clean(root: Path | None = None) -> None:
    contract = load_contract(root)
    violations = audit_contract(contract, root=root, check_live_files=True)
    if violations:
        raise RuntimeError(
            "CR002 source-sensitivity upstream-input scope violations:\n- "
            + "\n- ".join(violations)
        )


if __name__ == "__main__":
    assert_clean()
    print(json.dumps(mechanics_only_upstream_forcing_jet_witness(), indent=2))
