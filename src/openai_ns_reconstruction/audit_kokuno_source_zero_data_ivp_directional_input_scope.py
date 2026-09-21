"""Fail-closed CR002 audit for A2 #1044 zero-data IVP input retirement."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


CONTRACT_RELATIVE = Path(
    "configs/kokuno_source_zero_data_ivp_directional_input_scope.json"
)
PARENT_MODULE_RELATIVE = Path(
    "src/openai_ns_reconstruction/kokuno_source_zero_data_amplitude_ivp.py"
)
CONSTRAINTS_RELATIVE = Path("configs/constraints.json")
PROJECT_STATUS_RELATIVE = Path("project_status.json")

EXPECTED_PARENT_HEAD = "b3ac194b410073820d8160a470824034e08e9b5a"
EXPECTED_PARENT_BLOB = "1e13e81bea6c586938c91c610e0812f7283f131d"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_blob_sha1(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


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
    """Return all representation/provenance violations; empty means coherent."""

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
        != "CR002-KOKUNO-ZERO-DATA-IVP-DIRECTIONAL-INPUT-SCOPE-111",
        "contract id drift",
    )
    _append_if(violations, parent.get("pr") != 1044, "parent PR drift")
    _append_if(
        violations,
        parent.get("head") != EXPECTED_PARENT_HEAD,
        "parent head drift",
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
        "source_projected_amplitude_rhs_materialized",
        "source_zero_data_ivp_numerical_solver_materialized",
        "caller_supplies_mode_forcing_f_m",
        "caller_supplies_corrected_background",
        "source_phase_jet_regenerated_analytically_along_v",
    }
    required_false = {
        "caller_supplies_mode_state_t_m",
        "caller_can_tune_integrator_resolution",
        "source_exact_duhamel_frame_B_materialized",
        "source_exact_duhamel_propagator_Vm_materialized",
        "source_amplitude_directional_jet_materialized",
        "source_amplitude_mode_provider_complete",
        "actual_corrected_background_provider_materialized",
        "source_support_cutoff_runtime_materialized",
        "self_contained_velocity_xyzt_provider",
        "source_to_runtime_parameter_map_complete",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "paper_exact",
    }
    for key in sorted(required_true):
        _append_if(violations, truth.get(key) is not True, f"{key} must remain true")
    for key in sorted(required_false):
        _append_if(violations, truth.get(key) is not False, f"{key} must remain false")

    _append_if(
        violations,
        retirement.get("retired_caller_inputs") != ["t_m"],
        "only t_m is retired from caller supply by #1044",
    )
    _append_if(
        violations,
        retirement.get("still_caller_supplied_output_determining_inputs")
        != ["corrected_background_jet", "forcing_m(v)"],
        "remaining caller-supplied output-determining inputs drift",
    )
    _append_if(
        violations,
        retirement.get("still_unmaterialized_amplitude_directional_outputs")
        != ["D_r t_m", "D_z t_m"],
        "unmaterialized amplitude directional outputs drift",
    )
    _append_if(
        violations,
        retirement.get("retiring_t_m_does_not_retire_directional_jets") is not True,
        "t_m retirement must not imply directional-jet retirement",
    )
    _append_if(
        violations,
        retirement.get("retiring_t_m_does_not_make_source_route_self_contained")
        is not True,
        "t_m retirement must not imply a self-contained source route",
    )

    required_boundary_true = {
        "generated_zero_data_mode_state_is_not_complete_directional_jet",
        "fixed_rk4_refinement_evidence_is_not_exact_source_duhamel_propagator_evidence",
        "caller_free_t_m_is_not_caller_free_source_realization",
        "source_ivp_provider_is_not_unified_cartesian_velocity_provider",
        "source_route_incompleteness_does_not_downgrade_canonical_eq45_delivery",
    }
    for key in sorted(required_boundary_true):
        _append_if(
            violations,
            boundary.get(key) is not True,
            f"representation boundary {key} must remain true",
        )

    for key, value in forbidden.items():
        _append_if(
            violations,
            value is not False,
            f"forbidden promotion {key} must remain false",
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
            violations,
            cr001.get(key) != expected,
            f"canonical CR001 field {key} drift",
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
            violations,
            eq45.get(key) != expected,
            f"canonical Eq45 state {key} drift",
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
    if any("required live file missing" in item for item in violations):
        return violations

    _append_if(
        violations,
        _git_blob_sha1(parent_path.read_bytes()) != EXPECTED_PARENT_BLOB,
        "exact #1044 parent module blob mismatch",
    )
    _append_if(
        violations,
        _git_blob_sha1(constraints_path.read_bytes()) != EXPECTED_CONSTRAINTS_BLOB,
        "canonical configs/constraints.json blob mismatch",
    )

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    domain = constraints.get("domain", {})
    validation = constraints.get("validation", {})
    thresholds = validation.get("thresholds", {})
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
        constraints.get("forcing", {}).get("mode")
        != "restricted_two_parameter_family",
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


def mechanics_only_partial_input_retirement_witness() -> dict[str, Any]:
    """Toy logic: zero-data state generation does not close all source inputs/jets."""

    pulse_v = 0.5
    state_from_forcing_one = pulse_v * 1.0
    state_from_forcing_two = pulse_v * 2.0
    if state_from_forcing_one == state_from_forcing_two:
        raise RuntimeError("forcing-dependence witness unexpectedly collapsed")

    r0 = 1.0
    z0 = -0.25
    state0 = 0.7

    def amp_a(r: float, z: float) -> float:
        return state0 + 1.0 * (r - r0) + 2.0 * (z - z0)

    def amp_b(r: float, z: float) -> float:
        return state0 + 3.0 * (r - r0) - 1.0 * (z - z0)

    if amp_a(r0, z0) != amp_b(r0, z0):
        raise RuntimeError("directional-jet witness must share the same point value")

    h = 1e-6
    dr_a = (amp_a(r0 + h, z0) - amp_a(r0 - h, z0)) / (2.0 * h)
    dr_b = (amp_b(r0 + h, z0) - amp_b(r0 - h, z0)) / (2.0 * h)
    dz_a = (amp_a(r0, z0 + h) - amp_a(r0, z0 - h)) / (2.0 * h)
    dz_b = (amp_b(r0, z0 + h) - amp_b(r0, z0 - h)) / (2.0 * h)
    if dr_a == dr_b or dz_a == dz_b:
        raise RuntimeError("directional-jet witness unexpectedly collapsed")

    return {
        "classification": "autonomous mechanics only; not Kokuno/OpenAI candidate data",
        "zero_data_toy_ivp": "dt/dv=f(v), t(0)=0",
        "pulse_v": pulse_v,
        "state_for_constant_forcing_1": state_from_forcing_one,
        "state_for_constant_forcing_2": state_from_forcing_two,
        "caller_forcing_remains_output_determining": True,
        "shared_point_state": state0,
        "directional_jet_a": [dr_a, dz_a],
        "directional_jet_b": [dr_b, dz_b],
        "same_point_state_does_not_determine_directional_jet": True,
    }


def assert_clean(root: Path | None = None) -> None:
    contract = load_contract(root)
    violations = audit_contract(contract, root=root, check_live_files=True)
    if violations:
        raise RuntimeError(
            "CR002 zero-data IVP directional-input scope violations:\n- "
            + "\n- ".join(violations)
        )


if __name__ == "__main__":
    assert_clean()
    print(json.dumps(mechanics_only_partial_input_retirement_witness(), indent=2))
