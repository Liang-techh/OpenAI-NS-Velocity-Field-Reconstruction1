"""Fail-closed CR002 audit for the remaining source-amplitude delivery inputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


CONTRACT_RELATIVE = Path("configs/kokuno_source_amplitude_input_scope.json")
PARENT_MODULE_RELATIVE = Path(
    "src/openai_ns_reconstruction/kokuno_source_phase_coefficient_jet.py"
)
CONSTRAINTS_RELATIVE = Path("configs/constraints.json")
PROJECT_STATUS_RELATIVE = Path("project_status.json")

EXPECTED_PARENT_HEAD = "0f4006bfd159281d96a7dc90a3ab50f0566f13e1"
EXPECTED_PARENT_BLOB = "d3c9684d47cec5dfabf713db01ffab45ea9461c4"
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
    """Return every scope/provenance violation; empty means the contract is coherent."""

    violations: list[str] = []
    parent = contract.get("exact_parent", {})
    truth = contract.get("parent_interface_truth", {})
    remaining = contract.get("remaining_caller_inputs", {})
    forbidden = contract.get("forbidden_promotions", {})
    cr001 = contract.get("canonical_cr001", {})
    eq45 = contract.get("canonical_eq45_state_must_remain_independent", {})

    _append_if(
        violations,
        contract.get("contract_id") != "CR002-KOKUNO-SOURCE-AMPLITUDE-INPUT-SCOPE-107",
        "contract id drift",
    )
    _append_if(violations, parent.get("pr") != 1021, "parent PR drift")
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

    required_true = {
        "source_phase_vector_jet_executable",
        "source_coefficient_directional_jets_derived_analytically",
        "source_complete_curl_consumes_derived_coefficient_jets",
        "caller_supplies_background_phase_derivative_jet",
        "caller_supplies_amplitude_directional_jet",
    }
    required_false = {
        "caller_supplies_coefficient_directional_jets",
        "source_amplitude_ode_materialized",
        "source_support_cutoff_runtime_materialized",
        "self_contained_velocity_xyzt_provider",
        "source_to_runtime_parameter_map_complete",
        "current_runtime_complete_curl_source_equivalence_verified",
        "paper_exact",
        "complete_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
    }
    for key in sorted(required_true):
        _append_if(violations, truth.get(key) is not True, f"{key} must remain true")
    for key in sorted(required_false):
        _append_if(violations, truth.get(key) is not False, f"{key} must remain false")

    background = remaining.get("background_phase_derivative_jet", {})
    amplitude = remaining.get("amplitude_directional_jet", {})
    _append_if(
        violations,
        background.get("required") is not True,
        "background phase derivative jet must remain an explicit caller input",
    )
    _append_if(
        violations,
        amplitude.get("required") is not True,
        "amplitude directional jet must remain an explicit caller input",
    )
    _append_if(
        violations,
        background.get("fields")
        != [
            "f",
            "g",
            "f_r",
            "g_r",
            "f_z",
            "g_z",
            "f_rr",
            "g_rr",
            "f_rz",
            "g_rz",
            "f_zz",
            "g_zz",
        ],
        "background derivative-jet field set drift",
    )
    _append_if(
        violations,
        amplitude.get("fields") != ["t_m", "dr_t_m", "dz_t_m"],
        "amplitude directional-jet field set drift",
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
    if violations and any("missing" in item for item in violations):
        return violations

    _append_if(
        violations,
        _git_blob_sha1(parent_path.read_bytes()) != EXPECTED_PARENT_BLOB,
        "exact #1021 parent module blob mismatch",
    )
    _append_if(
        violations,
        _git_blob_sha1(constraints_path.read_bytes()) != EXPECTED_CONSTRAINTS_BLOB,
        "canonical configs/constraints.json blob mismatch",
    )

    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    _append_if(violations, constraints.get("nu") != 0.01, "live nu drift")
    _append_if(
        violations,
        constraints.get("domain", {}).get("physical") != "R^3",
        "live physical domain drift",
    )
    _append_if(
        violations,
        constraints.get("domain", {}).get("evaluation_box")
        != [[-2, 2], [-2, 2], [-2, 2]],
        "live evaluation box drift",
    )
    _append_if(
        violations,
        constraints.get("domain", {}).get("support") != "r < 2 and abs(z) < 2",
        "live support drift",
    )
    _append_if(
        violations,
        constraints.get("domain", {}).get("time_interval") != [0.25, 0.75],
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
        constraints.get("validation", {}).get("seed") != 914027,
        "live validation seed drift",
    )
    _append_if(
        violations,
        constraints.get("validation", {}).get("held_out_points") != 4096,
        "live held-out count drift",
    )
    thresholds = constraints.get("validation", {}).get("thresholds", {})
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


def mechanics_only_amplitude_dependence_witness() -> dict[str, Any]:
    """Toy proof that caller-supplied transverse amplitude remains output-determining."""

    n = (1.0, 0.0, 1.0)
    t_a = (0.0, 1.0, 0.0)
    t_b = (1.0, 0.0, -1.0)

    def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
        return sum(x * y for x, y in zip(a, b))

    def cross(
        a: tuple[float, float, float], b: tuple[float, float, float]
    ) -> tuple[float, float, float]:
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    if dot(n, t_a) != 0.0 or dot(n, t_b) != 0.0:
        raise RuntimeError("mechanics-only amplitudes must be transverse")
    norm_sq = dot(n, n)
    c_a = tuple(1j * value / norm_sq for value in cross(n, t_a))
    c_b = tuple(1j * value / norm_sq for value in cross(n, t_b))
    if c_a == c_b:
        raise RuntimeError("amplitude-dependence witness unexpectedly collapsed")
    return {
        "classification": "autonomous mechanics only; not Kokuno/OpenAI candidate data",
        "same_phase_vector": n,
        "transverse_amplitude_a": t_a,
        "transverse_amplitude_b": t_b,
        "coefficient_a": tuple(str(value) for value in c_a),
        "coefficient_b": tuple(str(value) for value in c_b),
        "outputs_differ": True,
    }


def assert_clean(root: Path | None = None) -> None:
    contract = load_contract(root)
    violations = audit_contract(contract, root=root, check_live_files=True)
    if violations:
        raise RuntimeError(
            "CR002 source-amplitude input scope violations:\n- "
            + "\n- ".join(violations)
        )


if __name__ == "__main__":
    assert_clean()
    print(json.dumps(mechanics_only_amplitude_dependence_witness(), indent=2))
