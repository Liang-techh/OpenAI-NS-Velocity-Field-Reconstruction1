import json
from pathlib import Path

from openai_ns_reconstruction.constrained_constraint_governance import audit_constraint_family


def _base():
    return {
        "schema_version": 1,
        "status": "preregistered_not_validated",
        "experiment_id": "base",
        "source_map": "docs/CONSTRAINT_SOURCES.md",
        "units": "dimensionless",
        "equation": "NS",
        "nu": 0.01,
        "domain": {
            "physical": "R3",
            "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
            "support": "compact",
            "boundary": "zero extension",
            "time_interval": [0.25, 0.75],
            "initial_condition": "nonzero",
            "pressure_gauge": "p=0 outside",
        },
        "structure": {
            "h": 0.005,
            "tau": "1-t",
            "radial_scale": "sqrt(tau)",
            "axial_scale": "tau^D",
            "leading_swirl_axial_scale": "tau^-A",
            "symmetry": "axisymmetric with swirl",
            "core_probe": "probe",
            "core_sign_requirements": ["u_r < 0", "u_theta > 0", "u_z > 0"],
        },
        "forcing": {
            "mode": "restricted",
            "restriction": "No residual-dependent basis or pointwise free force.",
        },
        "nontriviality": {
            "reference_time": 0.25,
            "reference_energy": 1.0,
        },
        "optimization": {
            "seed": 11,
            "sampling": "train",
            "interior_points": 32,
            "maximum_function_evaluations": 50,
            "loss_weights": {"pde": 1.0},
            "hard_constraints": ["divergence"],
            "candidate_parameter_bounds": {
                "amplitude": "derived, not optimized; reject outside [1e-4,100]"
            },
        },
        "validation": {
            "seed": 22,
            "held_out_points": 64,
            "times": [0.25, 0.5, 0.75],
            "operator": "independent",
            "derivative_steps": [0.02, 0.01, 0.005],
            "quadrature_orders_per_axis": [24, 48, 96],
            "vary_one_at_a_time": True,
            "norms": ["max", "L2"],
            "thresholds": {"pde": 0.001, "divergence": 1e-5},
            "failure_policy": "retain failed results",
            "time_derivatives": "one-sided endpoints",
        },
        "excluded_claims": ["blow-up proof"],
    }


def _write(path: Path, data):
    path.write_text(json.dumps(data), encoding="utf-8")


def test_audit_accepts_parameter_growth_but_preserves_preregistered_contract(tmp_path):
    first = _base()
    second = json.loads(json.dumps(first))
    second["experiment_id"] = "extended"
    second["optimization"]["candidate_parameter_bounds"]["new_basis"] = [-1, 1]
    second["parent_experiment"] = "base"
    second["change_reason"] = "add one bounded basis direction"

    p1 = tmp_path / "constraints.json"
    p2 = tmp_path / "constraints_extended.json"
    _write(p1, first)
    _write(p2, second)

    report = audit_constraint_family([p1, p2])
    assert report.protected_contract_pass
    assert report.semantic_errors == ()
    assert report.identity_collisions == {}
    assert report.governance_pass


def test_audit_catches_threshold_force_seed_and_amplitude_regressions(tmp_path):
    first = _base()
    broken = json.loads(json.dumps(first))
    broken["experiment_id"] = "broken"
    broken["validation"]["thresholds"]["pde"] = 0.1
    broken["validation"]["seed"] = broken["optimization"]["seed"]
    broken["forcing"]["restriction"] = "fit any pointwise residual"
    broken["optimization"]["candidate_parameter_bounds"]["amplitude"] = [-100, 100]

    p1 = tmp_path / "constraints.json"
    p2 = tmp_path / "constraints_broken.json"
    _write(p1, first)
    _write(p2, broken)

    report = audit_constraint_family([p1, p2])
    assert not report.protected_contract_pass
    assert any("validation.thresholds differs" in item for item in report.protected_mismatches)
    assert not report.training_validation_seeds_separate
    assert any("residual-defined/free forcing" in item for item in report.semantic_errors)
    assert any("amplitude-collapse guard" in item for item in report.semantic_errors)
    assert not report.governance_pass


def test_audit_surfaces_duplicate_experiment_ids_without_hiding_contract_state(tmp_path):
    first = _base()
    second = json.loads(json.dumps(first))
    second["experiment_id"] = first["experiment_id"]

    p1 = tmp_path / "constraints.json"
    p2 = tmp_path / "constraints_duplicate.json"
    _write(p1, first)
    _write(p2, second)

    report = audit_constraint_family([p1, p2])
    assert report.protected_contract_pass
    assert report.semantic_errors == ()
    assert report.identity_collisions == {
        "base": ("constraints.json", "constraints_duplicate.json")
    }
    assert not report.governance_pass


def test_repository_constraint_family_keeps_protected_contract_if_present():
    root = Path(__file__).resolve().parents[1]
    paths = sorted((root / "configs").glob("constraints*.json"))
    if not paths:
        return

    report = audit_constraint_family(paths)
    assert report.protected_contract_pass
    assert report.semantic_errors == ()
    assert report.training_validation_seeds_separate
    assert report.governance_pass == (not report.identity_collisions)
