"""Fail-closed CR002 audit for the source-specific symmetric stress-to-force map.

This module is governance only. It distinguishes the public corrected-reader
representation (two scalar stresses -> three cylindrical force components) from
current-lineage materialization. It does not construct pressure, forcing,
velocity corrections, residuals, or source-exact coefficients.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Mapping

SCHEMA = "cr002-kokuno-symmetric-stress-force-component-scope-v1"
TASK_ID = "CR002-KOKUNO-SYMMETRIC-STRESS-FORCE-COMPONENT-SCOPE-100"
CONTRACT_RELATIVE = Path("configs/kokuno_symmetric_stress_force_component_scope.json")

EXPECTED_BASE_HEAD = "dae4797fbd270a4f676b113f5ca3f7801d71a47e"
EXPECTED_A3_982_HEAD = "dd11348e8eecf298656d7da58aee4526c6102306"
EXPECTED_CR002_979_HEAD = "cef54441091d470efd137796909596c4018d5dea"
EXPECTED_READER_HEAD = "143f6773feb424ad9ed3a8d116653200f20346b7"
EXPECTED_READER_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"


def _root() -> Path:
    return Path(__file__).resolve().parents[2]


def load_contract(path: Path | None = None) -> dict[str, Any]:
    target = path or (_root() / CONTRACT_RELATIVE)
    return json.loads(target.read_text(encoding="utf-8"))


def mechanics_force(r: float, z: float) -> dict[str, float]:
    """Autonomous witness for sigma_1=r*z, sigma_2=r^2 under R33."""
    if r <= 0.0:
        raise ValueError("mechanics witness requires r>0")
    sigma_1 = r * z
    sigma_2 = r * r
    dz_sigma_1 = r
    dr_sigma_1 = z
    dr_sigma_2 = 2.0 * r
    return {
        "radial": dz_sigma_1,
        "tangential": dr_sigma_2 + 2.0 * sigma_2 / r,
        "axial": dr_sigma_1 + sigma_1 / r,
    }


def audit_contract(contract: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    require(contract.get("schema") == SCHEMA, "schema drift")
    require(contract.get("task_id") == TASK_ID, "task id drift")

    base = contract.get("base", {})
    require(base.get("a5_pr") == 984, "A5 base PR drift")
    require(base.get("head") == EXPECTED_BASE_HEAD, "A5 base head drift")
    require(
        base.get("active_constrained_integration_branch") == "codex/cr001-constraints",
        "active integration branch drift",
    )

    upstream = contract.get("upstream", {})
    a3 = upstream.get("agent3_pr_982", {})
    require(a3.get("head") == EXPECTED_A3_982_HEAD, "A3 #982 identity drift")
    require(a3.get("candidate_residual_evidence") is False, "radial force laundered into residual evidence")
    require(a3.get("complete_ns_defect") is False, "radial force laundered into complete NS defect")
    require(a3.get("authorized_as_correction_target") is False, "radial force prematurely authorized")
    require(
        upstream.get("cr002_pr_979", {}).get("head") == EXPECTED_CR002_979_HEAD,
        "CR002 #979 identity drift",
    )
    require(
        upstream.get("agent5_pr_984", {}).get("records_agent3_982_as_sibling_only") is True,
        "A5 #984 sibling-only scope lost",
    )

    provenance = contract.get("provenance", {})
    require(
        set(provenance.keys()) == {
            "user_requirement",
            "public_source_fact",
            "autonomous_design",
            "pending_unknown",
        },
        "four-way provenance partition drift",
    )
    public = provenance.get("public_source_fact", {})
    require(public.get("repository") == "KokunoYumeto/yang-mills-interacting-workbench", "reader repo drift")
    require(public.get("head") == EXPECTED_READER_HEAD, "reader head drift")
    require(public.get("blob") == EXPECTED_READER_BLOB, "reader blob drift")
    require(public.get("equations") == ["R32", "R33", "R34"], "reader equation scope drift")
    require(public.get("scalar_stresses") == ["sigma_2", "sigma_1"], "source scalar stress identity drift")
    require(public.get("reader_states_full_tensor_to_force_map") is True, "source full tensor-to-force statement lost")
    require(public.get("independent_third_scalar_radial_stress_present") is False, "invented third source scalar stress")
    require(
        public.get("force_components") == {
            "radial": "partial_z sigma_1",
            "tangential": "(partial_r + 2/r) sigma_2",
            "axial": "(partial_r + 1/r) sigma_1",
        },
        "source force component map drift",
    )
    require(
        "not an independent re-verification of the original OpenAI paper" in public.get("scope", ""),
        "corrected reader promoted to independently verified paper-exact fact",
    )

    autonomous = provenance.get("autonomous_design", [])
    require(any("mechanics-only" in item for item in autonomous), "mechanics witness not classified autonomous")
    pending = provenance.get("pending_unknown", [])
    require(any("complete three-component div(T) force vector" in item for item in pending), "full force vector no longer pending")
    require(any("complete NS defect" in item for item in pending), "complete NS defect no longer pending")
    require(any("Post-X_R" in item for item in pending), "global velocity completion no longer pending")
    require(any("paper exactness" in item for item in pending), "paper exactness no longer pending")

    resolution = contract.get("representation_resolution", {})
    require(resolution.get("generic_three_vector_requires_three_observed_components") is True, "generic witness truth was erased")
    require(resolution.get("generic_vector_witness_from_cr002_979_remains_logically_valid") is True, "CR002 #979 generic witness invalidated")
    require(resolution.get("generic_vector_witness_is_not_the_registered_source_tensor_representation") is True, "generic witness laundered into source tensor")
    require(resolution.get("source_specific_stress_scalar_channel_count") == 2, "source scalar stress count drift")
    require(resolution.get("source_specific_force_component_count") == 3, "source force component count drift")
    require(resolution.get("source_specific_independent_third_radial_stress_channel_required") is False, "redundant third radial stress reintroduced")
    require(resolution.get("source_specific_radial_force_obtained_from_axial_stress_cross_derivative") is True, "radial cross-derivative representation lost")
    require(resolution.get("source_specific_radial_force_formula") == "partial_z sigma_1", "radial force formula drift")
    require(resolution.get("cr002_979_alternative_representation_requirement_satisfied_by_public_reader") is True, "source-specific resolution of #979 fork lost")
    require(resolution.get("cr002_979_must_not_be_used_to_require_an_independent_third_radial_moment_inverse") is True, "#979 misused to require redundant inverse")
    require(resolution.get("agent3_982_current_radial_force_materialized") is True, "A3 #982 radial-force materialization lost")
    require(resolution.get("agent3_982_current_full_three_component_force_vector_materialized") is False, "radial force promoted to full force vector")
    require(resolution.get("agent3_982_complete_ns_defect_materialized") is False, "radial force promoted to complete defect")
    require(resolution.get("agent3_982_authorized_correction_target") is False, "radial force prematurely authorized")
    require(resolution.get("agent3_982_cartesian_correction_velocity_materialized") is False, "radial force promoted to correction velocity")

    witness = contract.get("mechanics_only_witness", {})
    require(
        witness.get("classification") == "autonomous_mechanics_only_not_kokuno_or_openai_data",
        "mechanics witness provenance laundering",
    )
    probe = witness.get("probe", {})
    try:
        observed = mechanics_force(float(probe["r"]), float(probe["z"]))
    except Exception as exc:
        errors.append(f"invalid mechanics witness probe: {exc}")
    else:
        expected = witness.get("derived_force", {})
        for key in ("radial", "tangential", "axial"):
            require(abs(observed[key] - float(expected.get(key, float("nan")))) <= 1e-12, f"mechanics witness {key} drift")
        require(all(abs(observed[key]) > 0.0 for key in observed), "mechanics witness became component-trivial")

    states = contract.get("truth_states", {})
    require(states.get("canonical_eq45_velocity_export_ready") is True, "Kokuno scope incorrectly downgraded Eq45 delivery")
    for key in (
        "canonical_eq45_visualization_ready",
        "canonical_eq45_visual_correspondence_verified",
        "canonical_eq45_pde_validated",
        "kokuno_velocity_export_ready",
        "kokuno_visual_correspondence_verified",
        "kokuno_pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        require(states.get(key) is False, f"forbidden truth-state promotion: {key}")

    cr001 = contract.get("cr001_snapshot", {})
    require(cr001.get("constraints_blob") == EXPECTED_CONSTRAINTS_BLOB, "canonical constraints identity drift")
    require(cr001.get("nu") == 0.01, "nu changed")
    require(cr001.get("physical_domain") == "R^3", "physical domain changed")
    require(cr001.get("evaluation_box") == [[-2, 2], [-2, 2], [-2, 2]], "evaluation box changed")
    require(cr001.get("support") == "r < 2 and abs(z) < 2", "support changed")
    require(cr001.get("time_interval") == [0.25, 0.75], "time interval changed")
    require(cr001.get("forcing_mode") == "restricted_two_parameter_family", "forcing family changed")
    require(cr001.get("reference_energy") == 1.0, "reference energy changed")
    require(cr001.get("reference_energy_abs_tolerance") == 0.001, "energy tolerance changed")
    require(cr001.get("validation_seed") == 914027, "validation seed changed")
    require(cr001.get("held_out_points") == 4096, "held-out count changed")
    require(cr001.get("derivative_steps") == [0.02, 0.01, 0.005], "derivative ladder changed")
    require(cr001.get("quadrature_orders_per_axis") == [24, 48, 96], "quadrature ladder changed")
    require(cr001.get("momentum_max_gate") == 0.001, "momentum max gate changed")
    require(cr001.get("momentum_l2_gate") == 0.001, "momentum L2 gate changed")
    require(cr001.get("divergence_max_gate") == 1e-5, "divergence max gate changed")
    require(cr001.get("divergence_l2_gate") == 1e-5, "divergence L2 gate changed")
    require(cr001.get("residual_defined_pointwise_free_force_allowed") is False, "free residual forcing enabled")
    require(cr001.get("candidate_collapse_allowed") is False, "candidate collapse enabled")
    require(cr001.get("post_hoc_threshold_relaxation_allowed") is False, "post-hoc threshold relaxation enabled")

    required_forbidden = {
        "two scalar source stresses -> current complete three-component force receipt",
        "current radial force -> complete NS defect",
        "current radial force -> authorized correction target",
        "current radial force -> Cartesian correction velocity",
        "public corrected-reader formula -> independently verified original-paper exactness",
        "scoped force/stress consistency -> PDE validation",
        "callable or visualizable velocity -> exact OpenAI field",
    }
    require(required_forbidden.issubset(set(contract.get("forbidden_promotions", []))), "forbidden-promotion firewall weakened")

    return errors


def assert_contract(contract: Mapping[str, Any] | None = None) -> None:
    payload = contract if contract is not None else load_contract()
    errors = audit_contract(payload)
    if errors:
        raise AssertionError("\n".join(errors))


def mutated(contract: Mapping[str, Any], path: tuple[str, ...], value: Any) -> dict[str, Any]:
    out = copy.deepcopy(dict(contract))
    cursor: dict[str, Any] = out
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    return out


if __name__ == "__main__":
    assert_contract()
    print("CR002 symmetric stress-to-force component scope: PASS")
