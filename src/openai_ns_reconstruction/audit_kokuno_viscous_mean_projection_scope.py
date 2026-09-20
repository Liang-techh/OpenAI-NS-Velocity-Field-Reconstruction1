"""Fail-closed CR002 audit for the Kokuno strict-inner viscous-mean projection scope.

Agent-3 PR #841 projects the Cartesian vector viscous term ``-nu*Delta u``
onto a rotating cylindrical basis and then takes an angular mean.  That object
must not be silently reinterpreted as a componentwise scalar Laplacian applied
to cylindrical mean components: radial/tangential vector-Laplacian components
carry geometric basis terms.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
from typing import Any

SCHEMA = "cr002-kokuno-viscous-mean-projection-scope-v1"
CONTRACT_PATH = "configs/kokuno_viscous_mean_projection_scope.json"
UPSTREAM_PATH = (
    "src/openai_ns_reconstruction/"
    "kokuno_inner_leading_oscillatory_viscous_mean.py"
)
UPSTREAM_BLOB = "2801435cd1f968194cd819bf19d32468b537f71d"

_REQUIRED_CLASSIFICATIONS = {
    "user_requirement",
    "public_source_fact",
    "autonomous_design",
    "pending_unknown",
}
_REQUIRED_UPSTREAM_CONSTANTS = {
    "PARENT_AGENT3_PR": 831,
    "AGENT1_PR": 839,
    "AGENT2_PR": 708,
    "AGENT1_LAPLACIAN_STEP": 1.0e-3,
    "AGENT2_LAPLACIAN_STEP": 4.5e-3,
    "REPOSITORY_VISCOSITY": 0.01,
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object at {path}")
    return value


def _git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()  # noqa: S324 - Git identity


def _literal_module_constants(source: str) -> dict[str, Any]:
    tree = ast.parse(source)
    values: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        if isinstance(node, ast.Assign):
            if len(node.targets) != 1 or not isinstance(node.targets[0], ast.Name):
                continue
            name = node.targets[0].id
            value_node = node.value
        else:
            if not isinstance(node.target, ast.Name) or node.value is None:
                continue
            name = node.target.id
            value_node = node.value
        try:
            values[name] = ast.literal_eval(value_node)
        except (ValueError, TypeError):
            continue
    return values


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def manufactured_vector_laplacian_witness(radius: float = 0.7) -> dict[str, float]:
    """Show why scalar-component and vector-Laplacian semantics differ.

    For the smooth Cartesian field ``u=(x,y,0)``, the cylindrical radial
    component is ``u_r=r``.  The Cartesian vector Laplacian is exactly zero.
    A naive axisymmetric scalar Laplacian of ``u_r`` gives ``1/r``; the
    cylindrical vector-Laplacian geometric term ``-u_r/r^2`` cancels it.
    """
    if radius <= 0.0:
        raise ValueError("radius must be positive for this off-axis witness")
    u_r = radius
    scalar_component_laplacian = 1.0 / radius
    geometric_basis_correction = -u_r / (radius * radius)
    vector_radial_laplacian = scalar_component_laplacian + geometric_basis_correction
    return {
        "radius": radius,
        "u_r": u_r,
        "naive_scalar_radial_laplacian": scalar_component_laplacian,
        "geometric_basis_correction": geometric_basis_correction,
        "vector_radial_laplacian": vector_radial_laplacian,
        "cartesian_vector_laplacian_radial_projection": 0.0,
    }


def _validate_contract(
    contract: dict[str, Any],
    cr001: dict[str, Any],
    upstream_constants: dict[str, Any],
    upstream_source: str,
) -> None:
    _require(contract.get("schema") == SCHEMA, "scope contract schema drift")
    _require(
        set(contract.get("classification", {})) == _REQUIRED_CLASSIFICATIONS,
        "four-way provenance classification drift",
    )

    upstream = contract["upstream"]
    _require(upstream["pr"] == 841, "upstream PR drift")
    _require(
        upstream["head"] == "d8b0e792d2a263bb97754d99ee6be998611c90c3",
        "upstream head drift",
    )
    _require(upstream["module"] == UPSTREAM_PATH, "upstream module drift")
    _require(upstream["module_blob"] == UPSTREAM_BLOB, "upstream blob drift")

    for name, expected in _REQUIRED_UPSTREAM_CONSTANTS.items():
        _require(upstream_constants.get(name) == expected, f"upstream constant drift: {name}")

    required_source_fragments = (
        "inner_viscous = -REPOSITORY_VISCOSITY * inner_laplacian",
        "oscillatory_viscous = -REPOSITORY_VISCOSITY * oscillatory_laplacian",
        "mean_total = _project_cartesian_ring_to_cylindrical_mean(total_viscous, c, s)",
        "for order in ANGULAR_ORDERS",
        '"complete_ns_defect": False',
    )
    for fragment in required_source_fragments:
        _require(fragment in upstream_source, f"upstream operator semantics drift: {fragment}")

    semantics = contract["upstream_operator_semantics"]
    _require(semantics["nu"] == 0.01, "viscosity semantics drift")
    _require(semantics["inner_laplacian_step"] == 1.0e-3, "inner step drift")
    _require(semantics["oscillatory_laplacian_step"] == 4.5e-3, "oscillatory step drift")
    _require(semantics["angular_orders"] == [32, 64, 128], "angular ladder drift")
    for false_key in (
        "naive_componentwise_scalar_laplacian_of_mean_components",
        "operator_on_mean_equivalence_established_here",
        "radial_and_axial_derivatives_of_mean_components_computed_here",
        "complete_ns_defect",
        "directly_comparable_to_st006_full_residual",
    ):
        _require(semantics[false_key] is False, f"operator-scope laundering: {false_key}")

    boundary = contract["cylindrical_vector_laplacian_boundary"]
    _require("- r^-2" in boundary["axisymmetric_m0_radial_component"], "radial geometric term lost")
    _require("- r^-2" in boundary["axisymmetric_m0_tangential_component"], "tangential geometric term lost")
    _require("r^-2" not in boundary["axisymmetric_m0_axial_component"], "spurious axial geometric term")
    required_nonclaims = {
        "naive scalar-Laplacian equivalence for radial or tangential mean components",
        "full operator-on-mean closure",
        "global/corrected-candidate viscous residual",
        "complete momentum residual",
    }
    _require(
        required_nonclaims.issubset(set(boundary["not_established_by_841"])),
        "vector-Laplacian promotion boundary weakened",
    )

    witness = manufactured_vector_laplacian_witness()
    _require(
        abs(witness["vector_radial_laplacian"]) <= 1.0e-15,
        "manufactured vector-Laplacian witness drift",
    )
    _require(
        witness["naive_scalar_radial_laplacian"] > 1.0,
        "manufactured naive scalar witness unexpectedly vanished",
    )
    _require(
        witness["naive_scalar_radial_laplacian"]
        != witness["cartesian_vector_laplacian_radial_projection"],
        "manufactured witness no longer distinguishes operator semantics",
    )

    canonical = contract["canonical_cr001"]
    validation = cr001["validation"]
    forcing = cr001["forcing"]
    nontriviality = cr001["nontriviality"]
    domain = cr001["domain"]

    _require(cr001["nu"] == canonical["nu"] == 0.01, "nu drift")
    _require(domain["physical"] == canonical["physical_domain"] == "R^3", "physical domain drift")
    _require(domain["evaluation_box"] == canonical["evaluation_box"], "evaluation box drift")
    _require(domain["support"] == canonical["support"], "support drift")
    _require(domain["time_interval"] == canonical["time_interval"], "time interval drift")
    _require(
        forcing["mode"] == canonical["forcing_mode"] == "restricted_two_parameter_family",
        "forcing mode drift",
    )
    _require(
        "No residual-dependent basis or pointwise free force" in forcing["restriction"],
        "free-force prohibition missing",
    )
    _require(canonical["free_residual_defined_force_allowed"] is False, "free residual force enabled")
    _require(nontriviality["reference_energy"] == canonical["reference_energy"] == 1.0, "energy target drift")
    _require(
        nontriviality["reference_energy_abs_tolerance"]
        == canonical["reference_energy_abs_tolerance"]
        == 1.0e-3,
        "energy tolerance drift",
    )
    _require("reject collapsed candidates" in nontriviality["enforcement"], "amplitude collapse guard missing")
    _require(canonical["reject_amplitude_collapse"] is True, "amplitude collapse allowed")

    _require(validation["seed"] == canonical["validation_seed"] == 914027, "validation seed drift")
    _require(validation["held_out_points"] == canonical["held_out_points"] == 4096, "held-out count drift")
    _require(validation["times"] == canonical["validation_times"], "validation times drift")
    _require(validation["derivative_steps"] == canonical["derivative_steps"], "derivative ladder drift")
    _require(
        validation["quadrature_orders_per_axis"] == canonical["quadrature_orders_per_axis"],
        "quadrature ladder drift",
    )
    thresholds = validation["thresholds"]
    _require(thresholds["pde_residual_max"] == canonical["momentum_max_threshold"] == 1.0e-3, "momentum max threshold drift")
    _require(thresholds["pde_residual_L2"] == canonical["momentum_l2_threshold"] == 1.0e-3, "momentum L2 threshold drift")
    _require(thresholds["divergence_max"] == canonical["divergence_max_threshold"] == 1.0e-5, "divergence max threshold drift")
    _require(thresholds["divergence_L2"] == canonical["divergence_l2_threshold"] == 1.0e-5, "divergence L2 threshold drift")
    _require("changing thresholds requires a new experiment version" in validation["failure_policy"], "post-hoc threshold guard missing")
    _require(canonical["posthoc_threshold_relaxation_allowed"] is False, "post-hoc relaxation enabled")

    states = contract["states"]
    _require(states["viscous_mean_projection_registered"] is True, "upstream projection registration lost")
    _require(states["upstream_841_ci_passed"] is None, "queued/unresolved CI pre-promoted")
    for false_key in (
        "operator_on_mean_equivalence_established",
        "canonical_cr001_momentum_assessed_by_this_contract",
        "st006_comparable",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
    ):
        _require(states[false_key] is False, f"truth-state promotion: {false_key}")

    promotion = contract["promotion_boundary"]
    required_promotion_nonclaims = {
        "naive_componentwise_scalar_laplacian_of_cylindrical_mean",
        "operator_on_mean_equivalence",
        "complete_global_kokuno_candidate",
        "canonical_cr001_momentum_acceptance",
        "st006_same_protocol_improvement",
        "pde_validated",
        "visual_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    }
    _require(
        required_promotion_nonclaims.issubset(
            set(promotion["future_841_pass_does_not_establish"])
        ),
        "promotion boundary weakened",
    )
    _require(
        promotion["upstream_failure_invalidates_existing_callable_velocity_delivery"] is False,
        "scoped upstream failure incorrectly blocks callable velocity delivery",
    )


def audit(
    contract_path: Path | None = None,
    cr001_path: Path | None = None,
    upstream_path: Path | None = None,
) -> dict[str, Any]:
    root = _repo_root()
    contract_path = contract_path or root / CONTRACT_PATH
    cr001_path = cr001_path or root / "configs/constraints.json"
    upstream_path = upstream_path or root / UPSTREAM_PATH

    contract = _load_json(contract_path)
    cr001 = _load_json(cr001_path)
    upstream_bytes = upstream_path.read_bytes()
    actual_blob = _git_blob_sha(upstream_bytes)
    _require(actual_blob == UPSTREAM_BLOB, f"upstream source blob drift: {actual_blob}")
    upstream_source = upstream_bytes.decode("utf-8")
    upstream_constants = _literal_module_constants(upstream_source)
    _validate_contract(contract, cr001, upstream_constants, upstream_source)
    return {
        "schema": SCHEMA,
        "status": "pass",
        "upstream_blob": actual_blob,
        "scope": "cartesian_vector_viscous_mean_projection_not_naive_scalar_laplacian_of_cylindrical_mean",
        "manufactured_witness": manufactured_vector_laplacian_witness(),
        "pde_validated": False,
    }


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
