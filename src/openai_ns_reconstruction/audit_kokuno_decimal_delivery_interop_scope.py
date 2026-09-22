"""CR002 audit for precision-qualified Kokuno velocity delivery interoperability.

This module is governance-only. It does not alter candidate mathematics, pressure,
forcing, residuals, validation thresholds, or scientific readiness states.

The exact parent PR #1179 closes the split-channel representation seam with a
fixed-precision Decimal total and decimal-string JSON stage export. This auditor
keeps a separate boundary between that precision-qualified stage surface and any
future binary64/MATLAB-native/global delivery claim.
"""

from __future__ import annotations

import ast
import json
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from pathlib import Path
from typing import Any, Mapping

CONTRACT_PATH = Path("configs/kokuno_decimal_delivery_interop_scope.json")
PARENT_SOURCE_PATH = Path(
    "src/openai_ns_reconstruction/"
    "kokuno_pa16_current_cartesian_relative_swirl_decimal_composed.py"
)
CANONICAL_CONSTRAINTS_PATH = Path("configs/constraints.json")

EXPECTED_PARENT_HEAD = "8c5c5b6d55a68de8285ed1dd0ebb55128f3078a4"
EXPECTED_PARENT_SOURCE_BLOB = "78247e3f9e21ddc19f68e109d8fe1a069ba0360f"
EXPECTED_CONSTRAINTS_BLOB = "6c559e42895a606e2ef025ade4cb448966d75814"
EXPECTED_DECIMAL_DIGITS = 96
EXPECTED_ROUNDING = "ROUND_HALF_EVEN"


class GovernanceError(ValueError):
    """Raised when a representation or truth-boundary contract drifts."""


def _literal_assignment(tree: ast.AST, name: str) -> Any:
    for node in getattr(tree, "body", []):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == name:
                    return ast.literal_eval(node.value)
    raise GovernanceError(f"missing module assignment {name}")


def _method_source_requirements(tree: ast.AST) -> None:
    methods: dict[str, ast.FunctionDef] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            methods[node.name] = node

    for required in (
        "velocity",
        "velocity_decimal_split",
        "export_velocity_decimal_json",
        "configuration",
        "save_configuration",
        "load_configuration",
    ):
        if required not in methods:
            raise GovernanceError(f"parent missing required method {required}")

    velocity_dump = ast.dump(methods["velocity"], include_attributes=False)
    if "_compose_decimal" not in velocity_dump:
        raise GovernanceError("parent velocity no longer composes through _compose_decimal")
    if "object" not in velocity_dump:
        raise GovernanceError("parent velocity no longer materializes an object/Decimal array")

    export_method = methods["export_velocity_decimal_json"]
    exact_string_export = False
    for node in ast.walk(export_method):
        if not isinstance(node, ast.Dict):
            continue
        for key, value in zip(node.keys, node.values):
            if not (
                isinstance(key, ast.Constant)
                and key.value == "values_row_major"
                and isinstance(value, ast.ListComp)
                and isinstance(value.elt, ast.Call)
                and isinstance(value.elt.func, ast.Name)
                and value.elt.func.id == "str"
            ):
                continue
            exact_string_export = True
    if not exact_string_export:
        raise GovernanceError("parent export no longer preserves Decimal values as strings")


def audit_parent_source_text(source_text: str, contract: Mapping[str, Any]) -> dict[str, Any]:
    """Fail closed if the exact-parent representation/truth seam drifts."""
    tree = ast.parse(source_text)

    if _literal_assignment(tree, "PARENT_EXACT_HEAD") != (
        "657818dd83e119e6091bd2490d3804c04c3ef723"
    ):
        raise GovernanceError("unexpected #1179 internal parent identity")
    if _literal_assignment(tree, "DECIMAL_DIGITS") != EXPECTED_DECIMAL_DIGITS:
        raise GovernanceError("Decimal precision drift")
    if _literal_assignment(tree, "PRECISION_MARGIN_DIGITS") != 24:
        raise GovernanceError("Decimal precision-margin drift")

    numerical = _literal_assignment(tree, "_NUMERICAL_REALIZATION")
    truth = _literal_assignment(tree, "_TRUTH_UPDATES")

    if "96" not in numerical.get("total_representation", ""):
        raise GovernanceError("fixed 96-digit representation declaration missing")
    if "Decimal.from_float" not in numerical.get("binary64_embedding", ""):
        raise GovernanceError("binary64 embedding declaration drift")
    if numerical.get("rounding") != EXPECTED_ROUNDING:
        raise GovernanceError("rounding declaration drift")

    expected_truth = contract["parent_required_truth"]
    for key, expected in expected_truth.items():
        if truth.get(key) is not expected:
            raise GovernanceError(
                f"parent truth drift for {key}: expected {expected!r}, got {truth.get(key)!r}"
            )

    _method_source_requirements(tree)

    return {
        "parent_head": contract["exact_parent"]["head"],
        "parent_source_blob": contract["exact_parent"]["source_blob"],
        "decimal_digits": EXPECTED_DECIMAL_DIGITS,
        "rounding": EXPECTED_ROUNDING,
        "binary64_total_velocity_export_ready": truth[
            "binary64_total_velocity_export_ready"
        ],
        "unified_global_cartesian_velocity_export_ready": truth[
            "unified_global_cartesian_velocity_export_ready"
        ],
        "pde_validated": truth["pde_validated"],
    }


def audit_canonical_constraints(
    constraints: Mapping[str, Any], contract: Mapping[str, Any]
) -> dict[str, Any]:
    """Confirm this delivery firewall does not mutate CR001's frozen gates."""
    expected = contract["canonical_constraints"]

    checks = {
        "nu": constraints["nu"],
        "physical_domain": constraints["domain"]["physical"],
        "evaluation_box": constraints["domain"]["evaluation_box"],
        "support": constraints["domain"]["support"],
        "time_interval": constraints["domain"]["time_interval"],
        "forcing_mode": constraints["forcing"]["mode"],
        "reference_energy": constraints["nontriviality"]["reference_energy"],
        "reference_energy_abs_tolerance": constraints["nontriviality"][
            "reference_energy_abs_tolerance"
        ],
        "validation_seed": constraints["validation"]["seed"],
        "held_out_points": constraints["validation"]["held_out_points"],
        "derivative_steps": constraints["validation"]["derivative_steps"],
        "quadrature_orders_per_axis": constraints["validation"][
            "quadrature_orders_per_axis"
        ],
        "momentum_max": constraints["validation"]["thresholds"]["pde_residual_max"],
        "momentum_L2": constraints["validation"]["thresholds"]["pde_residual_L2"],
        "divergence_max": constraints["validation"]["thresholds"]["divergence_max"],
        "divergence_L2": constraints["validation"]["thresholds"]["divergence_L2"],
    }
    for key, value in checks.items():
        if value != expected[key]:
            raise GovernanceError(
                f"canonical constraint drift for {key}: {value!r} != {expected[key]!r}"
            )

    restriction = constraints["forcing"]["restriction"].lower()
    if "no residual-dependent" not in restriction or "pointwise free force" not in restriction:
        raise GovernanceError("free-force prohibition drift")
    enforcement = constraints["nontriviality"]["enforcement"].lower()
    if "reject collapsed candidates" not in enforcement:
        raise GovernanceError("nontriviality/collapse prohibition drift")
    if "changing thresholds requires a new experiment version" not in constraints[
        "validation"
    ]["failure_policy"].lower():
        raise GovernanceError("threshold-relaxation policy drift")

    return checks


def autonomous_transport_witness() -> dict[str, Any]:
    """Mechanics-only witness: Decimal survives while binary64 downcast erases delta."""
    base = 1.0
    delta = 2.0 ** -60

    with localcontext() as ctx:
        ctx.prec = EXPECTED_DECIMAL_DIGITS
        ctx.rounding = ROUND_HALF_EVEN
        base_d = +Decimal.from_float(base)
        delta_d = +Decimal.from_float(delta)
        total_d = +(base_d + delta_d)
        replay_delta = +(total_d - base_d)

    decimal_string = str(total_d)
    parsed_string = Decimal(decimal_string)
    binary64_total = float(total_d)

    return {
        "classification": "autonomous_mechanics_only",
        "base_binary64": base,
        "delta_binary64": delta,
        "binary64_base_plus_delta_equals_base": (base + delta) == base,
        "decimal_total_differs_from_base": total_d != base_d,
        "decimal_total_minus_base_replays_delta": replay_delta == delta_d,
        "decimal_string_roundtrip_exact": parsed_string == total_d,
        "binary64_downcast_equals_base": binary64_total == base,
        "binary64_downcast_erases_nonzero_delta": (
            binary64_total == base and delta != 0.0
        ),
    }


def audit_delivery_representation_claim(
    claim: Mapping[str, Any], contract: Mapping[str, Any]
) -> None:
    """Reject evidence transfer across an un-replayed numeric representation change."""
    required = set(contract["delivery_representation_identity"]["required_fields"])
    missing = sorted(required - set(claim))
    if missing:
        raise GovernanceError(f"delivery representation identity missing: {missing}")

    parent = contract["delivery_representation_identity"]
    exact_parent_identity = {
        "scalar_type": parent["parent_scalar_type"],
        "precision_digits": parent["parent_precision_digits"],
        "rounding": parent["parent_rounding"],
        "binary64_embedding": parent["parent_binary64_embedding"],
        "export_encoding": parent["parent_export_encoding"],
    }
    changed = any(claim[key] != expected for key, expected in exact_parent_identity.items())

    survival_replayed = bool(claim.get("correction_survival_replayed", False))
    survival_observed = bool(claim.get("correction_survives_transport", False))
    save_load_replayed = bool(claim.get("save_load_replayed", False))

    if changed and not survival_replayed:
        raise GovernanceError(
            "changed delivery representation cannot inherit parent survival evidence"
        )
    if changed and not save_load_replayed:
        raise GovernanceError(
            "changed delivery representation requires save/load replay"
        )
    if claim.get("velocity_export_ready", False) and not survival_observed:
        raise GovernanceError(
            "velocity_export_ready cannot be promoted when the nonzero correction is erased"
        )
    if claim.get("binary64_total_velocity_export_ready", False) and (
        claim.get("scalar_type") == "binary64" and not survival_observed
    ):
        raise GovernanceError("binary64 export promotion erases the audited correction")

    for key in ("pde_validated", "paper_exact", "openai_field_identified", "blowup_proved"):
        if claim.get(key, False):
            raise GovernanceError(f"representation receipt cannot promote {key}")


def load_contract(path: str | Path = CONTRACT_PATH) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def run_repo_audit(root: str | Path = ".") -> dict[str, Any]:
    root = Path(root)
    contract = load_contract(root / CONTRACT_PATH)
    parent_text = (root / PARENT_SOURCE_PATH).read_text()
    constraints = json.loads((root / CANONICAL_CONSTRAINTS_PATH).read_text())

    parent_receipt = audit_parent_source_text(parent_text, contract)
    constraint_receipt = audit_canonical_constraints(constraints, contract)
    witness = autonomous_transport_witness()

    if not witness["decimal_total_differs_from_base"]:
        raise GovernanceError("mechanics witness failed: Decimal did not preserve delta")
    if not witness["binary64_downcast_erases_nonzero_delta"]:
        raise GovernanceError("mechanics witness failed: binary64 did not expose the seam")

    return {
        "schema": "cr002-kokuno-decimal-delivery-interop-audit-receipt-v1",
        "contract_schema": contract["schema"],
        "exact_parent_head": EXPECTED_PARENT_HEAD,
        "exact_parent_source_blob": EXPECTED_PARENT_SOURCE_BLOB,
        "canonical_constraints_blob": EXPECTED_CONSTRAINTS_BLOB,
        "parent": parent_receipt,
        "constraints": constraint_receipt,
        "mechanics_witness": witness,
        "scientific_state_promoted": False,
    }
