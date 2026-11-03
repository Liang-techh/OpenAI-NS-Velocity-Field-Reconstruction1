"""Bridge the saved v3 Green proof to the corrected parent provenance.

This leaf deliberately consumes *mappings*, not source files.  The v3 source
and AST proof was already run once by the owning runner.  The only task here
is to admit its saved ``actual_result`` against the direct pins in the cutoff
report.  In particular, the first-raw source is pinned by the parent's direct
``implementation_hashes`` mapping; it is not looked up in ``input_hashes``.

The acceptance remains conditional: it preserves the saved function-level
Green proof and its fixed-x/fixed-hb/nonzero-scale premises.  It does not
re-prove a live owner, a complete ``Kchi`` resolvent, a numerical residual, or
any ledger/value claim.  No AST parser, source import, source hash, or old
leaf is used in this module.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
from pathlib import Path
from typing import Any


SCHEMA = "n1_green_scale_provenance_bridge_v1"
V3_SCHEMA = "n1_green_scale_commutator_v3_result_v1"

V3_RESULT_PATH = Path("tmp/n1_green_scale_commutator_v3/result.json")
PARENT_REPORT_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
)
PARENT_REPORT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"

VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
FIRST_RAW_SOURCE = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"

VOLTERRA_SHA256 = "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
REGULAR_SHA256 = "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555"
FIRST_RAW_SHA256 = "f13f29b4d0652e8ed3b3347f19bab7eb1ea0e18fd4be0f50cf6d550a1aeb92c6"

ROWS = ("5", "6", "6", "4", "5", "6")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
SCALE_DIAGONAL = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)
CANONICAL_INPUT_COLUMNS = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")

EXPECTED_SCALE_RECORDS = [
    {
        "output_index": 1,
        "output_field": "F1",
        "input_row": "5",
        "input_column": "F1x",
        "output_scale": "F0(Z)/epsilon",
        "input_scale": "F0(Z)/epsilon",
        "same_scale": True,
    },
    {
        "output_index": 2,
        "output_field": "Uz1",
        "input_row": "6",
        "input_column": "Uz1x",
        "output_scale": "epsilon",
        "input_scale": "epsilon",
        "same_scale": True,
    },
    {
        "output_index": 3,
        "output_field": "K1",
        "input_row": "6",
        "input_column": "Uz1x",
        "output_scale": "epsilon",
        "input_scale": "epsilon",
        "same_scale": True,
    },
    {
        "output_index": 4,
        "output_field": "P1",
        "input_row": "4",
        "input_column": "P1",
        "output_scale": "epsilon*Pstar^2",
        "input_scale": "epsilon*Pstar^2",
        "same_scale": True,
    },
    {
        "output_index": 5,
        "output_field": "d_x_F1",
        "input_row": "5",
        "input_column": "F1x",
        "output_scale": "F0(Z)/epsilon",
        "input_scale": "F0(Z)/epsilon",
        "same_scale": True,
    },
    {
        "output_index": 6,
        "output_field": "d_x_Uz1",
        "input_row": "6",
        "input_column": "Uz1x",
        "output_scale": "epsilon",
        "input_scale": "epsilon",
        "same_scale": True,
    },
]
EXPECTED_SCALE_MAPPING = {
    "rows": list(ROWS),
    "records": EXPECTED_SCALE_RECORDS,
    "all_output_input_scales_equal": True,
}

EXPECTED_PROOF_IDENTITY = "S^-1*G*S = G"
EXPECTED_ENTRYWISE_ARGUMENT = (
    "G[i,j] != 0 only for j=ROWS[i]-1, and SCALE_DIAGONAL[i]=SCALE_DIAGONAL[j]; "
    "hence S_i^-1 G[i,j] S_j=G[i,j]"
)
EXPECTED_A_ARGUMENT = (
    "conditional algebraic reuse: if A=-2*chi*E51, then S1=S5=F0(Z)/epsilon "
    "gives S^-1*A*S=A; the external A support is not reverified by this leaf"
)
EXPECTED_GREEN_KERNEL_Z_DERIVATIVE = (
    "G_Z=0 under fixed x, fixed hb, and Z-independent endpoints/phase bounds"
)
EXPECTED_CONDITIONAL_RESIDUAL = {
    "premise_same_function_owner": "the raw and hatted six-component functions are the same owner",
    "premise_Qhat": "Qhat=S^-1*Q*S, with B1*S_Z retained in the ordinary-Z connection",
    "value_identity": "d_hat=S^-1*d_raw",
    "fixed_x_ordinary_Z_identity": "d_hat_Z=S^-1*d_raw_Z-Gamma*d_hat, Gamma=S^-1*S_Z",
    "source_owner_established_here": False,
}
EXPECTED_ORDINARY_Z_CONNECTION = {
    "S_Z": "diag(F0_Z/epsilon,0,0,0,F0_Z/epsilon,0)",
    "Gamma": "S^-1*S_Z=diag(F0_Z/F0,0,0,0,F0_Z/F0,0)",
    "B0_transformed": "S^-1*(B0*S+B1*S_Z)",
    "B1_transformed": "S^-1*B1*S",
    "g_transformed": "S^-1*g",
    "B1_S_Z_retained": True,
    "A_and_B1_source_reverified_here": False,
}
EXPECTED_PARENT_RAW_SCALE_BINDING = {
    "parent_actual_first_raw_scale_conjugacy": True,
    "parent_raw_preconditioned_resolvent_factorization": True,
    "parent_normalized_Kchi_kernel_substituted_is_false": True,
    "source_AST_scale_rows_match": True,
    "compatibility_check": True,
}
PROOF_LIMITATION_FALSE_FLAGS = (
    "complete_Kchi_resolvent_provided",
    "raw_numerical_residual_evaluated",
    "gauge_owner_proved",
    "historical_error_removed",
    "historical_error_added",
    "interval_or_tuple_overlap_used_as_function_identity",
    "source_provider_executed",
    "full_B0_B1_operator_conjugacy_beyond_the_proved_connection",
    "A_and_B1_support_reverified_here",
)

EXPECTED_SOURCE_HASHES = {
    VOLTERRA_SOURCE: VOLTERRA_SHA256,
    REGULAR_SOURCE: REGULAR_SHA256,
    FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
}
EXPECTED_SOURCE_PATHS = {
    "volterra": VOLTERRA_SOURCE,
    "regular": REGULAR_SOURCE,
    "first_raw": FIRST_RAW_SOURCE,
}
EXPECTED_V3_INPUT_HASHES = {
    str(PARENT_REPORT_PATH).replace("\\", "/"): PARENT_REPORT_SHA256,
    VOLTERRA_SOURCE: VOLTERRA_SHA256,
    REGULAR_SOURCE: REGULAR_SHA256,
    FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
}
EXPECTED_PER_SOURCE = {
    VOLTERRA_SOURCE: True,
    REGULAR_SOURCE: True,
    FIRST_RAW_SOURCE: False,
}

# These are the only false checks in the archived v3 actual result.  They
# reflect the old parent-admission bug, not a source or AST defect.
V3_PARENT_HASH_EXEMPT_FALSE = frozenset(
    {
        "parent_report_input_hashes_present",
        "parent_report_input_hashes_match",
        "parent_report_passed",
    }
)

PARENT_FIRST_RAW_FLAGS = (
    "actual_first_raw_source_definition_bound",
    "actual_first_raw_function_bound",
    "actual_first_raw_equation_bound",
    "actual_first_raw_core_entry_bound",
    "actual_first_raw_ordinary_Z0_to_3_bound",
    "actual_first_raw_scale_conjugacy",
    "actual_first_raw_all_word_convergence_bound",
    "raw_all_word_convergence_certified",
    "raw_preconditioned_resolvent_factorization",
)
PARENT_FIRST_RAW_NESTED_FLAGS = (
    "parent",
    "source_pins",
    "graph",
    "equation",
    "scale",
    "core_entry",
    "ordinary_Z",
    "owner_definition",
    "raw_all_word_convergence",
    "passed",
)


class ProvenanceContractError(ValueError):
    """Raised only when a requested JSON mapping cannot be loaded."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not _is_mapping(value):
        raise ProvenanceContractError(f"{label} must be a JSON object mapping")
    return value


def _rooted_path(path: str | Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute() or candidate.is_file():
        return candidate
    module_root = Path(__file__).resolve().parents[2]
    rooted = module_root / candidate
    return rooted if rooted.is_file() else candidate


def _load_mapping(value: Any, default_path: Path, label: str) -> Mapping[str, Any]:
    """Load one JSON mapping, without interpreting any source file."""

    if _is_mapping(value):
        return value
    path = _rooted_path(default_path if value is None else value)
    try:
        loaded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError) as exc:
        raise ProvenanceContractError(f"{label} is not readable JSON: {path}") from exc
    return _mapping(loaded, label)


def _exact_mapping(value: Any) -> dict[str, Any] | None:
    return dict(value) if _is_mapping(value) else None


def _family(value: Any) -> dict[str, Any] | None:
    """Return a direct source-family object, never a transitive receipt."""

    if not _is_mapping(value):
        return None
    candidate = value.get("source_family")
    return _exact_mapping(candidate)


def _saved_family(v3: Mapping[str, Any], actual: Mapping[str, Any]) -> dict[str, Any] | None:
    # The archived v3 result has no family object.  These direct optional
    # locations are supported for future saved packets and make a synthetic
    # wrong-family mapping fail closed; the archived packet uses its exact
    # direct parent-report pin instead.
    for owner in (v3, actual):
        candidate = _family(owner)
        if candidate is not None:
            return candidate
        manifest = owner.get("source_manifest")
        candidate = _family(manifest)
        if candidate is not None:
            return candidate
    return None


def _source_family(parent: Mapping[str, Any]) -> dict[str, Any] | None:
    return _family(parent)


def _direct_basename_pin(parent: Mapping[str, Any], path: str, expected: str) -> bool:
    pins = parent.get("input_hashes")
    if not _is_mapping(pins):
        return False
    basename = path.rsplit("/", 1)[-1]
    # The contract is intentionally exact: only the direct basename entry is
    # admitted.  A recursive/transitive receipt or a merely similar path is
    # not provenance for this leaf.
    return pins.get(basename) == expected and isinstance(pins.get(basename), str)


def _direct_implementation_pin(parent: Mapping[str, Any]) -> bool:
    implementations = parent.get("implementation_hashes")
    if not _is_mapping(implementations):
        return False
    return implementations.get(FIRST_RAW_SOURCE) == FIRST_RAW_SHA256


def _parent_checks(parent: Mapping[str, Any]) -> dict[str, bool]:
    checks: dict[str, bool] = {
        "parent_checks_passed": (
            _is_mapping(parent.get("checks"))
            and parent["checks"].get("passed") is True
        ),
        "parent_input_hash_volterra_basename": _direct_basename_pin(
            parent, VOLTERRA_SOURCE, VOLTERRA_SHA256
        ),
        "parent_input_hash_regular_basename": _direct_basename_pin(
            parent, REGULAR_SOURCE, REGULAR_SHA256
        ),
        "parent_implementation_hash_first_raw": _direct_implementation_pin(parent),
    }
    first_raw = parent.get("actual_first_raw_n1_function")
    checks["parent_actual_first_raw_present"] = _is_mapping(first_raw)
    if not _is_mapping(first_raw):
        for key in PARENT_FIRST_RAW_FLAGS:
            checks[f"parent_{key}"] = False
        for key in PARENT_FIRST_RAW_NESTED_FLAGS:
            checks[f"parent_actual_first_raw_checks_{key}"] = False
        checks["parent_actual_first_raw_passed"] = False
        checks["parent_actual_first_raw_source_only"] = False
        checks["parent_actual_first_raw_normalized_Kchi_kernel_is_false"] = False
        checks["parent_actual_first_raw_status"] = False
        return checks

    checks["parent_actual_first_raw_passed"] = first_raw.get("passed") is True
    checks["parent_actual_first_raw_source_only"] = first_raw.get("source_only") is True
    checks["parent_actual_first_raw_normalized_Kchi_kernel_is_false"] = (
        first_raw.get("normalized_Kchi_kernel_substituted") is False
    )
    checks["parent_actual_first_raw_status"] = (
        first_raw.get("status") == "source_bound_exact_raw_volterra_function"
    )
    for key in PARENT_FIRST_RAW_FLAGS:
        checks[f"parent_{key}"] = first_raw.get(key) is True
    nested = first_raw.get("checks")
    for key in PARENT_FIRST_RAW_NESTED_FLAGS:
        checks[f"parent_actual_first_raw_checks_{key}"] = (
            _is_mapping(nested) and nested.get(key) is True
        )
    return checks


def _actual_checks(v3: Mapping[str, Any]) -> tuple[dict[str, bool], Mapping[str, Any] | None]:
    actual = v3.get("actual_result")
    checks: dict[str, bool] = {
        "v3_schema": v3.get("schema") == V3_SCHEMA,
        "v3_status_failed": v3.get("status") == "failed",
        "v3_actual_leaf_calls_one": v3.get("actual_leaf_calls") == 1,
        "v3_fixture_passed": (
            _is_mapping(v3.get("fixture")) and v3["fixture"].get("passed") is True
        ),
        "v3_actual_result_present": _is_mapping(actual),
    }
    if not _is_mapping(actual):
        checks.update(
            {
                "v3_actual_result_passed_false": False,
                "v3_source_executed_false": False,
                "v3_old_modules_not_called": False,
                "v3_numeric_evaluation_not_performed": False,
                "v3_input_hashes_mapping": False,
                "v3_input_hash_parent_report_exact": False,
                "v3_input_hash_volterra_exact": False,
                "v3_input_hash_regular_exact": False,
                "v3_input_hash_first_raw_exact": False,
                "v3_false_checks_exact": False,
                "v3_all_non_exempt_checks_true": False,
                "v3_source_hashes_exact": False,
                "v3_source_paths_exact": False,
                "v3_rows_exact": False,
                "v3_fields_exact": False,
                "v3_scale_diagonal_exact": False,
                "v3_parent_report_passed_checks_true": False,
                "v3_original_per_source_preserved": False,
                "v3_fixed_x_premise": False,
                "v3_fixed_hb_premise": False,
                "v3_function_owner_present": False,
                "v3_function_level_proof_identity_exact": False,
                "v3_function_level_proof_entrywise_argument_exact": False,
                "v3_function_level_proof_A_argument_exact": False,
                "v3_function_level_proof_green_kernel_Z_derivative_exact": False,
                "v3_scale_mapping_records_six_exact": False,
                "v3_scale_mapping_exact": False,
                "v3_function_level_proof_scale_rows_exact": False,
                "v3_function_level_proof_conditional_residual_exact": False,
                "v3_conditional_residual_owner_unestablished": False,
                "v3_ordinary_Z_connection_exact": False,
                "v3_ordinary_Z_Gamma_exact": False,
                "v3_ordinary_Z_B1_S_Z_retained": False,
                "v3_parent_raw_scale_binding_exact": False,
                "v3_limitations_required_false": False,
            }
        )
        return checks, None

    checks["v3_actual_result_passed_false"] = actual.get("passed") is False
    checks["v3_source_executed_false"] = actual.get("source_executed") is False
    checks["v3_old_modules_not_called"] = actual.get("old_modules_called") is False
    checks["v3_numeric_evaluation_not_performed"] = (
        actual.get("numeric_evaluation_performed") is False
    )
    v3_input_hashes = v3.get("input_hashes")
    checks["v3_input_hashes_mapping"] = _is_mapping(v3_input_hashes)
    checks["v3_input_hash_parent_report_exact"] = (
        _is_mapping(v3_input_hashes)
        and v3_input_hashes.get(str(PARENT_REPORT_PATH).replace("\\", "/"))
        == PARENT_REPORT_SHA256
    )
    checks["v3_input_hash_volterra_exact"] = (
        _is_mapping(v3_input_hashes)
        and v3_input_hashes.get(VOLTERRA_SOURCE) == VOLTERRA_SHA256
    )
    checks["v3_input_hash_regular_exact"] = (
        _is_mapping(v3_input_hashes)
        and v3_input_hashes.get(REGULAR_SOURCE) == REGULAR_SHA256
    )
    checks["v3_input_hash_first_raw_exact"] = (
        _is_mapping(v3_input_hashes)
        and v3_input_hashes.get(FIRST_RAW_SOURCE) == FIRST_RAW_SHA256
    )
    actual_checks = actual.get("checks")
    false_keys = {
        key for key, value in actual_checks.items()
        if value is False
    } if _is_mapping(actual_checks) else set()
    checks["v3_false_checks_exact"] = false_keys == set(V3_PARENT_HASH_EXEMPT_FALSE)
    checks["v3_all_non_exempt_checks_true"] = (
        _is_mapping(actual_checks)
        and all(
            value is True
            for key, value in actual_checks.items()
            if key not in V3_PARENT_HASH_EXEMPT_FALSE
        )
    )
    checks["v3_source_hashes_exact"] = actual.get("source_hashes") == EXPECTED_SOURCE_HASHES
    checks["v3_source_paths_exact"] = actual.get("source_paths") == EXPECTED_SOURCE_PATHS
    checks["v3_rows_exact"] = tuple(actual.get("rows", ())) == ROWS
    checks["v3_fields_exact"] = tuple(actual.get("fields", ())) == FIELDS
    checks["v3_scale_diagonal_exact"] = tuple(actual.get("scale_diagonal", ())) == SCALE_DIAGONAL
    parent_report = actual.get("parent_report")
    checks["v3_parent_report_passed_checks_true"] = (
        _is_mapping(parent_report) and parent_report.get("parent_checks_passed") is True
    )
    checks["v3_original_per_source_preserved"] = (
        _is_mapping(parent_report)
        and parent_report.get("per_source") == EXPECTED_PER_SOURCE
    )
    premises = actual.get("premises")
    checks["v3_fixed_x_premise"] = _is_mapping(premises) and premises.get("fixed_x") is True
    checks["v3_fixed_hb_premise"] = _is_mapping(premises) and premises.get("fixed_hb") is True
    owner = actual.get("function_owner")
    checks["v3_function_owner_present"] = (
        _is_mapping(owner)
        and owner.get("green_apply") == "<module>"
        and owner.get("expected") == "<module>"
    )
    proof = actual.get("function_level_proof")
    checks["v3_function_level_proof_identity_exact"] = (
        _is_mapping(proof) and proof.get("identity") == EXPECTED_PROOF_IDENTITY
    )
    checks["v3_function_level_proof_entrywise_argument_exact"] = (
        _is_mapping(proof)
        and proof.get("entrywise_argument") == EXPECTED_ENTRYWISE_ARGUMENT
    )
    checks["v3_function_level_proof_A_argument_exact"] = (
        _is_mapping(proof) and proof.get("A_argument") == EXPECTED_A_ARGUMENT
    )
    checks["v3_function_level_proof_green_kernel_Z_derivative_exact"] = (
        _is_mapping(proof)
        and proof.get("green_kernel_Z_derivative") == EXPECTED_GREEN_KERNEL_Z_DERIVATIVE
    )
    scale_mapping = actual.get("scale_mapping")
    records = scale_mapping.get("records") if _is_mapping(scale_mapping) else None
    checks["v3_scale_mapping_records_six_exact"] = (
        isinstance(records, list)
        and len(records) == 6
        and records == EXPECTED_SCALE_RECORDS
    )
    checks["v3_scale_mapping_exact"] = scale_mapping == EXPECTED_SCALE_MAPPING
    checks["v3_function_level_proof_scale_rows_exact"] = (
        _is_mapping(proof) and proof.get("scale_rows") == EXPECTED_SCALE_MAPPING
    )
    conditional_residual = (
        proof.get("conditional_raw_residual_transport")
        if _is_mapping(proof)
        else None
    )
    checks["v3_function_level_proof_conditional_residual_exact"] = (
        conditional_residual == EXPECTED_CONDITIONAL_RESIDUAL
    )
    checks["v3_conditional_residual_owner_unestablished"] = (
        _is_mapping(conditional_residual)
        and conditional_residual.get("source_owner_established_here") is False
    )
    ordinary_z_connection = proof.get("ordinary_Z_connection") if _is_mapping(proof) else None
    checks["v3_ordinary_Z_connection_exact"] = (
        ordinary_z_connection == EXPECTED_ORDINARY_Z_CONNECTION
    )
    checks["v3_ordinary_Z_Gamma_exact"] = (
        _is_mapping(ordinary_z_connection)
        and ordinary_z_connection.get("Gamma")
        == EXPECTED_ORDINARY_Z_CONNECTION["Gamma"]
    )
    checks["v3_ordinary_Z_B1_S_Z_retained"] = (
        _is_mapping(ordinary_z_connection)
        and ordinary_z_connection.get("B1_S_Z_retained") is True
    )
    parent_raw_scale_binding = proof.get("parent_raw_scale_binding") if _is_mapping(proof) else None
    checks["v3_parent_raw_scale_binding_exact"] = (
        parent_raw_scale_binding == EXPECTED_PARENT_RAW_SCALE_BINDING
    )
    limitations = actual.get("limitations")
    checks["v3_limitations_required_false"] = (
        _is_mapping(limitations)
        and all(limitations.get(key) is False for key in PROOF_LIMITATION_FALSE_FLAGS)
    )
    return checks, actual


def _family_check(
    v3: Mapping[str, Any], actual: Mapping[str, Any], parent: Mapping[str, Any]
) -> tuple[bool, str]:
    v3_input_hashes = v3.get("input_hashes")
    direct_parent_report_link = (
        _is_mapping(v3_input_hashes)
        and v3_input_hashes.get(str(PARENT_REPORT_PATH).replace("\\", "/"))
        == PARENT_REPORT_SHA256
    )
    saved_family = _saved_family(v3, actual)
    if saved_family is not None:
        parent_family = _source_family(parent)
        return (
            direct_parent_report_link
            and parent_family is not None
            and saved_family == parent_family,
            "same pinned parent report, with optional explicit family comparison",
        )
    # The archived v3 packet has no family object.  Its direct input_hashes
    # pin the exact cutoff report and all three source bytes, which is the
    # admissible family identity rather than an inferred label.
    return (
        direct_parent_report_link,
        "same pinned parent report, source family inherited from that report",
    )


def audit_saved_mappings(
    v3_result: Mapping[str, Any],
    parent_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Audit saved mappings without reading or parsing any source text."""

    v3 = _mapping(v3_result, "v3_result")
    parent = _mapping(parent_report, "parent_report")
    checks, actual = _actual_checks(v3)
    parent_checks = _parent_checks(parent)
    checks.update(parent_checks)
    if actual is None:
        checks["source_family_same"] = False
        family_basis = "actual_result is missing"
    else:
        checks["source_family_same"], family_basis = _family_check(v3, actual, parent)

    passed = bool(checks) and all(value is True for value in checks.values())
    actual_facts = {
        "source_hashes": dict(actual.get("source_hashes", {})) if actual else {},
        "source_paths": dict(actual.get("source_paths", {})) if actual else {},
        "rows": list(actual.get("rows", ())) if actual else [],
        "fields": list(actual.get("fields", ())) if actual else [],
        "scale_diagonal": list(actual.get("scale_diagonal", ())) if actual else [],
        "function_owner": dict(actual.get("function_owner", {})) if actual and _is_mapping(actual.get("function_owner")) else {},
        "scale_mapping": dict(actual.get("scale_mapping", {})) if actual and _is_mapping(actual.get("scale_mapping")) else {},
        "function_level_proof": dict(actual.get("function_level_proof", {})) if actual and _is_mapping(actual.get("function_level_proof")) else {},
        "limitations": dict(actual.get("limitations", {})) if actual and _is_mapping(actual.get("limitations")) else {},
    }
    original_per_source = {}
    if actual and _is_mapping(actual.get("parent_report")):
        original_per_source = dict(actual["parent_report"].get("per_source", {}))
    actual_source_hashes = actual.get("source_hashes") if actual else None
    saved_input_hashes = v3.get("input_hashes")
    corrected_per_source = {
        VOLTERRA_SOURCE: (
            checks.get("parent_input_hash_volterra_basename") is True
            and _is_mapping(actual_source_hashes)
            and actual_source_hashes.get(VOLTERRA_SOURCE) == VOLTERRA_SHA256
            and _is_mapping(saved_input_hashes)
            and saved_input_hashes.get(VOLTERRA_SOURCE) == VOLTERRA_SHA256
        ),
        REGULAR_SOURCE: (
            checks.get("parent_input_hash_regular_basename") is True
            and _is_mapping(actual_source_hashes)
            and actual_source_hashes.get(REGULAR_SOURCE) == REGULAR_SHA256
            and _is_mapping(saved_input_hashes)
            and saved_input_hashes.get(REGULAR_SOURCE) == REGULAR_SHA256
        ),
        FIRST_RAW_SOURCE: (
            checks.get("parent_implementation_hash_first_raw") is True
            and _is_mapping(actual_source_hashes)
            and actual_source_hashes.get(FIRST_RAW_SOURCE) == FIRST_RAW_SHA256
            and _is_mapping(saved_input_hashes)
            and saved_input_hashes.get(FIRST_RAW_SOURCE) == FIRST_RAW_SHA256
        ),
    }
    return {
        "schema": SCHEMA,
        "passed": passed,
        "checks": checks,
        "v3_status": v3.get("status"),
        "v3_actual_leaf_calls": v3.get("actual_leaf_calls"),
        "v3_original_per_source": original_per_source,
        "corrected_per_source": corrected_per_source,
        "family_basis": family_basis,
        "actual_result_source_facts": actual_facts,
        "actual_result_checks": (
            dict(actual.get("checks", {}))
            if actual and _is_mapping(actual.get("checks"))
            else {}
        ),
        "source_AST_recomputed": False,
        "source_text_consumed": False,
        "old_modules_called": False,
    }


def _acceptance(
    audit: Mapping[str, Any],
    v3: Mapping[str, Any],
    parent: Mapping[str, Any],
) -> dict[str, Any]:
    actual = v3["actual_result"]
    facts = audit["actual_result_source_facts"]
    parent_family = (
        dict(parent["source_family"])
        if _is_mapping(parent.get("source_family"))
        else {}
    )
    return {
        "kind": "corrected_provenance_acceptance",
        "basis": "same pinned parent report, source family inherited from that report",
        "source_facts": facts,
        "original_v3_status": v3.get("status"),
        "original_v3_actual_result_passed": actual.get("passed"),
        "original_v3_checks_false_only": sorted(V3_PARENT_HASH_EXEMPT_FALSE),
        "corrected_parent_pins": {
            "v3_input_hashes": {
                str(PARENT_REPORT_PATH).replace("\\", "/"): PARENT_REPORT_SHA256,
                VOLTERRA_SOURCE: VOLTERRA_SHA256,
                REGULAR_SOURCE: REGULAR_SHA256,
                FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
            },
            "input_hashes": {
                VOLTERRA_SOURCE.rsplit("/", 1)[-1]: VOLTERRA_SHA256,
                REGULAR_SOURCE.rsplit("/", 1)[-1]: REGULAR_SHA256,
            },
            "implementation_hashes": {FIRST_RAW_SOURCE: FIRST_RAW_SHA256},
        },
        "source_family": parent_family,
        "parent_false_reason": (
            "The archived parent admission checked first_raw against input_hashes; "
            "the authoritative direct first_raw pin is in implementation_hashes. "
            "No source or AST check failed."
        ),
        "conditional_identity": "S^-1*G*S = G",
        "conditional_residual_transport": (
            "If the same owner has Qhat=S^-1*Q*S, retaining B1*S_Z, "
            "then at fixed x d_hat=S^-1*d_raw and "
            "d_hat_Z=S^-1*d_raw_Z-Gamma*d_hat."
        ),
        "premises": {
            "fixed_x": True,
            "fixed_hb": True,
            "endpoint_independent_of_Z": True,
            "scale_entries_nonzero": True,
            "source_structure_reused_from_saved_v3": True,
        },
        "limitations": {
            "live_Kchi_owner_proved": False,
            "complete_Kchi_resolvent_proved": False,
            "raw_numerical_residual_evaluated": False,
            "value_or_ledger_claim": False,
            "source_reexecuted": False,
            "source_AST_recomputed": False,
            "historical_error_removed": False,
            "historical_error_added": False,
        },
    }


def build_provenance_bridge(
    v3_result: Mapping[str, Any] | str | Path | None = None,
    parent_report: Mapping[str, Any] | str | Path | None = None,
    *,
    result_path: str | Path = V3_RESULT_PATH,
    parent_path: str | Path = PARENT_REPORT_PATH,
) -> dict[str, Any]:
    """Return corrected acceptance from the saved v3 and parent mappings."""

    v3 = _load_mapping(v3_result, Path(result_path), "v3_result")
    parent = _load_mapping(parent_report, Path(parent_path), "parent_report")
    audit = audit_saved_mappings(v3, parent)
    result: dict[str, Any] = dict(audit)
    result["schema"] = SCHEMA
    result["status"] = "accepted" if audit["passed"] else "rejected"
    result["parent_report_checks_passed"] = (
        _is_mapping(parent.get("checks")) and parent["checks"].get("passed") is True
    )
    result["parent_report_source_family"] = (
        dict(parent["source_family"]) if _is_mapping(parent.get("source_family")) else {}
    )
    if audit["passed"]:
        result["acceptance"] = _acceptance(audit, v3, parent)
    else:
        result["acceptance"] = None
    return result


def accept_saved_result(*args: Any, **kwargs: Any) -> dict[str, Any]:
    return build_provenance_bridge(*args, **kwargs)


audit_saved_result = audit_saved_mappings
build = build_provenance_bridge
run = build_provenance_bridge
prove = build_provenance_bridge
verify = build_provenance_bridge


__all__ = [
    "SCHEMA",
    "V3_SCHEMA",
    "V3_RESULT_PATH",
    "PARENT_REPORT_PATH",
    "PARENT_REPORT_SHA256",
    "VOLTERRA_SOURCE",
    "REGULAR_SOURCE",
    "FIRST_RAW_SOURCE",
    "VOLTERRA_SHA256",
    "REGULAR_SHA256",
    "FIRST_RAW_SHA256",
    "EXPECTED_V3_INPUT_HASHES",
    "ROWS",
    "FIELDS",
    "SCALE_DIAGONAL",
    "V3_PARENT_HASH_EXEMPT_FALSE",
    "ProvenanceContractError",
    "audit_saved_mappings",
    "audit_saved_result",
    "build_provenance_bridge",
    "accept_saved_result",
    "build",
    "run",
    "prove",
    "verify",
]
