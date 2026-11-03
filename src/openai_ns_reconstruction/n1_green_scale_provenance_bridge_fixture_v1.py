"""Independent mapping-only fixture for the provenance bridge v1.

The fixture never opens the archived v3 result, the cutoff report, or any
source file.  It constructs a small synthetic mapping contract and exercises
the new bridge auditor with one accepted packet and six distinct rejection
packets.  The caller supplies the bridge auditor so this file does not import
or re-enter an earlier leaf.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any


SCHEMA = "n1_green_scale_provenance_bridge_fixture_v1"
V3_SCHEMA = "n1_green_scale_commutator_v3_result_v1"

VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
FIRST_RAW_SOURCE = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"
PARENT_REPORT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
)
VOLTERRA_BASENAME = VOLTERRA_SOURCE.rsplit("/", 1)[-1]
REGULAR_BASENAME = REGULAR_SOURCE.rsplit("/", 1)[-1]

VOLTERRA_SHA256 = "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
REGULAR_SHA256 = "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555"
FIRST_RAW_SHA256 = "f13f29b4d0652e8ed3b3347f19bab7eb1ea0e18fd4be0f50cf6d550a1aeb92c6"
PARENT_REPORT_SHA256 = "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856"

ROWS = ["5", "6", "6", "4", "5", "6"]
FIELDS = ["F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1"]
SCALE_DIAGONAL = [
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
]
SCALE_RECORDS = [
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
SCALE_MAPPING = {
    "rows": list(ROWS),
    "records": SCALE_RECORDS,
    "all_output_input_scales_equal": True,
}
CONDITIONAL_RESIDUAL = {
    "premise_same_function_owner": "the raw and hatted six-component functions are the same owner",
    "premise_Qhat": "Qhat=S^-1*Q*S, with B1*S_Z retained in the ordinary-Z connection",
    "value_identity": "d_hat=S^-1*d_raw",
    "fixed_x_ordinary_Z_identity": "d_hat_Z=S^-1*d_raw_Z-Gamma*d_hat, Gamma=S^-1*S_Z",
    "source_owner_established_here": False,
}
ORDINARY_Z_CONNECTION = {
    "S_Z": "diag(F0_Z/epsilon,0,0,0,F0_Z/epsilon,0)",
    "Gamma": "S^-1*S_Z=diag(F0_Z/F0,0,0,0,F0_Z/F0,0)",
    "B0_transformed": "S^-1*(B0*S+B1*S_Z)",
    "B1_transformed": "S^-1*B1*S",
    "g_transformed": "S^-1*g",
    "B1_S_Z_retained": True,
    "A_and_B1_source_reverified_here": False,
}
PARENT_RAW_SCALE_BINDING = {
    "parent_actual_first_raw_scale_conjugacy": True,
    "parent_raw_preconditioned_resolvent_factorization": True,
    "parent_normalized_Kchi_kernel_substituted_is_false": True,
    "source_AST_scale_rows_match": True,
    "compatibility_check": True,
}
LIMITATIONS = {
    "complete_Kchi_resolvent_provided": False,
    "raw_numerical_residual_evaluated": False,
    "gauge_owner_proved": False,
    "historical_error_removed": False,
    "historical_error_added": False,
    "interval_or_tuple_overlap_used_as_function_identity": False,
    "source_provider_executed": False,
    "full_B0_B1_operator_conjugacy_beyond_the_proved_connection": False,
    "A_and_B1_support_reverified_here": False,
}

FAMILY = {
    "actual_five_defect_family_sha256": "synthetic-family-a-five-defect",
    "implicit_source_sha256": "synthetic-family-a-implicit-source",
    "datum_enclosure_sha256": "synthetic-family-a-datum-enclosure",
}
OTHER_FAMILY = {
    "actual_five_defect_family_sha256": "synthetic-family-b-five-defect",
    "implicit_source_sha256": "synthetic-family-b-implicit-source",
    "datum_enclosure_sha256": "synthetic-family-b-datum-enclosure",
}


def _parent_mapping(family: Mapping[str, str]) -> dict[str, Any]:
    nested_checks = {
        key: True
        for key in (
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
    }
    first_raw: dict[str, Any] = {
        "passed": True,
        "source_only": True,
        "normalized_Kchi_kernel_substituted": False,
        "status": "source_bound_exact_raw_volterra_function",
        "checks": nested_checks,
    }
    for key in (
        "actual_first_raw_source_definition_bound",
        "actual_first_raw_function_bound",
        "actual_first_raw_equation_bound",
        "actual_first_raw_core_entry_bound",
        "actual_first_raw_ordinary_Z0_to_3_bound",
        "actual_first_raw_scale_conjugacy",
        "actual_first_raw_all_word_convergence_bound",
        "raw_all_word_convergence_certified",
        "raw_preconditioned_resolvent_factorization",
    ):
        first_raw[key] = True
    return {
        "source_family": dict(family),
        "checks": {"passed": True},
        "input_hashes": {
            VOLTERRA_BASENAME: VOLTERRA_SHA256,
            REGULAR_BASENAME: REGULAR_SHA256,
        },
        "implementation_hashes": {FIRST_RAW_SOURCE: FIRST_RAW_SHA256},
        "actual_first_raw_n1_function": first_raw,
    }


def _v3_mapping(family: Mapping[str, str]) -> dict[str, Any]:
    actual_checks = {
        "saved_source_and_ast_structure": True,
        "source_pin_volterra": True,
        "source_pin_regular": True,
        "source_pin_first_raw": True,
        "parent_report_checks_passed": True,
        "scale_mapping_by_literal_input_row": True,
        "parent_report_input_hashes_present": False,
        "parent_report_input_hashes_match": False,
        "parent_report_passed": False,
    }
    actual = {
        "passed": False,
        "checks": actual_checks,
        "source_executed": False,
        "old_modules_called": False,
        "numeric_evaluation_performed": False,
        "source_hashes": {
            VOLTERRA_SOURCE: VOLTERRA_SHA256,
            REGULAR_SOURCE: REGULAR_SHA256,
            FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
        },
        "source_paths": {
            "volterra": VOLTERRA_SOURCE,
            "regular": REGULAR_SOURCE,
            "first_raw": FIRST_RAW_SOURCE,
        },
        "rows": list(ROWS),
        "fields": list(FIELDS),
        "scale_diagonal": list(SCALE_DIAGONAL),
        "function_owner": {"green_apply": "<module>", "expected": "<module>"},
        "premises": {"fixed_x": True, "fixed_hb": True},
        "scale_mapping": deepcopy(SCALE_MAPPING),
        "function_level_proof": {
            "A_argument": (
                "conditional algebraic reuse: if A=-2*chi*E51, then S1=S5=F0(Z)/epsilon "
                "gives S^-1*A*S=A; the external A support is not reverified by this leaf"
            ),
            "entrywise_argument": (
                "G[i,j] != 0 only for j=ROWS[i]-1, and SCALE_DIAGONAL[i]=SCALE_DIAGONAL[j]; "
                "hence S_i^-1 G[i,j] S_j=G[i,j]"
            ),
            "identity": "S^-1*G*S = G",
            "green_kernel_Z_derivative": (
                "G_Z=0 under fixed x, fixed hb, and Z-independent endpoints/phase bounds"
            ),
            "conditional_raw_residual_transport": deepcopy(CONDITIONAL_RESIDUAL),
            "ordinary_Z_connection": deepcopy(ORDINARY_Z_CONNECTION),
            "parent_raw_scale_binding": deepcopy(PARENT_RAW_SCALE_BINDING),
            "scale_rows": deepcopy(SCALE_MAPPING),
        },
        "limitations": deepcopy(LIMITATIONS),
        "parent_report": {
            "parent_checks_passed": True,
            "per_source": {
                VOLTERRA_SOURCE: True,
                REGULAR_SOURCE: True,
                FIRST_RAW_SOURCE: False,
            },
        },
    }
    return {
        "schema": V3_SCHEMA,
        "status": "failed",
        "actual_leaf_calls": 1,
        "fixture": {"passed": True},
        "actual_result": actual,
        "input_hashes": {
            PARENT_REPORT_SOURCE: PARENT_REPORT_SHA256,
            VOLTERRA_SOURCE: VOLTERRA_SHA256,
            REGULAR_SOURCE: REGULAR_SHA256,
            FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
        },
        "source_family": dict(family),
    }


def _auditor(auditor: Any) -> Any:
    candidate = getattr(auditor, "audit_saved_mappings", None)
    if not callable(candidate) and callable(auditor):
        candidate = auditor
    if not callable(candidate):
        raise TypeError("auditor must expose audit_saved_mappings or be callable")
    return candidate


def _invoke(auditor: Any, v3: Mapping[str, Any], parent: Mapping[str, Any]) -> Mapping[str, Any] | None:
    try:
        result = _auditor(auditor)(v3, parent)
    except (TypeError, ValueError, KeyError, OSError):
        return None
    return result if isinstance(result, Mapping) else None


def _false_keys(result: Mapping[str, Any] | None) -> set[str]:
    if not isinstance(result, Mapping) or not isinstance(result.get("checks"), Mapping):
        return {"<no_checks>"}
    return {
        str(key)
        for key, value in result["checks"].items()
        if value is False
    }


def _rejected_for(result: Mapping[str, Any] | None, expected: set[str]) -> bool:
    return (
        isinstance(result, Mapping)
        and result.get("passed") is False
        and expected.issubset(_false_keys(result))
    )


def run_fixture(
    auditor: Any,
    repo: str | Path | None = None,
    parent: Mapping[str, Any] | str | Path | None = None,
    source_texts: Any = None,
    **_: Any,
) -> dict[str, Any]:
    """Run mapping-only positive and negative provenance cases."""

    del repo, parent, source_texts
    candidate = _auditor(auditor)
    good_v3 = _v3_mapping(FAMILY)
    good_parent = _parent_mapping(FAMILY)
    accepted = _invoke(candidate, good_v3, good_parent)

    wrong_implementation_parent = deepcopy(good_parent)
    wrong_implementation_parent["implementation_hashes"][FIRST_RAW_SOURCE] = "wrong-implementation-hash"
    wrong_implementation = _invoke(candidate, good_v3, wrong_implementation_parent)

    wrong_input_parent = deepcopy(good_parent)
    wrong_input_parent["input_hashes"][REGULAR_BASENAME] = "wrong-input-hash"
    wrong_input = _invoke(candidate, good_v3, wrong_input_parent)

    extra_error_v3 = deepcopy(good_v3)
    extra_error_v3["actual_result"]["checks"]["saved_unrelated_error"] = False
    extra_error = _invoke(candidate, extra_error_v3, good_parent)

    zero_calls_v3 = deepcopy(good_v3)
    zero_calls_v3["actual_leaf_calls"] = 0
    zero_calls = _invoke(candidate, zero_calls_v3, good_parent)

    other_family_v3 = _v3_mapping(OTHER_FAMILY)
    other_family_v3["input_hashes"][PARENT_REPORT_SOURCE] = "wrong-parent-report-pin"
    other_family = _invoke(candidate, other_family_v3, good_parent)

    wrong_identity_v3 = deepcopy(good_v3)
    wrong_identity_v3["actual_result"]["function_level_proof"]["identity"] = "wrong-proof-identity"
    wrong_identity = _invoke(candidate, wrong_identity_v3, good_parent)

    checks = {
        "correct_mapping_accepted": (
            isinstance(accepted, Mapping) and accepted.get("passed") is True
        ),
        "wrong_implementation_hash_rejected": _rejected_for(
            wrong_implementation, {"parent_implementation_hash_first_raw"}
        ),
        "wrong_input_hash_rejected": _rejected_for(
            wrong_input, {"parent_input_hash_regular_basename"}
        ),
        "extra_error_rejected": _rejected_for(
            extra_error,
            {"v3_false_checks_exact", "v3_all_non_exempt_checks_true"},
        ),
        "zero_actual_leaf_calls_rejected": _rejected_for(
            zero_calls, {"v3_actual_leaf_calls_one"}
        ),
        "different_family_rejected": _rejected_for(
            other_family,
            {"source_family_same", "v3_input_hash_parent_report_exact"},
        ),
        "wrong_proof_identity_rejected": _rejected_for(
            wrong_identity, {"v3_function_level_proof_identity_exact"}
        ),
    }
    false_by_case = {
        "wrong_implementation_hash": _false_keys(wrong_implementation),
        "wrong_input_hash": _false_keys(wrong_input),
        "extra_error": _false_keys(extra_error),
        "zero_actual_leaf_calls": _false_keys(zero_calls),
        "different_family": _false_keys(other_family),
        "wrong_proof_identity": _false_keys(wrong_identity),
    }
    synthetic_records = {
        "correct": "new_mapping_only_corrected_parent_pins",
        "wrong_implementation_hash": "new_mapping_only_first_raw_implementation_pin",
        "wrong_input_hash": "new_mapping_only_regular_input_pin",
        "extra_error": "new_mapping_only_unrelated_v3_false_check",
        "zero_actual_leaf_calls": "new_mapping_only_v3_call_count",
        "different_family": "new_mapping_only_family_identity",
        "wrong_proof_identity": "new_mapping_only_saved_function_proof_identity",
    }
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_texts_consumed": False,
        "source_AST_recomputed": False,
        "actual_source_called_here": False,
        "v3_actual_leaf_called_here": False,
        "v1_or_v2_fixture_called": False,
        "old_modules_called": False,
        "provider_executed": False,
        "numeric_residual_evaluated": False,
        "implementation_imported": False,
        "synthetic_record_count": len(synthetic_records),
        "false_checks_by_case": {
            key: sorted(value) for key, value in false_by_case.items()
        },
        "synthetic_records": synthetic_records,
    }


def run(
    auditor: Any,
    repo: str | Path | None = None,
    parent: Mapping[str, Any] | str | Path | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    return run_fixture(auditor, repo=repo, parent=parent, **kwargs)


__all__ = ["SCHEMA", "run_fixture", "run"]
