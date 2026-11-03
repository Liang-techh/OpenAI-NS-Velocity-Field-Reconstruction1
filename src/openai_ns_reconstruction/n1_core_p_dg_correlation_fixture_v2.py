"""Independent exact-rational fixture for the core P/Dg correlation leaf.

The fixture builds a small verified source-functional reference and then
removes one proof at a time. It never imports a producer or the leaf under
test. Negative cases must close with ``status == 'conditional'`` and retain
the saved Dg rows; they must not manufacture a P alias from source-family
labels or interval overlap.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "n1_core_p_dg_correlation_fixture_v2"
RECEIPT_PATHS = {
    "KCHI": "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz",
    "FEEDBACK": "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz",
    "COMP": "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz",
    "LOW": "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz",
    "LOW3": "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz",
}


def _value(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _interval(c: Any, lower: Fraction | int, upper: Fraction | int | None = None) -> Any:
    if upper is None:
        upper = lower
    lo = _value(c, lower)
    hi = _value(c, upper)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))
        except Exception:
            pass
    return c.mpf([lo, hi])


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_clone(item) for item in value)
    return value


def _term(
    c: Any,
    *,
    parent_id: str,
    z_label: str,
    index: int,
    lower: Fraction,
    upper: Fraction,
    powers: tuple[int, int, int, int],
    error_only: bool = False,
) -> dict[str, Any]:
    source_id = f"v2-g5:{z_label}:{index}"
    source_path = f"receipt://v2/Dg/J1/{z_label}/{index}"
    coefficient = _interval(c, lower, upper)
    return {
        "parent_id": parent_id,
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": list(powers),
        "formal_powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [source_path],
        "error_classes": ["fixture_v2_saved_physical_Dg"],
        "uncertainty": error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
    }


def _packet(c: Any) -> dict[str, Any]:
    return {
        "physical_level": "m15-v2",
        "J1": {
            "Z0": [
                _term(c, parent_id="v2-dg-z0-a", z_label="Z0", index=0, lower=Fraction(2, 7), upper=Fraction(5, 7), powers=(2, -1, 4, 1)),
                _term(c, parent_id="v2-dg-z0-b", z_label="Z0", index=1, lower=Fraction(-11, 4), upper=Fraction(-3, 2), powers=(3, 1, -2, 5)),
            ],
            "Z1": [
                _term(c, parent_id="v2-dg-z1-a", z_label="Z1", index=0, lower=Fraction(-9, 5), upper=Fraction(7, 6), powers=(-2, 4, 1, -1)),
                _term(c, parent_id="v2-dg-z1-b", z_label="Z1", index=1, lower=Fraction(-2, 9), upper=Fraction(2, 9), powers=(1, 2, 3, -4), error_only=True),
            ],
        },
    }


def _receipt(path: str, index: int) -> dict[str, Any]:
    return {
        "path": path,
        "raw_hash_authenticated": True,
        "decoded": True,
        "raw_sha256": f"{(index + 40):02x}" * 32,
        "status": "returned",
        "source_manifest": {
            "sha256": f"{(index + 50):02x}" * 32,
            "files": [{"path": f"fixture/v2/source/record-{index}", "sha256": f"{(index + 60):02x}" * 32}],
        },
        "decoded_output": {"status": "returned", "source_family": "fixture-family-v2"},
    }


def _evidence(packet: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    receipts = {label: _receipt(path, index) for index, (label, path) in enumerate(RECEIPT_PATHS.items(), 1)}
    owners: dict[str, Any] = {}
    for z_label in ("Z0", "Z1"):
        for index, term in enumerate(packet["J1"][z_label]):
            owners[term["parent_id"]] = {
                "parent_id": term["parent_id"], "owner": "Dg_physical", "passed": True, "proved": True,
                "source_ids": list(term["source_ids"]), "source_paths": list(term["source_paths"]),
                "parent_hash": f"v2-parent-{z_label}-{index}",
                "geometry": "actual-core-v2-first-collar-chart",
                "labels": {"family": "fixture-family-v2", "physical_level": "m15-v2", "ordinary_Z": z_label},
                "exact_cell_path": f"v2.result.output.physical_Dg.J1.{z_label}.{index}",
                "output_path": f"v2.result.output.projection.J1.{z_label}.{index}",
            }
    certificate = {
        "passed": True, "source_receipts": receipts, "owner_by_parent": owners,
        "Dg_source_refs": {"receipt": RECEIPT_PATHS["FEEDBACK"], "output": "v2.result.output.m15.Dg"},
        "P_source_refs": {"receipt": RECEIPT_PATHS["COMP"], "output": "v2.result.output.potential_compensation_m15"},
        "same_actual_g5": True, "same_chi_Z": True, "same_green_K0": True,
        "same_complete_0_2_kernel": True, "same_physical_m15_measure": True,
        "same_leading_F": True, "low3_W0_same_K0_g5": True,
        "kchi_Dg_minus2_chi_K0_K0g5": True, "ordinary_Z1_chi_Z_retained": True,
        "ordinary_Z1_g5_Z_retained": True, "same_geometry": True, "same_labels": True,
        "same_parent_hashes": True, "exact_source_id_paths": True,
        "exact_cell_paths": True, "exact_outputs": True, "P_sign": 2,
        "physical_m15_functional_scope": "whole function only",
    }
    return certificate, receipts


def _call(build: Any, c: Any, packet: Any, certificate: Any, receipts: Any) -> Mapping[str, Any]:
    function: Callable[..., Any] = getattr(build, "build", build)
    result = function(c, packet, certificate, receipts)
    if not isinstance(result, Mapping):
        raise TypeError("leaf build did not return a mapping")
    return result


def _negative_case(build: Any, c: Any, packet: Any, certificate: Any, receipts: Any, key: str, value: Any) -> bool:
    bad = _clone(certificate)
    bad[key] = value
    try:
        result = _call(build, c, packet, bad, bad.get("source_receipts", receipts))
    except Exception:
        return False
    return (
        result.get("status") == "conditional" and result.get("passed") is False
        and result.get("correlation_proved") is False and result.get("pair_count") == 0
        and len(result.get("Dg", {}).get("Z0", ())) == 2
        and len(result.get("Dg", {}).get("Z1", ())) == 2
    )


def _functional_reference() -> dict[str, bool]:
    # A nontrivial rational K0 kernel verifies both ordinary routes without
    # adding any interval boxes.
    def kernel(q: Fraction) -> Fraction:
        return q * (1 - q / 3)

    chi0, chi1 = Fraction(3, 8), Fraction(-5, 9)
    g0, g1 = Fraction(7, 6), Fraction(5, 11)
    q0, q1 = Fraction(2, 5), Fraction(3, 5)
    d0 = -2 * kernel(q0) * chi0 * g0
    p0 = +2 * kernel(q0) * chi0 * g0
    d1 = -2 * kernel(q1) * (chi1 * g0 + chi0 * g1)
    p1 = +2 * kernel(q1) * (chi1 * g0 + chi0 * g1)
    return {
        "nontrivial_kernel": kernel(q0) != kernel(q1),
        "Z0_functional_negative": d0 + p0 == 0,
        "Z1_chi_Z_and_g5_Z_route": d1 + p1 == 0,
    }


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    packet = _packet(c)
    certificate, receipts = _evidence(packet)
    checks: dict[str, bool] = {}
    try:
        result = _call(build, c, packet, certificate, receipts)
    except Exception as exc:
        return {
            "schema": SCHEMA, "passed": False,
            "checks": {"build_completed": False}, "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}", "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "n1_core_p_dg_correlation_v2"
    checks["conditional_status"] = result.get("status") == "conditional_on_caller_assertions" and result.get("passed") is False
    checks["caller_assertions_recorded"] = (
        result.get("caller_certificate_verified") is False
        and result.get("caller_certificate_schema_satisfied") is True
        and result.get("caller_assertions_satisfied") is True
        and result.get("conditional_on_caller_assertions") is True
    )
    checks["actual_source_not_claimed"] = (
        result.get("actual_source_verified") is False
        and result.get("correlation_proved") is False
        and result.get("identity_verified") is False
        and result.get("zero_claimed") is False
        and result.get("residual_exact_zero") is False
    )
    checks["pair_count"] = result.get("pair_count") == 0 and result.get("functional_pair_count") == 1
    checks["ledger_untouched"] = result.get("ledger_written") is False and result.get("no_ledger") is True
    checks["P1_not_m15"] = result.get("P1_is_m15") is False and result.get("P1_not_m15") is True
    checks["no_interval_arithmetic"] = (
        result.get("coefficient_arithmetic_performed") is False
        and result.get("linewise_cancellation_performed") is False
        and result.get("independent_boxed_interval_sums") is False
    )
    checks["receipt_paths_only"] = result.get("receipt_paths_only") is True

    pairs = result.get("functional_pairs", [])
    checks["one_whole_function_pair"] = len(pairs) == 1
    functional = pairs[0] if pairs else {}
    expression = functional.get("functional_expression", {})
    checks["common_g5_expression"] = expression.get("W0_F") == "K0(g5)"
    checks["correct_kernel_sign_expression"] = (
        expression.get("Dg_F") == "-2 K0(chi K0(g5))"
        and expression.get("P_W0_F") == "+2 K0(chi K0(g5))"
    )
    checks["axis_chi_form"] = "-2 chi K0(K0(g5))" in expression.get("axis_chi_constant_form", "")
    checks["ordinary_Z1_expression"] = "chi_Z" in expression.get("ordinary_Z1", "") and "g5_Z" in expression.get("ordinary_Z1", "")
    checks["separate_source_refs"] = bool(functional.get("Dg_source_refs")) and bool(functional.get("P_source_refs"))
    checks["conditional_symbolic_implication"] = (
        functional.get("symbolic_identity_under_assumptions") is True
        and functional.get("actual_source_verified") is False
        and functional.get("identity_verified") is False
        and functional.get("zero_claimed") is False
    )
    checks["no_sector_P_alias"] = functional.get("independent_P_sector_materialized") is False and not result.get("P", {}).get("Z0") and not result.get("P", {}).get("Z1")
    checks["error_only_preserved"] = result.get("Dg", {}).get("Z1", [])[1].get("error_only") is True and result.get("Dg", {}).get("Z1", [])[1].get("uncertainty") is True
    checks.update({f"functional_{key}": value for key, value in _functional_reference().items()})

    checks["different_source_rejected"] = _negative_case(build, c, packet, certificate, receipts, "same_actual_g5", False)
    checks["different_measure_rejected"] = _negative_case(build, c, packet, certificate, receipts, "same_physical_m15_measure", False)
    checks["different_chi_rejected"] = _negative_case(build, c, packet, certificate, receipts, "same_chi_Z", False)
    checks["wrong_P_sign_rejected"] = _negative_case(build, c, packet, certificate, receipts, "P_sign", -2)
    checks["ordinary_Z1_missing_rejected"] = _negative_case(build, c, packet, certificate, receipts, "ordinary_Z1_g5_Z_retained", False)
    checks["coarse_LOW3_rejected"] = _negative_case(build, c, packet, certificate, receipts, "low3_W0_same_K0_g5", False)

    missing_owner = _clone(certificate)
    first_parent = packet["J1"]["Z0"][0]["parent_id"]
    missing_owner["owner_by_parent"][first_parent].pop("exact_cell_path", None)
    checks["missing_exact_owner_path_rejected"] = _negative_case(build, c, packet, missing_owner, receipts, "same_geometry", True)

    wrong_receipt = _clone(certificate)
    wrong_receipt["source_receipts"]["LOW3"]["path"] = "experiments/root_st073/v2_wrong_receipt.json.gz"
    checks["wrong_receipt_path_rejected"] = _negative_case(build, c, packet, wrong_receipt, receipts, "same_labels", True)

    insufficient = {"source_family": "fixture-family-v2", "same_g5": True, "same_K0": True, "same_chi": True}
    try:
        insufficient_result = _call(build, c, packet, insufficient, receipts)
        checks["same_family_only_rejected"] = insufficient_result.get("status") == "conditional" and insufficient_result.get("pair_count") == 0
    except Exception:
        checks["same_family_only_rejected"] = False

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks_passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "old_code_called": False,
        "boxed_intervals_summed": False,
        "coefficient_arithmetic_performed": False,
        "symbolic_identity_scope": "conditional symbolic implication under caller assertions",
        "unconditioned_global_zero": False,
    }


run = run_fixture



