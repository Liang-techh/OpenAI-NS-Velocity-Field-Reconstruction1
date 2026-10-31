"""Fresh synthetic fixture for the lower-parameter binding leaf.

The fixture never imports or executes the implementation.  A caller supplies
the leaf callable and this module supplies small mappings plus rejection cases
for a hash mismatch, a one-bit endpoint perturbation, and a wrong producer
path.
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Mapping
from typing import Any


EXPECTED_NATIVE_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)


def _record(lower: int, upper: int) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": [0, lower, -12, 40],
        "upper_exact_mpf_tuple": [0, upper, -12, 40],
        "directed": True,
        "synthetic_only": True,
    }


def _inputs() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str]:
    family = {
        "actual_five_defect_family_sha256": "fixture-family",
        "implicit_source_sha256": "fixture-source",
        "datum_enclosure_sha256": "fixture-datum",
    }
    delta = _record(3, 4)
    logpstar = _record(5, 6)
    native = {
        "source_family": copy.deepcopy(family),
        "actual_micro_Q_increment_bounds": {
            "source_family": copy.deepcopy(family),
            "actual_delta": copy.deepcopy(delta),
        },
        "actual_micro_background_I1_numeric": {
            "source_family": copy.deepcopy(family),
            "physical_scale": {"logPstar": copy.deepcopy(logpstar)},
        },
    }
    corrected = {
        "schema": "corrected_reserved_m15_inputs_v1",
        "source_only": True,
        "old_212_terms_and_all_errors_unchanged": True,
        "corrected_reserved_m15_source": {
            "source_only": True,
            "source_executed": False,
            "existing_ledger_edited": False,
        },
        "corrected_terminal_source_inputs": {
            "schema": "corrected_reserved_m15_inputs_v1",
            "source_only": True,
            "native_sha256": EXPECTED_NATIVE_SHA256,
            "parent_sha256": "fixture-parent",
            "old_graphs_unchanged": True,
            "numerical_ancestor_executed": False,
            "numeric_endpoint_cell": ["0", ".5"],
            "physical_endpoint_source_rows": {"R": {}, "Pstar": {}},
            "physical_endpoint_source_proof": {
                "source_family": copy.deepcopy(family),
                "actual_V_y_rows_0_through_3_factorial_converted_once": True,
                "common_Pstar_applied_once": True,
            },
            "delta": copy.deepcopy(delta),
            "actual_logPstar_parameter": copy.deepcopy(logpstar),
        },
    }
    receipt = {
        "result": {
            "result": {
                "source_family": copy.deepcopy(family),
                "corrected_terminal_source_contract": {
                    "actual_delta": copy.deepcopy(delta),
                    "actual_logPstar_parameter": copy.deepcopy(logpstar),
                    "source_only": True,
                },
                "ledger_written": False,
                "old_errors_removed": False,
                "source_rows_are_not_point_values": True,
                "endpoint_function_values_selected": False,
                "leaf": {
                    "ordinary_Z1_is_product_rule": True,
                    "ordinary_Z1_terms_not_merged": True,
                },
            }
        }
    }
    source_text = '''
import copy

PREFIX = "fixture_"
PARENT = PREFIX + "outer_pressure_enclosure.json"
PARENT_HASH = "fixture-parent"
NATIVE = PREFIX + "native_directed_controls_v2.json"
NATIVE_HASH = "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"

def setup():
    native = read_pinned(BASE / NATIVE, NATIVE_HASH)
    source = dict(
        parent_sha256=PARENT_HASH,
        native_sha256=NATIVE_HASH,
        delta=copy.deepcopy(native["actual_micro_Q_increment_bounds"]["actual_delta"]),
        actual_logPstar_parameter=copy.deepcopy(
            native["actual_micro_background_I1_numeric"]["physical_scale"]["logPstar"]
        ),
    )
    return source
'''
    return native, corrected, receipt, source_text


def _rejects(build: Callable[..., Mapping[str, Any]], *args: Any) -> bool:
    try:
        build(*args)
    except Exception:
        return True
    return False


def run_fixture(build: Callable[..., Mapping[str, Any]]) -> dict[str, Any]:
    """Run synthetic-only acceptance and rejection checks against a leaf."""

    native, corrected, receipt, source_text = _inputs()
    checks: dict[str, bool] = {}
    try:
        result = build(native, corrected, receipt, source_text)
        checks["accepted_exact_records"] = result.get("passed") is True
        ast_result = result.get("source_ast_contract", {})
        checks["direct_ast_producers"] = (
            ast_result.get("delta_direct_deepcopy") is True
            and ast_result.get("actual_logPstar_parameter_direct_deepcopy") is True
        )
        checks["live_identity_not_claimed"] = (
            result.get("provenance_scope", {}).get("endpoint_live_owner_identity") is False
        )
        checks["z1_not_rederived"] = (
            result.get("ordinary_Z1_scope", {}).get("product_rule_rederived_here") is False
        )
    except Exception as exc:
        return {
            "passed": False,
            "synthetic_only": True,
            "real_sources_read": False,
            "implementation_imported": False,
            "error": repr(exc),
            "checks": {"accepted_exact_records": False},
        }

    bad_hash = copy.deepcopy(corrected)
    bad_hash["corrected_terminal_source_inputs"]["native_sha256"] = "0" * 64
    checks["rejects_native_hash_mismatch"] = _rejects(
        build, native, bad_hash, receipt, source_text
    )

    one_bit = copy.deepcopy(corrected)
    one_bit["corrected_terminal_source_inputs"]["delta"][
        "lower_exact_mpf_tuple"
    ][1] += 1
    checks["rejects_one_bit_endpoint_change"] = _rejects(
        build, native, one_bit, receipt, source_text
    )

    wrong_producer = source_text.replace(
        'native["actual_micro_Q_increment_bounds"]["actual_delta"]',
        'native["wrong_parameter"]["actual_delta"]',
    )
    checks["rejects_incorrect_direct_producer"] = _rejects(
        build, native, corrected, receipt, wrong_producer
    )
    wrong_native_pin = source_text.replace(
        "BASE / NATIVE, NATIVE_HASH",
        "BASE / WRONG_NATIVE, NATIVE_HASH",
    )
    checks["rejects_wrong_native_pin_path"] = _rejects(
        build, native, corrected, receipt, wrong_native_pin
    )
    return {
        "passed": all(checks.values()),
        "synthetic_only": True,
        "real_sources_read": False,
        "implementation_imported": False,
        "checks": checks,
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
