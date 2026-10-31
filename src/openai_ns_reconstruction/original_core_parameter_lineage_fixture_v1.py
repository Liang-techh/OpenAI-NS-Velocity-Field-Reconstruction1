"""Synthetic-only fixture for :mod:`original_core_parameter_lineage_v1`.

The fixture supplies independent source texts and small JSON-like mappings to
the caller's verifier.  It never imports or executes the implementation, reads
a report, opens a receipt, or evaluates a producer.  Rejection cases exercise
the direct read path, direct record path, one-bit binary changes, source-family
separation, and missing/mismatched source hashes.  Records with different
display metadata but identical binary endpoints are intentionally accepted.
"""

from __future__ import annotations

import copy
from collections.abc import Callable, Mapping
from typing import Any


Q_SOURCE_PATH = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
BACKGROUND_SOURCE_PATH = (
    "src/openai_ns_reconstruction/original_micro_background_I1_numeric.py"
)


def _record(lower: int, upper: int, display: str) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": [0, lower, -12, 40],
        "upper_exact_mpf_tuple": [0, upper, -12, 40],
        "directed": True,
        "lower_display": display + "-lower",
        "upper_display": display + "-upper",
    }


def _q_source() -> str:
    return '''
def build(root, parent):
    c = object()
    de = read(c, parent["actual_core_signed_direction"]["parameters"]["delta"])
    return dict(passed=True, actual_delta=record(de))
'''


def _background_source() -> str:
    return '''
def _physical_log_scale(c, parent, loghb):
    logP = read(c, parent["actual_core_signed_direction"]["parameters"]["logPstar"])
    return dict(logPstar=record(logP), loghb=record(loghb))

def build(root, parent):
    c = object()
    loghb = object()
    scale = _physical_log_scale(c, parent, loghb)
    return dict(passed=True, physical_scale=scale)
'''


def _inputs() -> tuple[dict[str, Any], dict[str, Any], str, str]:
    family = {
        "actual_five_defect_family_sha256": "fixture-family",
        "implicit_source_sha256": "fixture-source",
        "datum_enclosure_sha256": "fixture-datum",
    }
    delta = _record(3, 4, "core-delta")
    logpstar = _record(5, 6, "core-logPstar")
    q_source = _q_source()
    background_source = _background_source()
    native = {
        "schema": "native_directed_controls_v2_fixture",
        "source_family": copy.deepcopy(family),
        "actual_core_signed_direction": {
            "source_family": copy.deepcopy(family),
            "parameters": {
                "delta": copy.deepcopy(delta),
                "logPstar": copy.deepcopy(logpstar),
            },
        },
        "actual_micro_Q_increment_bounds": {
            "passed": True,
            "source_family": copy.deepcopy(family),
            "actual_delta": _record(3, 4, "native-actual-delta"),
        },
        "actual_micro_background_I1_numeric": {
            "passed": True,
            "source_family": copy.deepcopy(family),
            "physical_scale": {
                "logPstar": _record(5, 6, "native-actual-logPstar"),
            },
        },
    }
    # The hashes are declared by the same saved native packet that carries the
    # actual records.  The verifier computes them from the independently
    # supplied texts, without importing either producer.
    import hashlib
    native["implementation_hashes"] = {
        Q_SOURCE_PATH: hashlib.sha256(q_source.encode("utf-8")).hexdigest(),
        BACKGROUND_SOURCE_PATH: hashlib.sha256(
            background_source.encode("utf-8")
        ).hexdigest(),
    }
    previous = {
        "schema": "original_lower_parameter_binding_v1",
        "source_family": copy.deepcopy(family),
        "delta": {
            "native_record": _record(3, 4, "previous-native-delta"),
            "corrected_record": _record(3, 4, "previous-corrected-delta"),
            "exact_saved_record_match": True,
        },
        "actual_logPstar_parameter": {
            "native_record": _record(5, 6, "previous-native-logPstar"),
            "corrected_record": _record(5, 6, "previous-corrected-logPstar"),
            "exact_saved_record_match": True,
        },
    }
    source_provenance = {
        "native": {
            "path": (
                "experiments/root_st073/"
                "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
            ),
            "sha256": (
                "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
            ),
        },
        "corrected": {
            "path": (
                "experiments/root_st073/"
                "lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json"
            ),
            "sha256": (
                "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5"
            ),
        },
        "producer": {
            "path": "scripts/run_original_corrected_reserved_m15.py",
            "sha256": (
                "9b41d60ab277656734cc98aa3272123a1ce383da4337ba85982dd41142350d6d"
            ),
        },
    }
    # Exercise the published runner shape: the prior leaf and its provenance
    # sibling are nested at result.result.leaf/result.result.source_provenance,
    # rather than passing a live object or a producer module.
    receipt = {
        "result": {
            "result": {
                "leaf": previous,
                "source_provenance": source_provenance,
            }
        }
    }
    return native, receipt, q_source, background_source


def _check(value: Any, expected: bool) -> bool:
    return isinstance(value, Mapping) and value.get("passed") is expected


def _rejects(verifier: Callable[..., Mapping[str, Any]], *args: Any) -> bool:
    try:
        result = verifier(*args)
    except Exception:
        return True
    return isinstance(result, Mapping) and result.get("passed") is False


def run_fixture(verifier: Callable[..., Mapping[str, Any]]) -> dict[str, Any]:
    """Run synthetic acceptance and fail-closed mismatch checks."""

    native, receipt, q_source, background_source = _inputs()
    checks: dict[str, bool] = {}
    try:
        accepted = verifier(native, receipt, q_source, background_source)
        checks["accepts_exact_binary_lineage"] = _check(accepted, True)
        checks["status_proved"] = accepted.get("status") == "proved"
        checks["source_family_checked"] = (
            accepted.get("source_family", {}).get("passed") is True
        )
        checks["producer_ast_paths_checked"] = (
            accepted.get("source_AST", {})
            .get("Q_increment_producer", {})
            .get("direct_parameter_read") is True
            and accepted.get("source_AST", {})
            .get("background_scale_producer", {})
            .get("build_calls_helper") is True
        )
        checks["display_metadata_not_identity"] = (
            bool(accepted.get("comparisons", {})
                 .get("delta", {})
                 .get("native_core_vs_native_actual", {})
                 .get("nonbinary_metadata_differences"))
            and accepted.get("comparisons", {})
                .get("all_exact_binary_tuple_equal") is True
        )
        checks["overlap_not_used"] = (
            accepted.get("comparisons", {}).get("overlap_used_as_identity") is False
        )
        checks["delta_branch_not_reproved"] = (
            accepted.get("provenance_scope", {})
            .get("delta_branch_identity_checked") is False
        )
    except Exception as exc:
        return {
            "passed": False,
            "synthetic_only": True,
            "real_sources_read": False,
            "implementation_imported": False,
            "error": repr(exc),
            "checks": {"accepts_exact_binary_lineage": False},
        }

    wrong_read = q_source.replace(
        '["parameters"]["delta"]',
        '["parameters"]["logPstar"]',
    )
    wrong_read_result = verifier(native, receipt, wrong_read, background_source)
    checks["rejects_wrong_delta_read_path"] = (
        wrong_read_result.get("passed") is False
        and wrong_read_result.get("source_AST", {})
            .get("Q_increment_producer", {})
            .get("direct_parameter_read") is False
    )

    wrong_record = q_source.replace("record(de)", "record(logP)")
    wrong_record_result = verifier(native, receipt, wrong_record, background_source)
    checks["rejects_wrong_actual_record_binding"] = (
        wrong_record_result.get("passed") is False
        and wrong_record_result.get("source_AST", {})
            .get("Q_increment_producer", {})
            .get("direct_record") is False
    )

    one_bit = copy.deepcopy(native)
    one_bit["actual_micro_Q_increment_bounds"]["actual_delta"][
        "lower_exact_mpf_tuple"
    ][1] += 1
    one_bit_result = verifier(one_bit, receipt, q_source, background_source)
    checks["rejects_one_bit_binary_change"] = (
        one_bit_result.get("status") == "mismatch"
        and one_bit_result.get("passed") is False
        and one_bit_result.get("comparisons", {}).get("mismatch_count", 0) > 0
        and one_bit_result.get("comparisons", {})
            .get("mismatches", [{}])[0]
            .get("overlap_used_as_identity") is False
    )

    wrong_family = copy.deepcopy(native)
    wrong_family["actual_micro_background_I1_numeric"]["source_family"] = {
        "actual_five_defect_family_sha256": "other-family",
        "implicit_source_sha256": "other-source",
        "datum_enclosure_sha256": "other-datum",
    }
    wrong_family_result = verifier(wrong_family, receipt, q_source, background_source)
    checks["rejects_cross_source_family"] = (
        wrong_family_result.get("passed") is False
        and wrong_family_result.get("source_family", {}).get("passed") is False
    )

    wrong_hash = copy.deepcopy(native)
    wrong_hash["implementation_hashes"][Q_SOURCE_PATH] = "0" * 64
    wrong_hash_result = verifier(wrong_hash, receipt, q_source, background_source)
    checks["rejects_source_hash_mismatch"] = (
        wrong_hash_result.get("passed") is False
        and wrong_hash_result.get("source_hashes", {})
            .get("Q_increment_producer", {})
            .get("matches") is False
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
