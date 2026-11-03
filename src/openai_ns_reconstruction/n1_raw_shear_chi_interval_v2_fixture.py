"""Synthetic v2 fixture for the runner's source-family projection boundary.

The runner may decode a complete source-family receipt containing unrelated
metadata, including floating-point fields.  It must validate that receipt
first and pass only the validated ``{"source_family": family}`` projection to
the v1 leaf.  This fixture models that boundary and calls the supplied leaf
directly; it imports no producer or prior fixture.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SCHEMA = "n1_raw_shear_chi_interval_v2_fixture"
SYNTHETIC_GRAPH_HASH = "synthetic-raw-shear-chi-interval-v2"


def _rows(family: Any) -> dict[int, dict[str, Any]]:
    """Use the saved raw-program row shapes, with a synthetic family only."""

    return {
        1: {"node_id": 1, "operation": "exact_rational", "value": 1},
        25085: {"node_id": 25085, "operation": "negative", "argument": 1},
        8: {
            "node_id": 8,
            "operation": "current_original_source_parameter",
            "name": "hbB",
            "source_family": family,
        },
        35926: {"node_id": 35926, "operation": "sum", "arguments": [8, 25085]},
        35950: {"node_id": 35950, "operation": "product", "arguments": [35926, 35949]},
        35951: {"node_id": 35951, "operation": "sum", "arguments": [1, 35950]},
        35954: {"node_id": 35954, "operation": "product", "arguments": [35926, 35953]},
    }


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run one successful leaf call and two narrow negative controls."""

    F = Fraction
    names = (
        "api",
        "receipt_projection_excludes_unrelated_float",
        "leaf_call_numeric_intervals",
        "signed_s_and_Z_distinction",
        "source_owner_flags_unproved",
        "family_mismatch_rejected",
        "bad_row_rejected",
    )
    checks = {name: False for name in names}
    diagnostics: dict[str, Any] = {}
    required = ("compute", "EXPECTED_GRAPH_SHA256", "CHI_NODE_ID", "SIGMA_S_PRODUCT_NODE_ID")
    checks["api"] = all(hasattr(leaf, name) for name in required)
    if not checks["api"]:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": checks,
            "diagnostics": {"missing_api": [name for name in required if not hasattr(leaf, name)]},
            "synthetic_only": True,
        }

    F = Fraction
    family = {"synthetic_family": "raw_shear_chi_v2", "revision": 23}
    rows = _rows(family)
    hb_interval = [F(1, 4), F(1, 3)]
    sigma_interval = [F(1, 8), F(1, 4)]
    sigma_s_interval = [F(1, 16), F(1, 8)]
    hb_evidence = {
        "source_family": family,
        "strict_hb_admitted": True,
        "hb_sc_Z_independent": True,
        "hb_interval": hb_interval,
        "fixed_s": True,
        "fixed_hb": True,
        "fixed_s_hb": True,
        "s_derivative_not_Z_derivative": True,
        "ordinary_Z1": F(0),
    }
    z_condition = {
        "fixed_s": True,
        "fixed_hb": True,
        "fixed_s_hb": True,
        "s_derivative_not_Z_derivative": True,
        "ordinary_Z1": F(0),
    }

    # This models a complete validated receipt.  The unrelated float is not
    # handed to the leaf: the runner's exact projection is the only receipt
    # value consumed by the source-family check.
    complete_receipt = {
        "source_family": family,
        "unrelated_float_metadata": 0.125,
        "nested": {"status": "validated", "other_value": 7.25},
    }
    projected_receipt = {"source_family": complete_receipt["source_family"]}
    common = {
        "graph_hash": SYNTHETIC_GRAPH_HASH,
        "source_family": family,
        "source_family_receipt": projected_receipt,
        "z_condition": z_condition,
        "rows_authenticated": True,
        "sigma_bounds_authenticated": True,
        "synthetic": True,
    }

    try:
        # The one successful source test case calls the real leaf once.
        result = leaf.compute(rows, sigma_interval, sigma_s_interval, hb_evidence, **common)
        checks["receipt_projection_excludes_unrelated_float"] = (
            result.get("source_family_receipt_checked") is True
            and projected_receipt == {"source_family": family}
            and complete_receipt["unrelated_float_metadata"] == 0.125
        )
        checks["leaf_call_numeric_intervals"] = (
            result.get("chi_interval") == {"lower": F(13, 16), "upper": F(11, 12)}
            and result.get("chi_s_interval") == {"lower": F(-3, 32), "upper": F(-1, 24)}
            and result.get("arithmetic_call_count") == 4
        )
        checks["signed_s_and_Z_distinction"] = (
            result.get("ordinary_Z1", {}).get("value") == F(0)
            and result.get("coordinate") == "raw_s"
            and result.get("s_derivative_is_Z_derivative") is False
            and result.get("s_derivative_is_q_derivative") is False
            and result.get("chi_s_is_Z_derivative") is False
            and result.get("chi_s_is_q_derivative") is False
        )
        checks["source_owner_flags_unproved"] = (
            result.get("source_identity_verified") is False
            and result.get("owner_binding") is False
            and result.get("live_owner_verified") is False
            and result.get("orientation_identity_proved") is False
            and result.get("physical_m15") is False
        )
        wrong_family = {"source_family": {"synthetic_family": "other", "revision": 23}}
        checks["family_mismatch_rejected"] = _rejected(
            lambda: leaf.compute(
                rows,
                sigma_interval,
                sigma_s_interval,
                hb_evidence,
                **dict(common, source_family_receipt=wrong_family),
            )
        )
        forged_rows = dict(rows)
        forged_rows[35954] = {
            "node_id": 35954,
            "operation": "sum",
            "arguments": [35926, 35953],
        }
        checks["bad_row_rejected"] = _rejected(
            lambda: leaf.compute(forged_rows, sigma_interval, sigma_s_interval, hb_evidence, **common)
        )
        diagnostics = {
            "chi_interval": result.get("chi_interval"),
            "chi_s_interval": result.get("chi_s_interval"),
            "synthetic_graph_hash": SYNTHETIC_GRAPH_HASH,
            "source_family_projection_only": True,
        }
    except Exception as exc:
        diagnostics = {"exception": type(exc).__name__, "message": str(exc)}

    return {
        "schema": SCHEMA,
        "component": "raw_g5_lower_collar_shear_chi_interval_v2_fixture",
        "passed": all(checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "synthetic_only": True,
        "source_family_projection_only": True,
        "used_old_fixture": False,
        "used_old_module": False,
        "gzip_reads": False,
    }


__all__ = ["SCHEMA", "SYNTHETIC_GRAPH_HASH", "run_fixture"]
