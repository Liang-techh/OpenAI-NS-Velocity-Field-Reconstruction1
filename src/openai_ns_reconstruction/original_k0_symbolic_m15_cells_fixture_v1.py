"""Synthetic structural fixture for the symbolic K0 W1 density leaf.

The fixture creates only decoded JSON-shaped records.  ``run_fixture`` calls
the leaf supplied by its caller; it does not import a producer, read a
receipt, or perform interval arithmetic.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any


SCHEMA = "original_k0_symbolic_m15_cells_fixture_v1"


def _interval(seed: int) -> dict[str, Any]:
    lower, upper = 1 + seed, 2 + seed
    return {
        "lower_exact_mpf_tuple": [0, lower, -8, lower.bit_length()],
        "upper_exact_mpf_tuple": [0, upper, -8, upper.bit_length()],
        "directed": True,
    }


def _saved_term(index: int, order: str, field: str) -> dict[str, Any]:
    f_power = 2 if field == "F1" else 1
    lambda_power = -5 if order == "Z0" else -6
    original = {
        "coefficient": _interval(index + (0 if field == "F1" else 100)),
        "F0_power": f_power,
        "Lambda_power": lambda_power,
        "source_cell": index,
        "target_cell": index,
        "matrix_sector_identity": f"fixture-matrix-{field}-{index}",
        "W0_sector_identity": f"fixture-w0-{field}-{index}",
        "source_id": f"fixture-saved-{field}-{index}",
    }
    return {
        "coefficient": deepcopy(original["coefficient"]),
        "formal_source_powers": [0, f_power, 0, lambda_power],
        "original_term": original,
    }


def _saved_binding() -> dict[str, Any]:
    cells = []
    for index in range(8):
        rows = {
            order: {
                field: [_saved_term(index, order, field)]
                for field in ("F1", "Uz1")
            }
            for order in ("Z0", "Z1")
        }
        cells.append(
            {
                "index": index,
                "source_indices": list(range(index + 1)),
                "triangular_diagonal_source": index,
                "axis_endpoint_continuation": index == 0,
                "t": _interval(index),
                "x": _interval(index + 10),
                "rho": _interval(index + 20),
                "ordinary_Z": rows,
            }
        )
    return {
        "schema": "original_saved_k0_w1_cells_v1",
        "source_family": {
            "actual_five_defect_family_sha256": "a" * 64,
            "implicit_source_sha256": "b" * 64,
            "datum_enclosure_sha256": "c" * 64,
        },
        "source_logs": {"fixture": "synthetic decoded saved result"},
        "formal_source_basis": ["hb", "F0_anchor", "Pstar_squared", "Lambda"],
        "source_kernel_scope": "synthetic K0 only",
        "Kchi_minus_K0_included": False,
        "allword_tail_added": False,
        "pressure_added": False,
        "appended_to_cumulative_m15": False,
        "actual_W1_complete": False,
        "cells": cells,
    }


def _background() -> dict[str, Any]:
    return {
        "U_lead_Z0": {
            **_interval(200),
            "source_id": "fixture-background-U-Z0",
            "term_id": "fixture-background-U-term-Z0",
        },
        "U_lead_Z1": {
            **_interval(201),
            "source_id": "fixture-background-U-Z1",
            "term_id": "fixture-background-U-term-Z1",
        },
        "Phi_Z0": {
            **_interval(202),
            "source_id": "fixture-background-Phi-Z0",
            "term_id": "fixture-background-Phi-term-Z0",
        },
        "Phi_Z1": {
            **_interval(203),
            "source_id": "fixture-background-Phi-Z1",
            "term_id": "fixture-background-Phi-term-Z1",
        },
    }


def _call(leaf: Any, saved_binding: Any, saved_background: Any) -> Any:
    function = getattr(leaf, "build", leaf)
    if not callable(function):
        raise TypeError("leaf must expose build")
    return function(saved_binding, saved_background)


def _reject(leaf: Any, saved_binding: Any, saved_background: Any) -> bool:
    try:
        _call(leaf, saved_binding, saved_background)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run structural checks against a caller-supplied leaf."""

    saved_binding = _saved_binding()
    saved_background = _background()
    before_binding = deepcopy(saved_binding)
    before_background = deepcopy(saved_background)
    result = _call(leaf, saved_binding, saved_background)
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "original_k0_symbolic_m15_cells_v1",
        "eight_target_cells": len(result.get("cells", ())) == 8,
        "formula_channels": set(result.get("channels", ())) == {"axial", "swirl"},
        "order_channel_counts": result.get("term_counts") == {
            "axial": {"Z0": 8, "Z1": 16},
            "swirl": {"Z0": 8, "Z1": 24},
        },
        "saved_input_unmodified": saved_binding == before_binding and saved_background == before_background,
        "binary_coefficient_unchanged": result["cells"][3]["terms"]["axial"]["Z0"][0]["coefficient"]
        == before_binding["cells"][3]["ordinary_Z"]["Z0"]["Uz1"][0]["coefficient"],
        "background_interval_embedded": result["cells"][3]["terms"]["swirl"]["Z0"][0]["background_factor"]
        == before_background["Phi_Z0"],
        "rho_interval_embedded": result["cells"][3]["terms"]["swirl"]["Z0"][0]["rho_factor"]
        == before_binding["cells"][3]["rho"],
        "signed_factors": result["cells"][3]["terms"]["axial"]["Z0"][0]["signed_factor"] == 2
        and result["cells"][3]["terms"]["swirl"]["Z0"][0]["signed_factor"] == -2,
        "axial_lambda_external": result["cells"][3]["terms"]["axial"]["Z0"][0]["formal_source_powers"]
        == [0, 1, 0, -6, 0, 0],
        "swirl_scales_independent": result["cells"][3]["terms"]["swirl"]["Z0"][0]["formal_source_powers"]
        == [0, 2, 0, -7, 1, 0]
        and result["cells"][3]["terms"]["swirl"]["Z1"][1]["formal_source_powers"]
        == [0, 2, 0, -7, 1, 1],
        "ordinary_Z1_saved_rows_reused": result["ordinary_Z_orders"] == ["Z0", "Z1"]
        and result["ordinary_Z1_saved_rows_reused"] is True
        and result["saved_Z1_coefficients_or_powers_differentiated"] is False,
        "source_target_matrix_w0_retained": all(
            term.get("source_cell") == 3
            and term.get("target_cell") == 3
            and term.get("matrix_sector_identity")
            and term.get("W0_sector_identity")
            for term in result["cells"][3]["terms"]["axial"]["Z1"]
        ),
        "source_identity_unverified": result.get("source_identity_status") == "source identity unverified"
        and result.get("source_identity_verified") is False,
        "F0_scales_not_identified": result.get("F0_scale_identity", {}).get("identified") is False,
        "unsupported_terms_unknown": all(item.get("status") == "unknown" for item in result.get("unknown_terms", ()))
        and result.get("unknown_terms_are_not_zero") is True,
        "no_pressure_kchi_tail_ledger": result.get("N1p_included") is False
        and result.get("Kchi_included") is False
        and result.get("allword_tail_included") is False
        and result.get("appended_to_212_ledger") is False,
        "no_interval_operators": result.get("interval_arithmetic_performed") is False
        and result.get("integration_performed") is False,
    }

    bad = deepcopy(saved_background)
    del bad["U_lead_Z0"]["source_id"]
    del bad["U_lead_Z0"]["term_id"]
    checks["background_without_source_rejected"] = _reject(leaf, saved_binding, bad)

    bad = deepcopy(saved_binding)
    bad["cells"][1], bad["cells"][2] = bad["cells"][2], bad["cells"][1]
    checks["cell_order_rejected"] = _reject(leaf, bad, saved_background)

    bad = deepcopy(saved_binding)
    del bad["cells"][0]["ordinary_Z"]["Z1"]["Uz1"]
    checks["missing_required_saved_sector_rejected"] = _reject(leaf, bad, saved_background)

    bad = deepcopy(saved_binding)
    bad["cells"][0]["ordinary_Z"]["Z0"]["F1"][0]["coefficient"] = _interval(8)
    checks["changed_coefficient_rejected"] = _reject(leaf, bad, saved_background)

    bad = deepcopy(saved_binding)
    bad["cells"][0]["ordinary_Z"]["Z0"]["F1"][0]["formal_source_powers"][1] += 1
    checks["changed_source_power_rejected"] = _reject(leaf, bad, saved_background)

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_values_computed": False,
        "implementation_imported": False,
        "old_code_called": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
