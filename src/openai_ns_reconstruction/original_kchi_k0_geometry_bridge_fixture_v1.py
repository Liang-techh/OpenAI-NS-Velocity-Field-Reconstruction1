"""Synthetic checks for the Kchi/K0 coordinate bridge.

The fixture creates only small JSON-shaped coordinate packets.  It does not
read a receipt, import a producer, inspect coefficient terms, or represent an
original source.  Negative cases deliberately break one coordinate, sign,
index, source/target ownership, source log, or bounded-exponent requirement.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from typing import Any


SCHEMA = "original_kchi_k0_geometry_bridge_fixture_v1"


def _binary(value: Fraction) -> list[int]:
    """Encode a small dyadic Fraction in the saved four-integer shape."""

    if value == 0:
        return [0, 0, 0, 0]
    denominator = value.denominator
    if denominator <= 0 or denominator & (denominator - 1):
        raise ValueError("fixture coordinates must be dyadic")
    mantissa = abs(value.numerator)
    exponent = -(denominator.bit_length() - 1)
    return [
        1 if value < 0 else 0,
        mantissa,
        exponent,
        mantissa.bit_length(),
    ]


def _interval(lower: Fraction, upper: Fraction) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": _binary(lower),
        "upper_exact_mpf_tuple": _binary(upper),
        "directed": True,
    }


def _t(index: int) -> tuple[Fraction, Fraction]:
    return Fraction(index, 8), Fraction(index + 1, 8)


def _x(index: int) -> tuple[Fraction, Fraction]:
    lower, upper = _t(index)
    return 2 * lower, 2 * upper


def _rho(index: int) -> tuple[Fraction, Fraction]:
    lower, upper = _t(index)
    return 4 * lower * lower, 4 * upper * upper


def _packets() -> tuple[dict[str, Any], dict[str, Any]]:
    family = {"fixture_family": "same-source-label"}
    logs = {
        "hb": "fixture-hb",
        "F0_anchor": "fixture-F0-anchor",
        "Pstar_squared": "fixture-Pstar-squared",
        "Lambda": "fixture-Lambda",
    }
    kchi_cells: list[dict[str, Any]] = []
    k0_cells: list[dict[str, Any]] = []
    for index in range(8):
        x_lower, x_upper = _x(index)
        x_record = _interval(x_lower, x_upper)
        source_cells: list[dict[str, Any]] = []
        for source_index in range(index + 1):
            source_lower, source_upper = _x(source_index)
            source_cells.append(
                {
                    "source_index": source_index,
                    "target_index": index,
                    "triangular": source_index == index,
                    "source_x": _interval(source_lower, source_upper),
                    "target_x": deepcopy(x_record),
                }
            )
        kchi_cells.append(
            {
                "index": index,
                "target_cell": index,
                "target_index": index,
                "source_indices": list(range(index + 1)),
                "triangular_diagonal_source": index,
                "source_j_le_target_i": True,
                "x": x_record,
                "target_x": deepcopy(x_record),
                "source_cells": source_cells,
            }
        )

        t_lower, t_upper = _t(index)
        rho_lower, rho_upper = _rho(index)
        k0_cells.append(
            {
                "index": index,
                "target_cell": index,
                "source_indices": list(range(index + 1)),
                "triangular_diagonal_source": index,
                "t": _interval(t_lower, t_upper),
                "x": _interval(x_lower, x_upper),
                "rho": _interval(rho_lower, rho_upper),
            }
        )

    return (
        {
            "schema": "original_kchi_first_saved_cells_v2",
            "source_family": deepcopy(family),
            "source_logs": deepcopy(logs),
            "cells": kchi_cells,
        },
        {
            "schema": "original_k0_symbolic_m15_cells_v1",
            "source_family": deepcopy(family),
            "source_logs": deepcopy(logs),
            "cells": k0_cells,
        },
    )


def _call(leaf: Any, kchi: Any, k0: Any) -> Any:
    function = getattr(leaf, "build", leaf)
    if not callable(function):
        raise TypeError("leaf must expose build")
    return function(kchi, k0)


def _reject(leaf: Any, kchi: Any, k0: Any) -> bool:
    try:
        _call(leaf, kchi, k0)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run positive and fail-closed structural checks against a supplied leaf."""

    kchi, k0 = _packets()
    result = _call(leaf, kchi, k0)
    cell3 = result["cells"][3]
    checks: dict[str, bool] = {
        "schema": result.get("schema") == "original_kchi_k0_geometry_bridge_v1",
        "eight_cells": result.get("cell_count") == 8 and len(result.get("cells", ())) == 8,
        "x_t_rho_records_preserved": all(
            key in cell3["coordinates"]
            for key in ("kchi_x", "kchi_target_x", "k0_t", "k0_x", "k0_rho")
        ),
        "canonical_transform_flags": result.get("coordinate_correspondence_verified") is True
        and result.get("rho_equals_x_squared_verified") is True,
        "source_target_structure": cell3["kchi"]["source_indices"] == [0, 1, 2, 3]
        and cell3["kchi"]["source_cells"][3]["triangular"] is True
        and cell3["k0"]["source_indices"] == [0, 1, 2, 3],
        "cell3_exact_weights": cell3["weights"]["axial"]
        == {"numerator": 7, "denominator": 8}
        and cell3["weights"]["swirl"] == {"numerator": -175, "denominator": 256},
        "all_weights_predefined": all(
            item["weights"]["is_predefined_rational_weight"]
            and item["weights"]["is_integral"] is False
            for item in result["cells"]
        ),
        "same_source_labels_conditional": result.get("source_family_logs_equal") is True
        and result.get("source_identity_verified") is False
        and result.get("live_owner_verified") is False,
        "no_coefficients_or_integrals": result.get("coefficient_records_read") is False
        and result.get("integration_performed") is False
        and result.get("integrals_claimed") is False,
        "no_ledger_or_W1_claim": result.get("ledger_changed") is False
        and result.get("complete_W1") is False
        and result.get("complete_m15") is False,
    }

    bad_kchi, bad_k0 = _packets()
    bad_kchi["cells"][3]["x"] = _interval(Fraction(1), Fraction(5, 4))
    checks["wrong_kchi_x_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_kchi["cells"][4]["target_x"] = _interval(Fraction(1), Fraction(5, 4))
    checks["wrong_kchi_target_x_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_k0["cells"][4]["rho"] = _interval(Fraction(0), Fraction(1, 4))
    checks["wrong_k0_rho_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_k0["cells"][1]["x"] = _interval(Fraction(-1, 4), Fraction(0))
    checks["negative_coordinate_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_kchi["cells"][2]["index"] = 3
    checks["wrong_cell_index_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_kchi["cells"][4]["source_cells"][2]["triangular"] = True
    checks["wrong_triangular_source_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_k0["cells"][5]["source_indices"] = [0, 1, 2, 3, 5, 4]
    checks["wrong_source_indices_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_k0["source_logs"]["Lambda"] = "different-source-log"
    checks["source_log_mismatch_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_kchi["source_family"]["fixture_family"] = "different-source-family"
    checks["source_family_mismatch_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    bad_k0["cells"][0]["t"]["upper_exact_mpf_tuple"][2] = 257
    checks["huge_coordinate_exponent_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    bad_kchi, bad_k0 = _packets()
    del bad_kchi["cells"][0]["x"]["directed"]
    checks["missing_directed_marker_rejected"] = _reject(leaf, bad_kchi, bad_k0)

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "coefficient_records_computed": False,
        "integration_performed": False,
        "saved_data_read": False,
    }


run = run_fixture


__all__ = ["SCHEMA", "run_fixture", "run"]
