"""Independent synthetic checks for ``original_kchi_p_tail_m15_v1``.

The fixture creates eight target cells, five tail branches, signed/symmetric
error records, and two symbolic Phi wrappers.  It does not import a producer
or any prior mathematical helper.  The hand reference is rational and only
checks enclosure containment; it never treats an error bound as a field value.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from typing import Any


BRANCHES = ("G0_P_Rg", "D_P_Rg", "R_P_G0g", "R_P_Dg", "R_P_Rg")
ORDERS = ("Z0", "Z1")


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _record(c: Any, low: Fraction, high: Fraction, *, directed: bool = True) -> dict[str, Any]:
    left = _q(c, low)
    right = _q(c, high)
    result = {
        "lower_exact_mpf_tuple": list(left._mpi_[0]),
        "upper_exact_mpf_tuple": list(right._mpi_[1]),
    }
    if directed:
        result["directed"] = True
    return result


def _term(c: Any, target: int, source: int, order: str,
          branch: str, position: int, radius: Fraction,
          powers: tuple[int, int, int, int]) -> dict[str, Any]:
    source_id = ("tail:" + branch + ":" + str(target) + ":" + str(source)
                 + ":" + order + ":" + str(position))
    return {
        "coefficient": _record(c, -radius, radius),
        "powers": list(powers),
        "source_id": source_id,
        "source_ids": [source_id],
        "error_classes": ["synthetic_target_tail", branch],
        "error_only": True,
        "symmetric_error": True,
        "signed_field": False,
        "uncertainty": True,
        "target_cell": target,
        "source_cell": source,
    }


def _tail(c: Any) -> dict[str, Any]:
    branches: dict[str, Any] = {}
    for branch_index, branch in enumerate(BRANCHES):
        targets: list[dict[str, Any]] = []
        for target in range(8):
            by_source: dict[str, Any] = {}
            for source in range(target + 1):
                # Distinct E0/E1 radii and alternating source signs make the
                # three Z1 routes observable without ever summing them.
                e0_radius = Fraction(branch_index + target + source + 1, 1024)
                e1_radius = Fraction(branch_index + 2 * target + source + 2, 2048)
                e0 = _term(c, target, source, "Z0", branch, 0, e0_radius,
                           (1, 2, 3, -5 - branch_index))
                e1 = _term(c, target, source, "Z1", branch, 0, e1_radius,
                           (1, 2, 3, -5 - branch_index))
                by_source[str(source)] = {
                    "target_index": target,
                    "source_index": source,
                    "fields": {
                        "Z0": {"F1": [e0]},
                        "Z1": {"F1": [e1]},
                        # The consumer must not evaluate this field.
                        "d_x_F1": {"Z0": [{"not_to_be_read": True}]},
                    },
                }
            targets.append({
                "target_index": target,
                "target_radius_segment": [_record(c, Fraction(target, 4), Fraction(target + 1, 4)),
                                           _record(c, Fraction(target, 4), Fraction(target + 1, 4))],
                "by_source_cell": by_source,
            })
        branches[branch] = {
            "branch": branch,
            "error_only": True,
            "symmetric_error": True,
            "target_cells": targets,
        }
    return {
        "schema": "synthetic_kchi_target_tail_v1",
        "source_family": {"fixture_family": "synthetic-kchi-family"},
        "source_logs": {
            "hb": "hb-source", "F0_anchor": "f0-source",
            "Pstar_squared": "pstar-source", "Lambda": "lambda-source",
        },
        "branches": branches,
    }


def _geometry(c: Any, family: Any, logs: Any) -> dict[str, Any]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        x = _record(c, Fraction(index, 4), Fraction(index + 1, 4))
        t = _record(c, Fraction(index, 8), Fraction(index + 1, 8))
        rho = _record(c, Fraction(index * index, 16), Fraction((index + 1) * (index + 1), 16))
        weight = -16 * (Fraction(index + 1, 8)**4 - Fraction(index, 8)**4)
        cells.append({
            "index": index,
            "coordinates": {"kchi_x": x, "k0_t": t, "k0_rho": rho},
            "weights": {
                "swirl": {"numerator": weight.numerator, "denominator": weight.denominator},
            },
        })
    return {
        "schema": "synthetic_kchi_geometry_bridge_v2",
        "status": "coordinate_weight_verified",
        "coordinate_verified": True,
        "coordinate_correspondence_verified": True,
        "cutoff_owner_verified": False,
        "source_family": deepcopy(family),
        "source_logs": deepcopy(logs),
        "cells": cells,
    }


def _symbolic(c: Any, family: Any, logs: Any) -> dict[str, Any]:
    phi0 = {
        "interval": _record(c, Fraction(-1, 3), Fraction(1, 2)),
        "source_name": "synthetic-Phi0",
        "source_ids": ["phi-source-Z0"],
        "source_path": "synthetic/Phi_Z0",
    }
    phi1 = {
        "interval": _record(c, Fraction(-2, 5), Fraction(3, 7)),
        "source_name": "synthetic-Phi1",
        "source_ids": ["phi-source-Z1"],
        "source_path": "synthetic/Phi_Z1",
    }
    cells: list[dict[str, Any]] = []
    for index in range(8):
        cells.append({
            "index": index,
            "terms": {
                "swirl": {
                    "Z0": [{"background_factor_name": "Phi_Z0",
                             "background_factor": deepcopy(phi0)}],
                    "Z1": [{"background_factor_name": "Phi_Z1",
                             "background_factor": deepcopy(phi1)}],
                },
            },
        })
    return {
        "schema": "original_k0_symbolic_m15_cells_v1",
        "source_family": deepcopy(family),
        "source_logs": deepcopy(logs),
        "source_power_basis": ["hb", "F0_anchor", "Pstar_squared", "Lambda",
                                "F0_local", "L_local"],
        "cells": cells,
    }


def _connection(c: Any, family: Any, logs: Any) -> dict[str, Any]:
    # Deliberately omit ``directed`` on this historical log record.  The
    # consumer accepts the omission but rejects an explicit false marker.
    ell = _record(c, Fraction(-1, 4), Fraction(2, 3), directed=False)
    return {
        "source_family": deepcopy(family),
        "source_logs": deepcopy(logs),
        "result": {"unrelated_saved_swirl_payload": True},
        "F0_connection": {
            "ell_over_lambda": ell,
            "source_id": "synthetic-F0-connection",
            "source_path": "synthetic/F0_connection/ell_over_lambda",
        },
    }


def _fraction_endpoint(bits: Any) -> Fraction:
    sign, mantissa, exponent, _count = bits
    value = Fraction(int(mantissa), 1)
    if exponent >= 0:
        value *= 1 << int(exponent)
    else:
        value /= 1 << int(-exponent)
    return -value if sign else value


def _lower(value: Any) -> Fraction:
    return _fraction_endpoint(value._mpi_[0])


def _upper(value: Any) -> Fraction:
    return _fraction_endpoint(value._mpi_[1])


def _terms(result: dict[str, Any], branch: str, target: int = 3,
           source: int = 1) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    row = result["branches"][branch]["target_cells"][target]["by_source_cell"][str(source)]
    return row["fields"]["Z0"]["F1"], row["fields"]["Z1"]["F1"]


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent synthetic route, sign, and rejection checks."""

    family = {"fixture_family": "synthetic-kchi-family"}
    logs = {
        "hb": "hb-source", "F0_anchor": "f0-source",
        "Pstar_squared": "pstar-source", "Lambda": "lambda-source",
    }
    tail = _tail(c)
    geometry = _geometry(c, family, logs)
    symbolic = _symbolic(c, family, logs)
    connection = _connection(c, family, logs)
    result = leaf(c, tail, geometry, symbolic, connection)
    checks: dict[str, bool] = {}
    checks["five_branches"] = tuple(result.get("branch_names", ())) == BRANCHES
    checks["eight_cells_each"] = checks["five_branches"] and all(
        len(result["branches"][branch]["target_cells"]) == 8
        for branch in BRANCHES
    )
    checks["only_F1_output"] = checks["five_branches"] and all(
        set(source["fields"]) == {"Z0", "Z1"}
        and set(source["fields"]["Z0"]) == {"F1"}
        and set(source["fields"]["Z1"]) == {"F1"}
        for branch in BRANCHES
        for target in result["branches"][branch]["target_cells"]
        for source in target["by_source_cell"].values()
    )

    z0, z1 = _terms(result, BRANCHES[0])
    checks["Z0_one_route"] = len(z0) == 1 and z0[0]["route"] == "Phi0_E0"
    checks["Z1_three_routes"] = len(z1) == 3 and {
        term["route"] for term in z1
    } == {"Phi1_E0", "Phi0_E1", "Phi0_ell_E0"}
    checks["all_symmetric_error"] = all(
        term["error_only"] is True and term["symmetric_error"] is True
        and term["signed_field"] is False
        and term["uncertainty"] is True
        and _lower(term["coefficient"]) == -_upper(term["coefficient"])
        for branch in BRANCHES
        for target in result["branches"][branch]["target_cells"]
        for source in target["by_source_cell"].values()
        for order in ORDERS
        for term in source["fields"][order]["F1"]
    )
    checks["six_axis_powers"] = all(
        len(term["formal_source_powers"]) == 6
        and term["formal_source_powers"][4] == 1
        and term["formal_source_powers"][5] == 0
        and term["formal_source_powers"][3]
        == term["saved_formal_source_powers"][3] - (1 if term["route"] == "Phi0_ell_E0" else 2)
        for branch in BRANCHES
        for target in result["branches"][branch]["target_cells"]
        for source in target["by_source_cell"].values()
        for order in ORDERS
        for term in source["fields"][order]["F1"]
    )
    checks["ell_branch_is_conditional"] = (
        next(term for term in z1 if term["route"] == "Phi0_ell_E0")["conditional_identity"]
        == "L_local = Lambda * ell_over_lambda"
        and result["connection_provenance"]["connection_recomputed"] is False
        and result["F0_anchor_local_independent"] is True
    )
    checks["weight_applied_once"] = all(
        term["weight_applied_once"] is True
        and term["no_extra_minus_two"] is True
        and term["geometry_weight"]["absolute"] is True
        for branch in BRANCHES
        for target in result["branches"][branch]["target_cells"]
        for source in target["by_source_cell"].values()
        for order in ORDERS
        for term in source["fields"][order]["F1"]
    )
    checks["original_errors_retained"] = all(
        term["original_error"] == term["input_error"]
        and term["input_term_position"] == 0
        and term["input_source_ids"]
        for branch in BRANCHES
        for target in result["branches"][branch]["target_cells"]
        for source in target["by_source_cell"].values()
        for order in ORDERS
        for term in source["fields"][order]["F1"]
    )
    checks["input_path_uses_input_order"] = all(
        (".fields." + term["input_order"] + ".F1[") in term["input_path"]
        for term in z0 + z1
    )
    # Cell 3/source 1, first branch: |swirl|=175/256, |Phi0|=1/2, E0=5/1024.
    # The returned enclosure must contain this independently computed rational.
    hand_radius = Fraction(175, 256) * Fraction(1, 2) * Fraction(5, 1024)
    checks["independent_hand_rational_enclosure"] = (
        _lower(z0[0]["coefficient"]) <= -hand_radius
        and _upper(z0[0]["coefficient"]) >= hand_radius
    )
    checks["conditional_prior_flags"] = (
        result["source_owner_identity_proved"] is False
        and result["geometry"]["cutoff_same_source_owner_proved"] is False
        and result["pressure_included"] is False
        and result["collar_included"] is False
        and result["full_word_included"] is False
        and result["appended_to_ledger"] is False
    )

    bad_tail = deepcopy(tail)
    bad_tail["branches"][BRANCHES[0]]["target_cells"][0]["by_source_cell"]["0"]["fields"]["Z0"]["F1"][0]["error_only"] = False
    try:
        leaf(c, bad_tail, geometry, symbolic, connection)
    except (TypeError, ValueError, KeyError):
        checks["reject_not_error"] = True
    else:
        checks["reject_not_error"] = False

    bad_tail = deepcopy(tail)
    bad_tail["branches"][BRANCHES[0]]["target_cells"][0]["by_source_cell"]["0"]["fields"]["Z0"]["F1"][0]["coefficient"] = _record(c, Fraction(-1), Fraction(2))
    try:
        leaf(c, bad_tail, geometry, symbolic, connection)
    except (TypeError, ValueError, KeyError):
        checks["reject_asymmetric_error"] = True
    else:
        checks["reject_asymmetric_error"] = False

    bad_family = deepcopy(geometry)
    bad_family["source_family"]["fixture_family"] = "different"
    try:
        leaf(c, tail, bad_family, symbolic, connection)
    except (TypeError, ValueError, KeyError):
        checks["reject_family_mismatch"] = True
    else:
        checks["reject_family_mismatch"] = False

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "small_rational_inputs": True,
        "no_old_module_called": True,
        "no_saved_data_read": True,
        "limitations": {
            "pressure": False, "collar": False, "full_word": False,
            "ledger": False, "source_owner": False,
        },
    }


run = run_fixture


__all__ = ["run_fixture", "run"]
