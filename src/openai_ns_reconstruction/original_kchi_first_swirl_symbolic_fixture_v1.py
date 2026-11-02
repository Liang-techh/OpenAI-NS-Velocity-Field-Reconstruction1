"""Synthetic structural checks for the direct W0_chi swirl symbolic leaf."""

from __future__ import annotations

from copy import deepcopy
from typing import Any


SCHEMA = "original_kchi_first_swirl_symbolic_fixture_v1"


def _endpoint(seed: int = 1) -> list[int]:
    return [0, seed, -8, seed.bit_length()]


def _interval(seed: int = 1) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": _endpoint(seed),
        "upper_exact_mpf_tuple": _endpoint(seed + 1),
        "directed": True,
    }


def _finite_term(field: str, order: str, source: int, target: int, seed: int) -> dict[str, Any]:
    source_id = f"fixture:kchi:{field}:{order}:{source}:{target}"
    return {
        "coefficient": _interval(seed),
        "formal_source_powers": [0, 1, 0, -1],
        "powers": [0, 1, 0, -1],
        "source_ids": [source_id],
        "origin": {"source_ids": [source_id], "error_classes": ["fixture_source"]},
        "source_cell": source,
        "source_index": source,
        "source_cell_id": f"fixture-cell:{source}",
        "source_term_id": source_id,
        "source_term_index": 0,
        "target_cell": target,
        "term_role": "finite_j1_signed",
        "error_only": False,
        "symmetric_error": False,
        "signed_field": True,
        "kernel_order": 1,
        "principal_kernel": "(-2*chi) T K0",
    }


def _tail_term(field: str, order: int, seed: int) -> dict[str, Any]:
    source_id = f"fixture:tail:{field}:{order}"
    return {
        "coefficient": {
            "lower_exact_mpf_tuple": [1, seed, -8, seed.bit_length()],
            "upper_exact_mpf_tuple": [0, seed, -8, seed.bit_length()],
            "directed": True,
        },
        "powers": [0, 1, 0, -1],
        "source_ids": [source_id],
        "error_classes": ["Kchi_j_ge_2_kernel_tail", "fixture_tail"],
        "uncertainty": True,
        "source_order": order,
        "kernel_Z_order": order,
        "physical_source_width_bound": 2,
        "error_only": True,
        "symmetric_error": True,
        "signed_field": False,
    }


def _kchi_binding() -> dict[str, Any]:
    cells = []
    for target in range(1):
        source_cells = []
        for source in range(target + 1):
            ordinary = {
                order: {
                    "F1": [_finite_term("F1", order, source, target, 10 + order_index)],
                    "d_x_F1": [_finite_term("d_x_F1", order, source, target, 20 + order_index)],
                }
                for order_index, order in enumerate(("Z0", "Z1", "Z2"))
            }
            source_cells.append(
                {
                    "source_index": source,
                    "target_index": target,
                    "triangular": True,
                    "source_x": _interval(2),
                    "target_x": _interval(3),
                    "ordinary_Z": ordinary,
                }
            )
        cells.append(
            {
                "index": target,
                "target_cell": target,
                "source_indices": [0],
                "triangular_diagonal_source": 0,
                "x": _interval(3),
                "target_x": _interval(3),
                "rho": _interval(4),
                "source_cells": source_cells,
            }
        )
    tail = {
        field: {order: [_tail_term(field, order_index, 30 + order_index)] for order_index, order in enumerate(("Z0", "Z1", "Z2"))}
        for field in ("F1", "d_x_F1")
    }
    family = {"actual_five_defect_family_sha256": "a" * 64}
    logs = {"hb": _interval(5), "F0_anchor": _interval(6), "Pstar_squared": _interval(7), "Lambda": _interval(8)}
    return {
        "schema": "original_kchi_first_saved_cells_v2",
        "source_family": family,
        "source_logs": logs,
        "formal_source_basis": ["hb", "F0_anchor", "Pstar_squared", "Lambda"],
        "finite_j1_correction": True,
        "finite_correction_is_signed": True,
        "tail_is_error_only": True,
        "cells": cells,
        "symmetric_kernel_tail": tail,
    }


def _background(name: str, seed: int) -> dict[str, Any]:
    return {
        "interval": _interval(seed),
        "source_id": f"fixture:background:{name}",
        "source_path": "fixture/prior_k0_symbolic",
    }


def _k0_symbolic(kchi: dict[str, Any]) -> dict[str, Any]:
    family = deepcopy(kchi["source_family"])
    logs = deepcopy(kchi["source_logs"])
    phi0 = _background("Phi_Z0", 40)
    phi1 = _background("Phi_Z1", 42)
    old_term = {
        "background_factor_name": "Phi_Z0",
        "background_factor": phi0,
        "source_ids": ["fixture:old:w1"],
    }
    old_term_1 = {
        "background_factor_name": "Phi_Z1",
        "background_factor": phi1,
        "source_ids": ["fixture:old:w1:z1"],
    }
    return {
        "schema": "original_k0_symbolic_m15_cells_v1",
        "source_family": family,
        "source_logs": logs,
        "source_power_basis": ["hb", "F0_anchor", "Pstar_squared", "Lambda", "F0_local", "L_local"],
        "cells": [
            {
                "index": 0,
                "rho": _interval(4),
                "terms": {"swirl": {"Z0": [old_term], "Z1": [old_term_1, old_term]}},
            }
        ],
    }


def _call(leaf: Any, kchi: Any, old: Any) -> dict[str, Any]:
    function = getattr(leaf, "build", leaf)
    if not callable(function):
        raise TypeError("leaf must expose build")
    return function(kchi, old)


def _reject(leaf: Any, kchi: Any, old: Any) -> bool:
    try:
        _call(leaf, kchi, old)
    except (TypeError, ValueError, KeyError):
        return True
    return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run synthetic positive and negative structural checks."""

    kchi = _kchi_binding()
    old = _k0_symbolic(kchi)
    result = _call(leaf, kchi, old)
    z0 = result["cells"][0]["terms"]["Z0"]
    z1 = result["cells"][0]["terms"]["Z1"]
    branches = {term["branch"] for term in z1}
    checks = {
        "schema": result.get("schema") == "original_kchi_first_swirl_symbolic_v1",
        "one_plus_three_counts": len(z0) == 1 and len(z1) == 3,
        "z1_branch_names": branches == {
            "Z1_Phi1_F1_0",
            "Z1_L_local_Phi0_F1_0",
            "Z1_Phi0_F1_Z",
        },
        "z0_branch_name": z0[0]["branch"] == "Z0_Phi0_F1_0",
        "six_axis_power": z0[0]["formal_source_powers"] == [0, 1, 0, -3, 1, 0],
        "local_log_power": any(term["formal_source_powers"] == [0, 1, 0, -3, 1, 1] for term in z1),
        "signed_prefactor_and_physical_scales": all(
            term["signed_factor"] == -2
            and term["physical_scale_factors"] == ["epsilon^2=Lambda^-2", "F0_local", "rho", "-2"]
            for term in z0 + z1
        ),
        "ordinary_z1_background_branches": {
            term["branch"]: term["background_factor_name"] for term in z1
        } == {
            "Z1_Phi1_F1_0": "Phi_Z1",
            "Z1_L_local_Phi0_F1_0": "Phi_Z0",
            "Z1_Phi0_F1_Z": "Phi_Z0",
        },
        "coefficient_unmodified": z0[0]["coefficient"] == kchi["cells"][0]["source_cells"][0]["ordinary_Z"]["Z0"]["F1"][0]["coefficient"],
        "source_ids_preserved": z0[0]["source_ids"] == kchi["cells"][0]["source_cells"][0]["ordinary_Z"]["Z0"]["F1"][0]["source_ids"],
        "target_source_preserved": z0[0]["target_cell"] == 0 and z0[0]["source_cell"] == 0,
        "rho_same_cell": z0[0]["rho_factor"] == old["cells"][0]["rho"],
        "tail_separate": result["j_ge_2_tail_is_error_only"] and result["tail_is_not_merged_with_finite"],
        "tail_not_signed": all(
            term["signed_field"] is False and term["error_only"] is True
            for field in result["symmetric_kernel_tail"].values()
            for rows in field.values()
            for term in rows
        ),
        "dx_and_z2_reserved": result["d_x_F1_retained_for_feedback"] and result["Z2_retained_for_feedback"],
        "conditional_owner": result["source_identity_verified"] is False and result["owner_verified"] is False,
        "no_ledger_or_numeric_ops": result["appended_to_212_ledger"] is False and result["no_interval_arithmetic"] is True and result["no_integration"] is True,
    }

    bad_family = deepcopy(old)
    bad_family["source_family"]["actual_five_defect_family_sha256"] = "b" * 64
    checks["family_mismatch_rejected"] = _reject(leaf, kchi, bad_family)

    bad_tail = deepcopy(kchi)
    bad_tail["symmetric_kernel_tail"]["F1"]["Z0"][0]["signed_field"] = True
    checks["signed_tail_rejected"] = _reject(leaf, bad_tail, old)

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "source_values_computed": False,
        "interval_arithmetic_performed": False,
        "integration_performed": False,
        "implementation_imported": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
