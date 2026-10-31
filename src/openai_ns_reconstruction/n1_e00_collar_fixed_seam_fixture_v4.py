"""Independent synthetic fixture for the version-three fixed-E00 leaf.

All public rows deliberately use the explicit ``Z0`` and ``Z1`` spelling.
The fixture passes the leaf as an argument and has no project or older-leaf
imports.  Its sentinel is placed in ``d_x_F1/Z1``, which is consumed by the
F-from-Fx and Fx-from-Fx first-potential branches.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping

FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
CANONICAL = ("F", "U", "K", "P", "Fx", "Ux")
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
BRANCHES = ("F_from_Fc", "F_from_Fxc", "Fx_from_Fc", "Fx_from_Fxc")
FIRST_HB_POWER_ADDED = {"F_from_Fc": 2, "F_from_Fxc": 3, "Fx_from_Fc": 1, "Fx_from_Fxc": 2}
Q = Fraction
SENTINEL_CLASS = "fixture:v4:E00:sentinel_error"
SENTINEL_ID = "fixture:v4:E00:sentinel:d_x_F1:Z1"


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(c: Any, value: Fraction | int | Any, source_id: str, powers: tuple[int, int, int, int], *, path: str, error_only: bool = False, error_class: str = "fixture:v4:signed") -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    result = {"coefficient": coefficient, "powers": list(powers), "source_ids": [source_id], "source_error_classes": [error_class], "error_classes": [error_class], "uncertainty": error_only, "origin": {"operation": "independent_fixture_v4", "source_ids": [source_id], "source_error_classes": [error_class], "error_classes": [error_class], "term_path": path}, "term_path": path, "old_full_prefix_ref": "fixture:v4:old_full_prefix:E00"}
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _packet(c: Any, prefix: str, values: Mapping[str, tuple[Fraction, Fraction]], powers: Mapping[str, tuple[int, int, int, int]]) -> dict[str, Any]:
    packet: dict[str, Any] = {}
    for field, canonical in zip(FIELDS, CANONICAL):
        z0, z1 = values[canonical]
        pows = powers[canonical]
        packet[field] = {"Z0": [_term(c, z0, f"{prefix}:{field}:Z0", pows, path=f"{prefix}.fields.{field}.Z0")], "Z1": [_term(c, z1, f"{prefix}:{field}:Z1", pows, path=f"{prefix}.fields.{field}.Z1")]}
    return packet


def _zero_interval(value: Any) -> bool:
    bits = value._mpi_
    return bits[0][1] == 0 and bits[1][1] == 0


def _nonzero_interval(value: Any) -> bool:
    bits = value._mpi_
    return bits[0][1] != 0 or bits[1][1] != 0


def _same_interval(left: Any, right: Any) -> bool:
    return getattr(left, "_mpi_", None) == getattr(right, "_mpi_", None)


def _first_terms(result: Mapping[str, Any], branch: str, order: str, cell: int = 2) -> list[Mapping[str, Any]]:
    return result["q_cells"][cell]["first_potential"]["branches"][branch][order]


def _negative_guards(build: Any, c: Any, inlet: Mapping[str, Any], chi: Mapping[str, Any], atoms: Mapping[str, str], certificate: Mapping[str, Any], q_cells: Any) -> dict[str, bool]:
    checks: dict[str, bool] = {}
    missing_owner = dict(certificate)
    owners = dict(missing_owner["owners"])
    owners.pop("chi", None)
    missing_owner["owners"] = owners
    try:
        build(c, inlet=inlet, chi=chi, source_family="fixture:v4:E00", source_atoms=atoms, certificate=missing_owner, q_cells=q_cells)
    except (TypeError, ValueError):
        checks["missing_chi_owner_rejected"] = True
    else:
        checks["missing_chi_owner_rejected"] = False

    bad_id = {field: {order: list(rows) for order, rows in inlet[field].items()} for field in FIELDS}
    bad_id["F1"]["Z0"] = [dict(bad_id["F1"]["Z0"][0], source_ids=[])]
    try:
        build(c, inlet=bad_id, chi=chi, source_family="fixture:v4:E00", source_atoms=atoms, certificate=certificate, q_cells=q_cells)
    except (TypeError, ValueError):
        checks["missing_source_id_rejected"] = True
    else:
        checks["missing_source_id_rejected"] = False

    bad_power = {field: {order: list(rows) for order, rows in inlet[field].items()} for field in FIELDS}
    bad_power["K1"]["Z1"] = [dict(bad_power["K1"]["Z1"][0], powers=[4, 5, 6])]
    try:
        build(c, inlet=bad_power, chi=chi, source_family="fixture:v4:E00", source_atoms=atoms, certificate=certificate, q_cells=q_cells)
    except (TypeError, ValueError):
        checks["bad_formal_powers_rejected"] = True
    else:
        checks["bad_formal_powers_rejected"] = False

    try:
        build(c, inlet=inlet, chi=chi, source_family="fixture:v4:E00", source_atoms=atoms, certificate=certificate, q_cells=((-Q(1, 13), Q(1, 4)), (Q(1, 4), Q(1, 2)), (Q(1, 2), 1)))
    except (TypeError, ValueError):
        checks["negative_q_cell_rejected"] = True
    else:
        checks["negative_q_cell_rejected"] = False
    return checks


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run structural, public-key, source-lineage and interval checks."""
    build = leaf.build if hasattr(leaf, "build") else leaf
    if not callable(build):
        raise TypeError("leaf must be a build callable or expose build")

    signed_powers = {
        "F": (-11, 3, 2, 0),
        "U": (6, -4, 0, 0),
        "K": (4, 0, -5, 0),
        "P": (7, 4, 0, -4),
        "Fx": (-4, 0, 4, 3),
        "Ux": (5, 0, -5, 3),
    }
    signed_values = {
        "F": (Q(-31, 109), Q(37, 113)),
        "U": (Q(41, 127), Q(-43, 131)),
        "K": (Q(-47, 137), Q(53, 139)),
        "P": (Q(59, 149), Q(-61, 151)),
        "Fx": (Q(67, 157), Q(-71, 163)),
        "Ux": (Q(-73, 167), Q(79, 173)),
    }
    inlet = _packet(c, "fixture:v4:E00:signed", signed_values, signed_powers)
    inlet["K1"]["Z0"].append(_term(c, Q(-17, 179), "fixture:v4:E00:signed:K1:Z0:extra", signed_powers["K"], path="fixture.v4.E00.fields.K1.Z0.extra"))
    sentinel = _term(c, _box(c, Q(-1, 1013), Q(1, 1013)), SENTINEL_ID, signed_powers["Fx"], path="fixture.v4.E00.fields.d_x_F1.Z1.sentinel", error_only=True, error_class=SENTINEL_CLASS)
    inlet["d_x_F1"]["Z1"].append(sentinel)

    chi = {
        "Z0": [_term(c, Q(-83, 181), "fixture:v4:chi:Z0", (0, 0, 0, 0), path="fixture.v4.chi.Z0", error_class="fixture:v4:chi")],
        "Z1": [_term(c, Q(89, 191), "fixture:v4:chi:Z1", (0, 0, 0, 0), path="fixture.v4.chi.Z1", error_class="fixture:v4:chi")],
    }
    q_cells = ((0, Q(1, 4)), (Q(1, 4), Q(1, 2)), (Q(1, 2), 1))
    certificate = {
        "owners": {"q": "fixture:v4:actual:q", "x": "fixture:v4:actual:x", "rho": "fixture:v4:actual:rho", "chi": "fixture:v4:actual:chi"},
        "fixed_geometry": True,
        "ordinary_Z_geometry_fixed": True,
        "rho_radius_Z_independent": True,
        "hb": Q(2, 11),
        "sc": Q(11, 193),
        "q_cells": q_cells,
        "fixture_geometry_owner": "independent_fixed_E00_seam_v4",
    }
    atoms = {name: "fixture:v4:E00:" + name for name in SOURCE_POWER_BASIS}
    result = build(c, inlet=inlet, chi=chi, source_family="synthetic_fixed_E00_v4", source_atoms=atoms, certificate=certificate, q_cells=q_cells)

    checks: dict[str, bool] = {
        "schema": result.get("schema") == "n1_e00_collar_fixed_seam_v4",
        "fields": tuple(result.get("fields", ())) == FIELDS,
        "ordinary_Z_orders": tuple(result.get("ordinary_Z_orders", ())) == (0, 1),
        "three_dyadic_q_cells": len(result.get("q_cells", ())) == 3,
        "source_power_basis": tuple(result.get("source_power_basis", ())) == SOURCE_POWER_BASIS,
        "public_inlet_Z_keys": all(tuple(result["inlet"][field]) == ("Z0", "Z1") for field in FIELDS),
        "public_chi_Z_keys": tuple(result["chi_jets"]) == ("Z0", "Z1"),
        "fixed_target": result.get("geometry", {}).get("fixed_target") == "x=2,rho=4",
        "fixed_geometry": result.get("geometry", {}).get("ordinary_Z_geometry_fixed") is True and result.get("geometry", {}).get("ordinary_Z_geometry_differentiated") is False,
        "radial_derivative_preserved": result.get("physical_radial_derivative_preserved") is True,
        "source_terms_indexed": all(result["source_term_ledger"][field][order] for field in FIELDS for order in ("Z0", "Z1")),
        "source_paths_retained": all(entry["path_refs"] and entry["old_full_prefix_refs"] for field in FIELDS for order in ("Z0", "Z1") for entry in result["source_term_ledger"][field][order]),
        "all_input_source_ids_indexed": {item["source_ids"][0] for field in FIELDS for order in ("Z0", "Z1") for item in inlet[field][order]} <= {item["source_ids"][0] for field in FIELDS for order in ("Z0", "Z1") for item in result["source_term_ledger"][field][order]},
        "seam_componentwise_zero": result["seam"]["h0_identity"]["componentwise_zero"] is True,
        "seam_literal_zero": result["seam"]["h0_identity"]["literal_directed_arithmetic"] is True and result["seam"]["conditional_algebraic_seam_only"] is True,
        "seam_masses_zero": result["seam"]["all_potential_masses_zero"] is True,
        "seam_local_and_tail_zero": result["seam"]["local_forcing_exact_zero"] is True and result["seam"]["higher_potential_words_exact_zero"] is True,
        "no_old_W1_equality": result["operator_equality_to_old_W1_proved"] is False and result["old_W1_function_equal"] is False,
        "no_duplicate_claim": result["no_dup_across_unrelated_owners_proved"] is False and result["append_to_old_coarse_prefix"] is False,
        "K0_algebraic_only": result["K0_transfer_algebraic_only"] is True and result["legacy_coarse_W1_transfer_reused"] is False,
        "q_cell_K0_six_component_transfer": all(result["q_cells"][2]["k0_transfer"]["branches"][branch]["Z0"] and result["q_cells"][2]["k0_transfer"]["branches"][branch]["Z1"] for branch in ("F_from_Fc", "F_from_Fxc", "U_from_Uc", "U_from_Uxc", "K_from_Kc", "K_from_Uxc", "P_from_Pc", "Fx_from_Fxc", "Ux_from_Uxc")),
        "P1_separate_no_m15": "P1" in result["fields"] and result.get("m15_computed") is False,
        "j1_separate_j2_omitted": result["K0_j1_disjoint"] is True and result["j2_remainder_included"] is False,
    }

    seam_k0 = result["seam"]["k0_transfer"]["branches"]
    identity_routes = {"F_from_Fc": "F1", "U_from_Uc": "Uz1", "K_from_Kc": "K1", "P_from_Pc": "P1", "Fx_from_Fxc": "d_x_F1", "Ux_from_Uxc": "d_x_Uz1"}
    zero_routes = {"F_from_Fxc": "d_x_F1", "U_from_Uxc": "d_x_Uz1", "K_from_Uxc": "d_x_Uz1"}
    for branch, field in identity_routes.items():
        checks["identity_route_" + branch] = all(len(seam_k0[branch][order]) == len(result["inlet"][field][order]) and all(_same_interval(got["coefficient"], expected["coefficient"]) and got["identity_difference_exact_zero"] is True for got, expected in zip(seam_k0[branch][order], result["inlet"][field][order])) for order in ("Z0", "Z1"))
    for branch in zero_routes:
        checks["zero_route_" + branch] = all(_zero_interval(term["coefficient"]) and term["identity_difference_exact_zero"] is True for order in ("Z0", "Z1") for term in seam_k0[branch][order])
    checks["seam_radial_derivative_factors"] = all(seam_k0[branch]["Z0"][0]["physical_radial_derivative_preserved"] for branch in ("Fx_from_Fxc", "Ux_from_Uxc")) and all("q_derivative" not in data for data in result["seam"]["h0_identity"]["factors"].values())
    checks["seam_first_potential_zero"] = all(_zero_interval(term["coefficient"]) for branch in BRANCHES for order in ("Z0", "Z1") for term in result["seam"]["first_potential"]["branches"][branch][order])

    checks["nonzero_q_masses"] = all(_nonzero_interval(result["q_cells"][2]["potential_kernel_masses"][branch]["mass"]) for branch in BRANCHES)
    checks["first_potential_signed"] = any(term.get("signed_source_term") and not term.get("error_only") and _nonzero_interval(term["coefficient"]) for branch in BRANCHES for term in _first_terms(result, branch, "Z0"))
    checks["ordinary_Z1_product_rule"] = all({term.get("leibniz_piece") for branch in BRANCHES for term in _first_terms(result, branch, "Z1")} >= {"chi_Z*inlet_Z0", "chi*inlet_Z1"} for _ in (0,))
    sentinel_terms = [term for branch in ("F_from_Fxc", "Fx_from_Fxc") for term in _first_terms(result, branch, "Z1") if SENTINEL_ID in term.get("incoming_source_ids", ())]
    checks["sentinel_error_preserved"] = any(SENTINEL_CLASS in term.get("source_error_classes", ()) and SENTINEL_CLASS in term.get("error_classes", ()) and term.get("sentinel_error_propagated") is True and term.get("error_only") is True and term.get("symmetric_error") is True for term in sentinel_terms)
    checks["sentinel_source_id_preserved"] = any(SENTINEL_ID in term.get("incoming_source_ids", ()) and SENTINEL_ID in term.get("source_ids", ()) for term in sentinel_terms)
    checks["signed_negative_input_retained"] = result["inlet"]["F1"]["Z0"][0]["coefficient"]._mpi_[0][0] == 1 and any(term.get("signed_source_term") for term in _first_terms(result, "F_from_Fc", "Z0"))
    expected_sources = {"F_from_Fc": "F1", "F_from_Fxc": "d_x_F1", "Fx_from_Fc": "F1", "Fx_from_Fxc": "d_x_F1"}
    checks["formal_first_powers"] = all(next(term["powers"][0] - result["inlet"][expected_sources[branch]]["Z0"][0]["powers"][0] == FIRST_HB_POWER_ADDED[branch] for term in _first_terms(result, branch, "Z0")) for branch in BRANCHES)
    checks["j2_separate_not_claimed"] = all(result["q_cells"][2]["j2_tail"]["j2_remainder_included"] is False and result["q_cells"][2]["j2_tail"]["j2_bound_status"] == "not_computed_in_this_leaf" and not result["q_cells"][2]["j2_tail"]["branches"][branch][order] for branch in BRANCHES for order in ("Z0", "Z1"))
    checks["q_positive_local_omitted"] = all(result["q_cells"][index]["local_forcing"]["local_forcing_status"] == "not_computed" and result["q_cells"][index]["local_forcing"]["local_source_computed"] is False for index in (0, 1, 2))
    checks.update(_negative_guards(build, c, inlet, chi, atoms, certificate, q_cells))
    return {"passed": all(checks.values()), "checks": checks, "check_count": len(checks), "synthetic_only": True, "independent_exact_fraction_inputs": True, "directed_interval_arithmetic": True, "signed_six_fields": True, "sentinel_error_checked_exact_label": SENTINEL_CLASS, "public_Z_keys_checked": True, "seam_radial_derivative_factors_checked": True, "quadrature_used": False, "old_fixture_imported": False, "project_imports": False, "receipt_reads": False, "midpoint_values_used": False, "source_identity_proved": False, "operator_equality_to_old_W1_proved": False, "limits": {"seam": "conditional algebraic H0 identity only", "old_prefix": "no equality or no-dup claim", "j2": "not computed in this leaf", "local_forcing": "q>0 not computed; q=0 exact zero"}}


run = run_fixture

__all__ = ["FIELDS", "CANONICAL", "SOURCE_POWER_BASIS", "BRANCHES", "SENTINEL_CLASS", "run_fixture", "run"]
