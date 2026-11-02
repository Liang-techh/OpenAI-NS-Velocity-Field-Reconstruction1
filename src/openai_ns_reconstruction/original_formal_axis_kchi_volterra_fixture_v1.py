"""Independent rational checks for the formal axis-Kchi Volterra leaf.

The fixture owns every source term and provenance label.  It does not import
the leaf or any ancestor producer: ``run_fixture(leaf.build, c)`` receives the
callable explicitly.  Its return value contains only JSON-safe scalars,
strings, lists and mappings.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F", "Fx")
SENTINEL_ID = "fixture:formal-axis:Fx:Z1:inherited-error"
SENTINEL_CLASS = "fixture:formal-axis:inherited-error"

__all__ = ["SOURCE_POWER_BASIS", "FIELDS", "SENTINEL_ID", "SENTINEL_CLASS", "run_fixture", "run"]


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / value.denominator
    return c.mpf(value)


def _box(c: Any, lower: Fraction | int, upper: Fraction | int) -> Any:
    lo, hi = _q(c, lower), _q(c, upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _contains(box: Any, value: Any) -> bool:
    return box.a <= value.a and box.b >= value.b


def _term(
    c: Any,
    value: Fraction | int | Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    path: str,
    error_only: bool = False,
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else _q(c, value)
    result = {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [path],
        "source_error_classes": [SENTINEL_CLASS if error_only else f"{path}:signed"],
        "error_classes": [SENTINEL_CLASS if error_only else f"{path}:signed"],
        "term_path": path,
        "origin": {"operation": "independent_formal_axis_fixture", "source_ids": [source_id], "term_path": path},
        "error_only": error_only,
        "symmetric_error": error_only,
        "uncertainty": error_only,
    }
    return result


def _rows(c: Any, prefix: str, values: tuple[Any, Any], powers: tuple[int, int, int, int], *, error_z1: bool = False) -> dict[str, list[dict[str, Any]]]:
    return {
        "Z0": [_term(c, values[0], f"{prefix}:Z0", powers, path=f"{prefix}.Z0")],
        "Z1": [_term(c, _box(c, -abs(values[1]), abs(values[1])) if error_z1 else values[1], SENTINEL_ID if error_z1 else f"{prefix}:Z1", powers, path=f"{prefix}.Z1", error_only=error_z1)],
    }


def _inlet(c: Any) -> dict[str, Any]:
    return {
        "F": _rows(c, "fixture:formal-axis:Fc", (Fraction(3, 5), Fraction(-2, 7)), (1, 2, 3, -1)),
        "Fx": _rows(c, "fixture:formal-axis:Fxc", (Fraction(5, 11), Fraction(1, 13)), (2, -1, 0, 2), error_z1=True),
    }


def _forcing(c: Any) -> dict[str, Any]:
    return {
        "uniform_collar": True,
        "scope": "whole collar h5 uniform sector",
        "Z0": [_term(c, Fraction(2, 9), "fixture:formal-axis:h5:Z0", (0, 1, 1, 0), path="h5.Z0")],
        "Z1": [_term(c, Fraction(-1, 6), "fixture:formal-axis:h5:Z1", (1, 0, 2, -1), path="h5.Z1")],
    }


def _axis(c: Any, value: Fraction | int, derivative: Fraction | int, *, id_value: str = "fixture:formal-axis:Kchi") -> tuple[dict[str, Any], dict[str, Any]]:
    common = {
        "powers": [0, 0, 0, 0],
        "source_ids": [id_value],
        "source_paths": ["axis.Kchi.whole_Z"],
        "whole_Z": True,
        "scope": "same-source axis Kchi whole-Z interval",
    }
    chi = dict(common)
    chi.update({"coefficient": _q(c, value), "coefficient_interval": _q(c, value), "term_path": "axis.Kchi.Z0"})
    chi_z = dict(common)
    chi_z.update({"coefficient": _q(c, derivative), "coefficient_interval": _q(c, derivative), "term_path": "axis.Kchi.Z1"})
    return chi, chi_z


def _find(rows: Mapping[str, Any], *, base: str, word: int, piece: str | None = None) -> list[Any]:
    result = []
    for term in rows["Z0"] + rows["Z1"]:
        if term.get("base") == base and term.get("word_order") == word and (piece is None or term.get("leibniz_piece") == piece):
            result.append(term)
    return result


def _all_terms(result: Mapping[str, Any], key: str = "finite_prefix") -> list[Any]:
    terms: list[Any] = []
    for field in FIELDS:
        for label in ("Z0", "Z1"):
            terms.extend(result[key][field][label])
    return terms


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent analytic checks without importing a producer."""

    if hasattr(build, "build"):
        build = build.build
    checks: dict[str, bool] = {}
    inlet = _inlet(c)
    forcing = _forcing(c)
    sc = _q(c, Fraction(1, 2))
    chi_zero, chi_z = _axis(c, 0, Fraction(1, 10))

    zero = build(c, inlet, forcing, chi_zero, chi_z, sc, max_order=0)
    checks["schema_and_scope"] = (
        zero.get("schema") == "original_formal_axis_kchi_volterra_v1"
        and zero.get("fields") == ["F", "Fx"]
        and zero["scope"]["full_six_field"] is False
        and zero["scope"]["time_recursion"] is False
        and zero["scope"]["m15"] is False
        and zero["scope"]["ledger_written"] is False
        and zero["scope"]["same_owner_proof"] is False
    )
    checks["zero_power_convention"] = zero["tail_majorant"]["zero_power_convention"] == "0^0=1 when M=0 and N=0"
    checks["geometry_contract"] = (
        zero["geometry"]["x_formula"] == "x(q)=2*exp(hb*sc*q/8), 0<=q<=1"
        and zero["geometry"]["L_bound"] == "L <= hb*S"
        and zero["geometry"]["K0_bound"] == "0 <= K0(x,t) <= x-t"
        and zero["geometry"]["Kx_bound"] == "0 <= Kx(x,t) <= 1"
        and zero["conditions"]["hb_materialized"] is False
    )

    finite_zero = zero["finite_prefix"]
    fc_base = [term for term in finite_zero["F"]["Z0"] if term.get("base") == "Fc"]
    fxc_base = [term for term in finite_zero["F"]["Z0"] if term.get("base") == "Fxc"]
    forcing_base = [term for term in finite_zero["F"]["Z0"] if term.get("base") == "forcing"]
    checks["n0_F_bases"] = (
        len(fc_base) == len(fxc_base) == len(forcing_base) == 1
        and _contains(fc_base[0]["coefficient"], _q(c, Fraction(3, 5)))
        and fc_base[0]["powers"][0] == 1
        and fxc_base[0]["powers"][0] == 3
        and forcing_base[0]["powers"][0] == 2
    )
    fxc_fx = [term for term in finite_zero["Fx"]["Z0"] if term.get("base") == "Fxc"]
    forcing_fx = [term for term in finite_zero["Fx"]["Z0"] if term.get("base") == "forcing"]
    expected_fxc_fx = _q(c, Fraction(5, 11)) * c.exp(-3 * sc / 16)
    checks["n0_Fx_bases"] = (
        len(fxc_fx) == len(forcing_fx) == 1
        and _contains(fxc_fx[0]["coefficient"], expected_fxc_fx)
        and fxc_fx[0]["powers"][0] == 2
        and forcing_fx[0]["powers"][0] == 1
    )
    checks["chi_zero_n0_has_no_finite_words"] = all(term.get("word_order") == 0 for term in _all_terms(zero))

    zero_f_z_tail = zero["tail"]["F"]["Z1"]
    insertion_tail = [term for term in zero_f_z_tail if term.get("word_kind") == "chi_Z_insertion_tail"]
    checks["chi_zero_nonzero_chi_Z_tail"] = (
        len(insertion_tail) == 3
        and all(term.get("tail") is True and term.get("error_only") is True and term.get("symmetric_error") is True for term in insertion_tail)
        and all("fixture:formal-axis:Kchi" in term.get("source_ids", []) for term in insertion_tail)
        and all(term["coefficient"].b > 0 for term in insertion_tail)
    )
    inherited = [term for term in finite_zero["Fx"]["Z1"] if term.get("base") == "Fxc"]
    checks["inherited_error_separate"] = (
        len(inherited) == 1
        and inherited[0]["inherited_error_only"] is True
        and inherited[0]["term_kind"] == "inherited_error"
        and inherited[0]["signed_finite"] is False
        and any(term.get("term_kind") == "volterra_tail" for term in zero["tail"]["Fx"]["Z1"])
    )

    chi_const, chi_z_const = _axis(c, Fraction(1, 8), Fraction(1, 16))
    nonzero = build(c, inlet, forcing, chi_const, chi_z_const, sc, max_order=1)
    fc_n1 = [term for term in nonzero["finite_prefix"]["F"]["Z0"] if term.get("base") == "Fc" and term.get("word_order") == 1]
    expected_fc_n1 = _q(c, Fraction(3, 5)) * (-2) * _q(c, Fraction(1, 8)) * nonzero["geometry"]["S"] ** 2 / 2
    checks["constant_chi_first_word"] = len(fc_n1) == 1 and _contains(fc_n1[0]["coefficient"], expected_fc_n1) and fc_n1[0]["powers"][0] == 3
    fc_fx_n1 = [term for term in nonzero["finite_prefix"]["Fx"]["Z0"] if term.get("base") == "Fc" and term.get("word_order") == 1]
    expected_fc_fx_n1 = _q(c, Fraction(3, 5)) * (-2) * _q(c, Fraction(1, 8)) * (nonzero["geometry"]["S"])
    checks["constant_chi_Fx_first_word"] = len(fc_fx_n1) == 1 and _contains(fc_fx_n1[0]["coefficient"], expected_fc_fx_n1) and fc_fx_n1[0]["powers"][0] == 2
    z1_terms = nonzero["finite_prefix"]["F"]["Z1"]
    checks["ordinary_Z1_two_pieces"] = (
        any(term.get("base") == "Fc" and term.get("word_order") == 1 and term.get("leibniz_piece") == "chi^n*base_Z1" for term in z1_terms)
        and any(term.get("base") == "Fc" and term.get("word_order") == 1 and term.get("leibniz_piece") == "n*chi^(n-1)*chi_Z*base_Z0" for term in z1_terms)
    )
    checks["source_power_addition"] = (
        fc_n1[0]["powers"][1:] == [2, 3, -1]
        and fc_n1[0]["source_ids"] == ["fixture:formal-axis:Fc:Z0", "fixture:formal-axis:Kchi"]
    )
    checks["signed_finite_and_formal_hb"] = (
        fc_n1[0]["signed_finite"] is True
        and fc_n1[0]["error_only"] is False
        and fc_n1[0]["coefficient_contains_formal_hb"] is False
        and fc_fx_n1[0]["coefficient_contains_formal_hb"] is False
    )
    tail_fc = [term for term in nonzero["tail"]["F"]["Z0"] if term.get("base") == "Fc"][0]
    true_n2 = _q(c, Fraction(3, 5)) * (_q(c, Fraction(1, 8)) ** 2) * 4 * nonzero["geometry"]["S"] ** 4 / 24
    checks["tail_contains_simple_constant_n2"] = _contains(tail_fc["coefficient"], true_n2) and tail_fc["word_order"] == 2 and tail_fc["tail"] is True
    checks["tail_starts_after_prefix"] = nonzero["tail_starts_at_j"] == 2 and all(term["word_order"] == 2 for term in nonzero["tail"]["F"]["Z0"])
    checks["tail_ids_and_error_flow"] = (
        "fixture:formal-axis:Kchi" in tail_fc["source_ids"]
        and tail_fc["error_only"] is True
        and tail_fc["symmetric_error"] is True
        and tail_fc["inherited_error_only"] is False
    )
    empty_forcing = build(c, inlet, {"Z0": [], "Z1": []}, chi_zero, chi_z, sc, max_order=0)
    checks["empty_forcing_means_no_local_forcing"] = (
        empty_forcing["conditions"]["forcing_rows_supplied"] is False
        and empty_forcing["conditions"]["fake_zero_forcing_filled"] is False
        and all(term.get("base") != "forcing" for term in _all_terms(empty_forcing))
    )

    dedup = build(c, inlet, forcing, chi_const, chi_z_const, sc, max_order=3, min_order=2)
    dedup_terms = _all_terms(dedup, "finite_prefix")
    checks["min_order_two_deduplicates"] = (
        dedup["min_order"] == 2
        and dedup["finite_orders_emitted"] == [2, 3]
        and all(term.get("word_order") in {2, 3} for term in dedup_terms)
        and dedup["alternative_envelopes"]["new_tail_is_alternative_to_prior_j2"] is True
        and dedup["alternative_envelopes"]["combined"] is False
        and dedup["alternative_envelopes"]["never_add_old_and_new_tail"] is True
        and dedup["branches"]["Fc"]["prior_orders_not_recomputed"] == [0, 1]
    )
    checks["dedup_tail_starts_N_plus_one"] = dedup["tail_starts_at_j"] == 4 and all(term["word_order"] == 4 for term in dedup["tail"]["F"]["Z0"])

    def rejected(call: Any) -> bool:
        try:
            call()
        except (TypeError, ValueError, KeyError):
            return True
        return False

    missing_fx = {"F": inlet["F"]}
    checks["incomplete_inlet_rejected"] = rejected(lambda: build(c, missing_fx, forcing, chi_const, chi_z_const, sc, max_order=0))
    checks["nonpositive_sc_rejected"] = rejected(lambda: build(c, inlet, forcing, chi_const, chi_z_const, _q(c, 0), max_order=0))
    checks["negative_order_rejected"] = rejected(lambda: build(c, inlet, forcing, chi_const, chi_z_const, sc, max_order=-1))
    checks["resource_order_rejected"] = rejected(lambda: build(c, inlet, forcing, chi_const, chi_z_const, sc, max_order=33))
    bad_error = _inlet(c)
    bad_error["Fx"]["Z1"][0]["error_only"] = True
    bad_error["Fx"]["Z1"][0]["symmetric_error"] = False
    checks["bad_error_flags_rejected"] = rejected(lambda: build(c, bad_error, forcing, chi_const, chi_z_const, sc, max_order=0))
    asymmetric_error = _inlet(c)
    asymmetric_error["Fx"]["Z1"][0]["coefficient"] = _box(c, -1, 2)
    checks["asymmetric_error_rejected"] = rejected(lambda: build(c, asymmetric_error, forcing, chi_const, chi_z_const, sc, max_order=0))
    bad_axis = dict(chi_const)
    bad_axis.update(error_only=True, symmetric_error=True, coefficient=_box(c, -1, 1))
    checks["error_only_axis_rejected"] = rejected(lambda: build(c, inlet, forcing, bad_axis, chi_z_const, sc, max_order=0))
    bad_scope, bad_scope_z = _axis(c, Fraction(1, 8), Fraction(1, 16))
    bad_scope["scope"] = "lower-collar shear"
    checks["lower_collar_chi_rejected"] = rejected(lambda: build(c, inlet, forcing, bad_scope, bad_scope_z, sc, max_order=0))
    bad_id, bad_id_z = _axis(c, Fraction(1, 8), Fraction(1, 16), id_value="other-axis")
    checks["different_axis_source_rejected"] = rejected(lambda: build(c, inlet, forcing, chi_const, bad_id_z, sc, max_order=0))

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_rational_sources": True,
        "json_safe": True,
        "receipt_or_ancestor_read": False,
        "midpoint_values_used": False,
        "cases": ["chi=0, chi_Z!=0, N=0", "constant nonzero chi, N=1", "min_order=2 deduplication", "negative input guards"],
        "limits": {
            "hb": "formal 0<hb<1/2 assumption; no hb materialization",
            "scope": "F/Fx only; not six-field, time-recursive, m15, ledger, or same-owner proof",
            "tail": "new tail and any prior j>=2 tail are alternative envelopes, never summed",
        },
    }


run = run_fixture

