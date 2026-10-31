"""Independent synthetic fixture for the local-g5 physical-m15 leaf.

The fixture does not import the implementation.  The caller supplies the
``build`` callable and its directed interval context.  The synthetic source
records deliberately contain signed finite J1 intervals, symmetric J2 error
intervals, nontrivial Z1 rows, distinct source ids/paths, and nonzero source
widths.  Checks are enclosure checks; no interval midpoint is used.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Callable


__all__ = ["run_fixture", "run"]


def _is_interval(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _point_from_raw(c: Any, raw: Any) -> Any:
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((raw, raw))
        except Exception:
            pass
    for pair in ((raw, raw), [raw, raw]):
        try:
            return c.mpf(pair)
        except Exception:
            pass
    return c.mpf(raw)


def _endpoint(c: Any, value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return _point_from_raw(c, raw[1 if upper else 0])
    if hasattr(value, "a") and hasattr(value, "b"):
        candidate = value.b if upper else value.a
        nested_raw = getattr(candidate, "_mpi_", None)
        if nested_raw is not None:
            return _point_from_raw(c, nested_raw[1 if upper else 0])
        return c.mpf(candidate)
    return c.mpf(value)


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    if _is_interval(lower):
        lower = _endpoint(c, lower, False)
    if _is_interval(upper):
        upper = _endpoint(c, upper, True)
    return c.mpf([lower, upper])


def _same_interval(c: Any, left: Any, right: Any) -> bool:
    left_raw = getattr(left, "_mpi_", None)
    right_raw = getattr(right, "_mpi_", None)
    if left_raw is not None and right_raw is not None:
        return left_raw == right_raw
    return bool(
        _endpoint(c, left, False) == _endpoint(c, right, False)
        and _endpoint(c, left, True) == _endpoint(c, right, True)
    )


def _contains(c: Any, outer: Any, inner: Any) -> bool:
    return bool(
        _endpoint(c, outer, False) <= _endpoint(c, inner, False)
        and _endpoint(c, outer, True) >= _endpoint(c, inner, True)
    )


def _symmetric(c: Any, value: Any) -> bool:
    return bool(_endpoint(c, value, False) == -_endpoint(c, value, True))


def _width_positive(c: Any, value: Any) -> bool:
    return bool(_endpoint(c, value, True) > _endpoint(c, value, False))


def _sector(
    c: Any,
    coefficient: Any,
    powers: list[int],
    source_id: str,
    source_path: str,
    error_classes: list[str],
    *,
    error_only: bool,
    q_endpoint: Any | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "coefficient": coefficient,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [source_path],
        "error_classes": list(error_classes),
        "error_only": error_only,
    }
    if q_endpoint is not None:
        result["q_endpoint"] = q_endpoint
    return result


def _synthetic_inputs(c: Any) -> tuple[Any, Any, Any, Any, Any]:
    point_one = c.mpf(1)
    q_one = _interval(c, point_one, point_one)
    lead = {
        "Z0": [_sector(
            c, _interval(c, c.mpf(2), c.mpf(3)), [0, 1, 0, 0],
            "lead-F-Z0", "leading/F1/Z0", ["lead-z0"], error_only=False)],
        "Z1": [_sector(
            c, _interval(c, -c.mpf(1), -(c.mpf(1)/c.mpf(2))), [1, 0, 1, 0],
            "lead-F-Z1", "leading/F1/Z1", ["lead-z1"], error_only=False)],
    }
    # J1 has a negative Z0 endpoint and a positive Z1 endpoint.  The phase
    # hull therefore has a genuine q=0 endpoint and is not a point value.
    j1 = {
        "Z0": [_sector(
            c, _interval(c, -c.mpf(3), -c.mpf(1)), [0, 0, 0, 1],
            "j1-F1-Z0", "saved/local-g5/J1/Z0", ["j1-z0"],
            error_only=False, q_endpoint=q_one)],
        "Z1": [_sector(
            c, _interval(c, c.mpf(1), c.mpf(4)), [0, 1, 0, 1],
            "j1-F1-Z1", "saved/local-g5/J1/Z1", ["j1-z1"],
            error_only=False, q_endpoint=q_one)],
    }
    # J2 is already a whole-collar symmetric error.  It must remain symmetric
    # after the signed physical factor and source-product operation.
    j2 = {
        "Z0": [_sector(
            c, _interval(c, -c.mpf(2), c.mpf(2)), [0, 0, 1, 0],
            "j2-F1-Z0", "saved/local-g5/J2/Z0", ["j2-z0"], error_only=True)],
        "Z1": [_sector(
            c, _interval(c, -c.mpf(1), c.mpf(1)), [1, 0, 0, 1],
            "j2-F1-Z1", "saved/local-g5/J2/Z1", ["j2-z1"], error_only=True)],
    }
    sc = _interval(c, c.mpf(2), c.mpf(2))
    hb = _interval(c, (c.mpf(1)/c.mpf(4)), (c.mpf(1)/c.mpf(4)))
    certificate = {
        "hb": hb,
        "q_domain": _interval(c, c.mpf(0), c.mpf(1)),
        "source_chi_g_full_q_enclosure": True,
        "unified_exponent_bound": True,
        "geometry_Z_independent": True,
        "original_cutoff_mass_5_over_8": True,
        "chi_g_proof": "synthetic full-q directed enclosure",
        "exponent_proof": "synthetic one-envelope bound",
        "eta": _interval(c, c.mpf(5) / c.mpf(8), c.mpf(5) / c.mpf(8)),
    }
    return j1, j2, lead, sc, certificate


def run_fixture(build: Callable[..., Mapping[str, Any]], c: Any) -> dict[str, Any]:
    """Run independent synthetic checks against the supplied ``build``."""

    j1, j2, leading, sc, certificate = _synthetic_inputs(c)
    checks: dict[str, bool] = {}
    try:
        result = build(c, sc, j1, j2, leading, certificate)
    except Exception as exc:  # The report remains useful to the caller.
        return {
            "passed": False,
            "check_count": 1,
            "checks": {"build_completed": False},
            "error": repr(exc),
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["groups_are_separate"] = result.get("groups_separate") is True
    checks["no_chi_reapplication"] = result.get("chi_reapplied") is False
    checks["no_green_reapplication"] = result.get("green_reapplied") is False
    checks["no_ledger"] = result.get("ledger_212_written") is False
    checks["no_same_source_equality_claim"] = result.get("same_source_equality_claim") is False
    checks["source_only"] = result.get("source_only") is True
    checks["point_values_not_selected"] = result.get("point_values_selected") is False

    geometry = result.get("geometry", {})
    expected_exp = c.exp(sc / c.mpf(4))
    expected_exp_envelope = _interval(c, c.mpf(1), _endpoint(c, expected_exp, True))
    expected_rho2 = _interval(c, c.mpf(16), c.mpf(16) * _endpoint(c, expected_exp, True))
    expected_physical = -c.mpf(5) * sc * expected_exp_envelope
    checks["rho2_volume_bound"] = _same_interval(c, geometry.get("rho2_interval"), expected_rho2)
    checks["physical_factor_exact_scale"] = _same_interval(
        c, geometry.get("physical_factor"), expected_physical
    )
    checks["physical_factor_is_negative"] = bool(
        _endpoint(c, geometry.get("physical_factor"), True) < c.mpf(0)
    )
    checks["formal_power_shift"] = geometry.get("formal_power_shift") == [1, 0, 0, -2]

    for family, source in (("J1", j1), ("J2", j2)):
        group = result.get(family, {})
        checks[f"{family}_has_Z0"] = bool(group.get("Z0"))
        checks[f"{family}_has_Z1"] = len(group.get("Z1", ())) == 2
        for term in group.get("Z0", ()) + group.get("Z1", ()):
            checks[f"{family}_term_has_ids_{id(term)}"] = bool(term.get("source_ids"))
            checks[f"{family}_term_has_paths_{id(term)}"] = bool(term.get("source_paths"))
            checks[f"{family}_term_has_errors_{id(term)}"] = bool(term.get("error_classes"))
            checks[f"{family}_term_no_chi_{id(term)}"] = term.get("chi_reapplied") is False
            checks[f"{family}_term_no_green_{id(term)}"] = term.get("green_reapplied") is False
            checks[f"{family}_term_nonpoint_{id(term)}"] = _width_positive(c, term["coefficient"])
        # All expected source labels must survive at least once in the family.
        ids = [source[z][0]["source_ids"][0] for z in ("Z0", "Z1")]
        observed = [value for term in group.get("Z0", ()) + group.get("Z1", ())
                    for value in term.get("source_ids", ())]
        checks[f"{family}_source_ids_retained"] = all(item in observed for item in ids)
        paths = [source[z][0]["source_paths"][0] for z in ("Z0", "Z1")]
        observed_paths = [value for term in group.get("Z0", ()) + group.get("Z1", ())
                          for value in term.get("source_paths", ())]
        checks[f"{family}_source_paths_retained"] = all(item in observed_paths for item in paths)

    # J1 hulls contain the q=0 source, retain a finite signed endpoint, and
    # leave the positive-kernel phase as an interval rather than a sample.
    j1_terms = result["J1"]["Z0"] + result["J1"]["Z1"]
    checks["J1_q_endpoint_is_one"] = _same_interval(
        c, result["J1_contract"]["q_endpoint"], _interval(c, c.mpf(1), c.mpf(1))
    )
    checks["J1_q_zero_is_in_hull"] = result["J1_contract"].get("q_zero_is_in_hull") is True
    checks["J1_hull_flags"] = all(term.get("q_hull_applied") is True for term in j1_terms)
    checks["J1_signed_finite_flags"] = all(
        term.get("finite_signed_interval") is True
        and term.get("point_value_selected") is False
        for term in j1_terms
    )
    checks["J1_hulls_contain_zero"] = all(
        _endpoint(c, term["field_F1_coefficient"], False) <= c.mpf(0)
        <= _endpoint(c, term["field_F1_coefficient"], True)
        for term in j1_terms
    )
    checks["J1_signed_physical_nontrivial"] = any(
        bool(_endpoint(c, term["field_F1_coefficient"], True) == c.mpf(0))
        for term in j1_terms
    )

    # Every J2 output remains a symmetric error_only enclosure, despite the
    # negative physical factor and the signed leading-F intervals.
    j2_terms = result["J2"]["Z0"] + result["J2"]["Z1"]
    checks["J2_symmetric"] = all(
        term.get("error_only") is True
        and term.get("symmetric_error") is True
        and _symmetric(c, term["coefficient"])
        for term in j2_terms
    )
    checks["J2_not_hulled"] = all(term.get("q_hull_applied") is False for term in j2_terms)

    # Z1 must expose both product paths rather than one merged sector.
    j1_z1_kinds = {term.get("derivative_kind") for term in result["J1"]["Z1"]}
    j2_z1_kinds = {term.get("derivative_kind") for term in result["J2"]["Z1"]}
    checks["J1_Z1_all_paths"] = j1_z1_kinds == {
        "Z1_leading_Z1_times_field_Z0", "Z1_leading_Z0_times_field_Z1"
    }
    checks["J2_Z1_all_paths"] = j2_z1_kinds == {
        "Z1_leading_Z1_times_field_Z0", "Z1_leading_Z0_times_field_Z1"
    }
    checks["J1_Z0_power_scale"] = result["J1"]["Z0"][0]["powers"] == [1, 1, 0, -1]
    checks["J1_Z1_power_scale"] = {
        tuple(term.get("powers", ())) for term in result["J1"]["Z1"]
    } == {(2, 0, 1, -1), (1, 2, 0, -1)}

    # Independent numerical cross-check: a whole-phase q^4 field with a
    # nonconstant rho and a symmetric linear synthetic cutoff of mass 5/8.
    import mpmath as mp
    mp.mp.dps = 80
    q = lambda x: mp.mpf(x)
    weight = lambda x: 1 if x <= q("0.5") else 3 - 4*x
    rho2 = lambda x: 16*mp.exp(x/4)
    exact = -2*q("0.25")*(q(2)/4)*q("2.5")*(-q(2))*mp.quad(lambda x: weight(x)*rho2(x)*x**4, [0,q("0.5"),q("0.75")])
    interval = result["J1"]["Z0"][0]["coefficient"]
    lo, hi = mp.make_mpf(interval._mpi_[0])*q("0.25"), mp.make_mpf(interval._mpi_[1])*q("0.25")
    checks["independent_nonconstant_phase_integral"] = lo <= exact <= hi
    closed_certificate = dict(certificate, hb=_interval(c,c.mpf(0),c.mpf(1)/2), strict_hb_admitted=True)
    closed_result = build(c, sc, j1, j2, leading, closed_certificate)
    checks["closed_hb_is_formal_enclosure"] = bool(closed_result["geometry"]["hb"].a == c.mpf(0) and closed_result["geometry"]["hb"].b == c.mpf(1)/2)

    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "implementation_imported": False,
        "ancestor_called": False,
        "point_values_selected": False,
        "q_hull_is_conditional_on_certificate": True,
    }


run = run_fixture

