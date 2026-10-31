"""Independent synthetic checks for ``n1_local_second_potential_m15_v1``.

The fixture never imports the implementation.  A caller passes the build
callable and a directed interval context.  It checks source identity,
directed signs, q=0 metadata, ordinary-Z route separation, formal powers,
chi=0 with nonzero chi_Z, same-source P negatives, old-tail retention, and a
nonconstant positive-Green quadrature against the j=2 mass enclosure.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Callable

__all__ = ["run_fixture", "run"]


def _raw(x: Any, upper: bool) -> Any:
    raw = getattr(x, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    if hasattr(x, "a") and hasattr(x, "b"):
        endpoint = x.b if upper else x.a
        nested = getattr(endpoint, "_mpi_", None)
        return nested[1 if nested is not None else 0] if nested is not None else endpoint
    return x


def _point(c: Any, x: Any) -> Any:
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((x, x))
        except Exception:
            pass
    return c.mpf(x)


def _iv(c: Any, lo: Any, hi: Any) -> Any:
    a, b = _raw(lo, False), _raw(hi, True)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((a, b))
        except Exception:
            pass
    return c.mpf([_point(c, a), _point(c, b)])


def _lo(c: Any, x: Any) -> Any:
    return _point(c, _raw(x, False))


def _hi(c: Any, x: Any) -> Any:
    return _point(c, _raw(x, True))


def _same(c: Any, left: Any, right: Any) -> bool:
    return bool(_lo(c, left) == _lo(c, right) and _hi(c, left) == _hi(c, right))


def _sym(c: Any, x: Any) -> bool:
    return bool(_lo(c, x) == -_hi(c, x))


def _negate(c: Any, x: Any) -> Any:
    return _iv(c, -_hi(c, x), -_lo(c, x))


def _sector(c: Any, coefficient: Any, powers: list[int], source_id: str, path: str, error: str) -> dict[str, Any]:
    return {"coefficient": coefficient, "powers": list(powers), "source_ids": [source_id], "source_paths": [path], "error_classes": [error], "error_only": False}


def _inputs(c: Any, *, nonzero_chi: bool = False) -> tuple[Any, Any, Any, Any, Any, Any, Any, Any]:
    one = c.mpf(1)
    q_domain = _iv(c, c.mpf(0), one)
    chi = _sector(c, _iv(c, c.mpf(1) / 4, c.mpf(1) / 3) if nonzero_chi else _iv(c, c.mpf(0), c.mpf(0)), [0, 0, 0, 0], "chi-axis-zero" if not nonzero_chi else "chi-axis-positive", "synthetic/chi/Z0", "chi-axis")
    chi_z = _sector(c, _iv(c, c.mpf(1) / 5, c.mpf(1) / 4), [0, 0, 0, 0], "chi-axis-Z1", "synthetic/chi/Z1", "chi-Z1")
    g5 = {
        "Z0": [_sector(c, _iv(c, c.mpf(1), c.mpf(2)), [0, 1, 0, 0], "g5-Z0", "synthetic/g5/Z0/whole-phase", "g5-value")],
        "Z1": [_sector(c, _iv(c, c.mpf(3) / 2, c.mpf(2)), [1, 0, 1, 0], "g5-Z1", "synthetic/g5/Z1/whole-phase", "g5-Z1")],
    }
    leading = {
        "Z0": [_sector(c, _iv(c, c.mpf(1), c.mpf(3) / 2), [0, 0, 1, 0], "lead-Z0", "synthetic/leading/Z0", "leading")],
        "Z1": [_sector(c, _iv(c, c.mpf(1) / 2, c.mpf(1)), [0, 1, 0, 0], "lead-Z1", "synthetic/leading/Z1", "leading-Z1")],
    }
    hb = _iv(c, c.mpf(1) / 4, c.mpf(1) / 4)
    sc = _iv(c, c.mpf(2), c.mpf(2))
    certificate = {
        "same_source": True,
        "chi_is_axis_potential": True,
        "geometry_Z_independent": True,
        "g5_whole_phase_enclosure": True,
        "P_sign": 2,
        "eta_mass": _iv(c, c.mpf(5) / 8, c.mpf(5) / 8),
        "q_domain": q_domain,
    }
    old_tail = {"schema": "old-j>=2", "source_ids": ["old-tail-id"], "error_classes": ["old-j2-error"], "retained": True}
    return chi, chi_z, g5, leading, hb, sc, certificate, old_tail


def _quadrature_check(c: Any, result: Mapping[str, Any], sc: Any, hb: Any) -> bool:
    """Compare a positive nonconstant Green integral with the mass cap.

    This is an independent numerical sanity check of the positive-kernel
    envelope, not a replay of the implementation's interval arithmetic.
    """

    try:
        import mpmath as mp
    except Exception:
        return False
    mp.mp.dps = 60
    def mp_endpoint(value: Any, upper: bool = False) -> Any:
        raw = getattr(value, "_mpi_", None)
        if raw is not None:
            return mp.make_mpf(raw[1 if upper else 0])
        return mp.mpf(value)

    sc_value = mp_endpoint(sc)
    hb_value = mp_endpoint(hb)
    B = (sc_value / 4) * mp.exp(sc_value / 16)
    green = mp.quad(lambda q: (hb_value * B * q) ** 6 * (1 + q) / mp.mpf(6 * 5 * 4 * 3 * 2 * 1), [0, 1])
    mass = result["geometry"]["j2_positive_mass_envelope"]
    return green <= mp_endpoint(mass, upper=True)


def run_fixture(build: Callable[..., Mapping[str, Any]], c: Any) -> dict[str, Any]:
    checks: dict[str, bool] = {}
    try:
        chi0, chi_z, g5, leading, hb, sc, certificate, old_tail = _inputs(c, nonzero_chi=False)
        result = build(c, hb, sc, chi0, chi_z, g5, leading, certificate, old_j2_tail=old_tail)
    except Exception as exc:
        return {"passed": False, "check_count": 1, "checks": {"build_completed": False}, "error": repr(exc), "synthetic_only": True, "implementation_imported": False}

    checks["build_completed"] = isinstance(result, Mapping)
    checks["three_groups_separate"] = set(result.get("groups", {})) == {"finite_j2", "P_correlated_negative", "j_ge_3"}
    checks["certificate_passed"] = result.get("caller_certificate_verified") is True
    checks["formal_basis"] = tuple(result.get("source_power_basis", ())) == ("hb", "F0_anchor", "Pstar_squared", "Lambda")
    checks["physical_sign"] = result.get("geometry", {}).get("physical_factor_negative") is True
    checks["hb_closed_coefficient_bound"] = result.get("geometry", {}).get("hb_coefficient_bound_only") is True
    checks["no_ledger"] = result.get("ledger_written") is False and result.get("no_ledger") is True
    checks["old_tail_retained"] = result.get("old_j2_tail", {}).get("retained") is True if isinstance(result.get("old_j2_tail"), Mapping) else False
    checks["old_tail_not_cancelled"] = result.get("old_j2_tail_cancelled") is False and result.get("old_j2_error_bound_cancelled") is False
    checks["q0_recorded"] = result.get("q0_exact_zero") is True
    checks["M_zero"] = _lo(c, result["M"]) == c.mpf(0) and _hi(c, result["M"]) == c.mpf(0)
    checks["Mz_nonzero"] = _hi(c, result["Mz"]) > c.mpf(0)
    checks["j3_no_M_division"] = result["geometry"].get("j3_no_M_division") is True

    tail_field = result["j_ge_3"]["field"]
    checks["chiZ_atom_retained_when_M_zero"] = any("chi-axis-Z1" in term.get("chi_Z_source_ids", []) for term in tail_field["Z1"])
    checks["j3_symmetric_error_only"] = all(term.get("error_only") is True and term.get("symmetric_error") is True and _sym(c, term["coefficient"]) for terms in tail_field.values() for term in terms)
    checks["j3_has_no_assigned_value"] = all(term.get("assigned_tail_value") is False and term.get("signed_tail_value") is False for terms in tail_field.values() for term in terms)
    checks["j3_power"] = all(term.get("powers", [None])[0] == 8 for terms in tail_field.values() for term in terms)

    finite_field = result["finite_j2"]["field"]
    checks["finite_Z0_power"] = finite_field["Z0"][0]["powers"] == [6, 1, 1, 0]
    routes = {term.get("route") for term in finite_field["Z1"]}
    checks["finite_Z1_three_paths"] = routes == {"chi2_gZ", "chiZ_outer_g0", "chiZ_inner_g0"}
    checks["chiZ_insertions_distinct"] = sum(term.get("route") == "chiZ_outer_g0" for term in finite_field["Z1"]) == 1 and sum(term.get("route") == "chiZ_inner_g0" for term in finite_field["Z1"]) == 1
    checks["finite_q0_metadata"] = all(term.get("q0_exact_zero") is True for terms in finite_field.values() for term in terms)
    checks["finite_signed_uncertainty_fields"] = all(term.get("error_only") is False for terms in finite_field.values() for term in terms)

    nonzero_chi, nonzero_chi_z, g5_nz, leading_nz, hb_nz, sc_nz, cert_nz, _ = _inputs(c, nonzero_chi=True)
    try:
        nonzero = build(c, hb_nz, sc_nz, nonzero_chi, nonzero_chi_z, g5_nz, leading_nz, cert_nz)
        checks["nontrivial_build"] = True
        signed = nonzero["finite_j2"]["Z0"][0]["coefficient"]
        checks["directed_negative_physical_j2"] = _lo(c, signed) < c.mpf(0) and _hi(c, signed) <= c.mpf(0)
        checks["physical_j2_power"] = nonzero["finite_j2"]["Z0"][0]["powers"] == [7, 1, 1, -2]
        checks["j3_physical_power"] = nonzero["j_ge_3"]["Z0"][0]["powers"] == [9, 1, 1, -2]
        checks["positive_green_quadrature"] = _quadrature_check(c, nonzero, sc_nz, hb_nz)
        pairs = nonzero["P_correlated_pairs"]
        checks["P_pairs_present"] = bool(pairs)
        checks["P_exact_directed_negative"] = all(_same(c, pair["P"]["coefficient"], _negate(c, pair["D2g"]["coefficient"])) for pair in pairs)
        checks["P_does_not_box_sum"] = all(pair["symbolic_sum"]["independent_boxed_interval_sum_performed"] is False for pair in pairs)
        checks["P_powers_match"] = all(pair["P"]["powers"] == pair["D2g"]["powers"] for pair in pairs)
    except Exception as exc:
        checks["nontrivial_build"] = False
        checks["nontrivial_error_absent"] = False
        return {"passed": False, "check_count": len(checks), "checks": checks, "error": repr(exc), "synthetic_only": True, "implementation_imported": False}

    return {"passed": all(checks.values()), "check_count": len(checks), "checks": checks, "synthetic_only": True, "implementation_imported": False, "ancestor_called": False, "quadrature_is_sanity_check": True}


run = run_fixture
