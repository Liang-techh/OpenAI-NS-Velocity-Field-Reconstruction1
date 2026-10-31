"""Independent synthetic checks for the local outer-D physical m15 leaf.

The fixture owns all source sectors and performs an independent high-
precision q integral.  Importing this file is inert; the root runner is the
only caller of :func:`run_fixture`.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")


def _q(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _interval(c: Any, left: Any, right: Any) -> Any:
    lo = left if hasattr(left, "_mpi_") else c.mpf(left)
    hi = right if hasattr(right, "_mpi_") else c.mpf(right)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _term(
    c: Any,
    value: Any,
    source_id: str,
    powers: tuple[int, int, int, int],
    *,
    path: tuple[str, ...],
    error_only: bool = False,
    error_class: str = "synthetic_source",
) -> dict[str, Any]:
    coefficient = value if hasattr(value, "_mpi_") else c.mpf(value)
    result = {
        "coefficient": coefficient,
        "powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [path],
        "error_classes": [error_class],
        "uncertainty": error_only,
        "origin": {
            "source_ids": [source_id],
            "source_paths": [path],
            "error_classes": [error_class],
            "uncertainty": error_only,
            "error_only": error_only,
        },
    }
    if error_only:
        result.update({"error_only": True, "symmetric_error": True})
    return result


def _point_or_interval(c: Any, left: Fraction | int, right: Fraction | int | None = None) -> Any:
    if right is None:
        right = left
    return _interval(c, _q(c, left), _q(c, right))


def _evaluate_terms(c: Any, terms: list[dict[str, Any]], bases: tuple[Any, ...]) -> Any:
    total = c.mpf(0)
    for term in terms:
        value = term["coefficient"]
        for base, power in zip(bases, term["powers"]):
            value *= base ** power
        total += value
    return total


def _contains(interval: Any, value: Any) -> bool:
    from mpmath import mp

    target = mp.mpf(value)
    return mp.make_mpf(interval._mpi_[0]) <= target <= mp.make_mpf(interval._mpi_[1])


def _source_packet(c: Any) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    w0 = _term(
        c,
        _point_or_interval(c, Fraction(109, 100), Fraction(111, 100)),
        "synthetic:W1:F1:Z0",
        (0, 0, 0, 0),
        path=("synthetic", "low_word", "0", "F1", "saved"),
        error_class="synthetic_W1",
    )
    w1 = _term(
        c,
        _point_or_interval(c, Fraction(79, 100), Fraction(81, 100)),
        "synthetic:W1:F1:Z1",
        (0, 0, 0, 0),
        path=("synthetic", "low_word", "1", "F1", "saved"),
        error_class="synthetic_W1_Z1",
    )
    w1_error = _term(
        c,
        _interval(c, _q(c, Fraction(-1, 100)), _q(c, Fraction(1, 100))),
        "synthetic:W1:F1:Z1:error",
        (0, 0, 0, 0),
        path=("synthetic", "low_word", "1", "F1", "error"),
        error_only=True,
        error_class="synthetic_W1_error",
    )
    low_word = {
        "local_signed_sectors": {
            "0": {"F1": [w0]},
            "1": {"F1": [w1, w1_error]},
        },
        "whole_collar": True,
    }

    chi0 = _term(
        c,
        _point_or_interval(c, Fraction(69, 100), Fraction(71, 100)),
        "synthetic:axis_chi:Z0",
        (0, 0, 1, 0),
        path=("synthetic", "axis_chi", "Z0"),
        error_class="synthetic_chi",
    )
    chiz = _term(
        c,
        _point_or_interval(c, Fraction(-41, 100), Fraction(-39, 100)),
        "synthetic:axis_chi:Z1",
        (0, 0, 1, 0),
        path=("synthetic", "axis_chi", "Z1"),
        error_class="synthetic_chi_Z",
    )
    axis_chi = {"Z0": [chi0], "Z1": [chiz]}

    leading0 = _term(
        c,
        _point_or_interval(c, Fraction(119, 100), Fraction(121, 100)),
        "synthetic:leading_F:Z0",
        (0, 1, 0, 0),
        path=("synthetic", "leading_F", "Z0"),
        error_class="synthetic_leading",
    )
    leading1 = _term(
        c,
        _point_or_interval(c, Fraction(-61, 100), Fraction(-59, 100)),
        "synthetic:leading_F:Z1",
        (0, 1, 0, 0),
        path=("synthetic", "leading_F", "Z1"),
        error_class="synthetic_leading_Z",
    )
    leading = {"F1": {"Z0": [leading0], "Z1": [leading1]}}
    return low_word, axis_chi, leading


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run independent interval, provenance, and numerical checks."""

    import mpmath as mp

    sc = _q(c, Fraction(1, 4096))
    low_word, axis_chi, leading = _source_packet(c)
    atoms = {name: "synthetic:" + name for name in SOURCE_POWER_BASIS}
    certificate = {
        "hb": _interval(c, _q(c, 0), _q(c, Fraction(1, 2))),
        "strict_hb_admitted": True,
        "geometry_Z_independent": True,
        "q_all_phase": True,
        "original_cutoff_mass_5_over_8": True,
        "axis_chi_is_not_phase_cutoff": True,
        "fixture_owner": "independent_local_outer_dq_v1",
    }
    result = leaf.build(
        c,
        sc=sc,
        low_word=low_word,
        axis_chi=axis_chi,
        leading_F=leading,
        source_family="synthetic:raw_Q:C0_B1",
        source_atoms=atoms,
        certificate=certificate,
    )

    bases = (c.mpf(1) / 3, c.mpf(2), c.mpf(3), c.mpf(5))
    d0 = _evaluate_terms(c, result["D"]["F1"]["Z0"], bases)
    d1 = _evaluate_terms(c, result["D"]["F1"]["Z1"], bases)
    p0 = _evaluate_terms(c, result["swirl"]["Z0"], bases)
    p1 = _evaluate_terms(c, result["swirl"]["Z1"], bases)

    mp.mp.dps = 90
    hb_mp = mp.mpf(1) / 3
    sc_mp = mp.mpf(1) / 4096
    f0_mp, pstar_mp, lam_mp = mp.mpf(2), mp.mpf(3), mp.mpf(5)
    lead0_mp = mp.mpf(6) / 5 * f0_mp
    lead1_mp = -mp.mpf(3) / 5 * f0_mp
    chi0_mp = mp.mpf(7) / 10 * pstar_mp
    chiz_mp = -mp.mpf(2) / 5 * pstar_mp
    w0_mp = mp.mpf(11) / 10
    w1_mp = mp.mpf(4) / 5

    def kernel(q: Any) -> Any:
        d = hb_mp * sc_mp * q / 4
        return mp.cosh(d) - 1

    def eta(q: Any) -> Any:
        return mp.mpf(1) / 2 + q / 4

    def radial(q: Any) -> Any:
        return mp.exp(hb_mp * sc_mp * q / 2)

    def d0_value(q: Any) -> Any:
        return -2 * chi0_mp * w0_mp * kernel(q)

    def d1_value(q: Any) -> Any:
        return -2 * (chiz_mp * w0_mp + chi0_mp * w1_mp) * kernel(q)

    measure = -8 * hb_mp * sc_mp / lam_mp**2
    reference_d0 = mp.quad(lambda q: d0_value(q), [0, 1])
    reference_d1 = mp.quad(lambda q: d1_value(q), [0, mp.mpf("0.5"), 1])
    reference_p0 = mp.quad(lambda q: measure * eta(q) * radial(q) * lead0_mp * d0_value(q), [0, mp.mpf("0.5"), 1])
    reference_p1 = mp.quad(
        lambda q: measure * eta(q) * radial(q) * (lead1_mp * d0_value(q) + lead0_mp * d1_value(q)),
        [0, mp.mpf("0.5"), 1],
    )

    checks: dict[str, bool] = {}
    checks["schema"] = result.get("schema") == "n1_local_outer_dq_m15_v1"
    checks["only_F"] = result.get("fields") == ["F1"] and result.get("only_F_output") is True
    checks["raw_Q_partial"] = result.get("raw_Q_partial") is True and result.get("raw_Q_components") == ["C0", "B1"]
    checks["kernel_mass_starts_zero"] = result["kernel"]["coefficient_mass"].a == 0
    checks["kernel_upper_is_bound"] = result["kernel"]["upper_is_bound_not_value"] is True and result["kernel"]["coefficient_mass"].b > 0
    checks["q_zero_is_retained"] = all(
        term["coefficient"].a <= 0 <= term["coefficient"].b
        for rows in (result["D"]["F1"], result["swirl"])
        for terms in rows.values()
        for term in terms
    )
    checks["D_signed_negative"] = d0.b <= 0
    checks["physical_signed_positive_for_Z0"] = p0.a >= 0
    checks["D_direct_integral"] = _contains(d0, reference_d0)
    checks["D_Z1_direct_integral"] = _contains(d1, reference_d1)
    checks["physical_radial_integral"] = _contains(p0, reference_p0) and _contains(p1, reference_p1)
    checks["nontrivial_kernel"] = bool(kernel(mp.mpf(1)) > kernel(mp.mpf("0.5")) > 0)
    checks["cutoff_mass_certificate"] = result["geometry"]["cutoff_mass"] == c.mpf(5) / c.mpf(8)
    checks["closed_hb_proved"] = result["geometry"]["closed_hb_enclosure_allowed"] is True
    checks["formal_total_powers"] = all(
        term["powers"] == [3, 1, 1, -2]
        for rows in result["swirl"].values()
        for term in rows
    )
    checks["formal_D_hb_power"] = all(
        term["powers"][0] == 2
        for rows in result["D"]["F1"].values()
        for term in rows
    )
    checks["three_Z1_paths"] = {term["leibniz_piece"] for term in result["swirl"]["Z1"]} >= {
        "leading_Z*chi*W1_0", "leading*chi_Z*W1_0", "leading*chi*W1_1"
    }
    checks["chi_Z_and_leading_Z_retained"] = any(
        term["leibniz_piece"] == "leading*chi_Z*W1_0" and "synthetic:axis_chi:Z1" in term["source_ids"]
        for term in result["swirl"]["Z1"]
    ) and any(term["leibniz_piece"] == "leading_Z*chi*W1_0" for term in result["swirl"]["Z1"])
    checks["source_paths_retained"] = all(
        term["source_paths"] for rows in result["swirl"].values() for term in rows
    )
    checks["error_only_symmetric"] = all(
        term["coefficient"].a == -term["coefficient"].b
        for rows in (result["D"]["F1"], result["swirl"])
        for terms in rows.values()
        for term in terms
        if term.get("error_only")
    )
    checks["no_P_or_pressure"] = (
        result["P_reapplied"] is False
        and result["pressure_included"] is False
        and result["P1_included"] is False
        and result["physical_P_compensation_included"] is False
    )
    checks["no_reentry"] = result["green_reapplied"] is False and result["W1_reapplied"] is False
    checks["no_ledger_or_cancellation"] = result["ledger_written"] is False and result["same_source_cancellation_claimed"] is False
    checks["axis_chi_not_phase_cutoff"] = result["axis_chi_is_not_phase_cutoff"] is True

    weak = dict(certificate)
    weak.pop("strict_hb_admitted")
    try:
        leaf.build(
            c,
            sc=sc,
            low_word=low_word,
            axis_chi=axis_chi,
            leading_F=leading,
            source_family="synthetic:raw_Q:C0_B1",
            source_atoms=atoms,
            certificate=weak,
        )
    except (TypeError, ValueError):
        checks["closed_hb_requires_proof"] = True
    else:
        checks["closed_hb_requires_proof"] = False

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "independent_kernel_and_radial_integrals": True,
        "upper_mass_is_not_value": True,
        "receipt_reads": False,
        "old_builder_called": False,
        "automatic_execution": False,
    }


run = run_fixture

__all__ = ["SOURCE_POWER_BASIS", "run_fixture", "run"]
