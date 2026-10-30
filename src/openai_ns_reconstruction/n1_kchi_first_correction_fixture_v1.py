"""Independent synthetic checks for :mod:`n1_kchi_first_correction_v1`.

The fixture supplies eight source cells whose coefficients enclose simple
polynomials in the source radius.  It checks the signed ``-2*chi`` factor,
the ordinary Z0/Z1/Z2 Leibniz rows, formal powers and source IDs, source
ordering, the diagonal triangular cell, the removable axis limits, and the
last-cell endpoint.  It is synthetic only: it does not read an archive or
claim an original n1 source identity.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


def _q(c: Any, value: Fraction) -> Any:
    return c.mpf(value.numerator) / value.denominator


def _cell(c: Any, index: int) -> Any:
    left = _q(c, Fraction(index, 4))
    right = _q(c, Fraction(index + 1, 4))
    return c.make_mpf((left._mpi_[0], right._mpi_[1]))


def _term(c: Any, coefficient: Any, source_id: str, *, powers=(0, 1, 0, -1)) -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": coefficient if hasattr(coefficient, "_mpi_") else c.mpf(coefficient),
        "origin": {"source_ids": [source_id], "error_classes": ["synthetic_source"]},
        "source_cell": source_id,
    }


def _interval(c: Any, lo: Fraction, hi: Fraction) -> Any:
    return c.make_mpf((_q(c, lo)._mpi_[0], _q(c, hi)._mpi_[1]))


def _source_cells(c: Any) -> list[dict[str, Any]]:
    cells: list[dict[str, Any]] = []
    for index in range(8):
        x = _cell(c, index)
        # h0=1+s, h1=2-s and h2=3+2s.  Each interval is a whole-cell
        # enclosure, so multiplying by a positive kernel mass is rigorous.
        lo, hi = Fraction(index, 4), Fraction(index + 1, 4)
        rows = {
            0: [_term(c, _interval(c, 1 + lo, 1 + hi), f"g5:{index}:Z0")],
            1: [_term(c, _interval(c, 2 - hi, 2 - lo), f"g5:{index}:Z1")],
            2: [_term(c, _interval(c, 3 + 2 * lo, 3 + 2 * hi), f"g5:{index}:Z2")],
        }
        cells.append({"x": x, "source_jets": rows})
    return cells


def _endpoint(value: Any) -> Fraction:
    sign, mantissa, exponent, _ = value
    result = Fraction(int(mantissa), 1)
    result = result * (1 << exponent) if exponent >= 0 else result / (1 << (-exponent))
    return -result if sign else result


def _contains(interval: Any, value: Fraction) -> bool:
    return _endpoint(interval._mpi_[0]) <= value <= _endpoint(interval._mpi_[1])


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run independent inclusion and bookkeeping checks against ``build``."""

    cells = _source_cells(c)
    chi = {
        0: [_term(c, c.mpf(1) / 2, "chi:Z0", powers=(0, 0, 0, 0))],
        1: [_term(c, c.mpf(3) / 4, "chi:Z1", powers=(0, 0, 0, 0))],
        2: [_term(c, -c.mpf(1) / 3, "chi:Z2", powers=(0, 0, 0, 0))],
    }
    result = build(c, cells, chi)
    targets = result.get("target_cells", [])
    checks: dict[str, bool] = {
        "eight_targets": len(targets) == 8,
        "source_j_le_target_i": all(
            target.get("source_indices") == list(range(index + 1))
            for index, target in enumerate(targets)
        ),
        "diagonal_triangular": all(
            target.get("triangular_diagonal_source") == index
            and target["W0_correction_by_source_cell"][str(index)]["triangular"]
            for index, target in enumerate(targets)
        ),
        "axis_continuation": bool(targets)
        and targets[0].get("axis_endpoint_continuation") is True,
        "formal_powers_retained": True,
        "negative_signed_correction": True,
    }

    # The axis triangular bare masses are x^4/192 and x^3/48.  The first
    # source interval is positive and chi_Z0 is positive, hence both signed
    # correction coefficients are nonpositive.
    axis = targets[0]["W0_correction_by_source_cell"]["0"]
    axis_f = axis["kernel_masses"]["F1"]
    axis_dx = axis["kernel_masses"]["d_x_F1"]
    checks["axis_mass_contains_upper_limit"] = _contains(axis_f, Fraction(1, 192 * 4**4))
    checks["axis_dx_mass_contains_upper_limit"] = _contains(axis_dx, Fraction(1, 48 * 4**3))

    for target in targets:
        for source_index in target["source_indices"]:
            by_source = target["W0_correction_by_source_cell"][str(source_index)]
            for field in ("F1", "d_x_F1"):
                for order in ("Z0", "Z1", "Z2"):
                    terms = by_source["correction"][field][order]
                    key = f"{field}_{order}_{source_index}"
                    checks[key + "_directed"] = bool(terms) and all(
                        hasattr(term["coefficient"], "_mpi_") for term in terms
                    )
                    checks[key + "_formal"] = bool(terms) and all(
                        tuple(term["powers"]) == (0, 1, 0, -1) for term in terms
                    )
                    checks[key + "_source_ids"] = bool(terms) and all(
                        term["source_cell"] == source_index
                        and term["source_index"] == source_index
                        and term["origin"]["source_ids"] for term in terms
                    )
                    if order == "Z2":
                        checks[key + "_middle_leibniz"] = any(
                            term["leibniz_multiplier"] == 2
                            and term["chi_order"] == 1
                            and term["source_order"] == 1
                            for term in terms
                        )
                    # The Z0 source and chi are positive, so the signed
                    # first correction is nonpositive.  Higher jets are
                    # intentionally mixed-sign in this fixture.
                    if order == "Z0":
                        checks[key + "_sign"] = all(term["coefficient"].b <= 0 for term in terms)

    # The exact endpoint x=2 is represented by the last target interval's
    # upper bound; a source cell is never taken from a future j>i.
    checks["last_target_endpoint"] = targets[-1]["x"].b == 2
    checks["no_future_source"] = all(
        max(target["source_indices"], default=-1) == index
        for index, target in enumerate(targets)
    )
    checks["no_higher_kchi"] = result.get("Kchi_higher_words_included") is False
    checks["no_allword_tail"] = result.get("allword_tail_added") is False
    checks["source_to_w0_only"] = result.get("source_to_W0_only") is True

    # Independent continuous integrals for nonconstant h(s), including the
    # moving diagonal endpoint.  The generated intervals enclose an entire
    # target cell, so x=1.9 need not be represented by a midpoint.
    import mpmath as mp
    mp.mp.dps = 70
    x = mp.mpf('1.9')
    for source_index, left, right in ((2, mp.mpf('0.5'), mp.mpf('0.75')),
                                      (7, mp.mpf('1.75'), x)):
        for field in ("F1", "d_x_F1"):
            def kernel(s):
                if field == "F1":
                    return s * (x*x - s**4/(x*x) - 4*s*s*mp.log(x/s)) / 16
                return s*x*(1-(s/x)**2)**2/8
            for order in range(3):
                def integrand(s):
                    g0=1+s; g1=2-s; g2=3+2*s
                    product=(g0/2, 3*g0/4+g1/2, -g0/3+3*g1/2+g2/2)[order]
                    return -2*kernel(s)*product
                reference=mp.quad(integrand,[left,right])
                terms=targets[7]["W0_correction_by_source_cell"][str(source_index)]["correction"][field]["Z"+str(order)]
                bound=sum((term["coefficient"] for term in terms),c.mpf(0))
                low=mp.make_mpf(bound._mpi_[0]); high=mp.make_mpf(bound._mpi_[1])
                checks[f"variable_source_{source_index}_{field}_Z{order}"]=bool(low<=reference<=high)

    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "fixture_is_not_original_source": True,
        "method": "independent polynomial source, exact axis limits, and directed cell-mass bookkeeping",
    }


run = run_fixture


__all__ = ["run_fixture", "run"]

