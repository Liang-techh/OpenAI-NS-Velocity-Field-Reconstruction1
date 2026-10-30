"""Independent high-precision checks for :mod:`n1_kchi_feedback_m15_v2`.

The fixture is synthetic only.  It supplies nonconstant affine source fields
on each target cell, evaluates the same physical integrals independently with
ordinary ``mpmath`` quadrature, and checks enclosure, signs, formal powers,
ordinary ``Z1``, the three branch split, and provenance.  It intentionally
does not read a receipt, invoke a producer, or turn a tail bound into a
signed field.
"""

from __future__ import annotations

from typing import Any


BRANCHES = ("W0_chi", "E01", "E10", "E11")
W1_BRANCHES = ("E01", "E10", "E11")


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    lo = c.mpf(lower)
    hi = c.mpf(upper)
    return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))


def _cell_t(c: Any, index: int) -> Any:
    return _interval(c, c.mpf(index) / 8, c.mpf(index + 1) / 8)


def _term(c: Any, value: Any, source: str, *, error_class: str | None = None, uncertainty: bool = False, **markers: Any) -> dict[str, Any]:
    origin: dict[str, Any] = {"source_ids": [source], "uncertainty": uncertainty}
    if error_class is not None:
        origin["error_classes"] = [error_class]
    else:
        origin["error_classes"] = []
    origin.update(markers)
    return {
        "powers": (0, 0, 0, 0),
        "coefficient": value,
        "origin": origin,
    }


def _affine_bounds(c: Any, function: Any, left: Any, right: Any) -> Any:
    from mpmath import mp

    values = (function(left), function(right))
    lower_index = min(
        range(2), key=lambda index: mp.make_mpf(values[index]._mpi_[0])
    )
    upper_index = max(
        range(2), key=lambda index: mp.make_mpf(values[index]._mpi_[1])
    )
    return c.make_mpf((values[lower_index]._mpi_[0], values[upper_index]._mpi_[1]))


def _inputs(c: Any) -> tuple[Any, Any, Any]:
    """Return profile, leading U, leading F without using the leaf."""

    profile: list[dict[str, Any]] = []
    for cell in range(8):
        t = _cell_t(c, cell)
        left = c.mpf(cell) / 8
        right = c.mpf(cell + 1) / 8
        # The three signed routes have different affine, nonconstant sources.
        branch_data: dict[str, Any] = {}
        branch_specs = {
            "W0_chi": (
                lambda q: c.mpf("0.65") + c.mpf("0.12") * q,
                lambda q: c.mpf("0.35") + c.mpf("0.08") * q,
                lambda q: c.mpf("0.06") + c.mpf("0.01") * q,
                lambda q: c.mpf("0.04") + c.mpf("0.02") * q,
                "chi",
            ),
            "E01": (
                lambda q: -c.mpf("0.30") - c.mpf("0.04") * q,
                lambda q: -c.mpf("0.18") - c.mpf("0.02") * q,
                lambda q: c.mpf("-0.03") + c.mpf("0.01") * q,
                lambda q: c.mpf("-0.02") + c.mpf("0.01") * q,
                "B0",
            ),
            "E10": (
                lambda q: c.mpf("0.22") + c.mpf("0.03") * q,
                lambda q: c.mpf("0.11") + c.mpf("0.02") * q,
                lambda q: c.mpf("0.02") + c.mpf("0.01") * q,
                lambda q: c.mpf("0.01") + c.mpf("0.01") * q,
                "B1",
            ),
            "E11": (
                lambda q: -c.mpf("0.16") - c.mpf("0.03") * q,
                lambda q: -c.mpf("0.09") - c.mpf("0.01") * q,
                lambda q: c.mpf("0.015") + c.mpf("0.006") * q,
                lambda q: c.mpf("0.008") + c.mpf("0.004") * q,
                "E11",
            ),
        }
        for branch, (f1, f1_z, uz1, uz1_z, tag) in branch_specs.items():
            f1_term = _term(
                c,
                _affine_bounds(c, f1, left, right),
                f"{tag}-F1-{cell}",
            )
            f1_z_term = _term(
                c,
                _affine_bounds(c, f1_z, left, right),
                f"{tag}-F1Z-{cell}",
            )
            uz1_term = _term(
                c,
                _affine_bounds(c, uz1, left, right),
                f"{tag}-Uz1-{cell}",
            )
            uz1_z_term = _term(
                c,
                _affine_bounds(c, uz1_z, left, right),
                f"{tag}-Uz1Z-{cell}",
            )
            payload: dict[str, Any] = {
                "Z0": {"F1": [f1_term]},
                "Z1": {"F1": [f1_z_term]},
            }
            # Deliberately omit E10 Uz1: absent fields are natural zero.
            if branch != "E10":
                payload["Z0"]["Uz1"] = [uz1_term]
                payload["Z1"]["Uz1"] = [uz1_z_term]
            if branch == "E01":
                payload["Z0"]["F1"] = {
                    "field": [f1_term],
                    "tail_bound": [
                        _term(
                            c,
                            _interval(c, 0, c.mpf("0.005")),
                            f"{tag}-tail-{cell}",
                            error_class="kchi_tail",
                            uncertainty=True,
                            tail_bound=True,
                            bound_type="tail_bound",
                        )
                    ],
                }
            if branch == "E10":
                payload["Z1"]["F1"] = {
                    "field": [f1_z_term],
                    "symmetric_error": [
                        _term(
                            c,
                            _interval(c, 0, c.mpf("0.004")),
                            f"{tag}-sym-{cell}",
                            error_class="symmetric_kchi",
                            uncertainty=True,
                        )
                    ],
                }
            branch_data[branch] = payload
        profile.append({"t": t, "rho": c.mpf(4) * t * t, "branches": branch_data})

    def leading_u(order: int) -> list[dict[str, Any]]:
        if order == 0:
            function = lambda q: c.mpf("1.10") + c.mpf("0.20") * q
            source = "Ulead0"
        else:
            function = lambda q: c.mpf("0.25") + c.mpf("0.03") * q
            source = "Ulead1"
        return [
            _term(
                c,
                _interval(c, function(c.mpf(0)), function(c.mpf(1))),
                source,
            )
        ]

    def leading_f(order: int) -> list[dict[str, Any]]:
        if order == 0:
            function = lambda q: c.mpf("1.35") + c.mpf("0.10") * q
            source = "Flead0"
        else:
            function = lambda q: c.mpf("-0.20") + c.mpf("0.02") * q
            source = "Flead1"
        return [
            _term(
                c,
                _interval(c, function(c.mpf(0)), function(c.mpf(1))),
                source,
            )
        ]

    return profile, {0: leading_u(0), 1: leading_u(1)}, {0: leading_f(0), 1: leading_f(1)}


def _mp_bounds(value: Any) -> tuple[Any, Any]:
    from mpmath import mp

    return mp.make_mpf(value._mpi_[0]), mp.make_mpf(value._mpi_[1])


def _sector_value(c: Any, terms: Any, bases: tuple[Any, Any, Any, Any]) -> Any:
    total = c.mpf(0)
    for term in terms:
        value = term["coefficient"]
        for base, power in zip(bases, term["powers"]):
            value *= base**power
        total += value
    return total


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Run v2 synthetic checks; callers choose whether to execute the fixture."""

    from mpmath import mp

    mp.dps = 90
    profile, leading_u, leading_f = _inputs(c)
    atoms = {name: f"synthetic_{name}" for name in ("hb", "F0_anchor", "Pstar_squared", "Lambda")}
    result = leaf.build(
        c,
        profile=profile,
        leading_U=leading_u,
        leading_F=leading_f,
        source_family="synthetic_kchi_feedback",
        source_atoms=atoms,
        certificate={"geometry_epsilon_Z_independent": True},
    )

    checks: dict[str, bool] = {}
    checks["three_branches"] = tuple(result.get("branch_order", ())) == BRANCHES
    checks["three_w1_routes"] = all(route in result.get("branch_order", ()) for route in W1_BRANCHES)
    checks["eight_targets"] = len(result.get("target_cells", ())) == 8
    checks["g4_not_used"] = result.get("g4_used") is False and result.get("pressure_included") is False
    checks["tail_is_separate"] = result.get("tail_bounds_as_signed_fields") is False
    checks["not_complete_m15"] = result.get("complete_m15_computed") is False

    # Independent variable-radius references.  The leading source is the
    # affine function, while each branch is integrated separately on each
    # target cell.  The output sectors are evaluated at Lambda=8.
    eps = mp.mpf(1) / 8
    bases = (mp.mpf(1), mp.mpf(1), mp.mpf(1), mp.mpf(8))

    lead_u = {
        0: lambda q: mp.mpf("1.10") + mp.mpf("0.20") * q,
        1: lambda q: mp.mpf("0.25") + mp.mpf("0.03") * q,
    }
    lead_f = {
        0: lambda q: mp.mpf("1.35") + mp.mpf("0.10") * q,
        1: lambda q: mp.mpf("-0.20") + mp.mpf("0.02") * q,
    }
    source_specs = {
        "W0_chi": (
            lambda q: mp.mpf("0.65") + mp.mpf("0.12") * q,
            lambda q: mp.mpf("0.35") + mp.mpf("0.08") * q,
            lambda q: mp.mpf("0.06") + mp.mpf("0.01") * q,
            lambda q: mp.mpf("0.04") + mp.mpf("0.02") * q,
        ),
        "E01": (
            lambda q: -mp.mpf("0.30") - mp.mpf("0.04") * q,
            lambda q: -mp.mpf("0.18") - mp.mpf("0.02") * q,
            lambda q: mp.mpf("-0.03") + mp.mpf("0.01") * q,
            lambda q: mp.mpf("-0.02") + mp.mpf("0.01") * q,
        ),
        "E10": (
            lambda q: mp.mpf("0.22") + mp.mpf("0.03") * q,
            lambda q: mp.mpf("0.11") + mp.mpf("0.02") * q,
            lambda q: mp.mpf("0.02") + mp.mpf("0.01") * q,
            lambda q: mp.mpf("0.01") + mp.mpf("0.01") * q,
        ),
        "E11": (
            lambda q: -mp.mpf("0.16") - mp.mpf("0.03") * q,
            lambda q: -mp.mpf("0.09") - mp.mpf("0.01") * q,
            lambda q: mp.mpf("0.015") + mp.mpf("0.006") * q,
            lambda q: mp.mpf("0.008") + mp.mpf("0.004") * q,
        ),
    }

    def reference(branch: str, order: int, channel: str) -> Any:
        f1, f1_z, uz1, uz1_z = source_specs[branch]
        total = mp.mpf(0)
        for cell in range(8):
            left = mp.mpf(cell) / 8
            right = mp.mpf(cell + 1) / 8

            def integrand(q: Any) -> Any:
                rho = 4 * q * q
                if channel == "axial":
                    word = uz1(q) if order == 0 else uz1_z(q)
                    if branch == "E10":
                        word = mp.mpf(0) if order == 0 else mp.mpf(0)
                    if order == 1:
                        word = lead_u[1](q) * uz1(q) + lead_u[0](q) * uz1_z(q)
                        if branch == "E10":
                            word = mp.mpf(0)
                    return 2 * eps * lead_u[0 if order == 0 else 0](q) * word * 8 * q if order == 0 else 2 * eps * word * 8 * q
                word = f1(q) if order == 0 else f1_z(q)
                if order == 1:
                    word = lead_f[1](q) * f1(q) + lead_f[0](q) * f1_z(q)
                return -2 * eps**2 * rho * (lead_f[0](q) * word if order == 0 else word) * 8 * q

            total += mp.quad(integrand, [left, (left + right) / 2, right])
        return total

    for branch in BRANCHES:
        for order in (0, 1):
            axial_terms = result["branches"][branch][order]["axial"]
            swirl_terms = result["branches"][branch][order]["swirl"]
            axial_value = _sector_value(c, axial_terms, tuple(c.mpf(base) for base in bases))
            swirl_value = _sector_value(c, swirl_terms, tuple(c.mpf(base) for base in bases))
            axial_ref = reference(branch, order, "axial")
            swirl_ref = reference(branch, order, "swirl")
            axial_lo, axial_hi = _mp_bounds(axial_value)
            swirl_lo, swirl_hi = _mp_bounds(swirl_value)
            checks[f"axial_encloses_{branch}_{order}"] = axial_lo <= axial_ref <= axial_hi
            checks[f"swirl_encloses_{branch}_{order}"] = swirl_lo <= swirl_ref <= swirl_hi
            checks[f"axial_units_{branch}_{order}"] = all(term["powers"][3] == -1 for term in axial_terms)
            checks[f"swirl_units_{branch}_{order}"] = all(term["powers"][3] == -2 for term in swirl_terms)
            checks[f"provenance_{branch}_{order}"] = all(
                term["origin"].get("source_ids") for term in axial_terms + swirl_terms
            )

    checks["ordinary_z1_product_rule"] = all(
        any(term["origin"].get("derivative_piece") == "leading1*word0" for term in result["branches"][branch][1]["swirl"])
        and any(term["origin"].get("derivative_piece") == "leading0*word1" for term in result["branches"][branch][1]["swirl"])
        for branch in BRANCHES
    )
    checks["missing_field_is_zero"] = not result["branches"]["E10"][0]["axial"]
    checks["symmetric_error_not_signed"] = bool(
        result["branches"]["E10"][1]["swirl_error"]
    ) and all(
        term["origin"].get("signed_field") is False
        for term in result["branches"]["E10"][1]["swirl_error"]
    )
    checks["tail_not_signed"] = all(
        term["origin"].get("bound_type") != "tail_bound"
        for term in result["branches"]["E01"][0]["swirl"]
    )
    checks["signs"] = (
        _mp_bounds(_sector_value(c, result["branches"]["W0_chi"][0]["swirl"], tuple(c.mpf(base) for base in bases)))[1] < 0
        and _mp_bounds(_sector_value(c, result["branches"]["E01"][0]["swirl"], tuple(c.mpf(base) for base in bases)))[0] > 0
    )
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "reference": "independent high-precision quadrature of nonconstant affine sources in rho=4*t^2",
        "fixture_is_not_actual_field": True,
    }


run_fixture = run

__all__ = ["BRANCHES", "run", "run_fixture"]
