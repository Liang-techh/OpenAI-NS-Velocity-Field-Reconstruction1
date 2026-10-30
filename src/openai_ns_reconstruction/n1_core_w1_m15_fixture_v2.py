"""Corrected independent fixture for the core ``W1`` m15 adapter.

This fixture is intentionally self-contained.  It builds eight variable-radius
target cells, calls the supplied leaf, and checks every axial/swirl ordinary
``Z0``/``Z1`` cell against independent high-precision ``mpmath`` integrals.
The swirl reference includes the formal ``F0_anchor**2`` factor represented by
the two input ``F`` powers; the v1 reference omitted that factor while its
sector evaluator retained it, causing every swirl enclosure check to fail.
"""

from __future__ import annotations

from typing import Any


def _term(value: Any, source_id: str, powers=(0, 0, 0, 0), *, error: bool = False) -> dict[str, Any]:
    return {
        "powers": list(powers),
        "coefficient": value,
        "origin": {
            "source_ids": [source_id],
            "error_classes": ["synthetic_w1_error"] if error else [],
            "uncertainty": error,
        },
    }


def _endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _contains(got: Any, reference: Any) -> bool:
    raw = getattr(got, "_mpi_", None)
    if raw is None:
        return bool(got == reference)
    from mpmath import mp

    lower = mp.make_mpf(raw[0])
    upper = mp.make_mpf(raw[1])
    return bool(lower <= reference <= upper)


def _evaluate(c: Any, terms: Any) -> Any:
    bases = (c.mpf(1) / c.mpf(3), c.mpf(2), c.mpf(3), c.mpf(2))
    total = c.mpf(0)
    for sector in terms:
        value = sector["coefficient"]
        for base, power in zip(bases, sector["powers"]):
            value *= base**power
        total += value
    return total


def _ids(terms: Any) -> set[Any]:
    ids: set[Any] = set()
    for sector in terms:
        origin = sector.get("origin", {})
        values = origin.get("source_ids", ()) if isinstance(origin, dict) else ()
        ids.update(values)
    return ids


def run(leaf: Any, c: Any) -> dict[str, Any]:
    """Run the corrected variable-radius enclosure against *leaf*."""

    import mpmath as mp

    leading_u = {
        0: [_term(c.mpf([2, 3]), "lead_U_Z0")],
        1: [_term(c.mpf(3), "lead_U_Z1", error=True)],
    }
    leading_f = {
        0: [_term(c.mpf([1, 2]), "lead_F_Z0", (0, 1, 0, 0))],
        1: [_term(c.mpf(2), "lead_F_Z1", (0, 1, 0, 0), error=True)],
    }

    profile = []
    expected_ids: set[str] = set()
    for index in range(8):
        left = c.mpf(index) / c.mpf(8)
        right = c.mpf(index + 1) / c.mpf(8)
        t = c.make_mpf((_endpoint(left, False), _endpoint(right, True)))
        rho = c.mpf(4) * t * t
        fields = {
            0: {
                "1": [_term(c.mpf(2) - t, f"W1_F0:{index}", (0, 1, 0, 0))],
                "2": [_term(c.mpf(1) + t, f"W1_U0:{index}")],
            },
            1: {
                "1": [_term(c.mpf(1), f"W1_F1:{index}", (0, 1, 0, 0), error=True)],
                "2": [_term(c.mpf(2), f"W1_U1:{index}", error=True)],
            },
        }
        for row in fields.values():
            for sectors in row.values():
                expected_ids.update(sectors[0]["origin"]["source_ids"])
        profile.append({"t": t, "rho": rho, "W0_ordinary_Z": fields})

    source_atoms = {
        name: f"synthetic_{name}"
        for name in ("hb", "F0_anchor", "Pstar_squared", "Lambda")
    }
    certificate = {
        "geometry_epsilon_Z_independent": True,
        "actual_W1_target_cells": True,
    }
    output = leaf.build(
        c,
        profile=profile,
        leading_U=leading_u,
        leading_F=leading_f,
        source_family="synthetic_core_W1",
        source_atoms=source_atoms,
        certificate=certificate,
    )

    mp.mp.dps = 90
    epsilon = mp.mpf(1) / 2
    f0_anchor = mp.mpf(2)
    checks: dict[str, bool] = {}
    for order in (0, 1):
        for index in range(8):
            left = mp.mpf(index) / 8
            right = mp.mpf(index + 1) / 8

            if order == 0:
                axial_product = lambda t: (2 + t) * (1 + t)
                swirl_product = lambda t: (1 + t) * (2 - t)
            else:
                axial_product = lambda t: 3 * (1 + t) + (2 + t) * 2
                swirl_product = lambda t: 2 * (2 - t) + (1 + t) * 1

            axial_reference = mp.quad(
                lambda t: epsilon * 16 * t * axial_product(t),
                [left, right],
            )
            swirl_reference = mp.quad(
                lambda t: -(f0_anchor**2) * (epsilon**2) * 64 * t**3 * swirl_product(t),
                [left, right],
            )

            axial_terms = output["axial"][order][index]
            swirl_terms = output["swirl"][order][index]
            checks[f"axial_{order}_{index}_reference"] = _contains(
                _evaluate(c, axial_terms), axial_reference
            )
            checks[f"swirl_{order}_{index}_reference"] = _contains(
                _evaluate(c, swirl_terms), swirl_reference
            )
            checks[f"axial_{order}_{index}_ids"] = bool(_ids(axial_terms))
            checks[f"swirl_{order}_{index}_ids"] = bool(_ids(swirl_terms))
            checks[f"axial_{order}_{index}_four_powers"] = all(
                len(term["powers"]) == 4 for term in axial_terms
            )
            checks[f"swirl_{order}_{index}_four_powers"] = all(
                len(term["powers"]) == 4 for term in swirl_terms
            )

    all_output_ids = set()
    for channel in ("axial", "swirl"):
        for order in (0, 1):
            for index in range(8):
                all_output_ids.update(_ids(output[channel][order][index]))
    checks["source_ids_from_profile_retained"] = expected_ids <= all_output_ids
    checks["no_pressure_key"] = "pressure" not in output
    checks["g4_not_used"] = output["g4_used"] is False
    checks["ordinary_Z1_rule"] = output["ordinary_Z1_product_rule"] is True
    checks["not_complete_m15"] = output["complete_m15_computed"] is False
    checks["W1_included"] = output["W1_included"] is True
    return {
        "passed": all(checks.values()),
        "check_count": len(checks),
        "checks": checks,
        "synthetic_only": True,
        "reference": "independent high-precision variable-radius Z0/Z1 integrals with formal F0_anchor^2",
    }


__all__ = ["run"]
