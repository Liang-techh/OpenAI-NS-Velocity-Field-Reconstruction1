"""Independent exact-rational fixture for the local P/Dg pair leaf.

The fixture imports no implementation or producer.  It is inert until the
caller invokes ``run_fixture``; interval checks only inspect directed
endpoints, and the zero assertion is a correlated symbolic assertion rather
than a sum of two independently boxed intervals.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "n1_local_p_dg_correlated_pair_fixture_v1"
__all__ = ["SCHEMA", "run_fixture", "run"]


def _value(c: Any, value: Fraction | int) -> Any:
    if isinstance(value, Fraction):
        return c.mpf(value.numerator) / c.mpf(value.denominator)
    return c.mpf(value)


def _interval(c: Any, lower: Fraction | int, upper: Fraction | int | None = None) -> Any:
    if upper is None:
        upper = lower
    lo = _value(c, lower)
    hi = _value(c, upper)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((lo._mpi_[0], hi._mpi_[1]))
        except Exception:
            pass
    return c.mpf([lo, hi])


def _endpoint(value: Any, upper: bool) -> Any:
    if hasattr(value, "a") and hasattr(value, "b"):
        return value.b if upper else value.a
    raw = getattr(value, "_mpi_", None)
    return raw[1 if upper else 0] if raw is not None else value


def _negative_of(source: Any, candidate: Any) -> bool:
    try:
        return bool(
            _endpoint(candidate, False) == -_endpoint(source, True)
            and _endpoint(candidate, True) == -_endpoint(source, False)
        )
    except Exception:
        return False


def _term(
    c: Any,
    *,
    parent_id: str,
    lower: Fraction | int,
    upper: Fraction | int,
    powers: tuple[int, int, int, int],
    source_id: str,
    source_path: str,
    error_class: str,
    error_only: bool = False,
) -> dict[str, Any]:
    coefficient = _interval(c, lower, upper)
    return {
        "parent_id": parent_id,
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": list(powers),
        "formal_powers": list(powers),
        "source_ids": [source_id],
        "source_paths": [source_path],
        "error_classes": [error_class],
        "uncertainty": error_only,
        "error_only": error_only,
        "symmetric_error": error_only,
    }


def _packet(c: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """Create positive, negative, straddling, and error-only input sectors."""

    z0_positive = _term(
        c, parent_id="dg-positive-z0", lower=Fraction(1, 3), upper=Fraction(2, 3),
        powers=(1, -2, 3, 0), source_id="saved:Dg:J1:Z0:positive",
        source_path="saved/physical_Dg/J1/Z0/positive", error_class="saved_positive",
    )
    z0_negative = _term(
        c, parent_id="dg-negative-z0", lower=Fraction(-5, 2), upper=Fraction(-7, 6),
        powers=(2, 0, -1, 4), source_id="saved:Dg:J1:Z0:negative",
        source_path="saved/physical_Dg/J1/Z0/negative", error_class="saved_negative",
    )
    z1_straddling = _term(
        c, parent_id="dg-straddling-z1", lower=Fraction(-3, 2), upper=Fraction(5, 4),
        powers=(-1, 3, 0, -2), source_id="saved:Dg:J1:Z1:straddling",
        source_path="saved/physical_Dg/J1/Z1/straddling", error_class="saved_straddling",
    )
    z1_error = _term(
        c, parent_id="dg-error-z1", lower=Fraction(-1, 7), upper=Fraction(1, 7),
        powers=(0, 1, 2, -3), source_id="saved:Dg:J1:Z1:error",
        source_path="saved/physical_Dg/J1/Z1/error", error_class="saved_error",
        error_only=True,
    )
    packet = {
        "physical_level": "m15",
        "J1": {"Z0": [z0_positive, z0_negative], "Z1": [z1_straddling, z1_error]},
        "ordinary_derivatives": "ordinary Z0/Z1",
    }
    owner_by_parent = {
        term["parent_id"]: {"parent_id": term["parent_id"], "owner": "Dg_physical", "proved": True, "passed": True}
        for row in (packet["J1"]["Z0"], packet["J1"]["Z1"])
        for term in row
    }
    certificate = {
        "passed": True,
        "owner_by_parent": owner_by_parent,
        "same_g5": True,
        "same_K0": True,
        "same_chi": True,
        "P_sign": 2,
        "physical_measure_same": True,
        "Z_independent_ordinary_derivative": True,
        "certificate_scope": "saved physical J1 Dg sectors only",
    }
    return packet, certificate


def _pair_rows(result: Mapping[str, Any]) -> Mapping[str, Any]:
    rows = result.get("pairs", result.get("correlated_pairs"))
    if not isinstance(rows, Mapping):
        raise TypeError("result does not expose correlated pairs")
    if isinstance(rows.get("J1"), Mapping):
        rows = rows["J1"]
    if not isinstance(rows.get("Z0"), (list, tuple)) or not isinstance(rows.get("Z1"), (list, tuple)):
        raise TypeError("result pairs do not expose ordinary Z0/Z1 rows")
    return rows


def _call(build: Any, c: Any, packet: Any, certificate: Any) -> Mapping[str, Any]:
    function: Callable[..., Any] = getattr(build, "build", build)
    result = function(c, packet, certificate)
    if not isinstance(result, Mapping):
        raise TypeError("leaf build did not return a mapping")
    return result


def _independent_rational_identity() -> dict[str, bool]:
    """Use a fresh nontrivial rational kernel to illustrate Z0 and Z1."""

    def kernel(q: Fraction) -> Fraction:
        return q * (1 - q / 2)

    chi0 = Fraction(2, 5)
    chi1 = Fraction(-3, 7)
    g0 = Fraction(5, 4)
    g1 = Fraction(4, 9)
    leading = Fraction(7, 5)
    measure = Fraction(11, 13)
    cutoff = Fraction(5, 8)
    physical = leading * measure * cutoff
    q0, q1 = Fraction(1, 3), Fraction(2, 3)
    d0 = -2 * kernel(q0) * chi0 * g0
    p0 = +2 * kernel(q0) * chi0 * g0
    d1 = -2 * kernel(q1) * (chi1 * g0 + chi0 * g1)
    p1 = +2 * kernel(q1) * (chi1 * g0 + chi0 * g1)
    return {
        "kernel_nontrivial": kernel(q0) != 0 and kernel(q0) != kernel(q1),
        "Z0_functional_zero": d0 + p0 == 0,
        "Z1_functional_zero": d1 + p1 == 0,
        "Z0_physical_zero": physical * (d0 + p0) == 0,
        "Z1_physical_zero": physical * (d1 + p1) == 0,
        "physical_factor_nonzero": physical != 0,
    }


def run_fixture(build: Any, c: Any) -> dict[str, Any]:
    """Run exact endpoint, provenance, scope, and negative-guard checks."""

    packet, certificate = _packet(c)
    checks: dict[str, bool] = {}
    try:
        result = _call(build, c, packet, certificate)
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"build_completed": False},
            "check_count": 1,
            "error": f"{type(exc).__name__}: {exc}",
            "synthetic_only": True,
            "implementation_imported": False,
        }

    checks["build_completed"] = isinstance(result, Mapping)
    checks["schema"] = result.get("schema") == "n1_local_p_dg_correlated_pair_v1"
    checks["ordinary_Z_orders"] = result.get("ordinary_Z_orders") == [0, 1]
    checks["ordinary_Z_labels"] = result.get("ordinary_Z_labels") == ["Z0", "Z1"]
    validation = result.get("certificate_validation", {})
    checks["certificate_passed_is_verified"] = (
        validation.get("passed") is True and result.get("caller_certificate_verified") is True
    )
    checks["certificate_flags_are_passed"] = all(
        validation.get(name) is True
        for name in ("same_g5", "same_K0", "same_chi", "physical_measure_same", "Z_independent_ordinary_derivative")
    ) and validation.get("P_sign") == 2
    checks["one_pair_per_saved_sector"] = result.get("pair_count") == 4

    try:
        rows = _pair_rows(result)
        checks["pair_rows_shape"] = len(rows["Z0"]) == 2 and len(rows["Z1"]) == 2
    except Exception:
        rows = {"Z0": (), "Z1": ()}
        checks["pair_rows_shape"] = False

    source_rows = packet["J1"]
    for z_label in ("Z0", "Z1"):
        z_order = 0 if z_label == "Z0" else 1
        for index, pair in enumerate(rows[z_label]):
            parent = source_rows[z_label][index]
            dg = pair.get("Dg")
            p = pair.get("P")
            prefix = f"{z_label}_{index}"
            checks[f"{prefix}_parent_id"] = pair.get("parent_id") == parent["parent_id"]
            checks[f"{prefix}_ordinary_order"] = pair.get("ordinary_Z_order") == z_order
            checks[f"{prefix}_directed_negative"] = (
                isinstance(dg, Mapping) and isinstance(p, Mapping)
                and _negative_of(dg["coefficient"], p["coefficient"])
            )
            checks[f"{prefix}_four_formal_powers"] = (
                pair.get("powers") == parent["powers"]
                and pair.get("formal_powers") == parent["formal_powers"]
                and len(pair.get("powers", ())) == 4
                and p.get("powers") == parent["powers"]
            )
            checks[f"{prefix}_source_ids_preserved"] = (
                dg.get("source_ids") == parent["source_ids"]
                and p.get("source_ids") == parent["source_ids"]
                and pair.get("parent_source_ids") == parent["source_ids"]
            )
            checks[f"{prefix}_source_paths_preserved"] = (
                dg.get("source_paths") == parent["source_paths"]
                and p.get("source_paths") == parent["source_paths"]
                and pair.get("parent_source_paths") == parent["source_paths"]
            )
            checks[f"{prefix}_error_classes_preserved"] = (
                dg.get("error_classes") == parent["error_classes"]
                and p.get("error_classes") == parent["error_classes"]
                and pair.get("parent_error_classes") == parent["error_classes"]
            )
            symbolic = pair.get("symbolic_sum", {})
            checks[f"{prefix}_symbolic_zero_correlated_only"] = (
                symbolic.get("exact_zero") is True
                and symbolic.get("scope") == "correlated pair only"
                and symbolic.get("independent_boxed_interval_sum_performed") is False
                and pair.get("residual", {}).get("unconditioned_global") is False
            )
            checks[f"{prefix}_ordinary_scope"] = "ordinary Z" in str(pair.get("ordinary_derivative_scope"))
            checks[f"{prefix}_no_overlap_identity"] = pair.get("source_interval_overlap_used_as_identity") is False

    if rows["Z1"]:
        error_pair = rows["Z1"][1]
        checks["error_only_preserved"] = (
            error_pair["Dg"].get("error_only") is True
            and error_pair["P"].get("error_only") is True
            and error_pair["Dg"].get("uncertainty") is True
            and error_pair["P"].get("uncertainty") is True
            and error_pair["Dg"].get("symmetric_error") is True
            and error_pair["P"].get("symmetric_error") is True
        )
    else:
        checks["error_only_preserved"] = False
    positive_pair = rows["Z0"][0] if rows["Z0"] else {}
    negative_pair = rows["Z0"][1] if len(rows["Z0"]) > 1 else {}
    straddle_pair = rows["Z1"][0] if rows["Z1"] else {}
    try:
        checks["positive_interval_sign_reversed"] = _endpoint(positive_pair["P"]["coefficient"], True) < 0
        checks["negative_interval_sign_reversed"] = _endpoint(negative_pair["P"]["coefficient"], False) > 0
        checks["straddling_interval_remains_straddling"] = (
            _endpoint(straddle_pair["P"]["coefficient"], False) < 0
            and _endpoint(straddle_pair["P"]["coefficient"], True) > 0
        )
    except Exception:
        checks["positive_interval_sign_reversed"] = False
        checks["negative_interval_sign_reversed"] = False
        checks["straddling_interval_remains_straddling"] = False

    checks["no_ledger"] = result.get("ledger_written") is False and result.get("no_ledger") is True
    checks["no_full_P_branch_claim"] = (
        result.get("full_P_branch") is False
        and result.get("complete_P_branch") is False
        and result.get("core_W0_included") is False
        and result.get("D_P_included") is False
        and result.get("W0_tail_included") is False
    )
    checks["global_zero_not_claimed"] = (
        result.get("residual_scope") == "correlated pair only"
        and result.get("unconditioned_global_zero") is False
    )
    checks["implementation_mirror_not_used"] = result.get("implementation_mirror_used") is False
    checks.update({f"independent_{key}": value for key, value in _independent_rational_identity().items()})

    # Wrong and missing owner entries must fail closed.
    wrong_certificate = dict(certificate)
    wrong_owners = dict(certificate["owner_by_parent"])
    wrong_entry = dict(next(iter(wrong_owners.values())))
    wrong_entry["owner"] = "wrong-owner"
    wrong_owners[next(iter(wrong_owners))] = wrong_entry
    wrong_certificate["owner_by_parent"] = wrong_owners
    try:
        _call(build, c, packet, wrong_certificate)
    except Exception:
        checks["wrong_owner_guard"] = True
    else:
        checks["wrong_owner_guard"] = False

    missing_certificate = dict(certificate)
    missing_owners = dict(certificate["owner_by_parent"])
    missing_owners.pop(next(iter(missing_owners)))
    missing_certificate["owner_by_parent"] = missing_owners
    try:
        _call(build, c, packet, missing_certificate)
    except Exception:
        checks["missing_owner_guard"] = True
    else:
        checks["missing_owner_guard"] = False

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks_passed": all(checks.values()),
        "checks": checks,
        "check_count": len(checks),
        "synthetic_only": True,
        "implementation_imported": False,
        "old_code_called": False,
        "boxed_intervals_summed": False,
        "symbolic_zero_scope": "correlated pair only",
        "unconditioned_global_zero": False,
    }


run = run_fixture

