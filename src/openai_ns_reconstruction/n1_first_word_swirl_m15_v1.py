"""Finite first-word physical swirl ``m15`` source sectors.

This module is a source-only adapter for the first collar.  It consumes the
caller-supplied whole-phase leading physical ``Flead`` sectors and the raw
physical first-word ``F1_W1`` sectors.  No producer is evaluated and no
archive or ledger is read.  The omitted all-word tail is deliberately not
added.

With ``s = sc*q/4`` and ``rho = 4*exp(hb*s)``, the certified cutoff identity
``int_0^1 eta(q)dq = 5/8`` (``eta(q) = 1-sigma(4*q-2)``) turns

    -2*epsilon**2*hb*(sc/4)*int rho**2*Flead*F1_W1*eta dq

into the signed outer factor ``-5*sc`` times the positive interval enclosure
for ``exp(2*hb*s)``.  Since ``0 < hb < 1/2`` and ``sc < 1/1024``, that factor
is enclosed by ``[1, 1/(1-sc/4)]``.  The explicit formal powers are one
``hb`` and two inverse powers of ``Lambda``.

Core and local branches remain independent.  Products are outward-added only
when formal powers, uncertainty, and error-class sets agree.  Source IDs and
error classes are retained in every aggregate; addition is not a cancellation
claim.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
BRANCHES = ("core", "local")

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "BRANCHES",
    "build",
]


def _is_interval(value: Any) -> bool:
    """Return whether *value* is a caller-directed interval."""

    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Build an interval from directed endpoints without scalar promotion."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    """Expose interval endpoints as exact scalar MPF values for guards."""

    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    from mpmath import mp

    return mp.make_mpf(raw[0]), mp.make_mpf(raw[1])


def _as_items(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if isinstance(value, Sequence):
        return list(value)
    return [value]


def _freeze(value: Any) -> Any:
    try:
        hash(value)
    except TypeError:
        return (type(value).__name__, repr(value))
    return (type(value).__name__, value)


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        key = _freeze(value)
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Read compact provenance while retaining source/error unions."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {
            "source_ids": [],
            "error_classes": [],
            "uncertainty": False,
            "path_count": 1,
            "source_count": 1,
        }
    if isinstance(value, Mapping):
        _seen.add(marker)
        ids = list(_as_items(value.get("source_ids", value.get("sources"))))
        classes = list(_as_items(value.get(
            "error_classes",
            value.get("uncertainty_class", value.get("error_class")),
        )))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(
            token in operation or token in error_kind
            for token in ("cauchy", "symmetric", "uncertainty")
        )
        nested: list[dict[str, Any]] = []
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested.extend(_origin_info(item, _seen=_seen) for item in child)
            else:
                nested.append(_origin_info(child, _seen=_seen))
        for info in nested:
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or bool(info["uncertainty"])
        raw_path = value.get("path_count", value.get("aggregate_path_count"))
        raw_source = value.get("source_count", value.get("aggregate_source_count"))
        try:
            path_count = max(1, int(raw_path)) if raw_path is not None else 1
        except (TypeError, ValueError):
            path_count = 1
        try:
            source_count = max(1, int(raw_source)) if raw_source is not None else path_count
        except (TypeError, ValueError):
            source_count = path_count
        if nested:
            path_count = max(path_count, max(int(item["path_count"]) for item in nested))
            source_count = max(source_count, max(int(item["source_count"]) for item in nested))
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
            "path_count": path_count,
            "source_count": source_count,
        }
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        nested = [_origin_info(item, _seen=_seen) for item in value]
        if nested:
            ids: list[Any] = []
            classes: list[Any] = []
            uncertain = False
            path_count = 1
            source_count = 1
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or bool(info["uncertainty"])
                path_count = max(path_count, int(info["path_count"]))
                source_count = max(source_count, int(info["source_count"]))
            return {
                "source_ids": _unique(ids),
                "error_classes": _unique(classes),
                "uncertainty": uncertain,
                "path_count": path_count,
                "source_count": source_count,
            }
    return {
        "source_ids": [],
        "error_classes": [],
        "uncertainty": False,
        "path_count": 1,
        "source_count": 1,
    }


def _origin_ref(value: Any) -> dict[str, Any]:
    info = _origin_info(value)
    return {
        "source_ids": list(info["source_ids"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "path_count": int(info["path_count"]),
        "source_count": int(info["source_count"]),
    }


def _new_origin(operation: str, *inputs: Any, **metadata: Any) -> dict[str, Any]:
    infos = [_origin_info(item) for item in inputs]
    ids: list[Any] = []
    classes: list[Any] = []
    uncertain = bool(metadata.get("uncertainty", False))
    if operation in {"product", "multiply"} or operation.endswith("_product"):
        path_count = 1
        source_count = 1
        for info in infos:
            path_count *= int(info["path_count"])
            source_count *= int(info["source_count"])
    elif operation in {"aggregate", "add", "sum"}:
        path_count = 0
        source_count = 0
        for info in infos:
            path_count += int(info["path_count"])
            source_count += int(info["source_count"])
    else:
        path_count = int(infos[0]["path_count"]) if infos else 1
        source_count = int(infos[0]["source_count"]) if infos else 1
    for info in infos:
        ids.extend(info["source_ids"])
        classes.extend(info["error_classes"])
        uncertain = uncertain or bool(info["uncertainty"])
    if "source_ids" in metadata:
        ids.extend(_as_items(metadata["source_ids"]))
    if "error_classes" in metadata:
        classes.extend(_as_items(metadata["error_classes"]))
    origin: dict[str, Any] = {
        "operation": operation,
        "source_ids": _unique(ids),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "path_count": max(1, int(path_count)),
        "source_count": max(1, int(source_count)),
        "source_id_count": len(_unique(ids)),
    }
    if inputs:
        origin["parents"] = tuple(_origin_ref(item) for item in inputs)
    for key, value in metadata.items():
        if key not in {"source_ids", "error_classes", "uncertainty"}:
            origin[key] = value
    return origin


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not isinstance(power, int) or isinstance(power, bool) for power in value)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(power) for power in value)  # type: ignore[return-value]


def _check_sector(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    required = {"powers", "coefficient", "origin"}
    if not required.issubset(value):
        raise ValueError(f"{label} requires powers, coefficient, and origin")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    powers = _check_powers(value["powers"], f"{label}.powers")
    info = _origin_info(value["origin"])
    if not info["source_ids"]:
        raise ValueError(f"{label}.origin must retain source_ids")
    return {"powers": powers, "coefficient": value["coefficient"], "origin": value["origin"]}


def _term_list(value: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_check_sector(item, f"{label}[{index}]") for index, item in enumerate(values)]


def _order_key(value: Any, label: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"boolean is not an ordinary-Z order in {label}")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(f"invalid ordinary-Z order {value!r} in {label}") from exc
    raise ValueError(f"ordinary-Z order required in {label}")


def _orders(value: Any, label: str) -> dict[int, Any]:
    if isinstance(value, Mapping):
        result: dict[int, Any] = {}
        for raw_key, item in value.items():
            key = _order_key(raw_key, label)
            if key in result:
                raise ValueError(f"duplicate ordinary-Z order {key} in {label}")
            result[key] = item
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        result = {index: item for index, item in enumerate(value)}
    else:
        raise TypeError(f"{label} must be an ordinary-Z order mapping")
    if set(result) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly Z0 and Z1")
    return result


def _normalize_leading(value: Any) -> dict[int, list[dict[str, Any]]]:
    raw = _orders(value, "leading")
    return {order: _term_list(raw[order], f"leading[{order}]") for order in ORDINARY_Z_ORDERS}


def _normalize_first_word(value: Any) -> dict[int, dict[str, list[dict[str, Any]]]]:
    raw = _orders(value, "first_word")
    result: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in ORDINARY_Z_ORDERS:
        branches = raw[order]
        if not isinstance(branches, Mapping) or set(branches) != set(BRANCHES):
            raise ValueError(f"first_word[{order}] must contain only core and local branches")
        result[order] = {
            branch: _term_list(branches[branch], f"first_word[{order}][{branch}]")
            for branch in BRANCHES
        }
    return result


def _validate_source_atoms(source_atoms: Any) -> None:
    """Validate names only; atom values are formal and never evaluated."""

    names: set[Any] = set()
    if isinstance(source_atoms, Mapping):
        names = set(source_atoms)
    elif isinstance(source_atoms, Sequence) and not isinstance(source_atoms, (str, bytes, bytearray)):
        for item in source_atoms:
            if isinstance(item, Mapping):
                name = item.get("name", item.get("atom", item.get("power")))
                if name is not None:
                    names.add(name)
            else:
                names.add(item)
    else:
        raise ValueError("source_atoms must name the four formal source atoms")
    if names != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda")


def _group_key(term: Mapping[str, Any]) -> tuple[Any, bool, tuple[Any, ...]]:
    info = _origin_info(term["origin"])
    classes = tuple(sorted((_freeze(item) for item in info["error_classes"]), key=repr))
    return tuple(term["powers"]), bool(info["uncertainty"]), classes


def _aggregate(c: Any, terms: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Outward-add a product stream using the directed grouping key."""

    # ``c`` is part of the source adapter API; interval addition is delegated
    # to each caller-owned coefficient object, preserving its directed context.
    del c
    groups: dict[tuple[Any, bool, tuple[Any, ...]], dict[str, Any]] = {}
    for term in terms:
        info = _origin_info(term["origin"])
        key = _group_key(term)
        group = groups.get(key)
        if group is None:
            group = {
                "powers": tuple(term["powers"]),
                "coefficient": term["coefficient"],
                "ids": {},
                "classes": {},
                "uncertainty": False,
                "path_count": 0,
                "source_count": 0,
                "merged_term_count": 0,
            }
            groups[key] = group
        else:
            group["coefficient"] = group["coefficient"] + term["coefficient"]
        for value in info["source_ids"]:
            group["ids"][_freeze(value)] = value
        for value in info["error_classes"]:
            group["classes"][_freeze(value)] = value
        group["uncertainty"] = group["uncertainty"] or bool(info["uncertainty"])
        group["path_count"] += int(info["path_count"])
        group["source_count"] += int(info["source_count"])
        group["merged_term_count"] += 1
    return [
        {
            "powers": group["powers"],
            "coefficient": group["coefficient"],
            "origin": {
                "operation": "aggregate",
                "source_ids": list(group["ids"].values()),
                "error_classes": list(group["classes"].values()),
                "uncertainty": group["uncertainty"],
                "path_count": group["path_count"],
                "source_count": group["source_count"],
                "merged_term_count": group["merged_term_count"],
                "symbolic_cancellation_claimed": False,
            },
        }
        for group in groups.values()
    ]


def _physical_product(
    leading: Mapping[str, Any],
    word: Mapping[str, Any],
    swirl_factor: Any,
    *,
    ordinary_order: int,
    branch: str,
    derivative_piece: str,
) -> dict[str, Any]:
    powers = tuple(a + b for a, b in zip(leading["powers"], word["powers"]))
    powers = (powers[0] + 1, powers[1], powers[2], powers[3] - 2)
    return {
        "powers": powers,
        "coefficient": leading["coefficient"] * word["coefficient"] * swirl_factor,
        "origin": _new_origin(
            "first_word_swirl_m15_product",
            leading,
            word,
            ordinary_Z_order=ordinary_order,
            branch=branch,
            derivative_piece=derivative_piece,
            external_formal_powers=(1, 0, 0, -2),
            coefficient_outer_factor="-5*sc",
            radial_factor="exp(2*hb*s), s=sc*q/4",
            eta_profile="1-sigma(4q-2)",
            eta_integral="5/8",
            physical_epsilon_scale="Lambda^-2",
            symbolic_cancellation_claimed=False,
        ),
    }


def _product_stream(
    leading: Sequence[Mapping[str, Any]],
    words: Sequence[Mapping[str, Any]],
    swirl_factor: Any,
    *,
    ordinary_order: int,
    branch: str,
    derivative_piece: str,
) -> Iterable[dict[str, Any]]:
    for left in leading:
        for right in words:
            yield _physical_product(
                left,
                right,
                swirl_factor,
                ordinary_order=ordinary_order,
                branch=branch,
                derivative_piece=derivative_piece,
            )


def build(
    c: Any,
    *,
    sc: Any,
    leading: Any,
    first_word: Any,
    source_family: Any,
    source_atoms: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Build finite first-word swirl ``m15`` sectors through ordinary ``Z1``.

    ``leading[0]`` and ``leading[1]`` are already-physical ``Flead`` and
    ``Flead_Z`` sectors.  ``first_word[order]`` contains separate raw
    physical ``F1_W1`` core and local sectors.  Ordinary ``Z1`` uses
    ``Flead_Z*F1_W1 + Flead*F1_W1_Z``.  Every product receives the signed
    factor ``-5*sc*[1,1/(1-sc/4)]`` and formal powers ``hb^1*Lambda^-2``.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    _validate_source_atoms(source_atoms)
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and cutoff certificates are required")
    required_certificates = (
        "hb_strictly_between_zero_and_half",
        "sc_phase_Z_independent",
        "hb_epsilon_Z_independent",
        "original_cutoff_mass_5_over_8",
    )
    missing = [name for name in required_certificates if certificate.get(name) is not True]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))

    sc_iv = c.mpf(sc)
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    one = c.mpf(1)
    sc_lower, sc_upper = _scalar_bounds(sc_iv)
    if sc_lower <= 0 or sc_upper >= 1 / 1024:
        raise ValueError("sc must satisfy 0 < sc < 1/1024")
    if sc_lower != sc_upper:
        raise ValueError("sc must be exact (a point interval)")

    leading_rows = _normalize_leading(leading)
    first_rows = _normalize_first_word(first_word)

    # The caller certifies 0 < hb < 1/2.  Thus 0 <= 2*hb*s <= sc/4 and the
    # positive radial factor is bounded by exp(x) <= 1/(1-x).
    denominator = one - sc_iv / c.mpf(4)
    denominator_lower, _ = _scalar_bounds(denominator)
    if denominator_lower <= 0:
        raise ValueError("invalid positive rho-squared denominator")
    exp_two_hb_s = _make_interval(c, one, (one / denominator).b)
    eta_integral = c.mpf(5) / c.mpf(8)
    outer_weight = c.mpf(-5) * sc_iv
    swirl_factor = outer_weight * exp_two_hb_s

    swirl: dict[int, dict[str, list[dict[str, Any]]]] = {}
    for order in ORDINARY_Z_ORDERS:
        swirl[order] = {}
        for branch in BRANCHES:
            if order == 0:
                stream = _product_stream(
                    leading_rows[0],
                    first_rows[0][branch],
                    swirl_factor,
                    ordinary_order=0,
                    branch=branch,
                    derivative_piece="Flead*F1_W1",
                )
            else:
                first_piece = _product_stream(
                    leading_rows[1],
                    first_rows[0][branch],
                    swirl_factor,
                    ordinary_order=1,
                    branch=branch,
                    derivative_piece="Flead_Z*F1_W1",
                )
                second_piece = _product_stream(
                    leading_rows[0],
                    first_rows[1][branch],
                    swirl_factor,
                    ordinary_order=1,
                    branch=branch,
                    derivative_piece="Flead*F1_W1_Z",
                )

                def combined() -> Iterable[dict[str, Any]]:
                    yield from first_piece
                    yield from second_piece

                stream = combined()
            swirl[order][branch] = _aggregate(c, stream)

    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "swirl": swirl,
        "coefficient_weight": outer_weight,
        "exp_two_hb_s_interval": exp_two_hb_s,
        "rho_squared_interval": c.mpf(16) * exp_two_hb_s,
        "rho_squared_exp_interval": exp_two_hb_s,
        "eta_profile": "1-sigma(4q-2)",
        "eta_integral": eta_integral,
        "eta_integral_exact": True,
        "cutoff_mass": eta_integral,
        "positive_cutoff_mass_exact": True,
        "formal_external_powers": (1, 0, 0, -2),
        "formal_hb_measure": True,
        "physical_epsilon_scale": "Lambda^-2",
        "leading_role": "physical Flead",
        "leading_is_entry_point": False,
        "first_word_role": "raw physical F1 W1",
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "aggregation_key": "(powers, uncertainty, error_classes)",
        "aggregation_is_outward_interval_addition_only": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "finite_first_word_only": True,
        "finite_first_word_swirl_m15_computed": True,
        "m15_computed": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "scope": "finite first-word physical leading first-collar swirl m15 only; core and local F1_W1 retained separately",
    }
