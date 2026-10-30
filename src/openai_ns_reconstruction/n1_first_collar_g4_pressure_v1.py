"""Raw physical ``g4`` pressure sectors on the complete first collar.

This is a deliberately small source adapter.  The caller supplies the exact
ordinary-``Z`` whole-phase enclosures for the raw physical source
``g4 = 2*epsilon*x*N1p``.  The adapter does not evaluate ``N1p`` (or any
other source), replay a producer, read a ledger, or consult the older W0/W1
builders.  It only attaches the leading first-collar pressure weight

    ``-epsilon*hb*sc*Integral(exp(3*hb*sc*q/8)*g4(q,Z), q=0..1)``.

The supplied ``g4`` coefficient is a uniform box on ``q in [0, 1]``.  Hence
the complete-collar integral is enclosed by multiplying that box by the
positive interval ``[1, 1/(1 - 3*sc/16)]`` and by the signed factor ``-sc``.
The physical epsilon is represented by the external formal power
``Lambda^-1`` and the measure contributes one external ``hb``.  Ordinary Z1
is a caller-supplied ``g4[1]`` row; no geometry differentiation is performed
here.

Each input term produces exactly one output sector.  Terms are intentionally
not merged: source IDs, uncertainty labels, and error classes therefore stay
attached to the identity of the original enclosure.  This module contains no
eta/velocity cutoff, endpoint term, tail, or total-m15 assembly.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
REQUIRED_CERTIFICATES = (
    "hb_strictly_between_zero_and_half",
    "sc_phase_Z_independent",
    "hb_epsilon_Z_independent",
    "g4_equals_two_epsilon_x_N1p",
)

__all__ = [
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "REQUIRED_CERTIFICATES",
    "build",
]


def _is_interval(value: Any) -> bool:
    """Return whether *value* is a caller-owned directed interval."""

    return hasattr(value, "_mpi_") and hasattr(value, "a") and hasattr(value, "b")


def _raw_endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct a directed interval without midpointing either endpoint."""

    return c.make_mpf((_raw_endpoint(lower, False), _raw_endpoint(upper, True)))


def _scalar_bounds(value: Any) -> tuple[Any, Any]:
    """Expose raw directed endpoints for guards only.

    ``mpmath.iv`` stores endpoint tuples on ``_mpi_``.  The source and all
    output coefficients remain in the caller's interval context; these scalar
    values are used only to check positivity and exactness of ``sc``.
    """

    raw = getattr(value, "_mpi_", None)
    if raw is None:
        return value, value
    # Importing the arithmetic library is not importing a repository program.
    # ``mp.make_mpf`` converts the low-level directed endpoints exactly enough
    # for the validation comparisons below.
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
    """Collect provenance without evaluating or combining a source term."""

    if _seen is None:
        _seen = set()
    marker = id(value)
    if marker in _seen:
        return {
            "source_ids": [],
            "error_classes": [],
            "uncertainty": False,
        }

    if isinstance(value, Mapping):
        _seen.add(marker)
        raw_ids = value.get("source_ids", value.get("sources"))
        if raw_ids is None:
            raw_ids = value.get("source_id", value.get("source"))
        ids = list(_as_items(raw_ids))
        raw_classes = value.get(
            "error_classes",
            value.get("uncertainty_class", value.get("error_class")),
        )
        classes = list(_as_items(raw_classes))
        operation = str(value.get("operation", "")).lower()
        error_kind = str(value.get("error_kind", value.get("kind", ""))).lower()
        uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
        uncertain = uncertain or any(
            token in operation or token in error_kind
            for token in ("cauchy", "symmetric", "uncertainty")
        )
        for key in ("origin", "metadata", "parents"):
            child = value.get(key)
            if child is None or child is value:
                continue
            if isinstance(child, Sequence) and not isinstance(child, (str, bytes, bytearray)):
                nested = [_origin_info(item, _seen=_seen) for item in child]
            else:
                nested = [_origin_info(child, _seen=_seen)]
            for info in nested:
                ids.extend(info["source_ids"])
                classes.extend(info["error_classes"])
                uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": bool(uncertain),
        }

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        ids: list[Any] = []
        classes: list[Any] = []
        uncertain = False
        for item in value:
            info = _origin_info(item, _seen=_seen)
            ids.extend(info["source_ids"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or bool(info["uncertainty"])
        return {
            "source_ids": _unique(ids),
            "error_classes": _unique(classes),
            "uncertainty": uncertain,
        }

    return {"source_ids": [], "error_classes": [], "uncertainty": False}


def _check_powers(value: Any, label: str) -> tuple[int, int, int, int]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes, bytearray))
        or len(value) != 4
        or any(not isinstance(power, int) or isinstance(power, bool) for power in value)
    ):
        raise ValueError(f"{label} must contain four integer formal powers")
    return tuple(int(power) for power in value)  # type: ignore[return-value]


def _check_term(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{label} must be a sector mapping")
    if not {"powers", "coefficient", "origin"}.issubset(value):
        raise ValueError(f"{label} requires powers, coefficient, and origin")
    if not _is_interval(value["coefficient"]):
        raise TypeError(f"{label}.coefficient must be a directed interval")
    powers = _check_powers(value["powers"], f"{label}.powers")
    info = _origin_info(value["origin"])
    ids = list(info["source_ids"])
    for key in ("input_id", "source_id"):
        if value.get(key) is not None:
            ids.extend(_as_items(value[key]))
    if not ids:
        raise ValueError(f"{label}.origin must retain source_ids")
    return {
        "powers": powers,
        "coefficient": value["coefficient"],
        "origin": value["origin"],
        "source_ids": _unique(ids),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "input_id": value.get("input_id"),
    }


def _term_list(value: Any, label: str) -> list[dict[str, Any]]:
    if isinstance(value, Mapping) and "powers" in value:
        values = [value]
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values:
        raise ValueError(f"{label} must not be empty")
    return [_check_term(item, f"{label}[{index}]") for index, item in enumerate(values)]


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
        raise TypeError(f"{label} must provide ordinary-Z order rows")
    if set(result) != set(ORDINARY_Z_ORDERS):
        raise ValueError(f"{label} must provide exactly Z0 and Z1")
    return result


def _normalize_g4(value: Any) -> dict[int, list[dict[str, Any]]]:
    rows = _orders(value, "g4")
    return {
        order: _term_list(rows[order], f"g4[{order}]")
        for order in ORDINARY_Z_ORDERS
    }


def _validate_source_atoms(source_atoms: Any) -> None:
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
        raise ValueError(
            "source_atoms must contain exactly hb, F0_anchor, Pstar_squared, Lambda"
        )


def _validate_certificate(certificate: Any) -> None:
    if not isinstance(certificate, Mapping):
        raise ValueError("caller geometry and g4 identity certificates are required")
    missing = [name for name in REQUIRED_CERTIFICATES if certificate.get(name) is not True]
    if missing:
        raise ValueError("missing caller certificate: " + ", ".join(missing))


def _prepare_output_term(
    term: Mapping[str, Any],
    signed_exp_factor: Any,
    *,
    ordinary_order: int,
    input_index: int,
) -> dict[str, Any]:
    """Attach the pressure factor to one input term, without aggregation."""

    powers = tuple(term["powers"])
    output_powers = (
        powers[0] + 1,
        powers[1],
        powers[2],
        powers[3] - 1,
    )
    origin_value = term["origin"]
    origin = dict(origin_value) if isinstance(origin_value, Mapping) else {}
    origin["source_ids"] = list(term["source_ids"])
    origin["error_classes"] = list(term["error_classes"])
    origin["uncertainty"] = bool(term["uncertainty"])
    origin.update(
        {
            "operation": "first_collar_g4_pressure_product",
            "ordinary_Z_order": ordinary_order,
            "g4_input_term_index": input_index,
            "external_formal_powers": (1, 0, 0, -1),
            "signed_outer_factor": "-sc",
            "positive_exponential_factor": "exp(3*hb*sc*q/8)",
            "uniform_q_box": True,
            "z1_geometry_differentiated": False,
            "symbolic_cancellation_claimed": False,
            "unmerged_input_term": True,
        }
    )
    return {
        "powers": output_powers,
        "coefficient": term["coefficient"] * signed_exp_factor,
        "origin": origin,
    }


def build(
    c: Any,
    *,
    sc: Any,
    g4: Any,
    source_family: Any,
    source_atoms: Any,
    certificate: Any,
) -> dict[str, Any]:
    """Build the leading physical ``N1p`` first-collar pressure sectors.

    ``g4[0]`` and ``g4[1]`` are the caller's uniform whole-phase boxes for
    raw physical ``g4`` and its ordinary-Z1 row.  The output contains one
    pressure sector for each input term, in input order.  Ordinary Z1 is
    intentionally consumed directly from ``g4[1]`` rather than obtained by
    differentiating any geometry in this module.
    """

    if not hasattr(c, "mpf") or not hasattr(c, "make_mpf"):
        raise TypeError("c must be the caller's MPIntervalContext")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")
    _validate_source_atoms(source_atoms)
    _validate_certificate(certificate)

    sc_iv = c.mpf(sc)
    if not _is_interval(sc_iv):
        raise TypeError("sc must be a directed interval")
    sc_lower, sc_upper = _scalar_bounds(sc_iv)
    if sc_lower <= 0:
        raise ValueError("sc must be strictly positive")
    if sc_lower != sc_upper:
        raise ValueError("sc must be exact (a point interval)")

    # 0 <= 3*hb*sc*q/8 <= 3*sc/16 under the caller-certified
    # 0 < hb < 1/2 and 0 <= q <= 1.  The denominator guard is the only
    # numerical restriction needed for the geometric exponential bound.
    one = c.mpf(1)
    denominator = one - c.mpf(3) * sc_iv / c.mpf(16)
    denominator_lower, _ = _scalar_bounds(denominator)
    if denominator_lower <= 0:
        raise ValueError("invalid positive exponential denominator")
    exp_interval = _make_interval(c, one, (one / denominator).b)
    signed_factor = c.mpf(-1) * sc_iv
    signed_exp_factor = signed_factor * exp_interval

    rows = _normalize_g4(g4)
    pressure: dict[int, list[dict[str, Any]]] = {}
    for order in ORDINARY_Z_ORDERS:
        pressure[order] = [
            _prepare_output_term(
                term,
                signed_exp_factor,
                ordinary_order=order,
                input_index=index,
            )
            for index, term in enumerate(rows[order])
        ]

    return {
        "source_family": source_family,
        "source_atoms": source_atoms,
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "pressure": pressure,
        "uniform_g4_boxes": True,
        "q_domain": "[0,1]",
        "coefficient_weight": signed_factor,
        "signed_outer_factor": signed_factor,
        "signed_factor": signed_factor,
        "exp_3hb_sc_q_over_8_interval": exp_interval,
        "positive_exponential_interval": exp_interval,
        "exp_interval": exp_interval,
        "signed_pressure_factor_interval": signed_exp_factor,
        "formal_external_powers": (1, 0, 0, -1),
        "formal_hb_measure": True,
        "physical_epsilon_scale": "Lambda^-1",
        "pressure_integrand": "-epsilon*hb*sc*exp(3*hb*sc*q/8)*g4(q,Z)",
        "raw_g4_identity": "g4=2*epsilon*x*N1p",
        "ordinary_Z1_uses_g4_Z1_directly": True,
        "z1_geometry_differentiated": False,
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "source_ids_union_retained": True,
        "error_classes_union_retained": True,
        "uncertainty_retained": True,
        "terms_merged": False,
        "input_terms_kept_one_to_one": True,
        "pressure_terms_unmerged": True,
        "source_terms_unmerged": True,
        "velocity_cutoff_applied": False,
        "eta_cutoff_applied": False,
        "eta_cutoff_used": False,
        "explicit_pressure_endpoint_added": False,
        "arbitrary_endpoint_cancellation_claimed": False,
        "allword_tail_added": False,
        "allword_tail_included": False,
        "tail_added": False,
        "m15_computed": False,
        "total_m15_computed": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "point_value_computed": False,
        "numeric_point_evaluation": False,
        "no_point_value": True,
        "w0_w1_consulted": False,
        "ledger_212_consulted": False,
        "scope": "leading raw physical N1p first-collar pressure m15 contribution only",
    }
