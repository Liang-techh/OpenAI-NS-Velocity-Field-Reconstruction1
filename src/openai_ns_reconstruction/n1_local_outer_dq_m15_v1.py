"""Conditional raw-Q local outer ``D`` first-potential physical m15 leaf.

The leaf consumes already saved, whole-collar ``F1`` sectors.  In
particular, ``low_word`` is the caller-owned local
``W1_F = G0(C0 W0 + B1 W0_Z)`` packet; it is never reconstructed here.  The
axis potential and its ordinary-Z derivative are also caller-owned sectors.

For ``x(q) = 2 exp(hb*sc*q/8)`` and ``d = hb*sc*q/4`` the positive K0 mass is
bounded on the whole phase by

``hb**2*(sc/4)**2/2 * [0, exp(sc/8)]``.

The ``hb**2`` is a formal power shift, not a numerical value.  The interval
``[0, exp(sc/8)]`` is a whole-phase upper enclosure and is deliberately not
used as a point value.  Ordinary Z1 is kept as the two Leibniz pieces
``chi_Z*W1_Z0 + chi*W1_Z1``.  The physical swirl projection then adds the
single factor ``hb*Lambda**-2`` and ``-5*sc*[1, exp(sc/4)]``.  Thus every
projected sector has the source powers plus ``(3, 0, 0, -2)``.

This is a raw-Q (C0/B1) partial only.  It does not re-enter Green, W1 or P,
does not apply the physical P compensation, does not include other inner
corrections or a higher potential tail, and does not write a ledger or claim
a complete m15.  Only the F channel is emitted; P1 is not pressure here.
The module has no producer/archive imports and performs no work on import.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_local_outer_dq_m15_v1"
COMPONENT = "local_outer_dq_m15"
SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
FIELDS = ("F1",)
KERNEL_POWER_SHIFT = (2, 0, 0, 0)
PHYSICAL_POWER_SHIFT = (1, 0, 0, -2)
TOTAL_POWER_SHIFT = (3, 0, 0, -2)
Z1_PATHS = (
    "leading_Z*chi*W1_0",
    "leading*chi_Z*W1_0",
    "leading*chi*W1_1",
)

_MISSING = object()
_FIELD_KEYS = ("F1", "f1", "field_F1", "swirl", "field", 1, "1")
_CHI_KEYS = ("chi", "axis_chi", "chi_jets", "ordinary_Z", "jets", "rows")


__all__ = [
    "SCHEMA",
    "COMPONENT",
    "SOURCE_POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "FIELDS",
    "KERNEL_POWER_SHIFT",
    "PHYSICAL_POWER_SHIFT",
    "TOTAL_POWER_SHIFT",
    "Z1_PATHS",
    "kernel_mass",
    "build",
    "build_local_outer_dq_m15",
    "evaluate",
]


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _items(value: Any) -> list[Any]:
    if value is None or value is _MISSING:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if _is_sequence(value):
        return list(value)
    return [value]


def _raw_endpoint(value: Any, upper: bool) -> Any:
    """Return one endpoint without rounding it through a midpoint."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    if hasattr(value, "a") and hasattr(value, "b"):
        candidate = value.b if upper else value.a
        nested = getattr(candidate, "_mpi_", None)
        if nested is not None:
            return nested[1 if upper else 0]
        return candidate
    return value


def _make_interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct from exact raw endpoints without lower-precision rounding."""
    def endpoint(value: Any, upper: bool) -> Any:
        raw = _raw_endpoint(value, upper)
        if isinstance(raw, tuple) and len(raw) == 4 and all(type(x) is int for x in raw):
            return raw
        return _raw_endpoint(c.mpf(raw), upper)
    return c.make_mpf((endpoint(lower, False), endpoint(upper, True)))


def _as_interval(c: Any, value: Any, label: str) -> Any:
    if hasattr(value, "_mpi_"):
        return _make_interval(c, value, value)
    if _is_mapping(value):
        for key in ("interval", "enclosure", "box", "value", "record"):
            if key in value:
                return _as_interval(c, value[key], label)
        lower = value.get("lower", value.get("lo", _MISSING))
        upper = value.get("upper", value.get("hi", _MISSING))
        if lower is not _MISSING and upper is not _MISSING:
            return _make_interval(c, lower, upper)
        raise TypeError(f"{label} is not interval-compatible")
    if _is_sequence(value) and len(value) == 2:
        return _make_interval(c, value[0], value[1])
    try:
        return c.mpf(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _ordinary_endpoint(value: Any, upper: bool) -> Any:
    """Convert an endpoint to an ordinary mpmath scalar for guards only."""

    from mpmath import mp

    raw = _raw_endpoint(value, upper)
    try:
        return mp.make_mpf(raw)
    except (TypeError, ValueError):
        return mp.mpf(raw)


def _lower(value: Any) -> Any:
    return _ordinary_endpoint(value, False)


def _upper(value: Any) -> Any:
    return _ordinary_endpoint(value, True)


def _zero(c: Any) -> Any:
    return c.mpf(0)


def _one(c: Any) -> Any:
    return c.mpf(1)


def _symmetric(c: Any, value: Any) -> Any:
    lo = c.make_mpf((_raw_endpoint(value, False), _raw_endpoint(value, False)))
    hi = c.make_mpf((_raw_endpoint(value, True), _raw_endpoint(value, True)))
    alo = -lo if lo < 0 else lo
    ahi = -hi if hi < 0 else hi
    radius = alo if alo >= ahi else ahi
    return _make_interval(c, -radius, radius)


def _singleton(c: Any, value: Any, *, upper: bool) -> Any:
    raw = _raw_endpoint(value, upper)
    if not (isinstance(raw, tuple) and len(raw) == 4 and all(type(x) is int for x in raw)):
        raw = _raw_endpoint(c.mpf(raw), upper)
    return c.make_mpf((raw, raw))


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    seen: set[Any] = set()
    for value in values:
        try:
            key = (type(value).__name__, value)
            hash(key)
        except TypeError:
            key = (type(value).__name__, repr(value))
        if key not in seen:
            seen.add(key)
            result.append(value)
    return result


def _lookup(mapping: Mapping[Any, Any], *keys: Any) -> Any:
    for key in keys:
        try:
            if key in mapping:
                return mapping[key]
        except (TypeError, ValueError):
            continue
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for key in keys:
        if str(key).lower() in lowered:
            return lowered[str(key).lower()]
    return _MISSING


def _origin_info(value: Any, *, _seen: set[int] | None = None) -> dict[str, Any]:
    """Collect source IDs, paths, classes, and uncertainty recursively."""

    if _seen is None:
        _seen = set()
    if id(value) in _seen:
        return {"source_ids": [], "source_paths": [], "error_classes": [],
                "uncertainty": False, "error_only": False}
    if not _is_mapping(value):
        if _is_sequence(value):
            merged = {"source_ids": [], "source_paths": [], "error_classes": [],
                      "uncertainty": False, "error_only": False}
            for item in value:
                info = _origin_info(item, _seen=_seen)
                for key in ("source_ids", "source_paths", "error_classes"):
                    merged[key].extend(info[key])
                merged["uncertainty"] = merged["uncertainty"] or info["uncertainty"]
                merged["error_only"] = merged["error_only"] or info["error_only"]
            for key in ("source_ids", "source_paths", "error_classes"):
                merged[key] = _unique(merged[key])
            return merged
        return {"source_ids": [], "source_paths": [], "error_classes": [],
                "uncertainty": False, "error_only": False}

    _seen.add(id(value))
    ids: list[Any] = []
    paths: list[Any] = []
    classes: list[Any] = []
    for key in ("source_ids", "sources", "source_id", "source", "input_id"):
        if key in value:
            ids.extend(_items(value[key]))
    for key in ("source_paths", "source_path", "paths", "path"):
        if key in value:
            paths.extend(_items(value[key]))
    for key in ("error_classes", "error_class", "uncertainty_class"):
        if key in value:
            classes.extend(_items(value[key]))
    operation = str(value.get("operation", "")).lower()
    kind = str(value.get("kind", value.get("error_kind", ""))).lower()
    error_only = bool(value.get("error_only", value.get("symmetric_error", False)))
    uncertain = bool(value.get("uncertainty", value.get("is_uncertainty", False)))
    uncertain = uncertain or error_only or any(
        token in operation or token in kind
        for token in ("error", "uncertainty", "symmetric", "tail", "remainder")
    )
    for key in ("origin", "metadata", "parents"):
        child = value.get(key)
        if child is None or child is value:
            continue
        for item in (_items(child) if _is_sequence(child) else [child]):
            info = _origin_info(item, _seen=_seen)
            ids.extend(info["source_ids"])
            paths.extend(info["source_paths"])
            classes.extend(info["error_classes"])
            uncertain = uncertain or info["uncertainty"]
            error_only = error_only or info["error_only"]
    return {
        "source_ids": _unique(ids),
        "source_paths": _unique(paths),
        "error_classes": _unique(classes),
        "uncertainty": bool(uncertain),
        "error_only": bool(error_only),
    }


def _powers(record: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    raw = _lookup(record, "powers", "formal_powers", "power_tuple")
    if raw is _MISSING:
        raw = {
            "hb": record.get("hb_power", 0),
            "F0_anchor": record.get("F0_power", record.get("F0_anchor_power", 0)),
            "Pstar_squared": record.get("Pstar_squared_power", record.get("Pstar_power", 0)),
            "Lambda": record.get("Lambda_power", record.get("lambda_power", 0)),
        }
    if _is_mapping(raw):
        values = []
        aliases = {
            "hb": ("hb", 0),
            "F0_anchor": ("F0_anchor", "F0", 1),
            "Pstar_squared": ("Pstar_squared", "Pstar^2", "Pstar2", 2),
            "Lambda": ("Lambda", "lambda", 3),
        }
        for name, names in aliases.items():
            item = _lookup(raw, *names)
            if item is _MISSING:
                raise ValueError(f"{label}.powers is missing {name}")
            values.append(item)
        raw = values
    if not _is_sequence(raw) or len(raw) != 4:
        raise ValueError(f"{label}.powers must contain four integer formal powers")
    if any(not isinstance(item, int) or isinstance(item, bool) for item in raw):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return tuple(int(item) for item in raw)  # type: ignore[return-value]


def _normalise_term(c: Any, value: Any, label: str) -> dict[str, Any]:
    if not _is_mapping(value):
        raise TypeError(f"{label} must be a sector mapping")
    coefficient = _lookup(value, "coefficient", "coefficient_interval", "value", "interval")
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    coefficient = _as_interval(c, coefficient, f"{label}.coefficient")
    powers = _powers(value, label)
    info = _origin_info(value)
    if not info["source_ids"]:
        raise ValueError(f"{label} must retain at least one source_id")
    result = dict(value)
    result.update({
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": list(powers),
        "formal_powers": list(powers),
        "source_ids": list(info["source_ids"]),
        "source_paths": list(info["source_paths"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": bool(info["uncertainty"]),
        "error_only": bool(info["error_only"]),
        "symmetric_error": bool(info["error_only"] or value.get("symmetric_error", False)),
    })
    return result


def _term_list(value: Any, label: str, c: Any) -> list[dict[str, Any]]:
    if _is_mapping(value):
        nested = _lookup(value, "sectors", "terms", "source_sectors", "records")
        if nested is not _MISSING:
            return _term_list(nested, label, c)
        if _lookup(value, "coefficient", "coefficient_interval", "value", "interval") is not _MISSING:
            return [_normalise_term(c, value, label)]
        values = list(value.values())
    elif _is_sequence(value):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector sequence")
    if not values or not all(_is_mapping(item) for item in values):
        raise TypeError(f"{label} must contain sector mappings")
    return [_normalise_term(c, item, f"{label}[{index}]") for index, item in enumerate(values)]


def _order_key(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int) and value in (0, 1):
        return value
    text = str(value).lower()
    if text in {"0", "z0", "ordinary_z0"}:
        return 0
    if text in {"1", "z1", "ordinary_z1"}:
        return 1
    return None


def _unwrap_field(value: Any, label: str, field: str = "F1") -> Any:
    current = value
    if _is_mapping(current):
        for key in ("local_signed_sectors", "signed_sectors", "ordinary_Z", "jets", "rows"):
            nested = _lookup(current, key)
            if nested is not _MISSING:
                current = nested
                break
        if _is_mapping(current):
            is_term = _lookup(current, "powers", "coefficient") is not _MISSING
            has_order = any(_order_key(key) is not None for key in current)
            if not is_term and not has_order:
                for key in _FIELD_KEYS:
                    nested = _lookup(current, key)
                    if nested is not _MISSING:
                        current = nested
                        break
    return current


def _rows(c: Any, value: Any, label: str, *, field: str = "F1") -> dict[int, list[dict[str, Any]]]:
    current = _unwrap_field(value, label, field)
    if _is_mapping(current):
        nested = _lookup(current, field, field.lower(), "field_" + field)
        if nested is not _MISSING and _lookup(current, "powers", "coefficient") is _MISSING:
            current = nested
    result: dict[int, Any] = {}
    if _is_mapping(current):
        for key, item in current.items():
            order = _order_key(key)
            if order is not None:
                if order in result:
                    raise ValueError(f"duplicate {label} ordinary-Z order {order}")
                result[order] = item
    elif _is_sequence(current) and len(current) == 2:
        result = {0: current[0], 1: current[1]}
    else:
        raise TypeError(f"{label} must provide Z0 and Z1 sectors")
    if set(result) != {0, 1}:
        raise ValueError(f"{label} must provide exactly Z0 and Z1 sectors")
    return {
        order: _term_list(
            _unwrap_field(result[order], f"{label}.Z{order}", field),
            f"{label}.Z{order}",
            c,
        )
        for order in (0, 1)
    }


def _chi_rows(c: Any, chi: Any, chi_z: Any) -> dict[int, list[dict[str, Any]]]:
    if chi is None:
        raise ValueError("axis chi is required")
    current = chi
    if _is_mapping(current):
        for key in _CHI_KEYS:
            nested = _lookup(current, key)
            if nested is not _MISSING and nested is not current:
                current = nested
                break
    has_order = _is_mapping(current) and any(_order_key(key) is not None for key in current)
    if chi_z is not None:
        if has_order:
            raise ValueError("axis chi already contains Z1; chi_Z was supplied twice")
        return {0: _term_list(current, "axis_chi.Z0", c),
                1: _term_list(chi_z, "axis_chi.Z1", c)}
    return _rows(c, current, "axis_chi")


def _source_atoms(value: Any) -> list[Any]:
    if value is None:
        return list(SOURCE_POWER_BASIS)
    if _is_mapping(value):
        names = set(value)
    elif _is_sequence(value):
        names = set()
        for item in value:
            if _is_mapping(item):
                name = item.get("name", item.get("atom", item.get("power")))
                if name is not None:
                    names.add(name)
            else:
                names.add(item)
    else:
        raise ValueError("source_atoms must contain the four formal source atoms")
    aliases = {
        "hb": "hb",
        "F0": "F0_anchor",
        "F0_anchor": "F0_anchor",
        "Pstar^2": "Pstar_squared",
        "Pstar2": "Pstar_squared",
        "Pstar_squared": "Pstar_squared",
        "Lambda": "Lambda",
        "lambda": "Lambda",
    }
    canonical = {aliases.get(str(name), name) for name in names}
    if canonical != set(SOURCE_POWER_BASIS):
        raise ValueError("source_atoms must contain the four formal source atoms")
    return list(value.keys()) if _is_mapping(value) else list(value)


def _flag(certificate: Mapping[Any, Any], *names: str) -> bool:
    candidates: list[Mapping[Any, Any]] = [certificate]
    for key in ("geometry", "source_certificate", "proof", "proofs", "uniformity"):
        nested = certificate.get(key)
        if _is_mapping(nested):
            candidates.append(nested)
    for current in candidates:
        for name in names:
            value = _lookup(current, name)
            if value is True:
                return True
            if _is_mapping(value) and value.get("proved") is True:
                return True
    return False


def _validate_certificate(c: Any, sc: Any, certificate: Any, hb_override: Any = None) -> tuple[Any, Any, dict[str, Any]]:
    if not _is_mapping(certificate):
        raise TypeError("caller geometry/source certificate is required")
    sc_i = _as_interval(c, sc, "sc")
    if _lower(sc_i) <= 0:
        raise ValueError("sc must be strictly positive")

    hb_raw = hb_override if hb_override is not None else _lookup(certificate, "hb", "hb_interval", "h_b")
    hb_i = None if hb_raw is _MISSING else _as_interval(c, hb_raw, "certificate.hb")
    strict_proof = _flag(
        certificate,
        "hb_strictly_between_zero_and_half",
        "strict_hb_admitted",
        "strict_hb_proof",
        "strict_hb_admission_proof",
    )
    if hb_i is not None:
        strict = _lower(hb_i) > 0 and _upper(hb_i) < 0.5
        closed = _lower(hb_i) >= 0 and _upper(hb_i) <= 0.5
        if not strict and not (closed and strict_proof):
            raise ValueError("hb must prove 0<hb<1/2; closed [0,1/2] requires strict proof")
    elif not strict_proof:
        raise ValueError("certificate must provide hb semantics or strict hb proof")

    z_independent = _flag(
        certificate,
        "geometry_Z_independent",
        "geometry_z_independent",
        "hb_sc_Z_independent",
        "rho_radius_Z_independent",
    ) or (
        _flag(certificate, "sc_phase_Z_independent", "sc_Z_independent")
        and _flag(certificate, "hb_epsilon_Z_independent", "hb_Z_independent")
    )
    if not z_independent:
        raise ValueError("caller certificate must prove geometry Z independence")

    q_all_phase = _flag(
        certificate,
        "q_all_phase",
        "whole_phase_q",
        "whole_collar_q",
        "full_q_source_enclosure",
        "source_full_q_enclosure",
    )
    q_domain = _lookup(certificate, "q_domain", "q_bounds", "source_q_domain")
    q_domain_i = None
    if q_domain is not _MISSING:
        q_domain_i = _as_interval(c, q_domain, "certificate.q_domain")
        q_all_phase = q_all_phase or (_lower(q_domain_i) <= 0 and _upper(q_domain_i) >= 1)
    if not q_all_phase:
        raise ValueError("caller certificate must cover all q in [0,1]")

    cutoff = _flag(
        certificate,
        "original_cutoff_mass_5_over_8",
        "cutoff_mass_5_over_8",
        "eta_mass_5_over_8",
    )
    cutoff_value = _lookup(certificate, "cutoff_mass", "eta_mass", "original_eta_mass")
    if cutoff_value is not _MISSING:
        cutoff_i = _as_interval(c, cutoff_value, "certificate.cutoff_mass")
        cutoff = cutoff or (_lower(cutoff_i) == 0.625 and _upper(cutoff_i) == 0.625)
    if not cutoff:
        raise ValueError("caller certificate must prove cutoff mass 5/8")

    explicit_phase_flag = _lookup(certificate, "axis_chi_is_not_phase_cutoff", "axis_chi_not_phase_cutoff")
    if explicit_phase_flag is False or _flag(certificate, "chi_is_phase_cutoff", "phase_cutoff_is_axis_chi"):
        raise ValueError("axis chi cannot be the phase cutoff")

    geometry = {
        "hb": hb_i,
        "sc": sc_i,
        "q_domain": q_domain_i if q_domain_i is not None else _make_interval(c, 0, 1),
        "hb_strict_semantics": "0<hb<1/2",
        "closed_hb_enclosure_allowed": bool(strict_proof and (hb_i is None or (_lower(hb_i) >= 0 and _upper(hb_i) <= 0.5))),
        "formal_hb_not_materialized": True,
        "geometry_Z_independent": True,
        "q_all_phase": True,
        "axis_chi_is_not_phase_cutoff": True,
        "cutoff_mass": c.mpf(5) / c.mpf(8),
        "cutoff_mass_exact": True,
    }
    return sc_i, hb_i, geometry


def _source_metadata(*terms: Mapping[str, Any]) -> dict[str, Any]:
    infos = [_origin_info(term) for term in terms]
    ids = _unique([item for info in infos for item in info["source_ids"]])
    paths = _unique([item for info in infos for item in info["source_paths"]])
    classes = _unique([item for info in infos for item in info["error_classes"]])
    error_only = any(info["error_only"] for info in infos)
    uncertainty = error_only or any(info["uncertainty"] for info in infos)
    factor_names = ("leading_F", "axis_chi", "W1_F")
    by_factor: dict[str, dict[str, list[Any]]] = {}
    for name, info in zip(factor_names[-len(infos):], infos):
        by_factor[name] = {
            "source_ids": list(info["source_ids"]),
            "source_paths": list(info["source_paths"]),
            "error_classes": list(info["error_classes"]),
        }
    if "axis_chi" in by_factor:
        by_factor["chi"] = dict(by_factor["axis_chi"])
    if "W1_F" in by_factor:
        by_factor["W1"] = dict(by_factor["W1_F"])
    return {
        "source_ids": ids,
        "source_paths": paths,
        "error_classes": classes,
        "error_only": error_only,
        "uncertainty": uncertainty,
        "by_factor": by_factor,
        "source_ids_by_factor": {key: value["source_ids"] for key, value in by_factor.items()},
        "source_paths_by_factor": {key: value["source_paths"] for key, value in by_factor.items()},
        "error_classes_by_factor": {key: value["error_classes"] for key, value in by_factor.items()},
    }


def _field_product(
    c: Any,
    chi_term: Mapping[str, Any],
    source_term: Mapping[str, Any],
    kernel_factor: Any,
    *,
    chi_order: int,
    source_order: int,
    z_path: str,
) -> dict[str, Any]:
    metadata = _source_metadata(chi_term, source_term)
    coefficient = -c.mpf(2) * chi_term["coefficient"] * source_term["coefficient"] * kernel_factor
    if metadata["error_only"]:
        coefficient = _symmetric(c, coefficient)
        if "symmetric_error" not in metadata["error_classes"]:
            metadata["error_classes"].append("symmetric_error")
    powers = [chi_term["powers"][index] + source_term["powers"][index] + KERNEL_POWER_SHIFT[index] for index in range(4)]
    ordinary_order = 0 if (chi_order, source_order) == (0, 0) else 1
    origin = {
        "operation": "local_outer_D_first_potential_field",
        "source_ids": list(metadata["source_ids"]),
        "source_paths": list(metadata["source_paths"]),
        "error_classes": list(metadata["error_classes"]),
        "uncertainty": metadata["uncertainty"],
        "error_only": metadata["error_only"],
        "symmetric_error": metadata["error_only"],
        "chi_order": chi_order,
        "source_order": source_order,
        "ordinary_Z_order": ordinary_order,
        "leibniz_piece": z_path,
        "source_ids_by_factor": metadata["source_ids_by_factor"],
        "source_paths_by_factor": metadata["source_paths_by_factor"],
        "error_classes_by_factor": metadata["error_classes_by_factor"],
        "formal_hb_power_added": KERNEL_POWER_SHIFT[0],
        "kernel_mass_is_upper_bound": True,
        "kernel_mass_upper_is_not_a_value": True,
        "kernel_mass_zero_at_q_zero": True,
        "green_reapplied": False,
        "W1_reapplied": False,
        "P_reapplied": False,
        "symbolic_cancellation_claimed": False,
    }
    return {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": powers,
        "formal_powers": powers,
        "power_basis": list(SOURCE_POWER_BASIS),
        "kernel_mass": kernel_factor,
        "kernel_mass_is_coefficient_only": True,
        "kernel_mass_is_upper_bound": True,
        "kernel_mass_upper_is_not_a_value": True,
        "kernel_mass_zero_at_q_zero": True,
        "field": "F1",
        "source_kind": "raw_Q_local_outer_D_W1_F",
        "ordinary_Z_order": ordinary_order,
        "chi_order": chi_order,
        "source_order": source_order,
        "leibniz_piece": z_path,
        "source_ids": list(metadata["source_ids"]),
        "source_paths": list(metadata["source_paths"]),
        "error_classes": list(metadata["error_classes"]),
        "uncertainty": metadata["uncertainty"],
        "error_only": metadata["error_only"],
        "symmetric_error": metadata["error_only"],
        "source_ids_by_factor": metadata["source_ids_by_factor"],
        "source_paths_by_factor": metadata["source_paths_by_factor"],
        "error_classes_by_factor": metadata["error_classes_by_factor"],
        "source_components": {"axis_chi": dict(chi_term), "W1_F": dict(source_term)},
        "origin": origin,
    }


def _project_product(
    c: Any,
    leading: Mapping[str, Any],
    field: Mapping[str, Any],
    physical_factor: Any,
    *,
    z_path: str,
) -> dict[str, Any]:
    source_components = field.get("source_components", {})
    chi_source = source_components.get("axis_chi", field)
    w1_source = source_components.get("W1_F", field)
    info = _source_metadata(leading, chi_source, w1_source)
    coefficient = physical_factor * leading["coefficient"] * field["coefficient"]
    if info["error_only"]:
        coefficient = _symmetric(c, coefficient)
        if "symmetric_error" not in info["error_classes"]:
            info["error_classes"].append("symmetric_error")
    powers = [leading["powers"][index] + field["powers"][index] + PHYSICAL_POWER_SHIFT[index] for index in range(4)]
    ordinary_order = 0 if z_path == "leading*chi*W1_0" else 1
    origin = {
        "operation": "local_outer_D_physical_swirl_projection",
        "source_ids": list(info["source_ids"]),
        "source_paths": list(info["source_paths"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": info["uncertainty"],
        "error_only": info["error_only"],
        "symmetric_error": info["error_only"],
        "leibniz_piece": z_path,
        "ordinary_Z_order": ordinary_order,
        "source_ids_by_factor": info["source_ids_by_factor"],
        "source_paths_by_factor": info["source_paths_by_factor"],
        "error_classes_by_factor": info["error_classes_by_factor"],
        "formal_physical_power_added": list(PHYSICAL_POWER_SHIFT),
        "physical_measure_attached_once": True,
        "chi_reapplied": False,
        "green_reapplied": False,
        "W1_reapplied": False,
        "P_reapplied": False,
        "physical_P_compensation_included": False,
        "symbolic_cancellation_claimed": False,
    }
    return {
        "coefficient": coefficient,
        "coefficient_interval": coefficient,
        "powers": powers,
        "formal_powers": powers,
        "power_basis": list(SOURCE_POWER_BASIS),
        "physical_factor": physical_factor,
        "physical_power_shift": list(PHYSICAL_POWER_SHIFT),
        "field": "F1",
        "channel": "physical_swirl_m15",
        "source_kind": "raw_Q_local_outer_D_physical_swirl",
        "ordinary_Z_order": ordinary_order,
        "leibniz_piece": z_path,
        "source_ids": list(info["source_ids"]),
        "source_paths": list(info["source_paths"]),
        "error_classes": list(info["error_classes"]),
        "uncertainty": info["uncertainty"],
        "error_only": info["error_only"],
        "symmetric_error": info["error_only"],
        "source_ids_by_factor": info["source_ids_by_factor"],
        "source_paths_by_factor": info["source_paths_by_factor"],
        "error_classes_by_factor": info["error_classes_by_factor"],
        "leading_F": dict(leading),
        "D_F1": dict(field),
        "origin": origin,
    }


def _rows_view(rows: Mapping[int, Sequence[Mapping[str, Any]]]) -> dict[str, list[dict[str, Any]]]:
    return {f"Z{order}": list(rows[order]) for order in (0, 1)}


def kernel_mass(c: Any, sc: Any) -> Any:
    """Return the coefficient-only whole-phase K0 mass bound.

    The returned interval omits the formal ``hb**2`` and starts at zero.  It
    is an enclosure, not a point evaluation at q=1.
    """

    sc_i = _as_interval(c, sc, "sc")
    if _lower(sc_i) <= 0:
        raise ValueError("sc must be strictly positive")
    upper_exp = c.exp(sc_i / c.mpf(8))
    envelope = _make_interval(c, _zero(c), _singleton(c, upper_exp, upper=True))
    return (sc_i / c.mpf(4)) ** 2 / c.mpf(2) * envelope


def build(
    c: Any,
    sc: Any = None,
    low_word: Any = None,
    axis_chi: Any = None,
    leading_F: Any = None,
    *,
    chi: Any = None,
    chi_Z: Any = None,
    axis_chi_Z: Any = None,
    source_family: Any = None,
    source_atoms: Any = None,
    certificate: Any = None,
    geometry_certificate: Any = None,
    hb: Any = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build the conditional raw-Q local outer D physical F/m15 sectors.

    ``low_word`` is normally the complete caller result's
    ``local_signed_sectors`` mapping.  ``axis_chi`` contains Z0/Z1 chi
    sectors, and ``leading_F`` contains physical leading-F Z0/Z1 sectors.
    All three inputs retain four integer formal powers and source metadata.
    """

    if not all(hasattr(c, key) for key in ("mpf", "make_mpf", "exp")):
        raise TypeError("c must be a directed interval context")
    if sc is None:
        sc = aliases.pop("scale", aliases.pop("s_c", None))
    if low_word is None:
        low_word = aliases.pop("local_signed_sectors", aliases.pop("W1_F", aliases.pop("w1_F", None)))
    if axis_chi is None:
        axis_chi = aliases.pop("chi_jets", aliases.pop("axis", None))
    if leading_F is None:
        leading_F = aliases.pop("leading", aliases.pop("F_lead", aliases.pop("Flead", None)))
    if chi is not None:
        if axis_chi is not None:
            raise ValueError("axis_chi and chi were supplied twice")
        axis_chi = chi
    if axis_chi_Z is not None:
        if chi_Z is not None:
            raise ValueError("chi_Z and axis_chi_Z were supplied twice")
        chi_Z = axis_chi_Z
    if geometry_certificate is not None:
        if certificate is not None and geometry_certificate != certificate:
            raise ValueError("certificate and geometry_certificate disagree")
        certificate = geometry_certificate
    if aliases:
        names = ", ".join(sorted(str(name) for name in aliases))
        raise TypeError(f"unexpected local outer D arguments: {names}")
    if sc is None or low_word is None or axis_chi is None or leading_F is None:
        raise TypeError("sc, low_word, axis_chi, and leading_F are required")
    if source_family is None or (isinstance(source_family, str) and not source_family.strip()):
        raise ValueError("source_family is required")

    atoms = _source_atoms(source_atoms)
    sc_i, hb_i, geometry = _validate_certificate(c, sc, certificate, hb_override=hb)
    words = _rows(c, low_word, "low_word", field="F1")
    chi_rows = _chi_rows(c, axis_chi, chi_Z)
    leading_rows = _rows(c, leading_F, "leading_F", field="F1")

    mass = kernel_mass(c, sc_i)
    exp_sc8 = _make_interval(c, _one(c), _singleton(c, c.exp(sc_i / c.mpf(8)), upper=True))
    exp_sc4 = _make_interval(c, _one(c), _singleton(c, c.exp(sc_i / c.mpf(4)), upper=True))
    physical_factor = -c.mpf(5) * sc_i * exp_sc4

    field_rows: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    # D W1_F, Z0 = -2 K0 chi W1_0.
    for chi_term in chi_rows[0]:
        for source_term in words[0]:
            field_rows[0].append(_field_product(
                c, chi_term, source_term, mass,
                chi_order=0, source_order=0, z_path="chi*W1_0",
            ))
    # D W1_F, Z1 = -2 K0 (chi_Z W1_0 + chi W1_1), unmerged.
    for chi_term in chi_rows[1]:
        for source_term in words[0]:
            field_rows[1].append(_field_product(
                c, chi_term, source_term, mass,
                chi_order=1, source_order=0, z_path="chi_Z*W1_0",
            ))
    for chi_term in chi_rows[0]:
        for source_term in words[1]:
            field_rows[1].append(_field_product(
                c, chi_term, source_term, mass,
                chi_order=0, source_order=1, z_path="chi*W1_1",
            ))

    projected_rows: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for leading_term in leading_rows[0]:
        for field_term in field_rows[0]:
            projected_rows[0].append(_project_product(
                c, leading_term, field_term, physical_factor,
                z_path="leading*chi*W1_0",
            ))
    for leading_term in leading_rows[1]:
        for field_term in field_rows[0]:
            projected_rows[1].append(_project_product(
                c, leading_term, field_term, physical_factor,
                z_path="leading_Z*chi*W1_0",
            ))
    for leading_term in leading_rows[0]:
        for field_term in field_rows[1]:
            path = field_term["leibniz_piece"]
            projected_rows[1].append(_project_product(
                c, leading_term, field_term, physical_factor,
                z_path="leading*" + path,
            ))

    field_view = _rows_view(field_rows)
    projected_view = _rows_view(projected_rows)
    chi_view = _rows_view(chi_rows)
    words_view = _rows_view(words)
    leading_view = _rows_view(leading_rows)
    geometry.update({
        "x_formula": "2*exp(hb*sc*q/8)",
        "d_formula": "hb*sc*q/4",
        "q_zero": _zero(c),
        "q_one": _make_interval(c, 1, 1),
        "exp_sc_over_8": exp_sc8,
        "exp_sc_over_4": exp_sc4,
        "physical_factor": physical_factor,
        "physical_factor_formula": "-5*sc*[1,exp(sc/4)]",
        "physical_factor_negative": _upper(physical_factor) < 0,
        "formal_kernel_power_shift": list(KERNEL_POWER_SHIFT),
        "formal_physical_power_shift": list(PHYSICAL_POWER_SHIFT),
        "formal_total_power_shift": list(TOTAL_POWER_SHIFT),
    })
    kernel = {
        "operator": "D W1_F(q)=-2*integral_2^x(q) K0(x,t)*chi(t,Z)*W1_F(t,Z)dt",
        "K0": "K0(x,t)=t*(1-t^2/x^2)/2",
        "integral_K0": "cosh(d)-1",
        "d": "hb*sc*q/4",
        "coefficient_mass": mass,
        "mass_formula": "(sc/4)^2/2*[0,exp(sc/8)] with formal hb^2",
        "formal_hb_power": 2,
        "positive_whole_q_mass": True,
        "zero_at_q_zero": True,
        "upper_is_bound_not_value": True,
        "q_domain": geometry["q_domain"],
    }

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "source_family": source_family,
        "source_atoms": source_atoms if source_atoms is not None else dict(zip(SOURCE_POWER_BASIS, SOURCE_POWER_BASIS)),
        "source_power_basis": list(SOURCE_POWER_BASIS),
        "fields": list(FIELDS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": ["Z0", "Z1"],
        "geometry_certificate": certificate,
        "certificate": certificate,
        "geometry": geometry,
        "kernel": kernel,
        "axis_chi": chi_view,
        "chi": chi_view,
        "low_word": words_view,
        "W1_F": words_view,
        "leading_F": leading_view,
        "field_intermediate": {"F1": field_view, **field_view},
        "field": {"F1": field_view, **field_view},
        "D": {"F1": field_view, **field_view},
        "D_field": field_view,
        "D_F1": field_view,
        "projected": {"F1": projected_view, **projected_view},
        "swirl": projected_view,
        "physical_swirl": projected_view,
        "m15": projected_view,
        "F": projected_view,
        "projection": {
            "F1": projected_view,
            "physical_factor": physical_factor,
            "physical_factor_formula": "-5*sc*[1,exp(sc/4)]",
            "formal_power_shift": list(PHYSICAL_POWER_SHIFT),
            "ordinary_Z1_paths": list(Z1_PATHS),
            "measure_attached_once": True,
        },
        "kernel_mass": mass,
        "kernel_masses": {"F1": mass},
        "formal_total_power_shift": list(TOTAL_POWER_SHIFT),
        "ordinary_Z1_product_rule": "leading_Z*chi*W1_0 + leading*chi_Z*W1_0 + leading*chi*W1_1",
        "ordinary_Z1_paths": list(Z1_PATHS),
        "z1_terms_are_not_merged": True,
        "source_ids_preserved": True,
        "source_paths_preserved": True,
        "error_classes_preserved": True,
        "uncertainty_preserved": True,
        "error_only_inputs_retained_and_symmetric": True,
        "kernel_mass_upper_not_used_as_value": True,
        "q_zero_is_included": True,
        "same_source_certificate_consumed": True,
        "geometry_assumptions_caller_certified": True,
        "geometry_certificate_independently_verified": False,
        "axis_chi_is_not_phase_cutoff": True,
        "raw_Q": True,
        "raw_Q_partial": True,
        "raw_Q_components": ["C0", "B1"],
        "W1_F_saved_input": True,
        "W1_F_recomputed": False,
        "green_reapplied": False,
        "W1_reapplied": False,
        "P_reapplied": False,
        "physical_P_compensation_included": False,
        "other_inner_corrections_included": False,
        "higher_potential_tail_included": False,
        "higher_tail_included": False,
        "pressure_included": False,
        "P1_included": False,
        "only_F_output": True,
        "ledger_written": False,
        "no_ledger": True,
        "same_source_cancellation_claimed": False,
        "symbolic_cancellation_claimed": False,
        "complete_m15_computed": False,
        "complete_m15": False,
        "m15_computed": False,
        "finite_raw_Q_first_potential_only": True,
        "scope": "conditional raw-Q local outer D first-potential physical F swirl m15 partial only",
    }


build_local_outer_dq_m15 = build
evaluate = build
