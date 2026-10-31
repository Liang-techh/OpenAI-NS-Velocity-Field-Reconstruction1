"""Correlated local-``P``/``Dg`` physical-sector pair leaf.

This is a deliberately small, data-only bridge.  The caller supplies saved
physical ``Dg`` sectors for the finite ``J1`` branch and a certificate saying
that the two expressions have the same owner, ``g5``, ``K0`` and ``chi``
inputs, the ``+2`` pressure sign, the same physical measure/cutoff, and a
``Z``-independent ordinary derivative.  For every saved sector this leaf
returns two *correlated* records:

``Dg``
    the caller-owned record, copied without changing its interval, powers, or
    uncertainty metadata;
``P``
    a compensating record with the exact directed interval ``-Dg`` and the
    same formal powers and provenance.

The zero is a function-level identity for one pair, not an interval operation:
the leaf never adds the two boxed intervals.  The result therefore does not
claim that independently enclosed intervals sum to zero, and it does not
write a ledger or manufacture the missing core ``W0``, ``D P``, or ``W0`` tail
branches.

Importing this file performs no work and imports no producer or legacy module.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_local_p_dg_correlated_pair_v1"
COMPONENT = "local_p_dg_correlated_pair"
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")

_MISSING = object()
_ROW_KEYS = (
    ("Z0", "z0", "ordinary_Z0", "ordinary_z0", 0, "0"),
    ("Z1", "z1", "ordinary_Z1", "ordinary_z1", 1, "1"),
)
_WRAPPER_KEYS = (
    "physical_Dg",
    "Dg_physical",
    "dg_physical",
    "Dg",
    "dG",
    "J1",
    "j1",
    "ordinary_Z",
    "ordinary_z",
    "rows",
    "sectors",
    "terms",
)

__all__ = [
    "SCHEMA",
    "COMPONENT",
    "POWER_BASIS",
    "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS",
    "build",
    "build_local_p_dg_correlated_pair",
    "enclose",
    "evaluate",
]


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _lookup(mapping: Mapping[Any, Any], *keys: Any) -> Any:
    """Look up exact keys first, then their lower-case spellings."""

    for key in keys:
        try:
            if key in mapping:
                return mapping[key]
        except (TypeError, ValueError):
            continue
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for key in keys:
        value = lowered.get(str(key).lower(), _MISSING)
        if value is not _MISSING:
            return value
    return _MISSING


def _as_items(value: Any) -> list[Any]:
    if value is None or value is _MISSING:
        return []
    if isinstance(value, (str, bytes, bytearray)):
        return [value]
    if _is_sequence(value):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    """Preserve metadata order while removing only exact duplicate labels."""

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


def _raw_endpoint(value: Any, upper: bool) -> Any:
    """Return an endpoint accepted by ``make_mpf`` without midpointing."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    if hasattr(value, "a") and hasattr(value, "b"):
        endpoint = value.b if upper else value.a
        nested = getattr(endpoint, "_mpi_", None)
        if nested is not None:
            return nested[1 if upper else 0]
        return endpoint
    return value


def _point(c: Any, raw: Any) -> Any:
    """Create a singleton interval/scalar from a raw endpoint."""

    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((raw, raw))
        except Exception:
            pass
    try:
        return c.mpf(raw)
    except Exception:
        return raw


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    """Construct an interval from directed endpoints."""

    low_raw = _raw_endpoint(lower, False)
    high_raw = _raw_endpoint(upper, True)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((low_raw, high_raw))
        except Exception:
            pass
    try:
        return c.mpf([_point(c, low_raw), _point(c, high_raw)])
    except Exception:
        return c.mpf((low_raw, high_raw))


def _as_interval(c: Any, value: Any, label: str) -> Any:
    if hasattr(value, "_mpi_"):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return _interval(c, value, value)
    if _is_mapping(value):
        nested = _lookup(value, "interval", "enclosure", "box", "record", "value")
        if nested is not _MISSING and nested is not value:
            return _as_interval(c, nested, label)
        lower = _lookup(value, "lower", "lo", "left")
        upper = _lookup(value, "upper", "hi", "right")
        if lower is not _MISSING and upper is not _MISSING:
            return _interval(c, lower, upper)
        raise TypeError(f"{label} is not interval-compatible")
    if _is_sequence(value) and len(value) == 2:
        return _interval(c, value[0], value[1])
    try:
        return c.mpf(value)
    except Exception as exc:
        raise TypeError(f"{label} is not interval-compatible") from exc


def _bounds(value: Any) -> tuple[Any, Any]:
    if hasattr(value, "_mpi_"):
        return _point_from_raw_like(value, False), _point_from_raw_like(value, True)
    if hasattr(value, "a") and hasattr(value, "b"):
        return value.a, value.b
    raise TypeError("directed interval required")


def _point_from_raw_like(value: Any, upper: bool) -> Any:
    raw = _raw_endpoint(value, upper)
    # ``value.a``/``value.b`` are preferable for comparisons; use a tiny
    # local conversion only when an object exposes raw endpoints alone.
    endpoint = value.b if upper else value.a
    return endpoint if endpoint is not None else raw


def _negate_interval(c: Any, value: Any) -> Any:
    """Return the directed negative with endpoints reversed exactly.

    Negation is performed on the two represented endpoints, never on the
    interval as a midpoint or a sampled value.
    """

    lower, upper = _bounds(value)
    return _interval(c, -upper, -lower)


def _is_interval_like(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _normalise_labels(record: Mapping[Any, Any], label: str) -> dict[str, Any]:
    origin = record.get("origin") if _is_mapping(record.get("origin")) else {}

    source_ids = _lookup(record, "source_ids", "source_id", "sources", "source")
    if source_ids is _MISSING:
        source_ids = _lookup(origin, "source_ids", "source_id", "sources", "source")
    ids = _as_items(source_ids)
    if not ids:
        raise ValueError(f"{label} must retain at least one source_id")

    source_paths = _lookup(record, "source_paths", "source_path", "paths", "path")
    if source_paths is _MISSING:
        source_paths = _lookup(origin, "source_paths", "source_path", "paths", "path")
    paths = [] if source_paths is _MISSING else _as_items(source_paths)

    error_classes = _lookup(record, "error_classes", "error_class", "uncertainty_class")
    if error_classes is _MISSING:
        error_classes = _lookup(origin, "error_classes", "error_class", "uncertainty_class")
    classes = [] if error_classes is _MISSING else _as_items(error_classes)

    error_only = record.get("error_only", record.get("symmetric_error", False))
    if not isinstance(error_only, bool):
        raise TypeError(f"{label}.error_only must be boolean")
    symmetric_error = bool(record.get("symmetric_error", error_only))
    if not isinstance(symmetric_error, bool):
        raise TypeError(f"{label}.symmetric_error must be boolean")
    uncertainty = bool(
        record.get("uncertainty", record.get("is_uncertainty", origin.get("uncertainty", False)))
        or error_only
    )
    return {
        "source_ids": list(ids),
        "source_paths": list(paths),
        "error_classes": list(classes),
        "error_only": bool(error_only),
        "symmetric_error": symmetric_error,
        "uncertainty": uncertainty,
    }


def _normalise_powers(record: Mapping[Any, Any], label: str) -> list[int]:
    raw = _lookup(record, "powers", "formal_powers", "power_tuple")
    if raw is _MISSING:
        raw = [
            _lookup(record, "hb_power", "hb"),
            _lookup(record, "F0_power", "F0_anchor_power", "F0"),
            _lookup(record, "Pstar_squared_power", "Pstar_power", "Pstar2_power", "Pstar_squared"),
            _lookup(record, "Lambda_power", "lambda_power", "Lambda"),
        ]
    if _is_mapping(raw):
        raw = [
            _lookup(raw, "hb", 0),
            _lookup(raw, "F0_anchor", "F0", 1),
            _lookup(raw, "Pstar_squared", "Pstar^2", "Pstar2", 2),
            _lookup(raw, "Lambda", "lambda", 3),
        ]
    if not _is_sequence(raw) or len(raw) != 4:
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    if any(item is _MISSING or not isinstance(item, int) or isinstance(item, bool) for item in raw):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return [int(item) for item in raw]


def _parent_id(record: Mapping[Any, Any], labels: Mapping[str, Any], label: str) -> Any:
    parent = _lookup(
        record,
        "parent_id",
        "parent_ID",
        "parent",
        "term_id",
        "sector_id",
        "record_id",
        "id",
    )
    if parent is not _MISSING and parent is not None:
        return parent
    # A saved source id is a stable parent pointer when the source packet did
    # not carry a separate term id.  It is not invented independently of the
    # caller's provenance.
    if len(labels["source_ids"]) == 1:
        return labels["source_ids"][0]
    raise ValueError(f"{label} must retain parent_id or one unambiguous source_id")


def _normalise_sector(c: Any, record: Any, label: str) -> dict[str, Any]:
    if not _is_mapping(record):
        raise TypeError(f"{label} must be a sector mapping")
    coefficient = _lookup(record, "coefficient", "coefficient_interval", "interval", "enclosure", "value")
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    coefficient = _as_interval(c, coefficient, f"{label}.coefficient")
    if not _is_interval_like(coefficient):
        raise TypeError(f"{label}.coefficient must be directed")
    labels = _normalise_labels(record, label)
    powers = _normalise_powers(record, label)
    parent_id = _parent_id(record, labels, label)
    result = dict(record)
    result.update(
        {
            "coefficient": coefficient,
            "coefficient_interval": coefficient,
            "powers": list(powers),
            "formal_powers": list(powers),
            "power_basis": list(POWER_BASIS),
            "source_ids": list(labels["source_ids"]),
            "source_paths": list(labels["source_paths"]),
            "error_classes": list(labels["error_classes"]),
            "error_only": labels["error_only"],
            "symmetric_error": labels["symmetric_error"],
            "uncertainty": labels["uncertainty"],
            "parent_id": parent_id,
        }
    )
    return result


def _unwrap_rows(value: Any) -> Any:
    current = value
    # A sector itself must not be unwrapped through a field called ``value``.
    for _ in range(8):
        if not _is_mapping(current):
            break
        if _lookup(current, "coefficient", "coefficient_interval", "interval") is not _MISSING:
            break
        nested = _MISSING
        for key in _WRAPPER_KEYS:
            candidate = _lookup(current, key)
            if candidate is not _MISSING and candidate is not current:
                nested = candidate
                break
        if nested is _MISSING:
            break
        current = nested
    return current


def _row_value(current: Any, keys: Sequence[Any]) -> Any:
    if _is_mapping(current):
        return _lookup(current, *keys)
    return _MISSING


def _sector_list(value: Any, label: str) -> list[Any]:
    if _is_mapping(value):
        nested = _lookup(value, "sectors", "terms", "source_sectors", "records")
        if nested is not _MISSING and nested is not value:
            return _sector_list(nested, label)
        if _lookup(value, "coefficient", "coefficient_interval", "interval", "value") is not _MISSING:
            return [value]
        values = list(value.values())
    elif _is_sequence(value):
        values = list(value)
    else:
        raise TypeError(f"{label} must be a sector list")
    if not values or not all(_is_mapping(item) for item in values):
        raise TypeError(f"{label} must contain sector mappings")
    return values


def _normalise_rows(c: Any, value: Any, label: str) -> dict[str, list[dict[str, Any]]]:
    current = _unwrap_rows(value)
    if _is_mapping(current):
        # Permit one more field wrapper after an outer J1/Dg wrapper.
        nested_field = _lookup(current, "Dg", "dG", "physical_Dg", "Dg_physical")
        if nested_field is not _MISSING and nested_field is not current:
            current = _unwrap_rows(nested_field)
    if _is_sequence(current) and len(current) == 2:
        raw_rows = {"Z0": current[0], "Z1": current[1]}
    elif _is_mapping(current):
        raw_rows = {}
        for canonical, keys in zip(ORDINARY_Z_LABELS, _ROW_KEYS):
            selected = _row_value(current, keys)
            if selected is _MISSING:
                raise KeyError(f"{label} is missing {canonical}")
            raw_rows[canonical] = selected
    else:
        raise TypeError(f"{label} must provide ordinary Z0 and Z1 rows")
    return {
        canonical: [
            _normalise_sector(c, record, f"{label}.{canonical}[{index}]")
            for index, record in enumerate(_sector_list(raw_rows[canonical], f"{label}.{canonical}"))
        ]
        for canonical in ORDINARY_Z_LABELS
    }


def _walk_mappings(value: Any, seen: set[int] | None = None) -> list[Mapping[Any, Any]]:
    if seen is None:
        seen = set()
    if not _is_mapping(value) or id(value) in seen:
        return []
    seen.add(id(value))
    result: list[Mapping[Any, Any]] = [value]
    for child in value.values():
        if _is_mapping(child):
            result.extend(_walk_mappings(child, seen))
        elif _is_sequence(child):
            for item in child:
                if _is_mapping(item):
                    result.extend(_walk_mappings(item, seen))
    return result


def _named_values(certificate: Mapping[Any, Any], names: Sequence[str]) -> list[tuple[str, Any]]:
    wanted = {str(name).lower() for name in names}
    found: list[tuple[str, Any]] = []
    for mapping in _walk_mappings(certificate):
        for key, value in mapping.items():
            if str(key).lower() in wanted:
                found.append((str(key).lower(), value))
    return found


def _proof_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if _is_mapping(value):
        flags: list[bool] = []
        for key in ("passed", "proved", "verified", "certified", "same", "valid"):
            selected = _lookup(value, key)
            if selected is not _MISSING and isinstance(selected, bool):
                flags.append(selected)
        return bool(flags) and all(flags)
    return False


def _required_flag(
    certificate: Mapping[Any, Any],
    names: Sequence[str],
    *,
    label: str,
) -> bool:
    values = _named_values(certificate, names)
    if not values:
        raise ValueError(f"certificate is missing passed {label} proof")
    if any(not _proof_value(value) for _key, value in values):
        raise ValueError(f"certificate {label} proof is not passed")
    return True


def _same_id(left: Any, right: Any) -> bool:
    return left == right or str(left) == str(right)


def _owner_truth(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower().replace("_", " ").replace("-", " ") in {
            "dg",
            "d g",
            "dg physical",
            "physical dg",
            "saved dg physical",
            "local p dg",
            "local p dg physical",
        }
    if not _is_mapping(value):
        return False
    owner = _lookup(value, "owner", "owner_kind", "owner_role", "owned_by")
    if owner is not _MISSING:
        if isinstance(owner, bool):
            if not owner:
                return False
        elif not _owner_truth(owner):
            return False
    flags: list[bool] = []
    for key in ("owner_proved", "owner_verified", "proved", "verified", "passed", "valid"):
        selected = _lookup(value, key)
        if selected is not _MISSING:
            if not isinstance(selected, bool):
                return False
            flags.append(selected)
    # An explicit recognised owner is itself a proof; otherwise at least one
    # passed/verified flag is required.
    return bool(owner is not _MISSING or flags) and all(flags)


def _owner_entries(certificate: Mapping[Any, Any]) -> Any:
    names = (
        "owner_by_parent",
        "owners_by_parent",
        "owner_certificate",
        "owner_proof",
        "owners",
        "per_sector_owner",
        "owner",
    )
    for mapping in _walk_mappings(certificate):
        selected = _lookup(mapping, *names)
        if selected is not _MISSING and selected is not certificate:
            return selected
    return _MISSING


def _owner_for(parent_id: Any, raw: Any) -> Any:
    if raw is _MISSING:
        raise ValueError("certificate is missing per-sector owner proof")
    if _is_mapping(raw):
        # A sequence-like records wrapper is accepted without treating an
        # arbitrary global ``owner=True`` as proof for every sector.
        for key in ("records", "entries", "sectors", "parents"):
            nested = _lookup(raw, key)
            if nested is not _MISSING and nested is not raw:
                return _owner_for(parent_id, nested)
        direct = _lookup(raw, parent_id)
        if direct is _MISSING:
            direct = _lookup(raw, str(parent_id))
        if direct is not _MISSING:
            if not _owner_truth(direct):
                raise ValueError(f"owner proof for parent {parent_id!r} is not passed")
            return direct
        # A single record can carry its own parent id.
        candidate_parent = _lookup(raw, "parent_id", "parent", "term_id", "sector_id", "id")
        if candidate_parent is not _MISSING:
            if not _same_id(candidate_parent, parent_id) or not _owner_truth(raw):
                raise ValueError(f"owner proof for parent {parent_id!r} is not passed")
            return raw
        # Do not accept owner=True as a global proof for multiple sectors.
        raise ValueError(f"certificate has no owner entry for parent {parent_id!r}")
    if _is_sequence(raw):
        matches: list[Any] = []
        for item in raw:
            if not _is_mapping(item):
                continue
            candidate_parent = _lookup(item, "parent_id", "parent", "term_id", "sector_id", "id")
            if candidate_parent is not _MISSING and _same_id(candidate_parent, parent_id):
                matches.append(item)
        if len(matches) != 1 or not _owner_truth(matches[0]):
            raise ValueError(f"owner proof for parent {parent_id!r} is not passed")
        return matches[0]
    raise ValueError(f"certificate has no owner entry for parent {parent_id!r}")


def _pressure_sign(certificate: Mapping[Any, Any]) -> int:
    names = (
        "P_sign",
        "p_sign",
        "pressure_sign",
        "P_sign_value",
        "p_sign_value",
        "P_sign_proved",
        "p_sign_proved",
        "pressure_sign_proved",
    )
    values = _named_values(certificate, names)
    if not values:
        raise ValueError("certificate is missing passed +2 P-sign proof")
    sign_seen = False
    proof_seen = False
    for key, value in values:
        if isinstance(value, bool):
            if not value:
                raise ValueError("certificate P-sign proof is not passed")
            proof_seen = True
            # A direct boolean P_sign is a proof of the named +2 convention;
            # a separate *_value can still disambiguate it below.
            if key in {"p_sign", "pressure_sign"}:
                sign_seen = True
            continue
        if _is_mapping(value):
            sign = _lookup(value, "sign", "value", "coefficient")
            if sign is not _MISSING:
                if isinstance(sign, bool) or str(sign).strip() not in {"2", "+2"}:
                    try:
                        if sign != 2:
                            raise ValueError("wrong P sign")
                    except Exception as exc:
                        raise ValueError("certificate P sign is not +2") from exc
                sign_seen = True
            if not _proof_value(value):
                raise ValueError("certificate P-sign proof is not passed")
            proof_seen = True
            continue
        if isinstance(value, int) and not isinstance(value, bool):
            if value != 2:
                raise ValueError("certificate P sign is not +2")
            sign_seen = True
            proof_seen = True
            continue
        if isinstance(value, str) and value.strip() in {"2", "+2"}:
            sign_seen = True
            proof_seen = True
            continue
        raise ValueError("certificate P sign is not +2")
    if not (sign_seen and proof_seen):
        raise ValueError("certificate is missing passed +2 P-sign proof")
    return 2


def _validate_certificate(certificate: Any, parent_ids: Sequence[Any]) -> dict[str, Any]:
    if not _is_mapping(certificate):
        raise TypeError("caller certificate is required")
    declared_passed = _lookup(certificate, "passed", "checks_passed", "certificate_passed")
    if declared_passed is not _MISSING and declared_passed is not True:
        raise ValueError("caller certificate passed flag is not true")

    # The combined proof name is accepted only in addition to the individual
    # names below; each operator is still explicitly checked and a false
    # individual entry cannot be hidden by a combined flag.
    same_g5 = _required_flag(
        certificate,
        ("same_g5", "same_g5_source", "g5_same_source", "same_g5_proved", "same_g5_K0_chi"),
        label="same g5",
    )
    same_K0 = _required_flag(
        certificate,
        ("same_K0", "same_k0", "same_K0_source", "K0_same_source", "same_g5_K0_chi"),
        label="same K0",
    )
    same_chi = _required_flag(
        certificate,
        ("same_chi", "same_chi_source", "chi_same_source", "same_chi_proved", "same_g5_K0_chi"),
        label="same chi",
    )
    physical_measure = _required_flag(
        certificate,
        (
            "physical_measure_same",
            "same_physical_measure",
            "same_measure",
            "physical_measure_proved",
            "same_physical_measure_cutoff",
        ),
        label="same physical measure/cutoff",
    )
    z_independent = _required_flag(
        certificate,
        (
            "Z_independent_ordinary_derivative",
            "z_independent_ordinary_derivative",
            "ordinary_Z_derivative_Z_independent",
            "ordinary_derivative_Z_independent",
            "Z_independent",
        ),
        label="Z-independent ordinary derivative",
    )
    sign = _pressure_sign(certificate)

    raw_owner = _owner_entries(certificate)
    owners: dict[str, Any] = {}
    for parent_id in parent_ids:
        owner = _owner_for(parent_id, raw_owner)
        owners[str(parent_id)] = owner
    return {
        "passed": bool(same_g5 and same_K0 and same_chi and physical_measure and z_independent and sign == 2),
        "same_g5": same_g5,
        "same_K0": same_K0,
        "same_chi": same_chi,
        "P_sign": sign,
        "physical_measure_same": physical_measure,
        "Z_independent_ordinary_derivative": z_independent,
        "owner_by_parent": owners,
        "owner_entries_verified": len(owners),
    }


def _pair(
    c: Any,
    parent: Mapping[str, Any],
    *,
    z_label: str,
    z_order: int,
    index: int,
) -> dict[str, Any]:
    parent_id = parent["parent_id"]
    compensation = dict(parent)
    compensation["coefficient"] = _negate_interval(c, parent["coefficient"])
    compensation["coefficient_interval"] = compensation["coefficient"]
    compensation["powers"] = list(parent["powers"])
    compensation["formal_powers"] = list(parent["formal_powers"])
    compensation["source_ids"] = list(parent["source_ids"])
    compensation["source_paths"] = list(parent["source_paths"])
    compensation["error_classes"] = list(parent["error_classes"])
    compensation["error_only"] = parent["error_only"]
    compensation["symmetric_error"] = parent["symmetric_error"]
    compensation["uncertainty"] = parent["uncertainty"]
    compensation.update(
        {
            "role": "P_compensation",
            "correlated_with_parent_id": parent_id,
            "generated_from": "Dg",
            "interval_operation": "exact_directed_negative",
            "independent_interval_sum_performed": False,
        }
    )
    dg = dict(parent)
    dg.update(
        {
            "role": "Dg_parent",
            "correlated_parent_id": parent_id,
            "interval_operation": "caller_saved_interval_unchanged",
            "independent_interval_sum_performed": False,
        }
    )
    pair_id = f"{parent_id}::P_compensation"
    source_pointer = {
        "parent_id": parent_id,
        "source_ids": list(parent["source_ids"]),
        "source_paths": list(parent["source_paths"]),
        "error_classes": list(parent["error_classes"]),
        "uncertainty": parent["uncertainty"],
        "error_only": parent["error_only"],
    }
    return {
        "pair_id": pair_id,
        "index": index,
        "parent_id": parent_id,
        "ordinary_Z": z_label,
        "ordinary_Z_order": z_order,
        "Dg": dg,
        "P": compensation,
        "Dg_parent": dg,
        "P_compensation": compensation,
        "formal_powers": list(parent["formal_powers"]),
        "powers": list(parent["powers"]),
        "source_pointers": {
            "parent": dict(source_pointer),
            "Dg": dict(source_pointer),
            "P": dict(source_pointer),
        },
        "parent_source_ids": list(parent["source_ids"]),
        "parent_source_paths": list(parent["source_paths"]),
        "parent_error_classes": list(parent["error_classes"]),
        "error_only": parent["error_only"],
        "symmetric_error": parent["symmetric_error"],
        "uncertainty": parent["uncertainty"],
        "functional_identity": "G0 P W0_local = +2 K0(chi K0 g5) = -Dg",
        "functional_identity_scope": "same underlying local W0_F=K0 g5 function, correlated pair only",
        "physical_swirl_identity": "same leadingF * physical measure * cutoff preserves the sign reversal",
        "ordinary_derivative_scope": "ordinary Z0 and Z1; Z-independent derivative certificate consumed",
        "symbolic_sum": {
            "expression": f"Dg({parent_id!r}) + P_compensation({parent_id!r})",
            "exact": "0",
            "exact_zero": True,
            "scope": "correlated pair only",
            "residual_exact_zero": True,
            "not_unconditioned_global": True,
            "not_interval_overlap": True,
            "independent_boxed_interval_sum_performed": False,
        },
        "residual": {
            "value": "0",
            "kind": "correlated_symbolic_identity",
            "scope": "correlated pair only",
            "unconditioned_global": False,
        },
        "source_interval_overlap_used_as_identity": False,
        "errors_preserved_without_combination": True,
        "parent_record_copied": True,
    }


def build(
    c: Any,
    dg_physical: Any = None,
    certificate: Any = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build correlated ``Dg``/``P`` pairs from saved physical ``Dg`` rows.

    ``dg_physical`` must contain exactly ordinary ``Z0`` and ``Z1`` ``J1``
    rows (wrappers such as ``{"J1": ...}`` and ``{"Dg_physical": ...}``
    are accepted).  ``certificate`` is mandatory and is checked against every
    parent sector before any pair is returned.
    """

    if dg_physical is None:
        for name in (
            "Dg_physical",
            "physical_Dg",
            "saved_Dg_physical",
            "dg",
            "Dg",
            "input_sectors",
            "sectors",
        ):
            if name in aliases:
                dg_physical = aliases.pop(name)
                break
    if certificate is None:
        for name in ("source_certificate", "owner_certificate", "proof", "proof_certificate"):
            if name in aliases:
                certificate = aliases.pop(name)
                break
    if aliases:
        names = ", ".join(sorted(str(name) for name in aliases))
        raise TypeError(f"unexpected local P/Dg arguments: {names}")
    if dg_physical is None:
        raise TypeError("saved physical Dg sectors are required")
    if certificate is None:
        raise TypeError("caller proof certificate is required")
    if not hasattr(c, "mpf"):
        raise TypeError("c must provide directed interval arithmetic")

    rows = _normalise_rows(c, dg_physical, "dg_physical")
    parent_ids = [term["parent_id"] for z in ORDINARY_Z_LABELS for term in rows[z]]
    if len(parent_ids) != len({(type(item).__name__, str(item)) for item in parent_ids}):
        raise ValueError("saved Dg sectors must have unique parent_id values")
    certificate_record = _validate_certificate(certificate, parent_ids)

    pair_rows: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    dg_rows: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    p_rows: dict[str, list[dict[str, Any]]] = {"Z0": [], "Z1": []}
    for z_order, z_label in enumerate(ORDINARY_Z_LABELS):
        for index, parent in enumerate(rows[z_label]):
            pair = _pair(c, parent, z_label=z_label, z_order=z_order, index=index)
            pair_rows[z_label].append(pair)
            dg_rows[z_label].append(pair["Dg"])
            p_rows[z_label].append(pair["P"])

    # Every ``passed`` field below derives from validated caller evidence or
    # from the exact one-to-one construction above.  It is not an unconditional
    # assertion about a larger P branch.
    pair_count = sum(len(pair_rows[z]) for z in ORDINARY_Z_LABELS)
    result = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "passed": bool(certificate_record["passed"] and pair_count == len(parent_ids)),
        "checks_passed": bool(certificate_record["passed"] and pair_count == len(parent_ids)),
        "source_power_basis": list(POWER_BASIS),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "certificate": certificate_record,
        "certificate_validation": certificate_record,
        "caller_certificate_verified": bool(certificate_record["passed"]),
        "J1": pair_rows,
        "J1_pairs": pair_rows,
        "pairs": pair_rows,
        "correlated_pairs": pair_rows,
        "pair_rows": pair_rows,
        "pairs_by_family": {"J1": pair_rows},
        "Dg": dg_rows,
        "physical_Dg": dg_rows,
        "saved_Dg_physical": dg_rows,
        "P": p_rows,
        "P_compensation": p_rows,
        "physical_P_compensation": p_rows,
        "ordinary_Z": pair_rows,
        "source_pointers": {
            z: [pair["source_pointers"] for pair in pair_rows[z]] for z in ORDINARY_Z_LABELS
        },
        "functional_identity_scope": "same underlying local W0_F=K0 g5 function; correlated J1 pair only",
        "physical_swirl_identity_scope": "same leadingF, physical measure, and cutoff; physical m15 pair only",
        "ordinary_derivative_scope": "ordinary Z0 and Z1; no nonordinary derivative claim",
        "residual_exact_zero": True,
        "residual_scope": "correlated pair only",
        "unconditioned_global_zero": False,
        "symbolic_sum_zero_correlated_only": True,
        "interval_overlap_is_not_identity": True,
        "independent_boxed_interval_sums": False,
        "old_independent_errors_preserved": True,
        "error_only_preserved": True,
        "ledger_written": False,
        "writes_ledger": False,
        "no_ledger": True,
        "full_P_branch": False,
        "complete_P_branch": False,
        "full_branch_claimed": False,
        "core_W0_included": False,
        "D_P_included": False,
        "W0_tail_included": False,
        "missing_branches": ["core W0", "D P", "W0 tail"],
        "source_only_correlated_pair": True,
        "ancestor_replayed": False,
        "implementation_mirror_used": False,
        "pair_count": pair_count,
    }
    return result


def build_local_p_dg_correlated_pair(c: Any, dg_physical: Any = None, certificate: Any = None, **aliases: Any) -> dict[str, Any]:
    """Named alias for :func:`build`."""

    return build(c, dg_physical, certificate, **aliases)


enclose = build
evaluate = build

