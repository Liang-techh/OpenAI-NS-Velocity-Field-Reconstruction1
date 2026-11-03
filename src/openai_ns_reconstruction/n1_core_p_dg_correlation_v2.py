"""Conditional source-functional correlation for the core physical ``P``/``Dg`` pair.

This module is a data-only adapter.  It consumes a caller-owned physical
``Dg`` sector packet and caller-owned, already authenticated/decoded source
receipt evidence.  It does not open receipts, hash files, run a producer, or
perform an integral.  In particular, a shared source-family label or
overlapping intervals is not an owner proof.

When the exact source-functional contract is present, one serializable
object-level functional pair is emitted with separate Dg/P receipt
references.  It says that the two *whole functions* have opposite physical-
m15 sign; it deliberately does not materialize a per-sector P record or
add/subtract two interval boxes.
When any required proof is absent or false, the result is ``conditional``
with the saved Dg rows retained and no P alias.

The intended identity is conditional and local only::

    W0_F = K0 g5
    Dg   = -2 K0(chi K0 g5)
    P W0 = +2 K0(chi K0 g5) = -Dg

The ordinary Z1 proof must explicitly retain both chi_Z and g5_Z routes.  A
P1 collar object is not an m15 object, and no ledger, 212-term, or six-term
record is touched here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "n1_core_p_dg_correlation_v2"
COMPONENT = "n1_core_p_dg_correlation_v2"
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDINARY_Z_ORDERS = (0, 1)
ORDINARY_Z_LABELS = ("Z0", "Z1")

# These are the exact receipt inputs expected from the frozen consumer.  The
# leaf only checks the caller's evidence against these paths; it never opens
# or decompresses them.
RECEIPT_PATHS = {
    "KCHI": "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz",
    "FEEDBACK": "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz",
    "COMP": "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz",
    "LOW": "experiments/root_st073/n1_core_word1_z2_v2_receipt.json.gz",
    "LOW3": "experiments/root_st073/n1_core_low_word_v3_receipt.json.gz",
}
REQUIRED_RECEIPTS = tuple(RECEIPT_PATHS)

_MISSING = object()
_ROW_KEYS = (
    ("Z0", "z0", "ordinary_Z0", "ordinary_z0", 0, "0"),
    ("Z1", "z1", "ordinary_Z1", "ordinary_z1", 1, "1"),
)
_WRAPPER_KEYS = (
    "physical_Dg", "Dg_physical", "saved_Dg_physical", "Dg", "dG", "J1",
    "j1", "ordinary_Z", "ordinary_z", "rows", "sectors", "terms",
    "source_sectors",
)

__all__ = [
    "SCHEMA", "COMPONENT", "POWER_BASIS", "ORDINARY_Z_ORDERS",
    "ORDINARY_Z_LABELS", "RECEIPT_PATHS", "REQUIRED_RECEIPTS", "build",
    "build_core_p_dg_correlation", "enclose", "evaluate",
]


class _EvidenceFailure(ValueError):
    """A proof failure that must close the correlation without raising."""


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray))


def _lookup(mapping: Mapping[Any, Any], *keys: Any) -> Any:
    """Lookup exact keys, then case-insensitive spellings."""

    for key in keys:
        try:
            if key in mapping:
                return mapping[key]
        except (TypeError, ValueError):
            pass
    lowered = {str(key).lower(): value for key, value in mapping.items()}
    for key in keys:
        if str(key).lower() in lowered:
            return lowered[str(key).lower()]
    return _MISSING


def _as_items(value: Any) -> list[Any]:
    if value is _MISSING or value is None:
        return []
    if _is_sequence(value):
        return list(value)
    return [value]


def _endpoint(value: Any, upper: bool) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return raw[1 if upper else 0]
    if hasattr(value, "a") and hasattr(value, "b"):
        return value.b if upper else value.a
    return value


def _interval(c: Any, lower: Any, upper: Any) -> Any:
    """Build a directed interval without midpointing a source enclosure."""

    lo = _endpoint(lower, False)
    hi = _endpoint(upper, True)
    if hasattr(c, "make_mpf"):
        try:
            return c.make_mpf((lo, hi))
        except Exception:
            pass
    if hasattr(c, "mpf"):
        return c.mpf([lo, hi])
    raise TypeError("c must provide directed interval arithmetic")


def _as_interval(c: Any, value: Any, label: str) -> Any:
    if hasattr(value, "_mpi_"):
        return value
    if hasattr(value, "a") and hasattr(value, "b"):
        return _interval(c, value.a, value.b)
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
    if hasattr(c, "mpf"):
        return c.mpf(value)
    raise TypeError(f"{label} is not interval-compatible")


def _interval_like(value: Any) -> bool:
    return hasattr(value, "_mpi_") or (hasattr(value, "a") and hasattr(value, "b"))


def _normalise_labels(record: Mapping[Any, Any], label: str) -> dict[str, Any]:
    origin = record.get("origin") if _is_mapping(record.get("origin")) else {}
    source_ids = _lookup(record, "source_ids", "source_id", "sources", "source")
    if source_ids is _MISSING:
        source_ids = _lookup(origin, "source_ids", "source_id", "sources", "source")
    ids = _as_items(source_ids)
    if not ids:
        raise ValueError(f"{label} must retain at least one exact source id")
    source_paths = _lookup(record, "source_paths", "source_path", "paths", "path")
    if source_paths is _MISSING:
        source_paths = _lookup(origin, "source_paths", "source_path", "paths", "path")
    paths = _as_items(source_paths)
    if not paths:
        raise ValueError(f"{label} must retain an exact source id path")
    error_classes = _lookup(record, "error_classes", "error_class", "uncertainty_class")
    if error_classes is _MISSING:
        error_classes = _lookup(origin, "error_classes", "error_class", "uncertainty_class")
    classes = _as_items(error_classes)
    error_only = record.get("error_only", record.get("symmetric_error", False))
    if not isinstance(error_only, bool):
        raise TypeError(f"{label}.error_only must be boolean")
    symmetric_error = record.get("symmetric_error", error_only)
    if not isinstance(symmetric_error, bool):
        raise TypeError(f"{label}.symmetric_error must be boolean")
    uncertainty = bool(
        record.get("uncertainty", record.get("is_uncertainty", origin.get("uncertainty", False)))
        or error_only
    )
    return {
        "source_ids": list(ids), "source_paths": list(paths),
        "error_classes": list(classes), "error_only": error_only,
        "symmetric_error": symmetric_error, "uncertainty": uncertainty,
    }


def _normalise_powers(record: Mapping[Any, Any], label: str) -> list[int]:
    powers = _lookup(record, "powers", "formal_powers", "power_tuple")
    if powers is _MISSING:
        powers = [
            _lookup(record, "hb_power", "hb"),
            _lookup(record, "F0_power", "F0_anchor_power", "F0"),
            _lookup(record, "Pstar_squared_power", "Pstar_power", "Pstar2_power", "Pstar_squared"),
            _lookup(record, "Lambda_power", "lambda_power", "Lambda"),
        ]
    if _is_mapping(powers):
        powers = [
            _lookup(powers, "hb", 0), _lookup(powers, "F0_anchor", "F0", 1),
            _lookup(powers, "Pstar_squared", "Pstar^2", "Pstar2", 2),
            _lookup(powers, "Lambda", "lambda", 3),
        ]
    if not _is_sequence(powers) or len(powers) != 4:
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    if any(item is _MISSING or isinstance(item, bool) or not isinstance(item, int) for item in powers):
        raise TypeError(f"{label}.powers must contain four integer formal powers")
    return [int(item) for item in powers]


def _parent_id(record: Mapping[Any, Any], labels: Mapping[str, Any], label: str) -> Any:
    parent = _lookup(record, "parent_id", "parent_ID", "parent", "term_id", "sector_id", "record_id", "id")
    if parent is not _MISSING and parent is not None:
        return parent
    if len(labels["source_ids"]) == 1:
        return labels["source_ids"][0]
    raise ValueError(f"{label} must retain parent_id or one unambiguous source id")


def _normalise_sector(c: Any, record: Any, label: str) -> dict[str, Any]:
    if not _is_mapping(record):
        raise TypeError(f"{label} must be a sector mapping")
    coefficient = _lookup(record, "coefficient", "coefficient_interval", "interval", "enclosure", "value")
    if coefficient is _MISSING:
        raise KeyError(f"{label} is missing coefficient")
    coefficient = _as_interval(c, coefficient, f"{label}.coefficient")
    if not _interval_like(coefficient):
        raise TypeError(f"{label}.coefficient must be directed")
    labels = _normalise_labels(record, label)
    powers = _normalise_powers(record, label)
    parent_id = _parent_id(record, labels, label)
    result = dict(record)
    result.update({
        "coefficient": coefficient, "coefficient_interval": coefficient,
        "powers": list(powers), "formal_powers": list(powers),
        "power_basis": list(POWER_BASIS), "source_ids": labels["source_ids"],
        "source_paths": labels["source_paths"], "error_classes": labels["error_classes"],
        "error_only": labels["error_only"], "symmetric_error": labels["symmetric_error"],
        "uncertainty": labels["uncertainty"], "parent_id": parent_id,
    })
    return result


def _unwrap_rows(value: Any) -> Any:
    current = value
    for _ in range(10):
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


def _normalise_rows(c: Any, value: Any) -> dict[str, list[dict[str, Any]]]:
    current = _unwrap_rows(value)
    if _is_mapping(current):
        nested = _lookup(current, "Dg", "dG", "physical_Dg", "Dg_physical")
        if nested is not _MISSING and nested is not current:
            current = _unwrap_rows(nested)
    if _is_sequence(current) and len(current) == 2:
        raw_rows = {"Z0": current[0], "Z1": current[1]}
    elif _is_mapping(current):
        raw_rows = {}
        for canonical, keys in zip(ORDINARY_Z_LABELS, _ROW_KEYS):
            selected = _lookup(current, *keys)
            if selected is _MISSING:
                raise KeyError(f"physical Dg is missing {canonical}")
            raw_rows[canonical] = selected
    else:
        raise TypeError("physical Dg must provide ordinary Z0 and Z1 rows")
    return {
        canonical: [
            _normalise_sector(c, row, f"physical_Dg.{canonical}[{index}]")
            for index, row in enumerate(_sector_list(raw_rows[canonical], f"physical_Dg.{canonical}"))
        ]
        for canonical in ORDINARY_Z_LABELS
    }


def _walk(value: Any, seen: set[int] | None = None) -> list[Mapping[Any, Any]]:
    if seen is None:
        seen = set()
    if not _is_mapping(value) or id(value) in seen:
        return []
    seen.add(id(value))
    result = [value]
    for child in value.values():
        if _is_mapping(child):
            result.extend(_walk(child, seen))
        elif _is_sequence(child):
            for item in child:
                if _is_mapping(item):
                    result.extend(_walk(item, seen))
    return result


def _named(certificate: Mapping[Any, Any], names: Sequence[str]) -> list[tuple[str, Any]]:
    wanted = {str(name).lower() for name in names}
    found: list[tuple[str, Any]] = []
    for mapping in _walk(certificate):
        for key, value in mapping.items():
            if str(key).lower() in wanted:
                found.append((str(key).lower(), value))
    return found


def _proof(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if _is_mapping(value):
        flags = []
        for key in ("passed", "proved", "verified", "certified", "same", "valid", "available"):
            selected = _lookup(value, key)
            if selected is not _MISSING:
                if not isinstance(selected, bool):
                    return False
                flags.append(selected)
        return bool(flags) and all(flags)
    return False


def _require_true(certificate: Mapping[Any, Any], names: Sequence[str], label: str) -> bool:
    values = _named(certificate, names)
    if not values:
        raise _EvidenceFailure(f"missing explicit proof: {label}")
    if any(not _proof(value) for _key, value in values):
        raise _EvidenceFailure(f"proof is false or malformed: {label}")
    return True


def _sign(certificate: Mapping[Any, Any]) -> int:
    values = _named(certificate, ("P_sign", "p_sign", "pressure_sign", "P_sign_value", "p_sign_value"))
    if not values:
        raise _EvidenceFailure("missing explicit +2 P-sign proof")
    seen = False
    for key, value in values:
        if isinstance(value, bool):
            if not value:
                raise _EvidenceFailure("P-sign proof is false")
            continue
        if _is_mapping(value):
            sign = _lookup(value, "sign", "value", "coefficient")
            if sign is _MISSING or str(sign).strip() not in {"2", "+2"}:
                raise _EvidenceFailure("P sign is not +2")
            if not _proof(value):
                raise _EvidenceFailure("P-sign proof is not passed")
            seen = True
            continue
        if isinstance(value, int) and not isinstance(value, bool):
            if value != 2:
                raise _EvidenceFailure("P sign is not +2")
            seen = True
            continue
        if isinstance(value, str) and value.strip() in {"2", "+2"}:
            seen = True
            continue
        raise _EvidenceFailure("P sign is not +2")
    if not seen:
        raise _EvidenceFailure("missing explicit +2 P-sign value")
    return 2


def _receipt_map(certificate: Mapping[Any, Any], source_receipts: Any) -> Mapping[Any, Any]:
    value = source_receipts
    if value is None:
        for name in ("source_receipts", "receipt_evidence", "source_evidence", "authenticated_sources", "decoded_receipts"):
            selected = _lookup(certificate, name)
            if selected is not _MISSING:
                value = selected
                break
    if not _is_mapping(value):
        raise _EvidenceFailure("caller must provide authenticated decoded source receipts")
    return value


def _receipt_entry(receipts: Mapping[Any, Any], label: str) -> Mapping[Any, Any]:
    expected = RECEIPT_PATHS[label]
    selected = _lookup(receipts, label, label.lower(), expected)
    if selected is _MISSING:
        # A path-keyed packet may use a nested ``receipts`` field.
        nested = _lookup(receipts, "receipts", "entries", "by_component")
        if _is_mapping(nested):
            selected = _lookup(nested, label, label.lower(), expected)
    if not _is_mapping(selected):
        raise _EvidenceFailure(f"missing decoded receipt evidence: {label}")
    path = _lookup(selected, "path", "receipt_path", "source_path", "exact_path")
    if path is not _MISSING and path != expected:
        raise _EvidenceFailure(f"receipt path mismatch for {label}")
    return selected


def _hash_ok(entry: Mapping[Any, Any], label: str) -> None:
    values = _named(entry, ("raw_hash_authenticated", "raw_hash_verified", "hash_authenticated", "hash_verified", "authenticated"))
    if not values or any(value is not True for _key, value in values):
        raise _EvidenceFailure(f"raw hash is not authenticated for {label}")
    decoded = _named(entry, ("decoded", "decoded_receipt", "decoded_ok", "receipt_decoded"))
    if not decoded or any(value is not True for _key, value in decoded):
        raise _EvidenceFailure(f"receipt is not marked decoded for {label}")
    digest = _lookup(entry, "raw_sha256", "sha256", "raw_hash", "receipt_sha256")
    if digest is _MISSING or not isinstance(digest, str) or not digest:
        raise _EvidenceFailure(f"raw hash value is missing for {label}")
    manifest = _lookup(entry, "source_manifest", "manifest", "source", "source_manifest_record")
    if not _is_mapping(manifest):
        raise _EvidenceFailure(f"source manifest is missing for {label}")
    manifest_hash = _lookup(manifest, "sha256", "manifest_sha256", "hash", "source_manifest_sha256")
    files = _lookup(manifest, "files", "source_files", "entries")
    if manifest_hash is _MISSING and not (_is_sequence(files) and len(files) > 0):
        raise _EvidenceFailure(f"source manifest hash/files are missing for {label}")
    output = _lookup(entry, "decoded_output", "output", "result_output", "result")
    if output is _MISSING or output is None:
        raise _EvidenceFailure(f"decoded output is missing for {label}")
    status = _lookup(entry, "status", "result_status")
    if status is not _MISSING and status not in ("returned", "verified", "ok"):
        raise _EvidenceFailure(f"receipt status is not returned for {label}")


def _owner_entries(certificate: Mapping[Any, Any]) -> Any:
    for name in ("owner_by_parent", "owners_by_parent", "per_sector_owner", "owner_certificate", "owner_proof"):
        value = _lookup(certificate, name)
        if value is not _MISSING:
            return value
    return _MISSING


def _owner_ok(value: Any) -> bool:
    if not _is_mapping(value):
        return False
    owner = _lookup(value, "owner", "owner_kind", "owner_role", "owned_by")
    if owner is _MISSING:
        return False
    owner_text = str(owner).lower().replace("_", " ").replace("-", " ")
    if owner_text not in {"dg physical", "physical dg", "kchi dg physical", "saved dg physical"}:
        return False
    passed = _lookup(value, "passed", "proved", "verified", "owner_verified")
    return passed is True


def _owner_for(raw: Any, parent_id: Any) -> Mapping[Any, Any]:
    if _is_mapping(raw):
        direct = _lookup(raw, parent_id, str(parent_id))
        if direct is _MISSING:
            nested = _lookup(raw, "records", "entries", "sectors", "parents")
            if nested is not _MISSING:
                return _owner_for(nested, parent_id)
            candidate = _lookup(raw, "parent_id", "parent", "term_id", "sector_id", "id")
            if candidate is not _MISSING and str(candidate) == str(parent_id):
                direct = raw
        if not _is_mapping(direct) or not _owner_ok(direct):
            raise _EvidenceFailure(f"owner proof is missing for parent {parent_id!r}")
        return direct
    if _is_sequence(raw):
        matches = []
        for item in raw:
            if _is_mapping(item):
                candidate = _lookup(item, "parent_id", "parent", "term_id", "sector_id", "id")
                if candidate is not _MISSING and str(candidate) == str(parent_id):
                    matches.append(item)
        if len(matches) != 1 or not _owner_ok(matches[0]):
            raise _EvidenceFailure(f"owner proof is missing for parent {parent_id!r}")
        return matches[0]
    raise _EvidenceFailure(f"owner proof is missing for parent {parent_id!r}")


def _owner_path_checks(owner: Mapping[Any, Any], sector: Mapping[str, Any], parent_id: Any) -> None:
    parent = _lookup(owner, "parent_id", "parent", "term_id", "sector_id", "id")
    if parent is not _MISSING and str(parent) != str(parent_id):
        raise _EvidenceFailure(f"owner parent mismatch for {parent_id!r}")
    ids = _as_items(_lookup(owner, "source_ids", "source_id", "exact_source_ids"))
    if ids and ids != list(sector["source_ids"]):
        raise _EvidenceFailure(f"owner source ids mismatch for {parent_id!r}")
    paths = _as_items(_lookup(owner, "source_paths", "source_path", "exact_source_id_path", "exact_source_paths"))
    if not paths:
        raise _EvidenceFailure(f"exact source id path is missing for {parent_id!r}")
    if paths != list(sector["source_paths"]):
        raise _EvidenceFailure(f"owner source path mismatch for {parent_id!r}")
    for names, label in (
        (("parent_hash", "same_parent_hash", "parent_source_hash"), "parent hash"),
        (("geometry", "geometry_label", "same_geometry"), "geometry"),
        (("labels", "source_labels", "exact_labels"), "labels"),
        (("exact_cell_path", "cell_path", "source_cell_path"), "cell path"),
        (("output_path", "exact_output_path", "output"), "output path"),
    ):
        value = _lookup(owner, *names)
        if value is _MISSING or value is None or value == "":
            raise _EvidenceFailure(f"{label} is missing for {parent_id!r}")


def _validate_evidence(
    certificate: Any,
    source_receipts: Any,
    rows: Mapping[str, list[Mapping[str, Any]]],
) -> dict[str, Any]:
    if not _is_mapping(certificate):
        raise _EvidenceFailure("caller correlation certificate is required")
    declared = _lookup(certificate, "passed", "checks_passed", "certificate_passed")
    if declared is not _MISSING and declared is not True:
        raise _EvidenceFailure("caller certificate is not marked passed")

    receipts = _receipt_map(certificate, source_receipts)
    receipt_meta = {}
    for label in REQUIRED_RECEIPTS:
        entry = _receipt_entry(receipts, label)
        _hash_ok(entry, label)
        receipt_meta[label] = {
            "path": RECEIPT_PATHS[label],
            "raw_hash": _lookup(entry, "raw_sha256", "sha256", "raw_hash", "receipt_sha256"),
            "manifest": _lookup(entry, "source_manifest", "manifest", "source", "source_manifest_record"),
        }

    # Every one of these is an explicit functional/owner assertion.  Plain
    # same-family labels are intentionally not accepted as substitutes.
    checks = {
        "same_actual_g5": _require_true(certificate, ("same_actual_g5", "same_g5_function", "same_g5_source_function", "g5_exact_owner_verified"), "same actual g5 function"),
        "same_chi_Z": _require_true(certificate, ("same_chi_Z", "same_chi_function", "same_chi_source_function", "chi_Z_same_source"), "same chi(Z) function"),
        "same_green_K0": _require_true(certificate, ("same_green_K0", "same_K0_function", "K0_exact_same_source", "same_zero_inlet_green"), "same Green K0"),
        "same_complete_kernel": _require_true(certificate, ("same_complete_0_2_kernel", "same_complete_[0,2]_kernel", "complete_kernel_0_2_same", "same_complete_kernel"), "same complete [0,2] kernel"),
        "same_physical_m15_measure": _require_true(certificate, ("same_physical_m15_measure", "same_m15_measure", "physical_m15_measure_same", "same_physical_measure_cutoff"), "same physical m15 measure"),
        "same_leading_F": _require_true(certificate, ("same_leading_F", "same_leading_F0", "leading_F_same_source", "same_leading_field"), "same leading F"),
        "low3_W0_same_K0_g5": _require_true(certificate, ("low3_W0_same_K0_g5", "low3_W0_exact_K0_g5", "W0_same_K0_g5", "low3_exact_function_owner"), "LOW3 W0 = same K0*g5"),
        "kchi_Dg_identity": _require_true(certificate, ("kchi_Dg_minus2_chi_K0_K0g5", "kchi_Dg_same_minus2_chi_K0_K0g5", "Dg_same_minus2_chi_K0_K0g5", "kchi_Dg_function_identity"), "KCHI Dg function identity"),
        "ordinary_Z1_chi_Z": _require_true(certificate, ("ordinary_Z1_chi_Z_retained", "ordinary_Z1_retains_chi_Z", "Z1_chi_Z_present"), "ordinary Z1 chi_Z route"),
        "ordinary_Z1_g5_Z": _require_true(certificate, ("ordinary_Z1_g5_Z_retained", "ordinary_Z1_retains_g5_Z", "Z1_g5_Z_present"), "ordinary Z1 g5_Z route"),
        "same_geometry": _require_true(certificate, ("same_geometry", "geometry_same_source", "physical_geometry_same"), "same geometry"),
        "same_labels": _require_true(certificate, ("same_labels", "source_labels_same", "physical_labels_same"), "same labels"),
        "same_parent_hashes": _require_true(certificate, ("same_parent_hashes", "parent_hashes_same", "all_parent_hashes_match"), "same parent hashes"),
        "exact_source_id_paths": _require_true(certificate, ("exact_source_id_paths", "source_id_paths_verified", "source_paths_exact"), "exact source id paths"),
        "exact_cell_paths": _require_true(certificate, ("exact_cell_paths", "cell_paths_verified", "source_cell_paths_exact"), "exact cell paths"),
        "exact_outputs": _require_true(certificate, ("exact_outputs", "outputs_verified", "source_outputs_exact"), "exact output paths/records"),
    }
    sign = _sign(certificate)

    raw_owner = _owner_entries(certificate)
    if raw_owner is _MISSING:
        raise _EvidenceFailure("per-sector Dg owner proof is required")
    owners: dict[str, Any] = {}
    for z_label in ORDINARY_Z_LABELS:
        for sector in rows[z_label]:
            parent_id = sector["parent_id"]
            owner = _owner_for(raw_owner, parent_id)
            _owner_path_checks(owner, sector, parent_id)
            owners[str(parent_id)] = dict(owner)
    dg_refs = _lookup(certificate, "Dg_source_refs", "dg_source_refs", "saved_Dg_source_refs")
    p_refs = _lookup(certificate, "P_source_refs", "p_source_refs", "saved_P_source_refs")
    if dg_refs is _MISSING or not dg_refs:
        raise _EvidenceFailure("saved Dg source references are required")
    if p_refs is _MISSING or not p_refs:
        raise _EvidenceFailure("independent P source references are required")
    return {
        # This validator establishes only that the caller supplied the
        # requested assertions and provenance records. It does not prove the
        # underlying live source functions.
        "passed": False, "caller_assertions_satisfied": True,
        "conditional_on_caller_assertions": True, "caller_certificate_schema_satisfied": True,
        "P_sign": sign, "checks": checks,
        "owner_by_parent": owners, "owner_entries_asserted": len(owners),
        "source_receipts": receipt_meta, "Dg_source_refs": dg_refs, "P_source_refs": p_refs,
    }


def _functional_pair(validation: Mapping[str, Any], rows: Mapping[str, list[Mapping[str, Any]]]) -> dict[str, Any]:
    """Create one object-level signed functional pair, never per-sector P rows."""

    dg_refs = validation["Dg_source_refs"]
    p_refs = validation["P_source_refs"]
    counts = {z: len(rows[z]) for z in ORDINARY_Z_LABELS}
    return {
        "schema": SCHEMA, "pair_id": "core_physical_m15_whole_function",
        "scope": "whole function physical m15 only", "Dg_source_refs": dg_refs,
        "P_source_refs": p_refs, "Dg_sector_counts": counts,
        "functional_expression": {
            "W0_F": "K0(g5)", "Dg_F": "-2 K0(chi K0(g5))",
            "P_W0_F": "+2 K0(chi K0(g5))",
            "axis_chi_constant_form": "Dg_F=-2 chi K0(K0(g5)); P_W0_F=+2 chi K0(K0(g5))",
            "ordinary_Z1": "-2 K0(chi_Z K0(g5) + chi K0(g5_Z))",
            "signed_relation": "P_W0_F = -Dg_F",
        },
        "sign_relation": -1, "symbolic_identity_under_assumptions": True,
        "symbolic_implication_under_assumptions": True,
        "actual_source_verified": False, "identity_verified": False,
        "zero_claimed": False,
        "coefficient_arithmetic_performed": False,
        "independent_boxed_interval_sum_performed": False,
        "linewise_cancellation_performed": False,
        "source_interval_overlap_used_as_identity": False,
        "saved_Dg_rows_retained": True, "independent_P_sector_materialized": False,
        "P_alias_created": False, "whole_function_sign_only": True,
        "old_errors_preserved": True,
    }


def _conditional_result(rows: Mapping[str, list[dict[str, Any]]], reasons: list[str], *, certificate: Any = None) -> dict[str, Any]:
    dg_rows = {z: [dict(item) for item in rows[z]] for z in ORDINARY_Z_LABELS}
    return {
        "schema": SCHEMA, "component": COMPONENT, "status": "conditional", "passed": False,
        "checks_passed": False, "correlation_proved": False,
        "actual_source_verified": False, "identity_verified": False,
        "zero_claimed": False, "residual_exact_zero": False,
        "conditional_reasons": list(reasons), "certificate": certificate,
        "certificate_validation": {"passed": False, "reasons": list(reasons)},
        "caller_certificate_verified": False, "source_power_basis": list(POWER_BASIS),
        "functional_expression": {
            "W0_F": "K0(g5)",
            "Dg_F": "-2 K0(chi K0(g5))",
            "P_W0_F": "+2 K0(chi K0(g5))",
            "signed_relation": "P_W0_F = -Dg_F only if exact common-function proof passes",
            "ordinary_Z1": "-2 K0(chi_Z K0(g5) + chi K0(g5_Z))",
            "common_g5_retained": True,
            "identity_verified": False,
            "zero_claimed": False,
        },
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
        "Dg": dg_rows, "physical_Dg": dg_rows, "saved_Dg_physical": dg_rows,
        "pairs": {"Z0": [], "Z1": []}, "correlated_pairs": {"Z0": [], "Z1": []},
        "functional_pairs": [],
        "P": {"Z0": [], "Z1": []}, "P_compensation": {"Z0": [], "Z1": []},
        "physical_P_compensation": {"Z0": [], "Z1": []},
        "source_only_correlated_pair": False, "unproved_source_function": True,
        "unconditioned_global_zero": False, "residual_scope": "none without proof",
        "interval_overlap_is_not_identity": True, "independent_boxed_interval_sums": False,
        "old_independent_errors_preserved": True, "error_only_preserved": True,
        "ledger_written": False, "writes_ledger": False, "no_ledger": True,
        "core_W0_included": False, "D_P_included": False, "W0_tail_included": False,
        "P1_is_m15": False, "P1_not_m15": True, "physical_m15_computed": False,
        "coefficient_arithmetic_performed": False, "linewise_cancellation_performed": False,
        "receipt_paths_only": True, "pair_count": 0,
    }


def build(
    c: Any,
    dg_physical: Any = None,
    certificate: Any = None,
    source_receipts: Any = None,
    **aliases: Any,
) -> dict[str, Any]:
    """Build a conditional physical source-functional correlation.

    ``dg_physical`` is the saved physical m15 Dg packet with ordinary Z0/Z1
    rows.  ``certificate`` and ``source_receipts`` are caller-owned evidence;
    they must carry explicit exact-owner/source-manifest/path/output proofs.
    Missing or false evidence returns a conditional result with no signed
    functional pair.
    """

    if dg_physical is None:
        for name in ("Dg_physical", "physical_Dg", "saved_Dg_physical", "Dg", "input_sectors", "sectors"):
            if name in aliases:
                dg_physical = aliases.pop(name)
                break
    if certificate is None:
        for name in ("source_certificate", "proof", "proof_certificate", "correlation_certificate", "evidence"):
            if name in aliases:
                certificate = aliases.pop(name)
                break
    if source_receipts is None:
        for name in ("receipt_evidence", "source_evidence", "authenticated_sources", "decoded_receipts"):
            if name in aliases:
                source_receipts = aliases.pop(name)
                break
    if aliases:
        names = ", ".join(sorted(str(name) for name in aliases))
        raise TypeError(f"unexpected correlation arguments: {names}")
    if dg_physical is None:
        raise TypeError("saved physical Dg sectors are required")
    if not hasattr(c, "mpf"):
        raise TypeError("c must provide directed interval arithmetic")

    rows = _normalise_rows(c, dg_physical)
    parent_ids = [row["parent_id"] for z in ORDINARY_Z_LABELS for row in rows[z]]
    if len(parent_ids) != len({(type(item).__name__, str(item)) for item in parent_ids}):
        raise ValueError("saved Dg sectors must have unique parent_id values")
    try:
        validation = _validate_evidence(certificate, source_receipts, rows)
    except _EvidenceFailure as exc:
        return _conditional_result(rows, [str(exc)], certificate=certificate)

    dg_rows: dict[str, list[dict[str, Any]]] = {
        z: [dict(item) for item in rows[z]] for z in ORDINARY_Z_LABELS
    }
    functional_pair = _functional_pair(validation, rows)
    pair_rows = {"whole_function": [functional_pair], "Z0": [], "Z1": []}
    p_rows = {"Z0": [], "Z1": []}
    return {
        "schema": SCHEMA, "component": COMPONENT,
        "status": "conditional_on_caller_assertions", "passed": False,
        "checks_passed": False, "conditional_contract_satisfied": True,
        "actual_source_verified": False, "correlation_proved": False,
        "source_power_basis": list(POWER_BASIS), "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_labels": list(ORDINARY_Z_LABELS), "certificate": validation,
        "certificate_validation": validation, "caller_certificate_verified": False,
        "caller_certificate_schema_satisfied": True, "caller_assertions_satisfied": True,
        "conditional_on_caller_assertions": True,
        "identity_verified": False, "zero_claimed": False,
        "functional_expression": {
            "W0_F": "K0(g5)",
            "Dg_F": "-2 K0(chi K0(g5))",
            "P_W0_F": "+2 K0(chi K0(g5))",
            "signed_relation": "P_W0_F = -Dg_F",
            "common_g5_retained": True,
            "identity_verified": False,
            "zero_claimed": False,
            "symbolic_implication_under_assumptions": True,
        },
        "J1": {"Z0": [], "Z1": []}, "J1_pairs": pair_rows, "pairs": pair_rows,
        "correlated_pairs": pair_rows, "functional_pairs": [functional_pair],
        "pair_rows": pair_rows, "pairs_by_family": {"whole_function": [functional_pair]},
        "Dg": dg_rows, "physical_Dg": dg_rows, "saved_Dg_physical": dg_rows,
        "P": p_rows, "P_compensation": p_rows, "physical_P_compensation": p_rows,
        "ordinary_Z": {"Z0": [], "Z1": []},
        "source_pointers": {"Dg": validation["Dg_source_refs"], "P": validation["P_source_refs"]},
        "functional_identity_scope": "conditional same actual g5, chi(Z), K0, complete [0,2] kernel; physical m15 pair only",
        "physical_swirl_identity_scope": "same leading F and physical m15 measure; whole-function sign alias only",
        "ordinary_derivative_scope": "ordinary Z0/Z1 with chi_Z and g5_Z retained",
        "residual_exact_zero": False,
        "residual_scope": "conditional symbolic implication under caller assertions",
        "unconditioned_global_zero": False, "symbolic_sum_zero_correlated_only": False,
        "interval_overlap_is_not_identity": True, "independent_boxed_interval_sums": False,
        "coefficient_arithmetic_performed": False, "linewise_cancellation_performed": False,
        "old_independent_errors_preserved": True, "error_only_preserved": True,
        "ledger_written": False, "writes_ledger": False, "no_ledger": True,
        "full_P_branch": False, "complete_P_branch": False, "full_branch_claimed": False,
        "core_W0_included": False, "D_P_included": False, "W0_tail_included": False,
        "missing_branches": ["core W0", "D P", "W0 tail"], "source_only_correlated_pair": True,
        "ancestor_replayed": False, "implementation_mirror_used": False,
        "P1_is_m15": False, "P1_not_m15": True, "physical_m15_computed": False,
        "receipt_paths_only": True, "pair_count": 0, "functional_pair_count": 1,
        "independent_P_sector_materialized": False,
    }


def build_core_p_dg_correlation(c: Any, dg_physical: Any = None, certificate: Any = None, source_receipts: Any = None, **aliases: Any) -> dict[str, Any]:
    return build(c, dg_physical, certificate, source_receipts, **aliases)


enclose = build
evaluate = build



