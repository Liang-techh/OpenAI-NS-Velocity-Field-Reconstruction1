"""Read-only source packet for the remaining historical F endpoint roles.

This leaf is intentionally narrower than the prior historical gap lookup.  It
authenticates the saved HIGH C3 archive and the saved F archive/check, reads
the HIGH ``endpoints.gap`` tensor as opaque ordinary-Z row metadata, and
looks up only the six endpoint fields that the prior leaf did not inspect.
The two previously inspected fields (and the earlier twenty-role result) are
carried from the authenticated prior result; their F rewrite map is not
queried again here.

No source graph is opened or traversed, no old helper or constructor is
imported, and no opaque id is interpreted as a value.  A missing saved map key
stays unresolved.  The packet is evidence about saved source records only; it
does not claim an actual rewritten F field, capture-free substitution,
ordinary-Z derivative transport, endpoint equality, numerical value, bound or
error reduction.
"""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "historical_F_full_endpoint_roles_v1"

PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
HIGH_ARCHIVE = Path(PREFIX + "n1_pulse_pressure_C3_functions.json.gz")
HIGH_CHECK = Path(PREFIX + "n1_pulse_pressure_C3_functions_check.json")
F_ARCHIVE = Path(PREFIX + "Rp_native_C5_pulse_source.json.gz")
F_CHECK = Path(PREFIX + "Rp_native_C5_pulse_source_check.json")

HIGH_ARCHIVE_SHA256 = (
    "bdfcc70f57d8cea981f508074690a7e940f884382c2452805f126a3f96bdb59d"
)
HIGH_CHECK_SHA256 = (
    "cbac4685fc9031a2ed5c5a7ef624b438e1a4ce3f32b85b569512aECB1d38feb8".lower()
)
F_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
F_CHECK_SHA256 = (
    "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15"
)

PRIOR_LOOKUP_REPORT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_historical_frame_gap_saved_rewrite_lookup.json"
)
PRIOR_LOOKUP_REPORT_SHA256 = (
    "f9e6eff7457afbe591323268655be122455de98b14bf320a235c80137be75b00"
)
PRIOR_LOOKUP_ARTIFACT = "tmp/historical_frame_gap_saved_rewrite_lookup_v1/result.json"
PRIOR_LOOKUP_RESULT_SHA256 = (
    "44da94d04c52e020799692049e8a6c58f1585a1ea4f794a0146ba933c71cfc4d"
)

# These are static source anchors, not imports.  Their bytes are authenticated
# when the packet is built, but their graph-producing code is never executed.
ENDPOINT_DEFINITION_SOURCE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions.py"
)
ENDPOINT_DEFINITION_SOURCE_SHA256 = (
    "282f03acfbde22be550cce46d083f33545c8d2d550bfaeea8d1cdd8e83ed531c"
)
HIGH_BINDING_SOURCE = Path(
    PREFIX + "n1_pulse_pressure_C3_functions.py"
)
HIGH_BINDING_SOURCE_SHA256 = (
    "3a41add18b7f9bfb63ab3c484bea82c321c514aeddd59f46a5c41ca2c9edaad2"
)
PHYSICAL_MEASURE_SOURCE = Path(PREFIX + "n1_pulse_pressure_functions.py")
PHYSICAL_MEASURE_SOURCE_SHA256 = (
    "8cc7442487dbbd3054e29211ea799a04d27c43fac3f060f74f5fadd929737594"
)
F_REWRITE_SOURCE = Path(PREFIX + "Rp_native_C5_pulse_source.py")
F_REWRITE_SOURCE_SHA256 = (
    "369e6404d885dd4548e4f6451f2d5a5ed4160c0abb79a4aca2849520ede03541"
)
F_RADIUS_SOURCE = Path(PREFIX + "n1_pulse_own_fields.py")
F_RADIUS_SOURCE_SHA256 = (
    "1400066d2da5373a95f6cd30edad19359256211d0358b454203a92fac53c5ced"
)

ENDPOINT_FIELDS = (
    "R_q",
    "R_squared_q",
    "R_q_squared",
    "R_squared_q_squared",
    "R_U_squared",
    "R_squared_U_squared",
    "q_plus_q_y",
    "R_q_y",
)
PREVIOUS_ENDPOINT_FIELDS = ("q_plus_q_y", "R_q_y")
NEW_ENDPOINT_FIELDS = tuple(
    field for field in ENDPOINT_FIELDS if field not in PREVIOUS_ENDPOINT_FIELDS
)
SIDES = ("lower", "upper")
ORDERS = (0, 1, 2, 3)
PROTECTED = frozenset({0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346})

HIGH_REQUIRED_RECEIPT_FLAGS = (
    "current_original_n1_actual_native_and_Gamma_C5_functions_installed",
    "current_original_n1_actual_leading_pulse_U_mean_C5_Q_C4_functions_installed",
    "current_original_n1_true_pulse_pressure_C3_functions_installed",
)
F_REQUIRED_RECEIPT_FLAGS = (
    "current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed",
    "current_original_actual_Rp_native_C5_all_eight_source_functions_consumed",
    "current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions",
)

INPUT_PINS = {
    str(HIGH_ARCHIVE): HIGH_ARCHIVE_SHA256,
    str(HIGH_CHECK): HIGH_CHECK_SHA256,
    str(F_ARCHIVE): F_ARCHIVE_SHA256,
    str(F_CHECK): F_CHECK_SHA256,
    str(PRIOR_LOOKUP_REPORT): PRIOR_LOOKUP_REPORT_SHA256,
}

STATIC_SOURCE_PINS = {
    str(ENDPOINT_DEFINITION_SOURCE): ENDPOINT_DEFINITION_SOURCE_SHA256,
    str(HIGH_BINDING_SOURCE): HIGH_BINDING_SOURCE_SHA256,
    str(PHYSICAL_MEASURE_SOURCE): PHYSICAL_MEASURE_SOURCE_SHA256,
    str(F_REWRITE_SOURCE): F_REWRITE_SOURCE_SHA256,
    str(F_RADIUS_SOURCE): F_RADIUS_SOURCE_SHA256,
}


class HistoricalFullEndpointRolesInputError(ValueError):
    """A frozen direct input is absent, malformed or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise HistoricalFullEndpointRolesInputError(message)


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _relative(path: Path) -> str:
    return path.as_posix()


def _safe_path(root: Path, relative: Path) -> Path:
    root = root.resolve()
    _require(not relative.is_absolute() and ".." not in relative.parts,
             "frozen input path must be relative and contained")
    path = (root / relative).resolve()
    _require(path.is_relative_to(root), "frozen input path escapes repository")
    return path


def _read_authenticated_bytes(root: Path, relative: Path, expected: str) -> bytes:
    path = _safe_path(root, relative)
    _require(path.is_file(), f"missing frozen input: {_relative(relative)}")
    raw = path.read_bytes()
    _require(_sha256(raw) == expected,
             f"raw SHA-256 mismatch: {_relative(relative)}")
    return raw


def _read_json_bytes(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise HistoricalFullEndpointRolesInputError(
            f"invalid JSON direct input: {label}"
        ) from exc
    _require(isinstance(value, Mapping), f"JSON input is not an object: {label}")
    return value


def _direct(value: Any, *keys: str) -> Any:
    current = value
    for key in keys:
        _require(isinstance(current, Mapping) and key in current,
                 "missing direct field: " + ".".join(keys))
        current = current[key]
    return current


def _read_receipt(
    raw: bytes,
    label: str,
    archive_name: str,
    archive_sha256: str,
    required_flags: tuple[str, ...],
) -> Mapping[str, Any]:
    receipt = _read_json_bytes(raw, label)
    _require(receipt.get("all_passed") is True,
             f"{label} receipt all_passed is not true")
    for flag in required_flags:
        _require(receipt.get(flag) is True,
                 f"{label} required receipt flag is absent: {flag}")
    input_hashes = receipt.get("input_hashes")
    _require(isinstance(input_hashes, Mapping),
             f"{label} receipt input_hashes is missing")
    _require(input_hashes.get(archive_name) == archive_sha256,
             f"{label} receipt own archive pin mismatch")
    return receipt


def _decode_archive(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
        value = json.loads(decoded.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise HistoricalFullEndpointRolesInputError(
            f"invalid authenticated gzip archive: {label}"
        ) from exc
    _require(isinstance(value, Mapping), f"archive is not an object: {label}")
    return value


def _authenticate_static_sources(root: Path) -> dict[str, str]:
    authenticated: dict[str, str] = {}
    for relative_name, expected in STATIC_SOURCE_PINS.items():
        raw = _read_authenticated_bytes(root, Path(relative_name), expected)
        authenticated[relative_name] = _sha256(raw)
    return authenticated


def _read_prior_result(root: Path) -> tuple[Mapping[str, Any], dict[str, Any]]:
    raw_report = _read_authenticated_bytes(
        root, PRIOR_LOOKUP_REPORT, PRIOR_LOOKUP_REPORT_SHA256
    )
    report = _read_json_bytes(raw_report, _relative(PRIOR_LOOKUP_REPORT))
    _require(
        report.get("component") == "historical_frame_gap_saved_rewrite_lookup_v1",
        "prior lookup report component changed",
    )
    artifacts = report.get("artifacts")
    _require(isinstance(artifacts, Mapping), "prior lookup report artifacts missing")
    artifact = artifacts.get(PRIOR_LOOKUP_ARTIFACT)
    _require(isinstance(artifact, Mapping), "prior lookup result artifact missing")
    encoded = artifact.get("base64")
    _require(isinstance(encoded, str), "prior lookup result base64 missing")
    try:
        result_raw = base64.b64decode(encoded, validate=True)
    except (TypeError, ValueError) as exc:
        raise HistoricalFullEndpointRolesInputError(
            "prior lookup result base64 is malformed"
        ) from exc
    _require(_sha256(result_raw) == PRIOR_LOOKUP_RESULT_SHA256,
             "prior lookup result SHA-256 mismatch")
    packet = _read_json_bytes(result_raw, PRIOR_LOOKUP_ARTIFACT)
    _require(packet.get("status") == "returned",
             "prior lookup result was not returned")
    result = packet.get("result")
    _require(isinstance(result, Mapping), "prior lookup result object missing")
    _require(result.get("schema") == "historical_frame_gap_saved_rewrite_lookup_v1",
             "prior lookup schema changed")
    _require(result.get("frame_archive_raw_sha256") == F_ARCHIVE_SHA256,
             "prior lookup did not use the pinned F archive")
    _require(
        result.get("map_owner")
        == "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
        "prior lookup map owner changed",
    )
    roles = result.get("roles")
    _require(isinstance(roles, Mapping), "prior lookup roles missing")
    for field in PREVIOUS_ENDPOINT_FIELDS:
        for side in SIDES:
            name = f"{side}_{field}"
            rows = roles.get(name)
            _require(isinstance(rows, list) and len(rows) == 4,
                     f"prior lookup role is not a four-row jet: {name}")
    _require(result.get("roles_missing_count") == 20,
             "prior lookup twenty-role scope changed")
    return result, {
        "report_path": _relative(PRIOR_LOOKUP_REPORT),
        "report_raw_sha256": PRIOR_LOOKUP_REPORT_SHA256,
        "result_artifact": PRIOR_LOOKUP_ARTIFACT,
        "result_raw_sha256": PRIOR_LOOKUP_RESULT_SHA256,
        "raw_hash_authenticated_before_parse": True,
        "saved_selected_endpoint_fields_reused": list(PREVIOUS_ENDPOINT_FIELDS),
        "direct_map_lookup_repeated": False,
    }


def _opaque_rows(value: Any, source_role: str) -> list[int]:
    _require(isinstance(value, list) and len(value) == 4,
             f"{source_role} must contain four ordinary-Z rows")
    _require(all(type(item) is int for item in value),
             f"{source_role} rows must be opaque integer handles")
    return list(value)


def _read_high_endpoint_rows(
    high_bindings: Mapping[str, Any],
) -> dict[str, dict[str, list[int]]]:
    endpoints = _direct(high_bindings, "endpoints", "gap")
    _require(isinstance(endpoints, list) and len(endpoints) == 2,
             "HIGH endpoints.gap must contain lower and upper endpoints")
    result: dict[str, dict[str, list[int]]] = {}
    for side, endpoint in zip(SIDES, endpoints):
        _require(isinstance(endpoint, Mapping),
                 f"HIGH endpoints.gap {side} endpoint is malformed")
        side_rows: dict[str, list[int]] = {}
        for field in ENDPOINT_FIELDS:
            pointer = (
                "$.C5_pulse_C3_pressure_bindings.endpoints.gap["
                + ("0" if side == "lower" else "1")
                + "]."
                + field
            )
            side_rows[field] = _opaque_rows(endpoint.get(field), pointer)
        result[side] = side_rows
    return result


def _prior_rows(prior: Mapping[str, Any], side: str, field: str) -> list[Mapping[str, Any]]:
    roles = _direct(prior, "roles")
    rows = roles.get(f"{side}_{field}")
    _require(isinstance(rows, list) and len(rows) == 4,
             f"prior role rows missing: {side}_{field}")
    _require(all(isinstance(row, Mapping) for row in rows),
             f"prior role records malformed: {side}_{field}")
    return list(rows)


def _reuse_prior_selected_rows(
    prior: Mapping[str, Any],
    high_rows: Mapping[str, Mapping[str, list[int]]],
) -> dict[str, Any]:
    reused: dict[str, Any] = {}
    for side in SIDES:
        for field in PREVIOUS_ENDPOINT_FIELDS:
            name = f"{side}_{field}"
            rows = _prior_rows(prior, side, field)
            originals = [row.get("original_ref") for row in rows]
            _require(originals == list(high_rows[side][field]),
                     f"prior selected rows disagree with HIGH metadata: {name}")
            _require(
                all(row.get("status") == "not_in_saved_rewrite_map" for row in rows),
                f"prior selected role status changed: {name}",
            )
            normalized_records = []
            for order, row in enumerate(rows):
                normalized = dict(row)
                normalized.update({
                    "field": field,
                    "side": side,
                    "ordinary_Z_order": order,
                    "source_namespace": "HIGH.C5_pulse_C3_pressure_bindings",
                    "source_owner": "HIGH:n1_pulse_pressure_C3_functions",
                    "map_namespace": "F.Rp_native_C5_pulse_source",
                    "map_owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
                    "mapped_ref": row.get("rewritten_ref"),
                    "mapped_handle": None,
                    "source_value_selected": False,
                    "actual_F_field_claimed": False,
                    "capture_free_substitution_proved": False,
                    "ordinary_Z_transport_proved": False,
                })
                normalized_records.append(normalized)
            reused[name] = {
                "source_namespace": "HIGH.C5_pulse_C3_pressure_bindings",
                "source_owner": "HIGH:n1_pulse_pressure_C3_functions",
                "source_role": (
                    "$.C5_pulse_C3_pressure_bindings.endpoints.gap."
                    + ("0" if side == "lower" else "1")
                    + "."
                    + field
                ),
                "ordinary_Z_orders": list(ORDERS),
                "records": normalized_records,
                "prior_records": [dict(row) for row in rows],
                "status": "reused_from_prior_saved_lookup",
                "direct_map_lookup_repeated": False,
                "source_value_selected": False,
                "actual_F_field_claimed": False,
            }
    return reused


def _lookup_new_field(
    original_ref: int,
    side: str,
    field: str,
    order: int,
    rewrite_map: Mapping[str, Any],
) -> dict[str, Any]:
    key = str(original_ref)
    source_role = (
        "$.C5_pulse_C3_pressure_bindings.endpoints.gap."
        + ("0" if side == "lower" else "1")
        + "."
        + field
        + f"[{order}]"
    )
    record: dict[str, Any] = {
        "field": field,
        "side": side,
        "ordinary_Z_order": order,
        "source_namespace": "HIGH.C5_pulse_C3_pressure_bindings",
        "source_owner": "HIGH:n1_pulse_pressure_C3_functions",
        "source_role": source_role,
        "original_ref": original_ref,
        "map_namespace": "F.Rp_native_C5_pulse_source",
        "map_owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
        "map_key": key,
        "source_value_selected": False,
        "actual_F_field_claimed": False,
        "capture_free_substitution_proved": False,
        "ordinary_Z_transport_proved": False,
    }
    if key not in rewrite_map:
        record.update({
            "status": "unresolved_missing_from_saved_F_rewrite_map",
            "mapped_ref": None,
            "mapped_handle": None,
            "missing_mapping_is_not_identity": True,
        })
        return record
    mapped = rewrite_map[key]
    _require(type(mapped) is int,
             f"saved F rewrite target is not an opaque integer: {source_role}")
    record.update({
        "status": "present_in_saved_F_rewrite_map",
        "mapped_ref": mapped,
        "mapped_handle": {
            "source_node_id": mapped,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
        },
        "missing_mapping_is_not_identity": False,
    })
    return record


def _full_source_anchors() -> dict[str, Any]:
    return {
        "eight_field_original_definition": {
            "path": _relative(ENDPOINT_DEFINITION_SOURCE),
            "sha256": ENDPOINT_DEFINITION_SOURCE_SHA256,
            "lines": "35-39",
            "fields": list(ENDPOINT_FIELDS),
            "formulas": {
                "R_q": "R*Q",
                "R_squared_q": "R^2*Q",
                "R_q_squared": "R*Q^2",
                "R_squared_q_squared": "R^2*Q^2",
                "R_U_squared": "R*U^2",
                "R_squared_U_squared": "R^2*U^2",
                "q_plus_q_y": "Q+Q_y",
                "R_q_y": "R*Q_y",
            },
            "function_value_selected": False,
        },
        "HIGH_endpoint_binding": {
            "path": _relative(HIGH_BINDING_SOURCE),
            "sha256": HIGH_BINDING_SOURCE_SHA256,
            "lines": "159-176",
            "definition": "C3 stores per-chart and global finite endpoint bindings",
            "physical_Pstar_radial_measures_and_finite_endpoints_retained": True,
            "function_value_selected": False,
        },
        "physical_measure": {
            "path": _relative(PHYSICAL_MEASURE_SOURCE),
            "sha256": PHYSICAL_MEASURE_SOURCE_SHA256,
            "lines": "27-63",
            "measures": {
                "main_gap_I0": "R/mu dxi",
                "main_gap_I1": "R^2/mu dxi",
                "end_I0": "R ds",
                "end_I1": "R^2 ds",
            },
            "endpoint_independent_of_Z": True,
            "physical_Pstar_not_reapplied": True,
            "function_value_selected": False,
        },
        "F_saved_rewrite_map": {
            "path": _relative(F_REWRITE_SOURCE),
            "sha256": F_REWRITE_SOURCE_SHA256,
            "lines": "90-116",
            "map_owner": "Rp_native_C5_pulse_source_bindings.rewritten_nodes",
            "operation": "saved direct original-ref key lookup only",
            "graph_traversal": False,
            "function_value_selected": False,
        },
        "historical_F_radius_template": {
            "path": _relative(F_RADIUS_SOURCE),
            "sha256": F_RADIUS_SOURCE_SHA256,
            "lines": "132-140",
            "formula": "R=exp(logRp+xi/mu)",
            "owner": "original own-field template",
            "rewritten_F_radius_proved": False,
            "function_value_selected": False,
        },
    }


def _negative_fixture() -> dict[str, Any]:
    """Describe unexecuted negative controls for future callers/tests."""

    return {
        "schema": "historical_F_full_endpoint_roles_negative_fixture_v1",
        "executed": False,
        "cases": {
            "missing_map_key_is_not_identity": {
                "original_ref": 123,
                "map": {},
                "expected_status": "unresolved_missing_from_saved_F_rewrite_map",
                "must_not_synthesize_mapped_ref": True,
            },
            "equal_opaque_handles_do_not_prove_function_equality": {
                "original_ref": 33982,
                "mapped_ref": 33982,
                "must_not_admit_actual_F_radius": True,
            },
            "previous_two_fields_are_not_relooked_up": {
                "fields": list(PREVIOUS_ENDPOINT_FIELDS),
                "expected_direct_map_lookup_repeated": False,
            },
            "short_ordinary_Z_jet_is_rejected": {
                "row_count": 3,
                "must_be_rejected": True,
            },
        },
    }


def build_historical_full_endpoint_roles(root: str | Path = ".") -> dict[str, Any]:
    """Build the source-only 48-role historical F endpoint packet."""

    repo_root = Path(root)

    # Every raw archive/check is hashed before parsing or gzip decoding.
    high_raw = _read_authenticated_bytes(
        repo_root, HIGH_ARCHIVE, HIGH_ARCHIVE_SHA256
    )
    high_check_raw = _read_authenticated_bytes(
        repo_root, HIGH_CHECK, HIGH_CHECK_SHA256
    )
    f_raw = _read_authenticated_bytes(repo_root, F_ARCHIVE, F_ARCHIVE_SHA256)
    f_check_raw = _read_authenticated_bytes(repo_root, F_CHECK, F_CHECK_SHA256)
    _read_authenticated_bytes(repo_root, PRIOR_LOOKUP_REPORT, PRIOR_LOOKUP_REPORT_SHA256)
    static_hashes = _authenticate_static_sources(repo_root)

    high_check = _read_receipt(
        high_check_raw,
        _relative(HIGH_CHECK),
        HIGH_ARCHIVE.name,
        HIGH_ARCHIVE_SHA256,
        HIGH_REQUIRED_RECEIPT_FLAGS,
    )
    f_check = _read_receipt(
        f_check_raw,
        _relative(F_CHECK),
        F_ARCHIVE.name,
        F_ARCHIVE_SHA256,
        F_REQUIRED_RECEIPT_FLAGS,
    )
    f_input_hashes = f_check.get("input_hashes")
    _require(isinstance(f_input_hashes, Mapping), "F receipt input_hashes missing")
    _require(f_input_hashes.get(HIGH_ARCHIVE.name) == HIGH_ARCHIVE_SHA256,
             "F receipt HIGH archive cross-pin mismatch")

    high_archive = _decode_archive(high_raw, _relative(HIGH_ARCHIVE))
    f_archive = _decode_archive(f_raw, _relative(F_ARCHIVE))
    high_family = high_archive.get("source_family")
    f_family = f_archive.get("source_family")
    _require(high_family == high_check.get("source_family"),
             "HIGH archive/check source family mismatch")
    _require(f_family == f_check.get("source_family"),
             "F archive/check source family mismatch")
    _require(high_family == f_family,
             "HIGH/F source family mismatch")

    high_bindings = _direct(high_archive, "C5_pulse_C3_pressure_bindings")
    f_bindings = _direct(f_archive, "Rp_native_C5_pulse_source_bindings")
    _direct(f_bindings, "inputs", HIGH_ARCHIVE.name)
    rewrite_map = _direct(f_bindings, "rewritten_nodes")
    _require(isinstance(rewrite_map, Mapping), "F saved rewrite map is missing")
    high_rows = _read_high_endpoint_rows(high_bindings)

    prior, prior_input = _read_prior_result(repo_root)
    prior_reuse = _reuse_prior_selected_rows(prior, high_rows)

    new_roles: list[dict[str, Any]] = []
    # Reuse the prior direct lookup records if a different endpoint field
    # shares a handle; never probe a protected handle or repeat a lookup.
    prior_by_ref = {}
    for value in _direct(prior, "roles").values():
        items = value if isinstance(value, list) else [value]
        for item in items:
            if isinstance(item, Mapping) and type(item.get("original_ref")) is int:
                ref = item["original_ref"]
                if ref in prior_by_ref:
                    old = prior_by_ref[ref]
                    _require((old.get("status"), old.get("rewritten_ref")) ==
                             (item.get("status"), item.get("rewritten_ref")),
                             "conflicting prior saved lookup records")
                prior_by_ref[ref] = item
    new_lookup_cache = {}
    for side in SIDES:
        for field in NEW_ENDPOINT_FIELDS:
            for order in ORDERS:
                ref = high_rows[side][field][order]
                if ref in prior_by_ref or ref in PROTECTED or ref in new_lookup_cache:
                    if ref in prior_by_ref:
                        status, previous = "prior_saved_lookup_reused", prior_by_ref[ref]
                    elif ref in PROTECTED:
                        status, previous = "protected_lookup_skipped", None
                    else:
                        status, previous = "new_lookup_reused", new_lookup_cache[ref]
                    new_roles.append({"field": field, "side": side, "ordinary_Z_order": order,
                        "original_ref": ref, "status": status,
                        "prior_saved_record": previous, "source_value_selected": False,
                        "actual_F_field_claimed": False, "capture_free_substitution_proved": False,
                        "ordinary_Z_transport_proved": False})
                else:
                    record = _lookup_new_field(ref, side, field, order, rewrite_map)
                    new_lookup_cache[ref] = record
                    new_roles.append(record)

    present_count = sum(
        record["status"] == "present_in_saved_F_rewrite_map"
        for record in new_roles
    )
    missing_count = sum(
        record["status"] == "unresolved_missing_from_saved_F_rewrite_map"
        for record in new_roles
    )
    reused_count = sum(record["status"] in {"prior_saved_lookup_reused", "protected_lookup_skipped", "new_lookup_reused"} for record in new_roles)
    _require(len(new_roles) == 48, "new endpoint role scope is not 48 records")
    _require(present_count + missing_count + reused_count == 48,
             "new endpoint role status count is incomplete")

    high_definition = _direct(high_bindings, "definition")
    _require(
        high_definition.get("physical_Pstar_radial_measures_and_finite_endpoints_retained")
        is True,
        "HIGH finite-endpoint/physical-measure definition flag is absent",
    )

    prior_roles = _direct(prior, "roles")
    return {
        "schema": SCHEMA,
        "status": "historical_F_full_endpoint_saved_map_lookup_source_only",
        "source_only": True,
        "raw_hashes_authenticated_before_gzip_decode": True,
        "gzip_decoded_after_raw_authentication": True,
        "source_graph_traversed": False,
        "old_helper_imported": False,
        "old_leaf_opened": False,
        "source_values_computed": False,
        "source_value_selected": False,
        "required_inputs": {
            "high_archive": {
                "path": _relative(HIGH_ARCHIVE),
                "raw_sha256": HIGH_ARCHIVE_SHA256,
                "archive_kind": "HIGH C3 source archive",
            },
            "high_check": {
                "path": _relative(HIGH_CHECK),
                "raw_sha256": HIGH_CHECK_SHA256,
                "required_flags": list(HIGH_REQUIRED_RECEIPT_FLAGS),
            },
            "f_archive": {
                "path": _relative(F_ARCHIVE),
                "raw_sha256": F_ARCHIVE_SHA256,
                "archive_kind": "F Rp-native C5 source archive",
            },
            "f_check": {
                "path": _relative(F_CHECK),
                "raw_sha256": F_CHECK_SHA256,
                "required_flags": list(F_REQUIRED_RECEIPT_FLAGS),
                "HIGH_archive_cross_pin": HIGH_ARCHIVE_SHA256,
            },
            "prior_saved_lookup": prior_input,
        },
        "receipt_authentication": {
            "HIGH_all_passed": high_check.get("all_passed") is True,
            "F_all_passed": f_check.get("all_passed") is True,
            "HIGH_source_family": high_family,
            "F_source_family": f_family,
            "same_source_family": high_family == f_family,
            "F_embedded_HIGH_input_present": True,
            "HIGH_required_flags": {
                flag: high_check.get(flag) for flag in HIGH_REQUIRED_RECEIPT_FLAGS
            },
            "F_required_flags": {
                flag: f_check.get(flag) for flag in F_REQUIRED_RECEIPT_FLAGS
            },
        },
        "source_anchors": _full_source_anchors(),
        "source_definition_metadata": {
            "HIGH_endpoint_definition": high_definition,
            "endpoint_fields": list(ENDPOINT_FIELDS),
            "ordinary_Z_orders": list(ORDERS),
            "endpoint_order": ["lower", "upper"],
            "declared_direction": "lower_to_upper",
            "orientation_identity_proved": False,
            "source_value_selected": False,
        },
        "authenticated_static_source_sha256": static_hashes,
        "HIGH_endpoint_rows": {
            side: {field: list(rows) for field, rows in fields.items()}
            for side, fields in high_rows.items()
        },
        "prior_saved_lookup_reuse": {
            "prior_result_sha256": PRIOR_LOOKUP_RESULT_SHA256,
            "map_owner": prior.get("map_owner"),
            "roles_present_count": prior.get("roles_present_count"),
            "roles_missing_count": prior.get("roles_missing_count"),
            "selected_endpoint_fields": list(PREVIOUS_ENDPOINT_FIELDS),
            "direct_map_lookup_repeated": False,
            "records": prior_reuse,
            "field_radius_and_gap_m15_prior_roles": {
                "field_radius": prior_roles.get("field_radius"),
                "gap_m15": prior_roles.get("gap_m15"),
            },
        },
        "new_unchecked_endpoint_roles": {
            "fields": list(NEW_ENDPOINT_FIELDS),
            "sides": list(SIDES),
            "ordinary_Z_orders": list(ORDERS),
            "total": len(new_roles),
            "present_in_saved_F_rewrite_map": present_count,
            "unresolved_missing_from_saved_F_rewrite_map": missing_count,
            "prior_saved_lookup_reused": sum(r["status"] == "prior_saved_lookup_reused" for r in new_roles),
            "protected_lookup_skipped": sum(r["status"] == "protected_lookup_skipped" for r in new_roles),
            "new_lookup_reused": sum(r["status"] == "new_lookup_reused" for r in new_roles),
            "records": new_roles,
        },
        "endpoint_direction_radius_and_measure": {
            "endpoint_order_preserved": ["lower", "upper"],
            "declared_direction": "lower_to_upper",
            "orientation_proved": False,
            "physical_measure_source": _relative(PHYSICAL_MEASURE_SOURCE),
            "physical_measure_source_lines": "27-63",
            "radius_template_prior_role": prior_roles.get("field_radius"),
            "radius_formula": "R=exp(logRp+xi/mu)",
            "historical_F_radius_owner_bound": False,
            "moving_endpoint_chain_rule_proved": False,
            "common_coordinate_mu_measure_proved": False,
        },
        "error_and_comparison_preservation": {
            "old_errors_removed": False,
            "error_transport_certified": False,
            "actual_error_reduction_proved": False,
            "independent_reference_offset_removed": False,
            "physical_correction_roots_removed": False,
            "all_prior_error_terms_retained": True,
            "formal_P_minus_F_decomposition": "(P-G)+(G-F)",
            "P_equals_F_proved": False,
            "G_equals_F_proved": False,
            "P_equals_G_proved": False,
        },
        "source_gates": {
            "HIGH_eight_field_tensor_read": True,
            "new_48_role_records_no_repeated_or_protected_lookups": True,
            "prior_20_roles_reused_without_lookup": True,
            "F_historical_radius_owner_bound": False,
            "F_endpoint_owner_bound": False,
            "full_eight_field_tensor_bound": False,
            "capture_free_inlet_substitution": False,
            "ordinary_Z0_to_Z3_transport": False,
            "same_source_endpoint_seams": False,
            "local_to_global_m15_identity": False,
            "common_coordinate_measure_orientation": False,
            "physical_scale_and_error_transport": False,
            "P_equals_F": False,
            "G_equals_F": False,
            "P_equals_G": False,
            "numeric_use_admissible": False,
        },
        "negative_fixture": _negative_fixture(),
        "admissible_for_numeric_use": False,
        "full_goal_complete": False,
    }


def recover(root: str | Path = ".") -> dict[str, Any]:
    """Alias for the source-only full endpoint role packet."""

    return build_historical_full_endpoint_roles(root)


def build_full_endpoint_roles(root: str | Path = ".") -> dict[str, Any]:
    """Short interface alias for callers that name the endpoint scope directly."""

    return build_historical_full_endpoint_roles(root)


def synthetic_negative_fixture() -> dict[str, Any]:
    """Expose, but do not execute, the in-memory negative fixture."""

    return _negative_fixture()


__all__ = [
    "ENDPOINT_FIELDS",
    "F_ARCHIVE",
    "F_ARCHIVE_SHA256",
    "F_CHECK",
    "F_CHECK_SHA256",
    "INPUT_PINS",
    "HistoricalFullEndpointRolesInputError",
    "HIGH_ARCHIVE",
    "HIGH_ARCHIVE_SHA256",
    "HIGH_CHECK",
    "HIGH_CHECK_SHA256",
    "NEW_ENDPOINT_FIELDS",
    "PREVIOUS_ENDPOINT_FIELDS",
    "PRIOR_LOOKUP_REPORT",
    "PRIOR_LOOKUP_REPORT_SHA256",
    "PRIOR_LOOKUP_RESULT_SHA256",
    "SCHEMA",
    "STATIC_SOURCE_PINS",
    "build_full_endpoint_roles",
    "build_historical_full_endpoint_roles",
    "recover",
    "synthetic_negative_fixture",
]
