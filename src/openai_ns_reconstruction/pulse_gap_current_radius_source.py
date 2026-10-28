"""Source-only consumer for the saved current-radius pulse-gap chart.

This leaf deliberately consumes saved source handles.  It does not import a
source constructor, select a field value, walk a graph, replay an ancestor,
or evaluate a pressure/moment expression.  The owner introduced here is
``G``: the current-radius chart reconstructed from the saved actual field
rows.  ``G`` is not renamed to historical ``F`` and is not identified with
the prior physical ``P`` source.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import copy
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "pulse_gap_current_radius_source_v1"
ARCHIVE_RELATIVE_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.json.gz"
)
RECEIPT_RELATIVE_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check.json"
)
PREVIOUS_RESULT_RELATIVE_PATH = Path("tmp/pulse_gap_actual_source_program_v1/result.json")
SOURCE_DEFINITION_RELATIVE_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py"
)
CHECKER_SOURCE_RELATIVE_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check.py"
)
FRAME_ARCHIVE_RELATIVE_PATH = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)

EXPECTED_ARCHIVE_SHA256 = "e3cd768f0508fd7a315a058077e68bbaf784198464a6b75d4cdd33de28037bf8"
EXPECTED_RECEIPT_SHA256 = "ea5eebb9483916932c9e6e43773d061e3ea8f490d284a9291e9fc1e38e147ee4"
EXPECTED_PREVIOUS_RESULT_SHA256 = "e611cdc7e049ce40e6667275b4e9706c3478044880b92c93baebf8ab8d06ebf2"
EXPECTED_FRAME_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"

FIELD_NAMES = ("Q", "Q_y", "m", "m_y", "U", "U_y")
DERIVATIVE_FIELD_NAMES = ("m_y", "m_yy", "Q", "Q_y", "Q_yy")
REQUIRED_RECEIPT_GATES = (
    "current_original_pulse_Uyy_C5_Qyy_C4_source_functions_installed",
    "current_original_pulse_radial_N1p_pressure_C3_source_functions_installed",
    "current_original_n1_partial_pulse_pressure_physical_y_FTC_C3_functions_installed",
)


class SourceProgramInputError(ValueError):
    """Raised internally when a required saved source packet is absent."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _path_string(path: Path) -> str:
    return path.as_posix()


def _read_json_after_raw_auth(
    path: Path, expected_sha256: str, label: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if not path.is_file():
        raise SourceProgramInputError(f"missing {label}: {_path_string(path)}")
    raw = path.read_bytes()
    actual_sha256 = _sha256(raw)
    if actual_sha256 != expected_sha256:
        raise SourceProgramInputError(
            f"{label} raw SHA-256 mismatch: {actual_sha256} != {expected_sha256}"
        )
    try:
        decoded = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceProgramInputError(f"invalid JSON in {label}") from exc
    if not isinstance(decoded, Mapping):
        raise SourceProgramInputError(f"{label} top level is not an object")
    return decoded, {
        "path": _path_string(path),
        "raw_sha256": actual_sha256,
        "raw_authenticated_before_decode": True,
    }


def _read_gzip_json_after_raw_auth(
    path: Path, expected_sha256: str, label: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if not path.is_file():
        raise SourceProgramInputError(f"missing {label}: {_path_string(path)}")
    raw = path.read_bytes()
    actual_sha256 = _sha256(raw)
    if actual_sha256 != expected_sha256:
        raise SourceProgramInputError(
            f"{label} raw SHA-256 mismatch: {actual_sha256} != {expected_sha256}"
        )
    try:
        # The raw-byte hash check is intentionally complete before this decode.
        decoded_bytes = gzip.decompress(raw)
        decoded = json.loads(decoded_bytes.decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SourceProgramInputError(f"invalid gzip JSON in {label}") from exc
    if not isinstance(decoded, Mapping):
        raise SourceProgramInputError(f"{label} top level is not an object")
    return decoded, {
        "path": _path_string(path),
        "raw_sha256": actual_sha256,
        "raw_authenticated_before_gzip_decode": True,
    }


def _require_mapping(value: Any, pointer: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise SourceProgramInputError(f"expected object at {pointer}")
    return value


def _require_sequence(value: Any, pointer: str, minimum: int = 1) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise SourceProgramInputError(f"expected row list at {pointer}")
    if len(value) < minimum:
        raise SourceProgramInputError(
            f"too few saved rows at {pointer}: {len(value)} < {minimum}"
        )
    return value


def _source_ref(value: Any, pointer: str) -> dict[str, Any]:
    # Integers in these encoded packets are opaque graph handles, never field
    # values.  bool is rejected because it is an int subclass in Python.
    if isinstance(value, bool) or not isinstance(value, int):
        raise SourceProgramInputError(f"non-literal source handle at {pointer}")
    return {
        "status": "direct_source_ref",
        "source_role": pointer,
        "source_node_id": value,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _direct_tree(value: Any, pointer: str) -> Any:
    """Copy only a saved literal/source-reference tree, without graph access."""

    if isinstance(value, Mapping):
        return {
            str(key): _direct_tree(item, f"{pointer}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, (str, bool)) or value is None:
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return _source_ref(value, pointer)
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return [
            _direct_tree(item, f"{pointer}[{index}]")
            for index, item in enumerate(value)
        ]
    raise SourceProgramInputError(f"unsupported saved packet at {pointer}")


def _direct_rows(value: Any, pointer: str, minimum: int = 3) -> dict[str, Any]:
    rows = _require_sequence(value, pointer, minimum)
    refs = [_source_ref(row, f"{pointer}[{index}]") for index, row in enumerate(rows)]
    return {
        "status": "direct_literal_row_list",
        "source_role": pointer,
        "ordinary_Z_orders": list(range(len(refs))),
        "ordinary_Z0_to_Z2": refs[:3],
        "all_saved_rows": refs,
        "source_values_selected": False,
    }


def _node_ids(row_packet: Mapping[str, Any]) -> list[int]:
    rows = row_packet.get("ordinary_Z0_to_Z2")
    if not isinstance(rows, Sequence):
        raise SourceProgramInputError("malformed direct row packet")
    result: list[int] = []
    for row in rows:
        record = _require_mapping(row, "row")
        node_id = record.get("source_node_id")
        if isinstance(node_id, bool) or not isinstance(node_id, int):
            raise SourceProgramInputError("malformed opaque source node id")
        result.append(node_id)
    return result


def _previous_ref(item: Any, pointer: str) -> dict[str, Any]:
    record = _require_mapping(item, pointer)
    source_node_id = record.get("source_node_id")
    if isinstance(source_node_id, bool) or not isinstance(source_node_id, int):
        raise SourceProgramInputError(f"missing previous source node at {pointer}")
    # Preserve the saved source metadata, but never promote its node id to a
    # value.  Unknown extra keys are not interpreted.
    result = dict(record)
    result["source_role"] = str(record.get("source_role", pointer))
    result["source_node_id"] = source_node_id
    result["node_id_is_not_a_numeric_value"] = True
    result["source_value_selected"] = False
    return result


def _previous_rows(value: Any, pointer: str) -> dict[str, Any]:
    bundle = _require_mapping(value, pointer)
    rows = _require_sequence(bundle.get("ordinary_Z0_to_Z2"), f"{pointer}.ordinary_Z0_to_Z2", 3)
    refs = [
        _previous_ref(item, f"{pointer}.ordinary_Z0_to_Z2[{index}]")
        for index, item in enumerate(rows[:3])
    ]
    return {
        "status": "previous_saved_direct_rows",
        "source_role": str(bundle.get("source_role", pointer)),
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z0_to_Z2": refs,
        "source_values_selected": False,
    }


def _require_receipt_gate(receipt: Mapping[str, Any], key: str) -> bool:
    return receipt.get(key) is True


def _receipt_summary(receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _require_mapping(receipt.get("checks"), "$.checks")
    generic = _require_mapping(
        checks.get("generic_recovery_FTC"), "$.checks.generic_recovery_FTC"
    )
    if receipt.get("all_passed") is not True:
        raise SourceProgramInputError("saved receipt all_passed is not true")
    if generic.get("passed") is not True:
        raise SourceProgramInputError("saved generic_recovery_FTC receipt is not passed")
    gates = {key: _require_receipt_gate(receipt, key) for key in REQUIRED_RECEIPT_GATES}
    if not all(gates.values()):
        raise SourceProgramInputError("one or more required C3 source gates are false")

    commutation_rows = generic.get(
        "pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation"
    )
    q_rows = generic.get(
        "actual_Q_lower_rows_directly_proved_and_higher_rows_from_pinned_ordinary_derivative_source"
    )
    mass_rows = generic.get("actual_mass_radial_recovery_rows")
    leading_q_rows = generic.get("actual_leading_Q_radial_recovery_rows")
    base_ftc = generic.get("actual_signed_pressure_FTC_rows")
    required_mass = {f"gap_{field}_Z{order}" for field in ("m_y", "m_yy") for order in range(3)}
    required_q = {f"gap_{field}_Z{order}" for field in ("Q", "Q_y", "Q_yy") for order in range(3)}
    saved_gap_recovery_rows = {
        "mass_and_gap_U_m_recovery_rows_present": isinstance(mass_rows, Sequence)
        and not isinstance(mass_rows, (str, bytes))
        and required_mass.issubset(mass_rows),
        "leading_Q_recovery_rows_present": isinstance(leading_q_rows, Sequence)
        and not isinstance(leading_q_rows, (str, bytes))
        and required_q.issubset(leading_q_rows),
    }
    return {
        "all_passed": True,
        "generic_recovery_FTC_passed": True,
        "required_source_gates": gates,
        "theorem_rows_consumed_without_recomputation": {
            "pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation": isinstance(
                commutation_rows, Sequence
            )
            and not isinstance(commutation_rows, (str, bytes))
            and "gap_m15_y_Z1" in commutation_rows
            and isinstance(base_ftc, Sequence)
            and "gap_m15_y_Z0" in base_ftc,
            "actual_Q_lower_rows_directly_proved_and_higher_rows_from_pinned_ordinary_derivative_source": isinstance(
                q_rows, Sequence
            )
            and not isinstance(q_rows, (str, bytes))
            and len(q_rows) > 0,
        },
        "saved_gap_recovery_rows": saved_gap_recovery_rows,
        "saved_gap_Q_recovery_theorem_consumed": all(saved_gap_recovery_rows.values()),
        "source_values_selected": False,
        "recomputed": False,
    }


def _verify_receipt_source_pin(
    receipt: Mapping[str, Any], source_definition_path: Path
) -> dict[str, Any]:
    if not source_definition_path.is_file():
        raise SourceProgramInputError(
            f"missing source definition: {_path_string(source_definition_path)}"
        )
    input_hashes = receipt.get("input_hashes")
    if not isinstance(input_hashes, Mapping):
        raise SourceProgramInputError("receipt input_hashes missing")
    normalized_name = source_definition_path.name.lower()
    matches = [
        str(value)
        for key, value in input_hashes.items()
        if str(key).replace("\\", "/").lower().endswith(normalized_name)
    ]
    if not matches:
        raise SourceProgramInputError("receipt does not pin the C3 source definition")
    actual = _sha256(source_definition_path.read_bytes())
    if actual not in matches:
        raise SourceProgramInputError(
            "C3 source definition SHA-256 disagrees with the accepted receipt"
        )
    return {
        "path": _path_string(source_definition_path),
        "sha256": actual,
        "receipt_pin_verified": True,
    }


def _required_mapping_item(mapping: Mapping[str, Any], key: str, pointer: str) -> Any:
    if key not in mapping:
        raise SourceProgramInputError(f"missing saved item at {pointer}.{key}")
    return mapping[key]


def _optional_direct_tree(mapping: Mapping[str, Any], key: str, pointer: str) -> Any:
    if key not in mapping:
        return {
            "status": "unresolved_missing_direct_source_ref",
            "source_role": pointer,
            "source_values_selected": False,
        }
    return _direct_tree(mapping[key], pointer)


def _extract_previous_side(
    previous: Mapping[str, Any], side_name: str
) -> dict[str, Any]:
    result = _require_mapping(previous.get("result"), "$.result")
    sides = _require_mapping(result.get("sides"), "$.result.sides")
    side = _require_mapping(sides.get(side_name), f"$.result.sides.{side_name}")
    if side_name == "F" and side.get("source_container_role") != "$.Rp_native_C5_pulse_source_bindings":
        raise SourceProgramInputError("previous F source container role changed")
    leading = _require_mapping(
        side.get("leading_fields_gap"), f"$.result.sides.{side_name}.leading_fields_gap"
    )
    rows = {
        field: _previous_rows(
            leading.get(field),
            f"$.result.sides.{side_name}.leading_fields_gap.{field}",
        )
        for field in FIELD_NAMES
    }
    if side_name == "F":
        for field, packet in rows.items():
            expected_role = f"$.Rp_native_C5_pulse_source_bindings.leading_fields.gap.{field}"
            if packet["source_role"] != expected_role:
                raise SourceProgramInputError(f"previous F source field role changed: {field}")
    preserved: dict[str, Any] = {
        "side": side_name,
        "leading_fields_gap": rows,
        "source_values_selected": False,
    }
    for optional_key in ("geometry", "coordinate_owner", "finite_endpoints", "gap_prerequisites"):
        if optional_key in side:
            preserved[optional_key] = side[optional_key]
    return preserved


def _extract_previous_corrections(previous: Mapping[str, Any]) -> Any:
    result = _require_mapping(previous.get("result"), "$.result")
    correction = result.get("physical_pulse_correction")
    if correction is None:
        raise SourceProgramInputError("previous physical pulse corrections are missing")
    return copy.deepcopy(correction)


def _extract_gap_binding(archive: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    binding = _require_mapping(
        archive.get("pulse_radial_pressure_C3_bindings"),
        "$.pulse_radial_pressure_C3_bindings",
    )
    charts = _require_mapping(binding.get("charts"), "$.pulse_radial_pressure_C3_bindings.charts")
    gap = _require_mapping(charts.get("gap"), "$.pulse_radial_pressure_C3_bindings.charts.gap")
    endpoints = _require_mapping(
        binding.get("endpoints"), "$.pulse_radial_pressure_C3_bindings.endpoints"
    )
    gap_endpoints = _require_mapping(endpoints.get("gap"), "...endpoints.gap")
    for side in ("left", "right"):
        if side not in gap_endpoints:
            raise SourceProgramInputError(f"missing saved gap endpoint: {side}")
    return binding, gap


def _extract_gap_owner(
    binding: Mapping[str, Any], gap: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    fields = _require_mapping(gap.get("fields"), "$.pulse_radial_pressure_C3_bindings.charts.gap.fields")
    rows = {
        field: _direct_rows(
            fields.get(field),
            f"$.pulse_radial_pressure_C3_bindings.charts.gap.fields.{field}",
            3,
        )
        for field in FIELD_NAMES
    }
    derivative_rows = {
        field: _direct_rows(
            fields.get(field),
            f"$.pulse_radial_pressure_C3_bindings.charts.gap.fields.{field}",
            3,
        )
        for field in DERIVATIVE_FIELD_NAMES
    }
    density = _require_mapping(
        gap.get("density"), "$.pulse_radial_pressure_C3_bindings.charts.gap.density"
    )
    m15_y = _direct_rows(
        density.get("m15_y"),
        "$.pulse_radial_pressure_C3_bindings.charts.gap.density.m15_y",
        1,
    )
    endpoints = _require_mapping(binding.get("endpoints"), "...endpoints")
    gap_endpoints = _require_mapping(endpoints.get("gap"), "...endpoints.gap")
    endpoint_packet = {
        side: _direct_tree(gap_endpoints[side], f"$.pulse_radial_pressure_C3_bindings.endpoints.gap.{side}")
        for side in ("left", "right")
    }
    selected = {
        "owner": "G",
        "source_role": "$.pulse_radial_pressure_C3_bindings.charts.gap",
        "coordinate": _direct_tree(
            _required_mapping_item(gap, "coordinate", "$.pulse_radial_pressure_C3_bindings.charts.gap"),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.coordinate",
        ),
        "fields": rows,
        "derivative_rows": derivative_rows,
        "m15_y": m15_y,
        "m15_y_Z0_present": len(m15_y["all_saved_rows"]) >= 1,
        "m15_y_higher_rows_present": len(m15_y["all_saved_rows"]) >= 2,
        "radius": _direct_tree(
            _required_mapping_item(gap, "radius", "$.pulse_radial_pressure_C3_bindings.charts.gap"),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.radius",
        ),
        "domain": _direct_tree(
            _required_mapping_item(gap, "domain", "$.pulse_radial_pressure_C3_bindings.charts.gap"),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.domain",
        ),
        "physical_y_Jacobian": _direct_tree(
            _required_mapping_item(
                gap, "physical_y_Jacobian", "$.pulse_radial_pressure_C3_bindings.charts.gap"
            ),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.physical_y_Jacobian",
        ),
        "signed_pressure_increment": _direct_tree(
            _required_mapping_item(
                gap, "signed_pressure_increment", "$.pulse_radial_pressure_C3_bindings.charts.gap"
            ),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.signed_pressure_increment",
        ),
        "density": _direct_tree(
            density,
            "$.pulse_radial_pressure_C3_bindings.charts.gap.density",
        ),
        "FTC": _direct_tree(
            gap.get("FTC"),
            "$.pulse_radial_pressure_C3_bindings.charts.gap.FTC",
        ),
        "endpoints": {
            "gap": endpoint_packet,
            "source_values_selected": False,
        },
        "source_values_selected": False,
        "current_radius_owner_bound": True,
    }
    # Keep the chart object only as the source-owner envelope; no chart graph
    # or physical-to-frame comparison is performed here.
    return selected, {"source_role": "G", "field_rows": rows}


def _compare_direct_field_arrays(
    g_fields: Mapping[str, Any], f_fields: Mapping[str, Any]
) -> dict[str, Any]:
    comparisons: dict[str, Any] = {}
    for field in FIELD_NAMES:
        g_packet = _require_mapping(g_fields.get(field), f"G.{field}")
        f_packet = _require_mapping(f_fields.get(field), f"F.{field}")
        g_ids = _node_ids(g_packet)
        f_ids = _node_ids(f_packet)
        comparisons[field] = {
            "G_Z0_to_Z2_source_node_ids": [
                {"source_node_id": item, "node_id_is_not_a_numeric_value": True}
                for item in g_ids
            ],
            "historical_F_Z0_to_Z2_source_node_ids": [
                {"source_node_id": item, "node_id_is_not_a_numeric_value": True}
                for item in f_ids
            ],
            "direct_literal_row_list_equal": g_ids == f_ids,
            "comparison_scope": "direct field arrays only; no graph traversal",
            "source_values_selected": False,
        }
    return {
        "fields": comparisons,
        "all_six_direct_field_arrays_equal": all(
            item["direct_literal_row_list_equal"] for item in comparisons.values()
        ),
        "historical_F_pressure_or_m15_equality_proved": False,
        "source_graph_traversed": False,
        "physical_to_frame_graph_compared": False,
    }


def _scope_audit() -> dict[str, Any]:
    return {
        "raw_gzip_authentication_required_before_decode": True,
        "source_graph_traversed": False,
        "physical_to_frame_graph_compared": False,
        "incident_pair_DAG_reentered": False,
        "root_mapping_performed": False,
        "solver_or_ancestor_replay": False,
        "source_constructors_called": 0,
        "source_checkers_called": 0,
        "quadratures": 0,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "source_values_selected": False,
        "common_units_or_equality_asserted": False,
        "M_P1_mixing": False,
        "historical_errors_deleted": False,
    }


def _unresolved(
    repo_root: Path,
    reason: str,
    unresolved_gates: list[str],
    provenance: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "incomplete_source_boundary_unresolved",
        "owner": "G",
        "source_only": True,
        "source_bound_G": False,
        "admissible_for_numeric_use": False,
        "historical_F_pressure_or_gap_m15_equality_proved": False,
        "current_P_source_equality_proved": False,
        "P_minus_F": {
            "exact_additive_decomposition": "(P-G)+(G-F)",
            "P_minus_G_evaluated": False,
            "G_minus_F_evaluated": False,
            "either_term_zero_proved": False,
        },
        "unresolved_gates": list(dict.fromkeys(unresolved_gates)),
        "diagnostic": reason,
        "provenance": dict(provenance or {}),
        "scope_audit": _scope_audit(),
    }


def build_gap_current_radius_source(repo_root: str | Path) -> dict[str, Any]:
    """Return a JSON-serializable, source-only current-radius ``G`` packet.

    Every persisted input is hash-authenticated before JSON/gzip decoding.
    Missing or inconsistent evidence returns an unresolved packet rather than
    fabricating a field, radius, endpoint, equality, unit, or numerical gap.
    """

    root = Path(repo_root)
    archive_path = root / ARCHIVE_RELATIVE_PATH
    receipt_path = root / RECEIPT_RELATIVE_PATH
    previous_path = root / PREVIOUS_RESULT_RELATIVE_PATH
    source_path = root / SOURCE_DEFINITION_RELATIVE_PATH
    checker_path = root / CHECKER_SOURCE_RELATIVE_PATH

    provenance: dict[str, Any] = {
        "archive_expected_sha256": EXPECTED_ARCHIVE_SHA256,
        "receipt_expected_sha256": EXPECTED_RECEIPT_SHA256,
        "previous_result_expected_sha256": EXPECTED_PREVIOUS_RESULT_SHA256,
        "archive_key": "pulse_radial_pressure_C3_bindings",
        "owner": "G",
    }
    try:
        archive, archive_meta = _read_gzip_json_after_raw_auth(
            archive_path, EXPECTED_ARCHIVE_SHA256, "current C3 gzip archive"
        )
        receipt, receipt_meta = _read_json_after_raw_auth(
            receipt_path, EXPECTED_RECEIPT_SHA256, "current C3 check receipt"
        )
        previous, previous_meta = _read_json_after_raw_auth(
            previous_path, EXPECTED_PREVIOUS_RESULT_SHA256, "previous saved F/P result"
        )
        archive_family = archive.get("source_family")
        receipt_family = receipt.get("source_family")
        if not isinstance(archive_family, Mapping) or not archive_family or archive_family != receipt_family:
            raise SourceProgramInputError("C3 archive and receipt source families differ")
        frame_path = root / FRAME_ARCHIVE_RELATIVE_PATH
        if not frame_path.is_file() or _sha256(frame_path.read_bytes()) != EXPECTED_FRAME_ARCHIVE_SHA256:
            raise SourceProgramInputError("frame archive raw SHA changed")
        prior_result = _require_mapping(previous.get("result"), "$.result")
        prior_inputs = _require_mapping(prior_result.get("inputs"), "$.result.inputs")
        prior_frame = _require_mapping(prior_inputs.get("frame_archive"), "$.result.inputs.frame_archive")
        receipt_pins = _require_mapping(receipt.get("input_hashes"), "$.input_hashes")
        frame_name = "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
        if prior_frame.get("raw_sha256") != EXPECTED_FRAME_ARCHIVE_SHA256 or receipt_pins.get(frame_name) != EXPECTED_FRAME_ARCHIVE_SHA256:
            raise SourceProgramInputError("C3 receipt and previous F packet do not pin the same frame archive")
        provenance.update(
            {
                "archive": archive_meta,
                "receipt": receipt_meta,
                "previous_result": previous_meta,
                "frame_archive_raw_sha256": EXPECTED_FRAME_ARCHIVE_SHA256,
                "C3_archive_receipt_source_family_equal": True,
                "source_definition": _verify_receipt_source_pin(receipt, source_path),
                "checker_source": _verify_receipt_source_pin(receipt, checker_path),
            }
        )
        receipt_summary = _receipt_summary(receipt)
        binding, gap = _extract_gap_binding(archive)
        g_packet, _ = _extract_gap_owner(binding, gap)
        p_prior = _extract_previous_side(previous, "P")
        f_prior = _extract_previous_side(previous, "F")
        corrections = _extract_previous_corrections(previous)
        comparison = _compare_direct_field_arrays(g_packet["fields"], f_prior["leading_fields_gap"])
    except SourceProgramInputError as exc:
        return _unresolved(root, str(exc), ["saved_input_evidence_complete"], provenance)

    theorem_rows = receipt_summary["theorem_rows_consumed_without_recomputation"]
    row_gate = all(
        len(g_packet["derivative_rows"][field]["all_saved_rows"]) >= 3
        for field in DERIVATIVE_FIELD_NAMES
    )
    m15_gate = bool(g_packet["m15_y_Z0_present"] and g_packet["m15_y_higher_rows_present"] and theorem_rows[
        "pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation"
    ])
    saved_gap_recovery = bool(receipt_summary["saved_gap_Q_recovery_theorem_consumed"])
    zero_fields = all(_node_ids(g_packet["fields"][key]) == [0, 0, 0] for key in ("U", "U_y"))
    source_bound = bool(
        comparison["all_six_direct_field_arrays_equal"]
        and zero_fields
        and row_gate
        and m15_gate
        and saved_gap_recovery
        and receipt_summary["all_passed"]
        and receipt_summary["generic_recovery_FTC_passed"]
    )
    unresolved_gates = []
    if not comparison["all_six_direct_field_arrays_equal"]:
        unresolved_gates.append("saved_Q_m_U_roles_match_previous_F_Z0_to_Z2")
    if not row_gate:
        unresolved_gates.append("gap_derivative_rows_m_y_m_yy_Q_Q_y_Q_yy_Z0_to_Z2")
    if not m15_gate:
        unresolved_gates.append("gap_m15_y_Z0_and_higher_rows_by_commutation_Z1")
    if not saved_gap_recovery:
        unresolved_gates.append("saved_checker_gap_U_m_Q_recovery_rows")
    if not zero_fields:
        unresolved_gates.append("gap_U_and_U_y_source_zero_roles")

    gap_linear_relation = {
        "scope": "current-radius G only",
        "Q_y_equals_negative_Q_under_saved_checker_premises_Z0_to_Z2": bool(source_bound),
        "Q_yy_equals_Q_under_saved_checker_premises_Z0_to_Z2": bool(source_bound),
        "unconditional_source_identity_proved": False,
        "status": (
            "analytic_inference_from_saved_checker_proof_and_linear_formula"
            if source_bound
            else "unproved"
        ),
        "basis": {
            "radial_G_direct_fields_match_previous_F_rows": comparison[
                "all_six_direct_field_arrays_equal"
            ],
            "accepted_checker_gap_recovery_rows_consumed": saved_gap_recovery,
            "ordinary_Z0_to_Z2_only": True,
            "delta_Z_independence_and_source_checker_premises_required": True,
            "new_graph_identity_or_checker_replay": False,
        },
        "historical_F_pressure_or_gap_m15_equality_proved": False,
        "current_P_source_equality_proved": False,
        "source_values_selected": False,
    }
    ordinary_z1 = {
        "G_higher_pressure_FTC_rows_consumed_from_saved_ordinary_commutation": bool(
            theorem_rows[
                "pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation"
            ]
        ),
        "G_to_historical_F_same_source_proved": False,
        "G_to_current_P_same_source_proved": False,
        "new_derivative_or_graph_replay": False,
        "scope": "saved G chart only",
    }
    source_parameters = {
        "G_mu": _optional_direct_tree(binding, "mu", "$.pulse_radial_pressure_C3_bindings.mu"),
        "G_delta": _optional_direct_tree(binding, "delta", "$.pulse_radial_pressure_C3_bindings.delta"),
        "G_coordinate": g_packet["coordinate"],
        "P_prior_Pstar": _optional_direct_tree(
            _require_mapping(p_prior.get("geometry"), "$.result.sides.P.geometry"),
            "Pstar",
            "$.result.sides.P.geometry.Pstar",
        )
        if "geometry" in p_prior
        else {
            "status": "unresolved_missing_P_prior_geometry",
            "source_values_selected": False,
        },
        "common_coordinate_proved": False,
        "common_mu_delta_Pstar_proved": False,
        "same_source_pressure_measure_proved": False,
    }

    return {
        "schema": SCHEMA,
        "status": "source_bound_current_radius_G" if source_bound else "incomplete_source_boundary_unresolved",
        "owner": "G",
        "source_only": True,
        "source_bound_G": source_bound,
        "admissible_for_numeric_use": False,
        "G_is_not_historical_F": True,
        "G_equals_historical_F_pressure_or_gap_m15": False,
        "G_equals_current_P_source": False,
        "historical_F_radius_owner_resolved": False,
        "historical_F_endpoint_owner_resolved": False,
        "receipt_validation": receipt_summary,
        "saved_source_parameters": source_parameters,
        "saved_gap_linear_relation": gap_linear_relation,
        "ordinary_Z1_same_source": ordinary_z1,
        "G_current_radius_source": g_packet,
        "historical_F_prior_source_refs": f_prior,
        "P_prior_source_refs": p_prior,
        "physical_pulse_correction": corrections,
        "direct_field_array_comparison": comparison,
        "residual_decomposition": {
            "P_minus_F": {
                "exact_additive_form": "(P-G)+(G-F)",
                "P_minus_G": {
                    "source_role": "prior P source minus current-radius G source",
                    "evaluated": False,
                    "zero_proved": False,
                },
                "G_minus_F": {
                    "source_role": "current-radius G source minus historical F source",
                    "direct_field_role_match_only": comparison["all_six_direct_field_arrays_equal"],
                    "pressure_or_gap_m15_equality_proved": False,
                    "evaluated": False,
                    "zero_proved": False,
                },
                "numeric_values_computed": False,
                "numeric_difference_evaluated": False,
            }
        },
        "unresolved_gates": unresolved_gates,
        "scope_audit": _scope_audit(),
        "source_values_selected": False,
    }

