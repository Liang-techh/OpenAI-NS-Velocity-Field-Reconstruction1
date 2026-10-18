"""Source-only adoption receipt for the restored pressure join.

This module consumes authenticated saved reports.  It deliberately does not
rehydrate a graph, import a producer, replay an owner, or evaluate a pressure
formula.  The adopted roots are therefore a source-linked contract, while the
old uncontracted packet remains attached verbatim.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any


__all__ = ["adopt", "compute"]


_ENDPOINT_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_post_restore_logu_C1.json"
)
_ENDPOINT_BYTES = 657774
_ENDPOINT_SHA256 = "1d778c10c272acafee88bd42178ad4ae8198d231de1f8b12d185a09da75f8072"
_BACKGROUND_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
)
_BACKGROUND_BYTES = 31700329
_BACKGROUND_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
_PARENT_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
)
_PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
_SIGNED_NATIVE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_signed_native_sources.json"
)
_SIGNED_NATIVE_BYTES = 566530
_SIGNED_NATIVE_SHA256 = "9887a353aeab1b24f301a5b6c240f83f58a5d54d001eaac6d16f932b87979561"
_SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}
_RAW_ROOTS = [186084, 186085]
_ADOPTED_ROOTS = [186096, 186099]
_RM_HANDLES = [183696, 183697]
_TARGETS = [171869, 171870, 178093, 178094]
_PREFIX_NODES = 185704
_PRESSURE_NODES = 186086
_CORRELATED_NODES = 186100
_PRESSURE_PREFIX_SHA256 = "2d6efeb13efe91dd95d4f24fa03eb4626c7685d04e5f29f2bdc6532e101994e3"
_PRESSURE_SHA256 = "dfbba64e1f106fa42c233af28a096f8dfad1e35580163c54cad9eb523d21271e"
_CORRELATED_SHA256 = "83b8594b5b775520cd38c9da1e6633e800ebf170dfb48b6b12eb34f903514fc0"
_RM_SOURCE_IDENTITY = (
    "postrestore pressure join retains -Rm*Qy; the RQy source handle carries "
    "coefficient -1; (Rm*Qy)/Wpost=Qy/Rm"
)


def _check(ok: bool, message: str, issues: list[str]) -> None:
    if not ok:
        issues.append(message)


def _relative(root: Path, name: str) -> Path:
    rel = Path(name)
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError(f"unsafe saved-source path: {name!r}")
    return root / rel


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _file_binding(root: Path, name: str, expected_bytes: int | None, expected_sha256: str) -> dict[str, Any]:
    data = _relative(root, name).read_bytes()
    actual_bytes = len(data)
    actual_sha256 = _digest(data)
    if expected_bytes is not None and actual_bytes != expected_bytes:
        raise ValueError(f"saved-source byte count mismatch for {name}: {actual_bytes} != {expected_bytes}")
    if actual_sha256 != expected_sha256:
        raise ValueError(f"saved-source SHA-256 mismatch for {name}")
    return {"path": name, "bytes": actual_bytes, "sha256": actual_sha256}


def _saved_json(root: Path, name: str, expected_bytes: int, expected_sha256: str) -> tuple[dict[str, Any], dict[str, Any]]:
    binding = _file_binding(root, name, expected_bytes, expected_sha256)
    parsed = json.loads(_relative(root, name).read_bytes().decode("utf-8"))
    if not isinstance(parsed, dict):
        raise ValueError(f"saved report is not an object: {name}")
    return parsed, binding


def _canonical_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return _digest(encoded)


def _coerce_endpoint(value: Mapping[str, Any] | bytes | bytearray | str) -> dict[str, Any]:
    if isinstance(value, Mapping):
        return dict(value)
    if isinstance(value, (bytes, bytearray)):
        value = bytes(value).decode("utf-8")
    if isinstance(value, str):
        parsed = json.loads(value)
        if isinstance(parsed, dict):
            return parsed
    raise TypeError("endpoint_report_raw must be a report object or JSON object bytes")


def _verify_report_ref(root: Path, ref: Mapping[str, Any], label: str) -> dict[str, Any]:
    name = ref.get("path")
    expected_bytes = ref.get("bytes")
    expected_sha256 = ref.get("sha256")
    if not isinstance(name, str) or not isinstance(expected_bytes, int) or not isinstance(expected_sha256, str):
        raise ValueError(f"malformed endpoint input binding: {label}")
    return _file_binding(root, name, expected_bytes, expected_sha256)


def _verify_source_pins(root: Path, packets: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    entries: dict[str, dict[str, Any]] = {}
    owners: dict[str, list[str]] = {}
    for packet_name, packet in packets.items():
        pins = packet.get("source_pins")
        if not isinstance(pins, Mapping):
            raise ValueError(f"missing source pin map: {packet_name}")
        for name, pin in sorted(pins.items()):
            if not isinstance(pin, Mapping) or pin.get("path") != name:
                raise ValueError(f"malformed source pin: {packet_name}:{name}")
            expected = pin.get("expected_sha256")
            if not isinstance(expected, str):
                raise ValueError(f"missing expected source SHA-256: {packet_name}:{name}")
            binding = _file_binding(root, name, None, expected)
            if pin.get("actual_sha256") != binding["sha256"] or pin.get("passed") is not True:
                raise ValueError(f"saved source pin receipt failed: {packet_name}:{name}")
            if pin.get("source_executed") is not False:
                raise ValueError(f"source execution flag is not false: {packet_name}:{name}")
            if name in entries and entries[name]["sha256"] != binding["sha256"]:
                raise ValueError(f"duplicate source pin has different bytes: {name}")
            entries[name] = binding
            owners.setdefault(name, []).append(packet_name)
    for name, binding in entries.items():
        binding["owners"] = sorted(owners[name])
    return {"count": len(entries), "files": [entries[name] for name in sorted(entries)]}


def _validate_headers(packets: Mapping[str, Mapping[str, Any]], issues: list[str]) -> None:
    for label, packet in packets.items():
        _check(packet.get("passed") is True, f"{label} packet is not passed", issues)
        _check(packet.get("source_family") == _SOURCE_FAMILY, f"{label} source family mismatch", issues)
        _check(packet.get("parent_report") == _PARENT_REL, f"{label} parent report mismatch", issues)
        _check(packet.get("parent_report_sha256") == _PARENT_SHA256, f"{label} parent report SHA mismatch", issues)
        _check(packet.get("source_executed") is False, f"{label} source execution was recorded", issues)
        _check(packet.get("accepted_ancestors_replayed") is False, f"{label} replayed an ancestor", issues)


def _validate_graphs(pressure: Mapping[str, Any], correlated: Mapping[str, Any], issues: list[str]) -> None:
    pg = pressure["graph_patch"]
    _check(pg.get("base_node_count") == _PREFIX_NODES, "pressure prefix node count mismatch", issues)
    _check(pg.get("new_node_start") == _PREFIX_NODES, "pressure append start mismatch", issues)
    _check(pg.get("appended_node_count") == 382, "pressure append count is not 382", issues)
    _check(isinstance(pg.get("appended_nodes"), list) and len(pg["appended_nodes"]) == 382, "pressure append payload count mismatch", issues)
    _check(pg.get("full_node_count") == _PRESSURE_NODES, "pressure full node count mismatch", issues)
    _check(pg.get("base_prefix_unchanged") is True, "pressure prefix is not marked unchanged", issues)
    _check(pg.get("base_canonical_nodes_sha256") == _PRESSURE_PREFIX_SHA256, "pressure prefix digest mismatch", issues)
    _check(pg.get("full_canonical_nodes_sha256") == _PRESSURE_SHA256, "pressure graph digest mismatch", issues)

    cg = correlated["graph_patch"]
    _check(cg.get("base_node_count") == _PRESSURE_NODES, "correlated base node count mismatch", issues)
    _check(cg.get("base_prefix_unchanged") is True, "correlated prefix is not marked unchanged", issues)
    _check(isinstance(cg.get("appended_nodes"), list) and len(cg["appended_nodes"]) == 14, "correlated append count is not 14", issues)
    _check(cg.get("full_node_count") == _CORRELATED_NODES, "correlated full node count mismatch", issues)
    _check(cg.get("base_canonical_nodes_sha256") == _PRESSURE_SHA256, "correlated base digest mismatch", issues)
    _check(cg.get("full_canonical_nodes_sha256") == _CORRELATED_SHA256, "correlated graph digest mismatch", issues)


def _validate_pressure(pressure: Mapping[str, Any], issues: list[str]) -> None:
    joined = pressure["joined_pressure_source"]
    _check(joined.get("passed") is True, "raw pressure join receipt failed", issues)
    _check(joined.get("joined_uncontracted") == _RAW_ROOTS, "raw pressure roots changed", issues)
    _check(joined.get("joined_contracted") is None, "raw pressure packet was contracted", issues)
    _check(joined.get("contraction_allowed") is False, "raw pressure contraction gate changed", issues)
    _check(joined.get("contraction_status") == "BLOCKED_UNPROVEN_RESTORATION_POSTRESTORE_SOURCE_CONGRUENCE", "raw pressure block status changed", issues)
    _check(pressure.get("complete_pressure_numeric") is False, "raw packet claims numeric completion", issues)
    _check(pressure.get("global_five_moment_closure") is False, "raw packet claims global closure", issues)
    stable = pressure["stable_chart_handles"]["postrestore"]
    _check(stable.get("upper_RQy_traces") == _RM_HANDLES, "pressure Rm namespace handles changed", issues)


def _validate_correlation(correlated: Mapping[str, Any], issues: list[str]) -> dict[str, str]:
    _check(correlated.get("full_modified_seam_proved") is True, "correlated seam proof is absent", issues)
    _check(correlated.get("composition_conditional") is False, "correlated composition is conditional", issues)
    _check(correlated.get("unresolved_gates") == [], "correlated source has unresolved gates", issues)
    _check(correlated.get("outer_Rm_RQy_retained") == _RM_HANDLES, "correlated Rm handles changed", issues)
    _check(correlated.get("joined_contracted") == _ADOPTED_ROOTS, "adopted contracted roots changed", issues)
    linear = correlated["exact_linearization"]
    _check(linear.get("passed") is True, "exact linearization receipt failed", issues)
    _check(linear.get("targets") == _TARGETS, "linearization target set changed", issues)
    _check(linear.get("only_sum_negative_expanded") is True, "linearization sign receipt changed", issues)
    _check(linear.get("conservative_integer_dependency_check") is True, "linearization integer check failed", issues)
    rows = linear.get("rows")
    _check(isinstance(rows, list) and len(rows) == 2, "linearization row count is not two", issues)
    if isinstance(rows, list) and len(rows) == 2:
        for order, row in enumerate(rows):
            _check(row.get("ordinary_Z_order") == order, f"linearization Z order {order} changed", issues)
            _check(row.get("reduced_handle") == _ADOPTED_ROOTS[order], f"linearization reduced root {order} changed", issues)
            coeffs = row.get("target_coefficients")
            _check(isinstance(coeffs, Mapping) and all(coeffs.get(str(target)) == 0 for target in _TARGETS), f"linearization target coefficients {order} are not zero", issues)
            deps = row.get("surviving_target_dependencies")
            _check(isinstance(deps, Mapping) and all(not deps.get(str(handle)) for handle in row.get("signed_atoms", {})), f"linearization surviving dependencies {order} are nonempty", issues)
            atoms = row.get("signed_atoms", {})
            _check(atoms.get(str(_RM_HANDLES[order])) == -1, f"linearization Rm coefficient {order} changed", issues)
    seams = correlated["seams"]
    expected = {
        "reference_to_restoration": ([171869, 171870], [185885, 185886]),
        "restoration_to_postrestore": ([178093, 178094], [186068, 186069]),
    }
    seam_digests: dict[str, str] = {}
    for name, (left, right) in expected.items():
        seam = seams.get(name, {})
        _check(seam.get("passed") is True, f"{name} seam receipt failed", issues)
        _check(seam.get("finite_N_Qy_identity_proven") is True, f"{name} finite-N Qy identity missing", issues)
        _check(seam.get("physical_Qy_C1_equal") is True, f"{name} physical Qy identity missing", issues)
        _check(seam.get("physical_RQy_C1_equal") is True, f"{name} physical RQy identity missing", issues)
        _check(seam.get("left_RQy_upper") == left, f"{name} left seam handles changed", issues)
        _check(seam.get("right_RQy_lower_wrapper") == right, f"{name} right seam handles changed", issues)
        _check(seam.get("ordinary_Z_orders") == [0, 1], f"{name} Z orders changed", issues)
        _check(seam.get("equality_kind") == "actual common smooth source functions on inherited Z interval", f"{name} equality kind changed", issues)
        seam_digests[name] = _canonical_digest(seam)
    return seam_digests


def _validate_endpoint(endpoint: Mapping[str, Any], pressure: Mapping[str, Any], signed_native: dict[str, Any], issues: list[str]) -> dict[str, Any]:
    _check(endpoint.get("passed") is True, "endpoint packet is not passed", issues)
    _check(endpoint.get("external_Rm_handles") == _RM_HANDLES, "endpoint Rm handles changed", issues)
    _check(endpoint.get("external_Rm_coefficient") == -1, "endpoint Rm coefficient changed", issues)
    boundary = endpoint.get("external_pressure_boundary", {})
    _check(boundary.get("quantity") == "-Rm*Qy", "endpoint boundary quantity changed", issues)
    _check(boundary.get("coefficient") == -1, "endpoint boundary coefficient changed", issues)
    _check(boundary.get("external_source_handles") == _RM_HANDLES, "endpoint boundary handles changed", issues)
    _check(boundary.get("signed_coefficients_by_Z") == {"0": -1, "1": -1}, "endpoint boundary signed coefficients changed", issues)
    _check(boundary.get("source_identity") == _RM_SOURCE_IDENTITY, "endpoint boundary source identity changed", issues)
    _check(signed_native.get("path") == _SIGNED_NATIVE_REL, "endpoint signed-native path changed", issues)
    _check(signed_native.get("bytes") == _SIGNED_NATIVE_BYTES, "endpoint signed-native byte count changed", issues)
    _check(signed_native.get("sha256") == _SIGNED_NATIVE_SHA256, "endpoint signed-native SHA changed", issues)
    _check(endpoint.get("input_reports", {}).get("signed_native_sources", {}).get("sha256") == signed_native["sha256"], "endpoint signed-native alias changed", issues)
    _check(endpoint.get("endpoint_C1", {}).get("passed") is True, "endpoint C1 receipt failed", issues)
    enclosure = endpoint.get("honest_normalized_C1_enclosure", {})
    _check(enclosure.get("directed") is True, "normalized C1 enclosure is not directed", issues)
    _check(enclosure.get("ordinary_Z_orders") == [0, 1], "normalized C1 Z orders changed", issues)
    _check(enclosure.get("enclosure_not_actual_zero_or_value") is True, "normalized C1 enclosure was collapsed", issues)
    _check(pressure["stable_chart_handles"]["postrestore"].get("upper_RQy_traces") == endpoint.get("external_Rm_handles"), "endpoint and pressure Rm namespaces differ", issues)
    return copy.deepcopy(enclosure)


def adopt(root: str | Path, endpoint_report_raw: Mapping[str, Any] | bytes | bytearray | str) -> dict[str, Any]:
    """Return a source-linked pressure adoption receipt.

    The endpoint object is required to be byte-identical to the authenticated
    saved endpoint report.  Semantic namespace failures become an explicit
    conditional receipt; no contracted root is emitted in that case.
    """

    repo = Path(root)
    endpoint, endpoint_binding = _saved_json(repo, _ENDPOINT_REL, _ENDPOINT_BYTES, _ENDPOINT_SHA256)
    supplied = _coerce_endpoint(endpoint_report_raw)
    if supplied != endpoint:
        raise ValueError("endpoint_report_raw does not match the authenticated saved endpoint report")
    background, background_binding = _saved_json(repo, _BACKGROUND_REL, _BACKGROUND_BYTES, _BACKGROUND_SHA256)
    pressure = background["actual_restoration_pressure_source"]
    correlated = background["actual_restoration_correlated_pressure"]
    # The source pin bytes are authenticated before any saved graph metadata is
    # admitted as an adoption receipt.  Archives are only hashed, never read.
    source_bindings = _verify_source_pins(repo, {"pressure": pressure, "correlated": correlated})
    issues: list[str] = []
    _validate_headers({"pressure": pressure, "correlated": correlated}, issues)
    _validate_graphs(pressure, correlated, issues)
    _validate_pressure(pressure, issues)
    seam_digests = _validate_correlation(correlated, issues)

    input_refs = endpoint["actual_post_restore_logu_C1"]["input_reports"]
    endpoint_inputs: dict[str, dict[str, Any]] = {}
    for label, ref in sorted(input_refs.items()):
        endpoint_inputs[label] = _verify_report_ref(repo, ref, label)
    signed_native = endpoint_inputs["signed_native_sources"]
    enclosure = _validate_endpoint(endpoint["actual_post_restore_logu_C1"], pressure, signed_native, issues)
    pin_by_name = {entry["path"]: entry["sha256"] for entry in source_bindings["files"]}
    _check(pressure.get("module_sha256") == pin_by_name.get("src/openai_ns_reconstruction/original_restoration_pressure_source.py"), "pressure module source hash is not pinned", issues)
    _check(pin_by_name.get("src/openai_ns_reconstruction/original_restoration_correlated_pressure.py") == "64f909cd6fc7be559d73b6924e9ac5efd896a4f187127ef843b06caacf75dbcc", "correlated module source hash is not pinned", issues)

    namespace = {
        "prefix": {"node_count": _PREFIX_NODES, "sha256": _PRESSURE_PREFIX_SHA256},
        "pressure_source": {"node_count": _PRESSURE_NODES, "sha256": _PRESSURE_SHA256},
        "correlated_pressure_source": {"node_count": _CORRELATED_NODES, "sha256": _CORRELATED_SHA256},
        "raw_uncontracted_roots": list(_RAW_ROOTS),
        "adopted_contracted_roots": list(_ADOPTED_ROOTS),
        "external_Rm_handles": list(_RM_HANDLES),
        "pressure_handle_path": "stable_chart_handles.postrestore.upper_RQy_traces",
        "correlated_handle_path": "outer_Rm_RQy_retained",
        "endpoint_boundary": {
            "report": _ENDPOINT_REL,
            "handles": list(_RM_HANDLES),
            "source_identity": _RM_SOURCE_IDENTITY,
            "source_graph": {"node_count": _PRESSURE_NODES, "sha256": _PRESSURE_SHA256},
            "handle_path": "actual_post_restore_logu_C1.external_pressure_boundary.external_source_handles",
        },
    }
    endpoint_packet = endpoint["actual_post_restore_logu_C1"]
    boundary = endpoint_packet["external_pressure_boundary"]
    namespace_identity_supported = not issues and boundary.get("external_source_handles") == _RM_HANDLES
    source_adoption = namespace_identity_supported
    status = "SOURCE_ADOPTION_READY_FOR_ROOT_REVIEW" if source_adoption else "CONDITIONAL_NAMESPACE_IDENTITY_UNVERIFIED"
    result: dict[str, Any] = {
        "schema": "original_restoration_pressure_adoption_v1",
        "passed": source_adoption,
        "status": status,
        "source_adoption": source_adoption,
        "namespace_identity_supported": namespace_identity_supported,
        "conditional_reasons": list(issues),
        "source_execution": False,
        "accepted_ancestors_replayed": False,
        "full_joined_pressure_numeric": False,
        "global_five_moment_closure": False,
        "complete_pressure_numeric": False,
        "full_goal_complete": False,
        "raw_uncontracted_pressure": {
            "roots": list(_RAW_ROOTS),
            "joined_contracted": None,
            "status": "BLOCKED_UNPROVEN_RESTORATION_POSTRESTORE_SOURCE_CONGRUENCE",
            "source_report": _BACKGROUND_REL,
            "source_report_sha256": background_binding["sha256"],
        },
        "adopted_pressure_join": {
            "roots": list(_ADOPTED_ROOTS) if source_adoption else None,
            "ordinary_Z_orders": [0, 1],
            "source_report": _BACKGROUND_REL,
            "source_namespace": namespace,
            "target_coefficients_zero": source_adoption,
            "external_Rm_retained": list(_RM_HANDLES),
        },
        "namespace": namespace,
        "seam_identity_digests": seam_digests,
        "external_Rm_normalized_C1": {
            "status": "BOUND",
            "source_report": _ENDPOINT_REL,
            "source_handles": list(_RM_HANDLES),
            "quantity": "-Rm*Qy",
            "coefficient": -1,
            "physical_scale": "Wpost=Rm^2",
            "normalization": "(Rm*Qy)/Wpost=Qy/Rm",
            "enclosure": enclosure,
        },
        "source_bindings": {
            "endpoint_report": endpoint_binding,
            "background_report": background_binding,
            "endpoint_inputs": endpoint_inputs,
            "source_pins": source_bindings,
        },
    }
    module_path = repo / "src/openai_ns_reconstruction/original_restoration_pressure_adoption.py"
    result["api_path"] = "src/openai_ns_reconstruction/original_restoration_pressure_adoption.py"
    result["api_sha256"] = _digest(module_path.read_bytes())
    return result


def compute(root: str | Path, endpoint_report_raw: Mapping[str, Any] | bytes | bytearray | str) -> dict[str, Any]:
    """Alias for :func:`adopt`, retained as the public compute entry point."""

    return adopt(root, endpoint_report_raw)
