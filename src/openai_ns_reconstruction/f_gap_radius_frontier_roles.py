"""Bind the saved radius-spine frontiers to saved source roles.

This leaf consumes only authenticated, already-produced JSON result packets.
It never opens an archive, walks a graph, indexes a protected node, imports a
source constructor, evaluates a source value, or calls an ancestor consumer.

The two correlations emitted here are deliberately role correlations:

* the opaque ``logRp`` child of the saved F spine is correlated with Gamma's
  saved imported and schedule ``logRp`` references;
* the opaque ``mu`` denominator of the saved F spine is correlated with the
  saved F ``geometry.mu`` reference, while the saved current G ``mu`` role is
  reported separately.

Equal opaque metadata is not function equality.  In particular, the protected
frontier labelled ``17`` is never dereferenced or used as a numerical value.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "f_gap_radius_frontier_roles_v1"
COMPONENT = SCHEMA

SPINE_RESULT = Path("tmp/f_gap_radius_spine_binding_v1/result.json")
HISTORICAL_RESULT = Path("tmp/historical_frame_owner_binding_v1/result.json")
GAMMA_RESULT = Path("tmp/frame_axis_gamma_roles_v1/result.json")
F_RESULT = Path("tmp/pulse_gap_current_radius_source_v1/result.json")

# The fifth packet is an optional, separately-authenticated conditional source
# witness.  It is not required for the four-result role binding.  When a
# caller's frozen manifest includes it, it is consumed only after its raw
# bytes have passed the fixed hash below.
OPTIONAL_LOGRP_BRIDGE_RESULT = Path(
    "tmp/original_native_logRp_source_bridge_v1/result.json"
)

SPINE_RESULT_SHA256 = (
    "3bbf29e3479632f5140c34931b9f752250ae770a9a1c2093819465c5e822854d"
)
HISTORICAL_RESULT_SHA256 = (
    "8c1acc7746ea6044bfcacc1ca8aa29dd97db83076e795f4af68e6779530b880d"
)
GAMMA_RESULT_SHA256 = (
    "2facb6d999cbbbb386cb8e53945cee571a2a7ab13688bd3c87507aa6cbb67a91"
)
F_RESULT_SHA256 = (
    "c6d170f6935e0dad3937e155fbbb7e67c6590ee19369291594d10479883fef05"
)
OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256 = (
    "41df6f56fc4a071de606372261a01425c445cd98962c76dbb46164e9699b621c"
)

REQUIRED_INPUTS = {
    SPINE_RESULT.as_posix(): SPINE_RESULT_SHA256,
    HISTORICAL_RESULT.as_posix(): HISTORICAL_RESULT_SHA256,
    GAMMA_RESULT.as_posix(): GAMMA_RESULT_SHA256,
    F_RESULT.as_posix(): F_RESULT_SHA256,
}
OPTIONAL_INPUTS = {
    OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix(): OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256,
}

SOURCE_PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_"
F_NAMESPACE = "$.Rp_native_C5_pulse_source_bindings"
G_NAMESPACE = "$.pulse_radial_pressure_C3_bindings"
GAMMA_NAMESPACE = "Gamma.function_bindings"
LOGRP_FRONTIER_LABEL = 31161
MU_FRONTIER_LABEL = 17


class FGapRadiusFrontierRolesInputError(ValueError):
    """A frozen result manifest or role packet is malformed or changed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FGapRadiusFrontierRolesInputError(message)


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canonical(value: Any) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _canonical_sha256(value: Any) -> str:
    return _sha256(_canonical(value))


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    value = value.lower()
    _require(
        all(character in "0123456789abcdef" for character in value),
        label + " SHA malformed",
    )
    return value


def _manifest_path(value: Any, root: Path) -> str:
    """Normalize a manifest path without permitting an escape from ``root``."""

    _require(isinstance(value, str) and value, "manifest path missing")
    text = value.replace("\\", "/")
    candidate = Path(text)
    # The pulse source checkpoint uses absolute Windows paths.  Convert those
    # to a repository-relative literal before the raw-byte pass.
    if candidate.is_absolute() or (
        len(text) >= 3 and text[1] == ":" and text[2] == "/"
    ):
        try:
            relative = candidate.resolve().relative_to(root.resolve())
        except ValueError as exc:
            raise FGapRadiusFrontierRolesInputError(
                "manifest path is outside repository: " + text
            ) from exc
        text = relative.as_posix()
    relative = Path(text)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes repository: " + text,
    )
    return relative.as_posix()


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes repository: " + name,
    )
    resolved = (root / relative).resolve()
    _require(
        resolved.is_relative_to(root.resolve()),
        "manifest path escapes repository: " + name,
    )
    return resolved


def _manifest_hashes(
    source: Mapping[str, Any],
    hashes: Mapping[str, str] | None,
    root: Path,
) -> dict[str, str]:
    entries: dict[str, str] = {}

    def add(name: Any, digest: Any, label: str) -> None:
        key = _manifest_path(name, root)
        value = _digest(digest, label + "." + key)
        prior = entries.get(key)
        _require(
            prior is None or prior == value,
            "manifest hash conflict: " + key,
        )
        entries[key] = value

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(
                isinstance(item, Mapping),
                "manifest files entry malformed: " + str(index),
            )
            add(item.get("path"), item.get("sha256"), "manifest.files[" + str(index) + "]")

    for field in (
        "sha256",
        "input_hashes",
        "input_sha256",
        "code_sha256",
        "hashes",
    ):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, digest in table.items():
                add(name, digest, "manifest." + field)

    if isinstance(hashes, Mapping):
        for name, digest in hashes.items():
            add(name, digest, "caller hashes")

    _require(entries, "caller frozen manifest has no raw-byte hashes")
    return entries


def _find_hash(entries: Mapping[str, str], name: str) -> str | None:
    if name in entries:
        return entries[name]
    candidates = [
        digest
        for key, digest in entries.items()
        if key.replace("\\", "/").endswith("/" + name.replace("\\", "/"))
    ]
    if len(candidates) == 1:
        return candidates[0]
    return None


def _authenticate(
    root: Path,
    source: Mapping[str, Any],
    hashes: Mapping[str, str] | None,
) -> tuple[dict[str, bytes], dict[str, str], dict[str, str], int]:
    """Authenticate every manifest byte before decoding any result JSON."""

    entries = _manifest_hashes(source, hashes, root)
    for name, expected in REQUIRED_INPUTS.items():
        _require(
            _find_hash(entries, name) == expected,
            "required result SHA changed or is missing: " + name,
        )
    optional_pinned = _find_hash(entries, OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix())
    if optional_pinned is not None:
        _require(
            optional_pinned == OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256,
            "optional logRp bridge SHA changed",
        )

    raw_by_name: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    # No JSON decode is reachable from this loop.  Every entry, including
    # source/checker files that the caller may include, is hashed first.
    for name, expected in entries.items():
        data = _safe(root, name).read_bytes()
        actual_sha = _sha256(data)
        _require(actual_sha == expected, "raw SHA changed: " + name)
        raw_by_name[name] = data
        actual[name] = actual_sha

    # A suffix alias can satisfy the required pin while retaining a different
    # manifest spelling.  Canonicalize the required literal path only after the
    # complete manifest pass and before parsing.
    for name, expected in REQUIRED_INPUTS.items():
        if name not in raw_by_name:
            data = _safe(root, name).read_bytes()
            actual_sha = _sha256(data)
            _require(actual_sha == expected, "required raw SHA changed: " + name)
            raw_by_name[name] = data
            actual[name] = actual_sha

    if optional_pinned is not None and OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix() not in raw_by_name:
        data = _safe(root, OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix()).read_bytes()
        actual_sha = _sha256(data)
        _require(
            actual_sha == OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256,
            "optional logRp bridge raw SHA changed",
        )
        raw_by_name[OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix()] = data
        actual[OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix()] = actual_sha

    return raw_by_name, actual, entries, len(entries)


def _json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FGapRadiusFrontierRolesInputError("invalid JSON: " + label) from exc
    _require(isinstance(value, Mapping), label + " is not an object")
    return value


def _prior_body(
    raw: bytes,
    label: str,
    schema: str | None,
) -> Mapping[str, Any]:
    packet = _json(raw, label)
    _require(packet.get("status") == "returned", label + " wrapper status changed")
    body = packet.get("result")
    _require(isinstance(body, Mapping), label + " result missing")
    if schema is not None:
        _require(body.get("schema") == schema, label + " result schema changed")
    return body


def _mapping(value: Any) -> Mapping[str, Any] | None:
    return value if isinstance(value, Mapping) else None


def _unresolved(role: str, reason: str) -> dict[str, Any]:
    return {
        "status": "unresolved_" + reason,
        "source_role": role,
        "source_node_id": None,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _saved_ref(
    value: Any,
    role: str,
    expected_source_role: str | None = None,
) -> dict[str, Any]:
    """Copy a saved opaque reference without interpreting its handle."""

    if type(value) is int:
        return {
            "status": "opaque_saved_role_ref",
            "source_role": role,
            "source_node_id": value,
            "node_id_is_not_a_numeric_value": True,
            "source_value_selected": False,
        }
    if not isinstance(value, Mapping):
        return _unresolved(role, "missing_or_malformed_role_ref")
    node = value.get("source_node_id")
    if type(node) is not int:
        return _unresolved(role, "missing_or_malformed_role_ref")
    saved_role = value.get("source_role")
    if expected_source_role is not None and saved_role != expected_source_role:
        return _unresolved(role, "saved_role_path_changed")
    return {
        "status": value.get("status") if isinstance(value.get("status"), str) else "direct_saved_role_ref",
        "source_role": saved_role if isinstance(saved_role, str) else role,
        "source_node_id": node,
        "node_id_is_not_a_numeric_value": True,
        "source_value_selected": False,
    }


def _same_opaque_metadata(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Compare saved handle metadata only; never dereference either handle."""

    left_id = left.get("source_node_id")
    right_id = right.get("source_node_id")
    return (
        type(left_id) is int
        and type(right_id) is int
        and left.get("node_id_is_not_a_numeric_value") is True
        and right.get("node_id_is_not_a_numeric_value") is True
        and left.get("source_value_selected") is False
        and right.get("source_value_selected") is False
        and left_id == right_id
    )


def _expected_frontier(
    reference: Mapping[str, Any], expected: int, role: str
) -> dict[str, Any]:
    node = reference.get("source_node_id")
    if type(node) is not int:
        return dict(reference)
    if node != expected:
        return _unresolved(role, "frontier_label_changed")
    return dict(reference)


def _nested_ref(
    container: Mapping[str, Any] | None,
    key: str,
    role: str,
    expected_source_role: str | None = None,
) -> dict[str, Any]:
    if container is None or key not in container:
        return _unresolved(role, "missing_role_field")
    return _saved_ref(container.get(key), role, expected_source_role)


def _source_family_view(
    spine: Mapping[str, Any],
    historical: Mapping[str, Any],
    gamma: Mapping[str, Any],
) -> dict[str, Any]:
    bindings = _mapping(spine.get("source_bindings"))
    f_family = _mapping(bindings.get("F_source_family")) if bindings else None
    old_family = _mapping(bindings.get("old_source_family")) if bindings else None
    receipts = _mapping(historical.get("receipts"))

    receipt_families: dict[str, Any] = {}
    if receipts is not None:
        for owner in ("preRp", "Gamma", "frame"):
            packet = _mapping(receipts.get(owner))
            family = _mapping(packet.get("source_family")) if packet else None
            receipt_families[owner] = dict(family) if family is not None else {
                "status": "unresolved_missing_source_family"
            }

    family_values = [
        value for value in receipt_families.values() if isinstance(value, Mapping) and "status" not in value
    ]
    receipts_equal = len(family_values) == 3 and all(
        value == family_values[0] for value in family_values
    )
    historical_family_sha = (
        _canonical_sha256(family_values[0]) if family_values and receipts_equal else None
    )
    f_metadata_sha = f_family.get("metadata_sha256") if f_family else None
    old_metadata_sha = old_family.get("metadata_sha256") if old_family else None
    spine_equal = (
        isinstance(f_metadata_sha, str)
        and isinstance(old_metadata_sha, str)
        and f_metadata_sha == old_metadata_sha
    )
    spine_historical_equal = (
        spine_equal
        and isinstance(historical_family_sha, str)
        and f_metadata_sha == historical_family_sha
    )
    return {
        "prefix": SOURCE_PREFIX,
        "spine_F": dict(f_family) if f_family is not None else {
            "status": "unresolved_missing_F_source_family"
        },
        "spine_old": dict(old_family) if old_family is not None else {
            "status": "unresolved_missing_old_source_family"
        },
        "historical_receipts": receipt_families,
        "metadata_sha256": {
            "spine_F": f_metadata_sha,
            "spine_old": old_metadata_sha,
            "historical_receipts": historical_family_sha,
        },
        "correlations": {
            "spine_F_equals_old_metadata": spine_equal,
            "spine_F_equals_historical_receipts": spine_historical_equal,
            "frame_axis_source_family_matches": gamma.get("source_family_matches") is True,
            "source_family_is_not_function_identity": True,
        },
        "status": (
            "saved_source_family_metadata_aligned"
            if spine_historical_equal and gamma.get("source_family_matches") is True
            else "saved_source_family_metadata_partially_resolved"
        ),
    }


def _owners_and_namespaces(
    spine: Mapping[str, Any],
    historical: Mapping[str, Any],
) -> dict[str, Any]:
    bindings = _mapping(spine.get("source_bindings"))
    archives = _mapping(historical.get("archive_paths"))
    return {
        "prefix": SOURCE_PREFIX,
        "namespaces": {
            "F": F_NAMESPACE,
            "G": G_NAMESPACE,
            "Gamma": GAMMA_NAMESPACE,
        },
        "archives": {
            "F": bindings.get("F_archive") if bindings else None,
            "Gamma": archives.get("Gamma") if archives else None,
            "frame": archives.get("frame") if archives else None,
        },
        "role_owners": {
            "logRp": {
                "spine_owner": "F.log_radius.other_child",
                "saved_owner": "Gamma imported/schedule logRp",
                "Gamma_imported_role": "Gamma.function_bindings.native_binding.imported_logRp",
                "Gamma_schedule_role": "Gamma.function_bindings.schedule.logRp",
            },
            "mu": {
                "spine_owner": "F.positive_quotient.denominator",
                "F_geometry_owner": F_NAMESPACE + ".geometry.mu",
                "G_owner": G_NAMESPACE + ".mu",
                "owner_distinction_is_preserved": True,
            },
        },
    }


def _formal_candidate() -> dict[str, Any]:
    return {
        "R_Z": {
            "expression": "R_Z = R * ((log Rp)_Z + xi_Z/mu - xi*mu_Z/mu^2)",
            "status": "formal_candidate_referenced",
            "source_value_selected": False,
            "numeric_evaluation": False,
            "source_function_identity_required": True,
        },
        "candidate_is_not_an_actual_F_formula": True,
    }


def _bridge_witness(
    body: Mapping[str, Any], spine_log: Mapping[str, Any]
) -> dict[str, Any]:
    proof = _mapping(body.get("proof"))
    gamma_witness = _mapping(body.get("Gamma_source_witness"))
    imported = _saved_ref(
        gamma_witness.get("imported_logRp") if gamma_witness else None,
        "original_native_logRp_source_bridge.Gamma_source_witness.imported_logRp",
    )
    schedule = _saved_ref(
        gamma_witness.get("schedule_logRp") if gamma_witness else None,
        "original_native_logRp_source_bridge.Gamma_source_witness.schedule_logRp",
    )
    old = _saved_ref(
        proof.get("old_logRp_node") if proof else None,
        "original_native_logRp_source_bridge.proof.old_logRp_node",
    )
    conditional = bool(
        body.get("source_identity_scope") == "absolute-radius log recipe only"
        and body.get("source_identity_is_global_alias") is False
        and proof is not None
        and proof.get("conditional") is True
        and proof.get("numeric_evaluation") is False
    )
    return {
        "status": "conditional_saved_source_witness" if conditional else "unresolved_bridge_scope",
        "source_identity_scope": body.get("source_identity_scope"),
        "conditional": proof.get("conditional") if proof else False,
        "numeric_evaluation": proof.get("numeric_evaluation") if proof else False,
        "old_logRp": old,
        "Gamma_imported_logRp": imported,
        "Gamma_schedule_logRp": schedule,
        "same_spine_to_old_opaque_metadata": _same_opaque_metadata(spine_log, old),
        "same_spine_to_Gamma_imported_opaque_metadata": _same_opaque_metadata(spine_log, imported),
        "role_reference_correspondence": True,
        "scope_independent_of_role_identity": True,
        "not_a_role_join": True,
        "role_join_proved": False,
        "function_identity_proved": False,
        "capture_free": False,
        "F_radius_owner_bound": False,
        "source_values_selected": False,
        "numeric_use_admissible": False,
    }


def build_frontier_roles(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Authenticate frozen result bytes and emit bounded role correlations."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    repo = Path(root).resolve()
    raw, actual, entries, manifest_entry_count = _authenticate(repo, source, hashes)

    # All four required JSON packets are decoded only after the complete
    # manifest raw-byte pass above.  No gzip or source archive is involved.
    spine_body = _prior_body(
        raw[SPINE_RESULT.as_posix()],
        SPINE_RESULT.as_posix(),
        None,
    )
    historical = _prior_body(
        raw[HISTORICAL_RESULT.as_posix()],
        HISTORICAL_RESULT.as_posix(),
        "historical_frame_owner_binding_v1",
    )
    gamma = _prior_body(
        raw[GAMMA_RESULT.as_posix()],
        GAMMA_RESULT.as_posix(),
        "frame_axis_gamma_roles_v1",
    )
    f_body = _prior_body(
        raw[F_RESULT.as_posix()],
        F_RESULT.as_posix(),
        "pulse_gap_current_radius_source_v1",
    )

    spine = _mapping(spine_body.get("radius_spine_binding")) or {}
    _require(
        spine.get("schema") == "f_gap_radius_spine_binding_v1",
        SPINE_RESULT.as_posix() + " nested result schema changed",
    )
    formula = _mapping(spine.get("formula_structure"))
    f_formula = _mapping(formula.get("F")) if formula else None
    spine_log = _expected_frontier(
        _nested_ref(
            f_formula,
            "log_radius_other_child",
            "F.log_radius.other_child",
            "F.log_radius.other_child",
        ),
        LOGRP_FRONTIER_LABEL,
        "F.log_radius.other_child",
    )
    spine_mu = _expected_frontier(
        _nested_ref(
            f_formula,
            "scale_denominator",
            "F.positive_quotient.denominator",
            "F.positive_quotient.denominator",
        ),
        MU_FRONTIER_LABEL,
        "F.positive_quotient.denominator",
    )

    historical_numbers = _mapping(historical.get("saved_numbers"))
    hist_imported = _expected_frontier(
        _nested_ref(
            historical_numbers,
            "Gamma_function_bindings_native_binding_imported_logRp",
            "Gamma.function_bindings.native_binding.imported_logRp",
        ),
        LOGRP_FRONTIER_LABEL,
        "Gamma.function_bindings.native_binding.imported_logRp",
    )
    hist_schedule = _expected_frontier(
        _nested_ref(
            historical_numbers,
            "Gamma_function_bindings_schedule_logRp",
            "Gamma.function_bindings.schedule.logRp",
        ),
        LOGRP_FRONTIER_LABEL,
        "Gamma.function_bindings.schedule.logRp",
    )
    gamma_native = _mapping(gamma.get("gamma_native_binding"))
    gamma_imported = _expected_frontier(
        _nested_ref(
            gamma_native,
            "imported_logRp",
            "Gamma.function_bindings.native_binding.imported_logRp",
        ),
        LOGRP_FRONTIER_LABEL,
        "Gamma.function_bindings.native_binding.imported_logRp",
    )

    f_prior = _mapping(f_body.get("historical_F_prior_source_refs"))
    f_geometry = _mapping(f_prior.get("geometry")) if f_prior else None
    f_mu = _expected_frontier(
        _nested_ref(
            f_geometry,
            "mu",
            F_NAMESPACE + ".geometry.mu",
            F_NAMESPACE + ".geometry.mu",
        ),
        MU_FRONTIER_LABEL,
        F_NAMESPACE + ".geometry.mu",
    )
    saved_parameters = _mapping(f_body.get("saved_source_parameters"))
    g_mu = _expected_frontier(
        _nested_ref(
            saved_parameters,
            "G_mu",
            G_NAMESPACE + ".mu",
        ),
        MU_FRONTIER_LABEL,
        G_NAMESPACE + ".mu",
    )

    log_checks = {
        "historical_check_imported_equals_schedule": (
            _mapping(historical.get("checks")) or {}
        ).get("Gamma_imported_logRp_matches_schedule_logRp") is True,
        "historical_binding_imported_equals_schedule": (
            _mapping(historical.get("binding")) or {}
        ).get("same_Gamma_imported_schedule_node") is True,
        "spine_to_Gamma_imported_opaque_metadata": _same_opaque_metadata(spine_log, hist_imported),
        "Gamma_imported_to_schedule_opaque_metadata": _same_opaque_metadata(hist_imported, hist_schedule),
        "Gamma_result_imported_matches_historical_opaque_metadata": _same_opaque_metadata(
            hist_imported, gamma_imported
        ),
    }
    log_complete = all(log_checks.values())

    mu_checks = {
        "spine_to_F_geometry_mu_opaque_metadata": _same_opaque_metadata(spine_mu, f_mu),
        "F_geometry_mu_to_G_mu_opaque_metadata": _same_opaque_metadata(f_mu, g_mu),
        "F_and_G_roles_are_distinct": f_mu.get("source_role") != g_mu.get("source_role"),
        "F_prior_side_is_F": f_prior is not None and f_prior.get("side") == "F",
        "historical_F_radius_owner_still_unresolved": f_body.get("historical_F_radius_owner_resolved") is False,
        "G_is_not_historical_F": f_body.get("G_is_not_historical_F") is True,
    }
    mu_complete = all(mu_checks.values())

    bridge_pinned = OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix() in raw
    if bridge_pinned:
        bridge_body = _prior_body(
            raw[OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix()],
            OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix(),
            "original_native_logRp_source_bridge_v1_unentered",
        )
        bridge = _bridge_witness(bridge_body, spine_log)
    else:
        bridge = {
            "status": "unresolved_not_in_frozen_manifest",
            "role_reference_correspondence": False,
            "scope_independent_of_role_identity": True,
            "not_a_role_join": True,
            "role_join_proved": False,
            "function_identity_proved": False,
            "capture_free": False,
            "F_radius_owner_bound": False,
            "source_values_selected": False,
            "numeric_use_admissible": False,
        }

    scope_flags = {
        "source_value_selected": False,
        "source_values_computed": False,
        "source_graph_traversed": False,
        "protected_root_indexed": False,
        "protected_frontier_dereferenced": False,
        "capture_free": False,
        "capture_free_entire_expression": False,
        "F_radius_owner_bound": False,
        "ordinary_Z_transport_proved": False,
        "m15_closure_proved": False,
        "m15_identity_proved": False,
        "old_errors_removed": False,
        "old_error_terms_removed": False,
        "error_deletion_proved": False,
        "error_transport_certified": False,
        "function_identity_proved": False,
        "same_source_function_proved": False,
        "numeric_use_admissible": False,
        "integrals_evaluated": False,
    }

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": (
            "saved_role_correspondence_complete"
            if log_complete and mu_complete
            else "saved_role_correspondence_unresolved"
        ),
        "source_only": True,
        "role_reference_correspondence_is_not_function_equivalence": True,
        "role_reference_correspondence": {
            "is_function_equivalence": False,
            "is_source_value_binding": False,
            "logRp": {
                "frontier": spine_log,
                "Gamma_imported": hist_imported,
                "Gamma_schedule": hist_schedule,
                "Gamma_result_imported": gamma_imported,
                "checks": log_checks,
                "complete": log_complete,
            },
            "mu": {
                "frontier": spine_mu,
                "F_geometry_mu": f_mu,
                "G_mu": g_mu,
                "checks": mu_checks,
                "complete": mu_complete,
            },
        },
        "source_family": _source_family_view(spine, historical, gamma),
        "owners_and_namespaces": _owners_and_namespaces(spine, historical),
        "protected_frontier_17": {
            "label": MU_FRONTIER_LABEL,
            "comparison": "saved opaque metadata only",
            "indexed": False,
            "dereferenced": False,
            "used_as_numeric_value": False,
            "source_value_selected": False,
        },
        "optional_logRp_source_witness": bridge,
        "formal_candidate": _formal_candidate(),
        "scope_flags": scope_flags,
        "raw_hash_authentication": {
            "all_manifest_raw_hashes_checked_before_result_json_read": True,
            "manifest_entry_count": manifest_entry_count,
            "authenticated_paths": sorted(actual),
            "required_results": {
                name: {
                    "path": name,
                    "expected_sha256": expected,
                    "actual_sha256": actual.get(name),
                    "json_decoded_after_raw_authentication": True,
                }
                for name, expected in REQUIRED_INPUTS.items()
            },
            "optional_logRp_bridge": {
                "path": OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix(),
                "expected_sha256": OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256,
                "actual_sha256": actual.get(OPTIONAL_LOGRP_BRIDGE_RESULT.as_posix()),
                "pinned": bridge_pinned,
            },
        },
        "consumed_result_schemas": {
            SPINE_RESULT.as_posix(): "f_gap_radius_spine_binding_v1_result_v1",
            HISTORICAL_RESULT.as_posix(): "historical_frame_owner_binding_v1",
            GAMMA_RESULT.as_posix(): "frame_axis_gamma_roles_v1",
            F_RESULT.as_posix(): "pulse_gap_current_radius_source_v1",
        },
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Run four in-memory guard checks without reading repository input."""

    log_spine = _saved_ref(31161, "F.log_radius.other_child")
    log_imported = _saved_ref(
        31161,
        "Gamma.function_bindings.native_binding.imported_logRp",
    )
    log_schedule = _saved_ref(
        31161,
        "Gamma.function_bindings.schedule.logRp",
    )
    f_mu = _saved_ref(17, F_NAMESPACE + ".geometry.mu")
    g_mu = _saved_ref(17, G_NAMESPACE + ".mu")
    missing = _unresolved(F_NAMESPACE + ".geometry.mu", "missing_role_field")
    checks = {
        "logRp_imported_schedule_correlate_as_opaque_metadata": (
            _same_opaque_metadata(log_spine, log_imported)
            and _same_opaque_metadata(log_imported, log_schedule)
        ),
        "protected_mu_is_not_indexed_or_dereferenced": (
            _same_opaque_metadata(f_mu, g_mu)
            and f_mu.get("node_id_is_not_a_numeric_value") is True
            and f_mu.get("source_value_selected") is False
        ),
        "F_and_G_mu_owner_roles_remain_distinct": (
            f_mu.get("source_role") != g_mu.get("source_role")
        ),
        "missing_role_stays_unresolved": (
            missing.get("status") == "unresolved_missing_role_field"
            and missing.get("source_value_selected") is False
        ),
    }
    return {
        "schema": SCHEMA + "_synthetic_fixture_v1",
        "component": COMPONENT,
        "status": "synthetic_guard_checks_passed" if all(checks.values()) else "synthetic_guard_checks_failed",
        "synthetic_only": True,
        "actual_source_read": False,
        "source_graph_traversed": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_use_admissible": False,
        "protected_root_indexed": False,
        "checks": checks,
        "checks_passed": bool(all(checks.values())),
        "role_reference_correspondence_is_not_function_equivalence": True,
    }


__all__ = [
    "COMPONENT",
    "F_NAMESPACE",
    "F_RESULT",
    "F_RESULT_SHA256",
    "FGapRadiusFrontierRolesInputError",
    "GAMMA_RESULT",
    "GAMMA_RESULT_SHA256",
    "G_NAMESPACE",
    "HISTORICAL_RESULT",
    "HISTORICAL_RESULT_SHA256",
    "LOGRP_FRONTIER_LABEL",
    "MU_FRONTIER_LABEL",
    "OPTIONAL_INPUTS",
    "OPTIONAL_LOGRP_BRIDGE_RESULT",
    "OPTIONAL_LOGRP_BRIDGE_RESULT_SHA256",
    "REQUIRED_INPUTS",
    "SCHEMA",
    "SPINE_RESULT",
    "SPINE_RESULT_SHA256",
    "build_frontier_roles",
    "verify_synthetic_fixture",
]
