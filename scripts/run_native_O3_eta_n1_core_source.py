"""Durable, unexecuted manager for the native O3 eta n1 core source leaf.

This file is intentionally a manager rather than a mathematical implementation.
``prepare`` authenticates frozen bytes, literals, hashes and AST contracts.  A
later root-only ``run`` records entry before the one import and one ``derive``
call, journals every strict two-argument callback before retaining it, and
saves the complete nested result before validation.  The kernel and capsule
digests are supplied by the frozen manifest at prepare time; no placeholder
digest is trusted here.
"""

from __future__ import annotations

import argparse
import ast
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import traceback
from typing import Any


REPO = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
WORK = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")

COMPONENT = "native_O3_eta_n1_core_source"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
MANIFEST_SCHEMA = f"{COMPONENT}_manifest_v1"
SOURCE_SCHEMA = f"{COMPONENT}_source_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_checkpoint_v1"
JOURNAL_SCHEMA = f"{COMPONENT}_journal_v1"
CAPSULE_SCHEMA = "native_O3_eta_n1_core_source_source_capsule_v1"

# The parent report was sealed before the current 119d branch commit.  Its
# own source lineage is cddf -> caac; the current leaf is based at 119d.
PARENT_REPORT = REPO / "experiments/root_st073/lei_ren_part1_paper_compliant_current_native_O3_eta_core_moment_functionals.json"
PARENT_REPORT_SHA256 = "32a498d1df35f8708ee28baaf22aa0fc0c203d3124be1a2484f7e728d7cce177"
PARENT_REPORT_BYTES = 48218881
PARENT_REPORT_INTERNAL_PARENT_COMMIT = "cddf394af87d3508f51718690d01a0ee4c4d349f"
PARENT_REPORT_GRANDPARENT_COMMIT = "caac8c0a55aa34dc7334f2a68311fc15ddb20a13"
CURRENT_HEAD_COMMIT = "119d81ff2e2303b90457dbec436859b5cd922944"
ACCEPTED_MATH_COMMIT = CURRENT_HEAD_COMMIT
PARENT_COMMIT = CURRENT_HEAD_COMMIT
PARENT_SOURCE_PIN_COUNT = 2459
PARENT_CHECKPOINT_COUNT = 348
PARENT_FAILED_ARTIFACT_COUNT = 211
PARENT_PROTECTED_FILE_COUNT = 16
PARENT_DIRTY_PATH_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
FROZEN_SOURCE_UNION_COUNT = 2464
PUBLISHED_SOURCE_UNION_COUNT = 2468
PUBLISHED_CHECKPOINT_COUNT = 352
FAILED_COUNT_IF_NEW_FAILURE = 218
PHYSICAL_BASELINE = "8f36f19e3c524b424a9286c46cc8e0b855f75307"

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "5469bb9c995bdb9a9d1d4646e7d09dcced76e69285f0b765a85483f3a25ef1a2"
MAIN_BRIEF_BYTES = 82318

KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_O3_eta_n1_core_source.py")
RUNNER_RELATIVE = Path("scripts/run_native_O3_eta_n1_core_source.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
INPUT_PATH = WORK / f"{COMPONENT}_input_2026_10_28.json"
MANIFEST_PATH = WORK / f"{COMPONENT}_manifest_2026_10_28.json"
BRIEF_PATH = WORK / f"{COMPONENT}_brief_2026_10_28.json"
CAPSULE_PATH = WORK / f"{COMPONENT}_source_capsule_2026_10_28.json"
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_PATH = CHECKPOINT_DIR / "source.json"
ENTERED_PATH = CHECKPOINT_DIR / "entered.json"
JOURNAL_PATH = CHECKPOINT_DIR / "journal.jsonl"
RESULT_PATH = CHECKPOINT_DIR / "result.json"
BRIEF_SOFT_LIMIT = 97_000

# These are deliberately absent until the root freezes the actual kernel,
# manifest, input and capsule.  All runtime pins come from their JSON specs.
KERNEL_SHA256: str | None = None
CAPSULE_SHA256: str | None = None

RHO_GRID = ("4",)
UNIFORM_RHO_GRID = ("0", "4")
ETA_GRID = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
UNIFORM_ETA_GRID = ("-1", "1")
PAYLOAD_RHO_LABELS = ("4", "uniform_rho")
PAYLOAD_ETA_LABELS = (*ETA_GRID, "uniform_eta")
NEW_MATH_FAMILIES = (
    "exact_view_registration",
    "mixed_core_jets",
    "gauge_third_jet",
    "n1_differential_operators",
    "canonical_rhs_assembly",
)
NEXT_SUCCESS_TASK = "source_native_O3_eta_n1_core_response_and_m15_binding"
GLOBAL_FALSE_KEYS = (
    "source_ftc_zero_identity",
    "source_C5_installation_claimed",
    "new_loop_domain_admission_claimed",
    "admissible_cone_admission_claimed",
    "cap_to_source_identity_proved",
    "quantitative_packet_error_reduction_proved",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "physical_layer_completed",
    "global_five_moment_closure",
    "time_order_recursion_claimed",
    "old_operators_replayed",
    "old_pointwise_S_recovered",
    "old_mass_relation_proved",
    "old212_adapter_identity",
    "old212_quantclosure",
    "native_old_density_identity",
    "native_old_density_closure",
    "old_errors_removed",
    "old_error_terms_removed",
    "full_goal_complete",
    "actual_C3_C5_incoming_source_identity_proved",
    "actual_C3_C5_control_source_identity_proved",
    "actual_C3_C5_profile_identity_proved",
    "actual_C3_cap_to_C5_source_bound_proved",
    "compact_volume_cap_same_source_bound_proved",
)

# Preserve this exact seven-rule policy in every checkpoint and metadata brief.
EXECUTION_RULES = [
    "Never rerun any entered ancestor mathematical/helper family, including failed leaves, old builders/providers/normalizers/checkers, source-value reducers, quadrature or recurrence. Static source/AST/literal JSON/hash inspection is allowed.",
    "Before a NEW mathematical/source-expression consumer: freeze code and input, persist source checkpoint and pending brief, persist entry before invocation, invoke once, save complete result before review. Preserve a failure and repair in a separate frozen leaf; never reenter.",
    "Authenticate raw gzip SHA BEFORE decoding. Never decode docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json or the huge Rm_loop_density_identity.json; raw SHA only.",
    "Never restart stopped session 17273. Retain nonzero incoming histories and ordinary Z1. Bounds, midpoints and interval overlap are not function values; remove old numerical errors only with a same-source packet theorem.",
    "Use current tree and GitHub; preserve all dirty and untracked files. Commit/push only exact accepted scoped files. Goal remains ACTIVE / INCOMPLETE.",
    "Children: luna_max_worker, gpt-5.6-luna/max fixed role, fork_turns:none, complete brief. Runtime model metadata unavailable means null; never claim it was verified.",
    "A worker function that checks mathematical source identities counts as a source consumer even when it consumes static receipts or is called a pure self-check. Workers return unexecuted code; root alone freezes code and inputs, records entry, invokes once and saves the complete result before review. AST parsing and metadata-only audits remain permitted. Do not execute or reconstruct the rejected memory draft result.",
]


def require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(char in "0123456789abcdefABCDEF" for char in value)


def grid_token(value: Any) -> str:
    aliases = {"-0.5": "-1/2", "-0.25": "-1/4", "0.0": "0", "0.25": "1/4", "0.5": "1/2", "1.0": "1", "4.0": "4"}
    text = str(value)
    return aliases.get(text, text)


def grid_tokens(value: Any) -> tuple[str, ...]:
    require(isinstance(value, (list, tuple)), "grid row must be a list")
    return tuple(grid_token(item) for item in value)


def wire(path: str | Path) -> str:
    return Path(path).as_posix()


def encode(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def encode_compact(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def write_new(path: Path, value: bytes | Any) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_replace(path: Path, value: bytes | Any) -> bytes:
    raw = value if isinstance(value, bytes) else encode(value)
    temporary = path.with_name(path.name + f".{COMPONENT}.replace")
    with temporary.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return raw


def candidates(value: str | Path) -> list[Path]:
    path = Path(value)
    if path.is_absolute():
        return [path]
    return list(dict.fromkeys((Path.cwd() / path, WORK / path, REPO / path, REPO / "experiments/root_st073" / path)))


def resolve_artifact(value: str | Path, expected: str | None = None) -> Path:
    pool = candidates(value)
    if expected is not None:
        require(valid_sha(expected), f"expected SHA malformed for {value}")
        for path in pool:
            if path.is_file() and sha_bytes(path.read_bytes()) == expected.lower():
                return path
    for path in pool:
        if path.is_file():
            return path
    raise FileNotFoundError(str(value))


def read_bytes(value: str | Path, expected: str | None = None) -> bytes:
    path = resolve_artifact(value, expected)
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), f"raw SHA changed: {wire(path)}")
    return raw


def read_json(value: str | Path, expected: str | None = None) -> tuple[Path, bytes, dict[str, Any]]:
    path = resolve_artifact(value, expected)
    raw = read_bytes(path, expected)
    parsed = json.loads(raw.decode("utf-8"))
    require(isinstance(parsed, dict), f"JSON object required: {wire(path)}")
    return path, raw, parsed


def artifact(path: str | Path, raw: bytes | None = None, *, embed: bool = False, expected: str | None = None) -> dict[str, Any]:
    resolved = resolve_artifact(path, expected)
    data = resolved.read_bytes() if raw is None else raw
    if expected is not None:
        require(sha_bytes(data) == expected.lower(), f"artifact SHA changed: {wire(resolved)}")
    item = {"path": wire(resolved), "sha256": sha_bytes(data), "bytes": len(data), "raw_bytes_exact": True, "raw_readonly": True}
    if embed:
        import base64
        item["raw_bytes_base64"] = base64.b64encode(data).decode("ascii")
    else:
        item["hash_only"] = True
    return item


def pointer(item: Mapping[str, Any]) -> dict[str, Any]:
    require(isinstance(item, Mapping), "artifact pointer must be a mapping")
    return {key: item[key] for key in ("path", "sha256", "bytes") if key in item}


def parse_ast(path: Path, raw: bytes, label: str) -> dict[str, Any]:
    try:
        tree = ast.parse(raw.decode("utf-8"), filename=str(path), mode="exec")
    except (SyntaxError, UnicodeDecodeError) as exc:
        raise ValueError(f"{label} AST parse failed: {exc}") from exc
    return {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "ast_parse": "ok",
        "node_count": sum(1 for _ in ast.walk(tree)),
        "top_level_definitions": [node.name for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))],
    }


def top_level_literals(raw: bytes, wanted: set[str], label: str) -> dict[str, Any]:
    try:
        tree = ast.parse(raw.decode("utf-8"), filename=label, mode="exec")
    except (SyntaxError, UnicodeDecodeError) as exc:
        raise ValueError(f"{label} AST parse failed: {exc}") from exc
    values: dict[str, Any] = {}
    for node in tree.body:
        targets: list[ast.expr] = []
        if isinstance(node, ast.Assign):
            targets = [target for target in node.targets if isinstance(target, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target]
        for target in targets:
            if target.id not in wanted:
                continue
            try:
                values[target.id] = ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError) as exc:
                raise ValueError(f"{label} literal {target.id} is not static") from exc
    return values


def choose_literal(values: Mapping[str, Any], primary: str, *aliases: str) -> Any:
    present = [name for name in (primary, *aliases) if name in values]
    require(present, f"kernel literal missing {primary}")
    first = values[present[0]]
    for name in present[1:]:
        require(values[name] == first, f"kernel literals disagree: {present[0]} and {name}")
    return first


def ast_kernel_contract(path: Path, raw: bytes) -> dict[str, Any]:
    # Only these assignments are literal AST contracts in the exact kernel.
    # EXPECTED_* are derived with dict(...), and API contains a SCHEMA Name;
    # neither is needed for this manager contract.
    values = top_level_literals(raw, {"SUCCESS_COUNTS", "SUCCESS_FLAGS", "CALLBACK_EVENT_COUNTS", "COUNT_EVENT_MAP", "COMPONENT"}, "kernel")
    counts = choose_literal(values, "SUCCESS_COUNTS")
    flags = choose_literal(values, "SUCCESS_FLAGS")
    callback = choose_literal(values, "CALLBACK_EVENT_COUNTS")
    event_map = values.get("COUNT_EVENT_MAP")
    require(isinstance(counts, Mapping) and isinstance(flags, Mapping) and isinstance(callback, Mapping), "kernel contract mappings missing")
    require(isinstance(event_map, Mapping), "kernel COUNT_EVENT_MAP missing")
    require("derive" in parse_ast(path, raw, "kernel")["top_level_definitions"], "kernel derive API missing")
    return {"SUCCESS_COUNTS": dict(counts), "SUCCESS_FLAGS": dict(flags), "CALLBACK_EVENT_COUNTS": dict(callback), "COUNT_EVENT_MAP": dict(event_map), "API": None, "component_literal": values.get("COMPONENT")}


def authenticate_spec(spec: Mapping[str, Any], label: str, expected_path: Path | None = None) -> tuple[Path, bytes, str]:
    require(isinstance(spec, Mapping), f"{label} spec missing")
    expected = spec.get("sha256")
    require(valid_sha(expected), f"{label} SHA missing or malformed")
    path = resolve_artifact(spec.get("path"), expected)
    raw = read_bytes(path, expected)
    if expected_path is not None:
        require(path.resolve() == expected_path.resolve(), f"{label} path changed")
    if "bytes" in spec:
        require(type(spec["bytes"]) is int and spec["bytes"] == len(raw), f"{label} byte count changed")
    return path, raw, sha_bytes(raw)


def hash_map(value: Any, label: str) -> dict[str, str]:
    require(isinstance(value, Mapping), f"{label} must be a mapping")
    result: dict[str, str] = {}
    for key, value in value.items():
        require(isinstance(key, str) and valid_sha(value), f"{label} contains malformed SHA")
        result[key] = value.lower()
    return result


def authenticate_hash_map(value: Any, label: str) -> dict[str, str]:
    result = hash_map(value, label)
    for path, expected in result.items():
        read_bytes(path, expected)
    return result


def load_parent() -> dict[str, Any]:
    path, raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    require(path.resolve() == PARENT_REPORT.resolve() and len(raw) == PARENT_REPORT_BYTES, "parent report raw pin changed")
    require(report.get("component") == "native_O3_eta_core_moment_functionals" and report.get("passed") is True and report.get("attempt_status") == "passed", "parent report status changed")
    require(report.get("parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT and report.get("prepublication_base_head") == PARENT_REPORT_INTERNAL_PARENT_COMMIT, "parent report internal lineage changed")
    require(report.get("parent_report_internal_parent_commit") == PARENT_REPORT_GRANDPARENT_COMMIT, "parent report grandparent changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT and report.get("source_union_count") == PARENT_SOURCE_PIN_COUNT, "parent source pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT and report.get("preserved_failed_artifact_count") == PARENT_FAILED_ARTIFACT_COUNT, "parent checkpoint count changed")
    require(report.get("parent_protected_file_count") == PARENT_PROTECTED_FILE_COUNT and report.get("parent_dirty_path_count") == PARENT_DIRTY_PATH_COUNT and report.get("global_false_count") == PARENT_FALSE_FLAG_COUNT, "parent caps changed")
    require(report.get("source_only") is True and report.get("metadata_only") is False and report.get("full_goal_complete") is False and report.get("physical_layer_completed") is False, "parent scope changed")
    require(report.get("numeric_evaluation") is True and report.get("current_leaf_numeric_evaluation") is True, "parent evaluation scope changed")
    require(report.get("new_flags") == report.get("expected_new_flags"), "parent flags changed")
    parent_flags = report.get("parent_false_flags")
    require(isinstance(parent_flags, list) and len(parent_flags) == PARENT_FALSE_FLAG_COUNT, "parent false flag list changed")
    require(report.get("parent_false_flag_values") == {name: False for name in parent_flags}, "parent false flag values changed")

    old = report.get("parent_report")
    require(isinstance(old, Mapping) and old.get("sha256") == "be8faeec5d47e595c34e8b1ce7f991472c3bf04adc17c372a384cdfa949702a7" and old.get("bytes") == 45809102 and "native_O3_anchor_second_Z_eta_chart.json" in str(old.get("path", "")), "second-Z parent alias changed")
    _, _, old_report = read_json(old["path"], old["sha256"])
    old_parent = old_report.get("parent_report")
    require(isinstance(old_parent, Mapping) and old_parent.get("sha256") == "4e926060bf944c752d74678cba8e13744abc2373b58c9a9916b9621f88dce072" and "native_O3_anchor_axial_Q_velocity.json" in str(old_parent.get("path", "")), "axial-Q ancestor alias changed")
    require("native_O3_anchor_first_picard_core.json" not in str(old_parent.get("path", "")), "stale first-Picard alias selected")

    maps = {key: authenticate_hash_map(report.get(key), f"parent {key}") for key in ("source_pins", "protected_file_hashes", "legacy_dirty_path_hashes", "preserved_failed_artifact_hashes_union")}
    require(len(maps["source_pins"]) == PARENT_SOURCE_PIN_COUNT and len(maps["protected_file_hashes"]) == PARENT_PROTECTED_FILE_COUNT and len(maps["legacy_dirty_path_hashes"]) == PARENT_DIRTY_PATH_COUNT and len(maps["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_ARTIFACT_COUNT, "parent preservation maps changed")
    retained = {key: report.get(key) for key in ("retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance", "preserved_error_cap_metadata")}
    construction = retained["retained_construction_data"]
    inventory = retained["retained_integral_term_inventory"]
    require(isinstance(construction, Mapping) and construction.get("cumulative_terms") == 212 and construction.get("reserved_contribution_registered_once") is True and construction.get("six_signed_channels_retained") is True, "212 cumulative packet changed")
    require(isinstance(inventory, Mapping) and inventory.get("reserved_contribution_registered_once") is True, "integral inventory changed")
    return {"report": report, "report_artifact": {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_bytes_exact": True, "raw_readonly": True, "hash_only": True}, **maps, **retained}


def capsule_status(capsule: Mapping[str, Any]) -> dict[str, Any]:
    """Authenticate capsule DATA shape without evaluating a coefficient."""
    require(capsule.get("schema") == CAPSULE_SCHEMA and capsule.get("component") == COMPONENT, "capsule identity changed")
    parameters = capsule.get("parameters")
    require(isinstance(parameters, Mapping), "capsule parameter contract changed")
    parameter_aliases = (("delta", "δ"), ("epsilon", "epsilon_core", "ε"), ("h",), ("sigma", "σ"), ("j",), ("s",), ("a",), ("U0a",), ("d_a",), ("L_a",), ("H1a",), ("H2a",))
    require(all(any(name in parameters for name in aliases) for aliases in parameter_aliases), "capsule parameter contract changed")
    source_atoms = capsule.get("source_atoms")
    require(isinstance(source_atoms, Mapping) and "F0_anchor" in source_atoms and source_atoms["F0_anchor"].get("materialized") is False, "formal F0 anchor changed")
    anchors = capsule.get("anchor_polynomials")
    require(isinstance(anchors, Mapping), "anchor polynomial mapping missing")
    require({"Phi1", "DU1", "PhiZ", "DUZ", "PhiZZ", "DUZZ"} <= set(anchors), "anchor polynomial keys changed")
    for name in ("Phi1", "DU1", "PhiZ", "DUZ", "PhiZZ", "DUZZ"):
        require(isinstance(anchors[name], list) and anchors[name], f"anchor polynomial {name} is not a degree interval list")
    norm = capsule.get("norm_data")
    require(isinstance(norm, Mapping) and {"E_norm", "Mphi", "Mpsi", "PhiZZ_error_existing", "DUZZ_error_existing"} <= set(norm), "norm data changed")
    mean_primitives = capsule.get("mean_primitives")
    require(isinstance(mean_primitives, Mapping) and set(PAYLOAD_ETA_LABELS) <= set(mean_primitives), "mean primitive eta rows changed")
    for key in PAYLOAD_ETA_LABELS:
        row = mean_primitives[key]
        require(isinstance(row, Mapping) and {"JD_candidate", "JDZ_candidate", "JD_enclosed", "JDZ_enclosed"} <= set(row), f"mean primitive envelope data changed: {key}")
    eta_rows = capsule.get("eta_rows")
    require(isinstance(eta_rows, Mapping) and set(PAYLOAD_ETA_LABELS) <= set(eta_rows), "eta row grid changed")
    for key in PAYLOAD_ETA_LABELS:
        row = eta_rows[key]
        require(isinstance(row, Mapping) and {"shift", "Ubase", "ratio", "ell1", "ell2_known"} <= set(row), f"eta row {key} changed")
        if key == "uniform_eta":
            require(row.get("ell2_known") is None or row.get("ell2_known") == [None] * 7, "uniform eta ell2 must remain unresolved")
        else:
            value = row.get("ell2_known")
            require(isinstance(value, (list, Mapping)), f"eta row {key} ell2 record changed")
    value_rows = capsule.get("value_rows", capsule.get("value_rows8x2DATA"))
    require(isinstance(value_rows, Mapping) and set(PAYLOAD_ETA_LABELS) <= set(value_rows), "value row grid changed")
    for key in PAYLOAD_ETA_LABELS:
        row = value_rows[key]
        require(isinstance(row, Mapping) and set(RHO_GRID) | {"uniform_rho"} <= set(row), f"value row {key} rho pair missing")
        for rho in (*RHO_GRID, "uniform_rho"):
            rho_row = row[rho]
            require(isinstance(rho_row, Mapping) and {"phi", "phi_minus_one", "DU"} <= set(rho_row), f"value row {key}.{rho} DATA fields changed")
            for name in ("phi", "phi_minus_one", "DU"):
                rec = rho_row[name]
                require(isinstance(rec, Mapping) and {"candidate", "error", "actual", "candidate_Z", "error_Z", "actual_Z"} <= set(rec), f"value row {key}.{rho}.{name} DATA changed")
    exact = capsule.get("exact_view_contract")
    require(isinstance(exact, Mapping), "exact-view contract missing")
    require(any(key in exact for key in ("old_ast_constructor", "oldASTconstructor", "old_AST_constructor", "old_core_checked_receipt")), "exact-view contract missing old AST/source receipt")
    require(any(key in exact for key in ("old_provider", "provider", "accepted_common_theorem")), "exact-view contract missing old provider/theorem")
    require(any(key in exact for key in ("theorem_refs", "theoremrefs", "theorems", "exact_field_definitions")), "exact-view contract missing theorem/field definitions")
    source_ids = capsule.get("source_ids")
    require(isinstance(source_ids, Mapping) and len(source_ids) == 10 and all(isinstance(key, str) and valid_sha(value) for key, value in source_ids.items()), "capsule source IDs changed")
    retained = capsule.get("retained_native_packet")
    require(isinstance(retained, Mapping), "retained native packet missing")
    require(retained.get("append_leading_moments_to_n1_ledger", retained.get("append_leading_moments_to_212_n1_ledger")) is False and retained.get("n1_source_values_present") is False and retained.get("old_core_to_eta_function_identity_proved") is False, "retained packet boundary changed")
    grid = capsule.get("grid", capsule.get("grid_data"))
    require(isinstance(grid, Mapping), "capsule grid missing")
    require(grid_tokens(grid.get("rho", ())) == RHO_GRID and grid_tokens(grid.get("uniform_rho", ())) == UNIFORM_RHO_GRID and grid_tokens(grid.get("eta", ())) == ETA_GRID and grid_tokens(grid.get("uniform_eta", ())) == UNIFORM_ETA_GRID, "capsule grid changed")
    require("tau" not in grid and "physical_tau" not in grid, "tau field is out of scope")
    return {"parameters": dict(parameters), "anchor_polynomials": dict(anchors), "norm_data": dict(norm), "eta_rows": dict(eta_rows), "value_rows": dict(value_rows), "exact_view_contract": dict(exact), "source_ids": dict(source_ids), "retained_native_packet": dict(retained), "grid": dict(grid)}


def capsule_from_spec(spec: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    path, raw, capsule_sha = authenticate_spec(spec, "capsule", CAPSULE_PATH)
    require(path.resolve() == CAPSULE_PATH.resolve(), "capsule path changed")
    try:
        capsule = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("capsule JSON parse failed") from exc
    require(isinstance(capsule, Mapping), "capsule JSON object required")
    value = capsule_status(capsule)
    canonical = json.dumps(capsule, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    canonical_sha = sha_bytes(canonical)
    return dict(capsule), {"path": wire(path), "sha256": capsule_sha, "bytes": len(raw), "canonical_sha256": canonical_sha, "schema": capsule["schema"], "component": capsule["component"], "status": value}


def load_active_brief() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    path, raw, brief = read_json(MAIN_BRIEF, MAIN_BRIEF_SHA256)
    require(path.resolve() == MAIN_BRIEF.resolve() and len(raw) == MAIN_BRIEF_BYTES, "active brief raw pin changed")
    require(brief.get("execution_rules") == EXECUTION_RULES and brief.get("execution_rules_preserved_exact") is True, "execution rules changed")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "goal state changed")
    accepted_progress = brief.get("accepted_progress")
    physical_baseline = brief.get("physical_baseline")
    require(isinstance(accepted_progress, Mapping), "accepted progress missing")
    if not isinstance(physical_baseline, Mapping):
        # The sealed parent brief keeps the baseline under accepted_progress;
        # carry that exact record rather than synthesizing a new value.
        physical_baseline = dict(accepted_progress)
    require(physical_baseline.get("commit") == PHYSICAL_BASELINE, "physical baseline changed")
    return raw, brief, {"accepted_progress": dict(accepted_progress), "physical_baseline": dict(physical_baseline), "execution_rules": list(EXECUTION_RULES)}


def encode_brief(value: Mapping[str, Any]) -> bytes:
    mutable = dict(value)
    raw = encode(mutable)
    if len(raw) < BRIEF_SOFT_LIMIT:
        return raw
    archive_path = WORK / f"{COMPONENT}_full_brief_metadata_archive_{sha_bytes(raw)[:16]}.json"
    archive = {"path": wire(archive_path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_bytes_exact": True, "raw_readonly": True}
    if archive_path.exists():
        require(archive_path.read_bytes() == raw, "brief archive collision")
    else:
        write_new(archive_path, raw)
    mutable["latest_brief_compaction_archive"] = {**archive, "metadata_only": True, "new_math_calls": 0, "all_original_sections_preserved": True}
    # Keep accepted progress, physical baseline, prior failed state and all cap
    # inventories in the active brief; only historical display aliases may be
    # projected through an exact archive pointer.
    protected = {"accepted_progress", "physical_baseline", "last_accepted_mathematical_state", "preserved_failed_artifact_hashes_union", "preserved_error_cap_metadata", "retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance", "execution_rules", "execution_rules_preserved_exact"}
    historical = [key for key in mutable if key not in protected and key.startswith(("latest_published_", "latest_accepted_", "previous_", "prior_", "unaccepted_"))]
    removed: list[str] = []
    raw = encode_compact(mutable)
    for key in sorted(historical, key=lambda name: len(encode(mutable[name])), reverse=True):
        if len(raw) < BRIEF_SOFT_LIMIT:
            break
        mutable.pop(key, None)
        removed.append(key)
        mutable["archived_historical_brief_sections"] = {**archive, "selectors": ["$." + name for name in removed], "preserved_exact": True, "metadata_only": True, "new_math_calls": 0}
        raw = encode_compact(mutable)
    require(len(raw) < BRIEF_SOFT_LIMIT, "updated brief exceeds 97 KB")
    return raw


def manifest_value(manifest: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in manifest:
            return manifest[name]
    return None


def load_manifest(value: str | Path) -> tuple[Path, bytes, dict[str, Any]]:
    path, raw, manifest = read_json(value)
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("component") == COMPONENT, "manifest identity changed")
    require(manifest_value(manifest, "expected_head", "current_head", "head") == CURRENT_HEAD_COMMIT, "manifest HEAD changed")
    parent_spec = manifest_value(manifest, "parent_report", "parent_report_pin")
    require(isinstance(parent_spec, Mapping) and parent_spec.get("sha256") == PARENT_REPORT_SHA256, "manifest parent report changed")
    require(str(parent_spec.get("path", "")).replace("\\", "/") == wire(PARENT_REPORT), "manifest parent report path changed")
    artifacts = manifest.get("artifacts", manifest)
    require(isinstance(artifacts, Mapping), "manifest artifacts missing")
    kernel_path, kernel_raw, kernel_sha = authenticate_spec(artifacts.get("kernel"), "manifest kernel", KERNEL_PATH)
    runner_path, runner_raw, runner_sha = authenticate_spec(artifacts.get("runner"), "manifest runner", RUNNER_PATH)
    require(kernel_sha != runner_sha or kernel_path.resolve() != runner_path.resolve(), "kernel and runner aliases are forbidden")
    kernel_ast = parse_ast(kernel_path, kernel_raw, "kernel")
    contract = ast_kernel_contract(kernel_path, kernel_raw)
    require(contract.get("component_literal") in (None, COMPONENT), "kernel component literal changed")
    require(contract["SUCCESS_COUNTS"].get("source_family_calls") == 1 and contract["SUCCESS_COUNTS"].get("new_math_family_calls") == 5 and contract["SUCCESS_COUNTS"].get("derive_calls") == 1 and contract["SUCCESS_COUNTS"].get("numeric_evaluation_calls") == 1 and contract["SUCCESS_COUNTS"].get("new_mathematical_proofs") == 1, "kernel lifecycle counts changed")
    for family in NEW_MATH_FAMILIES:
        require(contract["SUCCESS_COUNTS"].get(f"{family}_family_calls") == 1, f"kernel family count changed: {family}")
    callback_counts = {str(key): int(value) for key, value in contract["CALLBACK_EVENT_COUNTS"].items()}
    require(sum(callback_counts.values()) == 13 and callback_counts.get("proof_bundle") == 1, "callback contract must be 13 rows and one proof")

    wrappers = manifest.get("source_wrappers", manifest.get("wrappers"))
    if isinstance(wrappers, Mapping) and "capsule" not in wrappers and "source_capsule" in wrappers:
        wrappers = {"capsule": wrappers["source_capsule"]}
    require(isinstance(wrappers, Mapping) and "capsule" in wrappers, "manifest capsule wrapper missing")
    capsule, capsule_meta = capsule_from_spec(wrappers["capsule"])
    declared_counts = manifest_value(manifest, "expected_counts", "success_counts")
    declared_flags = manifest_value(manifest, "expected_new_flags", "new_flags", "expected_flags")
    declared_callback_raw = manifest_value(manifest, "callback_event_counts", "expected_callback_event_counts")
    require(isinstance(declared_counts, Mapping) and isinstance(declared_flags, Mapping) and isinstance(declared_callback_raw, Mapping), "manifest expected contract mappings missing")
    expected_counts = {str(key): int(value) for key, value in declared_counts.items()}
    expected_flags = dict(declared_flags)
    declared_callback = {str(key): int(value) for key, value in declared_callback_raw.items()}
    require(expected_counts == contract["SUCCESS_COUNTS"] and expected_flags == contract["SUCCESS_FLAGS"] and declared_callback == callback_counts, "manifest and kernel contracts differ")
    event_map = manifest_value(manifest, "count_event_map")
    require(isinstance(event_map, Mapping), "manifest event map missing")
    policy = manifest_value(manifest, "api_policy", "policy")
    require(isinstance(policy, Mapping) and policy.get("metadata_only") is False and policy.get("source_only") is True and policy.get("numeric_evaluation") is True and policy.get("full_goal_complete") is False and policy.get("physical_layer_completed") is False, "manifest policy changed")
    require(policy.get("actual_n1_solution_materialized", False) is False and policy.get("regular_n1_binding_installed", False) is False, "n1 solution scope changed")
    plan = manifest_value(manifest, "next_unexecuted_plan", "next_plan", "future_unexecuted_plan")
    require(isinstance(plan, Mapping) and plan.get("path") and valid_sha(plan.get("sha256")), "next plan pointer missing")
    return path, raw, {
        "manifest": {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw)},
        "kernel": {**kernel_ast, "sha256": kernel_sha},
        "runner": {"path": wire(runner_path), "sha256": runner_sha, "bytes": len(runner_raw), "ast_parse": "ok"},
        "kernel_contract": contract,
        "expected_counts": expected_counts,
        "expected_new_flags": expected_flags,
        "callback_event_counts": callback_counts,
        "count_event_map": {str(key): tuple(value) for key, value in event_map.items()},
        "source_wrappers": {"capsule": {"path": wire(CAPSULE_PATH), "sha256": capsule_meta["sha256"], "bytes": capsule_meta["bytes"], "canonical_sha256": capsule_meta["canonical_sha256"]}},
        "capsule": capsule,
        "capsule_meta": capsule_meta,
        "api_policy": dict(policy),
        "nexttask": NEXT_SUCCESS_TASK,
        "next_plan": dict(plan),
        "kernel_result_schema": manifest_value(manifest, "kernel_result_schema", "result_schema") or f"{COMPONENT}_kernel_result_v1",
    }


def merge_hashes(base: Mapping[str, str], additions: Mapping[str, str]) -> dict[str, str]:
    merged = dict(base)
    for path, expected in additions.items():
        if path in merged:
            require(merged[path] == expected, f"source hash collision changed: {path}")
        else:
            merged[path] = expected
    return merged


def build_source(parent: Mapping[str, Any], checked: Mapping[str, Any], input_path: Path, input_raw: bytes, input_value: Mapping[str, Any], brief_context: Mapping[str, Any], manifest_path: Path, manifest_raw: bytes) -> dict[str, Any]:
    additions = {
        checked["kernel"]["path"]: checked["kernel"]["sha256"],
        checked["runner"]["path"]: checked["runner"]["sha256"],
        checked["manifest"]["path"]: checked["manifest"]["sha256"],
        wire(input_path): sha_bytes(input_raw),
        checked["source_wrappers"]["capsule"]["path"]: checked["source_wrappers"]["capsule"]["sha256"],
    }
    source_union = merge_hashes(parent["source_pins"], additions)
    require(len(source_union) == FROZEN_SOURCE_UNION_COUNT, "frozen source union count changed")
    input_spec = {"path": wire(input_path), "sha256": sha_bytes(input_raw), "bytes": len(input_raw), "raw_bytes_exact": True, "raw_readonly": True}
    source = {
        "schema": SOURCE_SCHEMA,
        "component": COMPONENT,
        "phase": "prepared",
        "parent_commit": PARENT_COMMIT,
        "current_head_commit": CURRENT_HEAD_COMMIT,
        "accepted_math_commit": ACCEPTED_MATH_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "parent_report_grandparent_commit": PARENT_REPORT_GRANDPARENT_COMMIT,
        "parent_failure_preserved": True,
        "last_accepted_mathematical_state": {"commit": ACCEPTED_MATH_COMMIT, "component": "native_O3_eta_core_moment_functionals", "conditional": True, "source_only": True, "metadata_only": False, "full_goal_complete": False, "physical_layer_completed": False},
        "physical_baseline": brief_context["physical_baseline"],
        "input": dict(input_value),
        "input_spec": input_spec,
        "code2": {"kernel": checked["kernel"], "runner": checked["runner"]},
        "manifest": checked["manifest"],
        "source_wrappers": checked["source_wrappers"],
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "expected_counts": checked["expected_counts"],
        "expected_new_flags": checked["expected_new_flags"],
        "callback_event_counts": checked["callback_event_counts"],
        "callback_event_count": sum(checked["callback_event_counts"].values()),
        "count_event_map": checked["count_event_map"],
        "journal_event_count": sum(checked["callback_event_counts"].values()) + 2,
        "kernel_result_schema": checked["kernel_result_schema"],
        "nexttask": checked["nexttask"],
        "next_plan": checked["next_plan"],
        "new_math_families": list(NEW_MATH_FAMILIES),
        "rho_grid": list(RHO_GRID),
        "uniform_rho_grid": list(UNIFORM_RHO_GRID),
        "eta_grid": list(ETA_GRID),
        "uniform_eta_grid": list(UNIFORM_ETA_GRID),
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "numeric_values_computed": False,
        "source_exact_leading_view_conditional": True,
        "conditional_old_core_eta_view_registered": True,
        "unconditional_exact_view_admission": False,
        "regular_n1_solution_view_installed": False,
        "local_canonical_n1_rhs_C0_C1": True,
        "actual_n1_solution_materialized": False,
        "regular_n1_binding_installed": False,
        "matrix_jet_bound_only": True,
        "Z3_bound_only": True,
        "global_false_flags": {name: False for name in GLOBAL_FALSE_KEYS},
        "full_goal_complete": False,
        "physical_layer_completed": False,
        "runtime_model_metadata": None,
        "execution_rules": list(EXECUTION_RULES),
        "execution_rules_preserved_exact": True,
        "do_not_reenter": True,
    }
    return source


def update_brief(current_raw: bytes, current: Mapping[str, Any], patch: Mapping[str, Any]) -> dict[str, Any]:
    updated = dict(current)
    for key, value in patch.items():
        if key == "in_progress" and isinstance(updated.get(key), Mapping) and isinstance(value, Mapping):
            merged = dict(updated[key]); merged.update(dict(value)); updated[key] = merged
        else:
            updated[key] = value
    updated["execution_rules"] = list(EXECUTION_RULES)
    updated["execution_rules_preserved_exact"] = True
    new_raw = encode_brief(updated)
    write_replace(MAIN_BRIEF, new_raw)
    return {"path": wire(MAIN_BRIEF), "sha256": sha_bytes(new_raw), "bytes": len(new_raw), "preimage_sha256": sha_bytes(current_raw)}


def prepare(manifest_value_path: str) -> dict[str, Any]:
    manifest_path, manifest_raw, manifest = load_manifest(manifest_value_path)
    parent = load_parent()
    brief_raw, brief, brief_context = load_active_brief()
    input_path, input_raw, input = read_json(INPUT_PATH)
    require(input_path.resolve() == INPUT_PATH.resolve() and input.get("schema") == INPUT_SCHEMA and input.get("component") == COMPONENT, "input identity changed")
    require(input.get("current_head") == CURRENT_HEAD_COMMIT or input.get("parent_commit") == CURRENT_HEAD_COMMIT, "input HEAD changed")
    require(input.get("metadata_only") is False and input.get("source_only") is True and input.get("numeric_evaluation") is True and input.get("full_goal_complete") is False and input.get("physical_layer_completed") is False, "input scope changed")
    require(input.get("actual_n1_solution_materialized", False) is False and input.get("regular_n1_binding_installed", False) is False, "input n1 scope changed")
    input_capsule = input.get("capsule", input.get("source_capsule"))
    require(isinstance(input_capsule, Mapping) and input_capsule.get("sha256") == manifest["source_wrappers"]["capsule"]["sha256"], "input capsule pointer changed")
    preimage_path = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(brief_raw)[:12]}.json"
    require(not preimage_path.exists(), "main brief preimage already exists")
    source = build_source(parent, manifest, input_path, input_raw, input, brief_context, manifest_path, manifest_raw)
    preimage_raw = write_new(preimage_path, brief_raw)
    source["main_brief_preimage"] = {"path": wire(preimage_path), "sha256": sha_bytes(preimage_raw), "bytes": len(preimage_raw), "raw_bytes_exact": True, "raw_readonly": True}
    source_raw = write_new(SOURCE_PATH, source)
    pending = {
        "schema": f"{COMPONENT}_brief_v1",
        "component": COMPONENT,
        "phase": "prepared",
        "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw), "bytes": len(source_raw)},
        "input_spec": source["input_spec"],
        "code2": source["code2"],
        "manifest": source["manifest"],
        "source_wrappers": source["source_wrappers"],
        "expected_counts": source["expected_counts"],
        "expected_new_flags": source["expected_new_flags"],
        "callback_event_counts": source["callback_event_counts"],
        "count_event_map": source["count_event_map"],
        "main_brief_preimage": source["main_brief_preimage"],
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "numeric_values_computed": False,
        "actual_n1_solution_materialized": False,
        "regular_n1_binding_installed": False,
        "runtime_model_metadata": None,
        "execution_rules": list(EXECUTION_RULES),
        "execution_rules_preserved_exact": True,
        "no_run_performed": True,
    }
    pending_raw = write_new(BRIEF_PATH, pending)
    main_brief = update_brief(brief_raw, brief, {"phase": "prepared", "current_component": COMPONENT, "component": COMPONENT, "goal_active": True, "full_goal_complete": False, "next_task": NEXT_SUCCESS_TASK, "next_executable_task": NEXT_SUCCESS_TASK, "in_progress": {COMPONENT: {"phase": "prepared", "source": wire(SOURCE_PATH), "source_sha256": sha_bytes(source_raw), "input": wire(INPUT_PATH), "input_sha256": sha_bytes(input_raw), "metadata_only": False, "source_only": True, "numeric_evaluation": True, "numeric_values_computed": False, "do_not_reenter": True}}, "source_exact_leading_view_conditional": True, "local_canonical_n1_rhs_C0_C1": True, "actual_n1_solution_materialized": False, "regular_n1_binding_installed": False})
    return {"status": "prepared", "component": COMPONENT, "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw), "bytes": len(source_raw)}, "pending_brief": {"path": wire(BRIEF_PATH), "sha256": sha_bytes(pending_raw), "bytes": len(pending_raw)}, "main_brief": main_brief, "source_union_count": source["source_union_count"], "expected_counts": source["expected_counts"], "expected_new_flags": source["expected_new_flags"], "callback_event_counts": source["callback_event_counts"], "callback_event_count": source["callback_event_count"], "journal_event_count": source["journal_event_count"], "metadata_only": False, "source_only": True, "numeric_evaluation": True, "numeric_values_computed": False, "actual_n1_solution_materialized": False, "regular_n1_binding_installed": False, "runtime_model_metadata": None, "no_run_performed": True}


def load_source_for_run() -> tuple[bytes, dict[str, Any]]:
    source_path, source_raw, source = read_json(SOURCE_PATH)
    require(source_path.resolve() == SOURCE_PATH.resolve() and source.get("schema") == SOURCE_SCHEMA and source.get("component") == COMPONENT, "source checkpoint identity changed")
    require(source.get("parent_commit") == CURRENT_HEAD_COMMIT and source.get("current_head_commit") == CURRENT_HEAD_COMMIT and source.get("do_not_reenter") is True, "source lineage changed")
    require(source.get("parent_failure_preserved") is True and source.get("metadata_only") is False and source.get("source_only") is True and source.get("numeric_evaluation") is True and source.get("full_goal_complete") is False and source.get("physical_layer_completed") is False, "source scope changed")
    require(source.get("runtime_model_metadata") is None and source.get("execution_rules") == EXECUTION_RULES, "source execution policy changed")
    preimage = source.get("main_brief_preimage")
    require(isinstance(preimage, Mapping), "brief preimage missing")
    read_bytes(preimage["path"], preimage["sha256"])
    input_path, input_raw, input = read_json(source["input_spec"]["path"], source["input_spec"]["sha256"])
    require(input_path.resolve() == INPUT_PATH.resolve() and input == source.get("input"), "input literal changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})) == FROZEN_SOURCE_UNION_COUNT, "source union count changed")
    for path, expected in source["source_union_hashes"].items():
        read_bytes(path, expected)
    manifest_path, manifest_raw, manifest = read_json(source["manifest"]["path"], source["manifest"]["sha256"])
    _, _, checked = load_manifest(manifest_path)
    require(checked["expected_counts"] == source["expected_counts"] and checked["expected_new_flags"] == source["expected_new_flags"] and checked["callback_event_counts"] == source["callback_event_counts"], "manifest contract changed")
    capsule = source["source_wrappers"]["capsule"]
    _, _, capsule_value = read_json(capsule["path"], capsule["sha256"])
    capsule_status(capsule_value)
    return source_raw, source


def append_journal(row: Mapping[str, Any]) -> None:
    raw = encode_compact(dict(row))
    with JOURNAL_PATH.open("ab") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())


def journal_rows() -> list[dict[str, Any]]:
    if not JOURNAL_PATH.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in JOURNAL_PATH.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "journal row is not an object")
            rows.append(value)
    return rows


def import_kernel(path: Path) -> Any:
    require(path.resolve() != Path(__file__).resolve(), "kernel aliases runner")
    spec = importlib.util.spec_from_file_location(f"_{COMPONENT}_entry", path)
    require(spec is not None and spec.loader is not None, "kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = module
    require(sys.modules.get(spec.name) is module, "kernel module registration failed")
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "derive API missing")
    return module


def callback_rows_callback() -> tuple[Any, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []

    def emit(*args: Any) -> None:
        require(len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], Mapping), "emit must be exactly (event,row)")
        row = {"event": args[0], "row": dict(args[1])}
        append_journal({"schema": JOURNAL_SCHEMA, "event": "callback", **row})
        rows.append(row)

    return emit, rows


def count_prefix(source: Mapping[str, Any], rows: list[Mapping[str, Any]]) -> dict[str, int]:
    actual = {str(key): 0 for key in source["expected_counts"]}
    for row in rows:
        for target in source["count_event_map"].get(row.get("event"), ()):
            require(target in actual, f"unknown count target: {target}")
            actual[target] += 1
    return actual


def false_flags(source: Mapping[str, Any]) -> dict[str, bool]:
    return {str(name): False for name in source["expected_new_flags"]}


def result_payload(value: Mapping[str, Any]) -> Mapping[str, Any]:
    for candidate in (value, value.get("result") if isinstance(value.get("result"), Mapping) else {}):
        payload = candidate.get("payload") if isinstance(candidate, Mapping) else None
        if isinstance(payload, Mapping):
            return payload
    return {}


def mapping_record(value: Any, label: str, *, point_value: bool = False) -> Mapping[str, Any]:
    require(isinstance(value, Mapping), f"{label} record missing")
    if "point_value" in value:
        require(value.get("point_value") is point_value, f"{label} point-value scope changed")
    require(value.get("materialized", False) is False, f"{label} was materialized")
    return value


def validate_kernel_payload(value: Mapping[str, Any]) -> None:
    payload = result_payload(value)
    require(isinstance(payload.get("exact_view_registration"), Mapping), "exact-view registration missing")
    exact = payload["exact_view_registration"]
    require(exact.get("conditional") is True or exact.get("conditional_registration") is True or exact.get("source_exact_leading_view_conditional") is True or exact.get("conditional_leading_identity_registered") is True, "exact view is not conditional")
    require(exact.get("unconditional_admission", exact.get("admitted_unconditionally", False)) is False, "exact view was admitted unconditionally")
    require(exact.get("regular_n1_solution_view_installed", exact.get("regular_n1_solution_view", False)) is False, "regular n1 solution view was installed")
    require(exact.get("values_materialized", exact.get("n1_values_materialized", False)) is False, "exact view values were materialized")
    mixed = payload.get("mixed_jet_table")
    require((isinstance(mixed, list) or isinstance(mixed, Mapping)) and mixed, "mixed jet table missing")
    require(payload.get("mixed_jet_source_copy_count", 1) == 1, "mixed jet table is not a minimal one-copy table")
    gauge = payload.get("gauge_table")
    require((isinstance(gauge, list) or isinstance(gauge, Mapping)) and gauge, "gauge table missing")
    definitions = payload.get("definitions")
    proofpacket = payload.get("proofpacket", payload.get("proof_packet"))
    require(isinstance(definitions, Mapping) and isinstance(proofpacket, Mapping), "definitions/proofpacket missing")
    exact_definitions = definitions.get("exact_view")
    if not isinstance(exact_definitions, Mapping):
        exact_definitions = {}
    pressure_aliases = {"Paxis_scaled": ("Paxis_scaled", "Paxis_eta_exact"), "Ptotal_scaled": ("Ptotal_scaled", "Ptotal_eta_exact"), "oldPphysical": ("oldPphysical", "oldPphysical_eta_exact"), "scaled_pressure_identity": ("scaled_pressure_identity",), "p0_chart_is_envelope": ("p0_chart_is_envelope", "p0_chart_is_exact_function"), "Icandidate_is_envelope": ("Icandidate_is_envelope", "I_candidate_is_exact_function")}
    envelope_declaration = definitions.get("candidate_polynomials_are_envelopes") is True or exact_definitions.get("candidate_polynomials_are_envelopes") is True or proofpacket.get("candidate_polynomials_are_envelopes") is True
    for name, aliases in pressure_aliases.items():
        require(any(alias in definitions or alias in exact_definitions or alias in proofpacket for alias in aliases) or (name in ("p0_chart_is_envelope", "Icandidate_is_envelope") and envelope_declaration), f"pressure definition missing: {name}")
    require(envelope_declaration, "candidate pressure envelope declaration missing")
    rhs = payload.get("rhs_rows")
    require(isinstance(rhs, list) and len(rhs) == 48, "rhs table must contain 48 rows")
    expected = {(eta, rho, quantity) for eta in PAYLOAD_ETA_LABELS for rho in PAYLOAD_RHO_LABELS for quantity in ("g4", "g5", "g6")}
    actual = set()
    for index, row in enumerate(rhs):
        require(isinstance(row, Mapping), f"rhs row {index} is not a mapping")
        eta = str(row.get("eta")); rho = str(row.get("rho")); quantity = row.get("quantity")
        require((eta, rho, quantity) in expected, f"rhs selector changed at row {index}")
        actual.add((eta, rho, quantity))
        require(row.get("atom") in ("one", "F0_anchor"), f"rhs atom changed at row {index}")
        require((isinstance(row.get("sectors"), list) or isinstance(row.get("sectors"), Mapping)) and row["sectors"], f"rhs sectors missing at row {index}")
        for field in ("C0", "C1", "error_C0", "error_C1"):
            mapping_record(row.get(field), f"rhs {field} row {index}")
        for field in ("candidate_C0", "candidate_C1"):
            candidate = row.get(field)
            if candidate is not None:
                mapping_record(candidate, f"rhs {field} row {index}")
        if quantity == "g5":
            require(row.get("formal_positive_F0_coefficient", row.get("formal_positive", False)) is True or row.get("atom") == "F0_anchor", f"g5 row {index} lost formal F0 status")
        else:
            require(row.get("numeric_interval", True) is True, f"{quantity} row {index} is not numeric interval DATA")
    require(actual == expected, "rhs selectors are not exactly 8 eta by 2 rho by 3 quantities")
    summary = payload.get("summary") if isinstance(payload.get("summary"), Mapping) else {}
    require(summary.get("uniform_eta_direct_row") is True and summary.get("m3_candidate") is None and summary.get("m3_source_bound_only") is True and summary.get("uniform_rho") == "real [0,4] direct interval", "direct uniform eta/rho or Z3 boundary changed")
    require(proofpacket.get("m3_actual_jet_unresolved") is True, "third-Z signed value boundary changed")
    proof_domain = proofpacket.get("domain") if isinstance(proofpacket.get("domain"), Mapping) else {}
    require(proof_domain.get("extension_4_1_used") is False, "out-of-domain 4.1 extension changed")
    require(proofpacket.get("old_212_error_recertification") is False and proofpacket.get("leading_moments_appended_to_212") is False and proofpacket.get("actual_independent_pressure_target_retained") is True and proofpacket.get("reserved_own_swirl_retained_once") is True and proofpacket.get("no_N_materialization") is True and proofpacket.get("no_old_provider_calls") is True and proofpacket.get("no_checker_or_selftest") is True and proofpacket.get("no_global_closure") is True and proofpacket.get("no_time_recursion") is True and proofpacket.get("no_full_physical_field") is True and proofpacket.get("no_NS_claim") is True, "retained source packet boundary changed")
    bound_matrix = payload.get("matrix_jet_bound_only") is True or summary.get("m3_source_bound_only") is True
    bound_z3 = payload.get("Z3_bound_only") is True or summary.get("m3_source_bound_only") is True
    require(bound_matrix and bound_z3, "matrix/Z3 bound-only flags missing")
    local_rhs = payload.get("local_canonical_n1_rhs_C0_C1") is True or all(isinstance(row, Mapping) and "C0" in row and "C1" in row for row in rhs)
    source_conditional = payload.get("source_exact_leading_view_conditional") is True or payload.get("conditional") is True
    require(local_rhs and source_conditional, "local rhs/exact-view flags changed")
    require(payload.get("actual_n1_solution_materialized", False) is False and payload.get("regular_n1_binding_installed", False) is False, "regular n1 scope changed")
    for name in ("global_five_moment", "global_moment", "full_physical", "global_NS", "time_n_recursion", "full_goal"):
        if name in payload:
            require(payload[name] is False, f"global flag changed: {name}")
    require("moment_table" not in payload and "leading80_moment_table" not in payload and "tau" not in payload, "out-of-scope table or tau field present")


def kernel_failure(source: Mapping[str, Any], exc: BaseException, rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int) -> dict[str, Any]:
    failure = {"kind": "manager_import_or_derive_failure", "failure_stage": "derive" if imported else "import", "api_guard_diagnostic": type(exc).__name__, "error_type": type(exc).__name__, "error": repr(exc), "traceback": traceback.format_exc(), "kernel_imported": imported, "derive_entry_calls": entry_calls, "derive_calls": calls, "callback_prefix": list(rows), "callback_prefix_count": len(rows), "callback_prefix_counts": count_prefix(source, rows), "stopped": True}
    return {"schema": source["kernel_result_schema"], "component": COMPONENT, "passed": False, "failure": failure}


def wrap_result(source: Mapping[str, Any], kernel_result: Mapping[str, Any], imported: bool, entry_calls: int, calls: int, rows: list[dict[str, Any]], manager_failure: Mapping[str, Any] | None) -> dict[str, Any]:
    passed = kernel_result.get("passed") is True and manager_failure is None
    nested = dict(kernel_result)
    if not passed:
        nested["new_flags"] = false_flags(source)
        nested["proof_count"] = 0
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if passed else "failed",
        "passed": passed,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "numeric_values_computed": passed,
        "current_leaf_numeric_evaluation": True,
        "source_exact_leading_view_conditional": True,
        "conditional_old_core_eta_view_registered": True,
        "unconditional_exact_view_admission": False,
        "regular_n1_solution_view_installed": False,
        "local_canonical_n1_rhs_C0_C1": True,
        "actual_n1_solution_materialized": False,
        "regular_n1_binding_installed": False,
        "matrix_jet_bound_only": True,
        "Z3_bound_only": True,
        "full_goal_complete": False,
        "physical_layer_completed": False,
        "global_false_flags": {name: False for name in GLOBAL_FALSE_KEYS},
        "result": nested,
        "actual_counts": dict(kernel_result.get("actual_counts", count_prefix(source, rows))) if passed else count_prefix(source, rows),
        "new_flags": dict(kernel_result.get("new_flags", false_flags(source))) if passed else false_flags(source),
        "proof_count": int(kernel_result.get("proof_count", 0)) if passed else 0,
        "derive_entry_calls": entry_calls,
        "derive_calls": calls,
        "kernel_imported": imported,
        "callback_events": list(rows),
        "callback_event_count": len(rows),
        "callback_prefix_counts": count_prefix(source, rows),
        "failure": dict(manager_failure) if manager_failure else nested.get("failure"),
        "mathematical_proof_accepted": passed,
        "source_operator_implementation_accepted": passed,
        "math_source_identity_accepted": passed,
        "runtime_model_metadata": None,
        "source_checkpoint_sha256": sha_bytes(SOURCE_PATH.read_bytes()),
    }


def validate_saved(source: Mapping[str, Any], saved: Mapping[str, Any], rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int) -> dict[str, Any]:
    violations: list[str] = []
    if saved.get("schema") != RESULT_SCHEMA: violations.append("result schema")
    if saved.get("component") != COMPONENT: violations.append("result component")
    if saved.get("runtime_model_metadata") is not None: violations.append("runtime metadata")
    if saved.get("metadata_only") is not False or saved.get("source_only") is not True or saved.get("numeric_evaluation") is not True: violations.append("scope")
    if saved.get("actual_n1_solution_materialized") is not False or saved.get("regular_n1_binding_installed") is not False: violations.append("n1 materialization")
    allowed = set(source["callback_event_counts"])
    if any(row.get("event") not in allowed for row in rows): violations.append("callback event name")
    callback_counts = {name: sum(1 for row in rows if row.get("event") == name) for name in source["callback_event_counts"]}
    if saved.get("passed") is True:
        if saved.get("actual_counts") != source["expected_counts"]: violations.append("actual counts")
        if saved.get("new_flags") != source["expected_new_flags"]: violations.append("new flags")
        if saved.get("proof_count") != source["expected_counts"].get("new_mathematical_proofs", 1): violations.append("proof count")
        if not (imported and entry_calls == 1 and calls == 1): violations.append("lifecycle")
        if callback_counts != source["callback_event_counts"]: violations.append("callback counts")
        if len(journal_rows()) != source["journal_event_count"]: violations.append("journal count")
        try:
            validate_kernel_payload(saved["result"])
        except Exception as exc:
            violations.append(f"payload: {exc}")
    else:
        if saved.get("proof_count") != 0: violations.append("failed proof count")
        if len(journal_rows()) != 1 + int(imported) + len(rows): violations.append("failed journal prefix")
    return {"schema": f"{COMPONENT}_validation_v1", "component": COMPONENT, "passed": not violations and saved.get("passed") is True, "attempt_failed": saved.get("passed") is not True, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "violations": violations}


def run() -> dict[str, Any]:
    source_raw, source = load_source_for_run()
    require(not ENTERED_PATH.exists() and not JOURNAL_PATH.exists() and not RESULT_PATH.exists(), "run checkpoints already exist; reentry forbidden")
    entered = {"schema": ENTERED_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "kernel": source["code2"]["kernel"], "runner": source["code2"]["runner"], "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0, "sys_modules_registration_before_exec": True, "metadata_only": False, "source_only": True, "numeric_evaluation": True, "numeric_values_computed": False, "source_exact_leading_view_conditional": True, "local_canonical_n1_rhs_C0_C1": True, "actual_n1_solution_materialized": False, "regular_n1_binding_installed": False, "runtime_model_metadata": None, "entered_at_utc": now(), "do_not_reenter": True}
    entered_raw = write_new(ENTERED_PATH, entered)
    write_new(JOURNAL_PATH, b"")
    append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_entry_before_import", "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw), "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0})
    emit, rows = callback_rows_callback()
    imported = False; entry_calls = 0; calls = 0; manager_failure: Mapping[str, Any] | None = None
    try:
        module = import_kernel(KERNEL_PATH)
        imported = True; entry_calls = 1
        append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_derive_entry", "import_attempts": 1, "derive_entry_calls": 1, "derive_calls": 0, "sys_modules_registered": True})
        calls = 1
        capsule_ref = source["source_wrappers"]["capsule"]
        _, _, capsule = read_json(capsule_ref["path"], capsule_ref["sha256"])
        capsule_status(capsule)
        data = dict(source["input"])
        data["capsule"] = capsule
        returned = module.derive(data, emit)
        require(isinstance(returned, Mapping), "derive returned a non-mapping")
        kernel_result = dict(returned)
    except BaseException as exc:
        kernel_result = kernel_failure(source, exc, rows, imported, entry_calls, calls)
        manager_failure = kernel_result["failure"]
    result_value = wrap_result(source, kernel_result, imported, entry_calls, calls, rows, manager_failure)
    result_raw = write_new(RESULT_PATH, result_value)
    _, saved_raw, saved = read_json(RESULT_PATH, sha_bytes(result_raw))
    validation = validate_saved(source, saved, rows, imported, entry_calls, calls)
    return {"status": "run_complete", "component": COMPONENT, "attempt_status": saved.get("attempt_status"), "passed": bool(validation["passed"]), "metadata_only": False, "numeric_evaluation": True, "numeric_values_computed": bool(saved.get("passed")), "source_exact_leading_view_conditional": True, "local_canonical_n1_rhs_C0_C1": True, "actual_n1_solution_materialized": False, "regular_n1_binding_installed": False, "result": {"path": wire(RESULT_PATH), "sha256": sha_bytes(saved_raw), "bytes": len(saved_raw), "schema": saved.get("schema")}, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "validation": validation, "failure": saved.get("failure"), "runtime_model_metadata": None}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare"); prep.add_argument("--manifest", required=True)
    sub.add_parser("run")
    args = parser.parse_args(argv)
    try:
        receipt = prepare(args.manifest) if args.command == "prepare" else run()
    except Exception as exc:
        print(json.dumps({"status": "error", "component": COMPONENT, "error": f"{type(exc).__name__}: {exc}", "metadata_only": False, "runtime_model_metadata": None}, sort_keys=True))
        return 2
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt.get("passed", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
