"""One-shot metadata manager for the native O3 cylindrical pullback; lifecycle only, never the mathematical kernel."""
from __future__ import annotations
import argparse, base64, hashlib, importlib.util, json, os, sys
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence
REPO = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
WORK = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")
COMPONENT = "native_O3_cylindrical_pullback"
INPUT_SCHEMA = "native_O3_cylindrical_pullback_input_v1"
RESULT_SCHEMA = "original_native_O3_cylindrical_pullback"
PUBLICATION_SCHEMA = "original_native_O3_cylindrical_pullback_publication_v1"
ENTERED_SCHEMA = COMPONENT + "_entered_v1"
DIAGNOSTIC_SCHEMA = COMPONENT + "_diagnostic_v1"
PARENT_REPORT = REPO / "experiments/root_st073/lei_ren_part1_paper_compliant_current_native_history_radius_bundle.json"
PARENT_RESULT = REPO / "tmp/native_history_radius_bundle/result.json"
PARENT_REPORT_SHA256 = "7786f54ac15bf400cba235613ee5c2c89fd4a2556cb36eda3d1f1cfa58b43a93"
PARENT_RESULT_SHA256 = "b541342e47bed12ed9b7394facef236d03b4828b6436776739245715b93e0df9"
PARENT_COMPONENT = "native_history_radius_bundle"
PARENT_COMMIT = "c96d1ba093a82579743b5c82dc8341d63d1c2862"
CURRENT_HEAD_COMMIT = PARENT_COMMIT
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "ce34f9ad75a78714766bebca4d376cf2510e5499c054753fc5a0104f2cd58d5d"
MAIN_BRIEF_BYTES = 75965
KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_O3_cylindrical_pullback.py")
RUNNER_RELATIVE = Path("scripts/run_native_O3_cylindrical_pullback.py")
API_PATH = WORK / "native_O3_cylindrical_pullback_api_2026_10_27.json"
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
REPORT_RELATIVE = Path("experiments/root_st073/lei_ren_part1_paper_compliant_current_native_O3_cylindrical_pullback.json")
DOC_RELATIVE = Path("docs/CURRENT_NATIVE_O3_CYLINDRICAL_PULLBACK_2026_10_27.md")
HANDOFF_RELATIVES = tuple(Path(path) for path in ("docs/AGENT_TASKS.md", "docs/CROSS_DEVICE_HANDOFF.md", "docs/CURRENT_CHECKPOINT.md", "docs/NEXT_TASKS.md"))
SCOPED_RELATIVES = (Path(".gitattributes"), KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, DOC_RELATIVE, *HANDOFF_RELATIVES)
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
JOURNAL_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (SOURCE_CHECKPOINT, ENTERED_CHECKPOINT, RESULT_CHECKPOINT, JOURNAL_CHECKPOINT)
INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
FACT_SPECS = {"rm_trace": (WORK / "native_Rm_axial_dependency_trace_2026_10_27.json", "901ed6eca9505e20eb7b3819ddc6d87d8e162974f37132655c05679b4c028f3f", "native_Rm_axial_dependency_trace_v1"), "phase_trace": (WORK / "native_radial_phase_semantics_trace_2026_10_27.json", "8255f24dc44c6d0a47aad3c854b7a197cbf2fec809cbd965e715df3f5f3c49fb", "native_radial_phase_semantics_trace_v1"), "overlap_audit": (WORK / "native_radius_pullback_overlap_audit_2026_10_27.json", "d163d3997ad3b7b85d5dc5b4eba066a0d7b4ab398d3cc608b5c3377cf2ea7b13", "native_radius_pullback_overlap_audit_v1")}
PARENT_SOURCE_PIN_COUNT = 2265
PARENT_CHECKPOINT_COUNT = 264
PARENT_FAILED_COUNT = 183
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
PUBLISHED_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_COUNT + 7
EXPECTED_COUNTS = {"branch_frame_bindings": 2, "new_algebra_family_calls": 1, "new_mathematical_proofs": 1, "jacobian_identity_rows": 1, "inverse_coordinate_rows": 9, "cartesian_commutator_rows": 3, "scalar_laplacian_rows": 1, "branch_operator_manifests": 2, "policy_rejections": 3, "numeric_calls": 0, "source_value_provider_calls": 0, "old_math_calls": 0, "field_definition_replays": 0, "N_materializations": 0}
TRUE_FLAGS = ("conditional_native_chart_metric_identity_proved", "native_branch_total_derivative_policy_enforced", "native_cylindrical_volume_density_defined")
FALSE_FLAGS = ("actual_native_field_pullback_instantiated", "native_Rm_Z_zero_proved", "actual_physical_time_map_installed", "phase_identity_proved", "physical_layer_completed", "full_goal_complete", "caps_deleted", "numeric_closure", "temporal_recursion")
ALL_FLAGS = (*TRUE_FLAGS, *FALSE_FLAGS)
REPORT_SCOPE_ATTRIBUTES = ("parent_ledgers_and_false_flags", "retained_error_caps_and_construction_maps", "input_literal_and_static_facts", "kernel_runner_api_artifacts", "durable_entry_result_and_journal", "conditional_counts_and_flags", "independent_saved_byte_review", "conditional_geometry_documentation", "publication_read_first_and_scope_integrity")
EXECUTION_RULES = ["Never rerun any entered ancestor mathematical/helper family, including failed leaves, old builders/providers/normalizers/checkers, source-value reducers, quadrature or recurrence. Static source/AST/literal JSON/hash inspection is allowed.", "Before a NEW mathematical/source-expression consumer: freeze code and input, persist source checkpoint and pending brief, persist entry before invocation, invoke once, save complete result before review. Preserve a failure and repair in a separate frozen leaf; never reenter.", "Authenticate raw gzip SHA BEFORE decoding. Never decode docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json or the huge Rm_loop_density_identity.json; raw SHA only.", "Never restart stopped session 17273. Retain nonzero incoming histories and ordinary Z1. Bounds, midpoints and interval overlap are not function values; remove old numerical errors only with a same-source packet theorem.", "Use current tree and GitHub; preserve all dirty and untracked files. Commit/push only exact accepted scoped files. Goal remains ACTIVE / INCOMPLETE.", "Children: luna_max_worker, gpt-5.6-luna/max fixed role, fork_turns:none, complete brief. Runtime model metadata unavailable means null; never claim it was verified.", "A worker function that checks mathematical source identities counts as a source consumer even when it consumes static receipts or is called a pure self-check. Workers return unexecuted code; root alone freezes code and inputs, records entry, invokes once and saves the complete result before review. AST parsing and metadata-only audits remain permitted. Do not execute or reconstruct the rejected memory draft result."]
def require(condition: Any, message: str) -> None:
    if not condition:
        raise AssertionError(message)
def clone(value: Any) -> Any:
    return deepcopy(value)
def encode(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8") + b"\n"
def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()
def valid_sha(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(char in "0123456789abcdefABCDEF" for char in value)
def wire(path: str | Path) -> str:
    return Path(path).as_posix()
def resolve_path(value: str | Path, base: Path = REPO) -> Path:
    path = Path(value)
    return path if path.is_absolute() else next((candidate for candidate in (base / path, REPO / path, WORK / path) if candidate.exists()), base / path)
def read_bytes(path: str | Path, expected: str | None = None) -> bytes:
    path = Path(path)
    require(path.exists(), f"Missing artifact: {wire(path)}")
    raw = path.read_bytes()
    if expected is not None: require(sha_bytes(raw) == expected.lower(), f"Changed artifact: {wire(path)}")
    return raw
def read_json(path: str | Path, expected: str | None = None) -> tuple[bytes, dict[str, Any]]:
    raw = read_bytes(path, expected)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"JSON object required: {wire(path)}")
    return raw, value
def _write(path: str | Path, value: Any, mode: str) -> bytes:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open(mode) as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw
def write_new(path: str | Path, value: Any) -> bytes:
    return _write(path, value, "xb")
def write_replace(path: str | Path, value: Any) -> bytes:
    return _write(path, value, "wb")
def append_journal(event: Mapping[str, Any]) -> bytes:
    return _write(JOURNAL_CHECKPOINT, encode(dict(event)), "ab")
def journal_events() -> list[dict[str, Any]]:
    if not JOURNAL_CHECKPOINT.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in JOURNAL_CHECKPOINT.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "Journal rows must be JSON objects")
            rows.append(value)
    return rows
def artifact_from_bytes(path: str | Path, raw: bytes, *, embed: bool = True) -> dict[str, Any]:
    result: dict[str, Any] = {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_bytes_exact": True, "raw_readonly": True}
    if embed: result["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
    return result
def artifact(path: str | Path, expected: str | None = None, *, embed: bool = True) -> dict[str, Any]:
    return artifact_from_bytes(path, read_bytes(path, expected), embed=embed)
def merge_hashes(*maps: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for mapping in maps:
        require(isinstance(mapping, Mapping), "Hash map must be a mapping")
        for name, digest in mapping.items():
            path = wire(name)
            require(valid_sha(digest), f"Malformed SHA for {path}")
            digest = digest.lower()
            require(path not in result or result[path] == digest, f"Conflicting SHA for {path}")
            result[path] = digest
    return result
def brief_safe(value: Any) -> Any:
    if isinstance(value, dict): return {key: brief_safe(child) for key, child in value.items() if key not in {"raw_bytes_base64", "callback_event_raw_rows"}}
    if isinstance(value, list): return [brief_safe(child) for child in value]
    return value
def load_current_brief() -> tuple[bytes, dict[str, Any]]:
    raw, brief = read_json(MAIN_BRIEF, MAIN_BRIEF_SHA256)
    require(len(raw) == MAIN_BRIEF_BYTES, "Current main brief byte count changed"); require(brief.get("execution_rules") == EXECUTION_RULES, "The restored seven execution rules changed"); require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "Goal is not ACTIVE / INCOMPLETE"); require(brief.get("metadata_only") is True and brief.get("source_only") is True, "Current accepted baseline markers changed")
    return raw, brief
def brief_preimage(raw: bytes) -> dict[str, Any]:
    path = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(raw)[:12]}.json"
    return {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw), "preserved_exact_bytes": True, "raw_bytes_base64": base64.b64encode(raw).decode("ascii")}
def parent_report_view(report: Mapping[str, Any]) -> dict[str, Any]:
    keys = ("schema", "component", "passed", "metadata_only", "source_only", "physical_layer_completed", "full_goal_complete", "actual_counts", "new_flags", "expected_counts", "expected_new_flags", "parent_false_flags", "parent_false_flag_values", "next_executable_task", "runtime_model_metadata")
    return {key: clone(report[key]) for key in keys if key in report}
def load_parent() -> dict[str, Any]:
    report_raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    result_raw, result = read_json(PARENT_RESULT, PARENT_RESULT_SHA256)
    require(report.get("component") == PARENT_COMPONENT and report.get("passed") is True, "Parent report identity changed"); require(report.get("metadata_only") is True and report.get("source_only") is True and report.get("full_goal_complete") is False, "Parent report markers changed"); require(result.get("component") == PARENT_COMPONENT and result.get("passed") is True, "Parent result identity changed"); require(report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256, "Parent result pointer changed")
    required_maps = ("source_pins", "protected_file_hashes", "legacy_dirty_path_hashes", "preserved_failed_artifact_hashes_union", "parent_false_flags", "parent_false_flag_values", "preserved_error_cap_metadata", "retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance")
    for key in required_maps:
        require(key in report, f"Parent ledger missing: {key}")
    require(len(report["source_pins"]) == PARENT_SOURCE_PIN_COUNT and len(report["protected_file_hashes"]) == PARENT_PROTECTED_COUNT and len(report["legacy_dirty_path_hashes"]) == PARENT_DIRTY_COUNT and len(report["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_COUNT, "Parent ledger count changed")
    require(len(report["parent_false_flags"]) == PARENT_FALSE_FLAG_COUNT and report["parent_false_flag_values"] == {name: False for name in report["parent_false_flags"]}, "Parent false-flag ledger changed")
    for key in ("source_pins", "protected_file_hashes", "legacy_dirty_path_hashes", "preserved_failed_artifact_hashes_union"):
        require(all(valid_sha(value) for value in report[key].values()), f"Malformed parent {key}")
    return {"report": report, "result": result, "report_raw": report_raw, "result_raw": result_raw, "report_meta": artifact_from_bytes(PARENT_REPORT, report_raw, embed=False), "result_meta": artifact_from_bytes(PARENT_RESULT, result_raw, embed=False), "report_view": parent_report_view(report)}
def load_api(path: Path = API_PATH, expected: str | None = None) -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    raw, api = read_json(path, expected)
    require(api.get("schema") == "native_O3_cylindrical_pullback_api_v1", "Kernel API schema changed"); require(api.get("component") == COMPONENT, "Kernel API component changed"); require(api.get("artifact", {}).get("status") == "draft_unexecuted" and api["artifact"].get("runtime_model_metadata") is None, "Kernel API is not an inert draft contract"); require(api.get("input_contract", {}).get("schema") == INPUT_SCHEMA, "Kernel API input schema changed")
    result_contract = api.get("result_contract", {})
    require(result_contract.get("schema") == RESULT_SCHEMA and result_contract.get("component") == COMPONENT, "Kernel API result schema changed"); require(result_contract.get("metadata_only") is False and result_contract.get("source_only") is True, "Kernel API result markers changed")
    diagnostic = api.get("diagnostic", {})
    events = diagnostic.get("success_event_sequence")
    multiplicity = diagnostic.get("success_event_multiplicity")
    require(isinstance(events, list) and all(isinstance(item, str) for item in events), "Kernel callback event names missing"); require(isinstance(multiplicity, dict) and set(multiplicity) == set(events) and sum(multiplicity.values()) == 23, "Kernel callback event contract changed"); require(result_contract.get("success_counts_exact") == EXPECTED_COUNTS, "Kernel structural count contract changed"); require(tuple(result_contract.get("success_flags_true", ())) == TRUE_FLAGS and tuple(result_contract.get("explicit_false_flags", ())) == FALSE_FLAGS, "Kernel claim flag contract changed")
    plan = api.get("symbolic_call_plan", {})
    symbolic_expected = {"diff_api_calls": plan.get("diff_api_calls", {}).get("planned"), "determinant_api_calls": plan.get("determinant_api_calls"), "reduction_wrapper_calls": plan.get("reduction_wrapper_calls", {}).get("planned")}
    require(symbolic_expected == {"diff_api_calls": 70, "determinant_api_calls": 1, "reduction_wrapper_calls": 14}, "Kernel symbolic API count plan changed")
    source_contract = api.get("input_contract", {}).get("saved_sources", {})
    expected_source = {"accepted_result": PARENT_RESULT_SHA256, **{key: value[1] for key, value in FACT_SPECS.items()}}
    require(set(source_contract) >= set(expected_source), "Kernel source pins are incomplete")
    require({key: source_contract[key].get("raw_sha256") for key in expected_source} == expected_source, "Kernel source pins changed")
    return raw, api, {"event_names": list(events), "event_counts": dict(multiplicity), "symbolic_counts": symbolic_expected}
def load_facts() -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    values: dict[str, Any] = {}
    artifacts: dict[str, dict[str, Any]] = {}
    for name, (path, digest, schema) in FACT_SPECS.items():
        raw, value = read_json(path, digest)
        require(value.get("schema") == schema and value.get("runtime_model_metadata") is None, f"Static fact changed: {name}")
        values[name] = value
        artifacts[name] = artifact_from_bytes(path, raw, embed=True)
    require(values["rm_trace"].get("status") == "bounded_static_source_ast_literal_trace", "Rm trace status changed"); require(values["phase_trace"].get("status") == "passed_static_source_ast_trace", "Phase trace status changed"); require(values["overlap_audit"].get("audit_status") == "passed_static_overlap_audit", "Overlap audit status changed")
    return values, artifacts
def literal_wrapper(literal: Any, path: str | Path, raw: bytes) -> dict[str, Any]:
    return {"literal": clone(literal), "raw_sha256": sha_bytes(raw), "path": wire(path), "bytes": len(raw), "data_only": True}
def branch_manifests(parent_result: Mapping[str, Any]) -> list[dict[str, Any]]:
    manifests = parent_result.get("payload", {}).get("native_bundle_manifests")
    require(isinstance(manifests, list) and len(manifests) == 2, "Parent native branch manifests are missing")
    result: list[dict[str, Any]] = []
    for item in manifests:
        branch = item.get("branch")
        frame = item.get("frame")
        require(isinstance(item.get("manifest_id"), str) and item["manifest_id"], "Saved manifest_id is required")
        require(branch in {"background", "modified"} and isinstance(frame, dict), "Parent branch frame changed")
        require(frame.get("chart") == "O3_transition" and frame.get("coordinate_node") == 241536 and frame.get("radius_node") == 241787, "Native frame contract changed")
        result.append({"manifest_id": item["manifest_id"], "branch": branch, "frame": clone(frame), "radius_node": frame["radius_node"], "radius_definition": clone(item.get("radius_definition")), "radius_handle": {"owner": frame.get("native_Rm_owner", "self.Rm = self.leading.inlet['radius']")}})
    require({item["branch"] for item in result} == {"background", "modified"}, "Both branch manifests are required")
    return result
def build_input(parent: Mapping[str, Any], facts: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    saved_sources: dict[str, Any] = {"accepted_result": literal_wrapper(parent["result"], PARENT_RESULT, parent["result_raw"])}
    artifacts: dict[str, dict[str, Any]] = {}
    for name, value in facts.items():
        path, digest, _ = FACT_SPECS[name]
        raw = read_bytes(path, digest)
        saved_sources[name] = literal_wrapper(value, path, raw)
        artifacts[name] = artifact_from_bytes(path, raw, embed=True)
    input_value = {"schema": INPUT_SCHEMA, "component": COMPONENT, "parent_report_raw_sha256": PARENT_REPORT_SHA256, "parent_result_raw_sha256": PARENT_RESULT_SHA256, "saved_sources": saved_sources, "branch_manifests": branch_manifests(parent["result"]), "requested_context": {"held_phase_rows_as_total_x": True, "physical_measure_as_3d_volume": True, "rho_Z_zero_specialization": True}}; return input_value, artifacts
def update_main_brief(patch: Mapping[str, Any], expected_rules: Sequence[str] = EXECUTION_RULES) -> dict[str, Any]:
    raw, brief = read_json(MAIN_BRIEF)
    require(brief.get("execution_rules") == list(expected_rules), "Main brief execution_rules changed")
    preimage = brief_preimage(raw)
    preimage_path = Path(preimage["path"])
    if not preimage_path.exists():
        write_new(preimage_path, raw)
    for key, value in patch.items():
        if key == "in_progress" and isinstance(brief.get(key), dict) and isinstance(value, Mapping):
            merged = clone(brief[key])
            merged.update(clone(value))
            brief[key] = merged
        else:
            brief[key] = clone(value)
    brief["execution_rules"] = list(expected_rules)
    brief["execution_rules_preserved_exact"] = True
    new_raw = encode(brief_safe(brief))
    require(len(new_raw) < 100_000, "Updated main brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, new_raw)
    return {"path": wire(MAIN_BRIEF), "sha256": sha_bytes(new_raw), "bytes": len(new_raw), "preimage": {key: value for key, value in preimage.items() if key != "raw_bytes_base64"}}
def build_source(parent: Mapping[str, Any], input_value: Mapping[str, Any], input_raw: bytes, api_raw: bytes, api_contract: Mapping[str, Any], facts_artifacts: Mapping[str, Mapping[str, Any]], kernel_raw: bytes, runner_raw: bytes, brief_preimage_meta: Mapping[str, Any]) -> dict[str, Any]:
    additions = {wire(INPUT_SPEC): sha_bytes(input_raw), wire(API_PATH): sha_bytes(api_raw), wire(KERNEL_PATH): sha_bytes(kernel_raw), wire(RUNNER_PATH): sha_bytes(runner_raw)}
    additions.update({item["path"]: item["sha256"] for item in facts_artifacts.values()})
    source_union = merge_hashes(parent["report"]["source_pins"], additions)
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "parent_report_raw_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "source_pins": clone(parent["report"]["source_pins"]),
        "protected_file_hashes": clone(parent["report"]["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(parent["report"]["legacy_dirty_path_hashes"]),
        "preserved_failed_artifact_hashes_union": clone(parent["report"]["preserved_failed_artifact_hashes_union"]),
        "parent_false_flags": clone(parent["report"]["parent_false_flags"]),
        "parent_false_flag_values": clone(parent["report"]["parent_false_flag_values"]),
        "preserved_error_cap_metadata": clone(parent["report"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(parent["report"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(parent["report"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(parent["report"]["retained_cap_provenance"]),
        "input_literal": clone(dict(input_value)),
        "input_raw_sha256": sha_bytes(input_raw),
        "input_raw_bytes_base64": base64.b64encode(input_raw).decode("ascii"),
        "input_spec": {"path": wire(INPUT_SPEC), "sha256": sha_bytes(input_raw), "bytes": len(input_raw)},
        "api": artifact_from_bytes(API_PATH, api_raw, embed=True),
        "code2": {"kernel": artifact_from_bytes(KERNEL_PATH, kernel_raw, embed=True), "runner": artifact_from_bytes(RUNNER_PATH, runner_raw, embed=True)},
        "static_fact_artifacts": clone(dict(facts_artifacts)),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": {name: True for name in TRUE_FLAGS} | {name: False for name in FALSE_FLAGS},
        "callback_event_names": clone(api_contract["event_names"]),
        "callback_event_counts": clone(api_contract["event_counts"]),
        "symbolic_api_counts": clone(api_contract["symbolic_counts"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "main_brief_preimage": clone(dict(brief_preimage_meta)),
        "no_project_imports_before_run": True,
        "no_provider_calls": True,
    }
def load_source(*, require_empty_journal: bool = False) -> tuple[bytes, dict[str, Any]]:
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "Source checkpoint schema changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    input_raw = base64.b64decode(source["input_raw_bytes_base64"])
    require(sha_bytes(input_raw) == source["input_raw_sha256"] and json.loads(input_raw.decode("utf-8")) == source["input_literal"], "Frozen input literal changed")
    read_bytes(INPUT_SPEC, source["input_raw_sha256"])
    read_bytes(API_PATH, source["api"]["sha256"])
    for item in source["code2"].values():
        read_bytes(resolve_path(item["path"], REPO), item["sha256"])
    for item in source["static_fact_artifacts"].values():
        read_bytes(resolve_path(item["path"], WORK), item["sha256"])
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    if require_empty_journal:
        require(JOURNAL_CHECKPOINT.exists() and JOURNAL_CHECKPOINT.read_bytes() == b"", "Run requires an empty durable journal")
    return raw, source
def callback_summary(events: Sequence[Mapping[str, Any]], contract: Mapping[str, Any]) -> dict[str, Any]:
    require(len(events) >= 2, "Manager journal prefix is incomplete")
    require(events[0].get("event") == "manager_entry_started", "Manager entry prefix changed")
    require(events[1].get("event") in {"manager_kernel_imported", "manager_kernel_or_derive_failed"}, "Manager import prefix changed")
    rows = [dict(item) for item in events[2:]]
    counts: dict[str, int] = {}
    for row in rows:
        name = row.get("event")
        require(name in contract["event_names"], f"Unexpected kernel callback event: {name}")
        counts[name] = counts.get(name, 0) + 1
    return {"manager_prefix_event_count": 2, "callback_event_count": len(rows), "callback_event_counts": counts, "callback_event_raw_rows": rows, "journal_event_count": len(events), "journal_event_count_expected": 25}
def _zero_counts() -> dict[str, int]:
    return {key: 0 for key in EXPECTED_COUNTS}
def _false_flags() -> dict[str, bool]:
    return {name: False for name in ALL_FLAGS}
def _recovered_counts(rows: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    event_counts = {
        "branch_frame_binding": "branch_frame_bindings",
        "algebra_family_entry": "new_algebra_family_calls",
        "jacobian_identity_row": "jacobian_identity_rows",
        "inverse_coordinate_row": "inverse_coordinate_rows",
        "cartesian_commutator_row": "cartesian_commutator_rows",
        "scalar_laplacian_row": "scalar_laplacian_rows",
        "branch_operator_manifest": "branch_operator_manifests",
        "policy_rejection": "policy_rejections",
    }
    counts = _zero_counts()
    witnesses = []
    for row in rows:
        event = row.get("event")
        name = event_counts.get(event)
        if name:
            counts[name] += 1
        if event in {"jacobian_identity_row", "inverse_coordinate_row", "cartesian_commutator_row", "scalar_laplacian_row"}:
            witnesses.append(row)
    counts["new_mathematical_proofs"] = int(len(witnesses) == 14 and all(row.get("proved") is True for row in witnesses))
    return counts
def validate_result(result: Mapping[str, Any], source: Mapping[str, Any], api_contract: Mapping[str, Any], events: Sequence[Mapping[str, Any]], derive_entry_calls: int, derive_calls: int) -> dict[str, Any]:
    component = result.get("component_result") if isinstance(result.get("component_result"), Mapping) else result
    violations: list[dict[str, Any]] = []
    for key, expected in (("schema", RESULT_SCHEMA), ("component", COMPONENT), ("metadata_only", False), ("source_only", True), ("conditional", True), ("physical_layer_completed", False), ("full_goal_complete", False), ("runtime_model_metadata", None)):
        if component.get(key) != expected:
            violations.append({"field": key, "expected": expected, "actual": component.get(key)})
    counts = component.get("actual_counts")
    flags = component.get("new_flags")
    if not isinstance(counts, dict) or set(counts) != set(EXPECTED_COUNTS):
        violations.append({"field": "actual_counts", "reason": "exact structural count keys are required"})
    else:
        for key in EXPECTED_COUNTS:
            if not isinstance(counts[key], int) or isinstance(counts[key], bool) or counts[key] < 0:
                violations.append({"field": f"actual_counts.{key}", "reason": "count is not a nonnegative integer"})
    if not isinstance(flags, dict) or set(flags) != set(ALL_FLAGS) or not all(isinstance(flags.get(name), bool) for name in ALL_FLAGS):
        violations.append({"field": "new_flags", "reason": "exact claim flag keys are required"})
    elif component.get("passed") is True:
        if counts != EXPECTED_COUNTS or any(flags.get(name) is not True for name in TRUE_FLAGS) or any(flags.get(name) is not False for name in FALSE_FLAGS):
            violations.append({"field": "success_contract", "reason": "success counts or flags do not match the API"})
    elif isinstance(flags, dict) and any(flags.get(name) is not False for name in ALL_FLAGS):
        violations.append({"field": "failure_flags", "reason": "failed result flags must all be false"})
    payload = component.get("payload")
    if component.get("passed") is True:
        symbolic = payload.get("symbolic_call_counts") if isinstance(payload, dict) else None
        if not isinstance(symbolic, dict):
            violations.append({"field": "payload.symbolic_call_counts", "reason": "symbolic API counts missing"})
        else:
            for key, expected in api_contract["symbolic_counts"].items():
                if symbolic.get(key) != expected:
                    violations.append({"field": f"payload.symbolic_call_counts.{key}", "expected": expected, "actual": symbolic.get(key)})
    try:
        summary = callback_summary(events, api_contract)
    except AssertionError as exc:
        summary = {"manager_prefix_event_count": 0, "callback_event_count": 0, "callback_event_counts": {}, "callback_event_raw_rows": [], "journal_event_count": len(events)}
        violations.append({"field": "journal", "reason": str(exc)})
    if component.get("passed") is True:
        if derive_entry_calls != 1 or derive_calls != 1:
            violations.append({"field": "derive_calls", "reason": "successful entry must invoke derive exactly once"})
        if summary.get("callback_event_count") != 23 or summary.get("callback_event_counts") != api_contract["event_counts"]:
            violations.append({"field": "callback_events", "reason": "callback multiplicity does not match API"})
    return {"schema": COMPONENT + "_result_validation_v1", "component": COMPONENT, "passed": component.get("passed") is True and not violations, "explicit_passed": component.get("passed"), "derive_entry_calls": derive_entry_calls, "derive_calls": derive_calls, "actual_counts": clone(counts), "new_flags": clone(flags), "callback": summary, "violations": violations, "raw_result_preserved": True}
def import_kernel(path: Path) -> Any:
    require(path.resolve() != Path(__file__).resolve(), "Kernel path aliases metadata manager")
    module_name = "_native_O3_cylindrical_pullback_entry"
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, "Cannot load the frozen kernel")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "Frozen kernel has no derive API")
    return module
def prepare(args: argparse.Namespace) -> dict[str, Any]:
    require(not any(path.exists() for path in (INPUT_SPEC, BRIEF, PREPARE_RECEIPT, RUN_RECEIPT, PUBLISH_RECEIPT, *CHECKPOINTS)), "This component already has a durable attempt")
    parent = load_parent()
    current_raw, _ = load_current_brief()
    api_raw, _, api_contract = load_api(API_PATH, args.api_sha256.lower() if args.api_sha256 else None)
    facts, fact_artifacts = load_facts()
    kernel_raw = read_bytes(KERNEL_PATH, args.kernel_sha256.lower() if args.kernel_sha256 else None)
    runner_raw = read_bytes(RUNNER_PATH, args.runner_sha256.lower() if args.runner_sha256 else None)
    input_value, _ = build_input(parent, facts)
    input_raw = write_new(INPUT_SPEC, input_value)
    preimage = brief_preimage(current_raw)
    preimage_path = Path(preimage["path"])
    require(not preimage_path.exists(), "Main brief preimage already exists; metadata repair is external")
    write_new(preimage_path, current_raw)
    source = build_source(parent, input_value, input_raw, api_raw, api_contract, fact_artifacts, kernel_raw, runner_raw, preimage)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    journal_raw = write_new(JOURNAL_CHECKPOINT, b"")
    pending = {"schema": COMPONENT + "_brief_v1", "component": COMPONENT, "phase": "prepared", "source_checkpoint": {"path": wire(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)}, "input_spec": {"path": wire(INPUT_SPEC), "sha256": sha_bytes(input_raw)}, "api": artifact_from_bytes(API_PATH, api_raw, embed=False), "kernel": source["code2"]["kernel"], "runner": source["code2"]["runner"], "journal": {"path": wire(JOURNAL_CHECKPOINT), "sha256": sha_bytes(journal_raw)}, "metadata_only": True, "source_only": True, "do_not_reenter": True}
    pending_raw = write_new(BRIEF, pending)
    main_brief = update_main_brief({"phase": "prepared", "current_component": COMPONENT, "in_progress": {COMPONENT: {"brief": wire(BRIEF), "source": wire(SOURCE_CHECKPOINT), "source_sha256": sha_bytes(source_raw), "metadata_only": True, "source_only": True}}, "goal_active": True, "full_goal_complete": False, "next_task": "Run the frozen native O3 cylindrical pullback kernel exactly once after root entry."})
    receipt = {"schema": COMPONENT + "_prepare_receipt_v1", "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "input_spec_sha256": sha_bytes(input_raw), "journal_checkpoint_sha256": sha_bytes(journal_raw), "brief_sha256": sha_bytes(pending_raw), "main_brief": main_brief, "parent_report_sha256": PARENT_REPORT_SHA256, "parent_result_sha256": PARENT_RESULT_SHA256, "source_union_count": source["source_union_count"], "metadata_only": True, "source_only": True, "runtime_model_metadata": None}
    receipt_raw = write_new(PREPARE_RECEIPT, receipt)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw)}
def run(args: argparse.Namespace) -> dict[str, Any]:
    source_raw, source = load_source(require_empty_journal=True)
    require(PREPARE_RECEIPT.exists(), "Prepare receipt is missing")
    _, prepare_receipt = read_json(PREPARE_RECEIPT)
    require(prepare_receipt.get("source_checkpoint_sha256") == sha_bytes(source_raw), "Prepare receipt does not authenticate source")
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not RUN_RECEIPT.exists(), "Run has already been entered")
    entered = {"schema": ENTERED_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "kernel_path": wire(KERNEL_PATH), "kernel_sha256": source["code2"]["kernel"]["sha256"], "import_attempts": 1, "import_attempt_recorded_before_import": True, "sys_modules_registration_recorded_before_import": True, "derive_entry_calls": 0, "derive_entry_reserved": True, "kernel_imported": False, "runtime_model_metadata": None, "entered_at_utc": datetime.now(timezone.utc).isoformat(), "do_not_reenter": True}
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    append_journal({"schema": DIAGNOSTIC_SCHEMA, "event": "manager_entry_started", "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw), "import_attempts": 1, "derive_entry_calls": 0, "derive_entry_reserved": True})
    api_contract = {"event_names": source["callback_event_names"], "event_counts": source["callback_event_counts"], "symbolic_counts": source["symbolic_api_counts"]}
    kernel_imported = False
    derive_entry_calls = 0
    derive_calls = 0
    component_result: dict[str, Any] | None = None
    try:
        module = import_kernel(KERNEL_PATH)
        kernel_imported = True
        append_journal({"schema": DIAGNOSTIC_SCHEMA, "event": "manager_kernel_imported", "kernel_sha256": source["code2"]["kernel"]["sha256"], "sys_modules_registered": True})
        def callback(row: Mapping[str, Any]) -> None:
            require(isinstance(row, Mapping), "Kernel callback row must be a mapping")
            require(row.get("event") in source["callback_event_names"], "Kernel callback event is not in the frozen API")
            append_journal(dict(row))
        derive_entry_calls = 1
        derive_calls = 1
        component_result = module.derive(clone(source["input_literal"]), callback)
        require(isinstance(component_result, dict), "Kernel derive did not return a JSON mapping")
    except BaseException as exc:
        append_journal({"schema": DIAGNOSTIC_SCHEMA, "event": "manager_kernel_or_derive_failed", "error_type": type(exc).__name__, "error": repr(exc), "kernel_imported": kernel_imported, "derive_calls": derive_calls})
        rows = journal_events()[2:]
        existing = component_result if isinstance(component_result, dict) else {}
        actual_counts = existing.get("actual_counts") if isinstance(existing.get("actual_counts"), dict) and set(existing["actual_counts"]) == set(EXPECTED_COUNTS) else _recovered_counts(rows)
        component_result = clone(existing) if existing else {"schema": RESULT_SCHEMA, "component": COMPONENT, "runtime_model_metadata": None, "metadata_only": False, "source_only": True, "conditional": True, "physical_layer_completed": False, "full_goal_complete": False, "payload": {"recovered_callback_prefix": clone(rows)}, "diagnostics": {"event_count": len(rows), "callback_failure_stops_later_work": True}}
        component_result.update({"passed": False, "schema": RESULT_SCHEMA, "component": COMPONENT, "runtime_model_metadata": None, "metadata_only": False, "source_only": True, "conditional": True, "physical_layer_completed": False, "full_goal_complete": False, "actual_counts": clone(actual_counts), "new_flags": _false_flags(), "failure": {"kind": "manager_import_or_derive_failure", "error_type": type(exc).__name__, "error": repr(exc), "stopped": True}})
    outer = clone(component_result)
    outer["component_result"] = clone(component_result)
    outer["component_result_exact"] = True
    outer["derive_entry_calls"] = derive_entry_calls
    outer["derive_calls"] = derive_calls
    outer["kernel_imported"] = kernel_imported
    outer["execution_metadata"] = {"derive_entry_calls": derive_entry_calls, "derive_calls": derive_calls, "kernel_imported": kernel_imported, "import_attempts": 1, "old_math_reexecuted": False, "runner_executed": False}
    outer["passed"] = bool(component_result.get("passed") is True)
    outer["metadata_only"] = False
    outer["source_only"] = True
    outer["new_flags"] = clone(component_result.get("new_flags", _false_flags())) if outer["passed"] else _false_flags()
    outer["actual_counts"] = clone(component_result.get("actual_counts", _zero_counts()))
    result_raw = write_new(RESULT_CHECKPOINT, outer)
    validation = validate_result(outer, source, api_contract, journal_events(), derive_entry_calls, derive_calls)
    receipt = {"schema": COMPONENT + "_run_receipt_v1", "component": COMPONENT, "entered_checkpoint_sha256": sha_bytes(entered_raw), "result_checkpoint_sha256": sha_bytes(result_raw), "journal_checkpoint_sha256": sha_bytes(JOURNAL_CHECKPOINT.read_bytes()), "derive_entry_calls": derive_entry_calls, "derive_calls": derive_calls, "kernel_imported": kernel_imported, "durable_result_saved_before_validation": True, "callback_prefix_recovered_from_journal": True, "validation": validation, "attempt_status": "passed" if validation["passed"] else "failed", "metadata_only": False, "source_only": True, "runtime_model_metadata": None, "do_not_reenter": True}
    receipt_raw = write_new(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}
def load_review(path: str, digest: str) -> dict[str, Any]:
    require(valid_sha(digest), "Saved review SHA is malformed")
    review_path = resolve_path(path, WORK)
    raw, value = read_json(review_path, digest.lower())
    require(value.get("passed") is True and value.get("independent_review") is True, "Saved review is not independent and passed")
    require(value.get("math_source_identity_accepted") is True and value.get("mathematical_proof_accepted") is True, "Saved review did not accept the conditional mathematical source identity")
    scope = value.get("conditional_geometry_scope") is True or value.get("scope") in {"conditional_native_geometry", "conditional_geometry", "native_cylindrical_chart_metric"} or value.get("review_scope") in {"conditional_native_geometry", "conditional_geometry", "native_cylindrical_chart_metric"}
    require(scope, "Saved review is not explicitly scoped to conditional geometry")
    return {"path": wire(review_path), "sha256": sha_bytes(raw), "bytes": len(raw), "schema": value.get("schema"), "passed": True, "independent_review": True, "math_source_identity_accepted": True, "mathematical_proof_accepted": True, "conditional_geometry_scope": True, "raw_bytes_exact": True, "raw_readonly": True, "outside_science_pins": True, "raw_bytes_base64": base64.b64encode(raw).decode("ascii"), "value": value}
def attrs_candidate() -> bytes:
    path = REPO / ".gitattributes"
    old = path.read_bytes() if path.exists() else b""
    text = old.decode("utf-8")
    lines = [f"/{KERNEL_RELATIVE.as_posix()} -text", f"/{RUNNER_RELATIVE.as_posix()} -text", f"/{REPORT_RELATIVE.as_posix()} -text", f"/{DOC_RELATIVE.as_posix()} -text"]
    existing = text.splitlines()
    require(all(existing.count(line) <= 1 for line in lines), "Duplicate scope attribute already exists")
    for line in lines:
        if line not in existing:
            text += ("\n" if text and not text.endswith("\n") else "") + line + "\n"
            existing.append(line)
    return text.encode("utf-8")
def document_bytes(validation: Mapping[str, Any], result: Mapping[str, Any]) -> bytes:
    counts = result.get("actual_counts")
    flags = result.get("new_flags")
    component = result.get("component_result") if isinstance(result.get("component_result"), Mapping) else {}
    symbolic = component.get("payload", {}).get("symbolic_call_counts", {}) if isinstance(component.get("payload"), Mapping) else {}
    lines = [
        "# Current native O3 cylindrical pullback (2026-10-27)",
        "",
        f"Component: `{COMPONENT}`.",
        f"Attempt status: `{'passed' if validation.get('passed') else 'failed'}`; result is conditional geometry, source-only, and not physical completion.",
        "",
        f"Actual structural counts: `{counts}`.",
        f"Actual generic symbolic API counts: `{symbolic}`.",
        f"New flags: `{flags}`.",
        "",
        "The proved geometry is conditional on the universal cylindrical embedding `R = rho(Z) * exp(c + x)` with Jacobian determinant `R^2`, nine inverse-coordinate rows, three Cartesian commutators, and one scalar-Laplacian identity. Native `Rm` remains the unresolved symbolic `rho(Z)` handle; no `D_Z Rm = 0` claim is made.",
        "",
        "Mixed native-field derivatives are not supplied by the saved packets. The phase recipe remains opaque, fractional-part differentiation and phase identity are not claimed, and no physical-time map, physical volume-L2 result, or field-value pullback is installed.",
        "",
        "All inherited source pins, failed artifacts, error caps, construction data, integral inventories, and cap provenance remain retained. No caps are deleted and no global closure or temporal recursion is claimed.",
        "",
        "Next action: source_native_Rm_Ra_phase_parameter_Z_independence_certificate in a fresh scoped task; preserve unknown native Rm/Z and phase/time until that certificate; do not reenter this attempt.",
        "",
    ]
    return "\n".join(lines).encode("utf-8")
def publish(args: argparse.Namespace) -> dict[str, Any]:
    source_raw, source = load_source()
    report_path = REPO / REPORT_RELATIVE
    require(RUN_RECEIPT.exists() and not report_path.exists() and not PUBLISH_RECEIPT.exists(), "Publication already exists or run receipt is missing")
    _, run_receipt = read_json(RUN_RECEIPT)
    result_raw, result = read_json(RESULT_CHECKPOINT, run_receipt.get("result_checkpoint_sha256"))
    api_raw, _, api_contract = load_api(API_PATH, source["api"]["sha256"])
    validation = validate_result(result, source, api_contract, journal_events(), run_receipt.get("derive_entry_calls"), run_receipt.get("derive_calls"))
    review = load_review(args.review, args.review_sha256)
    kernel_artifact = artifact(KERNEL_PATH, source["code2"]["kernel"]["sha256"], embed=True)
    runner_artifact = artifact(RUNNER_PATH, source["code2"]["runner"]["sha256"], embed=True)
    input_artifact = artifact(INPUT_SPEC, source["input_spec"]["sha256"], embed=True)
    api_artifact = artifact_from_bytes(API_PATH, api_raw, embed=True)
    facts_artifacts = {name: artifact(path, digest, embed=True) for name, (path, digest, _) in FACT_SPECS.items()}
    checkpoint_artifacts = [artifact(path, embed=True) for path in CHECKPOINTS]
    helper_artifact = artifact(Path(__file__), embed=True)
    require(not (REPO / DOC_RELATIVE).exists(), "Publication document already exists; metadata repair is external")
    handoff_old: dict[str, bytes] = {}
    handoff_preimages: dict[str, dict[str, Any]] = {}
    handoff_artifacts: dict[str, dict[str, Any]] = {}
    prefix = f"\n\n## {COMPONENT} publication\n\nReport: `{REPORT_RELATIVE.as_posix()}`.\nStatus: `{'passed' if validation['passed'] else 'failed'}`; conditional geometry only; do not reenter.\n"
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        require(f"## {COMPONENT} publication" not in old.decode("utf-8", errors="ignore"), f"Existing handoff already contains {COMPONENT} publication")
        preimage_path = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(old)[:12]}.md"
        require(not preimage_path.exists(), "Handoff preimage already exists; metadata repair is external")
        new = prefix.encode("utf-8") + old
        handoff_old[wire(relative)] = old
        handoff_preimages[wire(relative)] = artifact_from_bytes(preimage_path, old, embed=False)
        handoff_artifacts[wire(relative)] = artifact_from_bytes(path, new, embed=False)
    doc_raw = document_bytes(validation, result)
    attrs_raw = attrs_candidate()
    current_brief_raw, current_brief = read_json(MAIN_BRIEF)
    require(current_brief.get("execution_rules") == EXECUTION_RULES and current_brief.get("goal_active") is True and current_brief.get("full_goal_complete") is False, "Main brief state changed before publication")
    publication_brief_preimage = brief_preimage(current_brief_raw)
    checkpoints = {item["path"]: item["sha256"] for item in checkpoint_artifacts}
    published_pins = merge_hashes(source["source_union_hashes"], checkpoints)
    failed = clone(source["preserved_failed_artifact_hashes_union"])
    if not validation["passed"]:
        failed = merge_hashes(failed, {item["path"]: item["sha256"] for item in (kernel_artifact, runner_artifact, input_artifact, *checkpoint_artifacts)})
    report = {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if validation["passed"] else "failed",
        "passed": validation["passed"],
        "metadata_only": False,
        "publication_metadata_only": True,
        "source_only": True,
        "conditional_geometry_scope": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "prepublication_base_head": CURRENT_HEAD_COMMIT,
        "parent_report": source["parent_report"],
        "parent_result": source["parent_result"],
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_raw_sha256": PARENT_RESULT_SHA256,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "checkpoint_count": PUBLISHED_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "preserved_failed_artifact_count": len(failed),
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "source_pin_count": len(published_pins),
        "source_union_count": len(published_pins),
        "source_pins": published_pins,
        "raw_sha_pins": published_pins,
        "checkpoint_hashes": checkpoints,
        "preserved_failed_artifact_hashes_union": failed,
        "protected_file_hashes": clone(source["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(source["legacy_dirty_path_hashes"]),
        "parent_false_flags": clone(source["parent_false_flags"]),
        "parent_false_flag_values": clone(source["parent_false_flag_values"]),
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "preserved_error_cap_metadata": clone(source["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(source["retained_construction_data"]),
        "retained_integral_term_inventory": clone(source["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(source["retained_cap_provenance"]),
        "input": clone(source["input_literal"]),
        "input_spec": input_artifact,
        "code2": {"kernel": kernel_artifact, "runner": runner_artifact},
        "api_artifact": api_artifact,
        "static_fact_artifacts": facts_artifacts,
        "checkpoint_artifacts": checkpoint_artifacts,
        "current_review": review,
        "component_brief": artifact(BRIEF, embed=True),
        "main_brief_preimage": publication_brief_preimage,
        "metadata_helper_artifact": helper_artifact,
        "document_artifact": artifact_from_bytes(REPO / DOC_RELATIVE, doc_raw, embed=True),
        "handoff_preimages": handoff_preimages,
        "handoff_artifacts": handoff_artifacts,
        "gitattributes_required": {"paths": [wire(REPO / path) for path in (KERNEL_RELATIVE, RUNNER_RELATIVE, REPORT_RELATIVE, DOC_RELATIVE)], "rule": "-text", "apply_before_freeze": True},
        "full_result": clone(result),
        "full_result_raw_sha256": sha_bytes(result_raw),
        "result_checkpoint": {"path": wire(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw), "bytes": len(result_raw)},
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": checkpoints[wire(ENTERED_CHECKPOINT)],
        "diagnostic_checkpoint_sha256": checkpoints[wire(JOURNAL_CHECKPOINT)],
        "actual_counts": validation["actual_counts"],
        "new_flags": validation["new_flags"],
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": {name: True for name in TRUE_FLAGS} | {name: False for name in FALSE_FLAGS},
        "callback_event_counts": clone(source["callback_event_counts"]),
        "callback_event_count": 23,
        "manager_prefix_event_count": 2,
        "journal_event_count": 25,
        "durable_result_validation": validation,
        "report_scope_attributes": list(REPORT_SCOPE_ATTRIBUTES),
        "read_first": [wire(REPO / "AGENTS.md"), wire(REPO / DOC_RELATIVE), wire(REPO / REPORT_RELATIVE)],
        "current_doc_sha256": sha_bytes(doc_raw),
        "current_document_sha256": sha_bytes(doc_raw),
        "runtime_model_metadata": None,
        "do_not_reenter": True,
        "next_executable_task": "source_native_Rm_Ra_phase_parameter_Z_independence_certificate",
        "lossless_artifacts": [kernel_artifact, runner_artifact, input_artifact, api_artifact, *facts_artifacts.values(), *checkpoint_artifacts, review, helper_artifact],
        "raw_bytes_lossless": True,
    }
    report_raw = write_new(report_path, report)
    write_new(REPO / DOC_RELATIVE, doc_raw)
    for relative, old in handoff_old.items():
        path = REPO / relative
        write_new(Path(handoff_preimages[relative]["path"]), old)
        write_replace(path, prefix.encode("utf-8") + old)
    write_replace(REPO / ".gitattributes", attrs_raw)
    brief = update_main_brief({"phase": "published_passed" if validation["passed"] else "published_failed", "current_component": COMPONENT, "report": wire(report_path), "report_sha256": sha_bytes(report_raw), "current_doc_sha256": sha_bytes(doc_raw), "current_document_sha256": sha_bytes(doc_raw), "read_first": report["read_first"], "goal_active": True, "full_goal_complete": False, "metadata_only": current_brief.get("metadata_only", True), "source_only": True, "next_task": report["next_executable_task"], "do_not_reenter": True})
    receipt = {"schema": COMPONENT + "_publication_receipt_v1", "component": COMPONENT, "report": wire(report_path), "report_sha256": sha_bytes(report_raw), "attempt_status": report["attempt_status"], "passed": report["passed"], "source_checkpoint_sha256": sha_bytes(source_raw), "result_checkpoint_sha256": sha_bytes(result_raw), "current_review": review, "main_resume_brief": brief, "checkpoint_count": PUBLISHED_CHECKPOINT_COUNT, "source_pin_count": len(published_pins), "preserved_failed_artifact_count": len(failed), "runtime_model_metadata": None, "do_not_reenter": True}
    receipt_raw = write_new(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}
def parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    prepare_parser = sub.add_parser("prepare")
    prepare_parser.add_argument("--api-sha256")
    prepare_parser.add_argument("--kernel-sha256")
    prepare_parser.add_argument("--runner-sha256")
    sub.add_parser("run")
    publish_parser = sub.add_parser("publish")
    publish_parser.add_argument("--review", required=True)
    publish_parser.add_argument("--review-sha256", required=True)
    return parser
def main(argv: Sequence[str] | None = None) -> dict[str, Any]:
    args = parser().parse_args(argv)
    value = {"prepare": prepare, "run": run, "publish": publish}[args.action](args)
    print(json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":")))
    return value
if __name__ == "__main__":
    main()
