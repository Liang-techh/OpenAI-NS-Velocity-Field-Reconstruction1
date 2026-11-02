"""One-entry raw P physical-m15 endpoint-sign classification consumer."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "n1_p_sector_signs_v3"
PREVIOUS_COMPONENT = "n1_p_compensation_cell_sum_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

P_RECEIPT = "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz"
PRIOR_RECEIPT = "experiments/root_st073/n1_w0_chi_ungrouped_v3_receipt.json.gz"
FAILURE_V1 = "experiments/root_st073/n1_p_compensation_cell_sum_v1_failure.json.gz"
FAILURE_V2 = "experiments/root_st073/n1_p_compensation_cell_sum_v2_failure.json.gz"
FEEDBACK_RECEIPT = "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz"
LEAF = "src/openai_ns_reconstruction/n1_p_sector_signs_v3.py"
FIXTURE = "src/openai_ns_reconstruction/n1_p_sector_signs_fixture_v3.py"
RUNNER = "scripts/run_n1_p_sector_signs_v3_once.py"
REQUIRED = (LEAF, FIXTURE, RUNNER, P_RECEIPT, PRIOR_RECEIPT, FAILURE_V1, FAILURE_V2)
RAW_GZIP_INPUTS = frozenset((P_RECEIPT, PRIOR_RECEIPT, FAILURE_V1, FAILURE_V2))
RAW_SHA256 = {
    P_RECEIPT: "c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf",
    PRIOR_RECEIPT: "3123a436957522fd3de3926dccacf6b463c28587a50951df6377337ffbab0821",
    FAILURE_V1: "012e2d410a74ee16eb6eaf6859bac00884a29514c55a9330006ce09d54e421a1",
    FAILURE_V2: "29e6215eb2137317585592acd7aa638cc84f212a1d7654331a247fe6bf172efa",
}
RESULT_SHA256 = {
    P_RECEIPT: "341a09693c0fc46682a6561c37df468ca536e62c59ee7dc944ad9841e16a53b5",
    PRIOR_RECEIPT: "d0f9f8d466828ea9cb02505f99d21500c6cab1966eeb113af22ea6941e5b4ee3",
}
FEEDBACK_RAW_SHA256 = "0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa"
BRANCHES = ("G0_P_W0", "G0_P_DeltaW0", "D_P_W0", "D_P_DeltaW0")
ORDERS, CELLS = (0, 1), tuple(range(8))
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FORBIDDEN_NAMES = {"rm_loop_density_identity.json", "lei_ren_part1_paper_compliant_current_original_rm_loop_density_identity.json"}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def _reject_float_tree(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(label + " contains a host floating-point value")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _reject_float_tree(key, label + ".key")
            _reject_float_tree(item, label)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            _reject_float_tree(item, label)


def encode(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return {"lower_exact_mpf_tuple": list(raw[0]), "upper_exact_mpf_tuple": list(raw[1]), "directed": True}
    if isinstance(value, Mapping):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [encode(item) for item in sorted(value, key=repr)]
    if isinstance(value, float):
        raise TypeError("cannot serialize a host floating-point value")
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def save(path: Path, value: object) -> str:
    raw = (json.dumps(encode(value), sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def _canonical_result_sha(result: Mapping[str, Any]) -> str:
    raw = (json.dumps(result, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")
    return sha(raw)


def _forbidden(path: Path) -> bool:
    return path.name.lower() in FORBIDDEN_NAMES or "continuation" in path.name.lower()


def _source_rows(receipt: Mapping[str, Any], label: str) -> dict[str, str]:
    source = receipt.get("source")
    require(isinstance(source, Mapping), label + " source manifest missing")
    rows = source.get("files")
    require(isinstance(rows, list), label + " source files missing")
    result: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), label + " source row")
        path, digest = row.get("path"), row.get("sha256")
        require(isinstance(path, str) and isinstance(digest, str) and path not in result, label + " source path")
        result[path] = digest
    return result


def preflight() -> tuple[dict[str, str], str, str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw, source_raw = BRIEF.read_bytes(), SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(str(checkpoint.get("path", ""))).resolve() == SOURCE.resolve() and checkpoint.get("sha256") == sha(source_raw), "source checkpoint")
    component, phase = brief.get("component"), brief.get("phase")
    require((component == COMPONENT and phase in (COMPONENT + "_pending", "component_pending")) or (component == PREVIOUS_COMPONENT and phase == PREVIOUS_COMPONENT + "_pending"), "pending predecessor brief guard")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False and brief.get("one_entry_only") is True and brief.get("do_not_reenter") is True, "brief lifecycle guard")
    preimage_path_value, preimage_sha = checkpoint.get("preimage_path"), checkpoint.get("preimage_sha256")
    if preimage_path_value is None and isinstance(brief.get("previous_pending_preimage"), Mapping):
        previous = brief["previous_pending_preimage"]
        preimage_path_value, preimage_sha = previous.get("path"), previous.get("sha256")
    require(isinstance(preimage_path_value, str) and isinstance(preimage_sha, str), "pending preimage pointer")
    preimage_path = Path(preimage_path_value)
    require(not _forbidden(preimage_path), "forbidden preimage path")
    preimage_raw = preimage_path.read_bytes()
    require(sha(preimage_raw) == preimage_sha, "pending preimage SHA")
    preimage = json.loads(preimage_raw)
    require(isinstance(preimage, Mapping) and preimage.get("component") == PREVIOUS_COMPONENT and preimage.get("phase") in (PREVIOUS_COMPONENT + "_pending", "component_pending"), "published predecessor preimage guard")
    manifest = json.loads(source_raw)
    require(isinstance(manifest, Mapping) and manifest.get("component") == COMPONENT and manifest.get("previous_component") == PREVIOUS_COMPONENT and manifest.get("one_entry_only") is True, "source manifest guard")
    policy = manifest.get("source_policy", {})
    require(isinstance(policy, Mapping), "source policy mapping")
    for key in ("authenticate_all_raw_gzip_bytes_before_decode", "fixture_before_saved_receipt_consumption", "previous_receipts_decode_only_after_entry", "complete_result_and_failures_are_saved", "no_ancestor_import_or_reconstruction", "source_bind_after_freeze"):
        require(policy.get(key) is True, "source policy guard " + key)
    rows = manifest.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "source manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "source manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes and not _forbidden(Path(path)), "unexpected source file")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed frozen source file " + path)
        if path in RAW_GZIP_INPUTS:
            require(row.get("raw_sha256_before_decode") is True and digest == RAW_SHA256[path], "changed pinned raw receipt " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "incomplete source manifest")
    return hashes, sha(source_raw), sha(brief_raw), str(preimage_sha)


def _module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module unavailable: " + path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _leaf_callable(module: Any) -> Callable[..., Any]:
    candidate = getattr(module, "build", None)
    require(callable(candidate), "sector-sign leaf must expose build")
    return candidate


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    require(isinstance(value, Mapping), label + " must be a mapping")
    return value


def _exact_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    require(isinstance(value, (list, tuple)) and len(value) == 4 and all(type(item) is int for item in value), label + " exact endpoint tuple")
    return tuple(value)  # type: ignore[return-value]


def _interval_token(value: Any, label: str) -> dict[str, Any]:
    require(isinstance(value, Mapping), label + " coefficient interval")
    return {"lower_exact_mpf_tuple": list(_exact_tuple(value.get("lower_exact_mpf_tuple"), label + ".lower")), "upper_exact_mpf_tuple": list(_exact_tuple(value.get("upper_exact_mpf_tuple"), label + ".upper")), "directed": True}


def _origin(term: Mapping[str, Any], label: str, order: int, cell: int) -> dict[str, Any]:
    raw = term.get("origin")
    require(isinstance(raw, Mapping), label + ".origin mapping")
    required = {"source_ids", "error_classes", "uncertainty", "target_cell", "lead_order", "word_order"}
    require(required.issubset(raw), label + ".origin fields")
    ids, classes = raw["source_ids"], raw["error_classes"]
    require(isinstance(ids, Sequence) and not isinstance(ids, (str, bytes, bytearray)) and bool(ids) and all(isinstance(item, str) and item for item in ids), label + ".origin.source_ids")
    require(isinstance(classes, Sequence) and not isinstance(classes, (str, bytes, bytearray)) and all(isinstance(item, str) for item in classes), label + ".origin.error_classes")
    require(isinstance(raw["uncertainty"], bool) and type(raw["target_cell"]) is int and raw["target_cell"] == cell, label + ".origin metadata")
    lead, word = raw["lead_order"], raw["word_order"]
    require(type(lead) is int and type(word) is int and lead in (0, 1) and word in (0, 1) and lead + word == order, label + ".origin derivative/order")
    copied = dict(raw)
    copied["source_ids"], copied["error_classes"] = list(ids), list(classes)
    return copied


def _potential_terms(output: Mapping[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    potential = output.get("potential_compensation_m15")
    require(isinstance(potential, Mapping) and set(potential) == set(BRANCHES), "P branch roster")
    packet: dict[str, Any] = {}
    provenance: list[dict[str, Any]] = []
    for branch in BRANCHES:
        branch_raw = potential[branch]
        require(isinstance(branch_raw, Mapping), branch + " branch mapping")
        packet[branch] = {}
        for order in ORDERS:
            order_key = str(order)
            row = branch_raw.get(order_key, branch_raw.get(order))
            require(isinstance(row, Mapping) and set(str(key) for key in row) == {str(cell) for cell in CELLS}, branch + " cell roster")
            packet[branch][order_key] = {}
            for cell in CELLS:
                cell_key = str(cell)
                raw_terms = row.get(cell_key, row.get(cell))
                require(isinstance(raw_terms, list), branch + " Z" + str(order) + " cell" + str(cell) + " terms")
                annotated: list[dict[str, Any]] = []
                for index, raw_term in enumerate(raw_terms):
                    require(isinstance(raw_term, Mapping), "P term mapping")
                    label = f"P.{branch}.Z{order}.cell{cell}[{index}]"
                    _reject_float_tree(raw_term, label)
                    powers = raw_term.get("powers")
                    require(isinstance(powers, Sequence) and not isinstance(powers, (str, bytes, bytearray)) and len(powers) == 4 and all(type(power) is int for power in powers), label + ".powers")
                    coefficient = _interval_token(raw_term.get("coefficient"), label + ".coefficient")
                    origin = _origin(raw_term, label, order, cell)
                    error_only = raw_term.get("error_only", origin.get("error_only", False))
                    require(isinstance(error_only, bool), label + ".error_only")
                    source_path = ["result", "output", "potential_compensation_m15", branch, order_key, cell_key, index]
                    term = dict(raw_term)
                    term.update({"branch": branch, "ordinary_Z_order": order, "target_cell": cell, "lead_order": origin["lead_order"], "word_order": origin["word_order"], "powers": list(powers), "coefficient": coefficient, "uncertainty": origin["uncertainty"], "error_only": error_only, "receipt": P_RECEIPT, "source_path": source_path, "input_index": index})
                    term_origin = dict(origin)
                    term_origin.update({"receipt": P_RECEIPT, "source_path": source_path})
                    term["origin"] = term_origin
                    annotated.append(term)
                    provenance.append({"branch": branch, "Z": order_key, "ordinary_Z_order": order, "cell": cell, "target_cell": cell, "lead_order": origin["lead_order"], "word_order": origin["word_order"], "powers": list(powers), "coefficient": coefficient, "receipt": P_RECEIPT, "source_path": source_path, "origin": term_origin, "uncertainty": origin["uncertainty"], "error_only": error_only, "input_index": index})
                packet[branch][order_key][cell_key] = annotated
    return packet, provenance


def _receipt(raw: Mapping[str, bytes], path: str, component: str) -> dict[str, Any]:
    require(sha(raw[path]) == RAW_SHA256[path], "pinned raw receipt changed " + path)
    envelope = json.loads(gzip.decompress(raw[path]))
    require(isinstance(envelope, Mapping), path + " receipt envelope")
    result = envelope.get("result")
    require(isinstance(result, Mapping) and result.get("status") == "returned" and result.get("component") == component, path + " status/component")
    require(_canonical_result_sha(result) == RESULT_SHA256[path], path + " saved result SHA")
    return dict(envelope)


def _source_atom(value: Any, label: str) -> dict[str, list[int]]:
    if isinstance(value, Mapping):
        lower, upper = value.get("lower_exact_mpf_tuple"), value.get("upper_exact_mpf_tuple")
        require(lower is not None and upper is not None, label + " exact source endpoints")
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        lower, upper = value
    else:
        raise TypeError(label + " source atom must expose exact endpoint tuples")
    return {"lower_exact_mpf_tuple": list(_exact_tuple(lower, label + ".lower")), "upper_exact_mpf_tuple": list(_exact_tuple(upper, label + ".upper"))}


def _source_scope(output: Mapping[str, Any], label: str) -> tuple[Any, Mapping[str, Any], dict[str, dict[str, list[int]]], Mapping[str, Any]]:
    family, logs, basis = output.get("source_family"), output.get("source_logs"), output.get("source_power_basis", POWER_BASIS)
    require(family is not None and isinstance(logs, Mapping) and isinstance(basis, Sequence) and not isinstance(basis, (str, bytes, bytearray)) and tuple(basis) == POWER_BASIS, label + " source scope")
    atoms: dict[str, dict[str, list[int]]] = {}
    for atom in POWER_BASIS:
        require(atom in logs, label + " missing source atom " + atom)
        atoms[atom] = _source_atom(logs[atom], label + "." + atom)
    index = output.get("input_index", output.get("source_index", {}))
    require(isinstance(index, Mapping), label + " input index mapping")
    return family, logs, atoms, index


def _same_source_atoms(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    return all(left[atom]["lower_exact_mpf_tuple"] == right[atom]["lower_exact_mpf_tuple"] and left[atom]["upper_exact_mpf_tuple"] == right[atom]["upper_exact_mpf_tuple"] for atom in POWER_BASIS)


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry " + path)
        raw[path] = data
    for path in RAW_GZIP_INPUTS:
        require(sha(raw[path]) == RAW_SHA256[path], "changed pinned raw gzip " + path)
    progress["raw_gzip_authenticated_before_decode"] = True
    leaf_module, fixture_module = _module(LEAF, COMPONENT + "_leaf"), _module(FIXTURE, COMPONENT + "_fixture")
    fixture_runner = getattr(fixture_module, "run_fixture", getattr(fixture_module, "run", None))
    require(callable(fixture_runner), "fixture must export run_fixture")
    progress["fixture_started"] = True
    proof = fixture_runner(leaf_module)
    progress["fixture"] = encode(proof)
    require(isinstance(proof, Mapping) and proof.get("passed") is True, "synthetic fixture failed")
    checks = proof.get("checks")
    if isinstance(checks, Mapping):
        require(all(value is True for value in checks.values()), "synthetic fixture has failed checks")
    progress["fixture_completed_before_receipt_decode"] = True
    progress["input_json_decode_started"] = True
    p_receipt, prior_receipt = _receipt(raw, P_RECEIPT, "n1_core_kchi_compensation_v2"), _receipt(raw, PRIOR_RECEIPT, "n1_w0_chi_ungrouped_v3")
    progress["p_receipt_decoded"], progress["prior_receipt_decoded"] = True, True
    p_output, prior_output = _as_mapping(p_receipt["result"].get("output"), "P output"), _as_mapping(prior_receipt["result"].get("output"), "prior output")
    p_family, p_logs, p_atoms, p_index = _source_scope(p_output, "P")
    prior_family, prior_logs, prior_atoms, prior_index = _source_scope(prior_output, "prior")
    require(p_family == prior_family and _same_source_atoms(p_atoms, prior_atoms), "source family/log endpoint mismatch")
    p_manifest, prior_manifest = _source_rows(p_receipt, "P"), _source_rows(prior_receipt, "prior")
    require(p_manifest.get(FEEDBACK_RECEIPT) == FEEDBACK_RAW_SHA256 and prior_manifest.get(FEEDBACK_RECEIPT) == FEEDBACK_RAW_SHA256, "common FEEDBACK source-manifest pin")
    for key in ("physical_scale_reapplied", "physical_factors_applied", "physical_weight_reapplied", "radial_weight_reapplied", "geometry_reapplied", "ledger_212_changed", "ledger_appended", "ledger_written"):
        if key in prior_output:
            require(prior_output[key] is False, "prior scope changed: " + key)
    packet, provenance = _potential_terms(p_output)
    prior_reference = {"receipt": PRIOR_RECEIPT, "raw_sha256": RAW_SHA256[PRIOR_RECEIPT], "result_sha256": RESULT_SHA256[PRIOR_RECEIPT], "output_path": ["result", "output"], "input_index_path": ["result", "output", "input_index"] if prior_index else None, "input_index_retained": bool(prior_index), "formal_parent_expanded": False, "prior_212_terms_reapplied": False, "prior_six_channels_reapplied": False, "prior_error_data_deleted": False}
    leaf = _leaf_callable(leaf_module)
    progress["actual_leaf_build_entered"] = True
    actual = leaf(packet, caller_obligation={"authenticated_branches": list(BRANCHES), "upstream_source_authenticated": True, "condition": "P finite branches are raw-authenticated exact endpoint records"})
    progress["actual_leaf_build_calls"] = 1
    return {"schema": COMPONENT, "method": "raw_exact_binary_interval_endpoint_sign_classification_and_integer_count", "actual_classification_result": actual, "input_term_provenance": provenance, "input_packet": packet, "p_receipt_reference": {"path": P_RECEIPT, "raw_sha256": RAW_SHA256[P_RECEIPT], "result_sha256": RESULT_SHA256[P_RECEIPT], "potential_compensation_m15_path": ["result", "output", "potential_compensation_m15"]}, "prior_receipt_reference": prior_reference, "source_family": encode(p_family), "source_logs": encode(p_logs), "prior_source_logs": encode(prior_logs), "source_power_basis": list(POWER_BASIS), "P_input_index": encode(p_index), "prior_input_index": encode(prior_index), "source_atoms_exact_endpoint_comparison": encode(p_atoms), "term_count": len(provenance), "coefficient_addition_performed": False, "coefficient_arithmetic_performed": False, "midpoint_evaluated": False, "function_value_evaluated": False, "endpoint_classification_only": True, "integer_counts_only": True, "physical_weight_reapplied": False, "physical_scale_reapplied": False, "physical_measure_reapplied": False, "radial_weight_reapplied": False, "P1_added": False, "P_is_P1": False, "P1_is_m15": False, "prior_212_terms_reapplied": False, "prior_six_channels_reapplied": False, "all_prior_error_data_retained": True, "gauge_common_owner_proved": False, "live_owner_claimed": False, "terminal_closure_proved": False, "physical_construction_claimed": False, "physical_field_claimed": False, "complete_m15_computed": False, "ancestor_reexecuted": False, "actual_leaf_build_calls": progress["actual_leaf_build_calls"], "fixture_leaf_build_calls": 1}


def _save_outcome(path: Path, outcome: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    try:
        return save(path, outcome), outcome
    except BaseException as exc:
        fallback = dict(outcome)
        fallback["output"] = None
        fallback["status"] = "failed"
        fallback["serialization_complete"] = False
        fallback["error"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc(), "prior_error": fallback.get("error")}
        return save(path, fallback), fallback


def main() -> int:
    hashes, source_sha256, brief_sha256, preimage_sha256 = preflight()
    DIR.mkdir(parents=True, exist_ok=True)
    entry_sha256 = save(ENTRY, {"component": COMPONENT, "previous_component": PREVIOUS_COMPONENT, "source_sha256": source_sha256, "brief_sha256": brief_sha256, "preimage_sha256": preimage_sha256, "input_hashes": hashes, "entered": True, "exclusive_entry_and_result": True, "raw_gzip_authenticated_before_decode": True, "fixture_before_saved_receipt_decode": True, "failures_raw_only_not_decompressed": True, "ancestor_reexecuted": False})
    progress: dict[str, Any] = {"raw_gzip_authenticated_before_decode": False, "failures_raw_only_not_decompressed": True, "fixture_started": False, "fixture_completed_before_receipt_decode": False, "input_json_decode_started": False, "p_receipt_decoded": False, "prior_receipt_decoded": False, "actual_leaf_build_entered": False, "actual_leaf_build_calls": 0, "fixture": None}
    outcome: dict[str, Any] = {"component": COMPONENT, "previous_component": PREVIOUS_COMPONENT, "source_sha256": source_sha256, "entry_sha256": entry_sha256, "status": "failed", "raw_gzip_authenticated_before_decode": True, "fixture_before_saved_receipt_decode": True, "failures_raw_only_not_decompressed": True, "progress": progress, "output": None, "serialization_complete": False}
    try:
        outcome["output"] = consume(hashes, progress)
        outcome["status"], outcome["error"] = "returned", None
    except BaseException as exc:
        outcome["error"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
    outcome["serialization_complete"] = True
    result_sha256, outcome = _save_outcome(RESULT, outcome)
    print(json.dumps({"component": COMPONENT, "status": outcome["status"], "entry_sha256": entry_sha256, "result_sha256": result_sha256, "error": (outcome.get("error") or {}).get("message")}, sort_keys=True))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())

