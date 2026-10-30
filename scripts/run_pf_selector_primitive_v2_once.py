"""Freeze and enter the P/F selector primitive-difference leaf once.

No historical producer or evaluator is imported. Every declared input is
authenticated as raw bytes before entry and before any archive decompression.
An exclusive entry receipt precedes all imports and source admission, and a
success or failure result is persisted before it is printed.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "pf_selector_primitive_v2"
PHASE = "pf_selector_primitive_v2_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_once(path: Path, obj: dict) -> str:
    data = json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    with path.open("xb") as f:
        f.write(data)
    return hashlib.sha256(data).hexdigest()


def load_new_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "cannot load new leaf")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def preflight() -> tuple[dict, dict[str, str], str]:
    require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    require(brief.get("component") == COMPONENT and brief.get("phase") == PHASE,
            "brief not pending for this component")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False
            and brief.get("do_not_reenter") is True, "brief goal flags mismatch")
    checkpoint = brief.get("source_checkpoint", {})
    require(Path(checkpoint.get("path", "")).resolve() == SOURCE.resolve()
            and checkpoint.get("sha256") == sha(SOURCE), "checkpoint mismatch")
    require(source.get("component") == COMPONENT and source.get("one_entry_only") is True
            and source.get("no_source_graph_traversal") is True
            and source.get("all_raw_inputs_authenticated_before_decompression") is True,
            "checkpoint safeguards missing")
    files = source.get("files")
    require(isinstance(files, list) and files, "empty input manifest")
    hashes: dict[str, str] = {}
    for item in files:
        relative = item["path"]
        candidate = Path(relative)
        require(not candidate.is_absolute() and ".." not in candidate.parts,
                "manifest path must stay inside repository")
        path = (ROOT / candidate).resolve()
        require(path.is_relative_to(ROOT.resolve()) and path.is_file(), "invalid input path")
        require(relative not in hashes and item["sha256"] == sha(path), "raw SHA mismatch")
        hashes[relative] = item["sha256"]
    required = {
        "src/openai_ns_reconstruction/pf_selector_primitive_delta_v2.py",
        "src/openai_ns_reconstruction/pf_selector_primitive_source_v2.py",
        "scripts/run_pf_selector_primitive_v2_once.py",
        "tmp/pf_selector_root_v1/result.json",
        "tmp/pf_selector_primitive_v1/result.json",
        "tmp/pf_selector_bridge_v1/result.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json",
    }
    require(required.issubset(hashes), "required frozen inputs missing")
    return source, hashes, sha(SOURCE)


def main() -> int:
    source, hashes, checkpoint_sha = preflight()
    entry_sha = write_once(ENTRY, {"component": COMPONENT, "phase": PHASE,
        "source_checkpoint_sha256": checkpoint_sha, "input_hashes": hashes,
        "consumer_calls_started": 0, "ancestor_math_calls": 0,
        "source_values_computed": False, "runtime_model_metadata": None})
    calls = []
    try:
        prior_failure = json.loads((ROOT / "tmp/pf_selector_primitive_v1/result.json").read_bytes())
        require(prior_failure.get("component") == "pf_selector_primitive_v1" and
                prior_failure.get("status") == "failed" and
                prior_failure.get("consumer_calls") == ["primitive_algebra_synthetic_fixture"] and
                prior_failure.get("source_values_computed") is False,
                "prior failed entry scope mismatch")
        algebra = load_new_module(ROOT / "src/openai_ns_reconstruction/pf_selector_primitive_delta_v2.py",
                                  "pf_selector_primitive_delta_v2_new_leaf")
        source_leaf = load_new_module(ROOT / "src/openai_ns_reconstruction/pf_selector_primitive_source_v2.py",
                                      "pf_selector_primitive_source_v2_new_leaf")
        calls.append("primitive_algebra_synthetic_fixture")
        fixture = algebra.verify_synthetic_fixture()
        require(fixture.get("checks_passed") is True, "primitive algebra synthetic fixture failed")
        calls.append("primitive_mapping_synthetic_controls")
        mapping_fixture = source_leaf.verify_mapping_synthetic()
        require(mapping_fixture.get("checks_passed") is True, "primitive mapping controls failed")
        calls.append("actual_primitive_role_binding")
        packet = source_leaf.build_primitive_source(ROOT, source)
        require(packet.get("source_values_computed") is False and packet.get("source_graph_traversed") is False
                and packet.get("source_function_difference_provider_bound") is False
                and packet.get("numeric_use_admissible") is False, "source packet scope invalid")
        output = {"component": COMPONENT, "phase": PHASE, "status": "returned",
                  "source_checkpoint_sha256": checkpoint_sha, "entry_sha256": entry_sha,
                  "input_hashes": hashes, "consumer_calls": calls,
                  "ancestor_math_calls": 0, "source_values_computed": False,
                  "certified": False, "runtime_model_metadata": None,
                  "result": {"algebra_fixture": fixture, "mapping_fixture": mapping_fixture,
                             "primitive_roles": packet}, "error": None}
    except Exception as exc:
        output = {"component": COMPONENT, "phase": PHASE, "status": "failed",
                  "source_checkpoint_sha256": checkpoint_sha, "entry_sha256": entry_sha,
                  "input_hashes": hashes, "consumer_calls": calls,
                  "ancestor_math_calls": 0, "source_values_computed": False,
                  "certified": False, "runtime_model_metadata": None,
                  "error": {"type": type(exc).__name__, "message": str(exc)},
                  "traceback": traceback.format_exc()}
    result_sha = write_once(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"],
                      "result_sha256": result_sha, "consumer_calls": calls}, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
