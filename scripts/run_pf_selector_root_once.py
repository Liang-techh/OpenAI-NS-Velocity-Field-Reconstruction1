"""Freeze and enter the non-common P/F selector root source leaf once.

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
COMPONENT = "pf_selector_root_v1"
PHASE = "pf_selector_root_pending"
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
        "src/openai_ns_reconstruction/pf_selector_correlated_root.py",
        "src/openai_ns_reconstruction/pf_selector_quadratic_source.py",
        "src/openai_ns_reconstruction/f_historical_full_endpoint_roles.py",
        "scripts/run_pf_selector_root_once.py",
        "tmp/pf_noncommon_moment_v1/result.json",
        "tmp/pf_selector_bridge_v1/result.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_native_physical_pulse_source.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check.json",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_historical_frame_gap_saved_rewrite_lookup.json",
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
        algebra = load_new_module(ROOT / "src/openai_ns_reconstruction/pf_selector_correlated_root.py",
                                  "pf_selector_correlated_root_new_leaf")
        source_leaf = load_new_module(ROOT / "src/openai_ns_reconstruction/pf_selector_quadratic_source.py",
                                      "pf_selector_quadratic_source_new_leaf")
        endpoint_leaf = load_new_module(ROOT / "src/openai_ns_reconstruction/f_historical_full_endpoint_roles.py",
                                        "f_historical_full_endpoint_roles_new_leaf")
        calls.append("root_algebra_synthetic_fixture")
        algebra_fixture = algebra.verify_synthetic_fixture()
        require(algebra_fixture.get("checks_passed") is True, "selector root synthetic fixture failed")
        calls.append("source_mapping_synthetic_fixture")
        source_fixture = source_leaf.verify_mapping_synthetic()
        require(source_fixture.get("checks_passed") is True, "source mapping synthetic fixture failed")
        calls.append("endpoint_mapping_synthetic_controls")
        absent = endpoint_leaf._lookup_new_field(999999991, "lower", "R_q", 0, {})
        present = endpoint_leaf._lookup_new_field(999999992, "upper", "R_q", 1, {"999999992": 999999993})
        endpoint_checks = {
            "missing_not_identity": absent.get("mapped_ref") is None and absent.get("missing_mapping_is_not_identity") is True,
            "present_entry_retained": present.get("mapped_ref") == 999999993 and present.get("actual_F_field_claimed") is False,
        }
        require(all(endpoint_checks.values()), "endpoint synthetic controls failed")
        calls.append("selector_quadratic_direct_roles")
        roles = source_leaf.build_source_packet(ROOT, source)
        require(roles.get("source_values_computed") is False and roles.get("source_graph_traversed") is False
                and roles.get("difference_provider_admissible_for_numeric_use") is False,
                "selector source packet flags invalid")
        calls.append("historical_endpoint_additional_roles")
        endpoints = endpoint_leaf.build_historical_full_endpoint_roles(ROOT)
        require(endpoints.get("source_values_computed") is False and
                endpoints.get("source_graph_traversed") is False and
                endpoints.get("admissible_for_numeric_use") is False,
                "historical endpoint packet flags invalid")
        output = {"component": COMPONENT, "phase": PHASE, "status": "returned",
                  "source_checkpoint_sha256": checkpoint_sha, "entry_sha256": entry_sha,
                  "input_hashes": hashes, "consumer_calls": calls,
                  "ancestor_math_calls": 0, "source_values_computed": False,
                  "certified": False, "runtime_model_metadata": None,
                  "result": {"algebra_fixture": algebra_fixture, "source_fixture": source_fixture,
                             "endpoint_synthetic_checks": endpoint_checks,
                             "source_roles": roles, "endpoints": endpoints}, "error": None}
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
